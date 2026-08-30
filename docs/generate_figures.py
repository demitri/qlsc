#!/usr/bin/env python
'''
Generate the figures used by the QLSC documentation.

Every figure in ``docs/_static/figures`` is produced by this script; none are
drawn by hand. The catalogues are synthetic and generated from a fixed random
seed, so the numbers quoted in the surrounding documentation pages match what
this script draws.

Usage::

	python generate_figures.py [--outdir _static/figures] [name ...]

Requires matplotlib and astropy in addition to qlsc's own dependencies; neither
is needed to build the documentation itself, only to regenerate these files.
'''

import argparse
import os
import sys

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
from matplotlib.patches import Circle, Ellipse, Polygon as MPLPolygon

from qlsc import QLSC, QLSCIndex, distance

# A colourblind-safe categorical set (Okabe & Ito); marker shape carries the
# same distinction so the figures do not rely on colour alone.
BLUE = "#0072B2"
VERMILLION = "#D55E00"
GREEN = "#009E73"
GREY = "#9a9a9a"

# A single perceptually uniform, colourblind-safe ramp for every density map.
DENSITY_CMAP = "cividis"

SEED = 20200422

SQ_DEG_PER_STERADIAN = (180. / np.pi) ** 2


# ---------------------------------------------------------------------------
# Pixel areas.
#
# QLSC pixels are only approximately equal in area: at any depth they vary by a
# factor of about 5.2, from the pixels at a cube VERTEX (the smallest) to those
# at a face centre (the largest). Pixel 0 is a cube-vertex pixel, so
# ``ipix_area(0, depth)`` is the SMALLEST area in the scheme and about 2.55
# times below the mean - which is what these figures used to divide their
# counts by, overstating every density by that factor while the surrounding
# prose said the mean had been used.
#
# There are two defensible scales and this file uses both, deliberately and by
# name: the mean, for anything characterising the scheme as a whole (a pixel
# side against depth), and the pixel's own area, for anything showing a real
# density. Neither is ``ipix_area(0, ...)``, and test_figure_pixel_areas.py
# fails if that spelling comes back.
# ---------------------------------------------------------------------------

def mean_pixel_area_sq_deg(q):
	'''
	The mean area of a pixel of a scheme, in square degrees: the sphere divided by
	the number of pixels. This is the scale that characterises a depth.
	'''
	return 4. * np.pi / q.nbins * SQ_DEG_PER_STERADIAN


def pixel_areas_sq_deg(q, ipix):
	'''
	The area of each of the given pixels, in square degrees.

	This is the scale that turns a count into a density, because the count came
	from that pixel and no other.
	'''
	return np.array([q.ipix_area(int(p), q.depth) for p in np.atleast_1d(ipix)]) \
		* SQ_DEG_PER_STERADIAN


def _style(ax, xlabel="right ascension (deg)", ylabel="declination (deg)"):
	ax.set_xlabel(xlabel)
	ax.set_ylabel(ylabel)
	ax.grid(True, color="#e6e6e6", linewidth=0.6, zorder=0)
	ax.set_axisbelow(True)
	for side in ("top", "right"):
		ax.spines[side].set_visible(False)
	for side in ("left", "bottom"):
		ax.spines[side].set_color("#b0b0b0")
	ax.tick_params(colors="#555555", labelsize=8)
	# right ascension labels are large numbers with small ranges; matplotlib's
	# default offset notation ("+2.143e2") is unreadable for them
	ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(5))
	ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter(useOffset=False))


def _legend(ax, **kwargs):
	'''A legend that stays readable when it sits over the data.'''
	kwargs.setdefault("loc", "lower left")
	kwargs.setdefault("fontsize", 8)
	legend = ax.legend(frameon=True, framealpha=0.88, edgecolor="none", **kwargs)
	legend.get_frame().set_facecolor("white")
	legend.set_zorder(10)
	return legend


def _sky_circle(ra, dec, radius_deg, **kwargs):
	'''
	A circle of the given angular radius drawn in (ra,dec) axes.

	The ra axis is compressed by cos(dec), so a circle on the sky is an ellipse
	in these coordinates; drawing a matplotlib Circle would misstate its size.
	'''
	sec = 1. / np.cos(np.radians(dec))
	return Ellipse((ra, dec), width=2 * radius_deg * sec, height=2 * radius_deg,
				   fill=False, **kwargs)


# ---------------------------------------------------------------------------
# 1. cone search around a transient  (docs/examples/cone_search.rst)
# ---------------------------------------------------------------------------

def transient_catalogue():
	rng = np.random.default_rng(seed=SEED)
	n_sky = 500_000
	sky_ra = rng.uniform(0., 360., n_sky)
	sky_dec = np.degrees(np.arcsin(rng.uniform(-1., 1., n_sky)))

	alert_ra, alert_dec = 214.36719, -12.83402
	n_field = 40_000
	r = 0.5 * np.sqrt(rng.uniform(0., 1., n_field))
	t = rng.uniform(0., 2. * np.pi, n_field)
	field_dec = alert_dec + r * np.sin(t)
	field_ra = alert_ra + r * np.cos(t) / np.cos(np.radians(alert_dec))

	ra = np.concatenate([sky_ra, field_ra])
	dec = np.concatenate([sky_dec, field_dec])
	return ra, dec, alert_ra, alert_dec


def figure_cone_search(path):
	ra, dec, alert_ra, alert_dec = transient_catalogue()
	index = QLSCIndex()
	index.add_points(ra=ra, dec=dec)

	radius_deg = 30. / 3600.
	matches = index.radial_query(ra=alert_ra, dec=alert_dec, radius=radius_deg,
								 return_distance=True)

	fig, (left, right) = plt.subplots(1, 2, figsize=(9.4, 4.3))

	# ---- the field ----
	sec = 1. / np.cos(np.radians(alert_dec))
	near = index.radial_query(ra=alert_ra, dec=alert_dec, radius=0.55)
	left.scatter(near[:, 0], near[:, 1], s=1.2, c=GREY, linewidths=0, zorder=2)
	left.add_patch(_sky_circle(alert_ra, alert_dec, 0.5, edgecolor=BLUE,
							   linewidth=1.6, zorder=3))
	left.plot([alert_ra], [alert_dec], marker="+", color=VERMILLION,
			  markersize=14, markeredgewidth=1.8, zorder=4)
	left.set_xlim(alert_ra + 0.62 * sec, alert_ra - 0.62 * sec)
	left.set_ylim(alert_dec - 0.62, alert_dec + 0.62)
	left.set_title(f"the field: {len(near):,} indexed sources", fontsize=10)
	_style(left)

	# ---- the cone, with the depth-14 pixels that cover it ----
	depth = 14
	q = QLSC(depth=depth)
	centre_pixel = q.ang2ipix(alert_ra, alert_dec)
	pixels = {int(centre_pixel)}
	frontier = set(pixels)
	for _ in range(3):
		new = set()
		for p in frontier:
			new.update(int(x) for x in q.ipix_neighbors(p))
		frontier = new - pixels
		pixels |= frontier
	for p in sorted(pixels):
		poly = q.ipix2polygon(ipix=p)
		right.add_patch(MPLPolygon(poly, closed=True, fill=False,
								   edgecolor="#cccccc", linewidth=0.6, zorder=1))

	window = index.radial_query(ra=alert_ra, dec=alert_dec, radius=90. / 3600.)
	right.scatter(window[:, 0], window[:, 1], s=26, c=GREY, marker="o",
				  linewidths=0, zorder=2, label="indexed source")
	right.scatter(matches["ra"], matches["dec"], s=52, marker="^", c=VERMILLION,
				  edgecolors="white", linewidths=0.8, zorder=4,
				  label=f"within 30\" ({len(matches)})")
	right.add_patch(_sky_circle(alert_ra, alert_dec, radius_deg, edgecolor=BLUE,
								linewidth=1.6, zorder=3))
	right.plot([alert_ra], [alert_dec], marker="+", color=BLUE, markersize=12,
			   markeredgewidth=1.6, zorder=5)
	span = 75. / 3600.
	right.set_xlim(alert_ra + span * sec, alert_ra - span * sec)
	right.set_ylim(alert_dec - span, alert_dec + span)
	right.set_title(f"the cone and the depth-{depth} pixels covering it", fontsize=10)
	_legend(right)
	_style(right)

	fig.tight_layout()
	fig.savefig(path, dpi=110, facecolor="white")
	plt.close(fig)


# ---------------------------------------------------------------------------
# 2. survey footprint  (docs/examples/survey_footprint.rst)
# ---------------------------------------------------------------------------

def figure_survey_footprint(path):
	rng = np.random.default_rng(seed=SEED)
	ra_c, dec_c, half = 150.0, 2.0, 2.5
	sec = 1. / np.cos(np.radians(dec_c))

	# a local patch, dense enough to see the boundary
	n = 6000
	ra = rng.uniform(ra_c - 4. * sec, ra_c + 4. * sec, n)
	dec = rng.uniform(dec_c - 4., dec_c + 4., n)
	index = QLSCIndex()
	index.add_points(ra=ra, dec=dec)

	tile = np.array([
		[ra_c - half * sec, dec_c - half],
		[ra_c + half * sec, dec_c - half],
		[ra_c + half * sec, dec_c + half],
		[ra_c - half * sec, dec_c + half]])
	notched = np.array([
		[ra_c - half * sec, dec_c - half],
		[ra_c + half * sec, dec_c - half],
		[ra_c + half * sec, dec_c + half],
		[ra_c, dec_c + half],
		[ra_c, dec_c],
		[ra_c - half * sec, dec_c]])

	fig, axes = plt.subplots(1, 2, figsize=(9.4, 4.4))
	for ax, poly, title in ((axes[0], tile, "convex tile"),
							(axes[1], notched, "concave footprint")):
		inside = index.polygon_query(polygon=poly)
		hit = set(zip(np.round(inside[:, 0], 9), np.round(inside[:, 1], 9)))
		is_in = np.array([(round(a, 9), round(b, 9)) in hit for a, b in zip(ra, dec)])
		ax.scatter(ra[~is_in], dec[~is_in], s=3, c=GREY, marker="o",
				   linewidths=0, zorder=2, label=f"outside ({(~is_in).sum()})")
		ax.scatter(ra[is_in], dec[is_in], s=7, c=VERMILLION, marker="^",
				   linewidths=0, zorder=3, label=f"inside ({is_in.sum()})")
		ax.add_patch(MPLPolygon(poly, closed=True, fill=False, edgecolor=BLUE,
								linewidth=1.8, zorder=4))
		ax.set_xlim(ra_c + 4.2 * sec, ra_c - 4.2 * sec)
		ax.set_ylim(dec_c - 4.2, dec_c + 4.2)
		ax.set_title(f"{title}: polygon_query", fontsize=10)
		_legend(ax)
		_style(ax)

	fig.tight_layout()
	fig.savefig(path, dpi=110, facecolor="white")
	plt.close(fig)


# ---------------------------------------------------------------------------
# 3. elliptical aperture  (docs/examples/galaxy_aperture.rst)
# ---------------------------------------------------------------------------

def galaxy_field():
	rng = np.random.default_rng(seed=SEED)
	gal_ra, gal_dec = 187.70593, 12.39112
	a, q_ratio, pa_deg = 0.25, 0.30, 55.0

	n_disc = 300
	t = rng.uniform(0., 2 * np.pi, n_disc)
	s = np.sqrt(rng.uniform(0., 1., n_disc))
	maj = a * s * np.cos(t)
	mnr = a * q_ratio * s * np.sin(t)
	pa = np.radians(pa_deg)
	disc_dec = gal_dec + maj * np.cos(pa) - mnr * np.sin(pa)
	disc_ra = gal_ra + (maj * np.sin(pa) + mnr * np.cos(pa)) / np.cos(np.radians(gal_dec))

	n_bg = 20_000
	b_r = 1.0 * np.sqrt(rng.uniform(0., 1., n_bg))
	b_t = rng.uniform(0., 2 * np.pi, n_bg)
	bg_dec = gal_dec + b_r * np.sin(b_t)
	bg_ra = gal_ra + b_r * np.cos(b_t) / np.cos(np.radians(gal_dec))

	return (gal_ra, gal_dec, a, q_ratio, pa_deg,
			np.concatenate([disc_ra, bg_ra]), np.concatenate([disc_dec, bg_dec]),
			n_disc)


def figure_galaxy_aperture(path):
	gal_ra, gal_dec, a, q_ratio, pa_deg, ra, dec, n_disc = galaxy_field()
	keys = [f"DISC{i:03d}" for i in range(n_disc)] + \
		   [f"BG{i:05d}" for i in range(len(ra) - n_disc)]
	index = QLSCIndex()
	index.add_points(ra=ra, dec=dec, keys=keys)

	ell = index.ellipse_query(ra=gal_ra, dec=gal_dec, semi_major_axis=a,
							  axis_ratio=q_ratio, position_angle=pa_deg,
							  return_key=True)
	cir = index.radial_query(ra=gal_ra, dec=gal_dec, radius=a, return_key=True)
	n_ell_disc = sum(1 for k in ell["key"] if k.startswith("DISC"))
	n_cir_disc = sum(1 for k in cir["key"] if k.startswith("DISC"))

	sec = 1. / np.cos(np.radians(gal_dec))
	fig, ax = plt.subplots(figsize=(6.4, 6.0))
	ax.scatter(ra[n_disc:], dec[n_disc:], s=3, c=GREY, marker="o",
			   linewidths=0, zorder=2, label="field")
	ax.scatter(ra[:n_disc], dec[:n_disc], s=12, c=VERMILLION, marker="^",
			   linewidths=0, zorder=3, label="disc population")
	ax.add_patch(_sky_circle(gal_ra, gal_dec, a, edgecolor=GREEN, linewidth=1.8,
							 linestyle="--", zorder=4,
							 label=f"circle: {len(cir)} objects, {n_cir_disc} from the disc"))
	ax.add_patch(Ellipse((gal_ra, gal_dec), width=2 * a * q_ratio * sec,
						 height=2 * a, angle=-pa_deg, fill=False,
						 edgecolor=BLUE, linewidth=1.8, zorder=5,
						 label=f"ellipse: {len(ell)} objects, {n_ell_disc} from the disc"))
	span = 0.34
	ax.set_xlim(gal_ra + span * sec, gal_ra - span * sec)
	ax.set_ylim(gal_dec - span, gal_dec + span)
	ax.set_title("an inclined disc: elliptical vs circular aperture", fontsize=10)
	_legend(ax, loc="lower center")
	_style(ax)

	fig.tight_layout()
	fig.savefig(path, dpi=110, facecolor="white")
	plt.close(fig)


# ---------------------------------------------------------------------------
# 4 & 5. sky density  (docs/examples/density_map.rst)
# ---------------------------------------------------------------------------

def star_catalogue():
	from astropy.coordinates import SkyCoord
	from astropy import units as u

	rng = np.random.default_rng(seed=SEED)
	n_disc, n_bulge, n_halo = 1_500_000, 300_000, 200_000

	l_disc = rng.uniform(0., 360., n_disc)
	b_disc = rng.laplace(0., 6., n_disc)
	keep = np.abs(b_disc) < 90.
	l_disc, b_disc = l_disc[keep], b_disc[keep]

	l_bulge = rng.normal(0., 12., n_bulge) % 360.
	b_bulge = rng.normal(0., 7., n_bulge)
	keep = np.abs(b_bulge) < 90.
	l_bulge, b_bulge = l_bulge[keep], b_bulge[keep]

	l_halo = rng.uniform(0., 360., n_halo)
	b_halo = np.degrees(np.arcsin(rng.uniform(-1., 1., n_halo)))

	icrs = SkyCoord(l=np.concatenate([l_disc, l_bulge, l_halo]) * u.deg,
					b=np.concatenate([b_disc, b_bulge, b_halo]) * u.deg,
					frame="galactic").icrs
	return icrs.ra.deg, icrs.dec.deg


def figure_sky_density(path, depth=5):
	ra, dec = star_catalogue()
	q = QLSC(depth=depth)
	occupied, counts = np.unique(q.ang2ipix(ra=ra, dec=dec), return_counts=True)
	count_of = dict(zip(occupied.tolist(), counts.tolist()))

	# sample the count map on a plate carree grid rather than drawing the pixel
	# polygons: the polygons on the polar faces span most of the ra axis.
	n_ra, n_dec = 1440, 720
	grid_ra = np.linspace(0., 360., n_ra, endpoint=False) + 360. / (2 * n_ra)
	grid_dec = np.degrees(np.arcsin(np.linspace(-1., 1., n_dec, endpoint=False)
									+ 1. / n_dec))
	mesh_ra, mesh_dec = np.meshgrid(grid_ra, grid_dec)
	ipix = q.ang2ipix(ra=mesh_ra.ravel().copy(), dec=mesh_dec.ravel().copy())

	# Each pixel's count is divided by THAT pixel's area, so the map is a real
	# density and not a count map with a single scale factor on it. Only the
	# occupied pixels need an area, and each is computed once however many grid
	# cells fall in it.
	areas = dict(zip(occupied.tolist(), pixel_areas_sq_deg(q, occupied).tolist()))
	density_of = {p: count_of[p] / areas[p] for p in count_of}

	density = np.array([density_of.get(int(p), 0.) for p in ipix]).reshape(mesh_dec.shape)
	mean_area = mean_pixel_area_sq_deg(q)

	fig, ax = plt.subplots(figsize=(9.4, 4.8))
	mesh = ax.pcolormesh(grid_ra, np.degrees(np.arcsin(
		np.linspace(-1., 1., n_dec, endpoint=False) + 1. / n_dec)),
		np.log10(np.clip(density, 0.5, None)), cmap=DENSITY_CMAP, shading="nearest")
	bar = fig.colorbar(mesh, ax=ax, pad=0.02)
	bar.set_label("log$_{10}$ sources per square degree", fontsize=9)
	bar.ax.tick_params(labelsize=8)
	ax.set_xlim(360., 0.)
	ax.set_ylim(-90., 90.)
	ax.set_title(f"source density at depth {depth} ({q.nbins:,} pixels, "
				 f"{mean_area:.3g} sq deg on average), each pixel's count "
				 f"divided by its own area", fontsize=10)
	_style(ax)
	ax.grid(False)
	ax.set_xticks([0, 60, 120, 180, 240, 300, 360])
	ax.set_yticks([-90, -60, -30, 0, 30, 60, 90])

	fig.tight_layout()
	fig.savefig(path, dpi=110, facecolor="white")
	plt.close(fig)


def figure_cube_faces(path, depth=5):
	ra, dec = star_catalogue()
	q = QLSC(depth=depth)
	occupied, counts = np.unique(q.ang2ipix(ra=ra, dec=dec), return_counts=True)
	facenum, x, y = q.ipix2xy(occupied)

	# as in figure_sky_density: each count divided by its own pixel's area
	densities = counts / pixel_areas_sq_deg(q, occupied)
	side = 2. / q.nside
	vmin = np.log10(max(densities.min(), 1e-3))
	vmax = np.log10(densities.max())

	# the conventional unfolded-cube layout: the four side faces in a row with
	# the top and bottom faces above and below face 1.
	slots = {0: (0, 1), 1: (1, 1), 2: (1, 2), 3: (1, 3), 4: (1, 0), 5: (2, 1)}
	fig, axes = plt.subplots(3, 4, figsize=(9.0, 6.4))
	for ax in axes.ravel():
		ax.set_visible(False)

	for face, (row, col) in slots.items():
		ax = axes[row][col]
		ax.set_visible(True)
		sel = facenum == face
		grid = np.full((q.nside, q.nside), np.nan)
		i = np.round((x[sel] + 1.) / side).astype(int)
		j = np.round((y[sel] + 1.) / side).astype(int)
		grid[j, i] = np.log10(densities[sel])
		mesh = ax.imshow(grid, origin="lower", extent=(-1, 1, -1, 1),
						 cmap=DENSITY_CMAP, vmin=vmin, vmax=vmax,
						 interpolation="nearest")
		ax.set_title(f"face {face}", fontsize=9)
		ax.set_xticks([-1, 0, 1])
		ax.set_yticks([-1, 0, 1])
		ax.tick_params(labelsize=7, colors="#555555")
		for side_name in ("top", "right", "left", "bottom"):
			ax.spines[side_name].set_color("#b0b0b0")

	fig.subplots_adjust(left=0.04, right=0.86, top=0.90, bottom=0.05,
						wspace=0.35, hspace=0.30)
	bar = fig.colorbar(mesh, ax=axes, fraction=0.03, pad=0.02)
	bar.set_label("log$_{10}$ sources per square degree", fontsize=9)
	bar.ax.tick_params(labelsize=8)
	fig.suptitle(f"the same density on the six cube faces, in face (x,y) "
				 f"coordinates (depth {depth})", fontsize=10)
	fig.savefig(path, dpi=110, facecolor="white")
	plt.close(fig)


# ---------------------------------------------------------------------------
# 6. pixel size against depth  (docs/user_guide/depth.rst)
# ---------------------------------------------------------------------------

def figure_pixel_scale(path):
	depths = np.arange(0, 31)
	# the MEAN pixel side, sqrt(4*pi/nbins), which is the number
	# docs/user_guide/depth.rst tabulates and the code block on that page prints;
	# ipix_area(0, d) is the smallest pixel in the scheme, about 1.6 times narrower
	sides = np.array([np.sqrt(4. * np.pi / QLSC(depth=int(d)).nbins) * 180. / np.pi
					  for d in depths])

	fig, ax = plt.subplots(figsize=(7.4, 4.4))
	ax.plot(depths, sides * 3600., color=BLUE, linewidth=2, marker="o",
			markersize=4, zorder=3)
	ax.set_yscale("log")
	for label, scale_arcsec, colour in [("1 degree", 3600., GREEN),
										("1 arcmin", 60., GREEN),
										("1 arcsec", 1., GREEN),
										("1 mas", 1e-3, GREEN)]:
		ax.axhline(scale_arcsec, color=colour, linewidth=1, linestyle=":", zorder=2)
		ax.annotate(label, xy=(30.4, scale_arcsec), fontsize=8, color="#555555",
					va="center", ha="left", annotation_clip=False)
	ax.set_xlim(-0.5, 30.5)
	ax.set_title("mean pixel side length against depth", fontsize=10)
	_style(ax, xlabel="depth", ylabel="mean pixel side (arcsec)")
	ax.set_xticks(np.arange(0, 31, 2))
	fig.tight_layout()
	fig.savefig(path, dpi=110, facecolor="white")
	plt.close(fig)


FIGURES = {
	"cone_search": figure_cone_search,
	"survey_footprint": figure_survey_footprint,
	"galaxy_aperture": figure_galaxy_aperture,
	"sky_density": figure_sky_density,
	"cube_faces": figure_cube_faces,
	"pixel_scale": figure_pixel_scale,
}


def main():
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("names", nargs="*", default=None,
						help="figures to generate; the default is all of them")
	parser.add_argument("--outdir", default=os.path.join(os.path.dirname(__file__),
														 "_static", "figures"))
	args = parser.parse_args()

	names = args.names or sorted(FIGURES)
	unknown = set(names) - set(FIGURES)
	if unknown:
		parser.error(f"unknown figure(s): {', '.join(sorted(unknown))}; "
					 f"choose from {', '.join(sorted(FIGURES))}")

	os.makedirs(args.outdir, exist_ok=True)
	for name in names:
		path = os.path.join(args.outdir, f"{name}.png")
		print(f"{name} -> {path}", flush=True)
		FIGURES[name](path)


if __name__ == "__main__":
	sys.exit(main())
