# Evidence for the CHANGELOG 1.2.0 entry on the q3c v2.0.6 ellipse rewrite.
#
# Not a test: a one-off differential measurement, kept so its numbers can be
# reproduced. Two steps, each run with the qlsc build it names on PYTHONPATH:
#
#   pip install --no-deps --no-build-isolation --target t205 <qlsc at 5630261>/source
#   pip install --no-deps --no-build-isolation --target t206 <qlsc at b760357>/source
#   PYTHONPATH=t205 python q3c_206_ellipse_differential.py sample e205.npz
#   PYTHONPATH=t206 python q3c_206_ellipse_differential.py sample e206.npz
#   python q3c_206_ellipse_differential.py adjudicate e205.npz e206.npz   # needs mpmath
#
# Result on r740, 2026-11-22: 800,000 answers, 38,364 differ; against the
# tangent plane conic in 50-digit mpmath, v2.0.6 right 38,364 times, v2.0.5 none.

import sys
import numpy as np

N_ELLIPSES, N_POINTS, SEED = 4000, 200, 20261122

def sample_geometry():
	"""Yield (ra0, dec0, a, axis_ratio, pa, ra[], dec[]) for the fixed random sample."""
	rng = np.random.default_rng(SEED)
	for _ in range(N_ELLIPSES):
		ra0 = rng.uniform(0, 360); dec0 = np.degrees(np.arcsin(rng.uniform(-1, 1)))
		a = 10 ** rng.uniform(-7, np.log10(80))	# 1e-7 .. 80 degrees
		ar = rng.uniform(0.01, 0.999); pa = rng.uniform(0, 360)
		# points within twice the semi-major axis, so many lie near the edge
		sep = a * rng.uniform(0, 2, N_POINTS); ang = rng.uniform(0, 2*np.pi, N_POINTS)
		d0 = np.radians(dec0); s = np.radians(sep)
		dec = np.arcsin(np.clip(np.sin(d0)*np.cos(s) + np.cos(d0)*np.sin(s)*np.cos(ang), -1, 1))
		dra = np.arctan2(np.sin(ang)*np.sin(s)*np.cos(d0), np.cos(s) - np.sin(d0)*np.sin(dec))
		yield ra0, dec0, a, ar, pa, (np.degrees(dra) + ra0) % 360, np.degrees(dec)

def sample(out):
	from qlsc import QLSC, q3c
	q = QLSC(depth=30)
	res, meta = [], []
	for ra0, dec0, a, ar, pa, ra, dec in sample_geometry():
		res.append([bool(q.point_in_ellipse(ra=float(x), dec=float(y), ra_cen=ra0, dec_cen=dec0,
											semi_major_axis=a, axis_ratio=ar, position_angle=pa))
					for x, y in zip(ra, dec)])
		meta.append((ra0, dec0, a, ar, pa))
	np.savez(out, res=np.array(res), meta=np.array(meta), version=q3c.version())

def adjudicate(old_file, new_file):
	import mpmath as mp
	mp.mp.dps = 50
	old, new = np.load(old_file), np.load(new_file)
	print("versions:", old["version"], new["version"])
	r_old, r_new = old["res"], new["res"]
	total = right_old = right_new = 0
	by_size = {}
	D = mp.pi / 180
	for i, (ra0, dec0, a, ar, pa, ra, dec) in enumerate(sample_geometry()):
		for j in np.where(r_old[i] != r_new[i])[0]:
			A0, De0, PA, T = mp.mpf(ra0)*D, mp.mpf(dec0)*D, mp.mpf(pa)*D, mp.tan(mp.mpf(a)*D)
			A1, De1 = mp.mpf(float(ra[j]))*D, mp.mpf(float(dec[j]))*D
			cosc = mp.sin(De0)*mp.sin(De1) + mp.cos(De0)*mp.cos(De1)*mp.cos(A1-A0)
			x = mp.cos(De1)*mp.sin(A1-A0)/cosc
			y = (mp.cos(De0)*mp.sin(De1) - mp.sin(De0)*mp.cos(De1)*mp.cos(A1-A0))/cosc
			u = x*mp.sin(PA) + y*mp.cos(PA); v = -x*mp.cos(PA) + y*mp.sin(PA)
			inside = cosc > 0 and (u/T)**2 + (v/(T*ar))**2 < 1
			total += 1; right_old += bool(r_old[i, j]) == inside; right_new += bool(r_new[i, j]) == inside
			b = by_size.setdefault(int(np.floor(np.log10(a))), [0, 0, 0])
			b[0] += 1; b[1] += bool(r_old[i, j]) == inside; b[2] += bool(r_new[i, j]) == inside
	print(f"answers: {r_old.size}, differing: {total}, old right: {right_old}, new right: {right_new}")
	for k in sorted(by_size):
		print(f"  a ~ 1e{k} deg: n={by_size[k][0]} old right={by_size[k][1]} new right={by_size[k][2]}")

if __name__ == "__main__":
	if sys.argv[1] == "sample":
		sample(sys.argv[2])
	elif sys.argv[1] == "adjudicate":
		adjudicate(sys.argv[2], sys.argv[3])
	else:
		sys.exit("usage: sample OUT.npz | adjudicate OLD.npz NEW.npz")
