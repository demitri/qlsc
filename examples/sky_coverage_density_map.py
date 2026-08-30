#!/usr/bin/env python

'''
This script demonstrates how to export the sky coverage of a QLSCIndex as JSON.

It is the qlsc half of the density-map demo: a catalog is indexed, the points are
counted per pixel at a coarse depth, and the result is written as a "sky-coverage"
version 1 document — the interchange format the astplot-js project's ``CoverageLayer``
reads, rendering the pixels as a filled choropleth in any projection.

The catalog used is ``1e5_wise.txt``, the 10,000 WISE sources that ship with these
examples. They cover the south polar cap from Dec -89.85 to -82.35 at every right
ascension, which makes the output a good test of a renderer as well as a demo: the
pixels wrap all the way around RA=0 and crowd together at the pole.

Usage:

	python sky_coverage_density_map.py                      # writes wise_coverage.json at depth 6
	python sky_coverage_density_map.py --depth 8 -o map.json

The document this writes can be handed straight to astplot-js:

	const doc = await (await fetch('wise_coverage.json')).json();
	plot.addLayer(new CoverageLayer(doc, { palette: 'viridis' }));
'''

import json
import argparse
import pathlib

import numpy as np

from qlsc import QLSC, QLSCIndex

INPUT_FILE = pathlib.Path(__file__).parent / "1e5_wise.txt"

def main():

	parser = argparse.ArgumentParser(description="Export the sky coverage of a catalog as sky-coverage v1 JSON.")
	parser.add_argument("--depth", type=int, default=6,
						help="the QLSC depth to count points at; larger is finer (default: 6)")
	parser.add_argument("-o", "--output", type=pathlib.Path, default=pathlib.Path("wise_coverage.json"),
						help="where to write the JSON document (default: wise_coverage.json)")
	args = parser.parse_args()

	# --- read the catalog -------------------------------------------------------------
	# columns: cntr, ra, dec (tab separated, one header line)
	catalog = np.loadtxt(INPUT_FILE, skiprows=1, delimiter="\t")
	ra, dec = catalog[:,1], catalog[:,2]
	print(f"Read {len(ra):,} sources from {INPUT_FILE.name}.")

	# --- index them -------------------------------------------------------------------
	# The index is built at depth 30, the finest available, so that the same index could
	# answer cone searches; the coverage below is counted at a much coarser depth. The
	# index does not have to be rebuilt to change the depth of the map.
	index = QLSCIndex(qlsc=QLSC(depth=30))
	index.add_points(ra=ra, dec=dec)
	print(f"Indexed {index.number_of_points:,} points (duplicate coordinates are discarded).")

	# --- export the coverage ----------------------------------------------------------
	# One GROUP BY inside SQLite does the counting; every pixel holding at least one point
	# becomes a polygon carrying its count. "meta" is the format's escape hatch for
	# provenance, which the renderer never interprets.
	document = index.coverage(depth=args.depth, meta={"catalog": INPUT_FILE.name})

	counts = [entry["value"] for entry in document["polygons"]]
	print(f"Depth {args.depth}: {len(counts):,} occupied pixels, "
		  f"{min(counts):.0f} to {max(counts):.0f} points each, {sum(counts):.0f} in total.")

	# The document is a plain dictionary of plain Python types, so this is the whole
	# hand-off: the standard library encoder, with no custom default= hook. allow_nan is
	# turned off because NaN and Infinity are not JSON and the format rejects them.
	with open(args.output, "w") as f:
		json.dump(document, f, allow_nan=False)

	print(f"Wrote {args.output} ({args.output.stat().st_size/1024:.0f} kB).")

	index.close()

if __name__ == "__main__":
	main()
