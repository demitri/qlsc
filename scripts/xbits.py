#!/usr/bin/env python

'''
A script to show how xbits and ybits are populated in Q3C.
'''


import numpy as np
np.set_printoptions(suppress=True) #, precision=10) # numpy "pretty print"

Q3C_INTERLEAVED_NBITS = 16
nbits = Q3C_INTERLEAVED_NBITS

xybits_size = 1 << nbits
xbits = np.zeros((xybits_size), dtype=np.int64)
ybits = np.zeros((xybits_size), dtype=np.int64)
xbits1 = np.zeros((xybits_size), dtype=np.int64)
ybits1 = np.zeros((xybits_size), dtype=np.int64)

# xbits & ybits
# -------------
xbits[0] = 0;
xbits[1] = 1;
ybits[0] = 0;
ybits[1] = 2;

m = 1
for i in range(2, xybits_size):

	k = i / m;
	if (k == 2):
		xbits[i] = xbits[i // 2] * 4;
		ybits[i] = 2 * xbits[i];
		m *= 2;
	else:
		xbits[i] = xbits[m] + xbits[i % m];
		ybits[i] = 2 * xbits[i];


# xbits1 & ybits1
# ---------------
xbits1[0] = 0;
xbits1[1] = 1;

m = 2
l = 2
for i in range(2, xybits_size):
	k = i // m;

	if k < 2:
		xbits1[i] = xbits1[i - m];
	else:
		if k == 4:
			xbits1[i] = xbits1[0];
			m *= 4;
			l *= 2;
		else:
			xbits1[i] = xbits1[i - 2 * m] + l;

ybits1[0] = 0
ybits1[1] = 0

m = 1
l = 1

for i in range(2, xybits_size):
	k = i // m;

	if k < 2:
		ybits1[i] = ybits1[i - m]
	else:
		if k == 4:
			ybits1[i] = ybits1[0]
			m *= 4
			l *= 2
		else:
			ybits1[i] = ybits1[i - 2 * m] + l


for x in ybits1:
	print(f"{x:10} = {x:064b}")