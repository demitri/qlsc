//
//  main.c
//  q3c_test
//
//  Created by Demitri Muna on 12/15/19.
//  Copyright © 2019 Demitri Muna. All rights reserved.
//

#define Q3C_VERSION "1.8.1"
#include <stdio.h>
#include <float.h> // ref: https://stackoverflow.com/questions/16839658/printf-width-specifier-to-maintain-precision-of-floating-point-value
#include "common.h"

int main(int argc, const char * argv[]) {

	struct q3c_prm *hprm = malloc(sizeof (struct q3c_prm));
	q3c_coord_t ra, dec; // (double)
	q3c_ipix_t ipix;
	
	q3c_ipix_t nside = 1073741824; // number of quadtree subdivisions
	
	init_q3c1(hprm, nside);
	
	ra = 134;
	dec = -66;
	
	//printf("%d", Q3C_I1);

	q3c_ang2ipix(hprm, ra, dec, &ipix);
	printf("ipix = %" PRId64 "\n", ipix); // int64_t format
	
	double pixel_area = q3c_pixarea(hprm, ipix, 15);
	printf("pixel area = %.10e\n", pixel_area);
	
	char facenum = q3c_get_facenum(ra, dec);
	printf("facenum for (%.2f, %.2f): %d\n", ra, dec, facenum);
	
	double ra_out, dec_out;
	q3c_ipix2ang(hprm, ipix, &ra_out, &dec_out);
	fprintf(stdout, "ipix2ang: %.3f, %.3f", ra_out, dec_out);
	
	printf("\n\n");
	free(hprm);
	
	// ----

	// dist calculations for testing
	{
		q3c_coord_t ra1, dec1, ra2, dec2, distance; // all degrees
		ra1 = 15;
		dec1 = 88;
		ra2 = 15;
		dec2 = -88;
		distance = q3c_dist(ra1, dec1, ra2, dec2);
		printf("dist: %.4f, %.4f, %.4f, %.4f, %.*f\n", ra1, dec1, ra2, dec2, DECIMAL_DIG, distance);
	}

	// sindist calculations for testing
	{
		q3c_coord_t ra1, dec1, ra2, dec2, distance; // all degrees
		ra1 = 45.45;
		dec1 = 48.8;
		ra2 = 127.89;
		dec2 = -27.7;
		distance = q3c_sindist(ra1, dec1, ra2, dec2);
		printf("sindist: %.4f, %.4f, %.4f, %.4f, %.*f\n", ra1, dec1, ra2, dec2, DECIMAL_DIG, distance);
	}
	
	{
		q3c_coord_t x = 1.01;	// -1 <= x <= 1
		q3c_coord_t y = 1;	// -1 <= y <= 1
		char facenum0 = 0;
		char result = q3c_xy2facenum(x, y, facenum0);
		printf("facenum for (x,y) = (%.4f, %.4f) at face_num0 %d: %hhd\n", x, y, facenum0, result);
	}
	
	printf("\n\n");
	return 0;

}
