/*
  permdist : Distribution of Permutation Tests by Streitberg and Roehmel
  Copyright (C) 2000  Torsten Hothorn <Torsten.Hothorn@rzmail.uni-erlangen.de>
    
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or (at
  your option) any later version.
            
  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License for more details.
                   
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA.
                          
  SYNOPSIS
               
  void cpermdist1(double *x, int *score_a, int *N)
  void cpermdist2(double *x, int *m, int *c, int *score_a, int *score_b, int *N)
                                             
  DESCRIPTION
                                                  
  cpermdist1	The density of the permutation distribution for paired 
  		observations.
  cpermdist2	The density of the permutation distribution for 
  		independent observations   	
                                                                  
*/

                                                                   
#include <R.h>
#include <R_ext/Mathlib.h>

/*
	N = m + n  <= 200 only 
*/

#define PERM_MAX_N 200

void cpermdist1(double *x, int *score_a, int *N)
{
	/*
		compute the joint permutation distribution of the 
		sum of the elements of score_a
		(usualy score_a = Wilcoxon scores
		which leads to the exact conditional distribution 
		in the paired two-sample situation).
	*/ 

	double *H; 
	int i, k, sum_a = 0, s_a = 0;
	double msum = 0.0;

	if (*N > PERM_MAX_N)
		error("N > %d in cpermdistr1", PERM_MAX_N); 
	
	for (i = 0; i < *N; i++) 
		sum_a += score_a[i];

	/*
		initialize H
	*/

	H = (double *) calloc(sum_a + 1, sizeof(double));
	if (!H)
		error("cpermdist1 allocation error %d", 1);
	for (i = 0; i <= sum_a; i++)
		H[i] = 0;
		
	/*
		start the algorithm with H[0] = 1
	*/
		
	H[0] = 1;
	
	for (k = 0; k < *N; k++) {
		s_a = s_a + score_a[k];
	
		for (i = s_a; i >= score_a[k]; i--)
			H[i] = H[i] + H[i - score_a[k]];

	}


	/* 
		get the values in row m and sum it up
	*/

	for (i = 0; i <= sum_a; i++)
		msum += H[i];	/* 0 is a possible realization */
	
	/*
		compute probabilities
		note: x holds the probabilities and can be read from R.
		[dpq] stuff is done in R

	*/ 
	
	for (i = 0; i <= sum_a; i++)
		x[i] = H[i]/msum;	/* 0 is a possible realization */
	
	/*
		free memory and exit
	*/
	
	free((void *) H);

}

void cpermdist2(double *x, int *m, int *c, int *score_a, int *score_b, int *N, int *lenx)
{
	/*
		compute the joint permutation distribution of the 
		sum of the first m elements of score_a and score_b
		(usualy score_a = rep(1, N) and 
			score_b = Wilcoxon (Ansari...) scores
		which leads to the exact conditional distribution 
		in the two-sample situation).
	*/ 

	double **H; 
	int i, j, k, z, sum_a = 0, sum_b = 0, s_a = 0, s_b = 0;
	double msum = 0.0;

	if (*N > PERM_MAX_N)
		error("N > %d in cpermdistr2", PERM_MAX_N); 
	
	for (i = 0; i < *N; i++) {
		sum_a += score_a[i];
		sum_b += score_b[i];
	}

	/*
		optimization according to Streitberg & Roehmel
	*/
	
	sum_a = imin2(sum_a, *m);
	sum_b = imin2(sum_b, *c);

	/*
		initialize H
	*/

	H = (double **) calloc(sum_a + 1, sizeof(double *));
	if (!H)
		error("cpermdist2 allocation error %d", 1);
	for (i = 0; i <= sum_a; i++) {
		H[i] = (double *) calloc(sum_b + 1, sizeof(double));
		if (!H)
			error("cpermdist2 allocation error %d", 2);
		for (j = 0; j <= sum_b; j++)
			H[i][j] = 0;
	}
		
	/*
		start the algorithm with H[0][0] = 1
	*/
		
	H[0][0] = 1;
	
	for (k = 0; k < *N; k++) {
		s_a = s_a + score_a[k];
		s_b = s_b + score_b[k];
	
	/*
		compute H up to row m and column c
		note: 
			sum_a = min(sum_a, *m)
			sum_b = min(sum_b, *c)
		in cpermdist
	*/
		
		for (i = imin2(*m, s_a); i >= score_a[k]; i--) {
			for (j = imin2(*c,s_b); j >= score_b[k]; j--) {
				H[i][j] = H[i][j] + H[i - score_a[k]][j - score_b[k]];
			}
		}
	}

	/*
		return the hole matrix H (not needed within exactRankTests)
	*/ 

	if (*lenx > sum_b)
	{
		z = 0;
		for (k = 0; k < *N; k++) {
			for (j = 0; j < sum_b; j++) {
				x[z] = H[k][j];
				z++;
			}
		}
	} else {	

		/* 
			get the values in row m and sum it up
		*/

		for (j = 0; j < sum_b; j++)
		{
			x[j] = H[*m][j+1];
			msum += x[j];
		}
	
		/*
			compute probabilities
			note: x holds the probabilities and can be read from R.
			[dpq] stuff is done in R
		*/ 
	
		for (j = 0; j < sum_b; j++)
			x[j] = x[j]/msum;
	
	}
		
	/*
		free memory and exit
	*/
	
	for (i = sum_a; i >= 0; i--) 
		if (H[i] != 0) free((void *) H[i]);
	free((void *) H);

}

