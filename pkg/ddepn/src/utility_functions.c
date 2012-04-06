//============================================================================
// Name        : utility_functions.c.c
// Author      : Christian Bender
// Version     :
// Copyright   : none
// Description : serial c-version of hmmsearch, using global theta estimation
//               adapted from function perform.hmmsearch_globalest in perform.hmmsearch.R
//
//============================================================================

#include <stdio.h>
#include <stdlib.h>  // for random numbers
#include <time.h>    // for seed for random numbers
#include <unistd.h>     // need for getpid() function
#include <math.h>    // for pow function
#include <string.h>  // for memcpy function
#include <R.h> // R definitions and macros, e.g. ISNA(x)
#include "utility_functions.h"

/*
 * estimate the normal density parameters
 */
void
estimate_theta(const double *X, int *GSsub, double *TH, const int N,
		const int T, const int R)
{
	int a_n, p_n, kk;
	double a, p, mua, mup, sda, sdp, tmp;
	double *params = malloc(4 * sizeof(double));
	double tmpval;

	// go through lines, i.e. proteins
	for (int i = 0; i != N; ++i)
	{
		a = 0.0;
		p = 0.0;
		a_n = 0;
		p_n = 0;
		// means; go through columns, i.e. timepoints
		for (int j = 0; j != T * R; ++j)
		{
			kk = i + j * N;
			tmpval = X[kk];
			if(!ISNA(tmpval)) {
                          if (GSsub[kk] == 1)
                          {
                                  a += tmpval;
                                  a_n++;
                          }
                          else
                          {
                                  p += tmpval;
                                  p_n++;
                          }
			}
		}

		mua, mup;
		if (a_n > 0)
			mua = a / a_n;
		else
			mua = HUGE_VAL;
		if (p_n > 0)
			mup = p / p_n;
		else
			mup = HUGE_VAL;

		// sds;
		a_n = 0;
		p_n = 0;
		a = 0.0;
		p = 0.0;
		for (int j = 0; j != T * R; ++j)
		{
			kk = i + j * N;
                        tmpval = X[kk];
                        if(!ISNA(tmpval)) {
                          if (GSsub[kk] == 1)
                          {
                                  tmp = tmpval - mua;
                                  a += (tmp * tmp);
                                  a_n++;
                          }
                          else
                          {
                                  tmp = tmpval - mup;
                                  p += (tmp * tmp);
                                  p_n++;
                          }
                        }
		}
		sda, sdp;
		a_n = a_n - 1;
		p_n = p_n - 1;

		if (a_n > 0 && !isnan(mua) && !isinf(mua))
			sda = sqrt(a / a_n);
		else
			sda = HUGE_VAL;
		if (p_n > 0 && !isnan(mup) && !isinf(mup))
			sdp = sqrt(p / p_n);
		else
			sdp = HUGE_VAL;
		params[0] = mua;
		params[1] = sda;
		params[2] = mup;
		params[3] = sdp;

		for (int j = 0; j != 4; ++j)
		{
			TH[i + j * N] = params[j];
		}

	} // end for
	free(params);
	return;
}
/*
 * Likelihood update
 */
double
calculate_likelihood(const double *X, int *GSsub, double *THsub, const int N,
		const int T, const int R)
{
	double dp, mu, sd, Lnew;
	int gp;
	Lnew = 0.0;
	for (int i = 0; i != N; ++i) // rows
	{
		// penalise the likelihood by the negative logged absolute difference in means
		double penalty = -1 * log(fabs(THsub[(int) (i + N * 0)] - THsub[(int) (i + N * 2)]));
		for (int j = 0; j != T * R; ++j) // columns
		{
			dp = X[i + j * N];
			gp = GSsub[i + j * N];
			if (gp == 1)
			{
				mu = THsub[(int) (i + N * 0)];
				sd = THsub[(int) (i + N * 1)];
			}
			else
			{
				mu = THsub[(int) (i + N * 2)];
				sd = THsub[(int) (i + N * 3)];
			}
			double dn = dnorm(dp, mu, sd);
			if (!isnan(dn) && !isinf(dn))
			{
				if (!isnan(penalty) && !isinf(penalty))
				{
					dn -= penalty;
				}
				Lnew += dn; // dnorm is logged data
			}
		}
	}
	return Lnew;
}

/*
 * normal density function
 */
double
dnorm(double dp, double mu, double sd)
{
	//double pi = 4 * atan(1);
	double factor = 1 / (sqrt(2 * M_PI) * sd);
	double exponent = -0.5 * (((dp - mu) * (dp - mu)) / (sd * sd));
	double density = factor * exp(exponent);
	density = log(density); // log to base 10
	return (density);
}
/*
 * print row sums and sum of a matrix
 * /
void
print_matrix_stats(double *x, int rows, int cols)
{
	double *rowsums = calloc(rows * cols, sizeof(double));
	double sum = 0.0;
	printf("Matrix stats (rowsums and sum): ");
	for (int i = 0; i != rows; ++i)
	{
		double rowsum = 0.0;
		for (int j = 0; j != cols; ++j)
		{
			rowsum += pow(2,x[i + j * rows]);
			sum += pow(2,x[i + j * rows]);
		}
		rowsums[i] = rowsum;
		printf("\n %f ", rowsum);
	}
	printf("\nSum: %f", sum);
	free(rowsums);
	return;
}
*/
/*
 * print a data matrix, passed as pointer
 * /
void
print_matrix(double *x, int rows, int cols)
{
	for (int i = 0; i != rows; ++i)
	{
		for (int j = 0; j != cols; ++j)
		{
			printf("%f ", x[i + j * rows]);
		}
		printf("\n");
	}
	printf("\n");
	return;
}

void
print_intmatrix(int *x, int rows, int cols)
{
	for (int i = 0; i != rows; ++i)
	{
		for (int j = 0; j != cols; ++j)
		{
			printf("%d ", x[i + j * rows]);
		}
		printf("\n");
	}
	printf("\n");
	return;
}
*/
void
normalise_rows(double *x, int rows, int cols, double *rowsums)
{
	for (int i = 0; i != rows; ++i)
	{
		for (int j = 0; j != cols; ++j)
		{
			x[i + j * rows] -= log2(rowsums[i]);
		}
	}
	return;
}

double
get_rowsums(double *x, int rows, int cols, double *rsums)
{
	double sum = 0.0;
	double rowsum;
	for (int i = 0; i != rows; ++i)
	{
		rowsum = 0.0;
		for (int j = 0; j != cols; ++j)
		{
			rowsum += pow(2,x[i + j * rows]);
			sum += pow(2,x[i + j * rows]);
		}
		rsums[i] = rowsum;
	}
	return sum;
}
void
updateA(double *A, int *maxemissionind, const int M, const int T)
{
	int Anum[M*M];
	for(int i=0; i!=M*M; ++i) { Anum[i] = 0; }
	int Adenom[M];
	for(int i=0; i!=M; ++i) { Adenom[i] = 0; }
	int pseudocount = 1;
	int pseudocountsum = M;
	int t, src, dst;

	// go along maxemissionind and count the transitions
	for (t = 0; t != (T - 1); t++)
	{
		src = maxemissionind[t]-1;
		dst = maxemissionind[t + 1]-1;
		int ind = dst*M + src;
		Anum[ind] = Anum[ind] + 1;
		Adenom[src] = Adenom[src] + 1;
	}
	// reset the transition probability to numocc / numtotaltransitionsfromsrc
	for (int i=0; i!=M; ++i) {
		for(int j=0; j!=M; ++j) {
			// found the transition
			if(Anum[j*M+i]!=0)
			{
				A[j*M+i] = log2(Anum[j*M+i] + pseudocount) - log2(Adenom[i]+pseudocountsum);
			}
		}
	}

	// record the rowSums for normalisation: transition matrix
	// is always normalised such that the rows sum to 1
	double *Asums = calloc(M * M, sizeof(double));
	double *rowsums = calloc(M, sizeof(double));
	double Asum = get_rowsums(A, M, M, rowsums);

	// now normalise
	normalise_rows(A, M, M, rowsums);

	free(rowsums);
	free(Asums);
}
/*
 * given a data vector and a state vector and the parameters theta,
 * get the emission probability p(x|gamma,th) = prod_i prod_r p(x_tir|gamma,th)
 */
double
getE(const double *X, int t, int *Gsub, int m, double *THsub, const int N,
		const int R)
{
	int startind = t * N * R; // start position in X
	int sel, gp;
	double dp, mu, sd, emission = 0;
	double dn1, dn2, mu1, mu2, sd1, sd2;

	// for each replicate
	for (int r = 0; r != R; ++r)
	{
		// for each protein
		for (int i = 0; i != N; ++i)
		{
			sel = (int) (startind + r * N + i);
			dp = X[sel]; // the measurement
			gp = Gsub[(int) m * N + i]; // state
			mu1 = THsub[(int) (i + N * 0)]; // mean active
			sd1 = THsub[(int) (i + N * 1)]; // standard deviation active
			mu2 = THsub[(int) (i + N * 2)]; // mean passive
			sd2 = THsub[(int) (i + N * 3)]; // standard deviation passive
			if (gp == 1)
			{
				dn1 = dnorm(dp, mu1, sd1);
				dn2 = dnorm(dp, mu2, sd2);
			}
			else
			{
				dn1 = dnorm(dp, mu2, sd2);
				dn2 = dnorm(dp, mu1, sd1);
			}
			// the emission probability log ratio
			double ret = dn1;
			if(isnan(ret))
				ret = 0.0;
			emission += ret;
		}
	}
	return (emission);
}

/*
 * get the basis 2 log
 */
double
log2(double x)
{
	double ret = log(x) / log(2);
	return (ret);
}

