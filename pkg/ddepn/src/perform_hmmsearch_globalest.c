//============================================================================
// Name        : perform_hmmsearch_globalest.c
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
#include "utility_functions.h"
#include "perform_hmmsearch_globalest.h"
/**
 * Performs HMM for network P and experiment in data matrix X
 * stimgrpsx: vector containing the lengths of each stimulus
 * numexperiments: number of experiments, i.e. length of stimuli
 *
 * stimidsx: unlist(stimuli)
 * stimgrpsx:  sapply(stimuli, length)
 * Xtt <- table(sub("_[0-9]+$","",colnames(X)))/T
 * R <- Xtt[match(names(Xtt), ord)]
 * numexperimentsx: length(R)
 */
void
perform_hmmsearch_globalest(int *Px, int *Nx, int *Tx, int *Rx, double *Xx, int *GSx,
		int *Gx, int *Glenx, double *THx, int *tpsx, int *stimidsx, int *stimgrpsx,
		int *numexperimentsx, double *Likx, int *hmmiterations, int *Msx)
{
	// create all objects that are needed. the const objects are not changed during
	// the run. the non - const objects are passed as reference and will be overridden
	// P is an array of length N * N, the adjacency matrix
	const int N = *Nx; // num of nodes
	const int T = *Tx; // num of timepoints
	const int hmmit = *hmmiterations;
	const double *X = Xx; // data matrix
	const int Glen = *Glenx;
	const int numexperiments = *numexperimentsx;
	double *TH = THx; // parameter matrix
	int ix = 0, i, j;
	double Lik = 0.0;
	*Likx = 0.0;

	// hmm search start
	Lik = hmmsearch(Px, N, T, Rx, Xx, GSx, Gx, Glen, THx, tpsx, stimidsx, stimgrpsx,
			numexperiments, hmmit, Msx);
	*Likx = Lik;
	return;
}

/**
 * For one network phi and one experiment extracted from the total data matrix X,
 * for one experiments states Gx use the viterbi algorithm
 * N: number of nodes
 * T: number of time points
 * R: number of replicates
 * X: data matrix (N x TxR)
 * GS: optim state matrix (N x TxR), initialised in R and given as argument
 * G: state matrix (N x sum_s(M_s)); for each stimulus s, there are M_s states
 * Glen:
 * TH: theta parameter matrix (N x 4)
 * tps: timepoint vector (T)
 * stimgrps: vector containing the stimulus indices
 * numexperimentsx: number of experiments (equals length of R and Ms)
 * hmmit: number of hmmiterations
 * Ms: number of system states for each experiment
 */
double hmmsearch(int *phi, const int N, const int T, const int *R,
		const double *X, int *GS,
		int *G, int Glen, double *TH,
		const int *tps,
		const int *stimids, const int *stimgrps,
		const int numexperiments, const int hmmit, int *Ms)
{
	// fixed for the moment, number of iterations in the em-algorithm
	int ncol_GS=0;
	int numstims, idstart=0, Gstart=0;

	for(int i=0; i!=numexperiments; ++i) {
		ncol_GS += T*R[i];
	}
	//  printf("~~~~~ \n ncol_GS: %d\n",ncol_GS);
	// init A, TH and L
	// sort of a hack for getting the - infinity value
	double temp = 1.0;
	double infinity = -1 * (temp / (temp - 1.0));
	int M = Glen/N; // total number of system states

	// allocate the transition probability matrix, for all experiments
	int A_sz = 0;
	for(int i=0; i!=numexperiments; ++i) {
		A_sz += pow(Ms[i],2);
	}

	// make a vector of doubles, holding the transition probabilities
	double *A = malloc(A_sz * sizeof(double));
	init_A(A, Ms, numexperiments); // init transition probabilities (sparse MxM matrix)

	/* main loop: for each iteration in hmmiterations:
	 * perform viterbi for each experiment separately
	 * get the parameters TH and updates for A and GS for all experiments combined
	 */
	int Mexp;
	int startA=0, startX=0, startG=0;
	double *Aexp = NULL;
	double *Xexp = NULL;
	int *Gexp = NULL;
	int *GSexp = NULL;
	int allR;
	double Lik = -1*infinity;
	double Likold = -1*infinity;
	int nLikEqual = 0;
	// keep track of the last 5 likelihood differences
	// if switching behaviour occurs, it can be seen here
	// take differences
	int K = 6;
	double diff, difftmp;
	//double *diffvec = calloc(K-1, sizeof(double));
	double *diffvec = malloc(K * sizeof(double));
	int *toswitch = calloc(N, sizeof(int));
	int nsw=0, maxsw=10;
	// new state matrix object
	for(int it=0; it!=hmmit; ++it) {
		allR=0;
		startA=0;
		startX=0;
		startG=0;
		// find switchable rows in TH
		nsw = find_switchable(TH, N, toswitch);
		if(nsw>0) {
			maxsw = min(nsw, maxsw);
		}
		/* here the experiment loop
		 * extract each experiment from X according to R and stimgrps and run the hmm
		 * indices of columns of X to be selected:
		 * expind is equivalent to the experiment index
		 * this defines the stimuli to take from stimgrps
		 */
		for(int expind=0; expind!=numexperiments; expind++) {
			//print_intmatrix(toswitch, 1, N);

			// if inconsistencies occur in the gamma matrix,
			// try switching the theta parameters upto
			// maxsw times to reduce the number of inconsitencies
			//maxsw = min(nsw, N);
			allR += R[expind];
			Mexp = Ms[expind];
			for(int swind=0; swind!=maxsw; ++swind) {

				// extract the sub-transition matrix:
				Aexp = realloc(Aexp, Mexp*Mexp * sizeof(double));
				extract_transitionmatrix(A, Aexp, Mexp, startA);

				// extract the sub-data matrix
				Xexp = realloc(Xexp, N*T*R[expind]*sizeof(double));
				extract_datamatrix(X, Xexp, N, T, R[expind], startX);

				// extract the sub-state matrix
				Gexp = realloc(Gexp, N*Mexp*sizeof(int));
				extract_statematrix(G, Gexp, N, Mexp, startG);

				// extract the sub-optimstate matrix
				GSexp = realloc(GSexp, N*T*R[expind]*sizeof(int));
				extract_statematrix(GS, GSexp, N, T*R[expind], startX);

				// run the viterbi
				viterbi(T, Mexp, Xexp, Gexp, TH, N, R[expind], Aexp, GSexp);

				// consitency check
				int inc = is_consistent(phi, GSexp, Gexp, N, T, R[expind]);

				if(inc>0 && nsw>0) {
					//printf("Inconsitensies in state series. Repeat HMM with modified thetaprime.\n");
					// select a row to switch randomly
					int rnum = rand() % nsw; // a random number between 1 and nsw
					int hit = 0, ii=0;
					for(ii=0; ii!=N; ++ii) {
						if(toswitch[ii]!=0) {
							if(hit==rnum) {
								break;
							}
							hit++;
						}
					}
					switch_theta_row(TH, ii, N);
					// remove node as possible switch node
					toswitch[ii] = 0;
					nsw--;
				} else {
					// update the GSnew matrix
					update_statematrix(GS, GSexp, startX, N, T, R[expind]);

					// update the new transition prob matrix
					update_transitionmatrix(A, Aexp, Mexp, startA);
					break;
				}
			} // switch loop end
			// increment the experiment start indices
			startA += pow(Mexp, 2);
			startX += N*T*R[expind];
			startG += N*Mexp;
		} // experiment loop end

		// M-step
		// update the theta matrix
		estimate_theta(X, GS, TH, N, T, allR);

		// calculate the new likelihood
		Lik = calculate_likelihood(X, GS, TH, N, T, allR); //Liktmp$L
		diff = fabs((fabs(Lik) - fabs(Likold)));

		// count number of equal differences in the last 10 likelihoods
		// if all 5 differences are equal, then stop
		int stopit=0, count = 0; // if all diffs are equal, stopit remains 1 and the loop is aborted
		// check if diffvec is filled
		if(it>=K) {
			// check elements in diffvec
			for(int k=0; k!=(K-1); ++k) {
				if(k>0) {
					if(fabs(diffvec[k]-difftmp)<=0.001) {
						count++;
						//stopit = 0;
					}
				}
				// remember the original value at current position
				difftmp = diffvec[k];
				// shift next value left one position
				if(k<K) {
					diffvec[k] = diffvec[k+1];
				} else {
					// update the new difference at last position
					diffvec[k] = diff;
				}
			}
			if(count==(K-1)) {
				stopit = 1;
				// make sure that the higher likelihood is taken
				// if switching occurs
				if(Likold>Lik) {
					stopit = 0;
				}
			}
		} else {
			// if not filled then add elements
			diffvec[it] = diff;
		}
		// another termination criterion: count if liklihood terms do not change
		if(Lik==Likold) {
			nLikEqual++;
		} else {
			nLikEqual = 0;
			Likold = Lik;
		}

		// abort baum-welch, if Lik does not change anymore
		if(nLikEqual>=10 || stopit==1) {
			break;
		}
	}
	free(toswitch);
	free(diffvec);
	free(GSexp);
	free(Gexp);
	free(Xexp);
	free(Aexp);
	free(A);
	return Lik;
}
/*
 * switch the active and passive states of a row rnum in TH
 */
void switch_theta_row(double *TH, int rnum, int N)
{
	double mua=TH[0*N + rnum];
	double sda=TH[1*N + rnum];
	double mup=TH[2*N + rnum];
	double sdp=TH[3*N + rnum];

	TH[0*N + rnum] = mup;
	TH[1*N + rnum] = sdp;
	TH[2*N + rnum] = mua;
	TH[3*N + rnum] = sda;
}
/*
 * inspect rows of TH
 * if no NA value is in a row or
 * if all values differ from 0.0
 * then switch of the active and passive
 * parameters is permitted
 * returns a vector of length N containing
 * 1 for each switchable row and 0 otherwise
 */
int find_switchable(double *TH, int N, int *toswitch)
{
	int nacnt, zerocnt, swcnt=0;
	for(int i=0; i!=N; ++i) {
		nacnt=0;
		zerocnt=0;
		for(int j=0; j!=4; ++j) {
			if(isfinite(TH[j*N + i])==0) {
				nacnt++;
			}
			if(TH[j*N + i]==0.0) {
				zerocnt++;
			}
		}
		if(nacnt>0 || zerocnt==4) {
			toswitch[i] = 0;
		} else {
			toswitch[i] = 1;
			swcnt++;
		}
	}
	return swcnt;
}
/* traverse all nodes
 * find inconsitency markers for each child depending on its parents
 *   inconsistency is found whenever:
 *   edge pa --> ch
 *        gammax[pa,t_i] == 1 # if the parent is active, then state switch from active to passive should not occur at the child
 *                      gammax[ch, t_i] == 1 & gammax[ch, t_(i+1)] == 0 not possible
 *        gammax[pa,t_i] == 0 # if the parent is passive, then state switch from passive to active should not occur at the child
 *                      gammax[ch, t_i] == 0 & gammax[ch, t_(i+1)] == 1 not possible
 *   edge pa --| ch
 *        gammax[pa,t_i] == 1 # if parent is active, then state switch from passive to active should not occur
 *                      gammax[ch, t_i] == 0 & gammax[ch, t_(i+1)] == 1 not possible
 *        gammax[pa,t_i] == 0
 *                      gammax[ch, t_i] == 1 & gammax[ch, t_(i+1)] == 0 not possible
 *
 */
int is_consistent(int *phi, int *GS, int *G, int N, int T, int R)
{
	int inc=0,
			cnt=0,
			pcnt=0;

	// initialise an array of parents, can be at most N elements
	int *parents = malloc(N*sizeof(int));

	// allocate consistency vector
	int *consvec = malloc((T*R-1)*sizeof(int));

	// for all nodes
	for(int i=0; i!=N; ++i) {

		// reset parents vector to 0
		memset(parents, 0, N*sizeof(int));

		// reset the consistency vector to all consistent
		for(int ci=0; ci!=(T*R-1); ++ci) {
			consvec[ci] = 1;
		}

		// check all parents
		pcnt = get_parents(phi, N, parents, i);

		if(pcnt>0) {
			// for each parent
			for(int j=0; j!=N; ++j) {
				if(parents[j]!=0) {
					checkC(phi, GS, i, j, N, T, R, consvec);
					// reset -14 value to 1
					for(int ci=0; ci!=(T*R-1); ci++) {
						if(consvec[ci]==-14) {
							consvec[ci] = 1;
						}
					}
				}
			} // parents end
			// increment the inconsistency count
			for(int ci=0; ci!=(T*R-1); ++ci) {
				if(consvec[ci]==0) {
					inc++;
				}
			}
		}

	}
	free(consvec);
	free(parents);
	return inc;
}
/*
 * check consistency of a child and on of its parents in a
 * state matrix GS.
 * returns a vector of integers with length equalling the number of
 * columns in GS-1. Values different from 0 indicate an inconsistency.
 */
void checkC(int *phi, int *GS, int ch, int pa, int N, int T, int R, int *consvec)
{
	// get the edge type
	int ed = phi[ch*N + pa], paind, chind, chind_next;
	char *eds = "ai";
	// each column
	for(int j=1; j!=(T*R); ++j) {
		paind = ((j-1))*N + pa;
		chind = ((j-1))*N + ch;
		chind_next = (j*N + ch);
		switch(eds[(ed-1)]){
		case 'a':
		// if parent is 0
			if(GS[paind]==0) {
				if(GS[chind]==0 && GS[chind_next]==0) { consvec[j-1] = min(consvec[j-1], 1); } // (0)
				if(GS[chind]==0 && GS[chind_next]==1) { consvec[j-1] = min(consvec[j-1], 1);} // (1) could be inconsistent, but is undeterminable, therefore set to 1
				if(GS[chind]==1 && GS[chind_next]==0) { consvec[j-1] = min(consvec[j-1], 1); } // (2)
				if(GS[chind]==1 && GS[chind_next]==1) { consvec[j-1] = min(consvec[j-1], 1); } // (3)
			} else {
				// or not pa == 0
				if(GS[chind]==0 && GS[chind_next]==0) { consvec[j-1] = min(consvec[j-1], 1); } // (4)
				if(GS[chind]==0 && GS[chind_next]==1) { consvec[j-1] = min(consvec[j-1], 1); } // (5)
				if(GS[chind]==1 && GS[chind_next]==0) { consvec[j-1] = min(consvec[j-1], 0); } // (6)
				if(GS[chind]==1 && GS[chind_next]==1) { consvec[j-1] = min(consvec[j-1], 1); } // (7)
			}
			break;
		case 'i':
			if(GS[((j-1)*N + pa)]==0) {
				if(GS[chind]==0 && GS[chind_next]==0) { consvec[j-1] = min(consvec[j-1], 1); } // (8)
				if(GS[chind]==0 && GS[chind_next]==1) { consvec[j-1] = min(consvec[j-1], 1); } // (9)
				if(GS[chind]==1 && GS[chind_next]==0) { consvec[j-1] = min(consvec[j-1], 1); } // (10)
				if(GS[chind]==1 && GS[chind_next]==1) { consvec[j-1] = min(consvec[j-1], 1); } // (11)
			} else {
				if(GS[chind]==0 && GS[chind_next]==0) { consvec[j-1] = min(consvec[j-1], 1); } // (12)
				if(GS[chind]==0 && GS[chind_next]==1) { consvec[j-1] = min(consvec[j-1], 0); } // (13)
				if(GS[chind]==1 && GS[chind_next]==0) { consvec[j-1] = min(consvec[j-1], -14); }  // (14) special case, can override case (6)
				if(GS[chind]==1 && GS[chind_next]==1) { consvec[j-1] = min(consvec[j-1], 0); } // (15)
			}
			break;
		}
	}
	return;
}

/*
 * get the parents of one leaf
 * store in vector parents
 */
int get_parents(int *phi, int N, int *parents, int leaf)
{
	int pcnt = 0;
	for(int i=0; i!=N; ++i) {
		// if any non-zero entry is found int the row, then
		// node is no leaf; jump back to row loop
		if(phi[leaf*N + i]!=0) {
			pcnt++;
			parents[i] = 1;
		}
	}
	return pcnt;
}
/*
 * get the leafs in a network phi
 * store in vector leafs
 * /
int get_leafs(int *phi, int N, int *leafs)
{
	int lcnt = N;
	for(int i=0; i!=N; ++i) {
		leafs[i] = 1;
		for(int j=0; j!=N; ++j) {
			// if any non-zero entry is found int the row, then
			// node is no leaf; jump back to row loop
			if(phi[j*N + i]!=0) {
				lcnt--;
				leafs[i] = 0;
				break;
			}
		}
	}
	return lcnt;
}*/

/*
 * update the optimised state matrix
 */
void update_statematrix(int *GS, int *GSexp, int start, int N, int T, int R)
{
	int ind = 0;
	int end = start + N*T*R;
	for(int i=start; i!=end; ++i) {
		GS[i] = GSexp[ind];
		ind++;
	}
}
/*
 * update the transition vector src using a sub-transition matrix
 * dst.
 * start from index up to index + nstates^2 in src and set each value
 * to the corresponding value in dst
 */
void update_transitionmatrix(double *src, double *dst, int nstates, int index)
{
	for(int i=0; i!=nstates; ++i) {
		for(int j=0; j!=nstates; ++j) {
			src[index] = dst[j*nstates + i];
			//dst[j*nstates + i] = src[index];
			index++;
		}
	}
}
/*
 * fill the viterbi matrix and identify the maximum likelihood path
 * through the system states
 */
void viterbi(int T, int M, double *X, int *G, double *TH, int N, int R, double *A, int *GS)
{
	// sort of a hack for getting the - infinity value
	double temp = 1.0;
	double infinity = -1 * (temp / (temp - 1.0));

	// the viterbi matrix as vector
	double *viterbi = malloc(M * T * R * sizeof(double));
	//for(int ind=0; ind!=M*T*R; ++ind) {
	//  viterbi[ind] = 0;
	//}

	// indices of most likely states
	int *maxemissionind = malloc(T * sizeof(int));

	double maxemission = infinity, viterbibefore, newem, maxviterbi;
	for (int t = 0; t != T; ++t)
	{ // all time points
		viterbibefore = maxemission;
		maxemission = infinity;
		for (int m = 0; m != M; ++m)
		{ // all states in gamma
			int ind = (int)(t*M+m);
			if (t == 0)
			{
				//initialise
				viterbi[ind] = -log2(M) + getE(X, t, G, m, TH, N, R);
			}
			else
			{
				int Aind, Vind;
				double tmp;
				maxviterbi = infinity;
				for(int m2=0; m2!=M; ++m2) {
					Aind = m*M+m2;
					Vind = (t-1)*M+m2;
					tmp = A[Aind] + viterbi[Vind];
					if(tmp>maxviterbi) {
						maxviterbi = tmp;
					}
				}
				viterbi[ind] = maxviterbi + getE(X, t, G, m, TH, N, R);
			}
			// remember the maximal value in this time point + the
			// state from which the data point was generated
			if (viterbi[ind] >= maxemission)
			{
				maxemission = viterbi[ind];
				maxemissionind[t] = m;
			}
		}
	}
	// fill the gammaprime matrix
	for (int t = 0; t != T; t++)
	{
		int state = maxemissionind[t];
		int startind = (int) t * N * R;

		for (int r = 0; r != R; r++)
		{
			int start = startind + r * N;
			for (int i = 0; i != N; i++)
			{
				GS[(int) (start + i)] = G[(int) (state * N + i)];
			}
		}
	}

	// M step: update transition probabilities for sub matrix
	updateA(A, maxemissionind, M, T);

	free(maxemissionind);
	free(viterbi);
	return;
}

void extract_datamatrix(const double *X, double *Xexp, int N, int T, int R, int start)
{
	int end = start + N*T*R;
	int ind = 0;
	for(int i=start; i!=end; ++i) {
		Xexp[ind] = X[i];
		ind++;
	}
}
void extract_statematrix(int *X, int *Xexp, int N, int cols, int start)
{
	int end = start + N*cols;
	int ind = 0;
	for(int i=start; i!=end; ++i) {
		Xexp[ind] = X[i];
		ind++;
	}
}

/*
 * extract a submatrix from the transition vector src
 * to destination matrix dst
 * starting from index up to index + nstates^2
 */
void extract_transitionmatrix(double *src, double *dst, int nstates, int index)
{
	for(int i=0; i!=nstates; ++i) {
		for(int j=0; j!=nstates; ++j) {
			dst[j*nstates + i] = src[index];
			index++;
		}
	}
}
/*
 * initialise the state matrix
 * use a vector representation, which saves a lot
 * of memory and time
 */
void init_A(double *A, int *Ms, int numexperiments)
{
	double rnum = 0.0;
	//double rnum;
	int startind = 0;
	for (int i = 0; i != numexperiments; ++i)
	{
		int Mi = Ms[i];
		int Misq = pow(Mi, 2);
		double *rowsums = calloc(Mi, sizeof(double)); // alloc space and init to 0
		for (int j = 0; j != Misq; ++j)
		{
			rnum = ((double) rand() / (double) RAND_MAX); // a random number
			A[(startind + j)] = rnum;
			rowsums[j % Mi] += A[(startind + j)];
		}
		// normalize transition matrix such that rowSums are 0
		for(int j = 0; j != Misq; ++j)
		{
			A[(startind + j)] = log2(A[(startind + j)]) - log2(rowsums[j % Mi]);
		}
		free(rowsums);
		// increment the index selecting the start of the
		// single experiment's transition matrix
		startind += Misq;
	}

}
/*
 * if Ms contains the limits of a number of bins,
 * find out, in which bin the number num falls.
 * numexperiments is the length of the vector Ms
 */
int find_bin(int *Ms, int num, int numexperiments)
{
	int bin = 0;
	int runM = 0;
	for(int i=0; i!=numexperiments; ++i)
	{
		runM += Ms[i];
		if(num<runM) {
			bin = i;
			break;
		}
	}
	return(bin);
}
/*
 * print the transition vector in matrix form
 */
void
print_transitions(double *A, int *Ms, int numexperiments, int Msum)
{
	int iind, jind;
	int Aind = 0;
	// build index matrix
	int *indmat = calloc(Msum*Msum, sizeof(int));
	// rows
	for(int j=0; j!=Msum; ++j)
	{
		jind = find_bin(Ms, j, numexperiments);
		// cols
		for(int i=0; i!=Msum; ++i)
		{
			iind = find_bin(Ms, i, numexperiments);
			if(jind==iind) {
				indmat[j*Msum + i] = Aind;
				Aind++;
			}
		}
	}
	// rows
	for(int i=0; i!=Msum; ++i)
	{
		iind = find_bin(Ms, i, numexperiments);
		// cols
		for(int j=0; j!=Msum; ++j)
		{
			jind = find_bin(Ms, j, numexperiments);
			if(jind!=iind) {
				printf("%f ", 0.0);
			} else {
				printf("%f ", A[indmat[j*Msum + i]]);
			}
		}
		printf("\n");
	}
	free(indmat);
}




