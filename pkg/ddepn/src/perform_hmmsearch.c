//============================================================================
// Name        : peform_hmmsearch.c
// Author      : Christian Bender
// Version     :
// Copyright   : none
// Description : serial c-version of hmmsearch, adapted from perform.hmmsearch.R
//
//============================================================================

#include <stdio.h>
#include <stdlib.h>  // for random numbers
#include <time.h>    // for seed for random numbers
#include <unistd.h>     // need for getpid() function
#include <math.h>    // for pow function
#include "utility_functions.h"
#include "perform_hmmsearch.h"

/**
 * Performs HMM for network P and experiment in data matrix X
*/
void
perform_hmmsearch(int *Px, int *Nx, int *Tx, int *Rx, double *Xx, int *GSx,
    int *Gx, int *Glenx, double *THx, int *tpsx, int *stimidsx, int *stimgrpsx,
    int *numexperimentsx, double *Likx, int *hmmiterations, int *Msx)
{
  // create all objects that are needed. the const objects are not changed during
  // the run. the non - const objects are passed as reference and will be overridden
  // P is an array of length n * n, the adjacency matrix
  const int N = *Nx;
  //const int T = *Tx;
  const int hmmit = *hmmiterations;
  const double *X = Xx;
  double *TH = THx;
  int ix = 0, i, j;
  double Lik = 0.0; //*Likx;
  *Likx = 0.0;
  // set the seed using the current time and the process ID
  // this ensures that every instance uses different seeds.
  //srand(time(NULL) + getpid());

  // hmm search start
  Lik = hmmsearch_singleest(Px, N, Tx, Rx, Xx, GSx, Gx, Glenx, THx, tpsx, stimidsx, stimgrpsx,
      numexperimentsx, hmmit, Msx);
  *Likx = Lik;
  //printf("Likx: %f\n",*Likx);
  //printf("\n\nSUCCESS\n");
  return;
}

/**
 * For one network phi and one experiment extracted from the total data matrix X,
 * for one experiments states Gx use the viterbi algorithm
*/
double hmmsearch_singleest(int *phi, const int N, const int *Tx, const int *R,
    const double *X, int *GS,
    int *G, int *Glen, double *TH,
    const int *tps,
    const int *stimids, const int *stimgrps,
    const int *numexperimentsx, const int hmmit, int *Ms)
{
  // fixed for the moment, number of iterations in the em-algorithm
  int lenofexp, allR=0, allT=0, T;
  int start = 0, end; // start and end positions in the data and gamma matrix
  int numstims, idstart=0, Gstart=0, M;

  // extract each experiment from X according to R and stimgrps and run the hmm
  // indices of columns of X to be selected:
  // rind is equivalent to the experiment index
  // this defines the stimuli to take from stimgrps
  for(int rind=0; rind!=*numexperimentsx; rind++) {
    T = Tx[rind];
    //printf("replicates: %d, T: %d\n", R[rind], T);
    allR += R[rind]; // count the total number of replicates
    allT += T; // count the total number of timepoints

    M = Ms[rind];
    // extract the stimuli for the rind'th experiment
    int lstids = stimgrps[rind];
    int *stids = (int*)malloc(lstids*sizeof (int));

    for(int st=0; st!=lstids; st++) {
      stids[st] = stimids[idstart+st];
    }
    //printf("lstids: %d  ", lstids);

    // remember the first stimulus from the next experiment
    idstart+=lstids;

    // how long is the current experiment, i.e. how many columns in the Gamma and X matrices
    lenofexp = R[rind] * T;

    // get the sub-data-matrix
    end=start+lenofexp;

    // run the viterbi algorithm
    // Gstart contains the position in the G matrix at which new states have to be inserted
    Gstart = runviterbi(phi, N, T, R[rind], X, GS, start, end, G, Glen,
        stids, lstids, Gstart, hmmit, M);
    start+=lenofexp;
    free(stids);
  }

  // number of replicates are the same in each experiment, since
  // pad-columns containing NAs were added
  allR = allR / *numexperimentsx;

  // now the GS and G matrices are ready, get the final likelihood and TH values
  // update the theta matrix
  estimate_theta(X, GS, TH, N, allT, allR);

  // calculate the new likelihood
  double Lik = calculate_likelihood(X, GS, TH, N, T, allR); //Liktmp$L

  return(Lik);
}

/*
 * run the viterbi algorithm for one experiment and one network
 */
int
runviterbi(int *phi, const int N, const int T, int R, const double *X, int *GS,
    int start, int end, int *G, int *Glen, int *stids,
    int lstids, int Gstart, const int hmmit, int M)
{
  // perform state propagation to obtain the gamma matrix
  int lenofmat = N * T * R, k, runningindex=0; //, pseudocountsum;   M=0,

  // create the system state matrix
  //int *Gsub = malloc(lenofmat * sizeof(int));
  int *Gsub = malloc(N*M * sizeof(int));
  //printf("Gstart: %d, G orig: \n", Gstart);
  //print_intmatrix(G, N, *Glen/N);
  for(int i=0; i!=N*M; ++i) {
      Gsub[i] = G[Gstart + i];
  }
  //printf("G subset: \n");
  //print_intmatrix(Gsub, N, M);

  int *GSsub = calloc(lenofmat, sizeof(int));

  // here the Gsub matrix is filled with the system states
  //M = propagate_effect(phi, stids, lstids, Gsub, N);
  //M = *Glen/N;
  //printf("BACK from propeffect: Found M=%d system states. Gstart: %d\n", M, Gstart);

  // free the unused memory in Gsub
  //Gsub = realloc(Gsub, (M*N)*sizeof(int));

  // get a theta matrix
  //double *THsub = malloc(4*N * sizeof(double));
  double *THsub = calloc(4*N, sizeof(double));
  estimate_theta(X, GSsub, THsub, N, T, R);

  // run the hmm
  hmm(X, Gsub, GSsub, THsub, N, T, R, M, hmmit);

  // start, end denote the start and end *columns*, get the indices:
  int GSstart = start * N;
  int GSend = end * N;
  // add to GS and G
  for(k=GSstart; k!=GSend; k++) {
      GS[k] = GSsub[runningindex];
      runningindex++;
  }
  int temp = (Gstart+M*N);
  /*
  // check if G is too small - shouldn't happen, though
  if(temp>*Glen) {
      //G = realloc(G, (Gstart+M*N) * sizeof(int)); //allocate more space if G wasn't big enough before
      printf("Warning: Too little memory allocated for G! temp: %d, Gstart %d, M: %d, N: %d\n", temp, Gstart, M, N);
  }
  */
  // printf("start: %d, end: %d, GSstart: %d, GSend: %d, Gstart: %d, temp: %d\n",start, end,GSstart, GSend, Gstart, temp);
  free(GSsub);
  free(Gsub);
  free(THsub);
  return(temp);
}

void
hmm(const double *X, int *Gsub, int *GSsub, double *THsub, const int N,
    const int T, const int R, int M, const int hmmit)
{
  // sort of a hack for getting the - infinity value
  double temp = 1.0;
  double infinity = -1 * (temp / (temp - 1.0));

  // initialise the state matrix
  initialise_GS(GSsub, Gsub, N, T, R, M);

  // initialise the transition probability matrix
  double *A = malloc(M * M * sizeof(double));
  initialise_A(A, M);

  // the viterbi matrix as vector
  double *viterbi = malloc(M * T * R * sizeof(double));

  // will hold the state indices with maximum likelihood
  int *maxemissionind = malloc(T * sizeof(int));

  double Lik = infinity, Lold;
  int diffold = -100, equally = 0, it, restarts = 0, tmp;
  double diff, diffsold[hmmit];
  for(int ddd=0; ddd!=hmmit; ++ddd) { diffsold[ddd] = ddd; } //HUGE_VAL; }

  // get initial likelihood
  //double Lnew = calculate_likelihood(X, GSsub, THsub, N, T, R);

  // main baum-welch loop
  //while (it <= hmmit)
  for(it=0; it!=hmmit; ++it)
    {
      //old likelihood
      Lold = Lik;

      // update the theta matrix
      estimate_theta(X, GSsub, THsub, N, T, R);

      // calculate the new likelihood
      Lik = calculate_likelihood(X, GSsub, THsub, N, T, R); //Liktmp$L

      if (isinf(Lold) != 1)
        {
          if (fabs((fabs(Lik) - fabs(Lold))) <= 1.0)
            break;
          diff = fabs((fabs(Lik) - fabs(Lold)));
          diffsold[it] = diff;
          // check for repeating patterns
          tmp = find_diff(diff, diffsold, it);
            if (tmp > 10)
              {
                //printf("REPEAT FOUND!!\n");
                // set equally to 2; the value of 2 here is used only because
                // of the next if statement, where the condition equally==3 is
                // supposed to be fulfilled after a repeating pattern was found
                equally = 2;
                diff = diffold;
              }
          if (diff == diffold)
            {
              equally++;
              //restart if switching behaviour occurs,
              //don't know where this comes from
              if (equally == 3)
                {
                  restarts++;
                  //only up to a number of restarts
                  if (restarts < 3)
                    {
                      //printf("RESTART!!\n");
                      // reset to some non-repeating dummy values
                      for(int ddd=0; ddd!=it; ++ddd) { diffsold[ddd] = ddd; }
                      // again initialise the state matrix
                      //initialise_A(A, M);
                      initialise_GS(GSsub, Gsub, N, T, R, M);
                      Lik = infinity;
                      diffold = -100;
                      equally = 0;
                      //restarts = 0;
                      continue;
                    }
                  else
                    {
                      break;
                    }
                }
            }
          diffold = diff;
        } // end if(isinf(Lold)!=1)

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
                  viterbi[ind] = -log2(M) + getE(X, t, Gsub, m, THsub, N, R);
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
                  viterbi[ind] = maxviterbi + getE(X, t, Gsub, m, THsub, N, R);
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
                  GSsub[(int) (start + i)] = Gsub[(int) (state * N + i)];
                }
            }
        }
      // M step: update transition probabilities
      updateA(A, maxemissionind, M, T);
    } // end while()
  free(maxemissionind);
  free(viterbi);
  free(A);
}

/*
 * initialisation of the state transition matrix
 * draw uniform probabilities
 */
void
initialise_A(double *A, int M)
{
  // set the seed using the current time and the process ID
  // this ensures that every instance uses different seeds.
 /* struct tm *current;
  time_t now;
  time(&now);
  current = localtime(&now);
  printf("the time is %i:%i:%i\n", current->tm_hour, current->tm_min, current->tm_sec);*/
  //srand(time(NULL) + getpid());
  double rnum, rnumsum = 0.0;
  for (int i = 0; i != M; ++i)
    {
      for (int j = 0; j != M; ++j)
        {
          rnum = ((double) rand() / (double) RAND_MAX);
          rnumsum += rnum;
          A[(i * M + j)] = rnum;
        }
    }
  //double *Asums = calloc(M * M, sizeof(double));
  double *rowsums = calloc(M, sizeof(double));
  double Asum = get_rowsums(A, M, M, rowsums);
  normalise_rows(A, M, M, rowsums);
  free(rowsums);
}

/*
 * look for value diff in vector diffsold and count
 * occurrences. maxit is the length of diffsold
 */
int
find_diff(double diff, double *diffsold, int maxit)
{
  int occ = 0;
  for (int i = 0; i != maxit; i++)
    {
      if (diffsold[i] == diff)
        occ++;
    }
  return (occ);
}

int
propagate_effect(int *phi, int *stids, int lstids, int *Gsub, const int N)
{
  int i, j, lenstims, binsum=0, index = N, goon = 1, it = 0;
  // allocate an array of integers of size 255
  // if this is too little, double the size of the array and reallocate and so on
  int arrsizeexp = 8;
  int arrsize = pow(2, arrsizeexp);
  int *binsums = (int*) calloc(arrsize, sizeof(int));

  //init_Gsub all stims are 1, rest 0
  for (i = 0; i != N; i++)
    {
      // check if protein is a stimulus
      if (is_stimulus(i, stids, lstids) == 1)
        {
          Gsub[i] = 1;
        }
      else
        {
          Gsub[i] = 0;
        }
      binsum += Gsub[i] * pow(2, i);
    }

  binsums[it++] = binsum;

  // create a new column for Gsub
  while (goon == 1)
    {
      int *newcol = (int*) malloc(N * sizeof(int));
      binsum = next_col(newcol, Gsub, N, stids, lstids, phi, index);
      //printf("new binsum: %d\n", binsum);

      // if the column was already found in the Gsub matrix, we're done
      if (check_binsum(binsums, binsum, it) == 1)
        {
          goon = 0;
        }
      else
        {
          // else add the columns id to the binsums vector and continue
          if (it > arrsize)
            {
              // extend the binsums array if not big enough
              arrsize = pow(2, ++arrsizeexp);
              binsums = (int *) realloc(binsums, arrsize * sizeof(int));
            }

          // and add the column id to binsums and the column to Gsub
          binsums[it] = binsum;
          for (i = 0; i != N; i++)
            {
              Gsub[(index + i)] = newcol[i];
            }
          it++;
        }
      index += N;
      free(newcol);

    }
  free(binsums);
  return (it);
}

/*
 * given parts of a matrix Gsub, create the next column in Gsub
 * N: number of nodes
 * stimids: indices of the stimuli nodes
 * phi: the adjacency matrix
 * index: the actual position in Gsub
 */
int
next_col(int *newcol, int *Gsub, int N, int *stids, int lstids, int *phi,
    int index)
{
  int i, j, tbefore, tnow, phiindex, binsum = 0;
  // for each node, get the parent nodes,
  // get the states of the parent nodes in the column
  // before and combine to the actual state
  for (i = 0; i != N; i++)
    {
      tnow = index + i;
      if (is_stimulus(i, stids, lstids) == 1)
        {
          newcol[i] = 1;
          binsum += (int)(newcol[i] * pow(2, i));
          continue;
        }
      else
        {
          newcol[i] = 0;
        }
      for (j = 0; j != N; j++)
        {
          if (i == j)
            continue;
          tbefore = index - N + j;
          // found an activation parent
          phiindex = i * N + j;
          if (phi[phiindex] == 1)
            {
              // parent was active at time before
              if (Gsub[tbefore] == 1)
                {
                  //Gsub[tnow] = 1;
                  newcol[i] = 1;
                }
            }
          // found an inhibition parent
          if (phi[phiindex] == 2)
            {
              // parent was active at time before => signal shut off
              // go to next node, since signal is off
              if (Gsub[tbefore] == 1)
                {
                  newcol[i] = 0;
                  break;
                }
            }
        }// second for
      // now newcol[i] is set, add to the binary sum for the column check
      binsum += (int)(newcol[i] * pow(2, i));
    }
  return (binsum);
}

/*
 * find binsum in array binsums, up to element it
 */
int
check_binsum(int *binsums, int binsum, int it)
{
  if (it == 0)
    {
      return (0);
    }
  for (int i = 0; i != it; ++i)
    {
      if (binsum == binsums[i])
        return (1);
    }
  return (0);
}

/*
 * check if a protein x is present in the sitmuli vectors
 */
int
is_stimulus(int x, int *stids, int lstids)
{
  for (int it = 0; it != lstids; it++)
    {
      if (stids[it] == x)
        {
          return (1);
        }
    }
  return (0);
}

/*
 * initialise the system state matrix that should be optimised
 * draw randomly M columns from Gsub in order of Gsub and
 * copy the element of these columns to GSsub
 */
void
initialise_GS(int *GSsub, int *Gsub, int N, int T, int R, int M)
{
  // set the seed using the current time and the process ID
  // this ensures that every instance uses different seeds.
  srand(time(NULL) + getpid());
  int rnum, i, gcolit, gcol, it = 0;
  int *rnums = (int*) malloc(T * R * sizeof(int));

  for (i = 0; i != T; ++i)
    {
      //sample columns from 1...M and add R times after each other
      rnum = ((int) rand() % M);
      for (it = i * R; it != (i + 1) * R; it++)
        {
          rnums[it] = rnum;
        }
      //printf("randnum: %d\n", rnum);
    }
  // sort the column sampling
  qsort(rnums, T * R, sizeof(int), compare_ints);

  // fill the GSsub matrix with the columns from Gsub
  gcolit = 0;
  for (i = 0; i != T * R * N; ++i)
    {
      if (i % N == 0 && i > 0)
        {
          gcolit++;
        }
      gcol = rnums[gcolit]; // select the sampled column from GSub
      GSsub[i] = Gsub[gcol * N + i % N];
    }
  free(rnums);
}

int
compare_doubles(const void *a, const void *b)
{
  const double *da = (const double *) a;
  const double *db = (const double *) b;

  return (*da > *db) - (*da < *db);
}

int
compare_ints(const void *a, const void *b)
{
  const int *ia = (const int *) a;
  const int *ib = (const int *) b;

  return (*ia > *ib) - (*ia < *ib);
}
