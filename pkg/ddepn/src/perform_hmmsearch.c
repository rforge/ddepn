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
#include "perform_hmmsearch.h"


/**
 * Performs HMM for network P and experiment in data matrix X
*/
void
perform_hmmsearch(int *Px, int *Nx, int *Tx, int *Rx, double *Xx, int *GSx,
    int *Gx, int *Glenx, double *THx, int *tpsx, int *stimidsx, int *stimgrpsx,
    int *numexperimentsx, double *Likx, int *hmmiterations)
{
  // create all objects that are needed. the const objects are not changed during
  // the run. the non - const objects are passed as reference and will be overridden
  // P is an array of length n * n, the adjacency matrix
  const int N = *Nx;
  const int T = *Tx;
  const int hmmit = *hmmiterations;
  const double *X = Xx;
  double *TH = THx;
  int ix = 0, i, j;
  double Lik = *Likx;

  // set the seed using the current time and the process ID
  // this ensures that every instance uses different seeds.
  //srand(time(NULL) + getpid());

  // hmm search start
  Lik = hmmsearch(Px, N, T, Rx, Xx, GSx, Gx, Glenx, THx, tpsx, stimidsx, stimgrpsx,
      numexperimentsx, hmmit);
  *Likx = Lik;
  //printf("\n\nSUCCESS\n");
  return;
}

/**
 * For one network phi and one experiment extracted from the total data matrix X,
 * for one experiments states Gx use the viterbi algorithm
*/
double hmmsearch(int *phi, const int N, const int T, const int *R,
    const double *X, int *GS,
    int *G, int *Glen, double *TH,
    const int *tps,
    const int *stimids, const int *stimgrps,
    const int *numexperimentsx, const int hmmit)
{
  // fixed for the moment, number of iterations in the em-algorithm
  int lenofexp, allR=0;
  int start = 0, end; // start and end positions in the data and gamma matrix
  int numstims, idstart=0, Gstart=0;

  // extract each experiment from X according to R and stimgrps and run the hmm
  // indices of columns of X to be selected:
  // rind is equivalent to the experiment index
  // this defines the stimuli to take from stimgrps
  for(int rind=0; rind!=*numexperimentsx; rind++) {
    //printf("replicates: %d\n", R[rind]);
    allR += R[rind]; // count the total number of replicates
    // extract the stimuli for the rind'th experiment
    int lstids = stimgrps[rind];
    int *stids = (int*)malloc(lstids*sizeof (int));

    for(int st=0; st!=lstids; st++) {
      //printf("STIMULUS ID: %d\n", stimids[idstart+st]);
      stids[st] = stimids[idstart+st];
    }
    //printf("lstids: %d  ", lstids);

    // remember the first stimulus from the next experiment
    idstart+=lstids;

    // how long is the current experiment, i.e. how many columns in the Gamma and X matrices
    lenofexp = R[rind] * T;

    // get the sub-data-matrix
    end=start+lenofexp;
    //printf("idstart %d; lenofexp: %d, end %d\n", idstart, lenofexp, end);
    // run the viterbi algorithm
    // Gstart contains the position in the G matrix at which new states have to be inserted
    Gstart = runviterbi(phi, N, T, R[rind], X, GS, start, end, G, Glen,
        stids, lstids, Gstart, hmmit);
    start+=lenofexp;
    free(stids);
  }

  // now the GS and G matrices are ready, get the final likelihood and TH values
  // update the theta matrix
  estimate_theta(X, GS, TH, N, T, allR);

  //printf("\t TH:\n");
  //print_matrix(TH, N, 4);

  // calculate the new likelihood
  double Lik = calculate_likelihood(X, GS, TH, N, T, allR); //Liktmp$L
  //printf("LIK::::: %f\n", Lik);
  return(Lik);
}

/*
 * run the viterbi algorithm for one experiment and one network
 */
int
runviterbi(int *phi, const int N, const int T, int R, const double *X, int *GS,
    int start, int end, int *G, int *Glen, int *stids,
    int lstids, int Gstart, const int hmmit)
{
  //printf("\nRUNVITERBI : start %d, end %d; N=%d; T=%d; R=%d\n", start, end, N,
  //    T, R);
  // perform state propagation to obtain the gamma matrix
  int lenofmat = N * T * R, M=0, k, runningindex=0; //, pseudocountsum;
  //printf("lenofmat: %d", lenofmat);

  // create the system state matrix
  int *Gsub = malloc(lenofmat * sizeof(int));
  int *GSsub = malloc(lenofmat * sizeof(int));

  // here the Gsub matrix is filled with the system states
  M = propagate_effect(phi, stids, lstids, Gsub, N);
  //pseudocountsum=M;
  //printf("BACK from propeffect: Found M=%d system states. Gstart: %d\n", M, Gstart);
  // free the unused memory in Gsub
  Gsub = realloc(Gsub, (M*N)*sizeof(int));
  //printf("\t Gsub:\n");
  //print_intmatrix(Gsub, N, M);

  // get a theta matrix
  double *THsub = malloc(4*N * sizeof(double));
  //TH = (double*)malloc(4*N * sizeof(double));

  // run the hmm
  hmm(X, Gsub, GSsub, THsub, N, T, R, M, hmmit);
  /*printf("RUNVITERBI:\n");
  print_intmatrix(GSsub, N, T*R);
  printf("RUNVITERBI GS:\n");
  print_intmatrix(GS, N, 2*T*R);
  */
  // start, end denote the start and end *columns*, get the indices:
  int GSstart = start * N;
  int GSend = end * N;
  // add to GS and G
  for(k=GSstart; k!=GSend; k++) {
      GS[k] = GSsub[runningindex];
      runningindex++;
  }
  //printf("RUNVITERBI GS:\n");
  //print_intmatrix(GS, N, 2*T*R);
  int temp = (Gstart+M*N);
  // allocate more space if G wasn't big enough before
  if(temp>*Glen) {
      //G = realloc(G, (Gstart+M*N) * sizeof(int));
      printf("Warning: Too little memory allocated for G!\n");
  }
  //    printf("start: %d, end: %d, GSstart: %d, GSend: %d, Gstart: %d, temp: %d\n",start, end,GSstart, GSend, Gstart, temp);
  runningindex = 0;
  for(k=Gstart; k!=temp; k++) {
      G[k] = Gsub[runningindex];
      runningindex++;
  }
  //printf("G, allM %d:\n",temp/N);
  //print_intmatrix(G, N, (temp/N));

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
  //for(int ind=0; ind!=M*T*R; ++ind) {
  //  viterbi[ind] = 0;
  //}
  int *maxemissionind = malloc(T * sizeof(int));

  double Lik = infinity, Lold;
  int diffold = -100, equally = 0, it = 0, restarts = 0, tmp;
  double diff, diffsold[hmmit];

  double Lnew = calculate_likelihood(X, GSsub, THsub, N, T, R);
  while (it <= hmmit)
    {
       it++;
      //old likelihood
      Lold = Lik;

      // update the theta matrix
      estimate_theta(X, GSsub, THsub, N, T, R);

      // calculate the new likelihood
      Lik = calculate_likelihood(X, GSsub, THsub, N, T, R); //Liktmp$L
      //printf("XXLik: %f ", fabs(Lik));
      //printf("XXLold: %f ", fabs(Lold));
      //printf("XXDIFF Lik Lold: %f\n", fabs((fabs(Lik) - fabs(Lold))));

      if (isinf(Lold) != 1)
        {
          if (fabs((fabs(Lik) - fabs(Lold))) <= 1.0)
            break;
          diff = fabs((fabs(Lik) - fabs(Lold)));
          diffsold[it] = diff;
          // check for repeating patterns
          tmp = find_diff(diff, diffsold, it);
          //printf("%d tmp: %d\n",it,tmp);
          //print_intmatrix(diffsold,1,it);
          //printf("%d DIFF: %d, EQUALLY: %d, RESTARTS: %d\n",it, diff, equally, restarts);
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
                      for(int ddd=0; ddd!=it; ++ddd) {
                        diffsold[ddd] = -1*it;
                      }
                      // again initialise the state matrix
                      //printf("A:\n");
                      //print_matrix(A, M, M);
                      //print_matrix_stats(A, M, M);
                      //srand(time(NULL) + getpid());
                      // alternative restart
                      //initialise_A(A, M);
                      initialise_GS(GSsub, Gsub, N, T, R, M);
                      //free(viterbi);
                      //free(maxemissionind);
                      //double *viterbi = malloc(M * T * R * sizeof(double));
                      //int *maxemissionind = malloc(T * sizeof(int));
                      //estimate_theta(X, GSsub, THsub, N, T, R);
                      /*it = 0;*/
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

      // fill the viterbi matrix
      //for(int ind=0; ind!=M*T*R; ++ind) {
      //  viterbi[ind] = 0;
      //}
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
              //printf("\n\t\tmaxemission: %f, viterbi: %f\n", maxemission, viterbi[ind]);
              if (viterbi[ind] >= maxemission)
                {
                  //printf("set maxemind\n");
                  maxemission = viterbi[ind];
                  maxemissionind[t] = m;
                }
            }
        }
     /* for(int xxx=0; xxx!=T; xxx++) {
          printf("%d ",maxemissionind[xxx]);
      }
*/
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
  //printf("ITERATIONS PERFORMED: %d / %d\n", it,hmmit);
  //print_intmatrix(GSsub, N, T*R);
  free(maxemissionind);
  free(viterbi);
  free(A);
}


void
updateA(double *A, int *maxemissionind, const int M, const int T)
{
  int *Anum = calloc(M * M, sizeof(int));
  int *Adenom = calloc(M, sizeof(int));
  //int *pseudocounts = calloc(M * M, sizeof(int));
  //int *pseudocountssum = calloc(M * M, sizeof(int));
  int t, src, dst;

  // go along maxemissionind and count the transitions
  for (t = 0; t != (T - 1); t++)
    {
      src = maxemissionind[t];
      dst = maxemissionind[t + 1];
      //printf("trans %d_%d\n",src,dst); //,Anum[ind]);
      int ind = dst*M + src;
      Anum[ind] = Anum[ind] + 1;
      //pseudocounts[ind] = 1;
      //pseudocountssum[ind] = M;
      Adenom[src] = Anum[src] + 1;
    }
  /*  printf("Anum:\n");
  print_intmatrix(Anum, M, M);
  printf("Adenom:\n");
  print_intmatrix(Adenom, 1, M);
  printf("A:\n");
  print_matrix(A, M, M);
  print_matrix_stats(A, M, M);
*/
  // reset the transition probability to numocc / numtotaltransitionsfromsrc
  double rowsum;
  for (t = 0; t != (M * M); t++)
    {
      rowsum = 0.0;
      // found the transition
      if (Anum[t] != 0)
        {
          A[t] = log2(Anum[t] + 1) - log2(Adenom[(t % M)] + M);
        }
      //Asum += pow(2, A[t]);
    }
  // is a wrong normalisation...
  double Asum = 0.0;
  // record the rowSums for normalisation: tranisition matrix
  // is always normalised such that the rows sum to 1
  double *Asums = calloc(M * M, sizeof(double));
  double *rowsums = calloc(M, sizeof(double));
  Asum = get_rowsums(A, M, M, rowsums);
  //printf("%%%%%%%%%%\nA before: \n");
  //print_matrix(A, M, M);
  //print_matrix_stats(A, M, M);
  // normalise the rows to sum to 1
  normalise_rows(A, M, M, rowsums);
  //printf("&&&&&&&&&&\nA after: \n");
  //print_matrix(A, M, M);
  //print_matrix_stats(A, M, M);
  //free(pseudocountssum);
  //free(pseudocounts);
  free(Anum);
  free(Adenom);
  free(Asums);
}
/*

  // now renormalise such that sum(A) = 1
  for (t = 0; t != (M * M); t++)
    {
      // found the transition
      A[t] = A[t] - log2(Asum);
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
  //double *rowsums = calloc(rows * cols, sizeof(double));
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
  double *Asums = calloc(M * M, sizeof(double));
  double *rowsums = calloc(M, sizeof(double));
  double Asum = get_rowsums(A, M, M, rowsums);

 /* printf("-------%%%%%%%%%%\nA before: \n");
   print_matrix(A, M, M);
   //print_matrix_stats(A, M, M);
   normalise_rows(A, M, M, rowsums);
   printf("--------&&&&&&&&&&\nA after: \n");
   print_matrix(A, M, M);
   //print_matrix_stats(A, M, M);
*/
  /* normalise to sum 1
  for (int i = 0; i != M; ++i)
    {
      for (int j = 0; j != M; ++j)
        {
          A[(i * M + j)] = log2(A[(i * M + j)] / rnumsum);
        }
    }
   */
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
          //double ret = dnorm(dp, mu, sd);
          //double ret = dn1 - dn2;
          double ret = dn1;
          if(isnan(ret))
            ret = 0.0;
          emission += ret;
        }
    }
  return (emission);
}

/*int
find_double_diff(double diff, double *diffsold, int maxit)
{
  int occ = 0;
  for (int i = 0; i != maxit; i++)
    {
      if (diffsold[i] == diff)
        occ++;
    }
  return (occ);
}
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
      //printf("Gsub[Gstart+i]: %d, binsum: %d\n", Gsub[i], (int)(Gsub[i] * pow(2, i)));
      binsum += Gsub[i] * pow(2, i);
    }

  binsums[it++] = binsum;
  //printf("binsums[it]: %d",binsums[(it-1)]);
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
          // print_intmatrix(Gsub, N, it + 2);
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
          //printf("newcol[i]: %d\n", newcol[i]);
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
      //printf("newcol[i]: %d\n", newcol[i]);
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
 * get the basis 2 log
 */
double
log2(double x)
{
  double ret = log(x) / log(2);
  return (ret);
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

/*
 * print row sums and sum of a matrix
 */
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
/*
 * print a data matrix, passed as pointer
 */
void
print_matrix(double *x, int rows, int cols)
{
  //printf("X:\n\n");
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
  //printf("X:\n\n");
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
          if (GSsub[kk] == 1)
            {
              a += X[kk];
              a_n++;
            }
          else
            {
              p += X[kk];
              p_n++;
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
          if (GSsub[kk] == 1)
            {
              tmp = X[kk] - mua;
              a += (tmp * tmp);
              a_n++;
            }
          else
            {
              tmp = X[kk] - mup;
              p += (tmp * tmp);
              p_n++;
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
/*  printf("\t GSsub:\n");
  print_intmatrix(GSsub, N, T*R);
  printf("N: %d, T: %d, R: %d",N,T,R);
  printf("\t TH:\n");
  print_matrix(THsub, N, 4);
*/  double dp, mu, sd, Lnew;
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
          //printf("%f ",dn);
          //if(i==4){
          //  printf("%f %f %f %f %d\n",dn,mu,sd, dp, gp);
          //}
          if (!isnan(dn) && !isinf(dn))
            {
              if (!isnan(penalty) && !isinf(penalty))
                {
                  dn -= penalty;
                }
              Lnew += dn; // dnorm is logged data
            }
        }
//      printf("\n");
    }
//  printf("Liklihood: %f", Lnew);
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
  //density = log(density) / log(2); // log to base 2
  density = log(density); // log to base 10
  //printf("dp: %f, mu: %f, sd: %f, density %f\n\n", dp, mu, sd, density);
  return (density);
}
