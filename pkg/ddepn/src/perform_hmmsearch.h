/* perform_hmmsearch.h */
#ifndef PERFORM_HMMSEARCH_H__
#define PERFORM_HMMSEARCH_H__

//#define HMMITERATIONS 100
//#define RAND_MAX 10000
#ifndef max
        #define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#ifndef min
        #define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif

double calculate_likelihood(const double *X, int *GSsub, double *THsub,
		const int N, const int T, const int R);

int check_binsum(int *binsums, int binsum, int it);

int compare_doubles(const void *a, const void *b);

int compare_ints(const void *a, const void *b);

double dnorm(double dp, double mu, double sd);

void estimate_theta(const double *X, int *GSsub, double *TH,
    const int N, const int T, const int R); //, int M);

int find_diff(double diff, double *diffsold, int maxit);

double getE(const double *X, int t, int *Gsub, int m, double *THsub,
     const int N, const int R);

void hmm(const double *X, int *Gsub, int *GSsub, double *TH, const int N, const int T,
    const int R, int M, const int hmmit);

/*
 * Perform Viterbi Training for a single network phi
 */
double hmmsearch(int *phi, const int N,
    const int T, const int *R,
    const double *X, int *GS,
    int *G, int *Glen, double *TH,
    const int *tps,
    const int *stimids, const int *stimgrps,
    const int *numexperimentsx, const int hmmit);

void initialise_A(double *A, int M);

void initialise_GS(int *GSsub, int *Gsub, const int N, const int T, const int R, int M);

int is_stimulus(int x, int *stids, int lstids);

double log2(double x);

int next_col(int *newcol, int *Gsub, const int N, int *stids, int lstids, int *phi, int index);

/*
 * Outer call to perform hmm-Viterbi training
 * For each network phi \in Px, perform the Viterbi
 */
void perform_hmmsearch(int *Px, int *Nx, int *Tx, int *Rx,
                       double *Xx, int *GSx,
                       int *Gx, int *Glenx, double *THx,
                       int *tpsx,
                       int *stimidsx, int *stimgrpsx,
                       int *numexperimentsx, double *Likx,
                       int *hmmiterations);

void print_matrix_stats(double *x, int rows, int cols);

double get_rowsums(double *x, int rows, int cols, double *rsums);
void normalise_rows(double *x, int rows, int cols, double *rowsums);

void print_matrix(double *x, int rows, int cols);

void print_intmatrix(int *x, int rows, int cols);

/*
 * Effect Propagation for one stimulus. Generates the
 * reachable states Gamma
 */
int propagate_effect(int *phi, int *stids, int lstids, int *Gsub, const int N);

/*
 * Run viterbi for one stimulus on a submatrix of data X
 */
int runviterbi(int *phi, const int N, const int T, int R,
    const double *X, int *GS, int start, int end,
    int *G, int *Glen,
    int *stids, int lstids,
    int Gstart, const int hmmit);

void updateA(double *A, int *maxemissionind, const int M, const int T);

#endif /* PERFORM_HMMSEARCH_H__ */
