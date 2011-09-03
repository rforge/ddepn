/* perform_hmmsearch_globalest.h */
#ifndef UTILITY_FUNCTIONS_H__
#define UTILITY_FUNCTIONS_H__

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

double dnorm(double dp, double mu, double sd);

void estimate_theta(const double *X, int *GSsub, double *TH,
		const int N, const int T, const int R); //, int M);

double getE(const double *X, int t, int *Gsub, int m, double *THsub,
		const int N, const int R);

//int get_leafs(int *phi, int N, int *leafs);

double get_rowsums(double *x, int rows, int cols, double *rsums);

double log2(double x);

void normalise_rows(double *x, int rows, int cols, double *rowsums);

void print_intmatrix(int *x, int rows, int cols);

void print_matrix(double *x, int rows, int cols);

void print_matrix_stats(double *x, int rows, int cols);

void updateA(double *A, int *maxemissionind, const int M, const int T);


/*
 * OLD functions, not used anymore
 *
int is_stimulus(int x, int *stids, int lstids);


int next_col(int *newcol, int *Gsub, const int N, int *stids, int lstids, int *phi, int index);


void initialise_A(double *A, int M);


int find_diff(double diff, double *diffsold, int maxit);


int check_binsum(int *binsums, int binsum, int it);

int compare_doubles(const void *a, const void *b);

int compare_ints(const void *a, const void *b);



void initialise_GS(int *GSsub, int *Gsub, const int N, const int T, const int R, int M);

int find_bin(int *Ms, int num, int numexperiments);

/ *
 * Effect Propagation for one stimulus. Generates the
 * reachable states Gamma
 * /
int propagate_effect(int *phi, int *stids, int lstids, int *Gsub, const int N);

/ *
 * Effect propagation wrapper
 * /
//void propagate_effect_wrap(int *phix, int *stidsx, int *lstidsx, int *Gsubx, int *Nx, int *Mx);

/ *
 * Run viterbi for one stimulus on a submatrix of data X
 * /
int runviterbi(int *phi, const int N, const int T, int R,
    const double *X, int *GS, int start, int end,
    int *G, int *Glen,
    int *stids, int lstids,
    int Gstart, const int hmmit);


void hmm(const double *X, int *Gsub, int *GSsub, double *TH, const int N, const int T,
    const int R, int M, const int hmmit);
 */
#endif /* UTILITY_FUNCTIONS_H__ */
