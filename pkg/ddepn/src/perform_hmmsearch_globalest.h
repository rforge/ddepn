/* perform_hmmsearch_globalest.h */
#ifndef PERFORM_HMMSEARCH_GLOBALEST_H__
#define PERFORM_HMMSEARCH_GLOBALEST_H__

//#define HMMITERATIONS 100
//#define RAND_MAX 10000
#ifndef max
#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#ifndef min
#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif

void checkC(int *phi, int *GS, int ch, int pa, int N, int T, int R, int *consvec);

void extract_transitionmatrix(double *src, double *dst, int nstates, int index);

void extract_datamatrix(const double *X, double *Xexp, int N, int T, int R, int start);

void extract_statematrix(int *X, int *Xexp, int N, int cols, int start);

int find_switchable(double *TH, int N, int *toswitch);

//int get_leafs(int *phi, int N, int *leafs);

int get_parents(int *phi, int N, int *parents, int leaf);

/*
 * Perform Viterbi Training for a single network phi
 */
double hmmsearch(int *phi, const int N,
		const int *Tx, const int *Rx,
		const double *X, int *GS,
		int *G, int Glen, double *TH,
		const int *tps,
		const int *stimids, const int *stimgrps,
		const int numexperimentsx, const int hmmit, int *Ms);

void init_A(double *A, int *Ms, int numexperiments);

int is_consistent(int *phi, int *GS, int *G, int N, int T, int R);


/*
 * Outer call to perform hmm-Viterbi training
 * For each network phi \in Px, perform the Viterbi
 */
void perform_hmmsearch_globalest(int *Px, int *Nx, int *Tx, int *Rx,
		double *Xx, int *GSx,
		int *Gx, int *Glenx, double *THx,
		int *tpsx,
		int *stimidsx, int *stimgrpsx,
		int *numexperimentsx, double *Likx,
		int *hmmiterations, int *Msx);

void print_intmatrix(int *x, int rows, int cols);

void print_matrix(double *x, int rows, int cols);

void print_matrix_stats(double *x, int rows, int cols);

/*void print_transitions(double *A, int *Ms, int numexperiments, int Msum);
*/
void switch_theta_row(double *TH, int rnum, int N);

void updateA(double *A, int *maxemissionind, const int M, const int T);

void update_statematrix(int *GS, int *GSexp, int start, int N, int T, int R);

void update_transitionmatrix(double *src, double *dst, int nstates, int index);

void viterbi(int T, int M, double *X, int *G, double *TH, int N, int R, double *A, int *GS);



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
#endif /* PERFORM_HMMSEARCH_GLOBALEST_H__ */
