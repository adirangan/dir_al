/* This holds a neuron state */
struct neuron
{
  double warning; /* -2 = spike that spikes, -1 = spike that doesn't, 0 = gizmo that doesn't, 1 = gizmo that spikes, 2 = bad *(n->V) */
  int row;
  int col;
  int pierow;
  int piecol;
  int inpierow;
  int inpiecol;
  double rad;
  int ang;
  int type;
  int sog;
  double *t2s;
  int sox;
  double lgnangle;
  double lgnangle_drift;
  double lgnphase;
  double lgndriftx;
  double lgndrifty;
  double sV1;
  double sV2;
  double sA1;
  double sA2;
  double sN1;
  double sN2;
  double sG1;
  double sG2;
  double g;
  double spiketime_guess;
  int spiketime_guess_flag;
  double inputrate;
  int spikeinput_flag;
  int spikeinput_multiplicity;
  long long int spikeinput_rseed;
  double spikeinput_time;
  double spikelast;
  double spiketime;
  double spikenext;
  int sparse_out;
  int sparse_in;
  double homogenizing_knob; /* between 0 (homogenous) to 1 (strictly sparse) */
  struct llitem *sparse_link;
  int lr_kernel_type;
  double *V;
  double *sA;
  double *sN;
  double *sG;
  double *VS;
};

/* This holds multiple neurons */
struct neuronarray
{
  int rows;
  int cols;
  double *Vra;
  double *sAra;
  double *sNra;
  double *sGra;
  double *VSra;
  double *t2sra;
  struct neuron ** N;
  int *t2stotal;
  int *natotal;
  double mIC;
  double mIS;
  double mEC;
  double mES;
};

/* here are cortex-dependent llists functions */
void lrupdate(struct neuronarray *,double,double);
int spikelast_compare(void *,void *);
int spiketime_compare(void *,void *);
int lgnangle_compare(void *,void *);
void llistunmake(struct llist **,struct llist **,struct llist *);
void llistremake(struct llist ***,struct llist ***,struct llist **);

/* Here are the neuron/neuronarray functions */
void makeshuffle(int *,int,int,int);
fftw_complex * makelrkernel(int);
void neuronsetang(struct neuronarray *);
void neuronmake(struct neuronarray *,int,int);
void neurontfree(struct neuron *);
struct neuron * nget(struct neuronarray *,int,int);
void nset(struct neuronarray *,int,int,void *);
struct neuronarray * neuronarraymake(int,int);
void neuronarraytfree(struct neuronarray *);
void setblocksize(double);;
void setinputrate(struct neuronarray *,double);
void spikeinput(struct neuronarray *,double,double);
int spikescan(struct neuronarray *,struct llist **,struct llist **,double,double);
double spikeguesseif(struct neuron *,double,double);
double spikeguess2(struct neuron *,double,double);
double spikeguess1(struct neuron *,double,double);
void spikesort(struct llist **,struct llist *);
void spikesort_tagchange(struct llist **,int,int,int,int,struct llist *,struct llist *);
void llprintf(struct llist *);
void llraprintf(struct llist **,int,int);
void lreconnect(struct llist *,struct llist *);
int distance(struct llist *,struct llist *);
int ilink(struct neuron *,double *,double *,double *);
int slink(struct neuron *,struct neuron *,double *,double *,double *,double *,double *,double *);
double ifrhs(double,double,double,double,double);
double eifrhs(double,double,double,double,double);
void slavecalc(double,double,double,double *,double *,double *);
void gizmointegrate2(struct neuron *,double,double,double,double,double);
void gizmointegrateeif(struct neuron *,double,double,double,double,double);
void clumpcorrect(struct llist *,struct llist *,double,double);
void spikecorrect(struct llist *,double,double);
void spikeconduct(struct llist *,double,double,int *,int *,int *,int *);
double quadrootfinder(double,double,double,double,double);
double cuberootfinder(double,double,double,double,double,double);
void gizmoconduct(struct llist **,struct llist **,double,double);
void gizmospiked(struct neuron *,double,double);
