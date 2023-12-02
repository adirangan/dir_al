/* This holds a lgn */
struct lgn
{
  char name[256]; /* name of file */
  int rows; /* rows in picture */
  int cols; /* columns in picture */
  double height; /* height in degrees */
  double width; /* width in degrees */
  double irad; /* inner radius (degrees) for lgn center, outer radius is 3 times inner radius */
  double lgnbaserate; /* lgn base background firing rate */
  double peakrate; /* peak firing rate (increase) of impulse response */
  fftw_plan plan_forward; /* forward plan for fast fft */
  fftw_plan plan_backward; /* backward plan for fast fft */
  double *pnm; /* holds picture data */
  double *onshape; /* spatial kernel for on-center cells - an extra 4 array elements hold stats */
  double *onin; /* summed intensity input to lgn on-center cells */
  double *onrate; /* used to determine lgn on-center firing rate */
  double *offrate; /* used to determine lgn off-center firing rate */
  double **angleshape; /* spatial kernel for gabor of a certain angle and phase - an extra 4 array elements (per phase,angle) hold stats */
  double **anglein; /* linear input to gabors of various angles - an extra 4 array elements (per phase,angle) hold stats */
  double **anglerate; /* used to determine various angular firing rates */
};

/* Here are the lgn functions */
struct lgn * lgnmake(char *,int,int,double,double,double,double,double);
void makegrating(struct lgn *,int);
void lgnworkhorse(struct lgn *);
double lgnonshapekernel(double,double,double,double,double);
int lgndump(struct lgn *);
int lgnread(struct lgn *);
void lgnconvolve(struct lgn *);
void lgnremake(struct lgn *,char *);
void lgnrefresh(struct lgn *,double,double);
void lgnevolve(struct lgn *,double,double);
void lgnresetrates(struct lgn *);
void lgntfree(struct lgn *);
int framefig(double);
