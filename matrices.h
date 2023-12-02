/* file_start */

/* This holds matrices in column major order */
struct matrix
{
  int rows;
  int cols;
  double *mtrx;
};

/* Here are the matrix functions */
struct matrix * mmake();
void minit0(struct matrix *,int,int);
void minit1(struct matrix *,int,int);
void minitr(struct matrix *,int,int);
void sentry(struct matrix *,int,int,double);
double gentry(struct matrix *,int,int);
void mcopy(struct matrix *,struct matrix *);
void mplugin(struct matrix *,int,int,int,int,struct matrix *,int,int);
void mplugout(struct matrix *,int,int,int,int,struct matrix *);
double mmax(struct matrix *);
double mmin(struct matrix *);
void mtrans(struct matrix *,struct matrix *);
void mtimesm(struct matrix *,struct matrix *,struct matrix *);
void mtimesm_block(struct matrix *,struct matrix *,struct matrix *,int,int);
void minteg_block(struct matrix *,struct matrix *,struct matrix *,int,int);
void mtimesd(struct matrix *,double,struct matrix *);
void mplusm(struct matrix *,struct matrix *,struct matrix *);
void mplusm_block(struct matrix *,struct matrix *,struct matrix *,int,int);
void msubtm(struct matrix *,struct matrix *,struct matrix *);
void mplusd(struct matrix *,double,struct matrix *);
void mfeval(struct matrix *, double (*)(double),struct matrix *);
void mvnoc(struct matrix *,struct matrix *,struct matrix *);
void mpint(struct matrix *,struct matrix *);
void mpder(struct matrix *,struct matrix *);
double mpeval(struct matrix *,double);
void mplu(struct matrix *,struct matrix *,struct matrix *);
void mplusolve(struct matrix *,struct matrix *, struct matrix *,struct matrix *);
void mplusolve_block(struct matrix *,struct matrix *,struct matrix *,struct matrix *,int,int);
void minv(struct matrix *,struct matrix *);
double mfrobnorm(struct matrix *);
void mqrh(struct matrix *,struct matrix *);
double mkget(struct matrix *);
double mkest(struct matrix *);
void mtfree(struct matrix *);
void mprintf(struct matrix *);
void mspy(struct matrix *,double);
void mplotf(struct matrix *,int,int);
void aprintf(struct matrix *);

/* file_end */

