/* measure the difference between an integrate-and-fire neuron and an exponential integrate-and-fire-neuron */
/* used for error analysis of single neuron integrate-and-fire solver */
#ifndef NOTGRAPH
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#endif /* NOTGRAPH */
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <fftw3.h>

#define ESCAPE_KEY 27
#define SPACE_KEY 32
#define ENTER_KEY 10
#define PI 3.141592653589793
#define maximum(A,B) ((A) > (B) ? (A) : (B))
#define minimum(A,B) ((A) < (B) ? (A) : (B))
#define periodize(A,B,C) ((A) < (B) ? (A) + (C) - (B) : ((A) >= (C) ? (A) - (C) + (B) : (A)))
#define crop(A,B,C) ((A) < (B) ? (B) : ((A) > (C) ? (C) : (A)))
#define rand01 ((double)rand()/(double)RAND_MAX)
#define entry(j,k,n) ((j) + (k)*(n)) /* This defines index for entry (j,k) in row major order (not 'fortran') n by m matrix */

#define RGB3 3 // 3 bytes of color info per pixel
#define RGBA 4 // 4 bytes of color+alpha info
#define UCHAR_MAX 255 // is this right?

#include "d_llists.h"
#include "matrices.h"

/* This holds a rho0 */
struct rho0
{
  int nv;
  double *midv;
  double *midv_rho;
  double *fluxv;
  double *coef_jumpin;
  double *coef_jumpout;
  double *rho;
  double sum;
  double firingrate;
  struct matrix *Ginv;
  struct matrix *Linv;
  struct matrix *G;
  struct matrix *L;
  struct matrix *L_base;
  double *Lprimerho;
  double *rho_tilde;
};

/* This holds a rho1 */
struct rho1
{
  int nv;
  int ng;
  int ns;
  double *midv;
  double *midg;
  double *midv_rho;
  double *midg_rho;
  double *fluxv;
  double *fluxg;
  double *coef_jumpin;
  double *coef_jumpout;
  double *rho;
  double sum;
  struct matrix *Ginv;
  struct matrix *Linv;
  struct matrix *G;
  struct matrix *L;
};

/* This holds a neuron */
struct neuron
{
  int types; /* number of different types */
/*   double fA; double nuA;  */
/*   double fN; double nuN;  */
/*   double fG; double nuG;  */
  struct llist *sA;
  struct llist *sN;
  struct llist *sG;
  struct llist **V; /* various instantaneous voltages */
  struct llist **spikenext; /* various spiketimes */
  struct llist **t; /* various times */
  double *sra; /* total spikes */
  double *Vra; /* mean voltages */
  double *tra; /* total time */
  double *hist; /* v,sA histogram, stored as [vbin + gbin*nbins + type*nbins*nbins] */
  int nbins;
};

/* Here are method functions */

struct rho0 *rho0make(int);
double rho0_vdot(double);
void rho0_getmidpoints(struct rho0 *,int);
void rho0_getjump(struct rho0 *);
void rho0_diffuse(struct rho0 *);
void rho0update(struct rho0 *,double,double);
void rho0_getrhotilde(struct rho0 *);
struct rho1 *rho1make(int,int);
double rho1_vdot(double,double);
double rho1_gdot(double,double);
void rho1_getmidpoints(struct rho1 *,int);
void rho1_getjump(struct rho1 *);
void rho1update(struct rho1 *,double,double);
struct neuron * neuronmake(int/*, double,double,double,double,double,double */);
void neurontfree(struct neuron *);
double quadrootfinder(double,double,double,double,double);
void slavecalc(double,double,double,double *,double *,double *);
int ifevolve(double,double,double,double,double,double,double,double,double *,double *);
int eifworkhorse(double,double,double,double,double,double,double *,double *);
int eifevolve(double,double,double,double,double,double,double,double,double *,double *);
void llistbanzai(struct llist *);
void neuronevolve(struct neuron *,double);
void setglobals();
void computestep();
#ifndef NOTGRAPH
int WritePPMFile(const char *,GLubyte *,int,int,int);
int DumpWindow(const char *,int,int);
void ftexto(float,float,float,char *);
void specialKeyPressed(int , int , int );
void keyPressed(unsigned char , int , int ); 
GLvoid glbox(double,double,double,double,double,double);
GLvoid Drawmenu(double,double,double);
GLvoid Drawneuron(double,double,double);
GLvoid DrawGLScene(GLvoid);
GLvoid ReSizeGLScene(GLsizei , GLsizei );
GLvoid InitGL(GLsizei, GLsizei);
#endif /* NOTGRAPH */

/* Here are the global functions */
double drand(int);
int processoptions(int, char **);

/* These are the real global variables */
int MEMDEBUG=0;
void * MEMDEBUGLIST=NULL;
#ifndef NOTGRAPH
int GLUTWINDOWNUMBER=0;
int GLOBAL_WINDOW_WIDTH=0;
int GLOBAL_WINDOW_HEIGHT=0;
#endif /* NOTGRAPH */

/* These are convenient global variables */
#ifndef NOTGRAPH
int GLOBAL_DRAW_FLAG=0;
GLfloat xdepth = 0.0f;
GLfloat ydepth = 0.0f;
GLfloat zdepth = 0.0f;
#endif /* NOTGRAPH */
int ON_MY_COMPUTER=1;
int GRAYSCALE=0;
int DRAW_FLAG=0;
double STD_VIEW=1.0;
double AMPA_RATE= 1250.0/1000.0;//200.0/1000.0;
double NMDA_RATE= 0;//200.0/1000.0;
double GABA_RATE= 0;
double AMPA_STRENGTH= 0.05;
double NMDA_STRENGTH= .01;
double GABA_STRENGTH= 10.0;
double CONDUCTANCE_LK=0.05;
double CONDUCTANCE_MAX=0.05;
double VOLTAGE_RESET=0;//-.375;
double VOLTAGE_LEAK=0;
double VOLTAGE_TAKEOFF=0.625;
double VOLTAGE_THRESHOLD_IF=1.0;
double VOLTAGE_THRESHOLD_EIF=4.375;
double VOLTAGE_EX=4.6666667;
double VOLTAGE_IN=-0.6666667;
double TAU_REF=2;//0.125;
double TAU_A=5;
double TAU_N=80;
double TAU_G=7;
double TAU_SPIKE=.0037;
double GHAT=0;
struct neuron *GLOBAL_NEURON=NULL;
int MAXLENGTH=1024*16;
double GLOBAL_time=0;
double GLOBAL_DT=1.0/1024.0;
double GLOBAL_TF=0;
int STEPS_PER_DRAW=0;
int FIDDLE_PARAMETER=0;
int RUN_DONE=0;
int GLOBAL_ARBOR_SIM=0;
int GLOBAL_NTYPES=1;
int GLOBAL_TRUETYPE=0;
int GLOBAL_RHO_VERSION=0;
struct rho1 *GLOBAL_rho1=NULL;
int GLOBAL_rho1_nv=8;
int GLOBAL_rho1_ng=12;
struct rho0 *GLOBAL_rho0=NULL;
int GLOBAL_rho0_nv=256;

#include "d_llists.c"
#include "matrices.c"

/* Here are the rho0 functions */

struct rho0 *rho0make(int rho0_nv)
{
  int verbose=0;
  int nv=0;
  int nv1=0,nv2=0;
  double sum=0;
  struct rho0 *r0=NULL;
  struct matrix *tmp1=NULL;
  r0 = (struct rho0 *) tcalloc(1,sizeof(struct rho0));
  r0->nv=rho0_nv;
  r0->midv = (double *) tcalloc(r0->nv+1,sizeof(double));
  r0->midv_rho = (double *) tcalloc(r0->nv+1,sizeof(double));
  r0->fluxv = (double *) tcalloc(r0->nv+1,sizeof(double));
  r0->coef_jumpin = (double *) tcalloc(r0->nv,sizeof(double));
  r0->coef_jumpout = (double *) tcalloc(r0->nv,sizeof(double));
  r0->rho = (double *) tcalloc(r0->nv,sizeof(double));
  rho0_getmidpoints(r0,1);
  r0->Ginv = mmake();
  minit0(r0->Ginv,r0->nv,r0->nv);
  for (nv1=0;nv1<r0->nv;nv1++){
    for (nv2=0;nv2<r0->nv;nv2++){ if (nv2==nv1){ r0->rho[nv2]=1;} else /* if (nv2!=nv1) */{ r0->rho[nv2]=0;}}
    rho0_getmidpoints(r0,0); 
    rho0_getjump(r0);
    for (nv2=0;nv2<r0->nv;nv2++){
      sentry(r0->Ginv,nv2,nv1,(r0->fluxv[nv2] - r0->fluxv[nv2+1]) + AMPA_RATE*(r0->coef_jumpin[nv2] - r0->coef_jumpout[nv2]));}}
  r0->Linv = mmake();
  minit1(r0->Linv,r0->nv,r0->nv);
  tmp1=mmake();mtimesd(r0->Ginv,GLOBAL_DT,tmp1);
  msubtm(r0->Linv,tmp1,r0->Linv);
  mtfree(tmp1);tfree(tmp1);
  r0->L=mmake();
  minv(r0->Linv,r0->L);
  r0->L_base=mmake(); mcopy(r0->L,r0->L_base);
  r0->Lprimerho = (double *) tcalloc(r0->nv,sizeof(double));
  r0->rho_tilde = (double *) tcalloc(r0->nv,sizeof(double));
/*   mprintf(r0->L); */
/*   mspy(r0->Ginv,0.0000001); */
/*   mspy(r0->Linv,0.0000001); */
/*   mspy(r0->L,0.0000001); */
/*   exit(0); */
  //for (nv=0;nv<r0->nv;nv++){ r0->rho[nv] = exp(-pow(nv-r0->nv/2,2));}
  for (nv=0;nv<r0->nv;nv++){ r0->rho[nv] = 1;}
  stats("double",r0->rho,r0->nv,NULL,NULL,&sum,NULL); sum *= r0->nv; ratimesequals(r0->rho,r0->nv,1.0/sum);
  if (verbose){ raprintf(r0->rho,"double",1,r0->nv,"rho: ");}
  stats("double",r0->rho,r0->nv,NULL,NULL,&(r0->sum),NULL); r0->sum *= r0->nv;
  rho0_getmidpoints(r0,1);
  rho0_getjump(r0);
  return r0;
}

double rho0_vdot(double v){ return -CONDUCTANCE_LK*(v-VOLTAGE_RESET);}

void rho0_getmidpoints(struct rho0 *r0,int setv_flag)
{
  /* obtains midpoint values for rho0 in between cells
     also obtains fluxes --- just for kicks. */
  int verbose=0;
  int nv=0;
  int tab_v1=0,tab_v2=0; double coef_v1=0,coef_v2=0;
  int tab_ll=0,tab_l=0,tab_r=0,tab_rr=0; double coef_ll=0,coef_l=0,coef_r=0,coef_rr=0;
  double hv=0;
  int stencil_number=2;
  hv = (VOLTAGE_THRESHOLD_IF-VOLTAGE_RESET)/(double)r0->nv;
  if (setv_flag){
    for (nv=0;nv<r0->nv+1;nv++){ r0->midv[nv] = VOLTAGE_RESET + (VOLTAGE_THRESHOLD_IF - VOLTAGE_RESET)*(double)nv/(double)r0->nv;}
    if (verbose){ raprintf(r0->midv,"double",1,r0->nv+1,"midv:");}}
  switch (stencil_number){
  case 0: /* 2 point symmetric stencil --  linear */
    for (nv=0;nv<r0->nv+1;nv++){
      if (nv==0){ tab_v1 = 0; tab_v2 = 1; coef_v1 = 1.5; coef_v2 = -0.5;}
      else if (nv==r0->nv){ tab_v1 = r0->nv-2; tab_v2 = r0->nv-1; coef_v1 = -0.5; coef_v2 = 1.5;}
      else /* if (nv>0 && nv<r0->nv) */{ tab_v1 = nv-1; tab_v2 = nv; coef_v1 = 0.5; coef_v2 = 0.5;}
      r0->midv_rho[nv] = coef_v1*r0->rho[tab_v1] + coef_v2*r0->rho[tab_v2];}
    if (verbose){ raprintf(r0->midv_rho,"double",1,r0->nv+1,"midv_rho: ");}
    break;
  case 1: /* 4 point symmetric stencil -- linear least squares*/
    for (nv=0;nv<r0->nv+1;nv++){
      if (nv==0){ tab_ll=0; tab_l=1; tab_r=2; tab_rr=3; coef_ll=17.0/20.0; coef_l=9.0/20.0; coef_r=1.0/20.0; coef_rr=-7.0/20.0;}
      else if (nv==1){ tab_ll=0; tab_l=1; tab_r=2; tab_rr=3; coef_ll=11.0/20.0; coef_l=7.0/20.0; coef_r=3.0/20.0; coef_rr=-1.0/20.0;}
      else if (nv>1 && nv<r0->nv-1){ tab_ll=nv-2; tab_l=nv-1; tab_r=nv; tab_rr=nv+1; coef_ll=5.0/20.0; coef_l=5.0/20.0; coef_r=5.0/20.0; coef_rr=5.0/20.0;}
      else if (nv==r0->nv-1){ tab_ll=r0->nv-4; tab_l=r0->nv-3; tab_r=r0->nv-2; tab_rr=r0->nv-1; coef_ll=-1.0/20.0; coef_l=3.0/20.0; coef_r=7.0/20.0; coef_rr=11.0/20.0;}
      else if (nv==r0->nv){ tab_ll=r0->nv-4; tab_l=r0->nv-3; tab_r=r0->nv-2; tab_rr=r0->nv-1; coef_ll=-7.0/20.0; coef_l=1.0/20.0; coef_r=9.0/20.0; coef_rr=17.0/20.0;}
      r0->midv_rho[nv] = coef_ll*r0->rho[tab_ll] + coef_l*r0->rho[tab_l] + coef_r*r0->rho[tab_r] + coef_rr*r0->rho[tab_rr];}
    if (verbose){ raprintf(r0->midv_rho,"double",1,r0->nv+1,"midv_rho: ");}
    break;
  case 2: /* 4 point symmetric stencil -- cubic */
    for (nv=0;nv<r0->nv+1;nv++){
      if (nv==0){ tab_ll=0; tab_l=1; tab_r=2; tab_rr=3; coef_ll=35.0/16.0; coef_l=-35.0/16.0; coef_r=21.0/16.0; coef_rr=-5.0/16.0;}
      else if (nv==1){ tab_ll=0; tab_l=1; tab_r=2; tab_rr=3; coef_ll=5.0/16.0; coef_l=15.0/16.0; coef_r=-5.0/16.0; coef_rr=1.0/16.0;}
      else if (nv>1 && nv<r0->nv-1){ tab_ll=nv-2; tab_l=nv-1; tab_r=nv; tab_rr=nv+1; coef_ll=-1.0/16.0; coef_l=9.0/16.0; coef_r=9.0/16.0; coef_rr=-1.0/16.0;}
      else if (nv==r0->nv-1){ tab_ll=r0->nv-4; tab_l=r0->nv-3; tab_r=r0->nv-2; tab_rr=r0->nv-1; coef_ll=1.0/16.0; coef_l=-5.0/16.0; coef_r=15.0/16.0; coef_rr=5.0/16.0;}
      else if (nv==r0->nv){ tab_ll=r0->nv-4; tab_l=r0->nv-3; tab_r=r0->nv-2; tab_rr=r0->nv-1; coef_ll=-5.0/16.0; coef_l=21.0/16.0; coef_r=-35.0/16.0; coef_rr=35.0/16.0;}
      r0->midv_rho[nv] = coef_ll*r0->rho[tab_ll] + coef_l*r0->rho[tab_l] + coef_r*r0->rho[tab_r] + coef_rr*r0->rho[tab_rr];}
    if (verbose){ raprintf(r0->midv_rho,"double",1,r0->nv+1,"midv_rho: ");}
    break;
  default: printf(" %% warning, incorrect stencil type in rho0_getmidpoints\n");}
  for (nv=0;nv<r0->nv+1;nv++){ r0->fluxv[nv] = rho0_vdot(r0->midv[nv])*r0->midv_rho[nv]/hv;}
  if (r0->fluxv[r0->nv]>0){ r0->fluxv[0] = r0->fluxv[r0->nv];}
  else /* if (r0->fluxv[r0->nv]<=0) */{ r0->fluxv[0] = 0; r0->fluxv[r0->nv] = 0;}
  if (verbose){ raprintf(r0->fluxv,"double",1,r0->nv+1,"fluxv: ");}
}

void rho0_getjump(struct rho0 *r0)
{
  int verbose=0;
  int nv=0;
  double jumpval = AMPA_STRENGTH;
  int tab_1=0,tab_2=0,tabp1=0,tabp2=0;
  double cell_back=0,cell_back_r=0,coef_1=0,coef_2=0;
  double hv=0;
  hv = (VOLTAGE_THRESHOLD_IF-VOLTAGE_RESET)/(double)r0->nv;
  cell_back = jumpval/hv;
  tab_1 = (int)floor(-cell_back); tab_2 = tab_1+1;
  cell_back_r = jumpval - hv*floor(jumpval/hv);
  coef_1 = cell_back_r/hv;
  coef_2 = 1.0 - coef_1;
  if (verbose>1){ printf(" %% jumpval %f, hv %f, cell_back %f, cell_back_r %f, tab_1 %d tab_2 %d coef_1 %f coef_2 %f\n",jumpval,hv,cell_back,cell_back_r,tab_1,tab_2,coef_1,coef_2);}
  r0->firingrate=0;
  for (nv=0;nv<r0->nv;nv++){ 
    tabp1 = periodize(nv+tab_1,0,r0->nv); tabp2 = periodize(nv+tab_2,0,r0->nv);
    r0->coef_jumpin[nv] = r0->rho[tabp1]*coef_1 + r0->rho[tabp2]*coef_2;
    r0->coef_jumpout[nv] = r0->rho[nv];
    if (tabp1>nv){ r0->firingrate += r0->rho[tabp1]*coef_1;}
    if (tabp2>nv){ r0->firingrate += r0->rho[tabp2]*coef_2;}}
  r0->firingrate *= AMPA_RATE*1024.0/1.0;
  if (verbose){ raprintf(r0->coef_jumpin,"double",1,r0->nv,"jumpin:"); raprintf(r0->coef_jumpout,"double",1,r0->nv,"jumpout:");}
}

void rho0_diffuse(struct rho0 *r0)
{
  /* hack job diffusion to stabilize computation */
  double *ra=NULL;
  int nv=0;
  double hv=0;
  double diffusion_coefficient=0;
  hv = (VOLTAGE_THRESHOLD_IF-VOLTAGE_RESET)/(double)r0->nv;
  ra = (double *) tcalloc(r0->nv,sizeof(double));
  for (nv=0;nv<r0->nv;nv++){ ra[nv] = 0.5*(nv-1 >= 0 ? r0->rho[(nv-1)] : r0->rho[(nv+0)]) + 0.5*(nv+1 < r0->nv ? r0->rho[(nv+1)] : r0->rho[(nv+0)]);}
  diffusion_coefficient = GLOBAL_DT;
  for (nv=0;nv<r0->nv;nv++){ r0->rho[nv] = (1-diffusion_coefficient)*r0->rho[nv] + diffusion_coefficient*ra[nv];}
  tfree(ra);ra=NULL;
}

void rho0update(struct rho0 *r0,double t,double dt)
{
  int verbose=0;
  int nv=0;
  double jumprate = AMPA_RATE;
  struct matrix *tmp1=mmake();
/*   for (nv=0;nv<r0->nv;nv++){  */
/*     r0->rho[nv] += (1-dt*jumprate)*((dt)*(r0->fluxv[nv] - r0->fluxv[nv+1])) + (dt*jumprate)*(r0->coef_jumpin[nv] - r0->coef_jumpout[nv]);} */
/*   //rho0_diffuse(r0); */
/*   rho0_getmidpoints(r0,0); */
/*   rho0_getjump(r0); */

  minit0(tmp1,r0->nv,1);
  raplugin(tmp1->mtrx,r0->nv,1,r0->rho,r0->nv,1,0,0);
  /* this line throws dynamics to the wind */
  if (GLOBAL_time<GLOBAL_DT*16){ mtimesm(r0->L,r0->L,r0->L);}
  mtimesm(r0->L,tmp1,tmp1);
  for (nv=0;nv<r0->nv;nv++){ r0->rho[nv] = gentry(tmp1,nv,0);}
  mtfree(tmp1);tfree(tmp1);
  rho0_getrhotilde(r0);
  stats("double",r0->rho,r0->nv,NULL,NULL,&(r0->sum),NULL); r0->sum *= r0->nv;
  if (verbose){ raprintf(r0->rho,"double",1,r0->nv,"rho: ");}
}

void rho0_getrhotilde(struct rho0 *r0)
{
  int nv=0,tab1=0,tab2=0;
  struct matrix *tmp1=mmake();
  for (nv=0;nv<r0->nv;nv++){ 
    tab1 = periodize(nv-1,0,r0->nv); tab2 = periodize(nv,0,r0->nv); 
    r0->Lprimerho[nv] = r0->rho[tab1]-r0->rho[tab2];}
  tmp1=mmake();
  minit0(tmp1,r0->nv,1);
  raplugin(tmp1->mtrx,r0->nv,1,r0->rho_tilde,r0->nv,1,0,0);
  mtimesm(r0->L_base,tmp1,tmp1);
  for (nv=0;nv<r0->nv;nv++){ r0->rho_tilde[nv] = gentry(tmp1,nv,0) + GLOBAL_DT*r0->Lprimerho[nv];}
  mtfree(tmp1);tfree(tmp1);
}

/* Here are the rho1 functions */

struct rho1 *rho1make(int rho1_nv,int rho1_ng)
{
  int verbose=0;
  int nv=0,ng=0;
  int nv1=0,ng1=0,ns1=0,nv2=0,ng2=0,ns2=0;
  double sum=0;
  struct rho1 *r1=NULL;
  struct matrix *tmp1=NULL;
  r1 = (struct rho1 *) tcalloc(1,sizeof(struct rho1));
  r1->nv=rho1_nv;
  r1->ng=rho1_ng;
  r1->ns=r1->nv*r1->ng;
  r1->midv = (double *) tcalloc(r1->nv+1,sizeof(double));
  r1->midg = (double *) tcalloc(r1->ng+1,sizeof(double));
  r1->midv_rho = (double *) tcalloc((r1->nv+1)*r1->ng,sizeof(double));
  r1->midg_rho = (double *) tcalloc(r1->nv*(r1->ng+1),sizeof(double));
  r1->fluxv = (double *) tcalloc((r1->nv+1)*r1->ng,sizeof(double));
  r1->fluxg = (double *) tcalloc(r1->nv*(r1->ng+1),sizeof(double));
  r1->coef_jumpin = (double *) tcalloc(r1->nv*r1->ng,sizeof(double));
  r1->coef_jumpout = (double *) tcalloc(r1->nv*r1->ng,sizeof(double));
  r1->rho = (double *) tcalloc(r1->nv*r1->ng,sizeof(double));
  rho1_getmidpoints(r1,1);
  r1->Ginv = mmake();
  minit0(r1->Ginv,r1->ns,r1->ns);
  for (nv1=0;nv1<r1->nv;nv1++){ for (ng1=0;ng1<r1->ng;ng1++){ 
    ns1 = nv1+ng1*r1->nv;
    for (nv2=0;nv2<r1->nv;nv2++){ for (ng2=0;ng2<r1->ng;ng2++){ 
      ns2 = nv2+ng2*r1->nv; if (ns2==ns1){ r1->rho[ns2]=1;} else /* if (ns2!=ns1) */{ r1->rho[ns2]=0;}}}
/*     printf("ns1 %d\n",ns1); for (nv2=0;nv2<r1->nv;nv2++){ for (ng2=0;ng2<r1->ng;ng2++){ ns2 = nv2+ng2*r1->nv; printf("%s",r1->rho[ns2]?"X":"O");} printf("\n");} printf("\n"); */
    rho1_getmidpoints(r1,0);
    rho1_getjump(r1);
/*     printf("ns1 %d\n",ns1); for (nv2=0;nv2<r1->nv+1;nv2++){ for (ng2=0;ng2<r1->ng;ng2++){ ns2 = nv2+ng2*(r1->nv+1); printf("%0.6f, ",r1->fluxv[ns2]);} printf("\n");}  */
/*     for (nv2=0;nv2<r1->nv;nv2++){ for (ng2=0;ng2<r1->ng+1;ng2++){ ns2 = nv2+ng2*r1->nv; printf("%0.6f, ",r1->fluxg[ns2]);} printf("\n");} printf("\n"); */
/*     printf("ns1 %d\n",ns1); for (nv2=0;nv2<r1->nv+1;nv2++){ for (ng2=0;ng2<r1->ng;ng2++){ ns2 = nv2+ng2*(r1->nv+1); printf("%s",r1->fluxv[ns2] ? "X":".");} printf("\n");}  */
/*     for (nv2=0;nv2<r1->nv;nv2++){ for (ng2=0;ng2<r1->ng+1;ng2++){ ns2 = nv2+ng2*r1->nv; printf("%s",r1->fluxg[ns2] ? "X":".");} printf("\n");} printf("\n"); */
/*     for (nv2=0;nv2<r1->nv;nv2++){ for (ng2=0;ng2<r1->ng;ng2++){ ns2 = nv2+ng2*r1->nv; printf("%s",r1->coef_jumpin[ns2] ? "X":".");} printf("    "); for (ng2=0;ng2<r1->ng;ng2++){ ns2 = nv2+ng2*r1->nv; printf("%s",r1->coef_jumpout[ns2] ? "X":".");} printf("\n");} printf("\n"); */
    for (nv2=0;nv2<r1->nv;nv2++){ for (ng2=0;ng2<r1->ng;ng2++){ 
      ns2 = nv2+ng2*r1->nv;
      sentry(r1->Ginv,ns2,ns1,(r1->fluxv[nv2+ng2*(r1->nv+1)] - r1->fluxv[(nv2+1)+ng2*(r1->nv+1)]) + (r1->fluxg[nv2+ng2*r1->nv] - r1->fluxg[nv2+(ng2+1)*r1->nv]) + AMPA_RATE*(r1->coef_jumpin[nv2+ng2*r1->nv] - r1->coef_jumpout[nv2+ng2*r1->nv]));}}
/*     printf("ns1 %d\n",ns1); for (nv2=0;nv2<r1->nv;nv2++){ for (ng2=0;ng2<r1->ng;ng2++){ ns2 = nv2+ng2*r1->nv; printf("%s",gentry(r1->Ginv,ns2,ns1) ? "X":".");} printf("\n");} */
}}
  r1->Linv = mmake();
  minit1(r1->Linv,r1->ns,r1->ns);
  tmp1=mmake();mtimesd(r1->Ginv,GLOBAL_DT,tmp1);
  msubtm(r1->Linv,tmp1,r1->Linv);
  mtfree(tmp1);tfree(tmp1);
  r1->L=mmake();
  minv(r1->Linv,r1->L);
/*   mspy(r1->Ginv,0.0000001); */
/*   mspy(r1->Linv,0.0000001); */
/*   mspy(r1->L,0.0000001); */
/*   exit(0); */
  //for (nv=0;nv<r1->nv;nv++){ for (ng=0;ng<r1->ng;ng++){ r1->rho[nv+ng*r1->nv] = exp(-pow(nv-r1->nv/2,2)-pow(ng-r1->ng/2,2));}}
  for (nv=0;nv<r1->nv;nv++){ for (ng=0;ng<r1->ng;ng++){ r1->rho[nv+ng*r1->nv] = exp(-pow(nv-r1->nv/10,2)-pow(ng-r1->ng/10,2));}}
  stats("double",r1->rho,r1->nv*r1->ng,NULL,NULL,&sum,NULL); sum *= r1->nv*r1->ng;
  for (nv=0;nv<r1->nv;nv++){ for (ng=0;ng<r1->ng;ng++){ r1->rho[nv+ng*r1->nv] /= sum;}}
  if (verbose){ raprintf(r1->rho,"double",r1->nv,r1->ng,"rho: ");}
  stats("double",r1->rho,r1->nv*r1->ng,NULL,NULL,&(r1->sum),NULL); r1->sum *= r1->nv*r1->ng;
  rho1_getmidpoints(r1,1);
  rho1_getjump(r1);
  return r1;
}

double rho1_vdot(double v,double g){ return -CONDUCTANCE_LK*(v-VOLTAGE_RESET) - g*(v-VOLTAGE_EX);}
double rho1_gdot(double v,double g){ return -g/TAU_A;}

void rho1_getmidpoints(struct rho1 *r1,int setvg_flag)
{
  /* obtains midpoint values for rho1 in between cells 
     also obtains fluxes --- just for kicks. */
  int verbose=0;
  int nv=0,ng=0;
  int tab_v1=0,tab_v2=0,tab_g1=0,tab_g2=0; double coef_v1=0,coef_v2=0,coef_g1=0,coef_g2=0;
  int tab_ll=0,tab_l=0,tab_r=0,tab_rr=0; double coef_ll=0,coef_l=0,coef_r=0,coef_rr=0;
  double hv=0,hg=0;
  int stencil_number=0;
  hv = (VOLTAGE_THRESHOLD_IF-VOLTAGE_RESET)/(double)r1->nv;
  hg = (CONDUCTANCE_MAX - 0)/(double)r1->ng;
  if (setvg_flag){
    for (nv=0;nv<r1->nv+1;nv++){ r1->midv[nv] = VOLTAGE_RESET + (VOLTAGE_THRESHOLD_IF - VOLTAGE_RESET)*(double)nv/(double)r1->nv;}
    for (ng=0;ng<r1->ng+1;ng++){ r1->midg[ng] = 0 + (CONDUCTANCE_MAX - 0)*(double)ng/(double)r1->ng;}
    if (verbose){ raprintf(r1->midv,"double",1,r1->nv+1,"midv:");raprintf(r1->midg,"double",1,r1->ng+1,"midg:");}}
  switch (stencil_number){
  case 0: /* 2 point symmetric stencil --  linear */
    for (nv=0;nv<r1->nv+1;nv++){ 
      if (nv==0){ tab_v1 = 0; tab_v2 = 1; coef_v1 = 1.5; coef_v2 = -0.5;}
      else if (nv==r1->nv){ tab_v1 = r1->nv-2; tab_v2 = r1->nv-1; coef_v1 = -0.5; coef_v2 = 1.5;}
      else /* if (nv>0 && nv<r1->nv) */{ tab_v1 = nv-1; tab_v2 = nv; coef_v1 = 0.5; coef_v2 = 0.5;}
      for (ng=0;ng<r1->ng;ng++){ 
	r1->midv_rho[nv + ng*(r1->nv+1)] = coef_v1*r1->rho[tab_v1 + ng*r1->nv] + coef_v2*r1->rho[tab_v2 + ng*r1->nv];}}  
    if (verbose){ raprintf(r1->midv_rho,"double",r1->nv+1,r1->ng,"midv_rho: ");}
    for (ng=0;ng<r1->ng+1;ng++){ 
      if (ng==0){ tab_g1 = 0; tab_g2 = 1; coef_g1 = 1.5; coef_g2 = -0.5;}
      else if (ng==r1->ng){ tab_g1 = r1->ng-2; tab_g2 = r1->ng-1; coef_g1 = -0.5; coef_g2 = 1.5;}
      else /* if (ng>0 && ng<r1->ng) */{ tab_g1 = ng-1; tab_g2 = ng; coef_g1 = 0.5; coef_g2 = 0.5;}
      for (nv=0;nv<r1->nv;nv++){
	r1->midg_rho[nv + ng*r1->nv] = coef_g1*r1->rho[nv + tab_g1*r1->nv] + coef_g2*r1->rho[nv + tab_g2*r1->nv];}}
    if (verbose){ raprintf(r1->midg_rho,"double",r1->nv,r1->ng+1,"midg_rho: ");}
    break;
  case 1: /* 4 point symmetric stencil -- linear least squares*/
    for (nv=0;nv<r1->nv+1;nv++){ 
      if (nv==0){ tab_ll=0; tab_l=1; tab_r=2; tab_rr=3; coef_ll=17.0/20.0; coef_l=9.0/20.0; coef_r=1.0/20.0; coef_rr=-7.0/20.0;}
      else if (nv==1){ tab_ll=0; tab_l=1; tab_r=2; tab_rr=3; coef_ll=11.0/20.0; coef_l=7.0/20.0; coef_r=3.0/20.0; coef_rr=-1.0/20.0;}
      else if (nv>1 && nv<r1->nv-1){ tab_ll=nv-2; tab_l=nv-1; tab_r=nv; tab_rr=nv+1; coef_ll=5.0/20.0; coef_l=5.0/20.0; coef_r=5.0/20.0; coef_rr=5.0/20.0;}
      else if (nv==r1->nv-1){ tab_ll=r1->nv-4; tab_l=r1->nv-3; tab_r=r1->nv-2; tab_rr=r1->nv-1; coef_ll=-1.0/20.0; coef_l=3.0/20.0; coef_r=7.0/20.0; coef_rr=11.0/20.0;}
      else if (nv==r1->nv){ tab_ll=r1->nv-4; tab_l=r1->nv-3; tab_r=r1->nv-2; tab_rr=r1->nv-1; coef_ll=-7.0/20.0; coef_l=1.0/20.0; coef_r=9.0/20.0; coef_rr=17.0/20.0;}
      for (ng=0;ng<r1->ng;ng++){ 
	r1->midv_rho[nv + ng*(r1->nv+1)] = coef_ll*r1->rho[tab_ll + ng*r1->nv] + coef_l*r1->rho[tab_l + ng*r1->nv] + coef_r*r1->rho[tab_r + ng*r1->nv] + coef_rr*r1->rho[tab_rr + ng*r1->nv];}}
    if (verbose){ raprintf(r1->midv_rho,"double",r1->nv+1,r1->ng,"midv_rho: ");}
    for (ng=0;ng<r1->ng+1;ng++){ 
      if (ng==0){ tab_ll=0; tab_l=1; tab_r=2; tab_rr=3; coef_ll=17.0/20.0; coef_l=9.0/20.0; coef_r=1.0/20.0; coef_rr=-7.0/20.0;}
      else if (ng==1){ tab_ll=0; tab_l=1; tab_r=2; tab_rr=3; coef_ll=11.0/20.0; coef_l=7.0/20.0; coef_r=3.0/20.0; coef_rr=-1.0/20.0;}
      else if (ng>1 && ng<r1->ng-1){ tab_ll=ng-2; tab_l=ng-1; tab_r=ng; tab_rr=ng+1; coef_ll=5.0/20.0; coef_l=5.0/20.0; coef_r=5.0/20.0; coef_rr=5.0/20.0;}
      else if (ng==r1->ng-1){ tab_ll=r1->ng-4; tab_l=r1->ng-3; tab_r=r1->ng-2; tab_rr=r1->ng-1; coef_ll=-1.0/20.0; coef_l=3.0/20.0; coef_r=7.0/20.0; coef_rr=11.0/20.0;}
      else if (ng==r1->ng){ tab_ll=r1->ng-4; tab_l=r1->ng-3; tab_r=r1->ng-2; tab_rr=r1->ng-1; coef_ll=-7.0/20.0; coef_l=1.0/20.0; coef_r=9.0/20.0; coef_rr=17.0/20.0;}
      for (nv=0;nv<r1->nv;nv++){
	r1->midg_rho[nv + ng*r1->nv] = coef_ll*r1->rho[nv + tab_ll*r1->nv] + coef_l*r1->rho[nv + tab_l*r1->nv] + coef_r*r1->rho[nv + tab_r*r1->nv] +  coef_rr*r1->rho[nv + tab_rr*r1->nv];}}
    if (verbose){ raprintf(r1->midg_rho,"double",r1->nv,r1->ng+1,"midg_rho: ");}
    break;
  case 2: /* 4 point symmetric stencil -- cubic */
    for (nv=0;nv<r1->nv+1;nv++){ 
      if (nv==0){ tab_ll=0; tab_l=1; tab_r=2; tab_rr=3; coef_ll=35.0/16.0; coef_l=-35.0/16.0; coef_r=21.0/16.0; coef_rr=-5.0/16.0;}
      else if (nv==1){ tab_ll=0; tab_l=1; tab_r=2; tab_rr=3; coef_ll=5.0/16.0; coef_l=15.0/16.0; coef_r=-5.0/16.0; coef_rr=1.0/16.0;}
      else if (nv>1 && nv<r1->nv-1){ tab_ll=nv-2; tab_l=nv-1; tab_r=nv; tab_rr=nv+1; coef_ll=-1.0/16.0; coef_l=9.0/16.0; coef_r=9.0/16.0; coef_rr=-1.0/16.0;}
      else if (nv==r1->nv-1){ tab_ll=r1->nv-4; tab_l=r1->nv-3; tab_r=r1->nv-2; tab_rr=r1->nv-1; coef_ll=1.0/16.0; coef_l=-5.0/16.0; coef_r=15.0/16.0; coef_rr=5.0/16.0;}
      else if (nv==r1->nv){ tab_ll=r1->nv-4; tab_l=r1->nv-3; tab_r=r1->nv-2; tab_rr=r1->nv-1; coef_ll=-5.0/16.0; coef_l=21.0/16.0; coef_r=-35.0/16.0; coef_rr=35.0/16.0;}
      for (ng=0;ng<r1->ng;ng++){ 
	r1->midv_rho[nv + ng*(r1->nv+1)] = coef_ll*r1->rho[tab_ll + ng*r1->nv] + coef_l*r1->rho[tab_l + ng*r1->nv] + coef_r*r1->rho[tab_r + ng*r1->nv] + coef_rr*r1->rho[tab_rr + ng*r1->nv];}}
    if (verbose){ raprintf(r1->midv_rho,"double",r1->nv+1,r1->ng,"midv_rho: ");}
    for (ng=0;ng<r1->ng+1;ng++){ 
      if (ng==0){ tab_ll=0; tab_l=1; tab_r=2; tab_rr=3; coef_ll=35.0/16.0; coef_l=-35.0/16.0; coef_r=21.0/16.0; coef_rr=-5.0/16.0;}
      else if (ng==1){ tab_ll=0; tab_l=1; tab_r=2; tab_rr=3; coef_ll=5.0/16.0; coef_l=15.0/16.0; coef_r=-5.0/16.0; coef_rr=1.0/16.0;}
      else if (ng>1 && ng<r1->ng-1){ tab_ll=ng-2; tab_l=ng-1; tab_r=ng; tab_rr=ng+1; coef_ll=-1.0/16.0; coef_l=9.0/16.0; coef_r=9.0/16.0; coef_rr=-1.0/16.0;}
      else if (ng==r1->ng-1){ tab_ll=r1->ng-4; tab_l=r1->ng-3; tab_r=r1->ng-2; tab_rr=r1->ng-1; coef_ll=1.0/16.0; coef_l=-5.0/16.0; coef_r=15.0/16.0; coef_rr=5.0/16.0;}
      else if (ng==r1->ng){ tab_ll=r1->ng-4; tab_l=r1->ng-3; tab_r=r1->ng-2; tab_rr=r1->ng-1; coef_ll=-5.0/16.0; coef_l=21.0/16.0; coef_r=-35.0/16.0; coef_rr=35.0/16.0;}
      for (nv=0;nv<r1->nv;nv++){
	r1->midg_rho[nv + ng*r1->nv] = coef_ll*r1->rho[nv + tab_ll*r1->nv] + coef_l*r1->rho[nv + tab_l*r1->nv] + coef_r*r1->rho[nv + tab_r*r1->nv] +  coef_rr*r1->rho[nv + tab_rr*r1->nv];}}
    if (verbose){ raprintf(r1->midg_rho,"double",r1->nv,r1->ng+1,"midg_rho: ");}
    break;
  default: printf(" %% warning, incorrect stencil type in rho1_getmidpoints\n");}
  for (nv=0;nv<r1->nv+1;nv++){ for (ng=0;ng<r1->ng;ng++){ 
    r1->fluxv[nv + ng*(r1->nv+1)] = rho1_vdot(r1->midv[nv],(r1->midg[ng]+r1->midg[ng+1])/2)*r1->midv_rho[nv+ng*(r1->nv+1)]/hv;}}
  for (ng=0;ng<r1->ng;ng++){ 
    if (r1->fluxv[r1->nv+ng*(r1->nv+1)]>0){ r1->fluxv[0 + ng*(r1->nv+1)] = r1->fluxv[r1->nv + ng*(r1->nv+1)];}
    else /* if (r1->fluxv[r1->nv+ng*(r1->nv+1)]<=0) */{ r1->fluxv[0 + ng*(r1->nv+1)] = 0; r1->fluxv[r1->nv + ng*(r1->nv+1)] = 0;}}  
  if (verbose){ raprintf(r1->fluxv,"double",r1->nv+1,r1->ng,"fluxv: ");}  
  for (nv=0;nv<r1->nv;nv++){ for (ng=0;ng<r1->ng+1;ng++){ 
    r1->fluxg[nv + ng*r1->nv] = rho1_gdot((r1->midv[nv]+r1->midv[nv+1])/1,r1->midg[ng])*r1->midg_rho[nv+ng*r1->nv]/hg;}}
  for (nv=0;nv<r1->nv;nv++){ r1->fluxg[nv + 0*r1->nv] = 0; r1->fluxg[nv + r1->ng*r1->nv] = 0;}
  if (verbose){ raprintf(r1->fluxg,"double",r1->nv,r1->ng+1,"fluxg: ");}  
}

void rho1_getjump(struct rho1 *r1)
{
  int verbose=0;
  int nv=0,ng=0;
  double jumpval = AMPA_STRENGTH;
  int tab_1=0,tab_2=0;
  double cell_back=0,cell_back_r=0,coef_1=0,coef_2=0;
  double hv=0,hg=0;
  hv = (VOLTAGE_THRESHOLD_IF-VOLTAGE_RESET)/(double)r1->nv;
  hg = (CONDUCTANCE_MAX - 0)/(double)r1->ng;
  cell_back = jumpval/hg;
  tab_1 = (int)floor(-cell_back); tab_2 = tab_1+1;
  cell_back_r = jumpval - hg*floor(jumpval/hg);
  coef_1 = cell_back_r/hg;
  coef_2 = 1.0 - coef_1;
  if (verbose>1){ printf(" %% jumpval %f, hg %f, cell_back %f, cell_back_r %f, tab_1 %d tab_2 %d coef_1 %f coef_2 %f\n",jumpval,hg,cell_back,cell_back_r,tab_1,tab_2,coef_1,coef_2);}
  for (nv=0;nv<r1->nv;nv++){ for (ng=0;ng<r1->ng-1;ng++){
    r1->coef_jumpin[nv + ng*r1->nv] = (ng+tab_1 >= 0)*r1->rho[nv+(ng+tab_1)*r1->nv]*coef_1 + (ng+tab_2 >= 0)*r1->rho[nv+(ng+tab_2)*r1->nv]*coef_2;
    r1->coef_jumpout[nv + ng*r1->nv] = r1->rho[nv+ng*r1->nv];}}
  for (nv=0;nv<r1->nv;nv++){ 
    r1->coef_jumpin[nv + (r1->ng-1)*r1->nv] = (r1->ng-1+tab_1 >= 0)*r1->rho[nv + (r1->ng-1+tab_1)*r1->nv]*coef_1;
    for (ng=r1->ng-1+tab_1+1;ng<r1->ng-1;ng++){ r1->coef_jumpin[nv + (r1->ng-1)*r1->nv] += r1->rho[nv + ng*r1->nv]*1.0;}
    r1->coef_jumpout[nv+(r1->ng-1)*r1->nv] = 0;}
  if (verbose){ 
    raprintf(r1->coef_jumpin,"double",r1->nv,r1->ng,"jumpin:");
    raprintf(r1->coef_jumpout,"double",r1->nv,r1->ng,"jumpout:");}
}

void rho1_diffuse(struct rho1 *r1)
{
  /* hack job diffusion to stabilize computation */
  double *ra=NULL;
  int nv=0,ng=0;
  double hv=0,hg=0;
  double diffusion_coefficient=0;
  hv = (VOLTAGE_THRESHOLD_IF-VOLTAGE_RESET)/(double)r1->nv;
  hg = (CONDUCTANCE_MAX - 0)/(double)r1->ng;
  ra = (double *) tcalloc(r1->nv*r1->ng,sizeof(double));
  for (nv=0;nv<r1->nv;nv++){ for (ng=0;ng<r1->ng;ng++){ ra[nv+ng*r1->nv] = 0.25*(nv-1 >= 0 ? r1->rho[(nv-1) + (ng+0)*r1->nv] : r1->rho[(nv+0) + (ng+0)*r1->nv]) + 0.25*(nv+1 < r1->nv ? r1->rho[(nv+1) + (ng+0)*r1->nv] : r1->rho[(nv+0) + (ng+0)*r1->nv]) + 0.25*(ng-1 >= 0 ? r1->rho[(nv+0) + (ng-1)*r1->nv] : r1->rho[(nv+0) + (ng+0)*r1->nv]) + 0.25*(ng+1 < r1->ng ? r1->rho[(nv+0) + (ng+1)*r1->nv] : r1->rho[(nv+0) + (ng+0)*r1->nv]);}}
  diffusion_coefficient = GLOBAL_DT;
  for (nv=0;nv<r1->nv;nv++){ for (ng=0;ng<r1->ng;ng++){ r1->rho[nv+ng*r1->nv] = (1-diffusion_coefficient)*r1->rho[nv+ng*r1->nv] + diffusion_coefficient*ra[nv+ng*r1->nv];}}
  tfree(ra);ra=NULL;
}

void rho1update(struct rho1 *r1,double t,double dt)
{
  int verbose=0;
  int nv=0,ng=0;
  double jumprate = AMPA_RATE;
/*   struct matrix *tmp1=NULL; */
  for (nv=0;nv<r1->nv;nv++){ for (ng=0;ng<r1->ng;ng++){
    r1->rho[nv+ng*r1->nv] += (1-dt*jumprate)*((dt)*(r1->fluxv[nv+ng*(r1->nv+1)] - r1->fluxv[(nv+1)+ng*(r1->nv+1)]) + (dt)*(r1->fluxg[nv+ng*r1->nv] - r1->fluxg[nv+(ng+1)*r1->nv])) + (dt*jumprate)*(r1->coef_jumpin[nv+ng*r1->nv] - r1->coef_jumpout[nv+ng*r1->nv]);}}
  //rho1_diffuse(r1);
  rho1_getmidpoints(r1,0);
  rho1_getjump(r1);
/*   tmp1=mmake(); */
/*   minit0(tmp1,r1->nv*r1->ng,1); */
/*   raplugin(tmp1->mtrx,r1->nv,r1->ng,r1->rho,r1->nv,r1->ng,0,0); */
/*   mtimesm(r1->L,tmp1,tmp1); */
/*   for (nv=0;nv<r1->nv;nv++){ for (ng=0;ng<r1->ng;ng++){ r1->rho[nv+ng*r1->nv] = gentry(tmp1,nv+ng*r1->nv,0);}} */
/*   mtfree(tmp1);tfree(tmp1);   */
  stats("double",r1->rho,r1->nv*r1->ng,NULL,NULL,&(r1->sum),NULL); r1->sum *= r1->nv*r1->ng;
  if (verbose){ raprintf(r1->rho,"double",r1->nv,r1->ng,"rho: ");}
}

/* Here are the neuron functions */

struct neuron * neuronmake(int types/* ,double fA,double nuA,double fN,double nuN,double fG,double nuG */)
{
  int nt=0;
  double *temp=NULL;
  struct neuron *n = (struct neuron *) malloc(sizeof(struct neuron));
  n->types = types;
/*   n->fA=fA; n->nuA=nuA; */
/*   n->fN=fN; n->nuN=nuN; */
/*   n->fG=fG; n->nuG=nuG; */
  n->sA = llistmake(); n->sN = llistmake(); n->sG = llistmake();
  temp=(double *)tcalloc(1,sizeof(double)); litemadd(n->sA,temp);
  temp=(double *)tcalloc(1,sizeof(double)); litemadd(n->sN,temp);
  temp=(double *)tcalloc(1,sizeof(double)); litemadd(n->sG,temp);
  n->V = (struct llist **) tcalloc(n->types,sizeof(struct llist *));
  n->spikenext = (struct llist **) tcalloc(n->types,sizeof(struct llist *));
  n->t = (struct llist **) tcalloc(n->types,sizeof(struct llist *));
  for (nt=0;nt<n->types;nt++){ 
    temp=(double *)tcalloc(1,sizeof(double));
    n->V[nt]=llistmake(); litemadd(n->V[nt],temp);
    temp=(double *)tcalloc(1,sizeof(double));  *temp = -1;
    n->spikenext[nt]=llistmake(); litemadd(n->spikenext[nt],temp);
    temp=(double *)tcalloc(1,sizeof(double)); 
    n->t[nt]=llistmake(); litemadd(n->t[nt],temp);}
  n->sra = (double *) tcalloc(n->types,sizeof(double));
  n->Vra = (double *) tcalloc(n->types,sizeof(double));
  n->tra = (double *) tcalloc(n->types,sizeof(double));
  n->nbins = 32;
  n->hist = (double *) tcalloc(n->nbins*n->nbins*n->types,sizeof(double));
  return n;
}

void neurontfree(struct neuron *n)
{
  int nt=0;
  llisttfree2(n->sA);llisttfree2(n->sN);llisttfree2(n->sG);
  for (nt=0;nt<n->types;nt++){
    llisttfree2(n->V[nt]); llisttfree2(n->spikenext[nt]); llisttfree2(n->t[nt]);}
  tfree(n->V);tfree(n->spikenext);tfree(n->t);
  tfree(n->sra);tfree(n->Vra);tfree(n->tra);
  tfree(n->hist);
  tfree(n); n=NULL;
}

double quadrootfinder(double VI, double k1, double VF, double VT, double DT)
{
  int i=0,imax=10;
  double a=0,b=0,c=0,e=0.00000000001;
  double r=0,dr=0;
  if (VI<= VT && VF >= VT){
    c = VI;
    b = k1;
    a = (VF - c - DT*b)/pow(DT,2);
    r = DT*(VT-VI)/(VF-VI);
    do{
      dr = -((a*r+b)*r+c-VT)/(2.0*a*r+b);
      r = r+dr;
      i++;}
    while (dr/r > e && i<imax);
    if (i>=imax || r<=0 || r>=DT){ 
      r = DT*(VT-VI)/(VF-VI);}}
  else{ r = 2*DT;}
  return r;
}

void slavecalc(double sa,double sn,double sg,double *gs,double *vs,double *vsprime)
{
  /* given conductances sa,sn,sg that decay exponentially, we compute the 
     slaving conductance gs
     slaving voltage vs
     slaving voltage prime vsprime */
  double a=0,b=0,c=0,d=0;
  a = CONDUCTANCE_LK + sa + sn + sg;
  b = -sa/TAU_A - sn/TAU_N - sg/TAU_G;
  c = CONDUCTANCE_LK*VOLTAGE_LEAK + (sa+sn)*VOLTAGE_EX + (sg)*VOLTAGE_IN;
  d = -(sa/TAU_A+sn/TAU_N)*VOLTAGE_EX - (sg/TAU_G)*VOLTAGE_IN;
  if (gs!=NULL){ *gs = a;}
  if (vs!=NULL){ *vs = c/a;}
  if (vsprime!=NULL){ *vsprime = (a*d - c*b)/(a*a);}
}

int ifevolve(double dtmax,double V,double spikenext,double sA,double sN,double sG,double t,double DT,double *V2,double *spikenext2)
{
  int spike_flag=0;
  double dt=DT,t2=t;
  double d=0,sA2=0,sN2=0,sG2=0,k1=0,eA=0,eN=0,eG=0,if1=0,if2=0,if3=0,vs1=0,vs2=0,vs3=0,vsp1=0,vsp2=0,vsp3=0,g=0,Vstart=0,Vfinish=0;
  double spiketime_guess=0;
  if (dt>dtmax){ dt=DT;while(dt>dtmax){dt/=2;}}
  while (t2<t+DT){
    if (spikenext >= t2+dt){
      d = dt;//t2+dt-t2;
      sA *= exp(-d/TAU_A); sN *= exp(-d/TAU_N); sG *= exp(-d/TAU_G);
      Vstart = VOLTAGE_RESET;
      Vfinish = VOLTAGE_RESET;}
    else{
      Vstart = V;
      if (spikenext >= t2 && spikenext < t2+dt){
	d = spikenext - t2;
	sA *= exp(-d/TAU_A); sN *= exp(-d/TAU_N); sG *= exp(-d/TAU_G);
	Vstart = VOLTAGE_RESET;
	Vfinish = VOLTAGE_RESET;}
      sA2 = 0; sN2 = 0; sG2 = 0;
      d = t2+dt - maximum(spikenext,t2);
      k1 = (CONDUCTANCE_LK+sA+sN+sG)*V + (VOLTAGE_LEAK*CONDUCTANCE_LK + VOLTAGE_EX*(sA+sN) + VOLTAGE_IN*(sG));
      eA = exp(-d/2/TAU_A); eN = exp(-d/2/TAU_N); eG = exp(-d/2/TAU_G);
      if1 = 0*CONDUCTANCE_LK+sA2+sN2+sG2;
      slavecalc(sA,sN,sG,NULL,&vs1,&vsp1);
      sA2 = TAU_A*sA*(1 - eA);
      sN2 = TAU_N*sN*(1 - eN);
      sG2 = TAU_G*sG*(1 - eG);
      sA *= eA;
      sN *= eN;
      sG *= eG;
      if2 = d/2*CONDUCTANCE_LK+sA2+sN2+sG2;
      slavecalc(sA,sN,sG,NULL,&vs2,&vsp2);
      sA2 += TAU_A*sA*(1 - eA);
      sN2 += TAU_N*sN*(1 - eN);
      sG2 += TAU_G*sG*(1 - eG);
      sA *= eA;
      sN *= eN;
      sG *= eG;
      if3 = d*CONDUCTANCE_LK+sA2+sN2+sG2;
      slavecalc(sA,sN,sG,NULL,&vs3,&vsp3);
      g = d/6.0*(exp(if1-if3)*vsp1 + 4*exp(if2-if3)*vsp2 + exp(if3-if3)*vsp3);
      Vfinish = ((Vstart-vs1)*exp(-if3) - g) + vs3;
      if (Vstart < VOLTAGE_THRESHOLD_IF && Vfinish >= VOLTAGE_THRESHOLD_IF){
	spiketime_guess = maximum(t2,spikenext) + quadrootfinder(Vstart,k1,Vfinish,VOLTAGE_THRESHOLD_IF,d);
	spike_flag = 1;}}
    V = Vfinish;//printf("%f\n",V);
    if (spiketime_guess > t2 && spiketime_guess <= t2+dt){ V=VOLTAGE_RESET; spikenext = spiketime_guess+TAU_REF;}
    t2 += dt;}
  *V2 = V; *spikenext2 = spikenext;
  return spike_flag;
}

void llistbanzai(struct llist *L)
{
  /* given a llist *L of double * items, kill the first element */
  struct litem *l=NULL;
  double *temp=NULL;
  if (L->length>MAXLENGTH){
    l=L->first;
    temp = (double *)l->item;
    L->first=l->child;
    litemtfree(l);
    tfree(temp);
    L->length -= 1;
    llistbanzai(L);}
}

void neuronevolve(struct neuron *n,double DT)
{
  /* evolves neuron */
  int nt=0,spike_flag=0;
  double V=0,spikenext=0,t=0,sA=0,sN=0,sG=0,V2=0,spikenext2=0,*temp=NULL;
  int vbin=0,sAbin=0;
  sA = *((double *) n->sA->last->item);
  sN = *((double *) n->sN->last->item);
  sG = *((double *) n->sG->last->item);
  for (nt=0;nt<n->types;nt++){
    V = *((double *) n->V[nt]->last->item);
    spikenext = *((double *) n->spikenext[nt]->last->item);
    t = *((double *) n->t[nt]->last->item);
    switch (nt){
    case 0: spike_flag = ifevolve(GLOBAL_DT,V,spikenext,sA,sN,sG,t,DT,&V2,&spikenext2); break;
    default: break;}
    n->Vra[nt] += V2*DT;
    n->sra[nt] += spike_flag;
    n->tra[nt] = t+DT;
    temp=(double *)tcalloc(1,sizeof(double));*temp=V2;litemadd(n->V[nt],temp); llistbanzai(n->V[nt]);
    temp=(double *)tcalloc(1,sizeof(double));*temp=spikenext2;litemadd(n->spikenext[nt],temp); llistbanzai(n->spikenext[nt]);
    temp=(double *)tcalloc(1,sizeof(double));*temp=t+DT;litemadd(n->t[nt],temp); llistbanzai(n->t[nt]);}
  sA *= exp(-DT/TAU_A); sN *= exp(-DT/TAU_N); sG *= exp(-DT/TAU_G);
  sA += (drand(0)<AMPA_RATE*DT ? AMPA_STRENGTH/TAU_A*(GLOBAL_ARBOR_SIM ? exp(-(pow(drand(1),2)+pow(drand(1),2))/.3)/PI/.3 : 1) : 0);
  sN += (drand(0)<NMDA_RATE*DT ? NMDA_STRENGTH/TAU_N*(GLOBAL_ARBOR_SIM ? exp(-(pow(drand(1),2)+pow(drand(1),2))/.3)/PI/.3 : 1) : 0);
  sG += (drand(0)<GABA_RATE*DT ? GABA_STRENGTH/TAU_G*(GLOBAL_ARBOR_SIM ? exp(-(pow(drand(1),2)+pow(drand(1),2))/.3)/PI/.3 : 1) : 0);
  temp=(double *)tcalloc(1,sizeof(double));*temp=sA;litemadd(n->sA,temp);
  temp=(double *)tcalloc(1,sizeof(double));*temp=sN;litemadd(n->sN,temp);
  temp=(double *)tcalloc(1,sizeof(double));*temp=sG;litemadd(n->sG,temp);
  llistbanzai(n->sA);
  llistbanzai(n->sN);
  llistbanzai(n->sG);
  for (nt=0;nt<n->types;nt++){
    vbin=crop((int)floor((double)n->nbins*(*(double *)n->V[nt]->last->item-VOLTAGE_RESET)/(VOLTAGE_THRESHOLD_IF-VOLTAGE_RESET)),0,n->nbins-1);
    sAbin=crop((int)floor((double)n->nbins*(*(double *)n->sA->last->item-0)/(2*GHAT-0)),0,n->nbins-1);
    n->hist[vbin + sAbin*n->nbins + nt*n->nbins*n->nbins] += DT;}
}

void setglobals()
{
#ifndef NOTGRAPH
  GLOBAL_WINDOW_WIDTH = 640;
  GLOBAL_WINDOW_HEIGHT = 480;
  xdepth = -.5f;
  ydepth = -1.75f;
  zdepth = -4.4f;
#endif /* NOTGRAPH */
  GHAT = CONDUCTANCE_LK*(VOLTAGE_THRESHOLD_IF-VOLTAGE_RESET)/(VOLTAGE_EX-VOLTAGE_THRESHOLD_IF);
  GLOBAL_NEURON = neuronmake(GLOBAL_NTYPES);
  switch (GLOBAL_RHO_VERSION){
  case 0: GLOBAL_rho0 = rho0make(GLOBAL_rho0_nv); break;
  case 1: GLOBAL_rho1 = rho1make(GLOBAL_rho1_nv,GLOBAL_rho1_ng); break;
  default: break;}
}

void computestep()
{
  switch(GLOBAL_RHO_VERSION){
  case 0: rho0update(GLOBAL_rho0,GLOBAL_time,GLOBAL_DT); break;
  case 1: rho1update(GLOBAL_rho1,GLOBAL_time,GLOBAL_DT); break;
  default: break;}
  GLOBAL_time += GLOBAL_DT;
}

#ifndef NOTGRAPH

void ftexto(float x, float y, float z, char *text)
{
  /* thanks to alex */
  char *p;
  glRasterPos3f(x,y,z);
  for (p = text; *p; p++){ glutBitmapCharacter(GLUT_BITMAP_8_BY_13, *p);}
}

GLvoid InitGL(GLsizei Width, GLsizei Height)
{
  /* A general OpenGL initialization function.  Sets all of the initial parameters. */
  glClearColor(0.0f, 0.0f, 0.0f, 0.0f);	// This Will Clear The Background Color To Black
  glClearDepth(1.0);				// Enables Clearing Of The Depth Buffer
  glDepthFunc(GL_LESS);			// The Type Of Depth Test To Do
  glEnable(GL_DEPTH_TEST);			// Enables Depth Testing
  glShadeModel(GL_SMOOTH);			// Enables Smooth Color Shading
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();				// Reset The Projection Matrix
  gluPerspective(45.0f,(GLfloat)Width/(GLfloat)Height,0.1f,100.0f);	// Calculate The Aspect Ratio Of The Window
  glMatrixMode(GL_MODELVIEW);
}

GLvoid ReSizeGLScene(GLsizei Width, GLsizei Height)
{
  /* The function called when our window is resized */
  if (Height==0){ Height=1;} /* don't divide by zero */
  glViewport(0, 0, Width, Height); /* reset viewport & perspective */
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(45.0f,(GLfloat)Width/(GLfloat)Height,0.1f,100.0f);
  glMatrixMode(GL_MODELVIEW);
}

GLvoid glbox(double x,double y,double r,double rcolor,double gcolor,double bcolor)
{
  glBegin(GL_QUADS);
  glColor3f(rcolor,gcolor,bcolor);
  glVertex3f(x+r/2,y+r/2,0);
  glVertex3f(x+r/2,y-r/2,0);
  glVertex3f(x-r/2,y-r/2,0);
  glVertex3f(x-r/2,y+r/2,0);
  glEnd();
}

GLvoid glrect(double x,double y,double rx,double ry,double rcolor,double gcolor,double bcolor)
{
  glBegin(GL_QUADS);
  glColor3f(rcolor,gcolor,bcolor);
  glVertex3f(x+rx/2,y+ry/2,0);
  glVertex3f(x+rx/2,y-ry/2,0);
  glVertex3f(x-rx/2,y-ry/2,0);
  glVertex3f(x-rx/2,y+ry/2,0);
  glEnd();
}

GLvoid glarrow(int nheads,double x1,double y1,double x2,double y2,double rmin,double rcolor,double gcolor,double bcolor)
{
  double theta = atan2(y2-y1,x2-x1),theta1=theta+PI-PI/8,theta2=theta+PI+PI/8;
  int nh=0;
  double x3=(x1+x2)/2.0,y3=(y1+y2)/2.0;
  double x4=0,y4=0;
  glBegin(GL_LINES);
  glColor3f(rcolor,gcolor,bcolor);
  glVertex3f(x1,y1,0); glVertex3f(x2,y2,0);
  for (nh=0;nh<nheads;nh++){
    x4 = x3 + (double)nh*rmin/5.0*cos(theta);
    y4 = y3 + (double)nh*rmin/5.0*sin(theta);
    glVertex3f(x4,y4,0);
    glVertex3f(x4+rmin*cos(theta1),y4+rmin*sin(theta1),0);
    glVertex3f(x4,y4,0);
    glVertex3f(x4+rmin*cos(theta2),y4+rmin*sin(theta2),0);}
  glEnd();
}

GLvoid Drawmenu(double side,double xoffset,double yoffset)
{
  int menupos=0;
  char text[32];
  double tcolor=0;
  struct neuron *n=GLOBAL_NEURON;
  int nt=0;
  menupos=0; tcolor = .25+.75*(FIDDLE_PARAMETER==menupos);glColor3f(tcolor,tcolor,tcolor);
  sprintf(text,"STEPS_PER_DRAW=%d",STEPS_PER_DRAW);
  ftexto(0*side+xoffset,menupos*side+yoffset,0,text);
  menupos=1; tcolor = .25+.75*(FIDDLE_PARAMETER==menupos);glColor3f(tcolor,tcolor,tcolor);
  sprintf(text,"GLOBAL_DT=%0.3f",GLOBAL_DT);
  ftexto(0*side+xoffset,menupos*side+yoffset,0,text);
  menupos=2; tcolor = .25+.75*(FIDDLE_PARAMETER==menupos);glColor3f(tcolor,tcolor,tcolor);
  sprintf(text,"TAU_A=%0.3f",TAU_A);
  ftexto(0*side+xoffset,menupos*side+yoffset,0,text);
  menupos=3; tcolor = .25+.75*(FIDDLE_PARAMETER==menupos);glColor3f(tcolor,tcolor,tcolor);
  sprintf(text,"AMPA_STRENGTH=%0.3f",AMPA_STRENGTH);
  ftexto(0*side+xoffset,menupos*side+yoffset,0,text);
  menupos=4; tcolor = .25+.75*(FIDDLE_PARAMETER==menupos);glColor3f(tcolor,tcolor,tcolor);
  sprintf(text,"AMPA_RATE=%0.3f",AMPA_RATE);
  ftexto(0*side+xoffset,menupos*side+yoffset,0,text);
  menupos=5; tcolor = .25+.75*(FIDDLE_PARAMETER==menupos);glColor3f(tcolor,tcolor,tcolor);
  sprintf(text,"DRAW_FLAG %d ",DRAW_FLAG);
  ftexto(0*side+xoffset,menupos*side+yoffset,0,text);
}

GLvoid Drawrho0(double side,double xoffset,double yoffset)
{
  int verbose=0;
  char text[128];
  struct rho0 *r0=GLOBAL_rho0;
  int nv=0;
  double xord=0,yord=0;
  double hv=0;
  hv = (VOLTAGE_THRESHOLD_IF-VOLTAGE_RESET)/(double)r0->nv;
  for (nv=0;nv<r0->nv;nv++){
    glBegin(GL_LINE_STRIP);glColor3f(1,1,1);
    xord = (double)(nv+0.0)/(double)(r0->nv);
    yord = 0*r0->rho[nv]/hv;
    glVertex3f(xord*side+xoffset,yord*side+yoffset,0);
    xord = (double)(nv+0.0)/(double)(r0->nv);
    yord = 1*r0->rho[nv]/hv;
    glVertex3f(xord*side+xoffset,yord*side+yoffset,0);
    xord = (double)(nv+1.0)/(double)(r0->nv);
    yord = 1*r0->rho[nv]/hv;
    glVertex3f(xord*side+xoffset,yord*side+yoffset,0);
    xord = (double)(nv+1.0)/(double)(r0->nv);
    yord = 0*r0->rho[nv]/hv;
    glVertex3f(xord*side+xoffset,yord*side+yoffset,0);
    xord = (double)(nv+0.0)/(double)(r0->nv);
    yord = 0*r0->rho[nv]/hv;
    glVertex3f(xord*side+xoffset,yord*side+yoffset,0);
    glEnd();}
  for (nv=0;nv<r0->nv;nv++){
    glBegin(GL_LINE_STRIP);glColor3f(0.5,0.76,0.6);
    xord = (double)(nv+0.0)/(double)(r0->nv);
    yord = 0*r0->rho_tilde[nv]/hv;
    glVertex3f(xord*side+xoffset+side,yord*side+yoffset,0);
    xord = (double)(nv+0.0)/(double)(r0->nv);
    yord = 1*r0->rho_tilde[nv]/hv;
    glVertex3f(xord*side+xoffset+side,yord*side+yoffset,0);
    xord = (double)(nv+1.0)/(double)(r0->nv);
    yord = 1*r0->rho_tilde[nv]/hv;
    glVertex3f(xord*side+xoffset+side,yord*side+yoffset,0);
    xord = (double)(nv+1.0)/(double)(r0->nv);
    yord = 0*r0->rho_tilde[nv]/hv;
    glVertex3f(xord*side+xoffset+side,yord*side+yoffset,0);
    xord = (double)(nv+0.0)/(double)(r0->nv);
    yord = 0*r0->rho_tilde[nv]/hv;
    glVertex3f(xord*side+xoffset+side,yord*side+yoffset,0);
    glEnd();}
  sprintf(text,"sum = %0.2f, firing rate %0.16f",r0->sum,r0->firingrate); glColor3f(1,1,1); ftexto(xoffset,yoffset-0.1,0,text);
}

GLvoid Drawrho1(double side,double xoffset,double yoffset)
{
  int verbose=0;
  char text[32];
  struct rho1 *r1=GLOBAL_rho1;
  int nv=0,ng=0;
  double xord=0,yord=0,rcolor=0,gcolor=0,bcolor=0,kcolor=0;
  double xord1=0,yord1=0,xord2=0,yord2=0;
  double gord=0,vord=0,gdot=0,vdot=0,norm=0,gdotn=0,vdotn=0;
  switch (abs(DRAW_FLAG)%3){
  case 0: /* draw flux arrows */
  for (nv=0;nv<r1->nv;nv++){ for (ng=0;ng<r1->ng;ng++){ 
    xord = (double)(ng+0.5)/(double)(r1->ng);
    yord = (double)(nv+0.5)/(double)(r1->nv);
    colorscale(7,r1->rho[nv+ng*r1->nv],4.0/(double)(r1->nv*r1->ng),0,&rcolor,&gcolor,&bcolor);
    glrect(xord*side+xoffset,yord*side+yoffset,0.5*side/(double)r1->ng,0.5*side/(double)r1->nv,rcolor,gcolor,bcolor);
    if (verbose){ sprintf(text,"%d,%d",nv,ng);glColor3f(1,0,1); ftexto(xord*side+xoffset,yord*side+yoffset,0,text);} }}
  for (nv=0;nv<r1->nv+1;nv++){ for (ng=0;ng<r1->ng;ng++){ 
    xord = (double)(ng+0.5)/(double)(r1->ng); yord = (double)(nv+0.0)/(double)(r1->nv);
    yord1 = yord-0.5/(double)(r1->nv); yord2 = yord+0.5/(double)(r1->nv); 
    xord = xord*side+xoffset; yord1 = yord1*side+yoffset; yord2 = yord2*side+yoffset;
    if (r1->fluxv[nv+ng*(r1->nv+1)]>0){ glarrow(1,xord,yord1,xord,yord2,0.125*side/(double)r1->nv,1,1,1);}
    else if (r1->fluxv[nv+ng*(r1->nv+1)]<0){ glarrow(1,xord,yord2,xord,yord1,0.125*side/(double)r1->nv,1,1,1);}}}
  for (nv=0;nv<r1->nv;nv++){ for (ng=0;ng<r1->ng+1;ng++){ 
    xord = (double)(ng+0.0)/(double)(r1->ng); yord = (double)(nv+0.5)/(double)(r1->nv);
    xord1 = xord-0.5/(double)(r1->ng); xord2 = xord+0.5/(double)(r1->ng); 
    yord = yord*side+yoffset; xord1 = xord1*side+xoffset; xord2 = xord2*side+xoffset;
    if (r1->fluxg[nv+ng*r1->nv]>0){ glarrow(1,xord1,yord,xord2,yord,0.125*side/(double)r1->ng,1,1,1);}
    else if (r1->fluxg[nv+ng*r1->nv]<0){ glarrow(1,xord2,yord,xord1,yord,0.125*side/(double)r1->ng,1,1,1);}}}
  break;
  case 1: /* only draw density */
  for (nv=0;nv<r1->nv;nv++){ for (ng=0;ng<r1->ng;ng++){ 
    xord = (double)(ng+0.5)/(double)(r1->ng);
    yord = (double)(nv+0.5)/(double)(r1->nv);
    colorscale(7,r1->rho[nv+ng*r1->nv],4.0/(double)(r1->nv*r1->ng),0,&rcolor,&gcolor,&bcolor);
    glrect(xord*side+xoffset,yord*side+yoffset,1.0*side/(double)r1->ng,1.0*side/(double)r1->nv,rcolor,gcolor,bcolor);}}
  break;
  case 2: /* only draw flux arrows */
    for (nv=0;nv<r1->nv;nv++){ for (ng=0;ng<r1->ng;ng++){
      xord = (double)(ng+0.5)/(double)r1->ng; yord = (double)(nv+0.5)/(double)r1->nv;
      vord = VOLTAGE_RESET + yord*(VOLTAGE_THRESHOLD_IF-VOLTAGE_RESET);
      gord = 0 + xord*(CONDUCTANCE_MAX-0);
      gdot = rho1_gdot(vord,gord); vdot = rho1_vdot(vord,gord);
      gdot = (gdot-0)/(CONDUCTANCE_MAX-0);
      vdot = (vdot-VOLTAGE_RESET)/(VOLTAGE_THRESHOLD_IF-VOLTAGE_RESET);
      norm = sqrt(pow(gdot,2)+pow(vdot,2)); 
      gdotn = gdot/norm/(double)r1->ng/sqrt(2); 
      vdotn = vdot/norm/(double)r1->nv/sqrt(2);
      glColor3f(1,1,1);
      glBegin(GL_LINES);
      glVertex3f((xord+gdotn)*side+xoffset,(yord+vdotn)*side+yoffset,0);
      glVertex3f((xord-gdotn)*side+xoffset,(yord-vdotn)*side+yoffset,0);
      glVertex3f((xord+gdotn)*side+xoffset,(yord+vdotn)*side+yoffset,0);
      glVertex3f((xord+0.5*gdotn-0.35*vdotn)*side+xoffset,(yord+0.5*vdotn+0.35*gdotn)*side+yoffset,0);
      glVertex3f((xord+gdotn)*side+xoffset,(yord+vdotn)*side+yoffset,0);
      glVertex3f((xord+0.5*gdotn+0.35*vdotn)*side+xoffset,(yord+0.5*vdotn-0.35*gdotn)*side+yoffset,0);
      glEnd();}}
    break;
  default: break;}
  sprintf(text,"sum = %f",r1->sum); glColor3f(1,1,1); ftexto(xoffset,yoffset-0.1,0,text);
}

GLvoid Drawrhovsa(double side,double xoffset,double yoffset)
{
  /* plots trajectory on v,sA plane */
  char text[64];
  int nt=0,nv=0,ng=0;
  double timetotal=1024*64,tlast=0;
  double t=0,V=0,sA=0;
  double xord=0,yord=0,rcolor=0,gcolor=0,bcolor=0,kcolor=0;
  double gord=0,vord=0,gdot=0,vdot=0,norm=0,gdotn=0,vdotn=0;
  struct neuron *n=GLOBAL_NEURON;
  struct litem *lt=NULL,*lV=NULL,*lA=NULL;
  double hmax=0,hmin=0,hmean=0,hstd=0;
  timetotal = GLOBAL_TF > 0 ? minimum(GLOBAL_TF,timetotal) : timetotal;
  for (nt=0;nt<n->types;nt++){
    stats("double",&(n->hist[0 + 0*n->nbins + nt*n->nbins*n->nbins]),n->nbins*(n->nbins-1),&hmax,&hmin,&hmean,&hstd);
    for (nv=0;nv<n->nbins;nv++){ for (ng=0;ng<n->nbins-1;ng++){
      if (n->hist[nv+ng*n->nbins+nt*n->nbins*n->nbins]>0){
	xord=(double)(ng/* +0.5 */)/(double)n->nbins;
	yord=(double)(nv/* +0.5 */)/(double)n->nbins;
	vord = VOLTAGE_RESET + yord*(VOLTAGE_THRESHOLD_IF-VOLTAGE_RESET);
	gord = 0 + xord*(2*GHAT-0);
	gdot = -gord/TAU_A;
	vdot = -CONDUCTANCE_LK*(vord-VOLTAGE_RESET) - gord*(vord-VOLTAGE_EX);
	gdot = (gdot-0)/(2*GHAT-0);
	vdot = (vdot-VOLTAGE_RESET)/(VOLTAGE_THRESHOLD_IF-VOLTAGE_RESET);
	norm = sqrt(pow(gdot,2)+pow(vdot,2)); 
	gdotn = gdot/norm/(double)n->nbins/sqrt(2); 
	vdotn = vdot/norm/(double)n->nbins/sqrt(2);
/* 	norm = gdotn; gdotn = -vdotn; vdotn = norm; */
	colorscale(7,log(n->hist[nv+ng*n->nbins+nt*n->nbins*n->nbins]),log(hmax),0,&rcolor,&gcolor,&bcolor);
	glColor3f(rcolor,gcolor,bcolor);
	glBegin(GL_LINES);
	glVertex3f((xord+gdotn)*side+xoffset-side,(yord+vdotn)*side+yoffset,0);
	glVertex3f((xord-gdotn)*side+xoffset-side,(yord-vdotn)*side+yoffset,0);
	glVertex3f((xord+gdotn)*side+xoffset-side,(yord+vdotn)*side+yoffset,0);
	glVertex3f((xord+0.5*gdotn-0.35*vdotn)*side+xoffset-side,(yord+0.5*vdotn+0.35*gdotn)*side+yoffset,0);
	glVertex3f((xord+gdotn)*side+xoffset-side,(yord+vdotn)*side+yoffset,0);
	glVertex3f((xord+0.5*gdotn+0.35*vdotn)*side+xoffset-side,(yord+0.5*vdotn-0.35*gdotn)*side+yoffset,0);
	glEnd();}}}
    sprintf(text,"max log(%d)=%f",(int)hmax,log(hmax));
    glColor3f(1,1,1);ftexto(xoffset,yoffset-0.1,0,text);    
    lt = n->t[nt]->last; lV = n->V[nt]->last; lA = n->sA->last;
    tlast = (lt==NULL ? 0 : *(double *)lt->item); t=tlast;
    if (0 && n->types>3){ colorscale(7,nt,n->types-1,0,&rcolor,&gcolor,&bcolor);}
    else{
      switch (nt){
      case 0: rcolor=1; gcolor=1; bcolor=1; break;
      case 1: rcolor=1; gcolor=0.5; bcolor=0; break;
      case 2: rcolor=0; gcolor=1; bcolor=1; break;
      default : colorscale(7,nt,n->types-1,0,&rcolor,&gcolor,&bcolor); break;}}
    glBegin(GL_POINTS);
    while (lt!=NULL && (tlast-t<timetotal)){
      t = *(double *)lt->item;
      V = *(double *)lV->item;
      sA = *(double *)lA->item;
      kcolor = 1.0-(tlast-t)/timetotal;
      glColor3f(rcolor*kcolor,gcolor*kcolor,bcolor*kcolor);
      xord = sA/(2*GHAT);
      yord = (V-VOLTAGE_RESET)/(VOLTAGE_THRESHOLD_IF-VOLTAGE_RESET);
      glVertex3f(xord*side+xoffset,yord*side+yoffset,0.0);
      lt = lt->parent;
      lV = lV->parent;
      lA = lA->parent;}
    glEnd();}
  glBegin(GL_LINE_LOOP);
  glColor3f(1,0,1);
  xord = 0; yord = 0; glVertex3f(xord*side+xoffset,yord*side+yoffset,0);
  xord = 0; yord = 1; glVertex3f(xord*side+xoffset,yord*side+yoffset,0);
  xord = 1; yord = 1; glVertex3f(xord*side+xoffset,yord*side+yoffset,0);
  xord = 1; yord = 0; glVertex3f(xord*side+xoffset,yord*side+yoffset,0);
  glEnd();
}

GLvoid Drawneuron(double side,double xoffset,double yoffset)
{
  //char text[32];
  int nt=0;
  double timetotal = GLOBAL_TF > 0 ? minimum(GLOBAL_TF,125) : 125,timescale = GLOBAL_TF > 0 ? minimum(GLOBAL_TF,50) : 50,tlast=0;
  double t=0,V=0,sA=0,sN=0,sG=0,gt=0,Vs=0;
  double xord=0,yord=0,rcolor=0,gcolor=0,bcolor=0;
  struct neuron *n=GLOBAL_NEURON;
  struct litem *lt=NULL,*lV=NULL,*lA=NULL,*lN=NULL,*lG=NULL;
  for (nt=0;nt<n->types;nt++){
    lt = n->t[nt]->last; lV = n->V[nt]->last; lA = n->sA->last; lN = n->sN->last; lG = n->sG->last;
    tlast = (lt==NULL ? 0 : *(double *)lt->item); t=tlast;
    if (0 && n->types>3){ colorscale(7,nt,n->types-1,0,&rcolor,&gcolor,&bcolor);}
    else{
      switch (nt){
      case 0: rcolor=1; gcolor=1; bcolor=1; break;
      case 1: rcolor=1; gcolor=0.5; bcolor=0; break;
      case 2: rcolor=0; gcolor=1; bcolor=1; break;
      default : colorscale(7,nt,n->types-1,0,&rcolor,&gcolor,&bcolor); break;}}
    while (lt!=NULL && (tlast-t<timetotal)){
      t = *(double *)lt->item;
      V = *(double *)lV->item;
      sA = *(double *)lA->item;
      sN = *(double *)lN->item;
      sG = *(double *)lG->item;
      gt = CONDUCTANCE_LK+sA+sN+sG;
      Vs = (CONDUCTANCE_LK*VOLTAGE_LEAK + (sA+sN)*VOLTAGE_EX + sG*VOLTAGE_IN)/gt;
      xord = (t-(tlast-timetotal))/timescale;
      yord = V;
      glbox(side*xord+xoffset,side*yord+yoffset,side/100,rcolor,gcolor,bcolor);
/*       yord = gt; */
/*       glbox(side*xord+xoffset,side*yord+yoffset,side/100,1,0,1); */
/*       yord = Vs; */
/*       glbox(side*xord+xoffset,side*yord+yoffset,side/100,1,1,0); */
      lt = lt->parent;
      lV = lV->parent;
      lA = lA->parent;
      lN = lN->parent;
      lG = lG->parent;}
    lt = n->t[nt]->last; lV = n->V[nt]->last; lA = n->sA->last; lN = n->sN->last; lG = n->sG->last;
    tlast = (lt==NULL ? 0 : *(double *)lt->item); t=tlast;
    glBegin(GL_LINE_STRIP);
    while (lt!=NULL && (tlast-t<timetotal)){
      t = *(double *)lt->item;
      V = *(double *)lV->item;
      sA = *(double *)lA->item;
      sN = *(double *)lN->item;
      sG = *(double *)lG->item;
      gt = CONDUCTANCE_LK+sA+sN+sG;
      Vs = (CONDUCTANCE_LK*VOLTAGE_LEAK + (sA+sN)*VOLTAGE_EX + sG*VOLTAGE_IN)/gt;
      xord = (t-(tlast-timetotal))/timescale;
      yord = V;
      glVertex3f(side*xord+xoffset,side*yord+yoffset,0);
      lt = lt->parent;
      lV = lV->parent;
      lA = lA->parent;
      lN = lN->parent;
      lG = lG->parent;}
    glEnd();}
  glBegin(GL_LINES);
  glColor3f(1,1,1);
  xord = 0; yord = 0; glVertex3f(xord*side+xoffset,yord*side+yoffset,0);
  xord = timetotal/timescale; yord = 0; glVertex3f(xord*side+xoffset,yord*side+yoffset,0);
  glColor3f(1,1,1);
  xord = 0; yord = VOLTAGE_THRESHOLD_IF; glVertex3f(xord*side+xoffset,yord*side+yoffset,0);
  xord = timetotal/timescale; yord = VOLTAGE_THRESHOLD_IF; glVertex3f(xord*side+xoffset,yord*side+yoffset,0);
/*   glColor3f(1,1,1); */
/*   xord = 0; yord = VOLTAGE_TAKEOFF; glVertex3f(xord*side+xoffset,yord*side+yoffset,0); */
/*   xord = timetotal/timescale; yord = VOLTAGE_TAKEOFF; glVertex3f(xord*side+xoffset,yord*side+yoffset,0); */
  glColor3f(1,1,1);
  xord = 0; yord = VOLTAGE_RESET; glVertex3f(xord*side+xoffset,yord*side+yoffset,0);
  xord = timetotal/timescale; yord = VOLTAGE_RESET; glVertex3f(xord*side+xoffset,yord*side+yoffset,0);
  glEnd();
}

GLvoid DrawGLScene(GLvoid)
{
  int index=0;
  for (index=0;index<STEPS_PER_DRAW;index++){ if (!RUN_DONE){ computestep();}}
  /* draw stuff */
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // Clear The Screen And The Depth Buffer
  glLoadIdentity(); // Reset The View
  glTranslatef(xdepth,ydepth,zdepth);
  switch (GLOBAL_RHO_VERSION){
  case 0: Drawrho0(1,0,0); break;
  case 1: Drawrho1(1,0,0); break;
  default: break;}
  Drawmenu(.1,-1,0);
  glutSwapBuffers(); /* since double buffered */
}

void keyPressed(unsigned char key, int x, int y) 
{
  /* The function called whenever a normal key is pressed. */
  usleep(100); /* still not sure why this is called */
  switch (key) {    
  case ESCAPE_KEY: /* stop everything */
    glutDestroyWindow(GLUTWINDOWNUMBER); /* shut down our window */
    printf("(x,y,z) = (%f,%f,%f)\n",xdepth,ydepth,zdepth);
    exit(EXIT_SUCCESS); /* exit the program...normal termination. */
    break; /* just in case */
  case SPACE_KEY: /* stop/restart temporarily */
    STEPS_PER_DRAW=!STEPS_PER_DRAW;if (!STEPS_PER_DRAW){ printf("PROCESS HALTED\n");}
    break;
  case 'w': FIDDLE_PARAMETER++; if (FIDDLE_PARAMETER > 22){ FIDDLE_PARAMETER = -1;} break;
  case 's': FIDDLE_PARAMETER--; if (FIDDLE_PARAMETER < -1){ FIDDLE_PARAMETER = 22;} break;
  case 'd': 
    switch (FIDDLE_PARAMETER){
    case 0: STEPS_PER_DRAW *= 2; if (STEPS_PER_DRAW==0){ STEPS_PER_DRAW=1;} break;
    case 1: GLOBAL_DT *= 2; break;
    case 2: TAU_A *= 1.1; break;
    case 3: AMPA_STRENGTH *= 1.1; break;
    case 4: AMPA_RATE *= 1.1; break;
    case 5: DRAW_FLAG += 1; break;
    default: break;}
    break;
  case 'a':
    switch (FIDDLE_PARAMETER){
    case 0: STEPS_PER_DRAW /= 2; break;
    case 1: GLOBAL_DT /= 2; break;
    case 2: TAU_A /= 1.1; break;
    case 3: AMPA_STRENGTH /= 1.1; break;
    case 4: AMPA_RATE /= 1.1; break;
    case 5: DRAW_FLAG -= 1; break;
    default: break;}
    break;
  default:
    break;}
}

void specialKeyPressed(int key, int x, int y) 
{
  /* The function called whenever a special key is pressed. */
  int mod=0;
  usleep(100);
  switch (key) {    
  case GLUT_KEY_PAGE_UP:
    zdepth /= 1.05;
    break;
  case GLUT_KEY_PAGE_DOWN:
    zdepth *= 1.05;
    break;
  case GLUT_KEY_UP:
    mod = glutGetModifiers();
    if (mod==GLUT_ACTIVE_SHIFT){}
    else{ ydepth+= -0.05*zdepth;}
    break;
  case GLUT_KEY_DOWN:
    mod = glutGetModifiers();
    if (mod==GLUT_ACTIVE_SHIFT){}
    else{ ydepth-= -0.05*zdepth;}
    break;
  case GLUT_KEY_LEFT:
    mod = glutGetModifiers();
    if (mod==GLUT_ACTIVE_SHIFT){}
    else{ xdepth-= -0.05*zdepth;}
    break;
  case GLUT_KEY_RIGHT:
    mod = glutGetModifiers();
    if (mod==GLUT_ACTIVE_SHIFT){}
    else{ xdepth+= -0.05*zdepth;}
    break;
  case GLUT_KEY_END:
    mod = glutGetModifiers();
    if (!STEPS_PER_DRAW){ computestep();}
    break;
  case GLUT_KEY_HOME:
    mod = glutGetModifiers();
    if (mod==GLUT_ACTIVE_SHIFT){}
    else { neurontfree(GLOBAL_NEURON); GLOBAL_NEURON = neuronmake(GLOBAL_NTYPES);}
    break;
  default:
    break;}
}

#endif /* NOTGRAPH */

/* Here are the global functions */

double drand(int xsym){ return (xsym ? ((double)rand()/(double)RAND_MAX - 0.5)*2 : (double)rand()/(double)RAND_MAX);}

int processoptions(int argc, char ** argv)
{
  int still_have_options = argc-1;
  while (still_have_options){
    switch(getopt(argc,argv,"d:D:f:F:")){
    case 'd': case 'D': GLOBAL_DT = atof(optarg); break;
    case 'f': case 'F': GLOBAL_TF = atof(optarg); break;
    default: exit(EXIT_SUCCESS); break;}
    still_have_options--;}
  return 1;
}

int main(int argc, char **argv) 
{
  setglobals(); /* initialize global variables */
  processoptions(argc,argv);
#ifndef NOTGRAPH
  glutInit(&argc, argv);  /* initialize glut */
  glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_ALPHA | GLUT_DEPTH); /* select display mode, 2xbuffer,rgba,alpha,depth_buffer */
  glutInitWindowSize(GLOBAL_WINDOW_WIDTH, GLOBAL_WINDOW_HEIGHT); /* small window */
  glutInitWindowPosition(0,0); /* starts at upper left corner of screen */
  GLUTWINDOWNUMBER = glutCreateWindow("solversforif_snx"); /* open window */
  glutDisplayFunc(&DrawGLScene); /* this does all the drawing for us */
  glutFullScreen(); GLOBAL_WINDOW_WIDTH=1024; GLOBAL_WINDOW_HEIGHT=768; /* if we wish to fullscreen... */
  glutIdleFunc(&DrawGLScene); /* continuously redraw */
  glutReshapeFunc(&ReSizeGLScene); /* resizing function */
  glutKeyboardFunc(&keyPressed);glutSpecialFunc(&specialKeyPressed); /* input functions */
  InitGL(GLOBAL_WINDOW_WIDTH,GLOBAL_WINDOW_HEIGHT); /* init window */
  glutMainLoop(); /* start processing */
#endif /* NOTGRAPH */
#ifdef NOTGRAPH
  while (!RUN_DONE){ computestep();}
#endif /* NOTGRAPH */
  return 1;
}
