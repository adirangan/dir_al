/* simple al */
#ifndef NOTGRAPH
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#endif /* NOTGRAPH */
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "blaswrap.h"
#include "f2c.h"
#include "clapack.h"
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
#include "al_ptreehandling.h"
#include "al_cortex.h"
#include "al_datahandling.h"
#include "al_inputhandling.h"
#include "al_ogletools.h"

/* Here are the method functions */
void setglobals();
void setdependencies();
void computestep();

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
int GLOBAL_DRAW_FLAG=1;
GLfloat xdepth = 0.0f;
GLfloat ydepth = 0.0f;
GLfloat zdepth = 0.0f;
GLfloat RXangle = 0.0f;
GLfloat RYangle = 0.0f;
GLfloat RXdangle = 0.0f;
GLfloat RYdangle = 0.0f;
#endif /* NOTGRAPH */
int ON_MY_COMPUTER=0;
char GLOBAL_STRING[32];
char GLOBAL_STRING_2[64];
int GLOBAL_RECORD_NUMBER=0;
int GLOBAL_SPIKEINPUT_RSEED=313;
int RUN_DONE=0;
struct neuronarray *GLOBAL_Nra=NULL;
int GLOBAL_NTYPES=2;
char **GLOBAL_TYPENAMES=NULL;
int GLOBAL_NVARS=8;
int GLOBAL_NCLUSTERS=1;
double *GLOBAL_CLUSTER_PRA_=NULL;
int GLOBAL_LINK_SPARSE_OR_DENSE=1;
int GLOBAL_CLUSTER_MAKE_OR_READ=1;
char GLOBAL_CLUSTER_READ_FILE[512];
double GLOBAL_REWEIGHT_STRENGTH=0.1;
char GLOBAL_REWEIGHT_READ_FILE[512];
int GLOBAL_REWEIGHT_MULTIPLY_OR_ADD=1;
int GLOBAL_NEURON_MODEL=1;
int GLOBAL_SNX_VERSION=0;
int *GLOBAL_SNX_NSTATES_=NULL;
char **GLOBAL_VARNAMES=NULL;
/* all global_indexing is -1 padded */
int GLOBAL_INDEXING_sra_LENGTH=0;
int *GLOBAL_INDEXING_REFILE_sra=NULL;
int *GLOBAL_INDEXING_CHECKOUT_sra=NULL;
int *GLOBAL_LENGTHRA=NULL;
double *VOLTAGE_=NULL;
double VOLTAGE_THRESHOLD_S=0;
double VOLTAGE_THRESHOLD_D=0;
double *CONDUCTANCE_=NULL;
double CONDUCTANCE_SD=0.05;
double CONDUCTANCE_DS=0.05;
double CURRENT_INJECTION_S=0.0;
double CURRENT_INJECTION_D=0.0;
double TAU_REF=2;
double *TAU_=NULL;
int *AUTAPSES_OFF=NULL;
double GLOBAL_CS_ORN_SCALE=1;
double *CS_ORN_=NULL;
double *CS_ORN_mainak_stim_=NULL;
double *CS_ORN_wilson_stim_=NULL;
double GLOBAL_CS_SCALE= 1;
double *GLOBAL_CS_PRETYPE_SCALE_= NULL;
double *GLOBAL_CS_POSTYPE_SCALE_= NULL;
double *GLOBAL_CS_SRA_SCALE_= NULL;
double *CS__=NULL;
double *SPARSE__=NULL;
double *SPARSEHIT__=NULL;
double *SPARSE_OGLI__=NULL;
double *SPARSEHIT_OGLI__=NULL;
double GLOBAL_P_FAIL_SCALE= 1;
double *P_FAIL__=NULL;
struct odor *GLOBAL_ODORra_START=NULL;
struct odor *GLOBAL_ODORra=NULL;
int GLOBAL_ODOR_BASE=2;
double INPUT_CONTRAST_START=0;
double INPUT_CONTRAST=0;
double INPUT_PULSE=1;
int ORN_BOTHER=0;
struct orn *GLOBAL_ORN=NULL;
double *GLOBAL_ORN_TIMES_=NULL;
double ORN_BACKRATE=0;
int SPIKETOL=0;
int GLOBAL_verbose=0,GLOBAL_dtok=0,GLOBAL_dtadapt=0;
double GLOBAL_time=0;
double GLOBAL_TI=0;
double GLOBAL_TF=0;
double GLOBAL_DT=0;
double GLOBAL_DTmax=0;
double GLOBAL_DTmin=0;
int FIDDLE_PARAMETER=0;
int FIDDLE_ROW=0;
int FIDDLE_COL=0;
int DRAW_FLAG=0;
int DRAW_FLAG2=0;
int STEPS_PER_DRAW=0;
double STD_VIEW=1;
double MEAN_VIEW=0;
int GRAYSCALE=0;
int GLOBAL_SPACE_SMOOTHER=0;
int SUPERGLOBAL_DRAW_FLAG=0;
double OUTPUT_DUMP_EVERY=0;
int POWER_BOTHER=0;
int GLOBAL_POWER_LENGTH=0;
struct power *GLOBAL_POWER=NULL;
int GLOBAL_POWER_INDEXING_NTYPE_LENGTH=0;
int *GLOBAL_POWER_INDEXING_NTYPE_CHECKOUT=NULL;
int *GLOBAL_POWER_INDEXING_NTYPE_REFILE=NULL;
int GLOBAL_POWER_INDEXING_NVAR_LENGTH=0;
int *GLOBAL_POWER_INDEXING_NVAR_CHECKOUT=NULL;
int *GLOBAL_POWER_INDEXING_NVAR_REFILE=NULL;
int GLOBAL_POWER_TRAJECTORY_LENGTH=1;
int GLOBAL_POWER_WINDOW_LENGTH=0;
int GLOBAL_POWER_WINDOW_UPDATE_EVERY=1;
double *GLOBAL_POWER_maxra_=NULL;
double *GLOBAL_POWER_minra_=NULL;
int GLOBAL_POWER_CYCLE_BOTHER=0;
int GLOBAL_POWER_CORRELATION_BOTHER=0;
int RHO_BOTHER=0;
struct rho *GLOBAL_RHO=NULL;
int GLOBAL_RHO_LENGTH=0;
int GLOBAL_RHO_INDEXING_NVAR_LENGTH=0;
int *GLOBAL_RHO_INDEXING_NVAR_CHECKOUT=NULL;
int *GLOBAL_RHO_INDEXING_NVAR_REFILE=NULL;
int *GLOBAL_RHO_NBINRA_=NULL;
int PTREE_BOTHER=0;
int GLOBAL_PTREE_BITBYBIT=0;
int GLOBAL_PTREE_BITBYBIT_CURRENT_NUMBER=0;
int GLOBAL_PTREE_ZZZ=1;
int GLOBAL_PTREE_NREGIONS=8;
int GLOBAL_PTREE_NLEGS=0;
int GLOBAL_PTREE_LEGTIME=0;
int GLOBAL_PTREE_EVENT_THRESHOLD=1;
double GLOBAL_PTREE_EVENT_WITHIN=2;
int GLOBAL_PTREE_REGION_TYPE=1;
int GLOBAL_PTREE_DUMP_TYPE=2;
struct ptree *GLOBAL_PTREE=NULL;
struct ptree *GLOBAL_INPUT_PTREE=NULL;
int CLUSTERDATA_BOTHER=0;
struct clusterdatara *GLOBAL_CDRA=NULL;
int SUITE_BOTHER=0;
int SUITE_NODORS=2;
struct odor *SUITE_ODORra_BACON=NULL;
double SUITE_CS_ORN_SCALE=1.1;
double SUITE_INPUT_CONTRAST=1.1;
int SUITE_NCONCENTRATIONS=1;
int SUITE_NBICUCULLINES=2;
int SUITE_NINSTANCES=1;
int SUITE_POWER_BOTHER=0;
int SUITE_PTREE_BOTHER=0;
double SUITE_4_REWEIGHT_STRENGTH=0;
int SUITE_7_CLEANUP=1;
int SUITE_8_CLEANUP=1;
double SUITE_8_wilson_axon_stim=1;
int SUITE_CLUSTERDATA_BOTHER=0;
int SUITE_NSECONDS=64;
int SUITE_DUMPEVERY=64;
int SUITE_BITBYBIT_RECORD=1;
int SUITE_BITBYBIT_REMOVE=1;
int SUITE_SINDEXMAX=0;
int SUITE_DINDEXMAX=0;
int SUITE_TINDEXMAX=1;
double SUITE_ONSET_TIME=0;
int LYAPUNOV_BOTHER=0;
int CAICOR_BOTHER=0;
struct caicor *GLOBAL_CAICOR=NULL;
int GLOBAL_CAICOR_NBINS=0;
int ISI_BOTHER=0;
struct isi *GLOBAL_ISI=NULL;
double GLOBAL_ISI_MAX=0;
double GLOBAL_ISI_MIN=0;
int GLOBAL_ISI_NBINS=0;
int LATTICE3D_BOTHER=0;
struct lattice3d *GLOBAL_LATTICE3D=NULL;
int GLOBAL_LATTICE3D_NBINS=0;
int GLOBAL_LATTICE3D_LENGTH=0;
int GLOBAL_LATTICE3D_DEPTH=0;
int GLOBAL_LATTICE3D_TET_VS_CUBE=0;
double *GLOBAL_LATTICE3D_CS__=NULL;
double *GLOBAL_LATTICE3D_CS_=NULL;
int HHLIB_BOTHER=0;
int GLOBAL_HHLIB_MAKE_OR_READ=0;
char GLOBAL_HHLIB_READ_FILE[512];
struct hhlib *GLOBAL_HHLIB=NULL;
double GLOBAL_HHLIB_TAU_SKIP=1.0;
int GLOBAL_HHLIB_INDEXING_TRIGGER=0;
int *GLOBAL_HHLIB_LOGFLAGRA=NULL;
int *GLOBAL_HHLIB_LIBUSEFLAGRA=NULL;
int GLOBAL_HHLIB_NBINS=0;
int GLOBAL_HHLIB_INDEXING_NVAR_LENGTH=0;
int *GLOBAL_HHLIB_INDEXING_NVAR_CHECKOUT=NULL;
int *GLOBAL_HHLIB_INDEXING_NVAR_REFILE=NULL;
double *GLOBAL_HHLIB_maxra_=NULL;
double *GLOBAL_HHLIB_minra_=NULL;
int GLOBAL_HHLIB_USE_AFTER=0;
int SNXDATA_BOTHER=0;
struct snxdata *GLOBAL_SNXDATA=NULL;
int GLOBAL_SNXDATA_LOOKBACK=1;
double GLOBAL_VESICLE_DEPLETION=0;

#include "al_inputhandling.c"
#include "d_llists.c"
#include "matrices.c"
#include "al_cortex.c"
#include "al_datahandling.c"
#include "al_ptreehandling.c"
#include "al_ogletools.c"

/* Here are the method functions */

void setglobals()
{
  int verbose=GLOBAL_verbose;
  int nt=0,nr=0,nc=0;
#ifndef NOTGRAPH
  GLOBAL_WINDOW_WIDTH = 640;
  GLOBAL_WINDOW_HEIGHT = 480;
  xdepth = -2.2f;
  ydepth = -1.2f;
  zdepth = -6.0f;
#endif /* NOTGRAPH */
  GLOBAL_time=GLOBAL_TI;
  GLOBAL_DT=GLOBAL_DTmax;
  GLOBAL_DTmin = GLOBAL_DTmax*pow(2,-10);
  GLOBAL_dtok=0;
  GLOBAL_ODORra = odormake(0,NULL,0,NULL,GLOBAL_ODORra_START);
  INPUT_CONTRAST=INPUT_CONTRAST_START;
  if (ORN_BOTHER){
    if (verbose){ printf(" %% making orn\n");}
    GLOBAL_ORN = ornmake(GLOBAL_ORN_TIMES_[0],GLOBAL_ORN_TIMES_[1],GLOBAL_ORN_TIMES_[2],GLOBAL_ORN_TIMES_[3]);}
  if (verbose){ printf(" %% making neuronarray\n");}
  GLOBAL_Nra = neuronarraymake(GLOBAL_NTYPES,GLOBAL_LENGTHRA,GLOBAL_NVARS,GLOBAL_INDEXING_sra_LENGTH,GLOBAL_NCLUSTERS);
/*   int nt=0,nt2=0,ns=0,tab=0; */
/*   for (nt=0;nt<GLOBAL_NTYPES;nt++){ for (nt2=0;nt2<GLOBAL_NTYPES;nt2++){ for (ns=0;ns<GLOBAL_INDEXING_sra_LENGTH;ns++){ */
/*     tab = nt + nt2*GLOBAL_NTYPES + ns*GLOBAL_NTYPES*GLOBAL_NTYPES; */
/*     if (GLOBAL_CS_PRETYPE_SCALE_ != NULL){ CS__[tab] *= GLOBAL_CS_PRETYPE_SCALE_[nt];} */
/*     if (GLOBAL_CS_POSTYPE_SCALE_ != NULL){ CS__[tab] *= GLOBAL_CS_POSTYPE_SCALE_[nt2];} */
/*     if (GLOBAL_CS_SRA_SCALE_ !=NULL){ CS__[tab] *= GLOBAL_CS_SRA_SCALE_[ns];} */
/*     CS__[tab] *= GLOBAL_CS_SCALE; */
/*     if (verbose){ printf(" %% CS__[%s-->%s (%s)]=%f\n",GLOBAL_TYPENAMES[nt],GLOBAL_TYPENAMES[nt2],GLOBAL_VARNAMES[GLOBAL_INDEXING_CHECKOUT_sra[ns]],CS__[tab]);}}}} */
  for (nt=0;nt<GLOBAL_NTYPES*GLOBAL_NTYPES*GLOBAL_INDEXING_sra_LENGTH;nt++){ P_FAIL__[nt] = maximum(0,minimum(1,P_FAIL__[nt]*GLOBAL_P_FAIL_SCALE));}
  if (verbose){ printf(" %% setting up data structures\n");}
  if (POWER_BOTHER){ GLOBAL_POWER = powermake(GLOBAL_Nra,GLOBAL_POWER_LENGTH,GLOBAL_POWER_INDEXING_NTYPE_LENGTH,GLOBAL_POWER_INDEXING_NTYPE_CHECKOUT,GLOBAL_POWER_INDEXING_NTYPE_REFILE,GLOBAL_POWER_INDEXING_NVAR_LENGTH,GLOBAL_POWER_INDEXING_NVAR_CHECKOUT,GLOBAL_POWER_INDEXING_NVAR_REFILE,GLOBAL_POWER_maxra_,GLOBAL_POWER_minra_,GLOBAL_POWER_CORRELATION_BOTHER);}
  if (RHO_BOTHER){ GLOBAL_RHO = rhomake(GLOBAL_Nra,GLOBAL_RHO_LENGTH,GLOBAL_RHO_INDEXING_NVAR_LENGTH,GLOBAL_RHO_INDEXING_NVAR_CHECKOUT,GLOBAL_RHO_INDEXING_NVAR_REFILE,GLOBAL_POWER_maxra_,GLOBAL_POWER_minra_,GLOBAL_RHO_NBINRA_);}
  if (PTREE_BOTHER){ GLOBAL_PTREE = ptreemake(GLOBAL_PTREE_NREGIONS,GLOBAL_PTREE_EVENT_WITHIN,GLOBAL_PTREE_EVENT_THRESHOLD,GLOBAL_PTREE_REGION_TYPE,GLOBAL_PTREE_NLEGS,GLOBAL_PTREE_LEGTIME);}
  if (CLUSTERDATA_BOTHER){ GLOBAL_CDRA = clusterdataramake(GLOBAL_Nra->gli,1,GLOBAL_POWER_LENGTH,GLOBAL_POWER_INDEXING_NVAR_LENGTH,GLOBAL_POWER_INDEXING_NVAR_CHECKOUT,GLOBAL_Nra->nvars,GLOBAL_POWER_INDEXING_NVAR_REFILE,GLOBAL_POWER_TRAJECTORY_LENGTH,GLOBAL_POWER_WINDOW_LENGTH,GLOBAL_POWER_WINDOW_UPDATE_EVERY,GLOBAL_POWER_maxra_,GLOBAL_POWER_minra_,GLOBAL_POWER_CYCLE_BOTHER,1,GLOBAL_PTREE_NLEGS,GLOBAL_PTREE_LEGTIME);}
  if (CAICOR_BOTHER){ GLOBAL_CAICOR = caicormake(GLOBAL_Nra,GLOBAL_POWER_INDEXING_NVAR_LENGTH,GLOBAL_POWER_INDEXING_NVAR_CHECKOUT,GLOBAL_Nra->nvars,GLOBAL_POWER_INDEXING_NVAR_REFILE,GLOBAL_POWER_maxra_,GLOBAL_POWER_minra_,GLOBAL_CAICOR_NBINS);}
  if (HHLIB_BOTHER){ if (GLOBAL_HHLIB_MAKE_OR_READ==0){ GLOBAL_HHLIB = hhlibread(GLOBAL_Nra,GLOBAL_HHLIB_READ_FILE); if (verbose){ for (nt=0;nt<GLOBAL_Nra->ntypes;nt++){ hhlib_project_printf(GLOBAL_HHLIB,nt,0);}}} else /* if (GLOBAL_HHLIB_MAKE_OR_READ!=0) */{ GLOBAL_HHLIB = hhlibmake(GLOBAL_Nra,GLOBAL_HHLIB_TAU_SKIP,GLOBAL_HHLIB_INDEXING_TRIGGER,GLOBAL_HHLIB_INDEXING_NVAR_LENGTH,GLOBAL_HHLIB_INDEXING_NVAR_CHECKOUT,GLOBAL_Nra->nvars,GLOBAL_HHLIB_INDEXING_NVAR_REFILE,GLOBAL_HHLIB_LOGFLAGRA,GLOBAL_HHLIB_LIBUSEFLAGRA,GLOBAL_HHLIB_maxra_,GLOBAL_HHLIB_minra_,GLOBAL_HHLIB_NBINS,GLOBAL_HHLIB_USE_AFTER);}}
  if (SNXDATA_BOTHER){ GLOBAL_SNXDATA = snxdatamake(GLOBAL_Nra,GLOBAL_SNXDATA_LOOKBACK);}
  if (ISI_BOTHER){ GLOBAL_ISI = isimake(GLOBAL_Nra,GLOBAL_ISI_NBINS,GLOBAL_ISI_MAX,GLOBAL_ISI_MIN);}
  if (LATTICE3D_BOTHER){
    if (GLOBAL_LATTICE3D_CS_!=NULL){
      if (GLOBAL_LATTICE3D_CS__==NULL){ GLOBAL_LATTICE3D_CS__=(double *)tcalloc(3*3,sizeof(double));}
      for (nr=0;nr<3;nr++){ for (nc=0;nc<3;nc++){
	GLOBAL_LATTICE3D_CS__[nr+nc*3] = (nr==nc?0:1)*GLOBAL_LATTICE3D_CS_[nr];}}}
    if (verbose){ raprintf(GLOBAL_LATTICE3D_CS__,"double",3,3,"global_latticd3d_cs__: ");}
    GLOBAL_LATTICE3D = threetree_test(GLOBAL_LATTICE3D_NBINS,GLOBAL_LATTICE3D_TET_VS_CUBE,GLOBAL_LATTICE3D_CS__); 
    exit(0);}
  if (verbose){ printf(" %% finished setglobals\n");}
}

void setdependencies()
{
  /* merely enforces global dependencies */
  int verbose=1;
  void *mention;
  char filename[2048];
  int continue_flag=0;
  FILE *fp=NULL;
  mention = &GLOBAL_QUADRATURE_N;
  mention = &GLOBAL_QUADRATURE_W;
  mention = &GLOBAL_QUADRATURE_Z;
  mention = &FIG_PREAMBLE_COLOR_0;
  mention = &FIG_PREAMBLE_COLOR_5;
  mention = &FIG_PREAMBLE_COLOR_7;
  if (GLOBAL_RECORD_NUMBER>0){ 
    if (verbose){ printf(" %% searching for global_record_number...\n");}
    do{
      sprintf(filename,"./%s_%d_key.log",GLOBAL_STRING,GLOBAL_RECORD_NUMBER);
      continue_flag = checktofind(filename); if (continue_flag){ GLOBAL_RECORD_NUMBER+=1;}}
    while (continue_flag);
    sprintf(filename,"./%s_%d_key.log",GLOBAL_STRING,GLOBAL_RECORD_NUMBER);
    if ((fp=fopen(filename,"w"))==NULL){ printf(" %% Warning! can't open %s in setdependencies\n",filename); fp=stdout;}
    fprintf(fp,"%s_%d_key.log",GLOBAL_STRING,GLOBAL_RECORD_NUMBER);
    if (fp!=stdout){ fclose(fp);fp=NULL;}
    if (verbose){ printf(" found %d\n",GLOBAL_RECORD_NUMBER);}
    srand(GLOBAL_RECORD_NUMBER);}
  sprintf(GLOBAL_STRING_2,"%s_%d_",GLOBAL_STRING,GLOBAL_RECORD_NUMBER); 
  if (verbose){ printf(" %s\n",GLOBAL_STRING_2);}
}

void computestep()
{
  int verbose=GLOBAL_verbose,dtok=0,pspikes=0,sspikes=0,gspikes=0;
  double t = GLOBAL_time,DT=GLOBAL_DT,DTmax=GLOBAL_DTmax;
  struct neuronarray *Nra=NULL;
  struct neuron *n=NULL;
  struct llist *LS=NULL;
  struct llist *LG=NULL;
  int nt=0,nr=0,nv=0;
  if (!RUN_DONE){
    dtok = GLOBAL_dtok; Nra = GLOBAL_Nra;
    if (verbose){ memprintf(0);}
    if (dtok > 4 && GLOBAL_dtadapt){ if (verbose){ printf("%% DT=%e has been successful for %d steps, doubling to make DT=%e\n",DT,dtok,2*DT);} DT=2*DT;}
    if (DT > DTmax){ if (verbose){ printf(" %% on second thought, DT large enough... setting DT=%e\n",DTmax);} DT=DTmax;}
    if (verbose){ printf("%% setting input spiketimes\n");}
    spikeinput(Nra,t,DT);
    if (GLOBAL_NEURON_MODEL==0 || GLOBAL_NEURON_MODEL==1 || GLOBAL_NEURON_MODEL==2 || GLOBAL_NEURON_MODEL==3){
      if (verbose){ printf("%% scanning for potential spikes\n");}
      do{
	if (verbose){ printf("%% \t unmaking previous llists\n");} 
	if (LS!=NULL){ llisttfree(LS); LS=NULL;} if (LG!=NULL){ llisttfree(LG); LG=NULL;}
	if (verbose){ printf("%% \t remaking llists\n");} 
	LS=llistmake();LG=llistmake();
	if (verbose){ printf("%% \t scanning...\n");}
	pspikes = spikescan(Nra,LS,LG,t,DT);
	if (verbose){ printf("%% \t found %d potential spikes\n",pspikes);}
	if (pspikes > SPIKETOL && GLOBAL_dtadapt){ if (verbose){ printf("%% \t categorically too many\n");}}
	else{ 
	  if (verbose){ printf("%% \t sorting ...\n");}
	  llistsort(LS->first,LS->last,LS->length,&spiketime_compare);
	  if (verbose){ printf("%% \t ...finished sorting\n");}}
	if (pspikes > SPIKETOL && GLOBAL_dtadapt){
	  if (verbose){ printf(" %%  pspikes=%d > %d, halving DT=%e to be DT/2=%e\n",pspikes,SPIKETOL,DT,DT/2);}
	  DT = DT/2;
	  dtok=0;}
	if (pspikes<=SPIKETOL || !GLOBAL_dtadapt || DT<=GLOBAL_DTmin){
	  if (verbose){ printf(" %% correcting spiketimes locally\n");}
	  for (nr=0;nr<pspikes;nr++){
	    llistsort(LS->first,LS->last,LS->length,&spiketime_compare); clumpcorrect(LS,LS,t,DT);}}
	llistsort(LS->first,LS->last,LS->length,&spiketime_compare); clumpcorrect(LG,LS,t,DT);
	gspikes = gizmospiked(LG,t,DT);
	if (gspikes > 0 && GLOBAL_dtadapt){
	  if (verbose){ printf(" %% gizmos fired %d times, halving DT=%e to DT/2=%e\n",gspikes,DT,DT/2);}
	  DT = DT/2;
	  dtok=0;}}
      while ((pspikes > SPIKETOL || gspikes>0) && (GLOBAL_dtadapt && DT>GLOBAL_DTmin));
      if (verbose){ printf(" %% step accepted with t=%0.3f,DT=%0.3f,pspikes=%d\n",t,DT,pspikes);}
      if (verbose){ printf(" %% evolving all neurons \n");}
      sspikes=0;spikeconduct(LS,t,DT,&sspikes,NULL);
      gizmoconduct(Nra,LS,LG,t,DT);
      if (LS!=NULL){ llisttfree(LS); LS=NULL;} if (LG!=NULL){ llisttfree(LG); LG=NULL;}}
    else if (GLOBAL_NEURON_MODEL==4){
      if (verbose){ printf(" %% evolving one spike at a time\n");}
      LS=llistmake(); LG=llistmake();
      for (nt=0;nt<Nra->ntypes;nt++){ for (nr=0;nr<Nra->lengthra[nt];nr++){ 
	n=nget(Nra,nt,nr);
	for (nv=0;nv<n->nvars;nv++){ n->vra[nv] = *(n->vpra[nv]);}
	litemadd(LS,nget(Nra,nt,nr)); litemadd(LG,nget(Nra,nt,nr));}}
      nllistevolve_cif(LS,LG,t,DT);
      spikeconduct(LS,t,DT,&sspikes,NULL);
      llisttfree(LS);LS=NULL; llisttfree(LG);LG=NULL;}
    else if (GLOBAL_NEURON_MODEL==5){
      if (verbose){ printf(" %% evolving mainak style\n");}
      mainak_continuous_evolve(Nra,t,DT);}
    else if (GLOBAL_NEURON_MODEL==6){
      if (verbose){ printf(" %% evolving wilson style\n");}
      wilson_evolve(Nra,t,DT);}
    else if (GLOBAL_NEURON_MODEL==7){
      if (verbose){ printf(" %% evolving snx style\n");}
      snx_evolve(Nra,t,DT);}
    else /* if GLOBAL_NEURON_MODEL not defined */{ printf(" %% warning! GLOBAL_NEURON_MODEL %d not defined\n",GLOBAL_NEURON_MODEL);}
    if (verbose){ printf(" %% updating data structures...");}
    if (ORN_BOTHER){ ornevolve(GLOBAL_ORN,t,DT);}
    if (POWER_BOTHER){ powerupdate(GLOBAL_POWER,t,DT);}
    if (RHO_BOTHER){ rhoupdate(GLOBAL_RHO,t,DT);}
    if (PTREE_BOTHER){ ptreeupdate(GLOBAL_PTREE,t,DT,1); if (GLOBAL_INPUT_PTREE!=NULL){ ptreeupdate(GLOBAL_INPUT_PTREE,t,DT,1);}}
    if (CLUSTERDATA_BOTHER){ clusterdataraupdate(GLOBAL_CDRA,t,DT);}
    if (CAICOR_BOTHER){ caicorupdate(GLOBAL_CAICOR,t,DT,sspikes);}
    if (HHLIB_BOTHER){ hhlibupdate(GLOBAL_HHLIB,t,DT);}
    if (SNXDATA_BOTHER){ snxdataupdate(GLOBAL_SNXDATA,t,DT);}
    if (verbose){ printf(" finished\n");}}
  if (verbose){ printf("\n %%%% finished for time interval [%f,%f] %%%% \n\n",t,t+DT);}
  if (SUITE_BOTHER){ system_monitor(Nra,t,DT);} 
  if (verbose){ memprintf(0);}
  t = t+DT;
  dtok++;
  GLOBAL_dtok=dtok;
  GLOBAL_time=t;
  GLOBAL_DT=DT;
  setinputrate(GLOBAL_Nra,GLOBAL_time);
  if (!SUITE_BOTHER && ((OUTPUT_DUMP_EVERY>0 && (int)floor(GLOBAL_time/OUTPUT_DUMP_EVERY)<(int)floor((GLOBAL_time+GLOBAL_DT)/OUTPUT_DUMP_EVERY)) || (GLOBAL_TF > 0 && GLOBAL_time > GLOBAL_TF))){
      printf("dumping output at time %f...",GLOBAL_time);
      if (PTREE_BOTHER){
	pnode_obs2dist_starter(0,1,0,NULL,GLOBAL_PTREE->postree,-1,0,1);
	if (verbose){ printf("dumping ptree... ");}
	ptreerate(GLOBAL_PTREE); 
	ptreedump_starter(GLOBAL_PTREE,NULL,GLOBAL_PTREE_DUMP_TYPE,0,0,0,+1,-1,0);
	if (GLOBAL_PTREE_BITBYBIT){ GLOBAL_PTREE_BITBYBIT_CURRENT_NUMBER += 1;}
	if (verbose){ printf("clearing ptree... ");}
	pnodeclear_starter(NULL,GLOBAL_PTREE->pretree); pnodeclear_starter(NULL,GLOBAL_PTREE->postree); GLOBAL_PTREE->total_time=0;
	if (verbose){ printf("\n");}}
      if (POWER_BOTHER){ powerdump(GLOBAL_POWER,NULL,0); powerdump(GLOBAL_POWER,NULL,1);}
      if (RHO_BOTHER){ rhodump(GLOBAL_RHO,NULL,0);}
      if (CLUSTERDATA_BOTHER){ clusterdataradump(GLOBAL_CDRA,NULL,0);}
      if (CAICOR_BOTHER){ caicordump(GLOBAL_CAICOR,NULL);}
      if (ISI_BOTHER){ isidump(GLOBAL_ISI,NULL);}
      if (HHLIB_BOTHER){ 
	if (verbose){ for (nt=0;nt<GLOBAL_HHLIB->Nra->ntypes;nt++){ hhlib_project_printf(GLOBAL_HHLIB,nt,0);}} 
	hhlibdump(GLOBAL_HHLIB,GLOBAL_HHLIB_READ_FILE);}
      if (SNXDATA_BOTHER){ snxdatadump(GLOBAL_SNXDATA);}
      printf("output dumped at time %f\n",GLOBAL_time);}
  if (GLOBAL_TF > 0 && GLOBAL_time > GLOBAL_TF){ cleanupoutput(); exit(EXIT_SUCCESS);}
}

/* Here is the main function */

int main(int argc, char **argv) 
{
  meminit(); /* initialize memory list for use */

/*   struct strobe **stra=NULL; */
/*   int ralength=4,length=18,window_length=1; */
/*   int nr=0,nt=0; */
/*   double *ra=NULL,*raproject=NULL; */
/*   stra = (struct strobe **) tcalloc(ralength,sizeof(struct strobe *)); */
/*   for (nr=0;nr<ralength;nr++){ */
/*     stra[nr]=strobemake(length,1,0,0,0,0); */
/*     for (nt=0;nt<length;nt++){ */
/*       stra[nr]->data[nt]=(nr+1)*nt+cos(4*PI*(nr+nt)/(double)ralength/(double)length)+sin(2*PI*(nr+nt)/(double)ralength/(double)length);}} */
/*   memprintf(0); */
/*   ra = (double *) tcalloc(ralength*length/window_length,sizeof(double)); */
/*   stra2ra(stra,ralength,window_length,ra); */
/*   memprintf(0); */
/*   raproject = (double *) tcalloc(ralength*length/window_length,sizeof(double)); */
/*   ra2pca(ra,length/window_length,ralength,raproject); */
/*   memprintf(0); */
/*   ON_MY_COMPUTER=1;pca2jpg(raproject,length/window_length,"./test"); */
/*   tfree(ra);ra=NULL; */
/*   tfree(raproject);raproject=NULL; */
/*   memprintf(0); */
/*   exit(0); */
  
/*   struct obsdisthist *odh1=NULL,*odh2=NULL; */
/*   FILE *fp=NULL; */
/*   memprintf(0); */
/*   odh1=obsdisthistmake(2); */
/*   obsdisthistupdate(odh1,0,11,1); obsdisthistupdate(odh1,0,13,1); obsdisthistupdate(odh1,0,15,1); obsdisthistupdate(odh1,0,17,1); */
/*   obsdisthistupdate(odh1,1,12,1); obsdisthistupdate(odh1,1,14,1); obsdisthistupdate(odh1,1,16,1); obsdisthistupdate(odh1,1,19,1); */
/*   obsdisthistappraise(odh1,4,0); */
/*   obsdisthistprintf(odh1,"odh1: "); */
/*   obsdisthistappraise(odh1,4,1); */
/*   obsdisthistprintf(odh1,"odh1: "); */
/*   fp=fopen("test.garbage","w"); */
/*   obsdisthist2file((void *)odh1,fp); */
/*   fclose(fp); */
/*   fp=fopen("test.garbage","r"); */
/*   odh2 = file2obsdisthist(fp); */
/*   obsdisthistprintf(odh2,"h2: "); */
/*   fclose(fp); */
/*   obsdisthisttfree(odh1);odh1=NULL; */
/*   obsdisthisttfree(odh2);odh2=NULL; */
/*   memprintf(0); */
/*   exit(0); */

/*   struct hist *h1=NULL,*h2=NULL; */
/*   FILE *fp=NULL; */
/*   memprintf(0); */
/*   h1=histmake(5,10,0); */
/*   histadd(h1,0.5,2);histadd(h1,1.5,1);histadd(h1,9.5,1);histadd(h1,5.5,10); */
/*   memprintf(0); */
/*   histprintf(h1,"h1: "); */
/*   fp=fopen("test.garbage","w"); */
/*   memprintf(0); */
/*   hist2file((void *)h1,fp); */
/*   memprintf(0); */
/*   fclose(fp); */
/*   fp=fopen("test.garbage","r"); */
/*   memprintf(0); */
/*   h2 = file2hist(fp); */
/*   histprintf(h2,"h2: "); */
/*   memprintf(0); */
/*   fclose(fp); */
/*   histtfree(h1);h1=NULL; */
/*   histtfree(h2);h2=NULL; */
/*   memprintf(0); */
/*   exit(0); */

  readinput(); /* get all the global variables from some input */
  setdependencies(); /* make sure a few things are kosher */    

/*   double ra1[6],ra2[6],*ra=NULL; */
/*   ra1[0]=1; ra1[2]=3; ra1[4]=5; */
/*   ra1[1]=2; ra1[3]=4; ra1[5]=6; */
/*   ra2[0]=-1; ra2[3]=4;  */
/*   ra2[1]=-2; ra2[4]=5;  */
/*   ra2[2]=-3; ra2[5]=6; */
/*   memprintf(0); */
/*   raprintf(ra1,"double",2,3,"ra1: "); raprintf(ra2,"double",3,2,"ra2: "); */
/*   memprintf(0); */
/*   ra = ra2ra_matrix_multiply(ra1,2,3,0,ra2,3,2,0); */
/*   raprintf(ra,"double",2,2,"ra1*ra2: "); tfree(ra);ra=NULL; */
/*   memprintf(0); */
/*   ra = ra2ra_matrix_multiply(ra1,2,3,1,ra2,3,2,1); */
/*   raprintf(ra,"double",3,3,"ra1'*ra2': ");tfree(ra);ra=NULL; */
/*   memprintf(0); */
/*   exit(0); */

/*   int length=3,depth=2,sdepth=1,pdepth=2; */
/*   int slength=(int)pow(2,length); */
/*   double *ra=NULL; */
/*   memprintf(0); */
/*   ra = binary_projection_s2s(length,depth,sdepth); */
/*   raprintf(ra,"double",pow(slength,sdepth),pow(slength,depth),"s2s: "); */
/*   tfree(ra);ra=NULL; */
/*   ra = binary_projection_s2p(0,length,depth,pdepth); */
/*   raprintf(ra,"double",pow(length,pdepth),pow(slength,depth),"s2p: "); */
/*   tfree(ra);ra=NULL; */
/*   memprintf(0); */
/*   exit(0); */
  
/*   int length=3,slength=0; */
/*   double *input=NULL,*synapse=NULL,*statetree=NULL,sum=0,*dspra=NULL; */
/*   slength=(int)pow(2,length); */
/*   input=(double *)tcalloc(length,sizeof(double)); */
/*   synapse=(double *)tcalloc(length*length,sizeof(double)); */
/*   input[0]=0.5;input[1]=0.5;input[2]=0.5; */
/*   synapse[0+0*length]=0.0;synapse[0+1*length]=0.1;synapse[0+2*length]=0.1; */
/*   synapse[1+0*length]=0.2;synapse[1+1*length]=0.0;synapse[1+2*length]=0.2; */
/*   synapse[2+0*length]=0.3;synapse[2+1*length]=0.3;synapse[2+2*length]=0.0; */
/*   statetree=(double *)tcalloc(slength*slength,sizeof(double)); */
/*   memprintf(0); */
/*   ideal_2_statetree(length,input,synapse,statetree); */
/*   raprintf(statetree,"double",slength,slength,"st: "); */
/*   stats("double",statetree,slength*slength,NULL,NULL,&sum,NULL);printf(" sum %f\n",sum*slength*slength); */
/*   memprintf(0); */
/*   dspra=ds_projection(slength); */
/*   exit(0); */

/*   int rows=5,cols=4,nrhs=8,nr=0,nc=0; */
/*   double *ra=NULL,*b=NULL,*x=NULL; */
/*   memprintf(0); */
/*   ra = (double *) tcalloc(rows*cols,sizeof(double)); */
/*   b = (double *) tcalloc(maximum(rows,cols)*nrhs,sizeof(double)); */
/*   x = (double *) tcalloc(maximum(rows,cols)*nrhs,sizeof(double)); */
/*   for (nr=0;nr<rows;nr++){ for (nc=0;nc<cols;nc++){ ra[nr+nc*rows]=sin(nr+nc+2)*(nc+2*nr+3)+cos(nr+nc+2)*(nr+2*nc+3);}} */
/*   for (nr=0;nr<rows;nr++){ for (nc=0;nc<nrhs;nc++){ b[nr+nc*rows]=(0.5+nr+1)*(nc+1);}} */
/*   ra2lss(ra,rows,cols,b,nrhs,x); */
/*   raprintf(ra,"double",rows,cols,"ra: "); */
/*   raprintf(b,"double",rows,nrhs,"b: "); */
/*   raprintf(x,"double",rows,nrhs,"x: "); */
/*   tfree(ra);ra=NULL; tfree(b);b=NULL; tfree(x);x=NULL;  */
/*   memprintf(0); */
/*   exit(0); */

/*   int rows=5,cols=5,nr=0,nc=0; */
/*   double *ra=NULL,*sv=NULL,*svl=NULL,*svr=NULL; */
/*   memprintf(0); */
/*   ra = (double *) tcalloc(rows*cols,sizeof(double)); */
/*   sv = (double *) tcalloc(minimum(rows,cols),sizeof(double)); */
/*   svl = (double *) tcalloc(rows*rows,sizeof(double)); */
/*   svr = (double *) tcalloc(cols*cols,sizeof(double)); */
/*   for (nr=0;nr<rows;nr++){ for (nc=0;nc<cols;nc++){ ra[nr+nc*rows]=sin(nr+nc+2)*(nc+2*nr+3)+cos(nr+nc+2)*(nr+2*nc+3);}} */
/*   ra2svd(ra,rows,cols,sv,svl,svr); */
/*   raprintf(sv,"double",1,minimum(rows,cols),"sv: "); */
/*   raprintf(svl,"double",rows,rows,"svl: "); */
/*   raprintf(svr,"double",cols,cols,"svr: "); */
/*   tfree(ra);ra=NULL; tfree(sv);sv=NULL; tfree(svl);svl=NULL; tfree(svr);svr=NULL; */
/*   memprintf(0); */
/*   exit(0); */

/*   int length=3; */
/*   double *ra=NULL,*eiwr=NULL,*eiwi=NULL,*eig=NULL; */
/*   ra = (double *) tcalloc(length*length,sizeof(double)); */
/*   eiwr = (double *) tcalloc(length,sizeof(double)); */
/*   eiwi = (double *) tcalloc(length,sizeof(double)); */
/*   eig = (double *) tcalloc(length*length,sizeof(double)); */
/*   ra[0 + 0*length] = 0.3; ra[0 + 1*length] = 0.2; ra[0 + 2*length] = 0.1; */
/*   ra[1 + 0*length] = 0.3; ra[1 + 1*length] = 0.5; ra[1 + 2*length] = 0.2; */
/*   ra[2 + 0*length] = 0.4; ra[2 + 1*length] = 0.3; ra[2 + 2*length] = 0.7; */
/*   ra2eig(ra,length,eiwr,eiwi,eig); */
/*   raprintf(eiwr,"double",1,length,"eiwr: "); */
/*   raprintf(eig,"double",length,length,"eig: "); */
/*   tfree(ra);ra=NULL; tfree(eiwr);eiwr=NULL; tfree(eiwi);eiwi=NULL; tfree(eig);eig=NULL; */
/*   exit(0); */

/*   int npoints=256; */
/*   GLOBAL_PTREE_REGION_TYPE=0; */
/*   synaptic_sphere_prediction_plot(npoints); */
/*   exit(0); */  

/*   int ndim=2,npoints=10000; */
/*   nsphere_index_generate(ndim,npoints); */
/*   nsphere_index_read(ndim,npoints,-1); */
/*   exit(0); */

/*   GLOBAL_PTREE_REGION_TYPE=0; */
/* /\*   granule_int2input_or_input2int(2,GLOBAL_ODORra_START,32,NULL,32,NULL); *\/ */
/* /\*   granule_local_plot(); *\/ */
/*   granule_plot(32,32); */
/* /\*   struct ptree *pW=NULL,*pW2=NULL; *\/ */
/* /\*   pW = granule_obtain_reweight("./ptree_i23j15_granule0record","./ptree_i22j16_granule0record",SUITE_DUMPEVERY); *\/ */
/* /\*   pnodeprintf(NULL,pW->postree,-1,0); *\/ */
/* /\*   ptreedump_starter(pW,"./ptree_reweight_23_15_22_16_granule0record",2,1,0,0,0.2,-0.2,1); *\/ */
/* /\*   pW2 = ptreadback("./ptree_reweight_23_15_22_16_granule0record",0); *\/ */
/* /\*   pnodeprintf(NULL,pW2->postree,-1,0); *\/ */
/*   exit(0); */

/*   char **fnamebase=NULL; */
/*   int nr1=0; */
/*   char *gs2=GLOBAL_STRING_2; */
/*   struct litem *l0=NULL; */
/*   struct pnode *pn=NULL; */
/*   struct obsdisthist *odh=NULL; */
/*   GLOBAL_PTREE_REGION_TYPE = 0; */
/*   GLOBAL_PTREE=ptree_readbacktemp("/home/rangan/Cstuff/v1/darena/dir_natfig2_results/ptree_dumptemp_ifengine_0_0_00record_time1048576",&file2obsdisthist); */
/*   //pnode_obs2dist_starter(0,0,(1024*SUITE_NSECONDS)/SUITE_DUMPEVERY,NULL,GLOBAL_PTREE->postree,-1,0,3); */
/* /\*   pnodeprintf(NULL,GLOBAL_PTREE->postree,0,0); *\/ */
/* /\*   pn = (struct pnode *)(GLOBAL_PTREE->postree->kidr->kidl->item); *\/ */
/* /\*   printf("\n\n%d, %d\n\n",(int)pn,pn->region->label); *\/ */
/* /\*   pnodeprintf(pn,pn->childllitem,0,0); *\/ */
/* /\*   pn = (struct pnode *)(pn->childllitem->kidr->kidr->kidr->kidr->kidr->kidr->kidr->item); *\/ */
/* /\*   printf("\n\n%d, %d\n\n",(int)pn,pn->region->label); *\/ */
/* /\*   pnodeprintf(pn,pn->childllitem,0,0); *\/ */
/* /\*   pn = (struct pnode *)(pn->childllitem->kidr->kidr->kidr->kidr->kidr->kidr->item); *\/ */
/* /\*   printf("\n\n%d, %d\n\n",(int)pn,pn->region->label); *\/ */
/* /\*   pnodeprintf(pn,pn->childllitem,0,0); *\/ */
/* /\*   pn = (struct pnode *)(pn->childllitem->item); *\/ */
/* /\*   printf("\n\n%d, %d\n\n",(int)pn,pn->region->label); *\/ */
/* /\*   pnodeprintf(pn,pn->childllitem,0,0); *\/ */
/* /\*   odh = (struct obsdisthist *) pn->temp; *\/ */
/* /\*   histdump(odh->hra[1],0,"ifengine_0_0_10record_odor1_hist_2","odor1",0); *\/ */
/* /\*   histdump(odh->hra[3],0,"ifengine_0_0_10record_odor3_hist_2","odor3",0); *\/ */
/*   pnode_obsdisthist_tofig_starter(GLOBAL_PTREE,NULL,4); */
/* /\*   fnamebase = (char **) tcalloc(SUITE_NODORS,sizeof(char *)); *\/ */
/* /\*   for (nr1=0;nr1<SUITE_NODORS;nr1++){ *\/ */
/* /\*     fnamebase[nr1] = (char *) tcalloc(256,sizeof(char)); *\/ */
/* /\*     sprintf(fnamebase[nr1],"ptree_odor%d_%srecord",nr1,gs2);} *\/ */
/* /\*   ptree_classify_starter(SUITE_NODORS,GLOBAL_PTREE_NLEGS+1,2,fnamebase,"temporary",0); *\/ */
/*   exit(0); */

  setglobals(); /* initialize global variables */

/*   char filename[1024]; */
/*   odorfprintf_full(stdout,GLOBAL_ODORra); */
/*   sprintf(filename,"odor_full_fwrite_test"); */
/*   odor_full_fwrite(filename,GLOBAL_ODORra); */
/*   struct odor *o=NULL; */
/*   o=odor_full_fread(filename); */
/*   odortfree(o); */
/*   exit(0); */
  
/*   ptree_test_mcpit();exit(0); */
/*   suite_7_power_process_helper("wilson_3_0_record");exit(0); */
/*   suite_8_power_process_helper("wilson5_0_record");exit(0); */

/* GLOBAL_verbose=3;  suite_7_power_process_helper("wilson6b_victor_0_0_record");exit(0); // I am not sure exactly what this does (Adi: 20230928) */

/*   int rows=4,cols=5; */
/*   int nr=0,nc=0; */
/*   double *ra=NULL,*ra2=NULL; */
/*   ra = (double *) tcalloc(rows*cols,sizeof(double)); */
/*   for (nr=0;nr<rows;nr++){ for (nc=0;nc<cols;nc++){ ra[nr+nc*rows]=nr+nc*rows;}} */
/*   raprintf(ra,"double",rows,cols,"ra: "); */
/*   ra2 = raplugout(ra,rows,cols,1,1,2,3); */
/*   raprintf(ra2,"double",2,3,"ra2: "); */
/*   exit(0); */

/*   char gs2[1024]; */
/*   sprintf(gs2,"%srecord",GLOBAL_STRING_2); */
/*   suite_8_power_process_helper2(gs2,2048,512); */
/*   exit(0); */

#ifndef NOTGRAPH
  if (SUPERGLOBAL_DRAW_FLAG){
    glutInit(&argc, argv);  /* initialize glut */
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_ALPHA | GLUT_DEPTH); /* select display mode, 2xbuffer,rgba,alpha,depth_buffer */
    glutInitWindowSize(GLOBAL_WINDOW_WIDTH, GLOBAL_WINDOW_HEIGHT); /* small window */
    glutInitWindowPosition(0,0); /* starts at upper left corner of screen */
    GLUTWINDOWNUMBER = glutCreateWindow("al"); /* open window */
    glutDisplayFunc(&DrawGLScene); /* this does all the drawing for us */
    glutFullScreen(); GLOBAL_WINDOW_WIDTH=1024; GLOBAL_WINDOW_HEIGHT=768; /* if we wish to fullscreen... */
    glutIdleFunc(&DrawGLScene); /* continuously redraw */
    glutReshapeFunc(&ReSizeGLScene); /* resizing function */
    glutKeyboardFunc(&keyPressed);glutSpecialFunc(&specialKeyPressed); /* input functions */
    InitGL(GLOBAL_WINDOW_WIDTH,GLOBAL_WINDOW_HEIGHT); /* init window */
    glutMainLoop(); /* start processing */}
  else{ while(!RUN_DONE){ computestep();}}
#endif /* NOTGRAPH */
#ifdef NOTGRAPH
  while (!RUN_DONE){ computestep();}
#endif /* NOTGRAPH */
  return 1;
}
