/* no retinotopy (yet), lgn through to v1, trying for minimalism */
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

#include "d_inputhandling.h"
#include "matrices.h"
#include "d_llists.h"
#include "d_lgn.h"
#include "d_cortex.h"
#include "d_datahandling.h"
#include "d_ogletools.h"

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
int GLOBAL_CLEANUP=1;
int NSLICES=0;
int NPHASES=0;
int CORTEX_BOTHER=0;
int AUTAPSES_OFF=1;
int LR_BOTHER=1;
int EIFvsIF=0;
int LGN_DETERMINISM=0;
int HANSHELLEY_FLAG=0;
double STARTING_HOMOGENIZATION=1;
int LGN_BOTHER=0;
int SUITE_BOTHER=0;
int SUITE_NSECONDS=64;
int SUITE_DUMPEVERY=64;
int SUITE_BITBYBIT_REMOVE=1;
int SUITE_SINDEXMAX=0;
int SUITE_DINDEXMAX=0;
int SUITE_TINDEXMAX=1;
int LGN_DUMP=0;
int LGN_TYPE_FLAG=0;
int LGN_DUMB=0;
int GLOBAL_ODOR=0;
int GLOBAL_ODOR_BASE=2;
int GLOBAL_ODOR_BACON=1;
int GLOBAL_FPS=0;
int GLOBAL_SPF=0;
int LGN_CELLS_PER_V1_CELL=0;
struct lgn *GLOBAL_LGN=NULL;
int ORDER=0;
int OBSCURED=0;
int DEPRESSION=0;
int SPIKETOL=0;
int SYSTEM_ROW_SIZE=0;
int SYSTEM_COL_SIZE=0;
int ARBOR_DIA=0;
int NARBORS_TALL=0;
int NARBORS_WIDE=0;
int PIE_ROW_DIA=0;
int PIE_COL_DIA=0;
int NPIEROWS=0;
int NPIECOLS=0;
double LR_DIST_DIA_INPIES=0;
double LR_ANGLE_DIA=0;
double LR_DIR_DIA=0;
double LR_RATE=0;
double LR_TO_AMPA=0;
double AXONAL_DELAY_PER_PIE;
double P_ES=0;
double P_EC=0;
double P_IS=0;
double P_IC=0;
double P_AMPA=1;
double P_NMDA=1;
double P_GABA=1;
double GABOR_DRIFT=0;
double LGNANGLE_DRIFT=0;
double AMPA_DIA=0;
double NMDA_DIA=0;
double GABA_DIA=0;
int GRATING_VS_LMI=0;
double CYCLE_LENGTH=0;
int INPUT_IGNORE_FLAG=0;
double GRATING_PULSE=0;
double GRATING_DRIFT=0;
double STIMULUS_ONSET_TIME=0;
double SQUARE_DURATION_TIME=0;
double SQUARE_DRAG_TIME=0;
double LINE_DELAY_TIME=0;
double OTHERLAYER_BACKRATE=0;
double OTHERLAYER_INPUTRATE=0;
double OTHERLAYER_STRENGTH=0;
double LGN_BACKRATE=0;
double LGN_STRENGTH=0;
double INPUT_CONTRAST=0;
double INPUT_SPACEK=0;
double INPUT_SPACEANGLE=0;
double INPUT_SPACEANGLE_BACON=0.1047; /* six degrees of separation */
double INPUT_SPACEPHASE=0;
double CS_ESES_A=0;
double CS_ESES_N=0;
double CS_ESIS_A=0;
double CS_ESIS_N=0;
double CS_ESEC_A=0;
double CS_ESEC_N=0;
double CS_ESIC_A=0;
double CS_ESIC_N=0;
double CS_ECES_A=0;
double CS_ECES_N=0;
double CS_ECIS_A=0;
double CS_ECIS_N=0;
double CS_ECEC_A=0;
double CS_ECEC_N=0;
double CS_ECIC_A=0;
double CS_ECIC_N=0;
double CS_ISES_G=0;
double CS_ISIS_G=0;
double CS_ISEC_G=0;
double CS_ISIC_G=0;
double CS_ICES_G=0;
double CS_ICIS_G=0;
double CS_ICEC_G=0;
double CS_ICIC_G=0;
double CS_ESLR=0;
double CS_ISLR=0;
double CS_ECLR=0;
double CS_ICLR=0;
double CS_LRES=0;
double CS_LRIS=0;
double CS_LREC=0;
double CS_LRIC=0;
int SPARSE_ESES=0;
int SPARSE_ESEC=0;
int SPARSE_ESIS=0;
int SPARSE_ESIC=0;
int SPARSE_ECES=0;
int SPARSE_ECEC=0;
int SPARSE_ECIS=0;
int SPARSE_ECIC=0;
int SPARSE_ISES=0;
int SPARSE_ISEC=0;
int SPARSE_ISIS=0;
int SPARSE_ISIC=0;
int SPARSE_ICES=0;
int SPARSE_ICEC=0;
int SPARSE_ICIS=0;
int SPARSE_ICIC=0;
int SPARSEHIT_ESES=0;
int SPARSEHIT_ESEC=0;
int SPARSEHIT_ESIS=0;
int SPARSEHIT_ESIC=0;
int SPARSEHIT_ECES=0;
int SPARSEHIT_ECEC=0;
int SPARSEHIT_ECIS=0;
int SPARSEHIT_ECIC=0;
int SPARSEHIT_ISES=0;
int SPARSEHIT_ISEC=0;
int SPARSEHIT_ISIS=0;
int SPARSEHIT_ISIC=0;
int SPARSEHIT_ICES=0;
int SPARSEHIT_ICEC=0;
int SPARSEHIT_ICIS=0;
int SPARSEHIT_ICIC=0;
double LR_PCONNECT=0;
double TAU_AMPA=0;
double TAU_NMDA=0;
double TAU_GABA=0;
double DEPRESS_pA=0;
double DEPRESS_pN=0;
double DEPRESS_pG=0;
double DEPRESS_pLR=0;
double TAU_DEPRESS_AMPA=0;
double TAU_DEPRESS_NMDA=0;
double TAU_DEPRESS_GABA=0;
double TAU_DEPRESS_LR=0;
double TAU_REF=0;
double VOLTAGE_REST=0;
double VOLTAGE_RESET=0;
double VOLTAGE_THRESHOLD=0;
double VOLTAGE_EX=0;
double VOLTAGE_IN=0;
double VOLTAGE_THRESHOLD_EIF=0;
double VOLTAGE_TAKEOFF=0;
double VOLTAGE_DELTAT=0;
double CONDUCTANCE_LK=0;
double CONDUCTANCE_EX_MAX=0;
double CONDUCTANCE_IN_MAX=0;
int GLOBAL_verbose=0,GLOBAL_dtok=0,GLOBAL_dtadapt=0;
double GLOBAL_time=0;
double GLOBAL_TI=0;
double GLOBAL_TF=0;
double GLOBAL_DT=0;
double GLOBAL_DTmax=0;
struct neuronarray *GLOBAL_Nra=NULL;
/* stored as [shuffled_ang_index + pie_row*PIE_ROW_DIA*PIE_COL_DIA + pie_col*PIE_ROW_DIA*PIE_COL_DIA*NPIEROWS] */
int *GLOBAL_SHUFFLES=NULL;
int NLRKERNELS=1;
/* stored as [pie_row + pie_col*NPIEROWS + ang_index*NPIEROWS*NPIECOLS] */
/* remember, row major order means that the `index' is (ang_index,pie_col,pie_row) */
fftw_complex **GLOBAL_LRKERNEL_X=NULL;
fftw_complex **GLOBAL_LRKERNEL_K=NULL;
fftw_complex **GLOBAL_LRSWAP_X=NULL;
fftw_complex **GLOBAL_LRSWAP_K=NULL;
fftw_plan *GLOBAL_FFTWPLAN_KERNEL_FORWARD=NULL;
fftw_plan *GLOBAL_FFTWPLAN_KERNEL_BACKWARD=NULL;
fftw_plan *GLOBAL_FFTWPLAN_SWAP_FORWARD=NULL;
fftw_plan *GLOBAL_FFTWPLAN_SWAP_BACKWARD=NULL;
int GLOBAL_BLOCK_ROW_DIA=0;
int GLOBAL_BLOCK_COL_DIA=0;
int GLOBAL_BLOCK_ROW_DIA_MIN=0;
int GLOBAL_BLOCK_COL_DIA_MIN=0;
int GLOBAL_NBLOCKS_TALL=0;
int GLOBAL_NBLOCKS_WIDE=0;
int FIDDLE_PARAMETER=0;
int FIDDLE_ROW=0;
int FIDDLE_COL=0;
int DRAW_FLAG=0;
int DRAW_FLAG2=0;
int STEPS_PER_DRAW=0;
double STD_VIEW=1;
double PTREE_VIEW=0;
int GRAYSCALE=0;
int GLOBAL_SPACE_SMOOTHER=0;
int SUPERGLOBAL_DRAW_FLAG=0;
int BIG_SYSTEM_FLAG=0;
double OUTPUT_DUMP_EVERY=0;
int MOVIE_BOTHER=0;
int FRAMEINDEX=0;
double MOVIE_LAST_TIME=0;
int PNM_BOTHER=0;
int PNM_REMOVE=0;
int PNM_FRAMEINDEX=0;
double PNM_MOVIE_LAST_TIME=0;
double MOVIE_START_RECORDING=0;
double MOVIE_TIME_REFRESH=0;
int RTC_BOTHER=0;
struct rtc *GLOBAL_RTC=NULL;
int GLOBAL_RTC_LENGTH=0;
int GLOBAL_RTC_FRAMELENGTH=0;
int GLOBAL_RTC_NANGLES=0;
int GLOBAL_RTC_NPHASES=0;
int STROBETRACE_BOTHER=0;
struct strobetrace *GLOBAL_STROBETRACE=NULL;
int GLOBAL_STROBETRACE_NANGLES=0;
int GLOBAL_STROBETRACE_CYCLE_BOTHER=0;
double GLOBAL_STROBETRACE_LENGTH=0;
double GLOBAL_STROBETRACE_AVALANCHE_UPDATE_EVERY=0;
int GLOBAL_STROBETRACE_AVALANCHE_LOGDMIN=0;
int GLOBAL_STROBETRACE_AVALANCHE_LOGDMAX=0;
int TUNINGCURVE_BOTHER=0;
struct tuningcurve *GLOBAL_TUNINGCURVE=NULL;
int GLOBAL_TUNINGCURVE_NANGLES=0;
int GLOBAL_TUNINGCURVE_NRADIUS=0;
int LMITRI_BOTHER=0;
struct lmitri *GLOBAL_LMITRI=NULL;
int GLOBAL_LMITRI_T0;
int GLOBAL_LMITRI_TIMELENGTH;
double GLOBAL_LMITRI_ROW_MAX;
double GLOBAL_LMITRI_ROW_MIN;
int PTREE_BOTHER=0;
int GLOBAL_PTREE_BITBYBIT=0;
int GLOBAL_PTREE_ZZZ=1;
struct ptree *GLOBAL_PTREE=NULL;
int GLOBAL_PTREE_NREGIONS=8;
int GLOBAL_PTREE_NLEGS=0;
int GLOBAL_PTREE_LEGTIME=0;
double P_REWIRE=0;
int SPARSE_ROW_DIA=4;
int SPARSE_COL_DIA=1;
int SPARSE_T2S_ORDERING=-1; /* non-negative values force t2s ordering */
int SPARSE_CONNECTION_TYPE=-1; /* negative values correspond to specific architectures as long as BIG_SYSTEM_FLAG<0 */
double SPARSE_CS_SCALE=1;
int GLOBAL_PTREE_EVENT_THRESHOLD=2;
double GLOBAL_PTREE_EVENT_WITHIN=5;
int GLOBAL_PTREE_REGION_TYPE=2;
int CLOSET_BOTHER=0;
struct closet *GLOBAL_CLOSET=NULL;
int YGGDRASIL_BOTHER=0;
int GLOBAL_YGGDRASIL_PPNREGIONS=0;
int GLOBAL_YGGDRASIL_PPNLEGS=0;
int GLOBAL_YGGDRASIL_PPLEGTIME=0;
int GLOBAL_YGGDRASIL_WEIGHT_MINIMUM=0;
struct yggdrasil *GLOBAL_YGGDRASIL=NULL;
int BONSAI_BOTHER=0;
struct bonsai *GLOBAL_BONSAI=NULL;
int HYDRA_BOTHER=0;
struct hydra *GLOBAL_HYDRA=NULL;
int GLOBAL_HYDRA_SWITCH_EVERY_TWO=0;
int GLOBAL_HYDRA_JUSTONTIME=0;
int GLOBAL_HYDRA_STAYONTIME=0;
int GLOBAL_HYDRA_DUMP_EVERY=0;
int LYAPUNOV_BOTHER=0;
struct lyapunov *GLOBAL_LYAPUNOV=NULL;
double GLOBAL_LYAPUNOV_UPDATE_EVERY=0;
double GLOBAL_LYAPUNOV_JIGGLE=0;
int POWER_BOTHER=0;
struct power *GLOBAL_POWER=NULL;
int GLOBAL_POWER_LENGTH=0;
int GLOBAL_POWER_HOWMANY=0;
int TAOF_BOTHER=0;
struct taof *GLOBAL_TAOF=NULL;
int GLOBAL_TAOF_LENGTH=0;
int GLOBAL_TAOF_STEP_EVERY=0;
int SEIDCORR_BOTHER=0;
struct seidcorr *GLOBAL_SEIDCORR=NULL;
int GLOBAL_SEIDCORR_SPACE_BIN_SIZE=0;
int GLOBAL_SEIDCORR_TIME_BIN_SIZE=0;
int GLOBAL_SEIDCORR_LENGTH=0;
double GLOBAL_SEIDCORR_TIME_START=0;

#include "d_inputhandling.c"
#include "matrices.c"
#include "d_llists.c"
#include "d_lgn.c"
#include "d_cortex.c"
#include "d_datahandling.c"
#include "d_ogletools.c"

/* Here are the method functions */

void setglobals()
{
  int verbose=GLOBAL_verbose;
  char text[64],command[256];
  int nr=0,nc=0,nl=0;
  FILE *fp=NULL;
#ifndef NOTGRAPH
  GLOBAL_WINDOW_WIDTH = 640;
  GLOBAL_WINDOW_HEIGHT = 480;
  xdepth = -2.2f;
  ydepth = -1.2f;
  zdepth = -6.0f;
#endif /* NOTGRAPH */
  GLOBAL_time=GLOBAL_TI;
  GLOBAL_DT=GLOBAL_DTmax;
  GLOBAL_dtok=0;
  GLOBAL_SPF=1024/GLOBAL_FPS;
  MOVIE_LAST_TIME=GLOBAL_TI;
  PNM_MOVIE_LAST_TIME=GLOBAL_TI;
  LGN_CELLS_PER_V1_CELL=16;
  if (verbose){ printf(" %% setting LGN_CELLS_PER_V1_CELL=%d\n",LGN_CELLS_PER_V1_CELL);}
  if (LGN_BOTHER){ 
    if (verbose){ printf(" %% trying to make LGN\n");}
    GLOBAL_LGN = lgnmake("blank_64x96",64,96,16.0,24.0,0.75,LGN_BACKRATE,3*LGN_BACKRATE);}
  if (CORTEX_BOTHER){
    if (verbose){ printf(" %% BIG_SYSTEM_FLAG=%d\n",BIG_SYSTEM_FLAG);}
    switch (BIG_SYSTEM_FLAG){ /* line, small, lmi or square */
    case -1: /* line */ NARBORS_WIDE=1; NARBORS_TALL=2; break;
    case 0: /* small */ NARBORS_WIDE=12; NARBORS_TALL=8; break;
    case 1: /* lmi */ NARBORS_WIDE=24; NARBORS_TALL=16; break;
    default: /* square */ NARBORS_WIDE=BIG_SYSTEM_FLAG; NARBORS_TALL=BIG_SYSTEM_FLAG; break;}
    if (BIG_SYSTEM_FLAG >= 0){
      if (verbose){ printf(" %% making big system\n");}    
      assert(minimum(NARBORS_WIDE,NARBORS_TALL) > 2);
      PIE_ROW_DIA=ARBOR_DIA*2;
      PIE_COL_DIA = ARBOR_DIA*2;
      NPIEROWS=NARBORS_TALL/2;
      NPIECOLS=NARBORS_WIDE/2;
      assert(PIE_ROW_DIA*NPIEROWS==ARBOR_DIA*NARBORS_TALL);
      assert(PIE_COL_DIA*NPIECOLS==ARBOR_DIA*NARBORS_WIDE);
      if ((NPIEROWS%2!=0 || NPIECOLS%2!=0) && GLOBAL_verbose){ printf(" %% warning, pies not correct, and longrange connections suspect\n");}
      SYSTEM_ROW_SIZE=NARBORS_TALL*ARBOR_DIA;
      SYSTEM_COL_SIZE=NARBORS_WIDE*ARBOR_DIA;
      AMPA_DIA=1.0*ARBOR_DIA;
      NMDA_DIA=0.5*ARBOR_DIA;
      GABA_DIA=1.0*ARBOR_DIA;
      if (HANSHELLEY_FLAG){ GLOBAL_BLOCK_ROW_DIA_MIN=SYSTEM_ROW_SIZE; GLOBAL_BLOCK_COL_DIA_MIN=SYSTEM_COL_SIZE;} 
      else{ GLOBAL_BLOCK_ROW_DIA_MIN=ARBOR_DIA; GLOBAL_BLOCK_COL_DIA_MIN=ARBOR_DIA;}
      GLOBAL_BLOCK_ROW_DIA=GLOBAL_BLOCK_ROW_DIA_MIN;
      GLOBAL_BLOCK_COL_DIA=GLOBAL_BLOCK_COL_DIA_MIN;
      GLOBAL_NBLOCKS_TALL=SYSTEM_ROW_SIZE/GLOBAL_BLOCK_ROW_DIA;
      GLOBAL_NBLOCKS_WIDE=SYSTEM_COL_SIZE/GLOBAL_BLOCK_COL_DIA;}
    else if (BIG_SYSTEM_FLAG<0){
      if (verbose){ printf(" %% making small system\n");}
      LR_BOTHER=0;
      PIE_ROW_DIA=ARBOR_DIA*2;
      PIE_COL_DIA = 1;
      NPIEROWS=NARBORS_TALL/2;
      NPIECOLS=1;
      assert(PIE_ROW_DIA*NPIEROWS==ARBOR_DIA*NARBORS_TALL);
      assert(PIE_COL_DIA*NPIECOLS==NARBORS_WIDE);
      SYSTEM_ROW_SIZE=NARBORS_TALL*ARBOR_DIA;
      SYSTEM_COL_SIZE=NARBORS_WIDE;
      AMPA_DIA=1.0*ARBOR_DIA;
      NMDA_DIA=0.5*ARBOR_DIA;
      GABA_DIA=1.0*ARBOR_DIA;
      if (HANSHELLEY_FLAG){ GLOBAL_BLOCK_ROW_DIA_MIN=SYSTEM_ROW_SIZE; GLOBAL_BLOCK_COL_DIA_MIN=SYSTEM_COL_SIZE;} 
      else{ GLOBAL_BLOCK_ROW_DIA_MIN=ARBOR_DIA; GLOBAL_BLOCK_COL_DIA_MIN=SYSTEM_COL_SIZE;}
      GLOBAL_BLOCK_ROW_DIA=GLOBAL_BLOCK_ROW_DIA_MIN;
      GLOBAL_BLOCK_COL_DIA=GLOBAL_BLOCK_COL_DIA_MIN;
      GLOBAL_NBLOCKS_TALL=SYSTEM_ROW_SIZE/GLOBAL_BLOCK_ROW_DIA;
      GLOBAL_NBLOCKS_WIDE=SYSTEM_COL_SIZE/GLOBAL_BLOCK_COL_DIA;}
    if (verbose){ printf(" %% setting up sparsity structure\n");}
    SPARSEHIT_ICIC = maximum(1,(int)floor(PI*pow(GABA_DIA,2)*P_IC/(P_ES+P_EC+P_IS+P_IC)));
    SPARSEHIT_ICIS = maximum(1,(int)floor(PI*pow(GABA_DIA,2)*P_IS/(P_ES+P_EC+P_IS+P_IC)));
    SPARSEHIT_ICEC = maximum(1,(int)floor(PI*pow(GABA_DIA,2)*P_EC/(P_ES+P_EC+P_IS+P_IC)));
    SPARSEHIT_ICES = maximum(1,(int)floor(PI*pow(GABA_DIA,2)*P_ES/(P_ES+P_EC+P_IS+P_IC)));
    SPARSEHIT_ISIC = maximum(1,(int)floor(PI*pow(GABA_DIA,2)*P_IC/(P_ES+P_EC+P_IS+P_IC)));
    SPARSEHIT_ISIS = maximum(1,(int)floor(PI*pow(GABA_DIA,2)*P_IS/(P_ES+P_EC+P_IS+P_IC)));
    SPARSEHIT_ISEC = maximum(1,(int)floor(PI*pow(GABA_DIA,2)*P_EC/(P_ES+P_EC+P_IS+P_IC)));
    SPARSEHIT_ISES = maximum(1,(int)floor(PI*pow(GABA_DIA,2)*P_ES/(P_ES+P_EC+P_IS+P_IC)));
    SPARSEHIT_ECIC = maximum(1,(int)floor(PI*pow(AMPA_DIA,2)*P_IC/(P_ES+P_EC+P_IS+P_IC)));
    SPARSEHIT_ECIS = maximum(1,(int)floor(PI*pow(AMPA_DIA,2)*P_IS/(P_ES+P_EC+P_IS+P_IC)));
    SPARSEHIT_ECEC = maximum(1,(int)floor(PI*pow(AMPA_DIA,2)*P_EC/(P_ES+P_EC+P_IS+P_IC)));
    SPARSEHIT_ECES = maximum(1,(int)floor(PI*pow(AMPA_DIA,2)*P_ES/(P_ES+P_EC+P_IS+P_IC)));
    SPARSEHIT_ESIC = maximum(1,(int)floor(PI*pow(AMPA_DIA,2)*P_IC/(P_ES+P_EC+P_IS+P_IC)));
    SPARSEHIT_ESIS = maximum(1,(int)floor(PI*pow(AMPA_DIA,2)*P_IS/(P_ES+P_EC+P_IS+P_IC)));
    SPARSEHIT_ESEC = maximum(1,(int)floor(PI*pow(AMPA_DIA,2)*P_EC/(P_ES+P_EC+P_IS+P_IC)));
    SPARSEHIT_ESES = maximum(1,(int)floor(PI*pow(AMPA_DIA,2)*P_ES/(P_ES+P_EC+P_IS+P_IC)));
    if (SPARSE_CS_SCALE!=1){ CS_ESES_A*=SPARSE_CS_SCALE;CS_ESIS_A*=SPARSE_CS_SCALE;CS_ESEC_A*=SPARSE_CS_SCALE;CS_ESIC_A*=SPARSE_CS_SCALE;CS_ECES_A*=SPARSE_CS_SCALE;CS_ECIS_A*=SPARSE_CS_SCALE;CS_ECEC_A*=SPARSE_CS_SCALE;CS_ECIC_A*=SPARSE_CS_SCALE;CS_ESES_N*=SPARSE_CS_SCALE;CS_ESIS_N*=SPARSE_CS_SCALE;CS_ESEC_N*=SPARSE_CS_SCALE;CS_ESIC_N*=SPARSE_CS_SCALE;CS_ECES_N*=SPARSE_CS_SCALE;CS_ECIS_N*=SPARSE_CS_SCALE;CS_ECEC_N*=SPARSE_CS_SCALE;CS_ECIC_N*=SPARSE_CS_SCALE;CS_ISES_G*=SPARSE_CS_SCALE;CS_ISIS_G*=SPARSE_CS_SCALE;CS_ISEC_G*=SPARSE_CS_SCALE;CS_ISIC_G*=SPARSE_CS_SCALE;CS_ICES_G*=SPARSE_CS_SCALE;CS_ICIS_G*=SPARSE_CS_SCALE;CS_ICEC_G*=SPARSE_CS_SCALE;CS_ICIC_G*=SPARSE_CS_SCALE;CS_LRES*=SPARSE_CS_SCALE;CS_LRIS*=SPARSE_CS_SCALE;CS_LREC*=SPARSE_CS_SCALE;CS_LRIC*=SPARSE_CS_SCALE;}
    DEPRESS_pLR=(LR_TO_AMPA)*DEPRESS_pA + (1-LR_TO_AMPA)*DEPRESS_pN;
    TAU_DEPRESS_LR=1.0/(LR_TO_AMPA*1.0/TAU_DEPRESS_AMPA + (1-LR_TO_AMPA)*1.0/TAU_DEPRESS_LR);
    CONDUCTANCE_EX_MAX=CONDUCTANCE_LK;
    CONDUCTANCE_IN_MAX=CONDUCTANCE_LK;
    GLOBAL_dtok=0;
    if (verbose){ printf(" %% making neuronarray with size %d,%d\n",SYSTEM_ROW_SIZE,SYSTEM_COL_SIZE);}
    GLOBAL_Nra=NULL;
    GLOBAL_Nra = neuronarraymake(SYSTEM_ROW_SIZE,SYSTEM_COL_SIZE);
    if (verbose){ printf(" %% making GLOBAL_SHUFFLES\n");}
    GLOBAL_SHUFFLES = (int *) tcalloc(NPIEROWS*NPIECOLS*PIE_ROW_DIA*PIE_COL_DIA,sizeof(int));
    for (nr=0;nr<NPIEROWS;nr++){ for (nc=0;nc<NPIECOLS;nc++){ makeshuffle(GLOBAL_SHUFFLES+nr*PIE_ROW_DIA*PIE_COL_DIA+nc*PIE_ROW_DIA*PIE_COL_DIA*NPIEROWS,PIE_ROW_DIA*PIE_COL_DIA,(int)floor(PIE_ROW_DIA*PIE_COL_DIA*LR_ANGLE_DIA/2/PI/2),nr+nc*NPIEROWS+GLOBAL_RECORD_NUMBER);}}
    if (verbose){ printf(" %% making GLOBAL_LRKERNEL\n");}
    GLOBAL_LRKERNEL_X = (fftw_complex **) tcalloc(NLRKERNELS,sizeof(fftw_complex *));
    GLOBAL_LRKERNEL_K = (fftw_complex **) tcalloc(NLRKERNELS,sizeof(fftw_complex *));
    GLOBAL_LRSWAP_X = (fftw_complex **) tcalloc(NLRKERNELS,sizeof(fftw_complex *));
    GLOBAL_LRSWAP_K = (fftw_complex **) tcalloc(NLRKERNELS,sizeof(fftw_complex *));
    for (nl=0;nl<NLRKERNELS;nl++){
      GLOBAL_LRKERNEL_X[nl] = makelrkernel(nl+GLOBAL_RECORD_NUMBER);
      GLOBAL_LRKERNEL_K[nl] = (fftw_complex *) tcalloc(NPIEROWS*NPIECOLS*PIE_ROW_DIA*PIE_COL_DIA*NLRKERNELS,sizeof(fftw_complex));
      GLOBAL_LRSWAP_X[nl] = (fftw_complex *) tcalloc(NPIEROWS*NPIECOLS*PIE_ROW_DIA*PIE_COL_DIA,sizeof(fftw_complex));
      GLOBAL_LRSWAP_K[nl] = (fftw_complex *) tcalloc(NPIEROWS*NPIECOLS*PIE_ROW_DIA*PIE_COL_DIA,sizeof(fftw_complex));}
    GLOBAL_FFTWPLAN_KERNEL_FORWARD = (fftw_plan *) tcalloc(NLRKERNELS,sizeof(fftw_plan));
    GLOBAL_FFTWPLAN_KERNEL_BACKWARD = (fftw_plan *) tcalloc(NLRKERNELS,sizeof(fftw_plan));
    GLOBAL_FFTWPLAN_SWAP_FORWARD = (fftw_plan *) tcalloc(NLRKERNELS,sizeof(fftw_plan));
    GLOBAL_FFTWPLAN_SWAP_BACKWARD = (fftw_plan *) tcalloc(NLRKERNELS,sizeof(fftw_plan));
    for (nl=0;nl<NLRKERNELS;nl++){
      GLOBAL_FFTWPLAN_KERNEL_FORWARD[nl] = fftw_plan_dft_3d(PIE_ROW_DIA*PIE_COL_DIA,NPIECOLS,NPIEROWS,GLOBAL_LRKERNEL_X[nl],GLOBAL_LRKERNEL_K[nl],-1,FFTW_ESTIMATE);
      GLOBAL_FFTWPLAN_KERNEL_BACKWARD[nl] = fftw_plan_dft_3d(PIE_ROW_DIA*PIE_COL_DIA,NPIECOLS,NPIEROWS,GLOBAL_LRKERNEL_X[nl],GLOBAL_LRKERNEL_K[nl],+1,FFTW_ESTIMATE);
      fftw_execute_dft(GLOBAL_FFTWPLAN_KERNEL_FORWARD[nl],GLOBAL_LRKERNEL_X[nl],GLOBAL_LRKERNEL_K[nl]);
      GLOBAL_FFTWPLAN_SWAP_FORWARD[nl] = fftw_plan_dft_3d(PIE_ROW_DIA*PIE_COL_DIA,NPIECOLS,NPIEROWS,GLOBAL_LRSWAP_X[nl],GLOBAL_LRSWAP_K[nl],-1,FFTW_ESTIMATE);
      GLOBAL_FFTWPLAN_SWAP_BACKWARD[nl] = fftw_plan_dft_3d(PIE_ROW_DIA*PIE_COL_DIA,NPIECOLS,NPIEROWS,GLOBAL_LRSWAP_X[nl],GLOBAL_LRSWAP_K[nl],+1,FFTW_ESTIMATE);}}
  if (verbose){ printf(" %% making datastructures\n");}
  if (RTC_BOTHER){ GLOBAL_RTC = rtcmake(GLOBAL_RTC_LENGTH,GLOBAL_RTC_FRAMELENGTH,GLOBAL_RTC_NANGLES,GLOBAL_RTC_NPHASES);}
  if (STROBETRACE_BOTHER){ GLOBAL_STROBETRACE = strobetracemake(nget(GLOBAL_Nra,0,0),GLOBAL_STROBETRACE_NANGLES,GLOBAL_STROBETRACE_LENGTH,1.0,GLOBAL_STROBETRACE_CYCLE_BOTHER,GLOBAL_STROBETRACE_AVALANCHE_UPDATE_EVERY,GLOBAL_STROBETRACE_AVALANCHE_LOGDMIN,GLOBAL_STROBETRACE_AVALANCHE_LOGDMAX);}
  if (TUNINGCURVE_BOTHER){ GLOBAL_TUNINGCURVE = tuningcurvemake(GLOBAL_TUNINGCURVE_NANGLES,GLOBAL_TUNINGCURVE_NRADIUS);}
  if (LMITRI_BOTHER){ GLOBAL_LMITRI = lmitrimake(GLOBAL_LMITRI_T0,GLOBAL_LMITRI_TIMELENGTH,GLOBAL_LMITRI_ROW_MAX,GLOBAL_LMITRI_ROW_MIN);}
  if (PTREE_BOTHER){ GLOBAL_PTREE = ptreemake(GLOBAL_PTREE_NREGIONS,GLOBAL_PTREE_EVENT_WITHIN,GLOBAL_PTREE_EVENT_THRESHOLD,GLOBAL_PTREE_REGION_TYPE,GLOBAL_PTREE_NLEGS,GLOBAL_PTREE_LEGTIME);}
  if (CLOSET_BOTHER){ GLOBAL_CLOSET = closetmake(0,SYSTEM_ROW_SIZE,0,SYSTEM_COL_SIZE,GLOBAL_PTREE_LEGTIME,0.1);}
  if (YGGDRASIL_BOTHER){ GLOBAL_YGGDRASIL = yggdrasilmake(GLOBAL_PTREE_NREGIONS,GLOBAL_PTREE_EVENT_WITHIN,GLOBAL_PTREE_EVENT_THRESHOLD,GLOBAL_PTREE_REGION_TYPE,GLOBAL_PTREE_NLEGS,GLOBAL_PTREE_LEGTIME,GLOBAL_YGGDRASIL_PPNREGIONS,GLOBAL_YGGDRASIL_PPNLEGS,GLOBAL_YGGDRASIL_PPLEGTIME,GLOBAL_YGGDRASIL_WEIGHT_MINIMUM);}
  if (BONSAI_BOTHER){ GLOBAL_BONSAI = bonsaimake(GLOBAL_Nra,GLOBAL_PTREE_NREGIONS,GLOBAL_PTREE_EVENT_WITHIN,GLOBAL_PTREE_EVENT_THRESHOLD,GLOBAL_PTREE_REGION_TYPE,GLOBAL_PTREE_NLEGS,GLOBAL_PTREE_LEGTIME);}
  if (HYDRA_BOTHER){ GLOBAL_HYDRA = hydramake(GLOBAL_HYDRA_JUSTONTIME,GLOBAL_HYDRA_STAYONTIME,GLOBAL_HYDRA_DUMP_EVERY,&ptree_inputswitchon,&ptree_inputswitchoff);}
  if (LYAPUNOV_BOTHER){ GLOBAL_LYAPUNOV = lyapunovmake(GLOBAL_Nra,GLOBAL_LYAPUNOV_UPDATE_EVERY,GLOBAL_LYAPUNOV_JIGGLE);}
  if (POWER_BOTHER){ GLOBAL_POWER = powermake(GLOBAL_Nra,GLOBAL_POWER_LENGTH,GLOBAL_POWER_HOWMANY);}
  if (TAOF_BOTHER){ GLOBAL_TAOF = taofmake(GLOBAL_Nra,GLOBAL_TAOF_LENGTH,GLOBAL_TAOF_STEP_EVERY);}
  if (SEIDCORR_BOTHER){ GLOBAL_SEIDCORR = seidcorrmake(GLOBAL_Nra,GLOBAL_SEIDCORR_SPACE_BIN_SIZE,GLOBAL_SEIDCORR_TIME_BIN_SIZE,GLOBAL_SEIDCORR_LENGTH,GLOBAL_SEIDCORR_TIME_START);}
  if (GLOBAL_RECORD_NUMBER>0){ 
    sprintf(text,"./%srecord_VS",GLOBAL_STRING_2); 
    if ((fp=fopen(text,"r"))!=NULL){ fclose(fp); sprintf(command,"rm %s;",text); system(command);}
    sprintf(text,"./%srecord_sN",GLOBAL_STRING_2); 
    if ((fp=fopen(text,"r"))!=NULL){ fclose(fp); sprintf(command,"rm %s;",text); system(command);}}
}

void setdependencies()
{
  /* merely enforces global dependencies */
  LR_BOTHER = LR_BOTHER && CORTEX_BOTHER;
  PNM_BOTHER =  PNM_BOTHER && (GLOBAL_RECORD_NUMBER>0);
  RTC_BOTHER = RTC_BOTHER && CORTEX_BOTHER && LGN_BOTHER;
  TUNINGCURVE_BOTHER = TUNINGCURVE_BOTHER && CORTEX_BOTHER && LGN_BOTHER;
  LMITRI_BOTHER = LMITRI_BOTHER && CORTEX_BOTHER;
  HANSHELLEY_FLAG = HANSHELLEY_FLAG && LGN_DETERMINISM && EIFvsIF<0 && BIG_SYSTEM_FLAG>1;
  STROBETRACE_BOTHER = STROBETRACE_BOTHER && CORTEX_BOTHER;
  GLOBAL_STROBETRACE_CYCLE_BOTHER = CORTEX_BOTHER && LGN_BOTHER && GRATING_VS_LMI==-1;
  GLOBAL_STROBETRACE_LENGTH = GLOBAL_STROBETRACE_CYCLE_BOTHER ? CYCLE_LENGTH : GLOBAL_STROBETRACE_LENGTH;
  PTREE_BOTHER = PTREE_BOTHER && CORTEX_BOTHER;
  CLOSET_BOTHER = CLOSET_BOTHER && CORTEX_BOTHER;
  YGGDRASIL_BOTHER = YGGDRASIL_BOTHER && CORTEX_BOTHER;
  BONSAI_BOTHER = BONSAI_BOTHER && CORTEX_BOTHER;
  HYDRA_BOTHER = HYDRA_BOTHER && CORTEX_BOTHER;
  LYAPUNOV_BOTHER = LYAPUNOV_BOTHER && CORTEX_BOTHER;
  POWER_BOTHER = POWER_BOTHER && CORTEX_BOTHER;
  TAOF_BOTHER = TAOF_BOTHER && CORTEX_BOTHER;
  SEIDCORR_BOTHER = SEIDCORR_BOTHER && CORTEX_BOTHER;
  //if (GLOBAL_RECORD_NUMBER>0){ GLOBAL_RECORD_NUMBER = findrecordnumber(); srand(GLOBAL_RECORD_NUMBER);} 
  sprintf(GLOBAL_STRING_2,"%s%d",GLOBAL_STRING,GLOBAL_RECORD_NUMBER);
}

void computestep()
{
  int ESr=0,ECr=0,ISr=0,ICr=0;
  int verbose=GLOBAL_verbose,dtok=0,pspikes=0,lllengthmax=0;
  double t = GLOBAL_time,DT=GLOBAL_DT,DTmax=GLOBAL_DTmax;
  struct neuronarray *Nra=NULL;
  struct llist **Sra=NULL;
  struct llist **Gra=NULL;
  struct llist *LL=NULL,*ll=NULL;
  struct litem *l=NULL;
  if (!RUN_DONE){
    if (CORTEX_BOTHER){
      dtok = GLOBAL_dtok; Nra = GLOBAL_Nra;
      if (verbose){ memprintf(0);}
      if (dtok > 4 && GLOBAL_dtadapt){ if (verbose){ printf("%% DT=%e has been successful for %d steps, doubling to make DT=%e\n",DT,dtok,2*DT);} DT=2*DT;}
      if (DT > DTmax){ if (verbose){ printf(" %% on second thought, DT large enough... setting DT=%e\n",DTmax);} DT=DTmax;}
      if (verbose){ printf("%% setting input spiketimes\n");}
      spikeinput(Nra,t,DT);
      if (verbose){ printf("%% scanning for potential spikes\n");}
      do{
	if (verbose){ printf("%% \t unmaking previous llists\n");}
	llistunmake(Sra,Gra,LL);
	if (verbose){ printf("%% \t setting blocksize\n");}
	setblocksize(DT);
	if (verbose){ printf("%% \t remaking llists\n");}
	llistremake(&Sra,&Gra,&LL);
	if (verbose){ printf("%% \t scanning...\n");}
	pspikes = spikescan(Nra,Sra,Gra,t,DT);
	if (verbose){ printf("%% \t found %d potential spikes\n",pspikes);}
	if (pspikes > SPIKETOL*NARBORS_TALL*NARBORS_WIDE/2.0 && GLOBAL_dtadapt){
	  if (verbose){ printf("%% \t categorically too many\n");}
	  lllengthmax = pspikes;}
	else{ 
	  if (verbose){ printf("%% \t sorting into local clusters...\n");}
	  spikesort(Sra,LL);
	  if (verbose){ printf("%% \t ...finished sorting into local clusters\n");}
	  lllengthmax=0;
	  l = LL->first;
	  while (l!=NULL){
	    ll = (struct llist *) l->item;
	    lllengthmax = maximum(lllengthmax,ll->length);
	    l = l->child;}}
	if (verbose){ printf("%% \n");}
	if (lllengthmax > SPIKETOL && GLOBAL_dtadapt){
	  if (verbose){ printf(" %%  pspikes=%d, lllengthmax=%d > %d, halving DT=%e to be DT/2=%e\n",pspikes,lllengthmax,SPIKETOL,DT,DT/2);}
	  DT = DT/2;
	  dtok=0;}}
      while (lllengthmax > SPIKETOL && GLOBAL_dtadapt);
      if (verbose){ printf(" %% step accepted with browdia=%d,bcoldia=%d,t=%0.3f,DT=%0.3f,pspikes=%d,lllengthmax=%d\n",GLOBAL_BLOCK_ROW_DIA,GLOBAL_BLOCK_COL_DIA,t,DT,pspikes,lllengthmax);}
      if (verbose){ printf(" %% correcting spiketimes locally and evolving spiking neurons \n");}
      spikecorrect(LL,t,DT); spikecorrect(LL,t,DT); 
      spikeconduct(LL,t,DT,&ESr,&ISr,&ECr,&ICr);
      Nra->mES = ESr*1000/DT/(double)Nra->t2stotal[3];
      Nra->mEC = ECr*1000/DT/(double)Nra->t2stotal[2];
      Nra->mIS = ISr*1000/DT/(double)Nra->t2stotal[1];
      Nra->mIC = ICr*1000/DT/(double)Nra->t2stotal[0];
      if (verbose){ printf(" %% %d ES, %d IS, %d EC, %d IC spiked... evolving nonspiking neurons \n",ESr,ISr,ECr,ICr);}
      gizmoconduct(Sra,Gra,t,DT);
      if (verbose){ printf(" %% NOT correcting the far-near-field conductances \n");}
      if (LR_BOTHER){
	if (verbose){ printf(" %% long range connections\n");}
	lrupdate(Nra,t,DT);}
      if (verbose){ printf(" %% freeing Sra,Gra,LL\n");}
      llistunmake(Sra,Gra,LL);
      if (verbose){ printf(" %% updating data structures...");}
      if (RTC_BOTHER){ rtcupdate(GLOBAL_RTC,t,DT);}
      if (STROBETRACE_BOTHER){ strobetraceupdate(GLOBAL_STROBETRACE,t,DT);}
      if (TUNINGCURVE_BOTHER){ tuningcurveupdate(GLOBAL_TUNINGCURVE,t,DT);}
      if (LMITRI_BOTHER){ lmitriupdate(GLOBAL_LMITRI,t,DT);}
      if (PTREE_BOTHER){ ptreeupdate(GLOBAL_PTREE,t,DT,1);}
      if (CLOSET_BOTHER){ ptreeupdate(GLOBAL_CLOSET->p,t,DT,1);}
      if (YGGDRASIL_BOTHER){ yggdrasilupdate(GLOBAL_YGGDRASIL,t,DT,1);}
      if (BONSAI_BOTHER){ bonsaiupdate(GLOBAL_BONSAI,t,DT,1);}
      if (HYDRA_BOTHER){ hydraupdate(GLOBAL_HYDRA,t,DT,1);}
      if (LYAPUNOV_BOTHER){ lyapunovupdate(GLOBAL_LYAPUNOV,t,DT);}
      if (POWER_BOTHER){ powerupdate(GLOBAL_POWER,t,DT);}
      if (TAOF_BOTHER){ taofupdate(GLOBAL_TAOF,t,DT);}
      if (SEIDCORR_BOTHER){ seidcorrupdate(GLOBAL_SEIDCORR,t,DT);}
      if (verbose){ printf(" finished\n");}}
    if (SUITE_BOTHER){ lgnswitch(GLOBAL_LGN,t,DT);} 
    if (LGN_BOTHER){ 
      lgnevolve(GLOBAL_LGN,t,DT); if (verbose){ printf(" %% updated lgn\n");}}
    if (verbose){ printf("\n %%%% finished for time interval [%f,%f] %%%% \n\n",t,t+DT);}
    if (verbose){ memprintf(0);}
    t = t+DT;
    dtok++;
    GLOBAL_dtok=dtok;
    GLOBAL_time=t;
    GLOBAL_DT=DT;
    if (CORTEX_BOTHER){ setinputrate(GLOBAL_Nra,GLOBAL_time);}
    if (!SUITE_BOTHER && ((OUTPUT_DUMP_EVERY>0 && (int)floor(GLOBAL_time/OUTPUT_DUMP_EVERY)<(int)floor((GLOBAL_time+GLOBAL_DT)/OUTPUT_DUMP_EVERY)) || (GLOBAL_TF > 0 && GLOBAL_time > GLOBAL_TF))){
      printf("dumping output at time %f...",GLOBAL_time);
      if (CORTEX_BOTHER && BIG_SYSTEM_FLAG<0){ connectionsdump(GLOBAL_Nra,1,"connections");}
      if (STROBETRACE_BOTHER){ strobetracedump(GLOBAL_STROBETRACE,0); strobetracedump(GLOBAL_STROBETRACE,1);}
      if (RTC_BOTHER){ rtcdump(GLOBAL_RTC,0); rtcdump(GLOBAL_RTC,1);}
      if (TUNINGCURVE_BOTHER){ tuningcurvedump(GLOBAL_TUNINGCURVE,1);}
      if (LMITRI_BOTHER){ lmitridump(GLOBAL_LMITRI,1);}
      if (PTREE_BOTHER){
	ptreerate(GLOBAL_PTREE);
	ptreedump_starter(GLOBAL_PTREE,NULL,2,0,0,0,+1,-1);
	if (GLOBAL_PTREE_BITBYBIT){
	  ptreereset(GLOBAL_PTREE);}}
      if (YGGDRASIL_BOTHER){
	ptreerate(GLOBAL_YGGDRASIL->p);
	ptreedump_starter(GLOBAL_YGGDRASIL->p,"yggdrasil_p",2,-1,0,0,+1,-1);
	if (GLOBAL_YGGDRASIL->pp!=NULL){ 
	  ptreerate(GLOBAL_YGGDRASIL->pp);
	  ptreedump_starter(GLOBAL_YGGDRASIL->pp,"yggdrasil_pp",2,-1,0,0,+1,-1);}}
      if (BONSAI_BOTHER){ bonsaidump(GLOBAL_BONSAI,"bonsai",0);bonsaidump(GLOBAL_BONSAI,"bonsai",1);bonsaidump(GLOBAL_BONSAI,"bonsai",2);}
      if (HYDRA_BOTHER){ hydradump(GLOBAL_HYDRA);}
      if (LYAPUNOV_BOTHER){ lyapunovdump(GLOBAL_LYAPUNOV,0); lyapunovdump(GLOBAL_LYAPUNOV,1);}
      if (POWER_BOTHER){ powerdump(GLOBAL_POWER,NULL,0); powerdump(GLOBAL_POWER,NULL,1);}
      if (TAOF_BOTHER){ taofdump(GLOBAL_TAOF,NULL,3);}
      printf("output dumped at time %f\n",GLOBAL_time);}
    if (GLOBAL_TF > 0 && GLOBAL_time > GLOBAL_TF){ if (GLOBAL_CLEANUP){ cleanupoutput();} exit(EXIT_SUCCESS);}}
}

/* Here is the main function */

int main(int argc, char **argv) 
{
  meminit(); /* initialize memory list for use */
  ping();
  readinput(); /* get all the global variables from some input */
  setdependencies(); /* make sure a few things are kosher */
  if (argc>1){ /* we have options */ 
    switch(getopt(argc,argv,"R:_:")){
    case 'R':
      GLOBAL_RECORD_NUMBER=atoi(optarg); printf(" %% GLOBAL_RECORD_NUMBER set to %d\n",GLOBAL_RECORD_NUMBER);
      break;
    case '_': 
      switch (atoi(optarg)){ 
      case 0: trialaverage(argc,argv); break; 
      case 1: ptree_trialaverage(argc,argv); break; 
      case 2: ptree_trialdistribution(argc,argv); break;
      case 3: seidcorr_compile(argc,argv); break;
      default: break;} 
      exit(EXIT_SUCCESS); break;
    default: printf(" need option type: either '-RGLOBAL_RECORD_NUMBER', or \n\t '-_0' (trialaverage) \n\t '-_1' (ptree_trialaverage) \n\t '-_2' (ptree_trialdistribution) \n\t '-_3' (seidcorr_compile)\n"); exit(EXIT_SUCCESS); break;}}
  setglobals(); /* initialize global variables */
/*   GLOBAL_PTREE = ptreadback("/home/rangan/Cstuff/v1/darena/dir_d_input_17/ptree_stayoff__disc17_1_00record_s0_all"); */
#ifndef NOTGRAPH
  if (SUPERGLOBAL_DRAW_FLAG){
    glutInit(&argc, argv);  /* initialize glut */
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_ALPHA | GLUT_DEPTH); /* select display mode, 2xbuffer,rgba,alpha,depth_buffer */
    glutInitWindowSize(GLOBAL_WINDOW_WIDTH, GLOBAL_WINDOW_HEIGHT); /* small window */
    glutInitWindowPosition(0,0); /* starts at upper left corner of screen */
    GLUTWINDOWNUMBER = glutCreateWindow("d_v1"); /* open window */
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
