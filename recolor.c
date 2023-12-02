#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
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

#define RGB3 3 // 3 bytes of color info per pixel
#define RGBA 4 // 4 bytes of color+alpha info
#define UCHAR_MAX 255 // is this right?

#include "d_llists.h"

/* These are the real global variables */
int MEMDEBUG=0;
void * MEMDEBUGLIST=NULL;

/* These are convenient global variables */
int GRAYSCALE=0;
double STD_VIEW=1;
int ON_MY_COMPUTER=1;

#include "d_llists.c"

int main(int argc, char **argv) 
{
  int verbose=1;
  int still_have_options = argc-1;
  char helpfile[1024],filename_base[512],filename[1024];
  int smoother=1;
  int colors=0,rows=0,cols=0;
  double max=0,min=0,mean=0,std=0;
  double *A=NULL,*B=NULL;
  double rcolor=0,gcolor=0,bcolor=0;
  meminit(); /* initialize memory list for use */
  sprintf(helpfile," need option type: -Ffilename_base [-Gsmoother] [-Sstdview]\n");
  if (still_have_options==0){ printf(helpfile); exit(EXIT_SUCCESS);}
  while (still_have_options){
    switch(getopt(argc,argv,"F:f:G:g:S:s:")){
    case 'F': case 'f': 
      sprintf(filename_base,"%s",optarg);
      break;
    case 'S': case 's': 
      STD_VIEW=atof(optarg); 
      break;
    case 'G': case'g':
      smoother=atoi(optarg);
      break;
    default: printf(helpfile); exit(EXIT_SUCCESS); break;}
    still_have_options--;}
  sprintf(filename,"%s.pnm",filename_base);
  A = ReadPNMfile(filename,1,1,&colors,&max,&min,&cols,&rows);
  B = spacesmear(A,rows,cols,smoother);
  stats("double",B,rows*cols,&max,&min,&mean,&std);
  sprintf(filename,"%s_adicolored_%d_%0.1f.pnm",filename_base,smoother,STD_VIEW);
  WritePNMfile_color(B,rows,cols,mean+STD_VIEW*std,mean-STD_VIEW*std,filename,7);
  tfree(A); tfree(B);
  return 1;
}
