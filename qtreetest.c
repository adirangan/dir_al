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

int ipow(int a,int b){ return (int)floor(pow(a,b));}

void simplexsample(double *ra,int length)
{
  int nr=0;
  double sum=0;
  for (nr=0;nr<length;nr++){ ra[nr] = -log(rand01); sum += ra[nr];}
  for (nr=0;nr<length;nr++){ ra[nr] /= sum;}
}

void rawsample(double *ra,int length)
{
  int nr=0;
  simplexsample(ra,length);
  for (nr=0;nr<length;nr++){ ra[nr] = ceil(ra[nr]*500);}
}

void num2regionhist(int num,int nregions,int nlegs,int *regionhist)
{
  /* assumes that regionhist has nlegs elements and num is less than pow(2*nregions,nlegs) */
  int verbose=0;
  int nr=0;
  if (verbose){ printf("given num %d, produced: ",num);}
  while (num>0 && nr<nlegs){
    regionhist[nr] = num%(2*nregions);
    num -= regionhist[nr];
    num /= 2*nregions;
    nr += 1;}
  if (verbose){ for (nr=0;nr<nlegs;nr++){ printf("%d ",regionhist[nr]);} printf("\n");}
}

void regionhist2char(int *regionhist,int nregions,int nlegs,char *chartemp)
{
  /* assumes that *regionhist has nlegs elements */
  int nl=0;
  sprintf(chartemp,"[");
  for (nl=nlegs-1;nl>=1;nl--){ 
    sprintf(chartemp,"%s%s%d|",chartemp,regionhist[nl]>=nregions ? "-" : "+",regionhist[nl]%nregions);}
  sprintf(chartemp,"%s%s%d]",chartemp,regionhist[0]>=nregions ? "-" : "+",regionhist[0]%nregions);
}

void num2state(int num,int nregions,int *state)
{
  /* assumes that *state has nregions elements and num is less than pow(2,nregions) */
  int verbose=0;
  int nr=0,base=2;
  for (nr=0;nr<nregions;nr++){
    state[nr] = (num%base)/(base/2);
    base *= 2;}
  if (verbose){
    printf("input num %d ",num);
    printf("output state ("); 
    for (nr=nregions-1;nr>=0;nr--){ printf("%d",state[nr]);} 
    printf(")\n");}
}

void num2statehist(int num,int nregions,int nlegs,int *statehist)
{
  /* assumes that *statehist has nregions*nlegs elements and num is less than pow(pow(2,nregions),nlegs) */
  int verbose=0;
  int nl=0,statebase=ipow(2,nregions),state=0;
  for (nl=0;nl<nlegs;nl++){
    state = (num%statebase)/(statebase/ipow(2,nregions));
    if (verbose){ printf("num %d, nl %d, state %d\n",num,nl,state);}
    num2state(state,nregions,&(statehist[0 + nl*nregions]));
    statebase *= ipow(2,nregions);}
}

void statehist2char(int *statehist,int nregions,int nlegs,char *chartemp)
{
  /* assumes that *statehist has nregions*nlegs elements */
  int nl=0,nr=0;
  sprintf(chartemp,"%s","");
  for (nl=nlegs-1;nl>=0;nl--){ 
    sprintf(chartemp,"%s(",chartemp); 
    for (nr=nregions-1;nr>=0;nr--){ sprintf(chartemp,"%s%d",chartemp,statehist[nr+nl*nregions]);}
    sprintf(chartemp,"%s), ",chartemp); }
}

void binarize(double *stree,int nregions,int nlegs,int index)
{
  /* uses binary indicator to collapse *stree over (the assumedly final) index<nlegs */
  int verbose=0;
  char chartemp[32];
  int nstates=ipow(2,nregions);
  int streepreindex = ipow(nstates,index);
  int preindex=0,nr=0,ns=0;
  double **qyes=NULL,**qnot=NULL,*qyessum=NULL,*qnotsum=NULL;
  if (verbose){ printf("[entered binarize] with nregions %d, nlegs %d, index %d\n",nregions,nlegs,index);}
  if (verbose){ printf(" passed stree:\n"); raprintf(stree,"double",nstates,ipow(nstates,index),"passed stree");}
  for (preindex=0;preindex<streepreindex;preindex++){
    if (verbose){ 
      printf(" preindex %d streepreindex %d, considering:\n",preindex,streepreindex); 
      for (ns=0;ns<nstates;ns++){ printf("%f ",stree[preindex + ns*streepreindex]);} printf("\n");}
    qyes = (double **) tcalloc(nregions,sizeof(double *));
    qnot = (double **) tcalloc(nregions,sizeof(double *));
    qyes[0] = (double *) tcalloc(ipow(2,nregions-1),sizeof(double));
    qnot[0] = (double *) tcalloc(ipow(2,nregions-1),sizeof(double));
    for (ns=0;ns<ipow(2,nregions-1);ns++){
      qyes[0][ns] = stree[preindex + (2*ns+1)*streepreindex];
      qnot[0][ns] = stree[preindex + (2*ns+0)*streepreindex];}
    for (nr=1;nr<nregions;nr++){
      qyes[nr] = (double *) tcalloc(ipow(2,nregions-1-nr),sizeof(double));
      qnot[nr] = (double *) tcalloc(ipow(2,nregions-1-nr),sizeof(double));
      for (ns=0;ns<ipow(2,nregions-1-nr);ns++){
	qyes[nr][ns] = qyes[nr-1][2*ns+1] + qnot[nr-1][2*ns+1];
	qnot[nr][ns] = qyes[nr-1][2*ns] + qnot[nr-1][2*ns];}}
    qyessum = (double *) tcalloc(nregions,sizeof(double));
    qnotsum = (double *) tcalloc(nregions,sizeof(double));
    for (nr=0;nr<nregions;nr++){
      for (ns=0;ns<ipow(2,nregions-1-nr);ns++){
	qyessum[nr] += qyes[nr][ns];
	qnotsum[nr] += qnot[nr][ns];}}
    if (verbose){
      printf(" produced\n"); 
      for (nr=0;nr<nregions;nr++){
	sprintf(chartemp,"qyes[%d]",nr); raprintf(qyes[nr],"double",1,ipow(2,nregions-1-nr),chartemp);
	sprintf(chartemp,"qnot[%d]",nr); raprintf(qnot[nr],"double",1,ipow(2,nregions-1-nr),chartemp);}
      raprintf(qyessum,"double",1,nregions,"qyessum");
      raprintf(qnotsum,"double",1,nregions,"qnotsum");}
    /* overwriting should be fine, since 2*nregions<nstates */
    for (nr=0;nr<nregions;nr++){
      stree[preindex + nr*streepreindex] = qyessum[nr];
      stree[preindex + (nr+nregions)*streepreindex] = qnotsum[nr];}
    for (nr=0;nr<nregions;nr++){ tfree(qyes[nr]); qyes[nr]=NULL; tfree(qnot[nr]); qnot[nr]=NULL;}
    tfree(qyes); qyes=NULL; tfree(qnot); qnot=NULL;
    tfree(qyessum); qyessum=NULL; tfree(qnotsum); qnotsum=NULL;}
}

int base2base(int num,int base1,int base2)
{
  /* stretch [region_1 + region_2*base1 + .. + region_nlegs*ipow(base1,nlegs-1)]
     towards [region_1 + region_2*base2 + .. + region_nlegs*ipow(base2,nlegs-1)] */
  int verbose=0;
  int *regions=NULL;
  int nlegs = num>0 ? (int)ceil(log(num)/log(base1)) : 1;
  int nl=0,nl2=0;
  int num2=0;
  if (verbose){ printf(" [entering base2base] with num %d base1 %d base2 %d\n",num,base1,base2);}
  regions = (int *) tcalloc(nlegs+1,sizeof(int));
  while (num>0){
    regions[nl] = num%(base1);
    num -= regions[nl];
    num /= (base1);
    nl += 1;}
  if (verbose){ printf(" found regions: "); for (nl2=0;nl2<nl;nl2++){ printf("%d ",regions[nl2]);} printf("\n");}
  for (nl2=0;nl2<nl;nl2++){
    num2 += regions[nl2]*ipow(base2,nl2);}
  tfree(regions); regions=NULL;
  if (verbose){ printf(" returning num2 %d\n",num2);}
  return num2;
}

void basesquish(double *ra,int base1,int nlegs,int base2)
{
  /* squish ra[region_1 + region_2*base1 + .. + region_nlegs*ipow(base1,nlegs-1)]
     making ra[region_1 + region_2*base2 + .. + region_nlegs*ipow(base2,nlegs-1)] */
  int nh=0;
  for (nh=0;nh<ipow(base2,nlegs);nh++){ ra[nh] = ra[base2base(nh,base2,base1)];}
}

double * stree2qtree(double *stree,int nregions,int nlegs)
{
  /* assumes that raw *stree has nstatehist elements, and raw *qtree is nregionhist elements*/
  int nstates=ipow(2,nregions);
  int nstatehist=ipow(nstates,nlegs);
  int nregionhist=ipow(2*nregions,nlegs);
  int nh=0,nh2=0,nl=0;
  double *qtree=NULL;
  qtree = tcalloc(nstatehist,sizeof(double));
  for (nh=0;nh<nstatehist;nh++){ qtree[nh]=stree[nh];}
  for (nl=0;nl<nlegs;nl++){
    for (nh=0;nh<ipow(2*nregions,nl);nh++){
      nh2 = base2base(nh,2*nregions,nstates);
      binarize(&(qtree[0 + nh2*ipow(nstates,nlegs-nl)]),nregions,nlegs,nlegs-1-nl);}}
  basesquish(qtree,nstates,nlegs,2*nregions);
  qtree = trealloc(qtree,nregionhist*sizeof(double));
  return qtree;
}

double * stree2qtree_slow(double *stree,int nregions,int nlegs)
{
  /* slowly generates raw *qtree of nregionhist elements */
  int nstates=ipow(2,nregions);
  int nstatehist=ipow(nstates,nlegs);
  int nregionhist=ipow(2*nregions,nlegs);
  int nr=0,nr2=0,nh=0,nh2=0,nl=0;
  double *qtree=NULL;
  int *regions=NULL,*states=NULL;
  int use_flag=0;
  qtree = tcalloc(nregionhist,sizeof(double));
  for (nr=0;nr<nregionhist;nr++){
    nr2=nr;nl=0;
    regions = (int *) tcalloc(nlegs,sizeof(int));
    while (nr2>0){
      regions[nl] = nr2%(2*nregions);
      nr2 -= regions[nl];
      nr2 /= 2*nregions;
      nl += 1;}
    for (nh=0;nh<nstatehist;nh++){
      nh2=nh;nl=0;
      states = (int *) tcalloc(nlegs,sizeof(int));
      while (nh2>0){
	states[nl] = nh2%nstates;
	nh2 -= states[nl];
	nh2 /= nstates;
	nl += 1;}
      use_flag=1;
      for (nl=0;nl<nlegs;nl++){
	if (regions[nl]>=nregions){ /* nonevent */
	  use_flag *= ((states[nl]%ipow(2,regions[nl]-nregions+1))/ipow(2,regions[nl]-nregions) == 0);}
	else if (regions[nl]<nregions){ /* event */
	  use_flag *= ((states[nl]%ipow(2,regions[nl]+1))/ipow(2,regions[nl]) == 1);}}
      if (use_flag){ qtree[nr] += stree[nh];}
      tfree(states); states=NULL;}
    tfree(regions); regions=NULL;}
  return qtree;
}

void recombinarize(double *qtree,int nregions,int nlegs,int index)
{
  /* recombines *qtree (assumed to be nstatehist elements) over (assumedly final) index<nlegs */
  int verbose=0;
  char chartemp[32];
  int nstates=ipow(2,nregions);
  int qtreepreindex = ipow(nstates,index);
  int preindex=0,nr=0,ns=0;
  double **stemp=NULL;
  if (verbose){ printf("[entered recombinarize] with nregions %d, nlegs %d, index %d\n",nregions,nlegs,index);}
  if (verbose){ printf(" passed qtree:\n"); raprintf(qtree,"double",nstates,ipow(nstates,index),"passed qtree");}
  for (preindex=0;preindex<qtreepreindex;preindex++){
    if (verbose){ 
      printf(" preindex %d qtreepreindex %d, considering:\n",preindex,qtreepreindex); 
      for (nr=0;nr<2*nregions;nr++){ printf("%f ",qtree[preindex + nr*qtreepreindex]);} printf("\n");}
    stemp = (double **) tcalloc(nregions,sizeof(double *));
    stemp[nregions-1] = (double *) tcalloc(ipow(2,nregions-(nregions-1)),sizeof(double));
    for (ns=0;ns<ipow(2,nregions-(nregions-1));ns++){
      stemp[nregions-1][ns] = qtree[preindex + (nregions-1 + (ns%2==0)*nregions)*qtreepreindex];}
    for (nr=nregions-2;nr>=0;nr--){
      stemp[nr] = (double *) tcalloc(ipow(2,nregions-nr),sizeof(double));
      for (ns=0;ns<ipow(2,nregions-nr);ns++){
	stemp[nr][ns] = stemp[nr+1][ns/2]*qtree[preindex + (nr + (ns%2==0)*nregions)*qtreepreindex];}}
    if (verbose){
      printf(" produced\n"); 
      for (nr=nregions-1;nr>=0;nr--){ sprintf(chartemp,"stemp[%d]",nr); raprintf(stemp[nr],"double",1,ipow(2,nregions-nr),chartemp);}}
    /* overwriting should be fine, since 2*nregions<nstates */
    for (ns=0;ns<nstates;ns++){
      qtree[preindex + ns*qtreepreindex] = stemp[0][ns];}
    for (nr=0;nr<nregions;nr++){ tfree(stemp[nr]); stemp[nr]=NULL;} tfree(stemp); stemp=NULL;}
}

double * qtree2stree(double *qtree,int nregions,int nlegs)
{
  /* assumes that normalized *qtree has nregionhist elements, and normalized *stree is nstatehist elements */
  int nstates=ipow(2,nregions);
  int nstatehist=ipow(nstates,nlegs);
  int nregionhist=ipow(2*nregions,nlegs);
  int nh=0,nh2=0,nl=0;
  double *stree=NULL;
  stree = tcalloc(nstatehist,sizeof(double));
  for (nh=0;nh<nregionhist;nh++){ nh2 = base2base(nh,2*nregions,nstates); stree[nh2]=qtree[nh];}
  for (nl=0;nl<nlegs;nl++){
    for (nh=0;nh<ipow(nstates,nl);nh++){
      recombinarize(&(stree[0 + nh*ipow(nstates,nlegs-nl)]),nregions,nlegs,nlegs-1-nl);}}
  return stree;
}

void streenormalize(double *stree,int nregions,int nlegs)
{
  /* normalizes *stree of nstatehist elements */
  int verbose=0;
  int nstates=ipow(2,nregions);
  int nstatehist=ipow(nstates,nlegs);
  int ns=0,nh=0;
  double sum=0;
  int *statehist=NULL;
  char *chartemp=NULL;
  if (verbose){ statehist = (int *) tcalloc(nregions*nlegs,sizeof(int)); chartemp = (char *) tcalloc(32,sizeof(char));}
  for (nh=0;nh<nstatehist/nstates;nh++){
    sum=0;
    if (verbose){ printf("summing legs: ");}
    for (ns=0;ns<nstates;ns++){
      sum += stree[ns + nh*nstates];
    if (verbose){ 	
      num2statehist(ns+nh*nstates,nregions,nlegs,statehist);statehist2char(statehist,nregions,nlegs,chartemp);
      printf("%s- ",chartemp);}}
    if (verbose){ printf(" sum=%f\n",sum);}
    for (ns=0;ns<nstates;ns++){
      stree[ns + nh*nstates] /= sum;}}
  if (verbose){ tfree(statehist); statehist=NULL; tfree(chartemp); chartemp=NULL;}
}

void qtreenormalize(double *qtree,int nregions,int nlegs)
{
  /* normalizes *qtree of nregionhist elements */
  int verbose=0;
  int nregionhist=ipow(2*nregions,nlegs);
  int nr=0,nh=0;
  double sum=0;
  int *regionhist=NULL;
  char *chartemp=NULL;
  if (verbose){ regionhist = (int *) tcalloc(nregions*nlegs,sizeof(int)); chartemp = (char *) tcalloc(32,sizeof(char));}
  for (nh=0;nh<nregionhist/(2*nregions);nh++){
    for (nr=0;nr<nregions;nr++){
      sum = qtree[nr+nh*(2*nregions)]+qtree[nr+nregions+nh*(2*nregions)];
      if (verbose){ 
	printf("summing legs: ");
	num2regionhist(nr+nh*(2*nregions),nregions,nlegs,regionhist);regionhist2char(regionhist,nregions,nlegs,chartemp);
	printf("%s- ",chartemp);
	num2regionhist(nr+nregions+nh*(2*nregions),nregions,nlegs,regionhist);regionhist2char(regionhist,nregions,nlegs,chartemp);
	printf("%s",chartemp);
	printf(" sum=%f\n",sum);}
      if (sum>0){
	qtree[nr+nh*(2*nregions)] /= sum;
	qtree[nr+nregions+nh*(2*nregions)] /= sum;}
      else{
	qtree[nr+nh*(2*nregions)] = 0.5;
	qtree[nr+nregions+nh*(2*nregions)] = 0.5;}}}
  if (verbose){ tfree(regionhist); regionhist=NULL; tfree(chartemp); chartemp=NULL;}
}

double radist(double *ra1,int length,double *ra2)
{
  /* frobnorm */
  double norm=0;
  int nr=0;
  for (nr=0;nr<length;nr++){ norm += pow(ra1[nr]-ra2[nr],2);}
  return sqrt(norm);
}

double rareldist(double *ra1,int length,double *ra2)
{
  /* average relative difference, assuming ra1 is correct */
  double norm=0;
  int nr=0;
  for (nr=0;nr<length;nr++){ norm += fabs((ra1[nr]-ra2[nr])/ra1[nr]);}
  return norm/length;
}

int main(int argc, char **argv) 
{
  int verbose=1;
  int iteration=0,iteration_max=0;
  int nregions=3,nlegs=2;
  int nstates=ipow(2,nregions);
  int nstatehist=ipow(nstates,nlegs);
  int nregionhist=ipow(2*nregions,nlegs);
  int nh=0;
  int *statehist=NULL;
  double *stree=NULL;
  double *qtree=NULL;
  double *sqtree=NULL;
  meminit(); /* initialize memory list for use */
  if (argc>1){ /* we have options */ 
    switch(getopt(argc,argv,"M:")){
    case 'M': 
      iteration_max=atoi(optarg); break;
    default: printf(" need option type:  [-Miteration_max]\n"); exit(EXIT_SUCCESS); break;}}
  do{
    statehist = (int *) tcalloc(nregions*nlegs,sizeof(int));
    /* stored as [state_1 + state_2*nstates + .. + state_nlegs*pow(nstates,nlegs-1)] */
    stree = (double *) tcalloc(nstatehist,sizeof(double));
    for (nh=0;nh<nstatehist/nstates;nh++){ rawsample(&(stree[0 + nh*nstates]),nstates);}
    if (verbose>1){ raprintf(stree,"double",nstates,nstatehist/nstates,"stree");}
    qtree = stree2qtree(stree,nregions,nlegs);
    if (verbose>1){ raprintf(qtree,"double",2*nregions,nregionhist/(2*nregions),"qtree");}
    streenormalize(stree,nregions,nlegs);
    if (verbose>1){ raprintf(stree,"double",nstates,nstatehist/nstates,"stree");}
    qtreenormalize(qtree,nregions,nlegs);
    if (verbose>1){ raprintf(qtree,"double",2*nregions,nregionhist/(2*nregions),"qtree");}
    sqtree = qtree2stree(qtree,nregions,nlegs); streenormalize(sqtree,nregions,nlegs);
    if (verbose>1){ raprintf(sqtree,"double",nstates,nstatehist/nstates,"sqtree");}
    if (verbose){ 
      printf("difference=%f, average relative difference=%f\n",radist(stree,nstatehist,sqtree),rareldist(stree,nstatehist,sqtree));}
    tfree(sqtree); sqtree=NULL;
    tfree(qtree); qtree=NULL;
    tfree(statehist); statehist=NULL;
    tfree(stree); stree=NULL;
    iteration+=1;}
  while (iteration<iteration_max);
  return 1;
}
