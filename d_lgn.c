/* /\* This holds a lgn *\/ */
/* struct lgn */
/* { */
/*   char name[64]; /\* name of file *\/ */
/*   int rows; /\* rows in picture *\/ */
/*   int cols; /\* columns in picture *\/ */
/*   double height; /\* height in degrees *\/ */
/*   double width; /\* width in degrees *\/ */
/*   double irad; /\* inner radius (degrees) for lgn center, outer radius is 3 times inner radius *\/ */
/*   double lgnbaserate; /\* lgn base background firing rate *\/ */
/*   double peakrate; /\* peak firing rate (increase) of impulse response *\/ */
/*   fftw_plan plan_forward; /\* forward plan for fast fft *\/ */
/*   fftw_plan plan_backward; /\* backward plan for fast fft *\/ */
/*   double *pnm; /\* holds picture data *\/ */
/*   double *onshape; /\* spatial kernel for on-center cells *\/ */
/*   double *onin; /\* summed intensity input to lgn on-center cells *\/ */
/*   double *onrate; /\* used to determine lgn on-center firing rate *\/ */
/*   double *offrate; /\* used to determine lgn off-center firing rate *\/ */
/*   double **angleshape; /\* spatial kernel for gabor of a certain angle and phase *\/ */
/*   double **anglein; /\* linear input to gabors of various angles *\/ */
/*   double **anglerate; /\* used to determine various angular firing rates *\/ */
/* }; */

/* /\* Here are the lgn functions *\/ */
/* double * ReadPNMfile(char *,int,int); */
/* struct lgn * lgnmake(char *,int,int,double,double,double,double,double); */
/* void makegrating(struct lgn *,int); */
/* void lgnworkhorse(struct lgn *); */
/* double lgnonshapekernel(double,double,double,double,double); */
/* int lgndump(struct lgn *); */
/* int lgnread(struct lgn *); */
/* void lgnconvolve(struct lgn *); */
/* void lgnremake(struct lgn *,char *); */
/* void lgnevolve(struct lgn *,double,double); */
/* void lgntfree(struct lgn *); */

/* Here are the lgn functions */

struct lgn * lgnmake(char *filename,int rows,int cols,double height,double width,double irad,double lgnbaserate,double peakrate)
{
  /* This reads in the pnm file *filename into p->pnm, and processes the image to yield the rates p->onin */
  int na=0,np=0;
  struct lgn *p=NULL;
  p = (struct lgn *) tmalloc(sizeof(struct lgn));
  p->rows = rows; p->cols = cols; p->width = width; p->height = height; p->irad = irad; p->lgnbaserate = lgnbaserate;
  p->peakrate = peakrate;
  p->plan_forward = NULL; p->plan_backward = NULL;
  p->pnm = (double *) tcalloc(p->rows*p->cols,sizeof(double));
  p->onin = (double *) tcalloc(p->rows*p->cols,sizeof(double));
  p->onshape = (double *) tcalloc(p->rows*p->cols+4,sizeof(double));
  p->onrate = (double *) tcalloc(5*p->rows*p->cols,sizeof(double));
  p->offrate = (double *) tcalloc(5*p->rows*p->cols,sizeof(double));
  p->angleshape = (double **) tcalloc(NSLICES*NPHASES,sizeof(double*));
  for (na=0;na<NSLICES;na++){ for (np=0;np<NPHASES;np++){
    p->angleshape[na + np*NSLICES] = (double *) tcalloc(p->rows*p->cols+4,sizeof(double));}}
  p->anglein = (double **) tcalloc(NSLICES*NPHASES,sizeof(double*));
  for (na=0;na<NSLICES;na++){ for (np=0;np<NPHASES;np++){
    p->anglein[na + np*NSLICES] = (double *) tcalloc(p->rows*p->cols+4,sizeof(double));}}
  p->anglerate = (double **) tcalloc(NSLICES*NPHASES,sizeof(double*));
  for (na=0;na<NSLICES;na++){ for (np=0;np<NPHASES;np++){
    p->anglerate[na + np*NSLICES] = (double *) tcalloc(5*p->rows*p->cols,sizeof(double));}}
  sprintf(p->name,filename); lgnremake(p,filename);
  return p;
}

void lgnremake(struct lgn *p,char *filename)
{
  /* this plugs a new picture into the lgn structure */
  int verbose=0;
  int frame=framefig(GLOBAL_time);
  char textbase[128],textfull[256];
  fftw_complex *tmp1=NULL,*tmp2=NULL;
  if (filename!=NULL){ sprintf(textbase,filename); sprintf(p->name,textbase);}
  else /* if (filename==NULL) */{ 
    if (GRATING_VS_LMI==-2){
      sprintf(textbase,"grating_%dx%d_angle%0.2f_phase%0.2f_k%0.2f",p->rows,p->cols,INPUT_SPACEANGLE,INPUT_SPACEPHASE,INPUT_SPACEK);
      sprintf(p->name,textbase);}
    else if (GRATING_PULSE && !GRATING_DRIFT){ 
      sprintf(textbase,"grating_pulse_%dx%d_frame_%d",p->rows,p->cols,frame); 
      sprintf(p->name,textbase);}
    else if (!GRATING_PULSE && GRATING_DRIFT){
      sprintf(textbase,"grating_%dx%d_angle%0.2f_phase%0.2f_k%0.2f",p->rows,p->cols,INPUT_SPACEANGLE,INPUT_SPACEPHASE,INPUT_SPACEK);
      sprintf(p->name,textbase);}
    else{
      sprintf(textbase,"grating_wave_%dx%d_frame_%d",p->rows,p->cols,frame); 
      sprintf(p->name,textbase);}}
  if (verbose){ printf("Trying to read %s.lgn\n",textbase);}
  if (lgnread(p)){ if (verbose){ printf("Read %s.lgn\n",textbase);}}
  else{
    if (verbose){ printf("Couldn't read %s.lgn, building...\n",textbase);}
    tmp1 = (fftw_complex *) fftw_malloc(p->rows*p->cols*sizeof(fftw_complex));
    tmp2 = (fftw_complex *) fftw_malloc(p->rows*p->cols*sizeof(fftw_complex));
    if (p->plan_forward != NULL){ fftw_destroy_plan(p->plan_forward); p->plan_forward=NULL;}
    if (p->plan_backward != NULL){ fftw_destroy_plan(p->plan_backward); p->plan_backward=NULL;}
    p->plan_forward = fftw_plan_dft_2d(p->rows,p->cols,tmp1,tmp2,-1,FFTW_MEASURE);
    p->plan_backward = fftw_plan_dft_2d(p->rows,p->cols,tmp1,tmp2,+1,FFTW_MEASURE);
    sprintf(textfull,"./%s.pnm",textbase);
    if (filename!=NULL){ tfree(p->pnm); p->pnm = ReadPNMfile(textfull,1,0,NULL,NULL,NULL,NULL,NULL);}
    else /* if (filename==NULL) */{ makegrating(p,frame);}
    lgnworkhorse(p);
    free(tmp1);free(tmp2);
    if (LGN_DUMP){ if (verbose){ printf("Writing %s.lgn\n",textbase);} lgndump(p);}}
}

void makegrating(struct lgn *p,int frame)
{
  /* this makes frame (out of GLOBAL_FPS) drifting grating stimulus for p->pnm */
  int nr=0,nc=0;
  double nr2=0,nc2=0;
  double rowsperdegree=0,colsperdegree=0;
  double weight=0;
  double angle = INPUT_SPACEANGLE,k1 = 2*PI*INPUT_SPACEK*cos(angle),k2 = 2*PI*INPUT_SPACEK*sin(angle);
  double pulsew = GRATING_PULSE;
  rowsperdegree = (double)p->rows/(double)p->height; colsperdegree = (double)p->cols/(double)p->width;
  for (nr=0;nr<p->rows;nr++){ for (nc=0;nc<p->cols;nc++){
    nr2 = (double)(periodize(nr,-p->rows/2,p->rows/2))/(double)rowsperdegree; 
    nc2 = (double)(periodize(nc,-p->cols/2,p->cols/2))/(double)colsperdegree;
    weight = cos((double)frame/(double)GLOBAL_FPS*2*PI*pulsew)*cos(k2*nr2 + k1*nc2 - (double)INPUT_SPACEPHASE);
    p->pnm[nr + nc*p->rows] = weight;}}
}

void lgnworkhorse(struct lgn *p)
{
  /* This computes all the lgn convolutions and stores them */
  int verbose=0;
  int nr=0,nc=0,na=0,np=0,plength=p->rows*p->cols;
  double nr2=0,nc2=0;
  double rowsperdegree=0,colsperdegree=0;
  double weight=0,da=0;
  double angle=0,phase=0,k=0,std=0,anglekernel_normalizer=0,k1=0,k2=0,cos2=0;
  rowsperdegree = (double)p->rows/(double)p->height; colsperdegree = (double)p->cols/(double)p->width;
  da = (double)1/(double)rowsperdegree/(double)colsperdegree;
  for (nr=0;nr<p->rows;nr++){ for (nc=0;nc<p->cols;nc++){
    nr2 = (double)(periodize(nr,-p->rows/2,p->rows/2))/(double)rowsperdegree; 
    nc2 = (double)(periodize(nc,-p->cols/2,p->cols/2))/(double)colsperdegree;
    //weight = exp(-(pow(nr2,2) + pow(nc2,2))/(p->irad/2))/(PI*p->irad/2) - exp(-(pow(nr2,2) + pow(nc2,2))/(3*p->irad/2))/(PI*3*p->irad/2);
    weight = lgnonshapekernel(p->irad,nc2,nr2,0,0);
    p->onshape[nr + nc*p->rows] = weight;}}
  stats("double",p->onshape,plength,p->onshape+plength+0,p->onshape+plength+1,p->onshape+plength+2,p->onshape+plength+3);
  for (nr=0;nr<p->rows;nr++){ for (nc=0;nc<p->cols;nc++){ p->onshape[nr+nc*p->rows] -= p->onshape[plength+2];}}
  stats("double",p->onshape,plength,p->onshape+plength+0,p->onshape+plength+1,p->onshape+plength+2,p->onshape+plength+3);
  for (np=0;np<NPHASES;np++){
    phase = 0 + (double)(np+0.5)/(double)NPHASES*PI;//-PI/2 + (double)(np+0.5)/(double)NPHASES*PI;
    k = 2*PI/4.0/p->irad;
    std = sqrt(3*p->irad/2);
    anglekernel_normalizer = PI*std*std /* *cos(phase)*exp(-k*k/std/std) */;
    for (na=0;na<NSLICES;na++){ 
      angle = -PI/2 + (double)(na+0.5)/(double)NSLICES*PI;
      k1 = k*cos(angle); k2 = k*sin(angle);
      for (nr=0;nr<p->rows;nr++){ for (nc=0;nc<p->cols;nc++){
	nr2 = (double)(periodize(nr,-p->rows/2,p->rows/2))/(double)rowsperdegree; 
	nc2 = (double)(periodize(nc,-p->cols/2,p->cols/2))/(double)colsperdegree;
	cos2 = cos(k2*nr2 + k1*nc2 + phase); //cos2 = cos2 > 0 ? +1 : cos2 < 0 ? -1 : 0;
	weight = cos2*exp(-(pow(nr2,2) + pow(nc2,2))/std/std)/anglekernel_normalizer;
	p->angleshape[na + np*NSLICES][nr + nc*p->rows] = weight;}}}} 
  for (np=0;np<NPHASES;np++){ for (na=0;na<NSLICES;na++){ 
    stats("double",p->angleshape[na+np*NSLICES],plength,p->angleshape[na+np*NSLICES]+plength+0,p->angleshape[na+np*NSLICES]+plength+1,p->angleshape[na+np*NSLICES]+plength+2,p->angleshape[na+np*NSLICES]+plength+3);}}
  for (np=0;np<NPHASES;np++){ for (na=0;na<NSLICES;na++){ for (nr=0;nr<p->rows;nr++){ for (nc=0;nc<p->cols;nc++){ 
    p->angleshape[na+np*NSLICES][nr+nc*p->rows] -= p->angleshape[na+np*NSLICES][plength+2];}}}}
  for (np=0;np<NPHASES;np++){ for (na=0;na<NSLICES;na++){ 
    stats("double",p->angleshape[na+np*NSLICES],plength,p->angleshape[na+np*NSLICES]+plength+0,p->angleshape[na+np*NSLICES]+plength+1,p->angleshape[na+np*NSLICES]+plength+2,p->angleshape[na+np*NSLICES]+plength+3);}}
  fftwconvolve(&p->plan_forward,&p->plan_backward,p->pnm,p->onshape,p->onin,p->rows,p->cols);
  for (na=0;na<NSLICES;na++){ for (np=0;np<NPHASES;np++){
    fftwconvolve(&p->plan_forward,&p->plan_backward,p->pnm,p->angleshape[na+np*NSLICES],p->anglein[na+np*NSLICES],p->rows,p->cols);}}
  for (nr=0;nr<p->rows;nr++){ for (nc=0;nc<p->cols;nc++){
    p->onin[nr+nc*p->rows] *= da;
    for (na=0;na<NSLICES;na++){ for (np=0;np<NPHASES;np++){
      p->anglein[na+np*NSLICES][nr+nc*p->rows] *= da;}}}}
  for (np=0;np<NPHASES;np++){ for (na=0;na<NSLICES;na++){ 
    stats("double",p->anglein[na+np*NSLICES],plength,p->anglein[na+np*NSLICES]+plength+0,p->anglein[na+np*NSLICES]+plength+1,p->anglein[na+np*NSLICES]+plength+2,p->anglein[na+np*NSLICES]+plength+3);}}
  if (verbose){
    printf("p->onshape,max=%0.2f,min=%0.2f,mean=%0.2f,std=%0.2f\n",p->onshape[plength+0],p->onshape[plength+1],p->onshape[plength+2],p->onshape[plength+3]);
    for (np=0;np<NPHASES;np++){ for (na=0;na<NSLICES;na++){ 
      printf("p->angle_%d_phase_%d_shape,max=%0.2f,min=%0.2f,mean=%0.2f,std=%0.2f\n",na,np,p->angleshape[na+np*NSLICES][plength+0],p->angleshape[na+np*NSLICES][plength+1],p->angleshape[na+np*NSLICES][plength+2],p->angleshape[na+np*NSLICES][plength+3]);
      printf("p->angle_%d_phase_%d_in,max=%0.2f,min=%0.2f,mean=%0.2f,std=%0.2f\n",na,np,p->anglein[na+np*NSLICES][plength+0],p->anglein[na+np*NSLICES][plength+1],p->anglein[na+np*NSLICES][plength+2],p->anglein[na+np*NSLICES][plength+3]);}}}
}

double lgnonshapekernel(double irad,double x,double y,double x0,double y0)
{
  /* returns the weight of lgn on-cell kernel centered at x0,y0 evaluated at x,y */
  double r2 = pow(x-x0,2)+pow(y-y0,2);
  return exp(-(r2)/(irad/2))/(PI*irad/2) - exp(-(r2)/(3*irad/2))/(PI*3*irad/2);
}

int lgndump(struct lgn *p)
{
  /* this dumps the lgn data structure to a file, writing everything in header order except for the name */
  char text[256];
  int ni=0;
  FILE * fp=NULL;
  sprintf(text,"./%s.lgn",p->name);
  if ((fp = fopen(text, "wb")) == NULL){ printf(" %% Warning: cannot open %s in lgndump\n", text); return 0;}
  fwrite(&p->rows,sizeof(int),1,fp);
  fwrite(&p->cols,sizeof(int),1,fp);
  fwrite(&p->height,sizeof(double),1,fp);
  fwrite(&p->width,sizeof(double),1,fp);
  fwrite(&p->irad,sizeof(double),1,fp);
  fwrite(&p->lgnbaserate,sizeof(double),1,fp);
  fwrite(&p->peakrate,sizeof(double),1,fp);
  //fwrite(&p->plan_forward,sizeof(fftw_plan),1,fp);
  //fwrite(&p->plan_backward,sizeof(fftw_plan),1,fp);
  fwrite(p->pnm,sizeof(double),p->rows*p->cols,fp);
  fwrite(p->onshape,sizeof(double),p->rows*p->cols+4,fp);
  fwrite(p->onin,sizeof(double),p->rows*p->cols,fp);
  //fwrite(p->onrate,sizeof(double),5*p->rows*p->cols,fp);
  //fwrite(p->offrate,sizeof(double),5*p->rows*p->cols,fp);
  for (ni=0;ni<NSLICES*NPHASES;ni++){ fwrite(p->angleshape[ni],sizeof(double),p->rows*p->cols+4,fp);}
  for (ni=0;ni<NSLICES*NPHASES;ni++){ fwrite(p->anglein[ni],sizeof(double),p->rows*p->cols+4,fp);}
  //for (ni=0;ni<NSLICES*NPHASES;ni++){ fwrite(p->anglerate[ni],sizeof(double),5*p->rows*p->cols,fp);}
  fclose(fp);
  return 1;
}

int lgnread(struct lgn *p)
{
  /* this reads the lgn data structure from a file, writing everything in header order except for the name */
  int verbose=0;
  char text[256];
  int ni=0;
  FILE * fp=NULL;
  sprintf(text,"./%s.lgn",p->name);
  if ((fp = fopen(text, "rb")) == NULL){ if (verbose){ printf(" %% cannot open %s in lgnread\n", text);} return 0;}
  fread(&p->rows,sizeof(int),1,fp);
  fread(&p->cols,sizeof(int),1,fp);
  fread(&p->height,sizeof(double),1,fp);
  fread(&p->width,sizeof(double),1,fp);
  fread(&p->irad,sizeof(double),1,fp);
  fread(&p->lgnbaserate,sizeof(double),1,fp);
  fread(&p->peakrate,sizeof(double),1,fp);
  //fread(&p->plan_forward,sizeof(fftw_plan),1,fp);
  //fread(&p->plan_backward,sizeof(fftw_plan),1,fp);ping();
  fread(p->pnm,sizeof(double),p->rows*p->cols,fp);
  fread(p->onshape,sizeof(double),p->rows*p->cols+4,fp);
  fread(p->onin,sizeof(double),p->rows*p->cols,fp);
  //fread(p->onrate,sizeof(double),5*p->rows*p->cols,fp);
  //fread(p->offrate,sizeof(double),5*p->rows*p->cols,fp);
  for (ni=0;ni<NSLICES*NPHASES;ni++){ fread(p->angleshape[ni],sizeof(double),p->rows*p->cols+4,fp);}
  for (ni=0;ni<NSLICES*NPHASES;ni++){ fread(p->anglein[ni],sizeof(double),p->rows*p->cols+4,fp);}
  //for (ni=0;ni<NSLICES*NPHASES;ni++){ fread(p->anglerate[ni],sizeof(double),5*p->rows*p->cols,fp);}
  fclose(fp);
  return 1;
}

void lgnconvolve(struct lgn *p)
{
  /* this convolves p->pnm to get p->onin */
  int verbose=0,ticker=0;
  int r=0,c=0,r2=0,c2=0,r3=0,c3=0;
  double rowsperdegree=0,colsperdegree=0;
  double weight=0,da=0,temp=0;
  rowsperdegree = (double)p->rows/(double)p->height; colsperdegree = (double)p->cols/(double)p->width;
  da = (double)1/(double)rowsperdegree/(double)colsperdegree;
  for (c=0;c<p->cols;c++){ for (r=0;r<p->rows;r++){ 
    if (verbose){ ticker++; if (ticker>100){ printf("%d%%\n",(int)floor(100*(double)(r+c*p->rows)/(double)(p->rows*p->cols))); ticker=0;}}
    temp=0;
    for (r2=-(int)floor(rowsperdegree*p->irad*3);r2<=(int)floor(rowsperdegree*p->irad*3);r2++){
      for (c2=-(int)floor(colsperdegree*p->irad*3);c2<=(int)floor(colsperdegree*p->irad*3);c2++){
	r3 = periodize(r+r2,0,p->rows); c3 = periodize(c+c2,0,p->cols);
	weight = (1/PI)*(exp(-(pow((double)r2/(double)rowsperdegree,2) + pow((double)c2/(double)colsperdegree,2))/(p->irad/2))/(p->irad/2) - exp(-(pow((double)r2/(double)rowsperdegree,2) + pow((double)c2/(double)colsperdegree,2))/(3*p->irad/2))/(3*p->irad/2));
	temp += weight*p->pnm[r3 + c3*p->rows];}}
    p->onin[r + c*p->rows] = temp*da;}}
  if (verbose){ printf("100%% -- done\n");}
}

void lgnrefresh(struct lgn *p,double t,double DT)
{
  /* this refreshes what's on the screen 
     given the key for GRATING_VS_LMI:
     -10: nothing at all
     -2: rtc input
     -1: pulsed grating
     0: full on grating
     1: line motion illusion
     2: drifting square
     3: growing bar
     4: blob motion illusion
     5: growing blob
     6: cute neuron motion illusion
     7: growing cute neuron
     8: seidemann's gabor
  */
  int frame1=0,frame2=0;
  double tcrit=0,tcycle=0;
  struct rtc *r=NULL;
  int nl=0,oldtab=0,maxback=0,dahumptab=0,gate_flag=0;
  if (INPUT_IGNORE_FLAG){
    if (strcmp(p->name,"blank_64x96")!=0){ lgnremake(p,"blank_64x96");}}
  else{
    if (GRATING_VS_LMI==-10){
      /* nothing at all */
      if (strcmp(p->name,"blank_64x96")!=0){ lgnremake(p,"blank_64x96");}}
    else if (GRATING_VS_LMI==-2){
      /* rtc input */
      if (t<=0 && t+DT>0){ lgnremake(p,"blank_64x96");}
      if (PTREE_BOTHER){ /* don't record unless */ 
	GLOBAL_PTREE->gate_flag = 0;
	if (RTC_BOTHER){
	  r = GLOBAL_RTC;
	  dahumptab = periodize((int)floor((r->nangles-1)*periodize(0 /* dahump angle */,0,PI)/PI),0,r->nangles-1);
	  gate_flag=0;
	  maxback = minimum(r->length,(int)floor((double)128/(double)r->update_every)); /* within 128ms of stimulus onset */
	  for (nl=0;nl<maxback;nl++){ oldtab = periodize(r->tab-nl,0,r->length); if (r->anglera[oldtab]==dahumptab){ gate_flag=1;}}
	  GLOBAL_PTREE->gate_flag = gate_flag;}}
      frame1 = framefig(t); frame2 = framefig(t+DT);
      if (frame2>frame1){ 
	if (rand01<1.0/(double)(RTC_BOTHER ? GLOBAL_RTC->nangles : NSLICES)){ lgnremake(p,"greyblank_64x96");}
	else{ INPUT_SPACEANGLE = rand01*PI; INPUT_SPACEPHASE = rand01*2*PI; lgnremake(p,NULL);}}}
    else if (GRATING_VS_LMI==-1){
      /* pulsed grating */
      if (t<=0 && t+DT>0){ lgnremake(p,"blank_64x96");}
      if (PTREE_BOTHER){ /* don't record unless... */ GLOBAL_PTREE->gate_flag = 0;}
      tcrit = STIMULUS_ONSET_TIME;
      if (t>tcrit){
	tcycle = (t-tcrit) - CYCLE_LENGTH*floor((t-tcrit)/CYCLE_LENGTH);
	/* now input single cycle */
	if (tcycle/CYCLE_LENGTH < 0.5){ 
	  if (PTREE_BOTHER){ /* ...within 128ms of stimulus onset */ GLOBAL_PTREE->gate_flag = (tcycle<128 ? 1 : 0);}
	  frame1=framefig(t);frame2=framefig(t+DT); if (frame2>frame1){ lgnremake(p,NULL);}}
	else{ if (strcmp(p->name,"blank_64x96")!=0){ lgnremake(p,"blank_64x96");}}}
      else{ if (strcmp(p->name,"blank_64x96")!=0){ lgnremake(p,"blank_64x96");}}}
    else if (GRATING_VS_LMI==0){
      /* full on grating */
      if (t<=0 && t+DT>0){ lgnremake(p,"blank_64x96");}
      tcrit = STIMULUS_ONSET_TIME;
      if (t>tcrit){
	frame1 = framefig(t);
	frame2 = framefig(t+DT);
	if (frame2>frame1){ 
	  INPUT_SPACEANGLE = INPUT_SPACEANGLE + 0;
	  INPUT_SPACEPHASE = periodize(INPUT_SPACEPHASE + 2*PI/GRATING_DRIFT /* *GLOBAL_SPF*1/128 */,0,2*PI);
	  lgnremake(p,NULL);}}
      else{ if (strcmp(p->name,"blank_64x96")!=0){ lgnremake(p,"blank_64x96");}}}
    else if (GRATING_VS_LMI==1){
      /* line motion illusion */
      if (t<=0 && t+DT>0){ lgnremake(p,"blank_64x96");}
      tcrit = STIMULUS_ONSET_TIME;
      if (t<=tcrit && t+DT>tcrit){ lgnremake(p,"square_64x96");}
      tcrit = STIMULUS_ONSET_TIME + SQUARE_DURATION_TIME;
      if (t<=tcrit && t+DT>tcrit){ lgnremake(p,"blank_64x96");}
      tcrit = STIMULUS_ONSET_TIME + SQUARE_DURATION_TIME + LINE_DELAY_TIME;
      if (t<=tcrit && t+DT>tcrit){ lgnremake(p,"line_64x96");}
      tcrit = STIMULUS_ONSET_TIME + SQUARE_DURATION_TIME + LINE_DELAY_TIME + SQUARE_DRAG_TIME;
      if (t<=tcrit && t+DT>tcrit){ lgnremake(p,"blank_64x96");}}
    else if (GRATING_VS_LMI==2){
      /* moving square */
      if (t<=0 && t+DT>0){ lgnremake(p,"blank_64x96");}
      tcrit = STIMULUS_ONSET_TIME;
      if (t<=tcrit && t+DT>tcrit){ lgnremake(p,"square_64x96");}
      tcrit = STIMULUS_ONSET_TIME + SQUARE_DURATION_TIME + SQUARE_DRAG_TIME*0.0/11.0;
      if (t<=tcrit && t+DT>tcrit){ lgnremake(p,"square1of8_64x96");}
      tcrit = STIMULUS_ONSET_TIME + SQUARE_DURATION_TIME + SQUARE_DRAG_TIME*1.0/11.0;
      if (t<=tcrit && t+DT>tcrit){ lgnremake(p,"square2of8_64x96");}
      tcrit = STIMULUS_ONSET_TIME + SQUARE_DURATION_TIME + SQUARE_DRAG_TIME*2.0/11.0;
      if (t<=tcrit && t+DT>tcrit){ lgnremake(p,"square3of8_64x96");}
      tcrit = STIMULUS_ONSET_TIME + SQUARE_DURATION_TIME + SQUARE_DRAG_TIME*3.0/11.0;
      if (t<=tcrit && t+DT>tcrit){ lgnremake(p,"square4of8_64x96");}
      tcrit = STIMULUS_ONSET_TIME + SQUARE_DURATION_TIME + SQUARE_DRAG_TIME*4.0/11.0;
      if (t<=tcrit && t+DT>tcrit){ lgnremake(p,"square5of8_64x96");}
      tcrit = STIMULUS_ONSET_TIME + SQUARE_DURATION_TIME + SQUARE_DRAG_TIME*5.0/11.0;
      if (t<=tcrit && t+DT>tcrit){ lgnremake(p,"square6of8_64x96");}
      tcrit = STIMULUS_ONSET_TIME + SQUARE_DURATION_TIME + SQUARE_DRAG_TIME*6.0/11.0;
      if (t<=tcrit && t+DT>tcrit){ lgnremake(p,"square7of8_64x96");}
      tcrit = STIMULUS_ONSET_TIME + SQUARE_DURATION_TIME + SQUARE_DRAG_TIME*7.0/11.0;
      if (t<=tcrit && t+DT>tcrit){ lgnremake(p,"square8of8_64x96");}
      tcrit = STIMULUS_ONSET_TIME + SQUARE_DURATION_TIME + SQUARE_DRAG_TIME*8.0/11.0;
      if (t<=tcrit && t+DT>tcrit){ lgnremake(p,"square9of8_64x96");}
      tcrit = STIMULUS_ONSET_TIME + SQUARE_DURATION_TIME + SQUARE_DRAG_TIME*9.0/11.0;
      if (t<=tcrit && t+DT>tcrit){ lgnremake(p,"square10of8_64x96");}
      tcrit = STIMULUS_ONSET_TIME + SQUARE_DURATION_TIME + SQUARE_DRAG_TIME*10.0/11.0;
      if (t<=tcrit && t+DT>tcrit){ lgnremake(p,"square11of8_64x96");}
      tcrit = STIMULUS_ONSET_TIME + SQUARE_DURATION_TIME + SQUARE_DRAG_TIME*11.0/11.0;
      if (t<=tcrit && t+DT>tcrit){ lgnremake(p,"blank_64x96");}}
    else if (GRATING_VS_LMI==3){
      /* growing bar */
      if (t<=0 && t+DT>0){ lgnremake(p,"blank_64x96");}
      tcrit = STIMULUS_ONSET_TIME;
      if (t<=tcrit && t+DT>tcrit){ lgnremake(p,"square_64x96");}
      tcrit = STIMULUS_ONSET_TIME + SQUARE_DURATION_TIME + SQUARE_DRAG_TIME*0.0/12.0;
      if (t<=tcrit && t+DT>tcrit){ lgnremake(p,"grobar1of8_64x96");}
      tcrit = STIMULUS_ONSET_TIME + SQUARE_DURATION_TIME + SQUARE_DRAG_TIME*1.0/12.0;
      if (t<=tcrit && t+DT>tcrit){ lgnremake(p,"grobar2of8_64x96");}
      tcrit = STIMULUS_ONSET_TIME + SQUARE_DURATION_TIME + SQUARE_DRAG_TIME*2.0/12.0;
      if (t<=tcrit && t+DT>tcrit){ lgnremake(p,"grobar3of8_64x96");}
      tcrit = STIMULUS_ONSET_TIME + SQUARE_DURATION_TIME + SQUARE_DRAG_TIME*3.0/12.0;
      if (t<=tcrit && t+DT>tcrit){ lgnremake(p,"grobar4of8_64x96");}
      tcrit = STIMULUS_ONSET_TIME + SQUARE_DURATION_TIME + SQUARE_DRAG_TIME*4.0/12.0;
      if (t<=tcrit && t+DT>tcrit){ lgnremake(p,"grobar5of8_64x96");}
      tcrit = STIMULUS_ONSET_TIME + SQUARE_DURATION_TIME + SQUARE_DRAG_TIME*5.0/12.0;
      if (t<=tcrit && t+DT>tcrit){ lgnremake(p,"grobar6of8_64x96");}
      tcrit = STIMULUS_ONSET_TIME + SQUARE_DURATION_TIME + SQUARE_DRAG_TIME*6.0/12.0;
      if (t<=tcrit && t+DT>tcrit){ lgnremake(p,"grobar7of8_64x96");}
      tcrit = STIMULUS_ONSET_TIME + SQUARE_DURATION_TIME + SQUARE_DRAG_TIME*7.0/12.0;
      if (t<=tcrit && t+DT>tcrit){ lgnremake(p,"grobar8of8_64x96");}
      tcrit = STIMULUS_ONSET_TIME + SQUARE_DURATION_TIME + SQUARE_DRAG_TIME*8.0/12.0;
      if (t<=tcrit && t+DT>tcrit){ lgnremake(p,"grobar9of8_64x96");}
      tcrit = STIMULUS_ONSET_TIME + SQUARE_DURATION_TIME + SQUARE_DRAG_TIME*9.0/12.0;
      if (t<=tcrit && t+DT>tcrit){ lgnremake(p,"grobar10of8_64x96");}
      tcrit = STIMULUS_ONSET_TIME + SQUARE_DURATION_TIME + SQUARE_DRAG_TIME*10.0/12.0;
      if (t<=tcrit && t+DT>tcrit){ lgnremake(p,"grobar11of8_64x96");}
      tcrit = STIMULUS_ONSET_TIME + SQUARE_DURATION_TIME + SQUARE_DRAG_TIME*11.0/12.0;
      if (t<=tcrit && t+DT>tcrit){ lgnremake(p,"grobar12of8_64x96");}
      tcrit = STIMULUS_ONSET_TIME + SQUARE_DURATION_TIME + SQUARE_DRAG_TIME*12.0/12.0;
      if (t<=tcrit && t+DT>tcrit){ lgnremake(p,"blank_64x96");}}
    else if (GRATING_VS_LMI==4){
      /* blob motion illusion */
      if (t<=0 && t+DT>0){ lgnremake(p,"blank_64x96");}
      tcrit = STIMULUS_ONSET_TIME;
      if (t<=tcrit && t+DT>tcrit){ lgnremake(p,"blob1of12_64x96");}
      tcrit = STIMULUS_ONSET_TIME + SQUARE_DURATION_TIME;
      if (t<=tcrit && t+DT>tcrit){ lgnremake(p,"blank_64x96");}
      tcrit = STIMULUS_ONSET_TIME + SQUARE_DURATION_TIME + LINE_DELAY_TIME;
      if (t<=tcrit && t+DT>tcrit){ lgnremake(p,"blob12of12_64x96");}
      tcrit = STIMULUS_ONSET_TIME + SQUARE_DURATION_TIME + LINE_DELAY_TIME + SQUARE_DRAG_TIME;
      if (t<=tcrit && t+DT>tcrit){ lgnremake(p,"blank_64x96");}}
    else if (GRATING_VS_LMI==5){
      /* growing blob */
      if (t<=0 && t+DT>0){ lgnremake(p,"blank_64x96");}
      tcrit = STIMULUS_ONSET_TIME;
      if (t<=tcrit && t+DT>tcrit){ lgnremake(p,"blob1of12_64x96");}
      tcrit = STIMULUS_ONSET_TIME + SQUARE_DURATION_TIME + SQUARE_DRAG_TIME*0.0/12.0;
      if (t<=tcrit && t+DT>tcrit){ lgnremake(p,"blob1of12_64x96");}
      tcrit = STIMULUS_ONSET_TIME + SQUARE_DURATION_TIME + SQUARE_DRAG_TIME*1.0/12.0;
      if (t<=tcrit && t+DT>tcrit){ lgnremake(p,"blob2of12_64x96");}
      tcrit = STIMULUS_ONSET_TIME + SQUARE_DURATION_TIME + SQUARE_DRAG_TIME*2.0/12.0;
      if (t<=tcrit && t+DT>tcrit){ lgnremake(p,"blob3of12_64x96");}
      tcrit = STIMULUS_ONSET_TIME + SQUARE_DURATION_TIME + SQUARE_DRAG_TIME*3.0/12.0;
      if (t<=tcrit && t+DT>tcrit){ lgnremake(p,"blob4of12_64x96");}
      tcrit = STIMULUS_ONSET_TIME + SQUARE_DURATION_TIME + SQUARE_DRAG_TIME*4.0/12.0;
      if (t<=tcrit && t+DT>tcrit){ lgnremake(p,"blob5of12_64x96");}
      tcrit = STIMULUS_ONSET_TIME + SQUARE_DURATION_TIME + SQUARE_DRAG_TIME*5.0/12.0;
      if (t<=tcrit && t+DT>tcrit){ lgnremake(p,"blob6of12_64x96");}
      tcrit = STIMULUS_ONSET_TIME + SQUARE_DURATION_TIME + SQUARE_DRAG_TIME*6.0/12.0;
      if (t<=tcrit && t+DT>tcrit){ lgnremake(p,"blob7of12_64x96");}
      tcrit = STIMULUS_ONSET_TIME + SQUARE_DURATION_TIME + SQUARE_DRAG_TIME*7.0/12.0;
      if (t<=tcrit && t+DT>tcrit){ lgnremake(p,"blob8of12_64x96");}
      tcrit = STIMULUS_ONSET_TIME + SQUARE_DURATION_TIME + SQUARE_DRAG_TIME*8.0/12.0;
      if (t<=tcrit && t+DT>tcrit){ lgnremake(p,"blob9of12_64x96");}
      tcrit = STIMULUS_ONSET_TIME + SQUARE_DURATION_TIME + SQUARE_DRAG_TIME*9.0/12.0;
      if (t<=tcrit && t+DT>tcrit){ lgnremake(p,"blob10of12_64x96");}
      tcrit = STIMULUS_ONSET_TIME + SQUARE_DURATION_TIME + SQUARE_DRAG_TIME*10.0/12.0;
      if (t<=tcrit && t+DT>tcrit){ lgnremake(p,"blob11of12_64x96");}
      tcrit = STIMULUS_ONSET_TIME + SQUARE_DURATION_TIME + SQUARE_DRAG_TIME*11.0/12.0;
      if (t<=tcrit && t+DT>tcrit){ lgnremake(p,"blob12of12_64x96");}
      tcrit = STIMULUS_ONSET_TIME + SQUARE_DURATION_TIME + SQUARE_DRAG_TIME*12.0/12.0;
      if (t<=tcrit && t+DT>tcrit){ lgnremake(p,"blank_64x96");}}
    else if (GRATING_VS_LMI==6){
      /* cute neuron motion illusion */
      if (t<=0 && t+DT>0){ lgnremake(p,"blank_64x96");}
      tcrit = STIMULUS_ONSET_TIME;
      if (t<=tcrit && t+DT>tcrit){ lgnremake(p,"neuron1of12_64x96");}
      tcrit = STIMULUS_ONSET_TIME + SQUARE_DURATION_TIME;
      if (t<=tcrit && t+DT>tcrit){ lgnremake(p,"blank_64x96");}
      tcrit = STIMULUS_ONSET_TIME + SQUARE_DURATION_TIME + LINE_DELAY_TIME;
      if (t<=tcrit && t+DT>tcrit){ lgnremake(p,"neuron12of12_64x96");}}
    else if (GRATING_VS_LMI==7){
      /* growing cute neuron */
      if (t<=0 && t+DT>0){ lgnremake(p,"blank_64x96");}
      tcrit = STIMULUS_ONSET_TIME;
      if (t<=tcrit && t+DT>tcrit){ lgnremake(p,"neuron1of12_64x96");}
      tcrit = STIMULUS_ONSET_TIME + SQUARE_DURATION_TIME + SQUARE_DRAG_TIME*0.0/12.0;
      if (t<=tcrit && t+DT>tcrit){ lgnremake(p,"neuron1of12_64x96");}
      tcrit = STIMULUS_ONSET_TIME + SQUARE_DURATION_TIME + SQUARE_DRAG_TIME*1.0/12.0;
      if (t<=tcrit && t+DT>tcrit){ lgnremake(p,"neuron2of12_64x96");}
      tcrit = STIMULUS_ONSET_TIME + SQUARE_DURATION_TIME + SQUARE_DRAG_TIME*2.0/12.0;
      if (t<=tcrit && t+DT>tcrit){ lgnremake(p,"neuron3of12_64x96");}
      tcrit = STIMULUS_ONSET_TIME + SQUARE_DURATION_TIME + SQUARE_DRAG_TIME*3.0/12.0;
      if (t<=tcrit && t+DT>tcrit){ lgnremake(p,"neuron4of12_64x96");}
      tcrit = STIMULUS_ONSET_TIME + SQUARE_DURATION_TIME + SQUARE_DRAG_TIME*4.0/12.0;
      if (t<=tcrit && t+DT>tcrit){ lgnremake(p,"neuron5of12_64x96");}
      tcrit = STIMULUS_ONSET_TIME + SQUARE_DURATION_TIME + SQUARE_DRAG_TIME*5.0/12.0;
      if (t<=tcrit && t+DT>tcrit){ lgnremake(p,"neuron6of12_64x96");}
      tcrit = STIMULUS_ONSET_TIME + SQUARE_DURATION_TIME + SQUARE_DRAG_TIME*6.0/12.0;
      if (t<=tcrit && t+DT>tcrit){ lgnremake(p,"neuron7of12_64x96");}
      tcrit = STIMULUS_ONSET_TIME + SQUARE_DURATION_TIME + SQUARE_DRAG_TIME*7.0/12.0;
      if (t<=tcrit && t+DT>tcrit){ lgnremake(p,"neuron8of12_64x96");}
      tcrit = STIMULUS_ONSET_TIME + SQUARE_DURATION_TIME + SQUARE_DRAG_TIME*8.0/12.0;
      if (t<=tcrit && t+DT>tcrit){ lgnremake(p,"neuron9of12_64x96");}
      tcrit = STIMULUS_ONSET_TIME + SQUARE_DURATION_TIME + SQUARE_DRAG_TIME*9.0/12.0;
      if (t<=tcrit && t+DT>tcrit){ lgnremake(p,"neuron10of12_64x96");}
      tcrit = STIMULUS_ONSET_TIME + SQUARE_DURATION_TIME + SQUARE_DRAG_TIME*10.0/12.0;
      if (t<=tcrit && t+DT>tcrit){ lgnremake(p,"neuron11of12_64x96");}
      tcrit = STIMULUS_ONSET_TIME + SQUARE_DURATION_TIME + SQUARE_DRAG_TIME*11.0/12.0;
      if (t<=tcrit && t+DT>tcrit){ lgnremake(p,"neuron12of12_64x96");}
      tcrit = STIMULUS_ONSET_TIME + SQUARE_DURATION_TIME + SQUARE_DRAG_TIME*12.0/12.0;
      if (t<=tcrit && t+DT>tcrit){ lgnremake(p,"blank_64x96");}}
    else if (GRATING_VS_LMI==8){
      /* seidemann's gabor */
      if (t<=0 && t+DT>0){ lgnremake(p,"blank_64x96");}
      tcrit = STIMULUS_ONSET_TIME;
      if (t<=tcrit && t+DT>tcrit){ lgnremake(p,"seidemann_gabor_64x96");}
      tcrit = STIMULUS_ONSET_TIME + SQUARE_DURATION_TIME;
      if (t<=tcrit && t+DT>tcrit){ lgnremake(p,"blank_64x96");}}}
}

void lgnevolve(struct lgn *p,double t,double DT)
{
  /* this evolves the lgn structure by a timestep DT */
  /* assuming the impulse response looks somewhat like:
     tau=30;sigma=31; t=0:.1:300; x=(exp(-t/tau)-exp(-t/sigma))/(1/sigma-1/tau)/sigma/tau; tau2=80;sigma2=tau*sigma/tau2; y=(exp(-t/tau2)-exp(-t/sigma2))/(1/sigma2-1/tau2)/sigma2/tau2; plot(t,x,t,y,t,x-y);
     anyway, we can estimate the peak by evaluating
     normalizer = (exp(-rise1/rise1)-exp(-rise1/rise2))/(1/rise2-1/rise1)/rise2/rise1 - (exp(-rise1/fall1)-exp(-rise1/fall2))/(1/fall2-1/fall1)/fall2/fall1;
     and then multiplying (x-y) via the factor p->peakrate/time_normalize */
  int verbose=0;
  int nr=0,nc=0,na=0,np=0,j=0,k=0;
  double rise1=0,rise2=0,fall1=0,fall2=0;
  double time_normalize=0;
  double e1=0,e2=0,e3=0,e4=0,tmp=0;
  if (verbose){ printf(" %% [entering lgnevolve] with t=%f,DT=%f\n",t,DT);}
  switch (LGN_TYPE_FLAG){
  case 0: /* monkey */ rise1=10;rise2=11;fall1=40;fall2=rise1*rise2/fall1;time_normalize=1.0/0.015*p->peakrate; break;
  case 1: /* cat */ rise1=15;rise2=16;fall1=64;fall2=rise1*rise2/fall1;time_normalize=1.0/0.011*p->peakrate; break;
  default: /* quick */ rise1=2;rise2=3;fall1=9;fall2=rise1*rise2/fall1;time_normalize=1.0/0.06*p->peakrate; break;}
  e1=exp(-DT/rise1);e2=exp(-DT/rise2);e3=exp(-DT/fall1);e4=exp(-DT/fall2);
  lgnrefresh(p,t,DT);
  for (nr=0;nr<p->rows;nr++){ for (nc=0;nc<p->cols;nc++){
    j = nr*5 + nc*p->rows*5;
    tmp = p->onin[nr+nc*p->rows];
    p->onrate[1+j] = e2*p->onrate[1+j] + 1.0/(rise1-rise2)*(e1-e2)*(p->onrate[0+j] - rise1*tmp) + (1-e2)*tmp;
    p->onrate[0+j] = e1*p->onrate[0+j] + rise1*(1-e1)*tmp;
    p->onrate[3+j] = e4*p->onrate[3+j] + 1.0/(fall1-fall2)*(e3-e4)*(p->onrate[2+j] - fall1*tmp) + (1-e4)*tmp;
    p->onrate[2+j] = e3*p->onrate[2+j] + fall1*(1-e3)*tmp;
    //p->onrate[4+j] = maximum(0,p->lgnbaserate + time_normalize*(p->onrate[1+j] - p->onrate[3+j]));
    p->onrate[4+j] = time_normalize*(p->onrate[1+j] - p->onrate[3+j]);
    tmp *= -1;
    p->offrate[1+j] = e2*p->offrate[1+j] + 1.0/(rise1-rise2)*(e1-e2)*(p->offrate[0+j] - rise1*tmp) + (1-e2)*tmp;
    p->offrate[0+j] = e1*p->offrate[0+j] + rise1*(1-e1)*tmp;
    p->offrate[3+j] = e4*p->offrate[3+j] + 1.0/(fall1-fall2)*(e3-e4)*(p->offrate[2+j] - fall1*tmp) + (1-e4)*tmp;
    p->offrate[2+j] = e3*p->offrate[2+j] + fall1*(1-e3)*tmp;
    //p->offrate[4+j] = maximum(0,p->lgnbaserate + time_normalize*(p->offrate[1+j] - p->offrate[3+j]));
    p->offrate[4+j] = time_normalize*(p->offrate[1+j] - p->offrate[3+j]);
    for (na=0;na<NSLICES;na++){ for (np=0;np<NPHASES;np++){
      j = nr*5 + nc*p->rows*5;
      k = na + np*NSLICES;
      tmp = p->anglein[k][nr+nc*p->rows];
      p->anglerate[k][1+j] = e2*p->anglerate[k][1+j] + 1.0/(rise1-rise2)*(e1-e2)*(p->anglerate[k][0+j] - rise1*tmp) + (1-e2)*tmp;
      p->anglerate[k][0+j] = e1*p->anglerate[k][0+j] + rise1*(1-e1)*tmp;
      p->anglerate[k][3+j] = e4*p->anglerate[k][3+j] + 1.0/(fall1-fall2)*(e3-e4)*(p->anglerate[k][2+j] - fall1*tmp) + (1-e4)*tmp;
      p->anglerate[k][2+j] = e3*p->anglerate[k][2+j] + fall1*(1-e3)*tmp;
      //p->anglerate[k][4+j] = maximum(0,p->lgnbaserate + time_normalize*(p->anglerate[k][1+j] - p->anglerate[k][3+j]));
      p->anglerate[k][4+j] = time_normalize*(p->anglerate[k][1+j] - p->anglerate[k][3+j]);}}}}
}

void lgnresetrates(struct lgn *p)
{
  /* simply sets the rates to 0 */
  int na=0,np=0;
  tfree(p->onrate); p->onrate = (double *) tcalloc(5*p->rows*p->cols,sizeof(double));
  tfree(p->offrate); p->offrate = (double *) tcalloc(5*p->rows*p->cols,sizeof(double));
  for (na=0;na<NSLICES;na++){ for (np=0;np<NPHASES;np++){ 
    tfree(p->anglerate[na+np*NSLICES]); p->anglerate[na+np*NSLICES] = (double *) tcalloc(5*p->rows*p->cols,sizeof(double));}}
}

/* void lgngabor(struct lgn *p) */
/* { */
/*   /\* this applies a simple gabor (7 adjacent circles) to the onrate and offrate to get the anglerate */
/*      ............ */
/*      ....@.O..... */
/*      ...@.O.@.... */
/*      ....O.@..... */
/*      ............ *\/ */
/*   int r=0,c=0,r2=0,c2=0,r3=0,c3=0,a=0; */
/*   double rowsperdegree=0,colsperdegree=0; */
/*   double angle=0,temp=0; */
/*   rowsperdegree = (double)p->rows/(double)p->height; colsperdegree = (double)p->cols/(double)p->width; */
/*   for (r=0;r<p->rows;r++){ for (c=0;c<p->cols;c++){ for (a=0;a<NSLICES;a++){ */
/*     angle = (double)a/(double)NSLICES*PI; */
/*     temp=0; */
/*     r2 = 0;c2 = 0; */
/*     r3 = periodize(r+r2,0,p->rows); c3 = periodize(c+c2,0,p->cols); */
/*     temp += p->onrate[4+r3*5+c3*p->rows*5]; */
/*     r2 = (int)floor(rowsperdegree*sin(angle+0.0*PI/3.0)*2*p->irad); c2 = (int)floor(colsperdegree*cos(angle+0.0*PI/3.0)*2*p->irad); */
/*     r3 = periodize(r+r2,0,p->rows); c3 = periodize(c+c2,0,p->cols); */
/*     temp += p->onrate[4+r3*5+c3*p->rows*5]; */
/*     r2 = (int)floor(rowsperdegree*sin(angle+1.0*PI/3.0)*2*p->irad); c2 = (int)floor(colsperdegree*cos(angle+1.0*PI/3.0)*2*p->irad); */
/*     r3 = periodize(r+r2,0,p->rows); c3 = periodize(c+c2,0,p->cols); */
/*     temp += -p->offrate[4+r3*5+c3*p->rows*5]; */
/*     r2 = (int)floor(rowsperdegree*sin(angle+2.0*PI/3.0)*2*p->irad); c2 = (int)floor(colsperdegree*cos(angle+2.0*PI/3.0)*2*p->irad); */
/*     r3 = periodize(r+r2,0,p->rows); c3 = periodize(c+c2,0,p->cols); */
/*     temp += -p->offrate[4+r3*5+c3*p->rows*5]; */
/*     r2 = (int)floor(rowsperdegree*sin(angle+3.0*PI/3.0)*2*p->irad); c2 = (int)floor(colsperdegree*cos(angle+3.0*PI/3.0)*2*p->irad); */
/*     r3 = periodize(r+r2,0,p->rows); c3 = periodize(c+c2,0,p->cols); */
/*     temp += p->onrate[4+r3*5+c3*p->rows*5]; */
/*     r2 = (int)floor(rowsperdegree*sin(angle+4.0*PI/3.0)*2*p->irad); c2 = (int)floor(colsperdegree*cos(angle+4.0*PI/3.0)*2*p->irad); */
/*     r3 = periodize(r+r2,0,p->rows); c3 = periodize(c+c2,0,p->cols); */
/*     temp += -p->offrate[4+r3*5+c3*p->rows*5]; */
/*     r2 = (int)floor(rowsperdegree*sin(angle+5.0*PI/3.0)*2*p->irad); c2 = (int)floor(colsperdegree*cos(angle+5.0*PI/3.0)*2*p->irad); */
/*     r3 = periodize(r+r2,0,p->rows); c3 = periodize(c+c2,0,p->cols); */
/*     temp += -p->offrate[4+r3*5+c3*p->rows*5]; */
/*     p->anglerate[a + r*NSLICES + c*NSLICES*p->rows] = temp/7.0;}}} */
/* } */

void lgntfree(struct lgn *p)
{
  int na=0,np=0;
  tfree(p->pnm);tfree(p->onin);tfree(p->onrate);tfree(p->offrate);
  for (na=0;na<NSLICES;na++){ for (np=0;np<NPHASES;np++){
    tfree(p->angleshape[na + np*NSLICES]);
    tfree(p->anglein[na + np*NSLICES]);
    tfree(p->anglerate[na + np*NSLICES]);}}
  tfree(p->angleshape);tfree(p->anglein);tfree(p->anglerate);
  if (p->plan_forward != NULL){ fftw_destroy_plan(p->plan_forward);}
  if (p->plan_backward != NULL){ fftw_destroy_plan(p->plan_backward);}
  tfree(p);p=NULL;
}

int framefig(double t){ return periodize((int)floor(GLOBAL_FPS*t/1024.0),0,GLOBAL_FPS)%GLOBAL_FPS;}
