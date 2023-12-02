/* int findrecordnumber(); */
/* int lgnpnmdump(struct lgn *,int,int,char *); */
/* int cortexdump(struct neuronarray *,char *); */
/* int cleanupoutput(); */

/* generically, dump_type satisfies:
   0: ascii
   1: colorscaled pnm or jpeg
   2: rawbits data dump
*/

int findrecordnumber()
{
  char text[64];
  int recordnumber=0,exit_flag=0,found_flag=0;
  int zero=0;
  FILE *fp=NULL;
  if (GLOBAL_RECORD_NUMBER>0){
    recordnumber=1;
    do{
      found_flag=0;
      if (!found_flag){
	sprintf(text,"./%s%drecord",GLOBAL_STRING,recordnumber);
	if ((fp=fopen(text,"r"))!=NULL){ fclose(fp); recordnumber++; exit_flag=0; found_flag=1;}}
      if (!found_flag){
	sprintf(text,"./%s%d_bundle.tar",GLOBAL_STRING,recordnumber);
	if ((fp=fopen(text,"r"))!=NULL){ fclose(fp); recordnumber++; exit_flag=0; found_flag=1;}}
      if (!found_flag){
	sprintf(text,"./%s%d_bundle.tar.gz",GLOBAL_STRING,recordnumber);
	if ((fp=fopen(text,"r"))!=NULL){ fclose(fp); recordnumber++; exit_flag=0; found_flag=1;}}}
    while (found_flag);
    sprintf(text,"./%s%drecord",GLOBAL_STRING,recordnumber);
    if ((fp=fopen(text,"w"))==NULL){ printf("findrecordnumber dump to stdout\n"); fp=stdout;} 
    fwrite(&SYSTEM_ROW_SIZE,sizeof(int),1,fp);
    fwrite(&SYSTEM_COL_SIZE,sizeof(int),1,fp);
    if (LGN_BOTHER){
      fwrite(&GLOBAL_LGN->rows,sizeof(int),1,fp);
      fwrite(&GLOBAL_LGN->cols,sizeof(int),1,fp);}
    else{
      fwrite(&zero,sizeof(int),1,fp);
      fwrite(&zero,sizeof(int),1,fp);}
    if (fp!=NULL){ fclose(fp);}}
  return recordnumber;
}

int lgnpnmdump(struct lgn *p,int na,int np)
{
  /* stores lgn data p->pnm and p->anglerate[na+np*NSLICES][4+nr*5+nc*p->rows*5] as records */
  int verbose=1;
  int output_flag=1;
  char text[64];
  int rows = p->rows,cols = p->cols;
  int nr=0,nc=0,nr2=0,nc2=0;
  double *pnmra=NULL,*angra=NULL;
  pnmra = (double *) tcalloc(rows*cols,sizeof(double));
  angra = (double *) tcalloc(rows*cols,sizeof(double));
  for (nr=0;nr<rows;nr++){ for (nc=0;nc<cols;nc++){
    nr2 = (int)floor((double)nr*((double)p->rows/(double)rows));
    nc2 = (int)floor((double)nc*((double)p->cols/(double)cols));
    pnmra[nr + nc*rows] = p->pnm[nr2 + nc2*p->rows];
    angra[nr + nc*rows] = p->anglerate[na+np*NSLICES][4+nr2*5+nc2*5*p->rows];}}
  if (GLOBAL_RECORD_NUMBER>0){
    if (verbose){ printf("recording lgn pnm and anglerate[%d+%d*NSLICES] arrays at time %0.1f+%0.1f (last time %0.1f+%0.1f)\n",na,np,GLOBAL_time,GLOBAL_DT,PNM_MOVIE_LAST_TIME,MOVIE_TIME_REFRESH);}
    sprintf(text,"./%srecord_LGNp",GLOBAL_STRING_2); output_flag *= rarecord(1,pnmra,rows*cols,text);
    sprintf(text,"./%srecord_LGNr",GLOBAL_STRING_2); output_flag *= rarecord(1,angra,rows*cols,text);}
  tfree(pnmra);tfree(angra);
  return output_flag;
}

int cortexdump(int datavsstring,struct neuronarray *Nra)
{
  /* stores data as record */
  int verbose=1;
  int output_flag=1;
  char text[64];
  if (GLOBAL_RECORD_NUMBER>0){
    if (verbose){ printf("recording VS and sN arrays at time %0.1f+%0.1f (last time %0.1f+%0.1f)\n",GLOBAL_time,GLOBAL_DT,PNM_MOVIE_LAST_TIME,MOVIE_TIME_REFRESH);}
    sprintf(text,"./%srecord_VS",GLOBAL_STRING_2); output_flag *= rarecord(datavsstring,Nra->VSra,Nra->rows*Nra->cols,text);
    sprintf(text,"./%srecord_sN",GLOBAL_STRING_2); output_flag *= rarecord(datavsstring,Nra->sNra,Nra->rows*Nra->cols,text);}
  return output_flag;
}

void connectionsdump(struct neuronarray *Nra,int dump_type,char *filename_base)
{
  /* assumes filename_base does NOT start with "./" */
  char filename[512];
  double *connections=NULL;
  int nr=0,nc=0;
  struct neuron *s=NULL,*n=NULL;
  double cA,cN,cG,mA,mN,mG;
  if (dump_type==1){
    connections = (double *) tcalloc((Nra->rows*Nra->cols)*(Nra->rows*Nra->cols+2),sizeof(double));
    for (nr=0;nr<Nra->rows*Nra->cols;nr++){ for (nc=0;nc<Nra->rows*Nra->cols;nc++){
      s = Nra->N[nr]; n = Nra->N[nc];
      if (slink(s,n,&cA,&cN,&cG,&mA,&mN,&mG)){ 
	if (*(s->t2s)<=1){ connections[nr+nc*Nra->rows*Nra->cols] = -cG;}
	else if (*(s->t2s)>1){ connections[nr+nc*Nra->rows*Nra->cols] = cA+cN;}}}}
    stats("double",connections,(Nra->rows*Nra->cols)*(Nra->rows*Nra->cols),&mA,&mG,NULL,NULL); mN=maximum(fabs(mA),fabs(mG));
    for (nr=0;nr<Nra->rows*Nra->cols;nr++){
      connections[nr+(Nra->rows*Nra->cols+1)*(Nra->rows*Nra->cols)] = -mN + 2*mN*(double)nr/(double)(Nra->rows*Nra->cols-1);}
    sprintf(filename,"./%s_%srecord.pnm",filename_base,GLOBAL_STRING_2);
    WritePNMfile_color(connections,Nra->rows*Nra->cols,Nra->rows*Nra->cols+2,mN,-mN,filename,7);
    tfree(connections);connections=NULL;}
}

int cleanupoutput()
{
  /* call this after calling everything else to clean up the output */
  int output_flag=1;
  char command[256];
  char *s=GLOBAL_STRING_2;
  printf(" %% trying to bundle everything into %s.tar.gz\n",s);
  sprintf(command,"nice -19 tar cvf ./%s_bundle.tar *%s*;",s,s); output_flag *= system(command);
  sprintf(command,"nice -19 rm *%srecord*;",s); output_flag *= system(command);
  sprintf(command,"nice -19 gzip ./%s_bundle.tar;",s); output_flag *= system(command);
  return output_flag;
}

int trialaverage(int argc,char **argv)
{
  int verbose=1;
  char *gs=GLOBAL_STRING;
  int remove_flag=PNM_REMOVE,transpose_flag=0,uselgnrows_flag=0,usecolor_flag=7;
  double framerate_indicator=0,expected_frame=0;
  char type[8];
  int framestart=0,frameend=0;
  double max=0,min=0;
  int smoother=0;
  int still_have_options=0;
  char text[64],framename[16],command[256];
  FILE *fp=NULL,**fpra=NULL;
  int rows=0,cols=0,lgnrows=0,lgncols=0,nrow=0,ncol=0;
  int recordnumber=0,exit_flag=0,nr=0,ni=0,nf=0,stillframesleft=0,readout=0,numberofframes=0;
  double *A=NULL,*B=NULL,*tempmean=NULL,*tempstdev=NULL,mean=0,stdev=0,*ra=NULL,*ra2=NULL,max2=0,min2=0;
  struct litem *l1=NULL,*l2=NULL;
  struct llist *meanL=llistmake(),*stdevL=llistmake();
  int maxperframe_flag=0;
  still_have_options = argc-2;
  sprintf(type,"VS");framestart=-1;frameend=-1;max=0;min=0;smoother=GLOBAL_SPACE_SMOOTHER;transpose_flag=0;
  while (still_have_options){
    switch(getopt(argc,argv,"VvNnPpAaS:s:E:e:M:m:O:R:r:Ttc:C:q:Q:")){
    case 'V': case 'v': sprintf(type,"VS"); break;
    case 'N': case 'n': sprintf(type,"sN"); break;
    case 'P': case 'p': sprintf(type,"LGNp"); uselgnrows_flag=1; usecolor_flag=0; break;
    case 'A': case 'a': sprintf(type,"LGNr"); uselgnrows_flag=1; usecolor_flag=0; break;
    case 'S': case 's': framestart = atoi(optarg); break;
    case 'E': case 'e': frameend = atoi(optarg); break;
    case 'M': max = atof(optarg); break;
    case 'm': min = atof(optarg); break;
    case 'O': smoother = atoi(optarg); break;
    case 'R': case 'r': remove_flag = atoi(optarg); break;
    case 'T': case 't': transpose_flag = 1; break;
    case 'C': case 'c': usecolor_flag = atoi(optarg); break;
    case 'Q': case 'q': framerate_indicator = atof(optarg); break;
    default: exit(EXIT_SUCCESS); break;}
    still_have_options--;}
  framerate_indicator = maximum(1,framerate_indicator);
  if (verbose){ printf("trialaverage with string %s,type %s,start %d,end %d,max %0.1e,min %0.1e,colorflag=%d,framerate=%0.2f;\n",gs,type,framestart,frameend,max,min,usecolor_flag,framerate_indicator);}
  sprintf(text,"./%sallrecord_%s",gs,type);
  if ((fp=fopen(text,"r"))==NULL){
    if (verbose){ printf("didn't find %s, making...\n",text);}
    recordnumber=1;
    do{
      sprintf(text,"./%s%drecord",gs,recordnumber);
      if (verbose){ printf("trying to read file %s... ",text);}
      if ((fp=fopen(text,"r"))==NULL){ printf("failed\n"); exit_flag=1;}
      else{ 
	printf("succeeded\n");
	fread(&rows,sizeof(int),1,fp); fread(&cols,sizeof(int),1,fp); 
	fread(&lgnrows,sizeof(int),1,fp); fread(&lgncols,sizeof(int),1,fp); 
	if (fp!=NULL){ fclose(fp);} recordnumber++; exit_flag=0;}}
    while (!exit_flag);
    if (verbose){ printf("read rows=%d,cols=%d,lgnrows=%d,lgncols=%d\n",rows,cols,lgnrows,lgncols);}
    if (uselgnrows_flag){ if (verbose){ printf("using lgnrows and lgncols\n");} rows = lgnrows; cols = lgncols;}
    recordnumber--;
    if (verbose){ printf("read %d records\n",recordnumber);}
    fpra = (FILE **) tmalloc(sizeof(FILE *)*recordnumber);
    for (nr=0;nr<recordnumber;nr++){ 
      sprintf(text,"./%s%drecord_%s",gs,nr+1,type); printf("setting file pointer %d to %s\n",nr,text); fpra[nr]=fopen(text,"r");}
    stillframesleft=recordnumber;nf=0;
    sprintf(text,"./%sallrecord_%s",gs,type);
    if ((fp=fopen(text,"w"))==NULL){ printf("error, couldn't create %s\n",text); exit(EXIT_FAILURE);}
    fwrite(&rows,sizeof(int),1,fp); fwrite(&cols,sizeof(int),1,fp);
    if (fp!=NULL){ fclose(fp);}
    while (stillframesleft>0){
      stillframesleft=recordnumber;
      A = (double *) tcalloc(rows*cols,sizeof(double));
      B = (double *) tcalloc(rows*cols,sizeof(double));
      if (verbose){ printf("averaging frame %d - records ",nf);}
      for (nr=0;nr<recordnumber;nr++){
	if (verbose){ printf("%d ",nr);}
	if (fpra[nr]!=NULL){ 
	  readout = fread(B,sizeof(double),rows*cols,fpra[nr]); 
	  if (readout==rows*cols){
	    for (ni=0;ni<rows*cols;ni++){ A[ni] += B[ni];}}
	else{ stillframesleft--;}}
	else{ stillframesleft--;}}
      if (verbose){ printf("\n");}
      tfree(B);
      if (stillframesleft>0){ 
	for (ni=0;ni<rows*cols;ni++){ A[ni]/=(double)stillframesleft;}
	sprintf(text,"./%sallrecord_%s",gs,type); rarecord(1,A,rows*cols,text);}
      tfree(A);
      nf++;}
    for (nr=0;nr<recordnumber;nr++){ if (fpra[nr]!=NULL){ fclose(fpra[nr]);}}
    tfree(fpra);}
  else{ 
    if (verbose){ printf("found old %s, keeping\n",text);}
    if (fp!=NULL){ fclose(fp);}}
  sprintf(text,"./%sallrecord_%s",gs,type);
  if ((fp=fopen(text,"r"))==NULL){ printf("error, %s should be here\n",text); exit(EXIT_FAILURE);} 
  fread(&rows,sizeof(int),1,fp); fread(&cols,sizeof(int),1,fp);
  B = (double *) tcalloc(rows*cols,sizeof(double)); 
  exit_flag=0;
  do{
    readout = fread(B,sizeof(double),rows*cols,fp);
    if (readout==rows*cols){
      tempmean = (double *) tmalloc(sizeof(double));
      tempstdev = (double *) tmalloc(sizeof(double));
      stats("double",B,rows*cols,NULL,NULL,tempmean,tempstdev);
      litemadd(meanL,tempmean);
      litemadd(stdevL,tempstdev);
      exit_flag=0;}
    else{ exit_flag=1;}}
  while (!exit_flag);
  tfree(B);
  if (fp!=NULL){ fclose(fp);}
  numberofframes=meanL->length; if (numberofframes==0){ printf("error, no frames read\n"); exit(EXIT_FAILURE);}
  if (framestart<0){ framestart=0;} if (frameend<0){ frameend=numberofframes;}
  framestart = maximum(0,minimum(framestart,frameend));
  frameend = minimum(numberofframes-1,maximum(framestart,frameend));
  tempmean = (double *) tcalloc(numberofframes,sizeof(double));
  tempstdev = (double *) tcalloc(numberofframes,sizeof(double));
  l1 = meanL->first; l2 = stdevL->first; nf=0;
  while (l1!=NULL && l2!=NULL){
    if (nf>=framestart && nf<=frameend){ tempmean[nf-framestart] = *(double *)l1->item; tempstdev[nf-framestart] = *(double *)l2->item;}
    l1 = l1->child; l2 = l2->child; nf++;}
  llisttfree2(meanL);llisttfree2(stdevL);
  stats("double",tempmean,numberofframes,NULL,NULL,&mean,NULL);
  stats("double",tempstdev,numberofframes,NULL,NULL,&stdev,NULL);
  tfree(tempmean);tfree(tempstdev);
  if (verbose){ printf("number of frames is %d, starting at %d, ending at %d\n",numberofframes,framestart,frameend);}
  if (verbose){ printf("found mean=%0.1e, stdev=%0.1e\n",mean,stdev);}
  if (max>min){ max2 = max; min2 = min; maxperframe_flag=0;}
  else if (max==min){ max2 = mean+STD_VIEW*stdev; min2 = mean-STD_VIEW*stdev; maxperframe_flag=0;} 
  else if (max<min){ max=0;min=0; maxperframe_flag=1;}
  if (verbose){ printf("using max=%0.1e min=%0.1e, with %s\n",max2,min2,maxperframe_flag ? "rescale per frame" : "same scale throughout");}
  sprintf(text,"./%sallrecord_%s",gs,type);
  if ((fp=fopen(text,"r"))==NULL){ printf("warning, cannot read %s in trialaverage",text);}
  A = (double *) tcalloc(rows*cols,sizeof(double));
  nf=0;exit_flag=0;expected_frame=0;
  do{
    readout = fread(A,sizeof(double),rows*cols,fp);
    if (readout==rows*cols){
      if (nf>=framestart && nf<=frameend && nf>=(int)rint(expected_frame)){
	num2frame(nf,framename);
	sprintf(text,"./%srecordframe_%s.%s.pnm",gs,type,framename);
	ra = spacesmear(A,rows,cols,smoother);
	if (maxperframe_flag){ stats("double",ra,rows*cols,NULL,NULL,&mean,&stdev); max2=mean+stdev*STD_VIEW; min2=mean-stdev*STD_VIEW;}
	if (transpose_flag){ 
	  ra2 = (double *) tcalloc(rows*cols,sizeof(double));
	  for (nrow=0;nrow<rows;nrow++){ for (ncol=0;ncol<cols;ncol++){ ra2[ncol+nrow*cols] = ra[nrow+ncol*rows];}}
	  WritePNMfile_color(ra2,cols,rows,max2,min2,text,usecolor_flag);
	  tfree(ra2);}
	else{ WritePNMfile_color(ra,rows,cols,max2,min2,text,usecolor_flag);}
	tfree(ra);
	if (remove_flag==0){ printf("converting..."); sprintf(command,"convert %s ./%s%s.eps;",text,type,framename); system(command);}
	do{ expected_frame += framerate_indicator;} while (expected_frame <= nf); expected_frame = minimum(frameend,expected_frame);
	printf("expected frame now %0.2f or %d\n",expected_frame,(int)rint(expected_frame)); }
      nf++; exit_flag=0; if (nf>frameend){ exit_flag=1;}}
    else{ exit_flag=1;}}
  while (!exit_flag);
  tfree(A);
  if (fp!=NULL){ fclose(fp);}
  if (framerate_indicator == 1){
    sprintf(text,"./%srecordframe_%s",gs,type);
    if (pnm2mpg(text,framestart,frameend)==0 && verbose){ printf("something went wrong in pnm2mpg()\n");}
    sprintf(command,"mv %srecordframe_%s.mpg %smovie_%s_s%d_e%d_m%0.4f_M%0.4f.mpg;",gs,type,gs,type,framestart,frameend,min2,max2); system(command);
    if (remove_flag){ sprintf(command,"rm *%srecordframe_%s*;",gs,type); system(command);}}
  return numberofframes;
}

/* Here are the rtc functions */

struct rtc * rtcmake(int length,int framelength,int nangles,int nphases)
{
  /* stores nangles-1 angles and one greyblank, and nphases-1 phases and one greyblank */
  int na=0;
  struct rtc *r=NULL;
  r = (struct rtc *) tmalloc(sizeof(struct rtc));
  r->length = length;
  r->nangles = nangles;
  r->nphases = nphases;
  r->update_every = framelength;
  r->total_time = r->length*r->update_every;
  r->update_last = GLOBAL_TI;
  r->tab = 0;
  r->anglera = (int *) tcalloc(r->length,sizeof(int)); for (na=0;na<r->length;na++){ r->anglera[na]=-2;} /* setup dummy variables */
  r->phasera = (int *) tcalloc(r->length,sizeof(int)); for (na=0;na<r->length;na++){ r->phasera[na]=-2;} /* setup dummy variables */
  /* stored as [ang_index + tab*NANGLES + t2s*length*NANGLES] */
  r->mra = (double *) tcalloc(r->length*r->nangles*4,sizeof(double));
  r->Vra = (double *) tcalloc(r->length*r->nangles*4,sizeof(double));
  r->sAra = (double *) tcalloc(r->length*r->nangles*4,sizeof(double));
  r->sNra = (double *) tcalloc(r->length*r->nangles*4,sizeof(double));
  r->sGra = (double *) tcalloc(r->length*r->nangles*4,sizeof(double));
  r->VSra = (double *) tcalloc(r->length*r->nangles*4,sizeof(double));
  /* stored as [ang_index + tab*NANGLES + t2s*length*NANGLES] */
  r->rtcmra = (double *) tcalloc(r->length*r->nangles*4,sizeof(double));
  r->rtcVra = (double *) tcalloc(r->length*r->nangles*4,sizeof(double));
  r->rtcsAra = (double *) tcalloc(r->length*r->nangles*4,sizeof(double));
  r->rtcsNra = (double *) tcalloc(r->length*r->nangles*4,sizeof(double));
  r->rtcsGra = (double *) tcalloc(r->length*r->nangles*4,sizeof(double));
  r->rtcVSra = (double *) tcalloc(r->length*r->nangles*4,sizeof(double));
  r->rtcra = (double *) tcalloc(r->length,sizeof(double));
  return r;
}

void rtctfree(struct rtc *r){ 
  tfree(r->anglera); tfree(r->phasera); 
  tfree(r->mra); tfree(r->Vra); tfree(r->sAra); tfree(r->sNra); tfree(r->sGra); tfree(r->VSra); 
  tfree(r->rtcmra); tfree(r->rtcVra); tfree(r->rtcsAra); tfree(r->rtcsNra); tfree(r->rtcsGra); tfree(r->rtcVSra); 
  tfree(r->rtcra);
  tfree(r); r=NULL;}

void rtcupdate(struct rtc *r,double t,double DT)
{
  /* slow and steady update of rtc *r */
  int verbose=0;
  int nr=0,nc=0,na=0,nl=0,nt2s=0,tab=0,oldtab=0,na2=0,tab2=0;
  struct neuronarray *Nra=GLOBAL_Nra;
  struct neuron *n=NULL;
  int nang = PIE_ROW_DIA*PIE_COL_DIA;
  struct lgn *p=GLOBAL_LGN;
  for (nr=0;nr<Nra->rows;nr++){ for (nc=0;nc<Nra->cols;nc++){
    n = nget(Nra,nr,nc);
    na = periodize((int)floor((r->nangles-1)*(double)n->ang/(double)nang),0,r->nangles-1);
    tab = na + r->tab*r->nangles + *(n->t2s)*r->length*r->nangles;
    r->mra[tab] += (n->spikelast == n->spiketime && n->spiketime>=t && n->spiketime<=t+DT);
    r->Vra[tab] += *(n->V)*DT;
    r->sAra[tab] += *(n->sA)*DT;
    r->sNra[tab] += *(n->sN)*DT;
    r->sGra[tab] += *(n->sG)*DT;
    r->VSra[tab] += *(n->VS)*DT;}}
  if (t >= r->update_last + r->update_every){
    for (nl=0;nl<r->length;nl++){
      oldtab = periodize(r->tab-nl,0,r->length);
      if (r->anglera[oldtab]>=0){ /* thus not dummy */
	for (na=0;na<r->nangles-1;na++){ 
	  na2 = periodize(r->anglera[oldtab]-na,0,r->nangles-1);
	  for (nt2s=0;nt2s<4;nt2s++){
	    tab = na + r->tab*r->nangles + nt2s*r->length*r->nangles;
	    tab2 = na2 + nl*r->nangles + nt2s*r->length*r->nangles;
	    r->rtcmra[tab2] += r->mra[tab];
	    r->rtcVra[tab2] += r->Vra[tab];
	    r->rtcsAra[tab2] += r->sAra[tab];
	    r->rtcsNra[tab2] += r->sNra[tab];
	    r->rtcsGra[tab2] += r->sGra[tab];
	    r->rtcVSra[tab2] += r->VSra[tab];}}
	r->rtcra[nl] += 1;}
      else if (r->anglera[oldtab]==-1){ /* thus greyblank */
	for (na=0;na<r->nangles-1;na++){
	  na2 = r->nangles-1;
	  /* fix later */
	  for (nt2s=0;nt2s<4;nt2s++){
	    tab = na + r->tab*r->nangles + nt2s*r->length*r->nangles;
	    tab2 = na2 + nl*r->nangles + nt2s*r->length*r->nangles;
	    r->rtcmra[tab2] += r->mra[tab]/(double)(r->nangles-1);
	    r->rtcVra[tab2] += r->Vra[tab]/(double)(r->nangles-1);
	    r->rtcsAra[tab2] += r->sAra[tab]/(double)(r->nangles-1);
	    r->rtcsNra[tab2] += r->sNra[tab]/(double)(r->nangles-1);
	    r->rtcsGra[tab2] += r->sGra[tab]/(double)(r->nangles-1);
	    r->rtcVSra[tab2] += r->VSra[tab]/(double)(r->nangles-1);}}
	r->rtcra[nl] += 1;}
      else /* if (r->anglera[oldtab]<=-2) */{ /* improper angle */
	if (verbose){ printf(" %% Warning, improper angle %d passed to rtcupdate\n",r->anglera[oldtab]);}}}
    r->tab += 1; r->tab = periodize(r->tab,0,r->length);
    for (na=0;na<r->nangles;na++){ for (nt2s=0;nt2s<4;nt2s++){
      tab = na + r->tab*r->nangles + nt2s*r->length*r->nangles;
      r->mra[tab] = 0; r->Vra[tab] = 0; r->sAra[tab] = 0; r->sNra[tab] = 0; r->sGra[tab] = 0; r->VSra[tab] = 0;}}
    r->update_last = maximum(t,r->update_last + r->update_every);}
  if (strncmp(p->name,"greyblank",9)==0){ r->anglera[r->tab] = -1; r->phasera[r->tab] = -1;}
  else if (strncmp(p->name,"grating",7)==0){ 
    r->anglera[r->tab] = periodize((int)floor((r->nangles-1)*periodize(INPUT_SPACEANGLE,0,PI)/PI),0,r->nangles-1);
    r->phasera[r->tab] = periodize((int)floor((r->nphases-1)*periodize(INPUT_SPACEPHASE,0,2*PI)/2/PI),0,r->nphases-1);}
  else if (strncmp(p->name,"blank",5)==0){ r->anglera[r->tab] = -1; r->phasera[r->tab] = -1;}
  else{ if (verbose){ printf(" %% Warning! improper file %s loaded into lgn for rtc\n",p->name); r->anglera[r->tab] = -2; r->phasera[r->tab] = -2;}}
}

void rtcdump(struct rtc *r,int dump_flag)
{
  /* dumps tuningcurve to text file */
  char filename[256];
  FILE *fp=NULL;
  int nt2s=0,na=0,nl=0,tab=0,tab2=0;
  double *ra=NULL;
  if (dump_flag==0){
    sprintf(filename,"./rtc_%srecord.m",GLOBAL_STRING_2);
    if ((fp=fopen(filename,"w"))==NULL){ printf("warning, cannot read %s in rtcdump",filename); fp=stdout;}
    fprintf(fp,"clear all;\n");
    for (na=0;na<r->nangles;na++){ for (nl=0;nl<r->length;nl++){ for (nt2s=0;nt2s<4;nt2s++){
      tab = na + nl*r->nangles + nt2s*r->length*r->nangles;
      tab2 = r->nangles-1 + nl*r->nangles + nt2s*r->length*r->nangles;
      if (na<r->nangles-1){ fprintf(fp,"rtcRmra(%d,%d,%d) = %0.16lf;\n",na+1,nl+1,nt2s+1,log(r->rtcmra[tab]/r->rtcmra[tab2]));}
      fprintf(fp,"rtcmra(%d,%d,%d) = %0.16lf;\n",na+1,nl+1,nt2s+1,r->rtcmra[tab]/maximum(1,r->rtcra[nl]));
      fprintf(fp,"rtcVra(%d,%d,%d) = %0.16lf;\n",na+1,nl+1,nt2s+1,r->rtcVra[tab]/maximum(1,r->rtcra[nl]));
      fprintf(fp,"rtcsAra(%d,%d,%d) = %0.16lf;\n",na+1,nl+1,nt2s+1,r->rtcsAra[tab]/maximum(1,r->rtcra[nl]));
      fprintf(fp,"rtcsNra(%d,%d,%d) = %0.16lf;\n",na+1,nl+1,nt2s+1,r->rtcsNra[tab]/maximum(1,r->rtcra[nl]));
      fprintf(fp,"rtcsGra(%d,%d,%d) = %0.16lf;\n",na+1,nl+1,nt2s+1,r->rtcsGra[tab]/maximum(1,r->rtcra[nl]));
      fprintf(fp,"rtcVSra(%d,%d,%d) = %0.16lf;\n",na+1,nl+1,nt2s+1,r->rtcVSra[tab]/maximum(1,r->rtcra[nl]));}}}
    fprintf(fp,";figure;clf;hold on;\n");
    for (nt2s=0;nt2s<4;nt2s++){
      fprintf(fp,"subplot(4,7,%d+1);hold on;imagesc(rtcRmra(:,:,%d));axis off;title('Rm t2s %d')\n",nt2s*7,nt2s+1,nt2s);
      fprintf(fp,"subplot(4,7,%d+2);hold on;imagesc(rtcmra(:,:,%d));axis off;title('m t2s %d')\n",nt2s*7,nt2s+1,nt2s);
      fprintf(fp,"subplot(4,7,%d+3);hold on;imagesc(rtcVra(:,:,%d));axis off;title('V t2s %d')\n",nt2s*7,nt2s+1,nt2s);
      fprintf(fp,"subplot(4,7,%d+4);hold on;imagesc(rtcsAra(:,:,%d));axis off;title('sA t2s %d')\n",nt2s*7,nt2s+1,nt2s);
      fprintf(fp,"subplot(4,7,%d+5);hold on;imagesc(rtcsNra(:,:,%d));axis off;title('sN t2s %d')\n",nt2s*7,nt2s+1,nt2s);
      fprintf(fp,"subplot(4,7,%d+6);hold on;imagesc(rtcsGra(:,:,%d));axis off;title('sG t2s %d')\n",nt2s*7,nt2s+1,nt2s);
      fprintf(fp,"subplot(4,7,%d+7);hold on;imagesc(rtcVSra(:,:,%d));axis off;title('VS t2s %d')\n",nt2s*7,nt2s+1,nt2s);}
    fprintf(fp,";\n");
    if (fp!=stdout){ fclose(fp);}}
  else if (dump_flag==1){ 
    for (nt2s=0;nt2s<4;nt2s++){
      sprintf(filename,"./rtc_%srecord_t2s%d_Rm.pnm",GLOBAL_STRING_2,nt2s);
      ra = (double *) tcalloc((r->nangles-1)*r->length,sizeof(double));
      for (na=0;na<r->nangles-1;na++){ for (nl=0;nl<r->length;nl++){
	tab = na + nl*r->nangles + nt2s*r->length*r->nangles;
	tab2 = r->nangles-1 + nl*r->nangles + nt2s*r->length*r->nangles;
	ra[na + nl*(r->nangles-1)] = log(r->rtcmra[tab]/r->rtcmra[tab2]);}}
      WritePNMfile_color(ra,r->nangles-1,r->length,0,0,filename,7);
      tfree(ra);ra=NULL;
      sprintf(filename,"./rtc_%srecord_t2s%d_m.pnm",GLOBAL_STRING_2,nt2s);
      ra = (double *) tcalloc(r->nangles*r->length,sizeof(double));
      for (na=0;na<r->nangles;na++){ for (nl=0;nl<r->length;nl++){
	tab = na + nl*r->nangles + nt2s*r->length*r->nangles;
	ra[na + nl*r->nangles] = r->rtcmra[tab]/maximum(1.0,r->rtcra[nl]);}}
      WritePNMfile_color(ra,r->nangles,r->length,0,0,filename,7);
      tfree(ra);ra=NULL;
      sprintf(filename,"./rtc_%srecord_t2s%d_V.pnm",GLOBAL_STRING_2,nt2s);
      ra = (double *) tcalloc(r->nangles*r->length,sizeof(double));
      for (na=0;na<r->nangles;na++){ for (nl=0;nl<r->length;nl++){
	tab = na + nl*r->nangles + nt2s*r->length*r->nangles;
	ra[na + nl*r->nangles] = r->rtcVra[tab]/maximum(1.0,r->rtcra[nl]);}}
      WritePNMfile_color(ra,r->nangles,r->length,0,0,filename,7);
      tfree(ra);ra=NULL;
      sprintf(filename,"./rtc_%srecord_t2s%d_sA.pnm",GLOBAL_STRING_2,nt2s);
      ra = (double *) tcalloc(r->nangles*r->length,sizeof(double));
      for (na=0;na<r->nangles;na++){ for (nl=0;nl<r->length;nl++){
	tab = na + nl*r->nangles + nt2s*r->length*r->nangles;
	ra[na + nl*r->nangles] = r->rtcsAra[tab]/maximum(1.0,r->rtcra[nl]);}}
      WritePNMfile_color(ra,r->nangles,r->length,0,0,filename,7);
      tfree(ra);ra=NULL;
      sprintf(filename,"./rtc_%srecord_t2s%d_sN.pnm",GLOBAL_STRING_2,nt2s);
      ra = (double *) tcalloc(r->nangles*r->length,sizeof(double));
      for (na=0;na<r->nangles;na++){ for (nl=0;nl<r->length;nl++){
	tab = na + nl*r->nangles + nt2s*r->length*r->nangles;
	ra[na + nl*r->nangles] = r->rtcsNra[tab]/maximum(1.0,r->rtcra[nl]);}}
      WritePNMfile_color(ra,r->nangles,r->length,0,0,filename,7);
      tfree(ra);ra=NULL;
      sprintf(filename,"./rtc_%srecord_t2s%d_sG.pnm",GLOBAL_STRING_2,nt2s);
      ra = (double *) tcalloc(r->nangles*r->length,sizeof(double));
      for (na=0;na<r->nangles;na++){ for (nl=0;nl<r->length;nl++){
	tab = na + nl*r->nangles + nt2s*r->length*r->nangles;
	ra[na + nl*r->nangles] = r->rtcsGra[tab]/maximum(1.0,r->rtcra[nl]);}}
      WritePNMfile_color(ra,r->nangles,r->length,0,0,filename,7);
      tfree(ra);ra=NULL;
      sprintf(filename,"./rtc_%srecord_t2s%d_VS.pnm",GLOBAL_STRING_2,nt2s);
      ra = (double *) tcalloc(r->nangles*r->length,sizeof(double));
      for (na=0;na<r->nangles;na++){ for (nl=0;nl<r->length;nl++){
	tab = na + nl*r->nangles + nt2s*r->length*r->nangles;
	ra[na + nl*r->nangles] = r->rtcVSra[tab]/maximum(1.0,r->rtcra[nl]);}}
      WritePNMfile_color(ra,r->nangles,r->length,0,0,filename,7);
      tfree(ra);ra=NULL;}}
}

/* Here are the strobe and strobetrace functions */

struct strobe * strobemake(int length,double update_timestep,int cycle_bother)
{
  struct strobe *st=NULL;
  st = (struct strobe *) tmalloc(sizeof(struct strobe));
  st->length = length;
  st->tab = 0;
  st->last_time = GLOBAL_TI;
  st->update_timestep = update_timestep;
  st->data = (double *) tcalloc(st->length,sizeof(double));
  st->cycle_bother = cycle_bother;
  st->cyclenum = 0;
  st->cycledata = NULL;
  if (st->cycle_bother){ st->cycledata = (double *) tcalloc(st->length,sizeof(double));}
  return st;
}

void strobeupdate(struct strobe *st,double t,double DT,double val)
{
  /* this averages, rather than strobes, which is better for spike statistics
     in addition, time-steps are resolved accurately */
  int i=0;
  double oldstep=0,newstep=0,nexttime=st->last_time+st->update_timestep;
  int oldtab=0,newtab=0;
  if (t+DT < nexttime){ 
    st->data[st->tab] += val*DT/st->update_timestep;}
  else /* if (t+DT >= nexttime) */{ 
    oldstep = minimum(st->update_timestep,maximum(0,nexttime-t)); newstep = minimum(st->update_timestep,maximum(0,t+DT-nexttime));
    oldtab = st->tab; 
    st->data[oldtab] += val*oldstep/st->update_timestep;
    newtab = st->tab+1; 
    if (newtab==st->length){ 
      newtab=0; 
      if (st->cycle_bother){ for (i=0;i<st->length;i++){ st->cycledata[i] += st->data[i];} st->cyclenum++;}}
    st->data[newtab] = val*newstep/st->update_timestep;
    st->tab = newtab;
    st->last_time = maximum(t+DT,st->last_time + st->update_timestep);}
}

void strobeupdate_sum(struct strobe *st,double t,double DT,double val)
{
  /* this sums, rather than strobes, which is only good for spike statistics */
  int i=0;
  st->data[st->tab] += val;
  if (t+DT >= st->last_time+st->update_timestep){
    st->tab++;
    if (st->tab >= st->length){ 
      st->tab=0; 
      if (st->cycle_bother){ for (i=0;i<st->length;i++){ st->cycledata[i] += st->data[i];} st->cyclenum++;}}
    st->data[st->tab] = 0;
    st->last_time = maximum(t+DT,st->last_time + st->update_timestep);}
}

void strobeupdate_bkp(struct strobe *st,double t,double DT,double val)
{
  /* this averages, rather than strobes, which is better for spike statistics */
  int i=0;
  st->data[st->tab] += val*DT/st->update_timestep;
  if (t+DT >= st->last_time+st->update_timestep){
    st->tab++;
    if (st->tab >= st->length){ 
      st->tab=0; 
      if (st->cycle_bother){ for (i=0;i<st->length;i++){ st->cycledata[i] += st->data[i];} st->cyclenum++;}}
    st->data[st->tab] = 0;
    st->last_time = maximum(t+DT,st->last_time + st->update_timestep);}
}

void strobeupdate_old(struct strobe *st,double t,double val)
{
  /* this strobes, as opposed to averaging, which is not good for discontinuous data like spike statistics */
  int i=0;
  if (t >= st->last_time+st->update_timestep){
    st->data[st->tab] = val; st->tab++;
    if (st->tab >= st->length){ 
      st->tab=0; 
      if (st->cycle_bother){ for (i=0;i<st->length;i++){ st->cycledata[i] += st->data[i];} st->cyclenum++;}}
    st->last_time = maximum(t,st->last_time + st->update_timestep);}
}

void stradump(struct strobe **stra,int ralength,int dump_type,char *filename_base)
{
  /* assumes filename_base does NOT start with "./" */
  int wrap_at=1024,rows=0,cols=0;
  char filename_base2[256],filename2[512];
  int sttab=stra[0]->tab,stlength=stra[0]->length,stcycle_bother=stra[0]->cycle_bother,stcyclenum=stra[0]->cyclenum;
  int na=0,nt=0,nt2=0;
  double *ra=NULL;
  FILE *fp=NULL;
  if (dump_type==0){
    sprintf(filename2,"%s.m",filename_base);
    if ((fp=fopen(filename2,"w"))==NULL){ printf(" %% Warning, couldn't open %s in stradump, writing to stdout\n",filename2); fp=stdout;}
    fprintf(fp,"clear %s_ra;\n",filename_base);
    for (na=0;na<ralength;na++){
      fprintf(fp,"%% array %d of %d;\n",na+1,ralength);
      for (nt=sttab;nt<sttab+stlength;nt++){
	nt2 = periodize(nt,0,stlength);
	fprintf(fp,"%s_ra(%d,%d)=%0.16f;\n",filename_base,na+1,nt-sttab+1,stra[na]->data[nt2]);}}
    fprintf(fp,"%% now plotting;\n");
    fprintf(fp,"figure;clf;hold on;imagesc(%s_ra);title('%s');hold off;\n",filename_base,filename_base);
    if (fp!=stdout){ fclose(fp);}
    if (stcycle_bother && stcyclenum>0){
      sprintf(filename2,"%s_cycle.m",filename_base);
      if ((fp=fopen(filename2,"w"))==NULL){ printf(" %% Warning, couldn't open %s in stradump, writing to stdout\n",filename2); fp=stdout;}
      fprintf(fp,"clear %s_cycle_ra;\n",filename_base);
	for (na=0;na<ralength;na++){
	  fprintf(fp,"%% array %d of %d;\n",na+1,ralength);
	  for (nt=0;nt<stlength;nt++){
	    fprintf(fp,"%s_cycle_ra(%d,%d)=%0.16f;\n",filename_base,na+1,nt+1,stra[na]->cycledata[nt]/stcyclenum);}}
      fprintf(fp,"%% now plotting;\n");
      fprintf(fp,"figure;clf;hold on;imagesc(%s_cycle_ra);title('%s_cycle');hold off;\n",filename_base,filename_base);
      if (fp!=stdout){ fclose(fp);}}}
  else if (dump_type==1){
    sprintf(filename_base2,"./%s",filename_base);
    if (wrap_at>0 && stlength>wrap_at){
      rows = (ralength+1)*((int)ceil((double)stlength/(double)wrap_at));
      cols = wrap_at;
      ra = (double *) tcalloc(rows*cols,sizeof(double));
      for (na=0;na<ralength;na++){ for (nt=sttab;nt<sttab+stlength;nt++){
	nt2 = periodize(nt,0,stlength);
	ra[na + ((nt-sttab)/wrap_at)*(ralength+1) + ((nt-sttab)%wrap_at)*rows] = stra[na]->data[nt2];}}
      if (ralength>1){ sprintf(filename2,"%s.pnm",filename_base2); WritePNMfile_color(ra,rows,cols,0,0,filename2,7);}
      ra2jpg(ra,"double",rows,cols,0,filename_base2,0);
      tfree(ra);}
    else{ 
      ra = (double *) tcalloc(ralength*stlength,sizeof(double));
      for (na=0;na<ralength;na++){ for (nt=sttab;nt<sttab+stlength;nt++){
	nt2 = periodize(nt,0,stlength);
	ra[na + (nt-sttab)*ralength] = stra[na]->data[nt2];}}
      if (ralength>1){ sprintf(filename2,"%s.pnm",filename_base2); WritePNMfile_color(ra,ralength,stlength,0,0,filename2,7);}
      ra2jpg(ra,"double",ralength,stlength,0,filename_base2,0);
      tfree(ra);}
    if (stcycle_bother && stcyclenum>0){
      if (wrap_at>0 && stlength>wrap_at){
	rows = (ralength+1)*((int)ceil((double)stlength/(double)wrap_at));
	cols = wrap_at;
	ra = (double *) tcalloc(rows*cols,sizeof(double));
	for (na=0;na<ralength;na++){ for (nt=0;nt<stlength;nt++){
	  ra[na + ((nt)/wrap_at)*(ralength+1) + ((nt)%wrap_at)*rows] = stra[na]->cycledata[nt]/stcyclenum;}}
	if (ralength>1){ sprintf(filename2,"%s_cycle.pnm",filename_base2); WritePNMfile_color(ra,rows,cols,0,0,filename2,7);}
	ra2jpg(ra,"double",rows,cols,0,filename_base2,0);
	tfree(ra);}
      else{ 
	ra = (double *) tcalloc(ralength*stlength,sizeof(double));
	for (na=0;na<ralength;na++){ for (nt=0;nt<stlength;nt++){
	  ra[na + (nt)*ralength] = stra[na]->cycledata[nt]/stcyclenum;}}
	if (ralength>1){ sprintf(filename2,"%s_cycle.pnm",filename_base2); WritePNMfile_color(ra,ralength,stlength,0,0,filename2,7);}
	ra2jpg(ra,"double",ralength,stlength,0,filename_base2,0);
	tfree(ra);}}}
}

void strobetfree(struct strobe *st){ tfree(st->data); if (st->cycle_bother){ tfree(st->cycledata);} tfree(st); st=NULL;}

void strobetraceupdate(struct strobetrace *st,double t,double DT)
{
  double *tempV=NULL,*tempA=NULL,*tempN=NULL,*tempG=NULL,*tempVS=NULL,*tempm=NULL,*temppm=NULL;
  int nr=0,nc=0,na=0,t2s=0,nr2=0,nc2=0,tab=0;
  struct neuron *n=NULL;
  struct neuronarray *Nra=GLOBAL_Nra;
  int rows=0,cols=0,s=0;
  int nang=PIE_ROW_DIA*PIE_COL_DIA;
  double sumsn=0,sumsnsn=0,sumvs=0,sumvsvs=0,sumsnvs=0;
  /* tack on */
  /* single neuron */
  strobeupdate_old(st->tst,t,t);
  strobeupdate(st->Vst,t,DT,*(st->n->V));
  strobeupdate(st->VSst,t,DT,*(st->n->VS));
  strobeupdate(st->sAst,t,DT,*(st->n->sA));
  strobeupdate(st->sNst,t,DT,*(st->n->sN));
  strobeupdate(st->sGst,t,DT,*(st->n->sG));
  /* angular populations */
  tempV = (double *) tcalloc(st->nangles,sizeof(double));
  tempA = (double *) tcalloc(st->nangles,sizeof(double));
  tempN = (double *) tcalloc(st->nangles,sizeof(double));
  tempG = (double *) tcalloc(st->nangles,sizeof(double));
  tempVS = (double *) tcalloc(st->nangles,sizeof(double));
  tempm = (double *) tcalloc(4,sizeof(double));
  for (nr=0;nr<Nra->rows;nr++){ for (nc=0;nc<Nra->cols;nc++){
    n = nget(Nra,nr,nc);
    na = periodize((int)floor(st->nangles*(double)n->ang/(double)nang),0,st->nangles);
    if (n->spiketime==n->spikelast && n->spiketime>=t && n->spiketime<=t+DT){ tempm[(int)*(n->t2s)] += 1;}
    tempV[na] += *(n->V);
    tempA[na] += *(n->sA);
    tempN[na] += *(n->sN);
    tempG[na] += *(n->sG);
    tempVS[na] += *(n->VS);}}
  for (na=0;na<st->nangles;na++){
    strobeupdate(st->Vstra[na],t,DT,tempV[na]/(double)st->natotal[na]);
    strobeupdate(st->VSstra[na],t,DT,tempVS[na]/(double)st->natotal[na]);
    strobeupdate(st->sAstra[na],t,DT,tempA[na]/(double)st->natotal[na]);
    strobeupdate(st->sNstra[na],t,DT,tempN[na]/(double)st->natotal[na]);
    strobeupdate(st->sGstra[na],t,DT,tempG[na]/(double)st->natotal[na]);}
  for (t2s=0;t2s<4;t2s++){ 
    strobeupdate(st->mstra[t2s],t,DT,tempm[t2s]/(double)st->t2stotal[t2s]);}
  tfree(tempV);tfree(tempA);tfree(tempN);tfree(tempG);tfree(tempVS);tfree(tempm);
  /* pcs similarity */
  rows = st->pcspiesacross*PIE_ROW_DIA;
  cols = st->pcspiestall*PIE_COL_DIA;
  s = maximum(0,GLOBAL_SPACE_SMOOTHER);
  tempVS = (double *) tcalloc((rows+2*s)*(cols+2*s),sizeof(double));  
  tempN = (double *) tcalloc((rows+2*s)*(cols+2*s),sizeof(double));  
  temppm = (double *) tcalloc(st->nangles,sizeof(double));
  for (nr=0;nr<rows+2*s;nr++){ for (nc=0;nc<cols+2*s;nc++){
    n = nget(Nra,periodize(st->n->row+nr-rows/2-s,0,Nra->rows),periodize(st->n->col+nc-cols/2-s,0,Nra->cols));
    if (/* ES,EC only */ n->type==+1 && n->spiketime==n->spikelast && n->spiketime>=t && n->spiketime<=t+DT){ 
      na = periodize((int)floor(st->nangles*(double)n->ang/(double)nang),0,st->nangles);
      temppm[na] += 1;}
    tempVS[nr+nc*(rows+2*s)] = *(n->VS);
    tempN[nr+nc*(rows+2*s)] = *(n->sN);}}
  tempV = spacesmear(tempVS,rows+2*s,cols+2*s,s);
  tfree(tempVS); tempVS = (double *) tcalloc(rows*cols,sizeof(double));
  for (nr=0;nr<rows;nr++){ for (nc=0;nc<cols;nc++){
    tempVS[nr+nc*rows] = tempV[(nr+s)+(nc+s)*(rows+2*s)];}}
  tfree(tempV);
  tempV = spacesmear(tempN,rows+2*s,cols+2*s,s);
  tfree(tempN); tempN = (double *) tcalloc(rows*cols,sizeof(double));
  for (nr=0;nr<rows;nr++){ for (nc=0;nc<cols;nc++){
    tempN[nr+nc*rows] = tempV[(nr+s)+(nc+s)*(rows+2*s)];}}
  tfree(tempV);
  for (na=0;na<st->nangles;na++){
    strobeupdate(st->patternstra_VS[na],t,DT,correlation(tempVS,&(st->pcsra[0+0*rows+na*rows*cols]),rows*cols));
    strobeupdate(st->patternstra_sN[na],t,DT,correlation(tempN,&(st->pcsra[0+0*rows+na*rows*cols]),rows*cols));
    strobeupdate(st->pmstra[na],t,DT,temppm[na]);}
  strobeupdate(st->patternc,t,DT,correlation(tempVS,tempN,rows*cols));
  tfree(tempVS);
  tfree(tempN);
  tfree(temppm);
  /* sN-VS temporal correlation */
  rows = SYSTEM_ROW_SIZE/NARBORS_TALL;
  cols = SYSTEM_COL_SIZE/NARBORS_WIDE;
  for (nr=0;nr<NARBORS_TALL;nr++){ for (nc=0;nc<NARBORS_WIDE;nc++){
    tempVS = (double *) tcalloc((rows+2*s)*(cols+2*s),sizeof(double));
    tempN = (double *) tcalloc((rows+2*s)*(cols+2*s),sizeof(double));
    for (nr2=0;nr2<rows;nr2++){ for (nc2=0;nc2<cols;nc2++){
      tempVS[nr2 + nc2*(rows+2*s)] = Nra->VSra[nr*rows+nr2 + (nc*cols+nc2)*SYSTEM_ROW_SIZE];
      tempN[nr2 + nc2*(rows+2*s)] = Nra->sNra[nr*rows+nr2 + (nc*cols+nc2)*SYSTEM_ROW_SIZE];}}
    tempV = spacesmear(tempVS,rows+2*s,cols+2*s,s);
    tfree(tempVS); tempVS=tempV;
    tempV = spacesmear(tempN,rows+2*s,cols+2*s,s);
    tfree(tempN); tempN=tempV;
    sumsn=0;sumsnsn=0;sumvs=0;sumvsvs=0;sumsnvs=0;
    for (nr2=0;nr2<rows;nr2++){ for (nc2=0;nc2<cols;nc2++){
      tab = nr2+nc2*rows;
      sumsn += tempN[tab];
      sumsnsn += tempN[tab]*tempN[tab];
      sumvs += tempVS[tab];
      sumvsvs += tempVS[tab]*tempVS[tab];
      sumsnvs += tempN[tab]*tempVS[tab];}}
    sumsn *= DT/((double)rows*cols); sumsnsn *= DT/((double)rows*cols); 
    sumvs *= DT/((double)rows*cols); sumvsvs *= DT/((double)rows*cols); 
    sumsnvs *= DT/((double)rows*cols);
    st->snvstc[0 + nr*5 + nc*NARBORS_TALL*5] = sumsn; st->snvstc[1 + nr*5 + nc*NARBORS_TALL*5] = sumsnsn;
    st->snvstc[2 + nr*5 + nc*NARBORS_TALL*5] = sumvs; st->snvstc[3 + nr*5 + nc*NARBORS_TALL*5] = sumvsvs;
    st->snvstc[4 + nr*5 + nc*NARBORS_TALL*5] = sumsnvs;
    tfree(tempVS);tfree(tempN);}}
  avalancheupdate(st->avalanche,t,DT);
  st->total_time += DT;
}

struct strobetrace * strobetracemake(struct neuron *n_center,int nangles,double timelength,double update_timestep,int cycle_bother,double avalanche_update_every,int avalanche_logdmin,int avalanche_logdmax)
{
  int verbose=0;
  int na=0,na2=0,nr=0,nc=0;
  int t2s=0;
  int nang=PIE_ROW_DIA*PIE_COL_DIA;
  struct strobetrace * st=NULL;
  struct neuron *n=NULL;
  struct neuronarray *Nra=GLOBAL_Nra;
  int i=0,j=0,rows=0,cols=0;
  double angle=0;
  if (verbose){ printf(" %% \n");}
  if (verbose){ printf(" %% [entering strobetracemake] with n (%d,%d), nangles %d, timelength %0.2f, update_timestep %0.2f cycle_bother %d avalanche_update_every %0.2f avalanche_logdmax %d\n",n_center->row,n_center->col,nangles,timelength,update_timestep,cycle_bother,avalanche_update_every,avalanche_logdmax);}
  st = (struct strobetrace *) tmalloc(sizeof(struct strobetrace));
  st->n = n_center;
  st->nangles = nangles;
  st->t2stotal = (int *) tcalloc(4,sizeof(int));
  st->natotal = (int *) tcalloc(st->nangles,sizeof(int));
  if (verbose){ printf(" %% making t2stotal,natotal\n");}
  for (nr=0;nr<Nra->rows;nr++){ for (nc=0;nc<Nra->cols;nc++){
    n = nget(Nra,nr,nc);
    na = periodize((int)floor(st->nangles*(double)n->ang/(double)nang),0,st->nangles);
    st->t2stotal[(int)*(n->t2s)] += 1;
    st->natotal[na] += 1;}}
  for (na=0;na<st->nangles;na++){ st->natotal[na] = maximum(1,st->natotal[na]);} 
  for (na=0;na<4;na++){ st->t2stotal[na] = maximum(1,st->t2stotal[na]);}
  st->timelength = maximum(0,timelength);
  st->update_timestep = minimum(1,update_timestep);
  st->length = (int)floor(st->timelength/st->update_timestep);
  st->cycle_bother = cycle_bother;
  if (verbose){ printf(" %% making tst,Vst,sAst,sNst,sGst,VSst\n");}
  st->tst = strobemake(st->length,st->update_timestep,st->cycle_bother);
  st->Vst = strobemake(st->length,st->update_timestep,st->cycle_bother);
  st->sAst = strobemake(st->length,st->update_timestep,st->cycle_bother);
  st->sNst = strobemake(st->length,st->update_timestep,st->cycle_bother);
  st->sGst = strobemake(st->length,st->update_timestep,st->cycle_bother);
  st->VSst = strobemake(st->length,st->update_timestep,st->cycle_bother);
  if (verbose){ printf(" %% making Vstra,sAstra,sNstra,sGstra,VSstra\n");}
  st->Vstra = (struct strobe **) tcalloc(st->nangles,sizeof(struct strobe *));
  st->sAstra = (struct strobe **) tcalloc(st->nangles,sizeof(struct strobe *));
  st->sNstra = (struct strobe **) tcalloc(st->nangles,sizeof(struct strobe *));
  st->sGstra = (struct strobe **) tcalloc(st->nangles,sizeof(struct strobe *));
  st->VSstra = (struct strobe **) tcalloc(st->nangles,sizeof(struct strobe *));
  for (na=0;na<st->nangles;na++){
    st->Vstra[na] = strobemake(st->length,st->update_timestep,st->cycle_bother);
    st->sAstra[na] = strobemake(st->length,st->update_timestep,st->cycle_bother);
    st->sNstra[na] = strobemake(st->length,st->update_timestep,st->cycle_bother);
    st->sGstra[na] = strobemake(st->length,st->update_timestep,st->cycle_bother);
    st->VSstra[na] = strobemake(st->length,st->update_timestep,st->cycle_bother);}
  if (verbose){ printf(" %% making mstra\n");}  
  st->mstra = (struct strobe **) tcalloc(4,sizeof(struct strobe *));
  for (t2s=0;t2s<4;t2s++){ st->mstra[t2s] = strobemake(st->length,st->update_timestep,st->cycle_bother);}
  st->pcspiesacross=minimum(4,NPIEROWS);
  st->pcspiestall=minimum(4,NPIECOLS);
  rows = st->pcspiesacross*PIE_ROW_DIA;
  cols = st->pcspiestall*PIE_COL_DIA;
  st->pcsra = (double *) tcalloc(st->nangles*rows*cols,sizeof(double));
  for (na=0;na<st->nangles;na++){ 
    angle = -PI/2 + (double)(na+0.5)*PI/st->nangles;
    for (i=0;i<rows;i++){ for (j=0;j<cols;j++){
      n = nget(Nra,periodize(st->n->row+i-rows/2,0,Nra->rows),periodize(st->n->col+j-cols/2,0,Nra->cols));
      na2 = periodize((int)floor(st->nangles*(double)n->ang/(double)nang),0,st->nangles);
      st->pcsra[i+j*rows+na*rows*cols] = exp(-pow(periodize((double)(na-na2),-st->nangles/2,st->nangles/2)/(st->nangles/4),2));}}}
  st->patternstra_VS = (struct strobe **) tcalloc(st->nangles,sizeof(struct strobe *));
  st->patternstra_sN = (struct strobe **) tcalloc(st->nangles,sizeof(struct strobe *));
  st->pmstra = (struct strobe **) tcalloc(st->nangles,sizeof(struct strobe *));
  for (na=0;na<st->nangles;na++){ 
    st->patternstra_VS[na] = strobemake(st->length,st->update_timestep,st->cycle_bother);
    st->patternstra_sN[na] = strobemake(st->length,st->update_timestep,st->cycle_bother);
    st->pmstra[na] = strobemake(st->length,st->update_timestep,st->cycle_bother);}
  st->patternc = strobemake(st->length,st->update_timestep,st->cycle_bother);
  /* stored as [ sum{sN} sum{sNsN} sum{VS} sum{VSVS} sum{sNVS} ... + arbor_row*5 + arbor_col*5*NARBORS_TALL] */
  st->snvstc = (double *) tcalloc(5*NARBORS_TALL*NARBORS_WIDE,sizeof(double));
  st->avalanche = avalanchemake(Nra,avalanche_update_every,avalanche_logdmin,avalanche_logdmax);
  st->total_time = 0;
  return st;
}

void strobetracedump(struct strobetrace *st,int dump_type)
{
  char gs2time[128],filename_base[256];
  int na=0,rows=0,cols=0;
  double *ra=NULL;
  int nr=0,nc=0;
  double sumsn=0,sumsnsn=0,sumvs=0,sumvsvs=0,sumsnvs=0,sumt=0;
  sprintf(gs2time,"%srecord_time%d",GLOBAL_STRING_2,(int)floor(GLOBAL_time));
  sprintf(filename_base,"strobetrace_%s_single_t",gs2time);
  stradump(&(st->tst),1,dump_type,filename_base);
  sprintf(filename_base,"strobetrace_%s_single_V",gs2time);
  stradump(&(st->Vst),1,dump_type,filename_base);
  sprintf(filename_base,"strobetrace_%s_single_sA",gs2time);
  stradump(&(st->sAst),1,dump_type,filename_base);
  sprintf(filename_base,"strobetrace_%s_single_sN",gs2time);
  stradump(&(st->sNst),1,dump_type,filename_base);
  sprintf(filename_base,"strobetrace_%s_single_sG",gs2time);
  stradump(&(st->sGst),1,dump_type,filename_base);
  sprintf(filename_base,"strobetrace_%s_m",gs2time);
  stradump(st->mstra,4,dump_type,filename_base);
  sprintf(filename_base,"strobetrace_%s_V",gs2time);
  stradump(st->Vstra,st->nangles,dump_type,filename_base);
  sprintf(filename_base,"strobetrace_%s_sA",gs2time);
  stradump(st->sAstra,st->nangles,dump_type,filename_base);
  sprintf(filename_base,"strobetrace_%s_sN",gs2time);
  stradump(st->sNstra,st->nangles,dump_type,filename_base);
  sprintf(filename_base,"strobetrace_%s_sG",gs2time);
  stradump(st->sGstra,st->nangles,dump_type,filename_base);
  sprintf(filename_base,"strobetrace_%s_VS",gs2time);
  stradump(st->VSstra,st->nangles,dump_type,filename_base);
  rows = st->pcspiesacross*PIE_ROW_DIA; cols = st->pcspiestall*PIE_COL_DIA;
  for (na=0;na<st->nangles;na++){
    sprintf(filename_base,"./strobetrace_%s_pcsra_%d.pnm",gs2time,na);
    WritePNMfile_color(&(st->pcsra[0+0*rows+na*rows*cols]),rows,cols,0,0,filename_base,7);}
  sprintf(filename_base,"strobetrace_%s_patternstra_VS",gs2time);
  stradump(st->patternstra_VS,st->nangles,dump_type,filename_base);
  sprintf(filename_base,"strobetrace_%s_patternstra_sN",gs2time);
  stradump(st->patternstra_sN,st->nangles,dump_type,filename_base);
  sprintf(filename_base,"strobetrace_%s_patternc",gs2time);
  stradump(&(st->patternc),1,dump_type,filename_base);
  sprintf(filename_base,"strobetrace_%s_pm",gs2time);
  stradump(st->pmstra,st->nangles,dump_type,filename_base);
  sprintf(filename_base,"./strobetrace_%s_snvstc.pnm",gs2time);
  rows = NARBORS_TALL; cols = NARBORS_WIDE;
  ra = (double *) tcalloc(rows*cols,sizeof(double));
  for (nr=0;nr<rows;nr++){ for (nc=0;nc<cols;nc++){
    sumsn = st->snvstc[0 + nr*5 + nc*rows*5]; sumsnsn = st->snvstc[1 + nr*5 + nc*rows*5];
    sumvs = st->snvstc[2 + nr*5 + nc*rows*5]; sumvsvs = st->snvstc[3 + nr*5 + nc*rows*5];
    sumsnvs = st->snvstc[4 + nr*5 + nc*rows*5]; sumt = st->total_time;
    ra[nr + nc*rows] = (sumt*sumsnvs - sumsn*sumvs)/sqrt((sumt*sumsnsn-sumsn*sumsn)*(sumt*sumvsvs-sumvs*sumvs));}}
  WritePNMfile_color(ra,rows,cols,+1,-1,filename_base,7);
  tfree(ra);
  sprintf(filename_base,"strobetrace_%s_avalanche",gs2time);
  avalanchedump(st->avalanche,dump_type,filename_base);
}

void strobetracetfree(struct strobetrace *st)
{
  int na=0;
  tfree(st->t2stotal);
  tfree(st->natotal);
  strobetfree(st->tst);
  strobetfree(st->Vst);
  strobetfree(st->VSst);
  strobetfree(st->sAst);
  strobetfree(st->sNst);
  strobetfree(st->sGst);
  for (na=0;na<st->nangles;na++){
    strobetfree(st->Vstra[na]);
    strobetfree(st->VSstra[na]);
    strobetfree(st->sAstra[na]);
    strobetfree(st->sNstra[na]);
    strobetfree(st->sGstra[na]);}
  tfree(st->Vstra);tfree(st->VSstra);tfree(st->sAstra);tfree(st->sNstra);tfree(st->sGstra);
  for (na=0;na<4;na++){
    strobetfree(st->mstra[na]);}
  tfree(st->mstra);
  for (na=0;na<st->nangles;na++){
    strobetfree(st->patternstra_VS[na]);
    strobetfree(st->patternstra_sN[na]);
    strobetfree(st->pmstra[na]);}
  strobetfree(st->patternc);
  tfree(st->patternstra_VS);
  tfree(st->patternstra_sN);
  tfree(st->pmstra);
  tfree(st->snvstc);
  avalanchetfree(st->avalanche);
  tfree(st); st=NULL;
}

/* Here are the steady state tuning curve functions */

struct tuningcurve * tuningcurvemake(int nangles,int nradius)
{
  struct tuningcurve *T=NULL;
  T = (struct tuningcurve *) tmalloc(sizeof(struct tuningcurve));
  T->nangles = nangles;
  T->nradius = nradius;
  T->mra = (double *) tcalloc(T->nangles*T->nradius*4,sizeof(double));
  T->Vra = (double *) tcalloc(T->nangles*T->nradius*4,sizeof(double));
  T->sAra = (double *) tcalloc(T->nangles*T->nradius*4,sizeof(double));
  T->sNra = (double *) tcalloc(T->nangles*T->nradius*4,sizeof(double));
  T->sGra = (double *) tcalloc(T->nangles*T->nradius*4,sizeof(double));
  T->VSra = (double *) tcalloc(T->nangles*T->nradius*4,sizeof(double));
  return T;
}

void tuningcurvetfree(struct tuningcurve *T)
{
  tfree(T->mra);tfree(T->Vra);tfree(T->sAra);tfree(T->sNra);tfree(T->sGra);tfree(T->VSra);tfree(T);T=NULL;
}

void tuningcurveupdate(struct tuningcurve *T,double t,double DT)
{
  struct neuronarray *Nra=GLOBAL_Nra;
  int nang=PIE_ROW_DIA*PIE_COL_DIA;
  int nr=0,nc=0,na=0,nd=0,tab=0;
  struct neuron *n=NULL;
  for (nr=0;nr<Nra->rows;nr++){ for (nc=0;nc<Nra->cols;nc++){
    n = nget(Nra,nr,nc);
    na = periodize((int)floor(T->nangles*(double)n->ang/(double)nang),0,T->nangles);
    nd = periodize((int)floor(T->nradius*(double)n->rad),0,T->nradius);
    tab = na + nd*T->nangles + *(n->t2s)*T->nangles*T->nradius;
    T->mra[tab] += (n->spiketime==n->spikelast && n->spiketime>=t && n->spiketime<=t+DT);
    T->Vra[tab] += *(n->V)*DT;
    T->sAra[tab] += *(n->sA)*DT;
    T->sNra[tab] += *(n->sN)*DT;
    T->sGra[tab] += *(n->sG)*DT;
    T->VSra[tab] += *(n->VS)*DT;}}
}

void tuningcurvedump(struct tuningcurve *T,int dump_type)
{
  /* dumps tuningcurve to text file */
  char filename[256],chartemp[16];
  FILE *fp=NULL;
  int na=0,nr=0,nt2s=0,tab=0;
  if (dump_type==0){
    sprintf(filename,"./tuningcurve_%srecord_data",GLOBAL_STRING_2);
    if ((fp=fopen(filename,"w"))==NULL){ printf("warning, cannot read %s in tuningcurvedump",filename); fp=stdout;}
    for (nt2s=0;nt2s<4;nt2s++){ for (nr=0;nr<T->nradius;nr++){ for (na=0;na<T->nangles;na++){ 
      tab = na + nr*T->nangles + nt2s*T->nangles*T->nradius;
      fprintf(fp,"%0.6lf %0.6lf %0.6lf %0.6lf %0.6lf %0.6lf\n",T->mra[tab],T->Vra[tab],T->sAra[tab],T->sNra[tab],T->sGra[tab],T->VSra[tab]);}}}
    if (fp!=stdout){ fclose(fp);}
    sprintf(filename,"./tuningcurve_%srecord.m",GLOBAL_STRING_2);
    if ((fp=fopen(filename,"w"))==NULL){ printf("warning, cannot read %s in tuningcurvedump",filename); fp=stdout;}
    fprintf(fp,"[mra,Vra,sAra,sNra,sGra,VSra] = textread('tuningcurve_%srecord_data','%%n%%n%%n%%n%%n%%n');\n",GLOBAL_STRING_2);
    fprintf(fp,"nangles=1:%d;\n",T->nangles);
    fprintf(fp,"nradius=%d;\n",T->nradius);
    fprintf(fp,"figure;clf;hold on;;\n");
    for (nr=0;nr<T->nradius;nr++){ for (nt2s=0;nt2s<4;nt2s++){
      switch(nt2s){ 
      case 0: sprintf(chartemp,"'c-'"); break;
      case 1: sprintf(chartemp,"'b-'"); break;
      case 2: sprintf(chartemp,"'m-'"); break;
      case 3: sprintf(chartemp,"'r-'"); break;
      default: sprintf(chartemp,"'k-'");break;}
      fprintf(fp,"subplot(nradius,6,%d);hold on;\n",1+nr*6);
      fprintf(fp,"ratemp=mra(%d+nangles);plot(nangles,ratemp/max(ratemp),%s);\n",nr*T->nangles + nt2s*T->nangles*T->nradius,chartemp);
      fprintf(fp,"axis([1 %d 0 1]);\n",T->nangles);
      fprintf(fp,"subplot(nradius,6,%d);hold on;\n",2+nr*6);
      fprintf(fp,"ratemp=Vra(%d+nangles);plot(nangles,(ratemp-min(ratemp))/(max(ratemp)-min(ratemp)),%s);\n",nr*T->nangles + nt2s*T->nangles*T->nradius,chartemp);
      fprintf(fp,"axis([1 %d 0 1]);\n",T->nangles);
      fprintf(fp,"subplot(nradius,6,%d);hold on;\n",3+nr*6);
      fprintf(fp,"ratemp=sAra(%d+nangles);plot(nangles,ratemp/max(ratemp),%s);\n",nr*T->nangles + nt2s*T->nangles*T->nradius,chartemp);
      fprintf(fp,"axis([1 %d 0 1]);\n",T->nangles);
      fprintf(fp,"subplot(nradius,6,%d);hold on;\n",4+nr*6);
      fprintf(fp,"ratemp=sNra(%d+nangles);plot(nangles,ratemp/max(ratemp),%s);\n",nr*T->nangles + nt2s*T->nangles*T->nradius,chartemp);
      fprintf(fp,"axis([1 %d 0 1]);\n",T->nangles);
      fprintf(fp,"subplot(nradius,6,%d);hold on;\n",5+nr*6);
      fprintf(fp,"ratemp=sGra(%d+nangles);plot(nangles,ratemp/max(ratemp),%s);\n",nr*T->nangles + nt2s*T->nangles*T->nradius,chartemp);
      fprintf(fp,"axis([1 %d 0 1]);\n",T->nangles);
      fprintf(fp,"subplot(nradius,6,%d);hold on;\n",6+nr*6);
      fprintf(fp,"ratemp=VSra(%d+nangles);plot(nangles,(ratemp-min(ratemp))/(max(ratemp)-min(ratemp)),%s);\n",nr*T->nangles + nt2s*T->nangles*T->nradius,chartemp);
      fprintf(fp,"axis([1 %d 0 1]);\n",T->nangles);}}
    if (fp!=stdout){fclose(fp);}}
  else if (dump_type==1){
    for (nt2s=0;nt2s<4;nt2s++){
      tab = 0 + 0*T->nangles + nt2s*T->nangles*T->nradius;
      sprintf(filename,"./tuningcurve_%srecord_t2s%d_m.pnm",GLOBAL_STRING_2,nt2s);
      WritePNMfile_color(&(T->mra[tab]),T->nangles,T->nradius,0,0,filename,7);
      sprintf(filename,"./tuningcurve_%srecord_t2s%d_m",GLOBAL_STRING_2,nt2s);
      ra2jpg(&(T->mra[tab]),"double",T->nangles,T->nradius,1,filename,0);
      sprintf(filename,"./tuningcurve_%srecord_t2s%d_V.pnm",GLOBAL_STRING_2,nt2s);
      WritePNMfile_color(&(T->Vra[tab]),T->nangles,T->nradius,0,0,filename,7);
      sprintf(filename,"./tuningcurve_%srecord_t2s%d_V",GLOBAL_STRING_2,nt2s);
      ra2jpg(&(T->Vra[tab]),"double",T->nangles,T->nradius,1,filename,0);
      sprintf(filename,"./tuningcurve_%srecord_t2s%d_sA.pnm",GLOBAL_STRING_2,nt2s);
      WritePNMfile_color(&(T->sAra[tab]),T->nangles,T->nradius,0,0,filename,7);
      sprintf(filename,"./tuningcurve_%srecord_t2s%d_sA",GLOBAL_STRING_2,nt2s);
      ra2jpg(&(T->sAra[tab]),"double",T->nangles,T->nradius,1,filename,0);
      sprintf(filename,"./tuningcurve_%srecord_t2s%d_sN.pnm",GLOBAL_STRING_2,nt2s);
      WritePNMfile_color(&(T->sNra[tab]),T->nangles,T->nradius,0,0,filename,7);
      sprintf(filename,"./tuningcurve_%srecord_t2s%d_sN",GLOBAL_STRING_2,nt2s);
      ra2jpg(&(T->sNra[tab]),"double",T->nangles,T->nradius,1,filename,0);
      sprintf(filename,"./tuningcurve_%srecord_t2s%d_sG.pnm",GLOBAL_STRING_2,nt2s);
      WritePNMfile_color(&(T->sGra[tab]),T->nangles,T->nradius,0,0,filename,7);
      sprintf(filename,"./tuningcurve_%srecord_t2s%d_sG",GLOBAL_STRING_2,nt2s);
      ra2jpg(&(T->sGra[tab]),"double",T->nangles,T->nradius,1,filename,0);
      sprintf(filename,"./tuningcurve_%srecord_t2s%d_VS.pnm",GLOBAL_STRING_2,nt2s);
      WritePNMfile_color(&(T->VSra[tab]),T->nangles,T->nradius,0,0,filename,7);
      sprintf(filename,"./tuningcurve_%srecord_t2s%d_VS",GLOBAL_STRING_2,nt2s);
      ra2jpg(&(T->VSra[tab]),"double",T->nangles,T->nradius,1,filename,0);}}
}

/* Here are the lmitri functions */

struct lmitri * lmitrimake(int time_start,int time_length,double row_max,double row_min)
{
  int verbose=0;
  struct lmitri *lt=NULL;
  struct neuronarray *Nra=GLOBAL_Nra;
  int nr=0,nc=0,nt2s=0;
  struct neuron *n=NULL;
  if (verbose){ printf(" %% [entering lmitrimake] time_start %d, time_length %d, row_max %f, row_min %f\n",time_start,time_length,row_max,row_min);}
  lt = (struct lmitri *) tmalloc(sizeof(struct lmitri));
  lt->time_start = time_start;
  lt->time_length = time_length;
  lt->space_length = Nra->cols;
  lt->rmax = crop((int) floor(row_max*(double)Nra->rows),0,Nra->rows-1);
  lt->rmin = crop((int) floor(row_min*(double)Nra->rows),0,Nra->rows-1);
  if (verbose){ printf(" %% rmax %d, rmin %d out of %d\n",lt->rmax,lt->rmin,Nra->rows);}
  lt->t2stotal = (double *) tcalloc(lt->space_length*5,sizeof(double));
  for (nc=0;nc<Nra->cols;nc++){ for(nr=lt->rmin;nr<=lt->rmax;nr++){ 
    n = nget(Nra,nr,nc); 
    lt->t2stotal[nc + (int)*(n->t2s)*lt->space_length] += 1;
    lt->t2stotal[nc + 4*lt->space_length] += 1;}}
  for (nc=0;nc<lt->space_length;nc++){ for (nt2s=0;nt2s<5;nt2s++){ 
    lt->t2stotal[nc+nt2s*lt->space_length] = maximum(1,lt->t2stotal[nc+nt2s*lt->space_length]);}}
  lt->nupdates = (double *) tcalloc(lt->time_length,sizeof(double));
  lt->mra = (double *) tcalloc(lt->space_length*lt->time_length*5,sizeof(double));
  lt->Vra = (double *) tcalloc(lt->space_length*lt->time_length*5,sizeof(double));
  lt->sAra = (double *) tcalloc(lt->space_length*lt->time_length*5,sizeof(double));
  lt->sNra = (double *) tcalloc(lt->space_length*lt->time_length*5,sizeof(double));
  lt->sGra = (double *) tcalloc(lt->space_length*lt->time_length*5,sizeof(double));
  lt->VSra = (double *) tcalloc(lt->space_length*lt->time_length*5,sizeof(double));
  if (verbose){ printf(" %% [exiting lmitrimake]\n");}
  return lt;
}

void lmitriupdate(struct lmitri *lt,double t,double DT)
{
  int verbose=0;
  int nr=0,nc=0,nt2s=0;
  double *summ=NULL,*sumV=NULL,*sumsA=NULL,*sumsN=NULL,*sumsG=NULL,*sumVS=NULL;
  struct neuronarray *Nra=GLOBAL_Nra;
  struct neuron *n=NULL;
  int timetab1=0,timetab2=0;
  double DT1=0,DT2=0;
  int tab=0;
  double norm_number=0;
  if (verbose){ printf(" %% [entering lmitriupdate] \n");}
  timetab1 = (int)floor(t-lt->time_start);
  timetab2 = (int)floor(t+DT-lt->time_start);
  if (timetab1==timetab2){ DT1=DT/2; DT2=DT/2;}
  else if (timetab2>timetab1){ DT2 = t+DT-lt->time_start-timetab2; DT1 = DT-DT2;}
  else{ printf(" %% Warning, timetab2>timetab1 in lmitriupdate\n");}
  if (verbose){ printf(" %% t=%0.2f, DT=%0.2f timetab1=%d, timetab2=%d, DT1=%0.2f, DT2=%0.2f\n",t,DT,timetab1,timetab2,DT1,DT2);}
  if (timetab2>=0 && timetab1<lt->time_length){
    for (nc=0;nc<Nra->cols;nc++){
      if (verbose>1){ printf(" %% column %d\n",nc);}
      summ = (double *) tcalloc(5,sizeof(double));
      sumV = (double *) tcalloc(5,sizeof(double));
      sumsA = (double *) tcalloc(5,sizeof(double));
      sumsN = (double *) tcalloc(5,sizeof(double));
      sumsG = (double *) tcalloc(5,sizeof(double));
      sumVS = (double *) tcalloc(5,sizeof(double));
      for (nr=lt->rmin;nr<=lt->rmax;nr++){ 
	if (verbose>1){ printf(" %% %% %% row %d\n",nr);}
	n=nget(Nra,nr,nc);
	if (n->spiketime==n->spikelast && n->spiketime>=t && n->spiketime<=t+DT){ summ[(int)*(n->t2s)] += 1; summ[4] += 1;}
	sumV[(int)*(n->t2s)] += *(n->V);
	sumsA[(int)*(n->t2s)] += *(n->sA);
	sumsN[(int)*(n->t2s)] += *(n->sN);
	sumsG[(int)*(n->t2s)] += *(n->sG);
	sumVS[(int)*(n->t2s)] += *(n->VS);
	sumV[4] += *(n->V);
	sumsA[4] += *(n->sA);
	sumsN[4] += *(n->sN);
	sumsG[4] += *(n->sG);
	sumVS[4] += *(n->VS);}
      for (nt2s=0;nt2s<5;nt2s++){
	if (verbose>1){ printf(" %% %% nt2s %d\n",nt2s);}
	norm_number = (double) lt->t2stotal[nc + nt2s*lt->space_length];
	if (timetab1 >= 0 && timetab1 < lt->time_length){
	  tab = nc + timetab1*lt->space_length + nt2s*lt->space_length*lt->time_length;
	  lt->mra[tab] += summ[nt2s]/norm_number;
	  lt->Vra[tab] += sumV[nt2s]*DT1/norm_number;
	  lt->sAra[tab] += sumsA[nt2s]*DT1/norm_number;
	  lt->sNra[tab] += sumsN[nt2s]*DT1/norm_number;
	  lt->sGra[tab] += sumsG[nt2s]*DT1/norm_number;
	  lt->VSra[tab] += sumVS[nt2s]*DT1/norm_number;}
	if (timetab2 >= 0 && timetab2 < lt->time_length){
	  tab = nc + timetab2*lt->space_length + nt2s*lt->space_length*lt->time_length;
	  lt->mra[tab] += summ[nt2s]/norm_number;
	  lt->Vra[tab] += sumV[nt2s]*DT2/norm_number;
	  lt->sAra[tab] += sumsA[nt2s]*DT2/norm_number;
	  lt->sNra[tab] += sumsN[nt2s]*DT2/norm_number;
	  lt->sGra[tab] += sumsG[nt2s]*DT2/norm_number;
	  lt->VSra[tab] += sumVS[nt2s]*DT2/norm_number;}}
      tfree(summ);tfree(sumV);tfree(sumsA);tfree(sumsN);tfree(sumsG);tfree(sumVS);}
    if (timetab1 >= 0 && timetab1 < lt->time_length){ lt->nupdates[timetab1] += DT1;}
    if (timetab2 >= 0 && timetab2 < lt->time_length){ lt->nupdates[timetab2] += DT2;}}
}

void lmitridump(struct lmitri *lt,int dump_type)
{
  char gs2time[128],filename[256];
  int nt2s=0,tab=0;
  sprintf(gs2time,"%srecord_time%d",GLOBAL_STRING_2,(int)floor(GLOBAL_time));
  if (dump_type==0){}
  else if (dump_type==1){
    for (nt2s=0;nt2s<5;nt2s++){
      tab = 0 + 0*lt->space_length + nt2s*lt->space_length*lt->time_length;
      sprintf(filename,"./lmitri_%s_t2s%d_m.pnm",gs2time,nt2s);
      WritePNMfile_color(&(lt->mra[tab]),lt->space_length,lt->time_length,0,0,filename,7);
      sprintf(filename,"./lmitri_%s_t2s%d_V.pnm",gs2time,nt2s);
      WritePNMfile_color(&(lt->Vra[tab]),lt->space_length,lt->time_length,0,0,filename,7);
      sprintf(filename,"./lmitri_%s_t2s%d_sA.pnm",gs2time,nt2s);
      WritePNMfile_color(&(lt->sAra[tab]),lt->space_length,lt->time_length,0,0,filename,7);
      sprintf(filename,"./lmitri_%s_t2s%d_sN.pnm",gs2time,nt2s);
      WritePNMfile_color(&(lt->sNra[tab]),lt->space_length,lt->time_length,0,0,filename,7);
      sprintf(filename,"./lmitri_%s_t2s%d_sG.pnm",gs2time,nt2s);
      WritePNMfile_color(&(lt->sGra[tab]),lt->space_length,lt->time_length,0,0,filename,7);
      sprintf(filename,"./lmitri_%s_t2s%d_VS.pnm",gs2time,nt2s);
      WritePNMfile_color(&(lt->VSra[tab]),lt->space_length,lt->time_length,0,0,filename,7);}}
}

void lmitritfree(struct lmitri *lt)
{
  tfree(lt->t2stotal); tfree(lt->nupdates);
  tfree(lt->mra); tfree(lt->Vra); tfree(lt->sAra); tfree(lt->sNra); tfree(lt->sGra); tfree(lt->VSra);
  tfree(lt); lt=NULL;
}

/* Here are the ptree functions anew */

void regionramake(struct ptree *p,double event_within,int event_threshold,int region_type)
{
  /* intended for orientation domains */
  struct neuronarray *Nra=GLOBAL_Nra;
  int nr=0,nc=0,rbin=0,pbin=0,tbin=0,abin=0;
  struct neuron *n=NULL;
  int nang=PIE_ROW_DIA*PIE_COL_DIA;
  for (nr=0;nr<p->nregions;nr++){ 
    p->regionra[nr] = (struct region *) tmalloc(sizeof(struct region)); 
    p->regionra[nr]->label = nr; /* critical to have label nr match location in regionra */
    p->regionra[nr]->last_event = GLOBAL_TI-rand01;
    p->regionra[nr]->event_within = event_within;
    p->regionra[nr]->event_threshold = event_threshold;
    p->regionra[nr]->neuronllist = llistmake();
    p->regionra[nr]->pn = NULL;}
  switch (region_type){
  case 0: /* do nothing */ break;
  case 1: /* individual neurons by index */
    for (nr=0;nr<p->nregions;nr++){
      litemadd(p->regionra[nr]->neuronllist,Nra->N[minimum(p->regionra[nr]->label,Nra->rows*Nra->cols-1)]);}
    break;
  case 2: /* angular subdivision excitatory neurons within a single pinwheel */
    for (nr=0;nr<PIE_ROW_DIA;nr++){ for (nc=0;nc<PIE_COL_DIA;nc++){
      n=nget(Nra,nr,nc);
      if (*(n->t2s)==2 || *(n->t2s)==3){
	rbin = periodize((int)floor(p->nregions*(double)n->ang/(double)nang),0,p->nregions);
	litemadd(p->regionra[rbin]->neuronllist,n);}}}
    break;
  case 3: /* angular and t2s subdivisions within a single pinwheel, requires nregions be a multiple of nt2s */
    for (nr=0;nr<PIE_ROW_DIA;nr++){ for (nc=0;nc<PIE_COL_DIA;nc++){
      n=nget(Nra,nr,nc);
      abin = periodize((int)floor((double)p->nregions/(double)4*(double)n->ang/(double)nang),0,p->nregions/4);
      tbin = (int)*(n->t2s);
      rbin = tbin + abin*2;
      litemadd(p->regionra[rbin]->neuronllist,n);}}
  case 4: /* angular and type subdivisions within two adjacent pinwheels, requires nregions to be a multiple of 2*ntypes=4 */
    for (nr=0;nr<PIE_ROW_DIA;nr++){ for (nc=0;nc<2*PIE_COL_DIA;nc++){
      n = nget(Nra,nr,nc);
      abin = periodize((int)floor((double)p->nregions/(double)4*(double)n->ang/(double)nang),0,p->nregions/4);
      tbin = (*(n->t2s) > 1);
      pbin = nc/PIE_COL_DIA;
      rbin = tbin + abin*2 + pbin*(p->nregions/2);
      litemadd(p->regionra[rbin]->neuronllist,n);}}
    for (abin=0;abin<p->nregions/4;abin++){ for (tbin=0;tbin<2;tbin++){ for (pbin=0;pbin<2;pbin++){
      rbin = tbin + abin*2 + pbin*(p->nregions/2);
      p->regionra[rbin]->event_threshold = (tbin==0?event_threshold:1);}}}
  default: break;}
}

void regionratfree(struct ptree *p)
{
  int nr=0;
  for (nr=0;nr<p->nregions;nr++){ 
    llisttfree(p->regionra[nr]->neuronllist); p->regionra[nr]->neuronllist=NULL;
    tfree(p->regionra[nr]); p->regionra[nr]=NULL;}
  tfree(p->regionra); p->regionra=NULL;
}

int region_has_event(struct region *r,double t,double DT)
{
  /* if we have threshold excitatory firing events since last_event, and these threshold are within within ms of each other */
  int verbose=0;
  int threshold=0,depth=0;
  double spikelast=0,within=0;
  struct neuron *n=NULL;
  struct llist *L=r->neuronllist;
  struct litem *l=NULL;
  int exit_flag=0;
  if (verbose){ printf(" %% [entering region_has_event] with region->label %d time %0.2f and DT %0.2f\n",r->label,t,DT);}
  threshold=minimum(L->length,r->event_threshold);
  within=r->event_within;
  if (verbose){ printf(" %% event_threshold %d event_within %0.2f\n",threshold,within);}
  if (L->length>0){ llistsort(L->first,L->last,L->length,&spikelast_compare);}
  exit_flag=0; depth=0;
  if (L->last!=NULL){
    l=L->last;
    n=(struct neuron *)L->last->item; spikelast = n->spikelast;
    if (n->spikelast>r->last_event && n->spikelast>=t && n->spikelast<=t+DT){
      if (verbose){ printf(" %% at final neuron (%d,%d,%d) with spikelast %f\n",n->row,n->col,n->ang,n->spikelast);}
      l=l->parent; depth=1;
      while (l!=NULL && !exit_flag){
	n=(struct neuron *)l->item;
	if (verbose){ printf(" %% at parent neuron (%d,%d,%d) with spikelast %f -- depth=%d\n",n->row,n->col,n->ang,n->spikelast,depth);}
	if (n->spikelast>r->last_event && (spikelast-n->spikelast)<within){ l=l->parent; depth+=1;}
	else{ exit_flag=1;}}}
    if (depth>=threshold){
      if (verbose){ printf(" %% depth=%d, event happened\n",depth);}
      r->last_event = spikelast;
      return 1;}
    else{ 
      if (verbose){ printf(" %% depth=%d, event didn't happened\n",depth);}
      return 0;}}
  return 0;
}

int region_has_event_old(struct region *r,double t,double DT)
{
  int threshold = r->neuronllist->length/2;
  struct litem *l=NULL;
  struct neuron *n=NULL;
  int nspikes = 0;
  double last_event = r->last_event;
  l=r->neuronllist->first;
  while (l!=NULL){
    n=(struct neuron *)l->item;
    if (n->spikelast==n->spiketime && n->spiketime>=t && n->spiketime <= t+DT){
      nspikes += 1;
      last_event = maximum(last_event,n->spiketime);}
    l=l->child;}
  if (nspikes>threshold){
    r->last_event = last_event;
    return 1;}
  return 0;
}

int yggdrasil_region_has_event(struct region *r,struct ptree *p,double t,double DT)
{
  /* checks to see if r->pn (as embedded in observed ptree *p) is a current event in p->eventra,
     assuming that p->eventra has already been updated for the current time-step */
  int verbose=0;
  struct pnode *pn=NULL;
  struct llist *rL=NULL;
  int tab_last=periodize(p->tab-1,0,p->length),tab_back=0,tab=0,nt=0;
  int events_found=0,event_found_local=0;
  rL = llistmake();pn=r->pn;
  while (pn!=NULL){
    litemadd(rL,pn->region);
    pn=pn->parent;}
  pn=r->pn;
  if (pn!=NULL){
    if (llitemaddorfind(0,p->eventra[tab_last],pn->region,&region2region_compare_label)!=NULL){ events_found=1; pn=pn->parent;}
    else{ events_found=0;}}
  else{ events_found=0;}
  tab_back=1;
  while (pn!=NULL && events_found){
    event_found_local=0;
    for (nt=0;nt<p->legtime;nt++){
      tab = periodize(tab_last-tab_back-nt,0,p->length);
      event_found_local += (llitemaddorfind(0,p->eventra[tab],pn->region,&region2region_compare_label)!=NULL);}
    if (event_found_local){ events_found += 1; pn=pn->parent;}
    else{ events_found=0;}
    tab_back += p->legtime;}
  if (verbose && events_found){ printf(" %% events_found %d, number required %d\n",events_found,rL->length);}
  llisttfree(rL);rL=NULL;
  if (events_found){ r->last_event = t-p->legtime; return events_found;}
  return 0;
}

struct pnode * pnodemake(struct pnode *parent,struct region *region,double weight,double relevance)
{
  struct pnode *p=NULL;
  p = (struct pnode *) tmalloc(sizeof(struct pnode));
  p->parent = parent;
  p->region = region;
  p->weight = weight;
  p->relevance = relevance;
  p->childllitem = llitemmake();
  p->broodsize = 1;
  p->temp = NULL;
  return p;
}

void pnodetfree(void *vp)
{
  /* recursively frees pnode *p and its children */
  int verbose=0;
  struct pnode *p = (struct pnode *) vp;
  if (verbose){ 
    printf(" trying to free pnode %d with parent %d, region %d, and childllitem...\n",(int)p,(int)p->parent,p->region->label);
    llitemprintf(p->childllitem,&pnodeprintf_label);}
  llitemtfree(p->childllitem,&pnodetfree); 
  p->childllitem=NULL;
  tfree(p);p=NULL;
}

struct ptree * ptreemake(int nregions,double event_within,int event_threshold,int region_type,int nlegs,int legtime)
{
  /* update to use region data structure */
  struct ptree *p=NULL;
  int nr=0;
  p = (struct ptree *) tmalloc(sizeof(struct ptree));
  p->nregions = nregions;
  p->regionra = (struct region **) tcalloc(p->nregions,sizeof(struct region*));
  regionramake(p,event_within,event_threshold,region_type);
  p->nlegs=nlegs;
  p->legtime=legtime;
  p->length = p->nlegs*p->legtime+1;
  p->eventra = (struct llitem **) tmalloc(p->length*sizeof(struct llitem *));
  for (nr=0;nr<p->length;nr++){ p->eventra[nr]=llitemmake();}
  p->tab=0;
  p->update_every=1;
  p->update_last=GLOBAL_TI;
  p->total_time=0;
  p->gate_flag=1;
  ptreemakepospre(p);
  p->wh = histmake(1,+1,-1);
  p->rh = histmake(1,+1,-1);
  return p;
}

void ptreemakepospre(struct ptree *p)
{
  /* assumes p->pretree and p->postree are both NULL */
  int nr=0;
  p->pretree = llitemmake();
  for (nr=0;nr<p->nregions;nr++){
    if (llitemaddorfind(0,p->pretree,p->regionra[nr],&region2region_compare_label)!=NULL){ printf(" %% duplicate region label in ptreemakepospre\n");}
    else{ llitemaddorfind(1,p->pretree,pnodemake(NULL,p->regionra[nr],0,0),&pnode2pnode_compare_label);}}
  llitembalance(p->pretree); p->pretree = llitemclimb(p->pretree);
  p->postree = llitemmake();
  for (nr=0;nr<p->nregions;nr++){
    if (llitemaddorfind(0,p->postree,p->regionra[nr],&region2region_compare_label)!=NULL){ printf(" %% duplicate region label in ptreemakepospre\n");}
    else{ llitemaddorfind(1,p->postree,pnodemake(NULL,p->regionra[nr],0,0),&pnode2pnode_compare_label);}}
  llitembalance(p->postree); p->postree = llitemclimb(p->postree);
}

void ptreereset(struct ptree *p)
{
  /* destroys postree and pretree, but retains regionra and eventra */
  p->total_time=0;
  llitemtfree(p->pretree,&pnodetfree); p->pretree=NULL;
  llitemtfree(p->postree,&pnodetfree); p->postree=NULL;
  ptreemakepospre(p);
}

void ptreetfree(struct ptree *p)
{
  int nr=0;
  llitemtfree(p->pretree,&pnodetfree); p->pretree=NULL;
  llitemtfree(p->postree,&pnodetfree); p->postree=NULL;
  for (nr=0;nr<p->length;nr++){ llitemtfree(p->eventra[nr],NULL); p->eventra[nr]=NULL;}
  tfree(p->eventra); p->eventra=NULL;
  regionratfree(p); p->regionra=NULL;
  histtfree(p->wh); p->wh=NULL;
  histtfree(p->rh); p->rh=NULL;
  tfree(p); p=NULL;
}

int region2region_compare_last_event(void *v1,void *v2)
{
  struct region *r1 = (struct region *)v1;
  struct region *r2 = (struct region *)v2;
  return (r1->last_event > r2->last_event ? 1 : r1->last_event < r2->last_event ? -1 : 0);
}

int region2region_compare_label(void *v1,void *v2)
{
  struct region *r1 = (struct region *)v1;
  struct region *r2 = (struct region *)v2;
  return (r1->label > r2->label ? 1 : r1->label < r2->label ? -1 : 0);
}

int region2pnode_compare_label(void *v1,void *v2)
{
  struct region *r1 = (struct region *)v1;
  struct pnode *p2 = (struct pnode *)v2;
  return (r1->label > p2->region->label ? 1 : r1->label < p2->region->label ? -1 : 0);
}

int pnode2pnode_compare_label(void *v1,void *v2)
{
  struct pnode *p1 = (struct pnode *)v1;
  struct pnode *p2 = (struct pnode *)v2;
  return (p1->region->label > p2->region->label ? 1 : p1->region->label < p2->region->label ? -1 : 0);
}

int pnode2pnode_compare_relevance(void *v1,void *v2)
{
  struct pnode *p1 = (struct pnode *)v1;
  struct pnode *p2 = (struct pnode *)v2;
  return (p1->relevance > p2->relevance ? 1 : p1->relevance < p2->relevance ? -1 : 0);
}

int pnode2pnode_compare_weight(void *v1,void *v2)
{
  struct pnode *p1 = (struct pnode *)v1;
  struct pnode *p2 = (struct pnode *)v2;
  return (p1->weight > p2->weight ? 1 : p1->weight < p2->weight ? -1 : 0);
}

void ptreeupdate_helper(struct ptree *p,double t,double DT,double weight,int verbose)
{
  /* call if t>= p->update-last + p->update_every */
  int nt=0,tab_back=0,exit_flag=0,tab=0;
  struct llist *LL=NULL,*L2=NULL;
  struct litem *l=NULL,*l2=NULL;
  struct llist *l0=NULL;
  struct region *r=NULL;
  if (verbose){ printf(" %% [calling ptreeupdate_helper] \n");}
  assert(t >= p->update_last + p->update_every);
  if (llitemlength(p->eventra[p->tab])>0){
    if (verbose){ printf(" %% have %d events, continuing\n",llitemlength(p->eventra[p->tab]));}
    if (verbose>1){ printf(" %% creating LL\n");}
    LL = llistmake();
    if (verbose>1){ printf(" %% creating l0\n");}
    l0 = llistmake();
    if (verbose>1){ printf(" %% growing l0 for tab %d\n",p->tab);}
    llistgrowllitem(l0,p->eventra[p->tab]);
    if (verbose>1){ printf(" %% pruning and sorting l0 \n");}
    llistprune(l0); llistsort(l0->first,l0->last,l0->length,&region2region_compare_last_event);
    if (verbose>1){ printf(" %% adding l0 to LL\n");}
    litemadd(LL,l0);
    if (verbose>1){ printf(" %% now adding further to LL\n");}
    tab_back=1;exit_flag=(LL->length==p->nlegs+1);
    while (!exit_flag){
      l0 = llistmake();
      for (nt=0;nt<p->legtime;nt++){
	tab = periodize(p->tab-tab_back-nt,0,p->length);
	if (verbose>5){ printf(" %% growing l0 for tab %d\n",tab);}
	llistgrowllitem(l0,p->eventra[tab]);}
      if (l0->length>0){
	if (verbose>5){ printf(" %% pruning and sorting l0\n");}
	llistprune(l0); llistsort(l0->first,l0->last,l0->length,&region2region_compare_last_event);
	if (verbose>5){ printf(" now adding further to LL\n");}
	litemadd(LL,l0);
	if (verbose>1){ printf(" %% old tab_back=%d . . .",tab_back);}
	tab_back += p->legtime;
	if (verbose>1){ printf(" new tab_back=%d\n",tab_back);}
	exit_flag = LL->length==p->nlegs+1;
	if (verbose>1){ printf(" exit_flag set to %d\n",exit_flag);}}
      else /* if (l0->length==0) */{
	if (verbose>5){ printf(" %% no events between tabs %d and %d\n",p->tab-tab_back,p->tab-tab_back-nt);}
	llisttfree(l0); l0=NULL;
	exit_flag = 1;}}
    if (verbose>5){
      printf(" %% at this point LL is composed of %d region llists:\n",LL->length);
      l = LL->first;
      while (l!=NULL){
	L2 = (struct llist *)l->item;
	printf(" %% %% llist %d of length %d\n",(int)L2,L2->length);
	l2 = L2->first;
	while (l2!=NULL){
	  r = (struct region *) l2->item;
	  printf(" %% %% %% region %d with label %d and last_event %0.3f\n",(int)r,r->label,r->last_event);
	  l2=l2->child;}
	l=l->child;}}
    if (verbose>6){ printf(" %% about to start pstrengthen_starter, ptree is:\n"); pnodeprintf(NULL,p->postree,-1,0);}
    pstrengthen_starter(p,LL,weight);
    if (verbose>6){ printf(" %% just finished pstrengthen_starter, ptree is:\n"); pnodeprintf(NULL,p->postree,-1,0);}}
  else{
    if (verbose){ printf(" %% have %d events, not updating\n",llitemlength(p->eventra[p->tab]));}}
  p->tab += 1; p->tab = periodize(p->tab,0,p->length);
  llitemtfree(p->eventra[p->tab],NULL); p->eventra[p->tab] = llitemmake();
  p->update_last = maximum(t,p->update_last + p->update_every);
}

void ptreeupdate(struct ptree *p,double t,double DT,double weight)
{
  /* updates ptree *p */
  int verbose=0;
  int nr=0;
  if (verbose){ printf(" %% [starting ptreeupdate] with p=%d, t=%0.3f, DT=%0.3f, weight=%0.3f\n",(int)p,t,DT,weight);}
  for (nr=0;nr<p->nregions;nr++){
    if (p->gate_flag && region_has_event(p->regionra[nr],t,DT) /* determine activity and set last_event_flag */){ 
      if (llitemaddorfind(0,p->eventra[p->tab],p->regionra[nr],&region2region_compare_label)!=NULL){ 
	if (verbose>1){ printf(" %% region %d active, but already in llitem\n",nr);}}
      else{
	if (verbose>1){ printf(" %% region %d active, adding to llitem\n",nr);}
	llitemaddorfind(1,p->eventra[p->tab],p->regionra[nr],&region2region_compare_label);}}}
  p->total_time += DT;
  if (t >= p->update_last + p->update_every){ ptreeupdate_helper(p,t,DT,weight,verbose);}
}

void pstrengthen_starter(struct ptree *p,struct llist *LL,double weight)
{
  /* starting recursive dropdown, and might as well free LL while we are at it */
  int verbose=0;
  struct litem *l=NULL,*l2=NULL;
  struct region *r=NULL;
  struct llist *L=NULL;
  if (verbose){ 
    printf(" %% [starting pstrengthen_starter] p=%d, LL->length=%d\n",(int)p,LL->length);
    printf(" %% LL is composed of llists of regions:\n");
    l = LL->first;
    while (l!=NULL){
      L = (struct llist *) l->item;
      printf(" %% %% length %d llist\n",L->length);
      l2 = L->first;
      while (l2!=NULL){
	r = (struct region *) l2->item;
	printf(" %% %% %% region %d\n",r->label);
	l2=l2->child;}
      l=l->child;}}
  while (LL->length>=1){
    pstrengthen_helper(1,p,NULL,p->postree,LL,LL->length,weight);
    pstrengthen_helper(0,p,NULL,p->pretree,LL,LL->length,weight);
    L = (struct llist *) LL->last->item;
    if (verbose){ printf(" %% trying to free final llist L=%d\n",(int)L);}
    litemminus(LL,L);
    llisttfree(L);L=NULL;
    if (verbose){ printf(" %% now LL->length=%d\n",LL->length);}}
  if (verbose){ printf(" %% now we have exhausted all the llists in LL (now of length %d)\n",LL->length);}
  llisttfree(LL); LL=NULL;
}

void pstrengthen_helper(int postorpre,struct ptree *p,struct pnode *parent,struct llitem *childllitem,struct llist *LL,int length,double weight)
{
  /* given a llitem of pnode2pnode_compare_label sorted pnodes, and a length-long llist *LL of region2region_compare_last_event
     sorted llists of regional events */
  int verbose=0;
  struct llist *L=NULL;
  struct litem *l=NULL,*l2=NULL;
  struct region *r=NULL,*r2=NULL;
  struct llitem *l0=NULL,*l02=NULL;
  struct pnode *pn=NULL,*pn2=NULL;
  int depth=0;
  if (verbose){ printf(" %% [starting pstrengthen_helper] with postorpre=%d, parent %d, childllitem %d of length %d, LL %d of length %d, length %d, weight %0.3f\n",postorpre,(int)parent,(int)childllitem,llitemlength(childllitem),(int)LL,LL->length,length,weight);}
  assert(LL->length>=length);
  if (parent!=NULL){ assert(parent->childllitem==childllitem);}
  if (postorpre==1){ /* postree - read length from front */
    l=LL->first; depth=1; while(l!=NULL && depth<length){ depth+=1; l=l->child;} L = (struct llist *) l->item;}
  else /* if (postorpre==0) */ { /* pretree - read length from back */
    l=LL->last; depth=1; while(l!=NULL && depth<length){ depth+=1; l=l->parent;} L = (struct llist *) l->item;}
  if (verbose){ printf(" %% read out llist L=%d of length %d\n",(int)L,L->length);}
  if (length==1){
    if (parent==NULL && LL->length==1){ /* singleton llist, must loop through causal events */
      if (postorpre==0){ /* pretree, starting at the end and heading back */
	if (verbose){ printf(" %% starting at the end and heading back\n");}
	l=L->last;
	while (l!=NULL){
	  r = (struct region *) l->item;
	  if (verbose>5){ printf(" %% at region %d\n",r->label);}
	  if ((l0=llitemaddorfind(0,childllitem,r,region2pnode_compare_label))!=NULL){
	    pn = (struct pnode *) l0->item;
	    pn->weight += weight;
	    if (verbose){ printf(" %% now at region %d out of %d\n",r->label,L->length);}
	    if (p->nlegs>0){
	      l2=l->parent;
	      while (l2!=NULL){
		r2 = (struct region *) l2->item;
		if ((l02=llitemaddorfind(0,pn->childllitem,r2,region2pnode_compare_label))!=NULL){
		  pn2 = (struct pnode *) l02->item;
		  pn2->weight += weight;}
		else{
		  pn2 = pnodemake(pn,r2,weight,0);
		  llitemaddorfind(1,pn->childllitem,pn2,pnode2pnode_compare_label);}
		l2 = l2->parent;}}}
	  else{ printf(" %% something's fishy in pstrengthen_helper\n");}
	  l = l->parent;}}
      else /* if (postorpre==1) */{ /* postree, starting at the beginning and heading forward */
	if (verbose){ printf(" %% starting at the beginning and heading forward\n");}
	l=L->first;
	while (l!=NULL){
	  r = (struct region *) l->item;
	  if (verbose>5){ printf(" %% at region %d\n",r->label);}
	  if ((l0=llitemaddorfind(0,childllitem,r,region2pnode_compare_label))!=NULL){
	    pn = (struct pnode *) l0->item;
	    pn->weight += weight;
	    if (verbose){ printf(" %% now at region %d out of %d\n",r->label,L->length);}
	    if (p->nlegs>0){
	      l2=l->child;
	      while (l2!=NULL){
		r2 = (struct region *) l2->item;
		if ((l02=llitemaddorfind(0,pn->childllitem,r2,region2pnode_compare_label))!=NULL){
		  pn2 = (struct pnode *) l02->item;
		  pn2->weight += weight;}
		else{
		  pn2 = pnodemake(pn,r2,weight,0);
		  llitemaddorfind(1,pn->childllitem,pn2,pnode2pnode_compare_label);}
		l2 = l2->child;}}}
	  else{ printf(" %% something's fishy in pstrengthen_helper\n");}
	  l = l->child;}}}
    else /* if (parent!=NULL || LL->length>1) */{ /* end of long llist, simply update */
      l=L->first;
      while (l!=NULL){
	r = (struct region *) l->item;
	if ((l0=llitemaddorfind(0,childllitem,r,region2pnode_compare_label))!=NULL){ /* already exists */
	  pn = (struct pnode *) l0->item;
	  pn->weight += weight;}
	else{ /* must make */
	  pn = pnodemake(parent,r,weight,0);
	  llitemaddorfind(1,childllitem,pn,pnode2pnode_compare_label);}
	l=l->child;}
      /* this is the most suspicious step */
      if (parent!=NULL && llitemlength(parent->childllitem)>2*parent->broodsize){
	llitembalance(parent->childllitem); parent->childllitem = llitemclimb(parent->childllitem); parent->broodsize *= 2;}}}
  else if (length>1){
    l=L->first;
    while (l!=NULL){
      r = (struct region *) l->item;
      if ((l0=llitemaddorfind(0,childllitem,r,region2pnode_compare_label))!=NULL){ /* descend */
	pn = (struct pnode *) l0->item;
	pstrengthen_helper(postorpre,p,pn,pn->childllitem,LL,length-1,weight);}
      else{ /* do nothing */ }
      l=l->child;}}
}

void pnodeprintf(struct pnode *parent,struct llitem *l0,int maxlevel,int level)
{
  /* recursively prints out the pnodes in llitem *l0 */
  char text[128];
  int nl=0;
  struct llist *L=llistmake();
  struct litem *l=NULL;
  struct pnode *pn2=NULL;
  llistgrowllitem(L,l0);
  if (level==0){ sprintf(text,"|__");}
  else{ sprintf(text,"|  "); for (nl=0;nl<level-1;nl++){ sprintf(text,"%s   ",text);} sprintf(text,"%s|__",text);}
  if (parent!=NULL){
    assert(parent->childllitem==l0);}
  //if (L->length>0){ printf(" %%%s %d->childllist of length %d contains:\n",text,(int)parent,L->length);}
  l = L->first;
  while (l!=NULL){
    pn2 = (struct pnode *) l->item;
    printf(" %%%s pnode %d->parent %d, region %d->label %d, weight %0.1f relevance %0.1f",text,(int)pn2,(int)pn2->parent,(int)pn2->region,(int)pn2->region->label,pn2->weight,pn2->relevance);
    if (llitemlength(pn2->childllitem)>0 && (maxlevel==-1 || level<maxlevel)){
      printf(", now descending into %d children...\n",llitemlength(pn2->childllitem));
      pnodeprintf(pn2,pn2->childllitem,maxlevel,level+1);}
    else{ printf("\n");}
    l=l->child;}
  llisttfree(L);L=NULL;
}

void pnodeprune_starter(int posorpre,struct ptree *p,struct pnode *parent,struct llitem *childllitem,struct llist *L)
{
  /* starts pruning process by descending p[re,os]tree and passing leaves (allong with [de,as]cendent llist *L) to pnodeprune_helper */
  int firstentry=0;
  struct pnode *pn=NULL;
  if (parent!=NULL && childllitem->parent==NULL){ /* first descent into parent->childllitem */ 
    assert(parent->childllitem==childllitem); firstentry=1;}
  if (firstentry){ litemadd(L,parent->region);}
  if (firstentry){ pnodeprune_helper(posorpre,p,parent,L);}
  if (childllitem->kidl!=NULL){ pnodeprune_starter(posorpre,p,parent,childllitem->kidl,L);}
  if (childllitem->item!=NULL){ pn=(struct pnode *)childllitem->item; pnodeprune_starter(posorpre,p,pn,pn->childllitem,L);}
  if (childllitem->kidr!=NULL){ pnodeprune_starter(posorpre,p,parent,childllitem->kidr,L);}
  if (firstentry){ llistkilllast(L);}
}

void pnodeprune_helper(int posorpre,struct ptree *p,struct pnode *pn,struct llist *L)
{
  /* given a llist of regions *L and a leaf node *pn (of depth L->length) */
  int verbose=0;
  struct llitem *childllitem1=NULL,*childllitem2=NULL,*childllitem3=NULL;
  struct litem *l=NULL;
  struct llitem *l0=NULL;
  struct region *r=NULL;
  struct pnode *pn2=NULL;
  int depth=0;
  double weight_expected=0,weight_first=0,weight_second=0,weight_third=0;
  double *ra1=NULL,*ra2=NULL,*ra3=NULL;
  double temp=0,relevance=0;
  if (verbose){
    if (posorpre){ childllitem1=p->postree;} else{ childllitem1=p->pretree;}
    l=L->first; depth=0; 
    while (l!=NULL){
      r = (struct region *) l->item;
      if ((l0=llitemaddorfind(0,childllitem1,r,region2pnode_compare_label))!=NULL){
	pn2 = (struct pnode *) l0->item;
	childllitem1 = pn2->childllitem;
	l=l->child;
	depth++;}
      else{ l=NULL;}}
    if (depth==L->length && pn2==pn){ printf(" %% passed correct node && lineage to pnodeprune_helper\n");}}
  if (L->length==1){
    weight_expected = pn->weight;
    pn->relevance = pnode_shear(pn->weight,weight_expected);}
  else if (L->length==2){ /* single leg */
    if (posorpre){ childllitem1=p->postree;} else{ childllitem1=p->pretree;}
    if (verbose){ printf(" %% %% only one leg\n");}
    r = (struct region *) L->first->item; weight_first=0;
    if ((l0=llitemaddorfind(0,childllitem1,r,region2pnode_compare_label))!=NULL){
      pn2 = (struct pnode *) l0->item;
      weight_first = pn2->weight;}
    r = (struct region *) L->last->item; weight_second=0;
    if ((l0=llitemaddorfind(0,childllitem1,r,region2pnode_compare_label))!=NULL){
      pn2 = (struct pnode *) l0->item;
      weight_second = pn2->weight;}
    weight_expected = weight_first*p->legtime*weight_second/p->total_time;
    pn->relevance = pnode_shear(pn->weight,weight_expected);}
  else if (L->length>2){ /* multiple legs */
    if (posorpre){ childllitem1=p->postree; childllitem2=p->pretree;} 
    else{ childllitem1=p->pretree; childllitem2=p->postree;}
    childllitem3=childllitem1;
    ra1 = (double *) tcalloc(L->length,sizeof(double));
    ra2 = (double *) tcalloc(L->length,sizeof(double));
    ra3 = (double *) tcalloc(L->length,sizeof(double));
    l=L->first; depth=0; 
    while (l!=NULL){
      r = (struct region *) l->item;
      if ((l0=llitemaddorfind(0,childllitem1,r,region2pnode_compare_label))!=NULL){
	pn2 = (struct pnode *) l0->item;
	ra1[depth] = pn2->weight;
	childllitem1 = pn2->childllitem;}
      else{ childllitem1=NULL;}
      l=l->child; depth++;}
    l=L->last; depth=0;
    while (l!=NULL){
      r = (struct region *) l->item;
      if ((l0=llitemaddorfind(0,childllitem2,r,region2pnode_compare_label))!=NULL){
	pn2 = (struct pnode *) l0->item;
	ra2[depth] = pn2->weight;
	childllitem2 = pn2->childllitem;}
      else{ childllitem2=NULL;}
      l=l->parent; depth++;}
    l=L->first; depth=0;
    while (l!=NULL){
      r = (struct region *) l->item;
      if ((l0=llitemaddorfind(0,childllitem3,r,region2pnode_compare_label))!=NULL){
	pn2 = (struct pnode *) l0->item;
	ra3[depth] = pn2->weight;}
      l=l->child;
      depth++;}
    relevance=0;
    for (depth=1;depth<L->length-1;depth++){
      weight_first = ra1[depth];
      weight_second = ra2[L->length-1-depth];
      weight_third = ra3[depth];
      weight_expected = weight_first*weight_second/weight_third;
      temp = pnode_shear(pn->weight,weight_expected);
      if (fabs(temp)>fabs(relevance)){ relevance=temp;}}
    pn->relevance = relevance;
    tfree(ra1);tfree(ra2);tfree(ra3);}
}

double pnode_shear(double pn_weight,double weight_expected)
{
  /* some sort of shear for pruning trees */
  double weight_tol=10;
  double relevance=0;
  if (pn_weight<weight_tol && weight_expected<weight_tol){ relevance=0;}
  else{ relevance = log(pn_weight/weight_expected);}
  if (!finite(relevance)){ relevance=0;}
  //if (!finite(relevance)){ if (pn_weight==0){ relevance=-10;} else if (weight_expected==0){ relevance=10;} else{ relevance=0;}}
  return relevance;
}

void pnodeZZZ_starter(int posorpre,struct ptree *p,struct pnode *parent,struct llitem *childllitem,double threshold,int depth)
{
  /* starts ZZZing process by descending p[re,os]tree and passing leaves (allong with [de,as]cendent llist *L) to pnodeZZZ_helper */
  int firstentry=0;
  struct pnode *pn=NULL;
  if (parent!=NULL && childllitem->parent==NULL){ /* first descent into parent->childllitem */ 
    assert(parent->childllitem==childllitem); firstentry=1;}
  if (firstentry && depth<p->nlegs){ pnodeZZZ_helper(posorpre,p,parent,threshold);}
  if (childllitem->kidl!=NULL){ pnodeZZZ_starter(posorpre,p,parent,childllitem->kidl,threshold,depth+1);}
  if (childllitem->item!=NULL){ pn=(struct pnode *)childllitem->item; pnodeZZZ_starter(posorpre,p,pn,pn->childllitem,threshold,depth+1);}
  if (childllitem->kidr!=NULL){ pnodeZZZ_starter(posorpre,p,parent,childllitem->kidr,threshold,depth+1);}
}

void pnodeZZZ_helper(int posorpre,struct ptree *p,struct pnode *pn,double threshold)
{
  /* checks pnode *pn for missing nodes given reference ptree *p */
  struct litem *l=NULL;
  struct llitem *l0=NULL;
  struct pnode *pn2=NULL;
  double expectance=0;
  struct llist *L=llistmake();
  if (posorpre){ llistgrowllitem(L,p->postree);} else{ llistgrowllitem(L,p->pretree);}
  l=L->first;
  while (l!=NULL){
    pn2 = (struct pnode *) l->item;
    if ((l0=llitemaddorfind(0,pn->childllitem,pn2->region,&region2pnode_compare_label))==NULL){
      expectance = pn->weight*p->legtime*pn2->weight/p->total_time;
      if (abs(expectance)>threshold){
	llitemaddorfind(1,pn->childllitem,pnodemake(pn,pn2->region,0,pnode_shear(0,expectance)),&pnode2pnode_compare_label);}}
    l=l->child;}
  llisttfree(L);L=NULL;
}

void pnodebalance_starter(struct pnode *parent,struct llitem *childllitem)
{
  int firstentry=0;
  struct pnode *pn=NULL;
  if (parent!=NULL && childllitem->parent==NULL){ /* first descent into parent->childllitem */ 
    assert(parent->childllitem==childllitem); firstentry=1;}
  if (childllitem->kidl!=NULL){ pnodebalance_starter(parent,childllitem->kidl);}
  if (childllitem->kidr!=NULL){ pnodebalance_starter(parent,childllitem->kidr);}
  if (childllitem->item!=NULL){ pn=(struct pnode *)childllitem->item; pnodebalance_starter(pn,pn->childllitem);}
  if (firstentry){ 
    if (parent!=NULL && llitemlength(parent->childllitem)>2*parent->broodsize){ 
      llitembalance(parent->childllitem); parent->childllitem = llitemclimb(parent->childllitem); 
      parent->broodsize = maximum(parent->broodsize*2,llitemlength(parent->childllitem));}}
}

void ptreerate(struct ptree *p)
{
  int verbose=0;
  struct llist *L=NULL;
  struct pnode *pn=NULL;
  if (verbose){ printf(" %% starting ptreerate, calling pnodeprune_starter\n");}
  L=llistmake(); pnodeprune_starter(1,p,NULL,p->postree,L); llisttfree(L); L=NULL;
  L=llistmake(); pnodeprune_starter(0,p,NULL,p->pretree,L); llisttfree(L); L=NULL;
  if (GLOBAL_PTREE_ZZZ){
    if (verbose){ printf(" %% finished calling pnodeprune_starter, calling pnodeZZZ_starter\n");}
    pnodeZZZ_starter(1,p,NULL,p->postree,1,0);
    pnodeZZZ_starter(0,p,NULL,p->pretree,1,0);
    if (verbose){ printf(" %% finished calling pnodeZZZ_starter, calling pnodebalance_starter\n");}}
  pnodebalance_starter(NULL,p->postree);
  pnodebalance_starter(NULL,p->pretree);
  L=llistmake(); pnode2llist_starter(NULL,p->postree,1,0,L); llistsort(L->first,L->last,L->length,pnode2pnode_compare_weight);
  pn = (struct pnode *) L->last->item;
  histtfree(p->wh); p->wh=histmake(128,pn->weight,0);
  llisttfree(L); L=NULL;
  histtfree(p->rh); p->rh=histmake(128,+4,-4);
  pnode2hist_starter(NULL,p->postree,-1,0,p->wh,p->rh);
}

struct llist * ptreextract_starter(struct ptree *p,double minworth,int nregions)
{
  /* extracts nregions most relevant regions from rated ptree *p, 
     simply by scanning through all legs and determining which 
     regions have the greatest worth = (relevance>threshold)*weight
     last_event is used to hold worth */
  int verbose=0;
  struct region **regionra=NULL;
  int nr=0,depth=0;
  struct llist *Lr=NULL,*Ll=NULL;
  struct litem *l=NULL;
  struct region *r=NULL;
  double *label=NULL;
  if (verbose){ printf(" %% [entering ptreextract_starter] with minworth %0.2f and nregions %d\n",minworth,nregions);}
  if (verbose){ printf(" %% making temporary regionra\n");}
  regionra = (struct region **) tcalloc(p->nregions,sizeof(struct region *));
  for (nr=0;nr<p->nregions;nr++){
    regionra[nr] = (struct region *) tcalloc(1,sizeof(struct region));
    regionra[nr]->label = p->regionra[nr]->label;
    regionra[nr]->neuronllist = NULL;}
  if (verbose){ printf(" %% calling ptreextract_helper\n");}
  ptreextract_helper(NULL,p->postree,regionra,minworth);
  if (verbose){ printf(" %% making and sorting regionllist\n");}
  Lr = llistmake(); for (nr=0;nr<p->nregions;nr++){ litemadd(Lr,regionra[nr]);} tfree(regionra); regionra=NULL;
  llistsort(Lr->first,Lr->last,Lr->length,&region2region_compare_last_event);
  if (verbose){ printf(" %% making labellist\n");}
  Ll = llistmake(); depth=0;
  l=Lr->last; 
  while (l!=NULL && depth<nregions){
    r = (struct region *) l->item;
    if (verbose){ printf(" %% %% adding label %d\n",r->label);}
    label = (double *) tcalloc(1,sizeof(double)); *label = (double)r->label; litemadd(Ll,label);
    l=l->parent; depth+=1;}
  if (verbose){ printf(" %% freeing temporary regions\n");}
  l=Lr->first;
  while (l!=NULL){
    r = (struct region *) l->item;
    tfree(r); r=NULL; l->item=NULL;
    l=l->child;}
  llisttfree(Lr);
  llistsort(Ll->first,Ll->last,Ll->length,&double_compare);
  return Ll;
}

void ptreextract_helper(struct pnode *parent,struct llitem *l0,struct region **regionra,double minworth)
{
  int firstentry=0;
  struct pnode *pn=NULL;
  if (parent!=NULL && l0->parent==NULL){ /* first descent into parent->childllitem */ 
    assert(parent->childllitem==l0); firstentry=1;}
  if (firstentry){ regionra[parent->region->label]->last_event += pnode_worth(parent->weight,parent->relevance,minworth);}
  if (l0->kidl!=NULL){ ptreextract_helper(parent,l0->kidl,regionra,minworth);}
  if (l0->item!=NULL){ pn=(struct pnode *)l0->item; ptreextract_helper(pn,pn->childllitem,regionra,minworth);}
  if (l0->kidr!=NULL){ ptreextract_helper(parent,l0->kidr,regionra,minworth);}
}

double pnode_worth(double weight,double relevance,double minworth){ return weight*(fabs(relevance)+minworth);}

void ptree2jpg_starter(struct ptree *p,char *filename_given,double wmax,double wmin,double rmax,double rmin)
{
  /* assuming filename_given starts with "./" */
  int remove_flag=0,jpg_flag=0,eps_flag=(ON_MY_COMPUTER ? 1 : 0);
  int coloringtype=7;
  char filename_w[256],filename2_w[512],filename_r[256],filename2_r[512],filename_wh[512],filename_rh[512];
  char filename_ws[256],filename2_ws[512],filename_rs[256],filename2_rs[512];
  FILE *fpw=NULL,*fpr=NULL,*fpws=NULL,*fprs=NULL;
  struct llist *L=NULL;
  char command[1024];
  double wmax2=0,wmin2=0,wmean=0,wstd=0,rmax2=0,rmin2=0,rmean=0,rstd=0;
  pnodestats_starter(NULL,p->postree,3,-1,0,NULL,&wmax2,&wmin2,&wmean,&wstd,&rmax2,&rmin2,&rmean,&rstd);
  if (wmax<wmin){ wmax = wmean + STD_VIEW*wstd; wmin = wmean - STD_VIEW*wstd;}
  else if (wmax==wmin){ wmax = wmax2; wmin = wmin2;}
  if (rmax<rmin){ rmax = rmean + STD_VIEW*rstd; rmin = rmean - STD_VIEW*rstd;}
  else if (rmax==rmin){ rmax = rmax2; rmin = rmin2;}
  if (filename_given==NULL){ sprintf(filename_w,"./ptree_%srecord_w",GLOBAL_STRING_2);}
  else{ sprintf(filename_w,"%s_w",filename_given);}
  sprintf(filename2_w,"%s.fig",filename_w);
  if ((fpw=fopen(filename2_w,"w"))==NULL){ printf("warning, cannot create %s in ptree2jpg_starter",filename2_w); fpw=stdout;}
  if (filename_given==NULL){ sprintf(filename_r,"./ptree_%srecord_r",GLOBAL_STRING_2);}
  else{ sprintf(filename_r,"%s_r",filename_given);}
  sprintf(filename2_r,"%s.fig",filename_r);
  if ((fpr=fopen(filename2_r,"w"))==NULL){ printf("warning, cannot create %s in ptree2jpg_starter",filename2_r); fpr=stdout;}
  if (filename_given==NULL){ sprintf(filename_ws,"./ptree_%srecord_ws",GLOBAL_STRING_2);}
  else{ sprintf(filename_ws,"%s_ws",filename_given);}
  sprintf(filename2_ws,"%s.fig",filename_ws);
  if ((fpws=fopen(filename2_ws,"w"))==NULL){ printf("warning, cannot create %s in ptree2jpg_starter",filename2_ws); fpws=stdout;}
  if (filename_given==NULL){ sprintf(filename_rs,"./ptree_%srecord_rs",GLOBAL_STRING_2);}
  else{ sprintf(filename_rs,"%s_rs",filename_given);}
  sprintf(filename2_rs,"%s.fig",filename_rs);
  if ((fprs=fopen(filename2_rs,"w"))==NULL){ printf("warning, cannot create %s in ptree2jpg_starter",filename2_rs); fprs=stdout;}
  if (filename_given==NULL){ sprintf(filename_wh,"./ptree_%srecord_wh",GLOBAL_STRING_2);}
  else{ sprintf(filename_wh,"%s_wh",filename_given);}
  if (filename_given==NULL){ sprintf(filename_rh,"./ptree_%srecord_rh",GLOBAL_STRING_2);}
  else{ sprintf(filename_rh,"%s_rh",filename_given);}
  histdump(p->wh,0,filename_wh," ",0); histdump(p->wh,1,filename_wh,"wh",0);
  histdump(p->rh,0,filename_rh," ",0); histdump(p->rh,1,filename_rh,"rh",0);
  fprintf(fpw,"%s",FIG_PREAMBLE); fprintf(fpr,"%s",FIG_PREAMBLE);
  fprintf(fpws,"%s",FIG_PREAMBLE); fprintf(fprs,"%s",FIG_PREAMBLE);
  switch (coloringtype){
  case 0: fprintf(fpw,"%s",FIG_PREAMBLE_COLOR_0); fprintf(fpr,"%s",FIG_PREAMBLE_COLOR_0); break;
  case 5: fprintf(fpw,"%s",FIG_PREAMBLE_COLOR_5); fprintf(fpr,"%s",FIG_PREAMBLE_COLOR_5); break;
  case 7: fprintf(fpw,"%s",FIG_PREAMBLE_COLOR_7); fprintf(fpr,"%s",FIG_PREAMBLE_COLOR_7); break;
  default: break;}
  fprintf(fpws,"%s",FIG_PREAMBLE_COLOR_7); fprintf(fprs,"%s",FIG_PREAMBLE_COLOR_7);
  ptree2jpg_helper(+1,p->nregions,wmax,wmin,rmax,rmin,NULL,NULL,fpw,fpr,NULL);
  ptree2jpg_helper(-1,p->nregions,wmax,wmin,rmax,rmin,NULL,NULL,fpw,fpr,NULL);
  L=llistmake();
  ptree2jpg_helper(+1,p->nregions,wmax,wmin,rmax,rmin,NULL,p->postree,fpw,fpr,L);
  llisttfree(L); L=NULL;
  L=llistmake();
  ptree2jpg_helper(-1,p->nregions,wmax,wmin,rmax,rmin,NULL,p->pretree,fpw,fpr,L);
  llisttfree(L); L=NULL;
  L=llistmake();
  ptree2star_helper(+1,p->nlegs,p->nregions,-1,wmax2,1,rmax,rmin,NULL,p->postree,fpws,fprs,L);
  ptree2star_helper(-1,p->nlegs,p->nregions,-1,wmax2,1,rmax,rmin,NULL,p->pretree,fpws,fprs,L);
  llisttfree(L); L=NULL;
  if (fpw!=stdout){ fclose(fpw);} if (fpr!=stdout){ fclose(fpr);}
  if (fpws!=stdout){ fclose(fpws);} if (fprs!=stdout){ fclose(fprs);}
  if (jpg_flag){ sprintf(command,"fig2dev -Ljpeg -q %d %s.fig %s.jpg;",/*quality*/5,filename_w,filename_w); system(command);}
  if (jpg_flag){ sprintf(command,"fig2dev -Ljpeg -q %d %s.fig %s.jpg;",/*quality*/5,filename_r,filename_r); system(command);}
  if (jpg_flag){ sprintf(command,"fig2dev -Ljpeg -q %d %s.fig %s.jpg;",/*quality*/5,filename_ws,filename_ws); system(command);}
  if (jpg_flag){ sprintf(command,"fig2dev -Ljpeg -q %d %s.fig %s.jpg;",/*quality*/5,filename_rs,filename_rs); system(command);}
  if (eps_flag){ sprintf(command,"fig2dev -Leps %s.fig %s.eps;",filename_w,filename_w); system(command);}
  if (eps_flag){ sprintf(command,"fig2dev -Leps %s.fig %s.eps;",filename_r,filename_r); system(command);}
  if (eps_flag){ sprintf(command,"fig2dev -Leps %s.fig %s.eps;",filename_ws,filename_ws); system(command);}
  if (eps_flag){ sprintf(command,"fig2dev -Leps %s.fig %s.eps;",filename_rs,filename_rs); system(command);}
  if (remove_flag){ sprintf(command,"rm %s.fig;",filename_w); system(command);}
  if (remove_flag){ sprintf(command,"rm %s.fig;",filename_r); system(command);}
  if (remove_flag){ sprintf(command,"rm %s.fig;",filename_ws); system(command);}
  if (remove_flag){ sprintf(command,"rm %s.fig;",filename_rs); system(command);}
}

void ptree2jpg_helper(int posorpre,int nregions,double wmax,double wmin,double rmax,double rmin,struct pnode *parent,struct llitem *l0,FILE *fpw,FILE *fpr,struct llist *L)
{
  /* runs through llitem *l0 of pnode items, the relevance is dumped up to depth 3
     any pnode order works, but we do reachfirst.
  */
  int firstentry=0;
  int nr=0,maxdia=10000;
  struct pnode *pn=NULL;
  double rcolor=0,gcolor=0,bcolor=0;
  int colorcode=0;
  struct region *r1=NULL,*r2=NULL,*r3=NULL,*r4=NULL;
  double rord = 1.0/(double)nregions,side=0,xside=0,yside=0;
  double xord=0,yord=0;
  if (l0==NULL){ /* draw skeleton */
    xord = 0.5 + 0.1*posorpre;
    yord = 1.25*(double)nregions/2.0 + 1*(double)(nregions-1.0)/2.0*rord;
    xside = rord; yside = 1.0;
    if (posorpre==+1){
      fprintf(fpw,"4 0 0 %d 0 %d %d 0.0000 4 %d %d %d %d %0.5f\\001\n",/*depth*/1,/*font*/12,/*point*/72,/*textheight*/(int)(24*13.5),/*textwidth*/(int)(24*9*18),/*xpos*/(int)floor(maxdia*(xord-2.5*rord)),/*ypos*/(int)floor(maxdia*(yord-(nregions+1)*rord/2)),wmax);
      fprintf(fpw,"4 0 0 %d 0 %d %d 0.0000 4 %d %d %d %d %0.5f\\001\n",/*depth*/1,/*font*/12,/*point*/72,/*textheight*/(int)(24*13.5),/*textwidth*/(int)(24*9*18),/*xpos*/(int)floor(maxdia*(xord-2.5*rord)),/*ypos*/(int)floor(maxdia*(yord+(nregions+1.5)*rord/2)),wmin);
      fprintf(fpr,"4 0 0 %d 0 %d %d 0.0000 4 %d %d %d %d %0.5f\\001\n",/*depth*/1,/*font*/12,/*point*/72,/*textheight*/(int)(24*13.5),/*textwidth*/(int)(24*9*18),/*xpos*/(int)floor(maxdia*(xord-2.5*rord)),/*ypos*/(int)floor(maxdia*(yord-(nregions+1)*rord/2)),rmax);
      fprintf(fpr,"4 0 0 %d 0 %d %d 0.0000 4 %d %d %d %d %0.5f\\001\n",/*depth*/1,/*font*/12,/*point*/72,/*textheight*/(int)(24*13.5),/*textwidth*/(int)(24*9*18),/*xpos*/(int)floor(maxdia*(xord-2.5*rord)),/*ypos*/(int)floor(maxdia*(yord+(nregions+1.5)*rord/2)),rmin);}
    fprintf(fpw,"2 2 0 4 0 7 900 0 -1 0.000 0 0 -1 0 0 5\n");
    fprintf(fpw,"\t %d %d %d %d %d %d %d %d %d %d\n",(int)floor(maxdia*(xord+xside/2)),(int)floor(maxdia*(yord+yside/2)),(int)floor(maxdia*(xord+xside/2)),(int)floor(maxdia*(yord-yside/2)),(int)floor(maxdia*(xord-xside/2)),(int)floor(maxdia*(yord-yside/2)),(int)floor(maxdia*(xord-xside/2)),(int)floor(maxdia*(yord+yside/2)),(int)floor(maxdia*(xord+xside/2)),(int)floor(maxdia*(yord+yside/2)));
    fprintf(fpr,"2 2 0 4 0 7 900 0 -1 0.000 0 0 -1 0 0 5\n");
    fprintf(fpr,"\t %d %d %d %d %d %d %d %d %d %d\n",(int)floor(maxdia*(xord+xside/2)),(int)floor(maxdia*(yord+yside/2)),(int)floor(maxdia*(xord+xside/2)),(int)floor(maxdia*(yord-yside/2)),(int)floor(maxdia*(xord-xside/2)),(int)floor(maxdia*(yord-yside/2)),(int)floor(maxdia*(xord-xside/2)),(int)floor(maxdia*(yord+yside/2)),(int)floor(maxdia*(xord+xside/2)),(int)floor(maxdia*(yord+yside/2)));
    for (nr=0;nr<nregions;nr++){
      xord = 0.5 + 0.75*posorpre;
      yord = 1.25*(double)nr + 1*(double)(nregions-1.0)/2.0*rord;
      xside = rord; yside = 1.0;
      fprintf(fpw,"2 2 0 4 0 7 900 0 -1 0.000 0 0 -1 0 0 5\n");
      fprintf(fpw,"\t %d %d %d %d %d %d %d %d %d %d\n",(int)floor(maxdia*(xord+xside/2)),(int)floor(maxdia*(yord+yside/2)),(int)floor(maxdia*(xord+xside/2)),(int)floor(maxdia*(yord-yside/2)),(int)floor(maxdia*(xord-xside/2)),(int)floor(maxdia*(yord-yside/2)),(int)floor(maxdia*(xord-xside/2)),(int)floor(maxdia*(yord+yside/2)),(int)floor(maxdia*(xord+xside/2)),(int)floor(maxdia*(yord+yside/2)));
      fprintf(fpr,"2 2 0 4 0 7 900 0 -1 0.000 0 0 -1 0 0 5\n");
      fprintf(fpr,"\t %d %d %d %d %d %d %d %d %d %d\n",(int)floor(maxdia*(xord+xside/2)),(int)floor(maxdia*(yord+yside/2)),(int)floor(maxdia*(xord+xside/2)),(int)floor(maxdia*(yord-yside/2)),(int)floor(maxdia*(xord-xside/2)),(int)floor(maxdia*(yord-yside/2)),(int)floor(maxdia*(xord-xside/2)),(int)floor(maxdia*(yord+yside/2)),(int)floor(maxdia*(xord+xside/2)),(int)floor(maxdia*(yord+yside/2)));}
    for (nr=0;nr<nregions;nr++){
      xord = 1.75*posorpre + (double)(nregions-1.0)/2.0*rord;
      yord = 1.25*(double)nr + 1*(double)(nregions-1.0)/2.0*rord;
      xside = 1.0; yside = 1.0;
      fprintf(fpw,"2 2 0 4 0 7 900 0 -1 0.000 0 0 -1 0 0 5\n");
      fprintf(fpw,"\t %d %d %d %d %d %d %d %d %d %d\n",(int)floor(maxdia*(xord+xside/2)),(int)floor(maxdia*(yord+yside/2)),(int)floor(maxdia*(xord+xside/2)),(int)floor(maxdia*(yord-yside/2)),(int)floor(maxdia*(xord-xside/2)),(int)floor(maxdia*(yord-yside/2)),(int)floor(maxdia*(xord-xside/2)),(int)floor(maxdia*(yord+yside/2)),(int)floor(maxdia*(xord+xside/2)),(int)floor(maxdia*(yord+yside/2)));
      fprintf(fpr,"2 2 0 4 0 7 900 0 -1 0.000 0 0 -1 0 0 5\n");
      fprintf(fpr,"\t %d %d %d %d %d %d %d %d %d %d\n",(int)floor(maxdia*(xord+xside/2)),(int)floor(maxdia*(yord+yside/2)),(int)floor(maxdia*(xord+xside/2)),(int)floor(maxdia*(yord-yside/2)),(int)floor(maxdia*(xord-xside/2)),(int)floor(maxdia*(yord-yside/2)),(int)floor(maxdia*(xord-xside/2)),(int)floor(maxdia*(yord+yside/2)),(int)floor(maxdia*(xord+xside/2)),(int)floor(maxdia*(yord+yside/2)));}
    for (nr=0;nr<nregions;nr++){
      xord = 3.0*posorpre + (double)(nregions-1.0)/2.0*rord;
      yord = 1.25*(double)nr + 1*(double)(nregions-1.0)/2.0*rord;
      xside = 1.0; yside = 1.0;
      fprintf(fpw,"2 2 0 4 0 7 900 0 -1 0.000 0 0 -1 0 0 5\n");
      fprintf(fpw,"\t %d %d %d %d %d %d %d %d %d %d\n",(int)floor(maxdia*(xord+xside/2)),(int)floor(maxdia*(yord+yside/2)),(int)floor(maxdia*(xord+xside/2)),(int)floor(maxdia*(yord-yside/2)),(int)floor(maxdia*(xord-xside/2)),(int)floor(maxdia*(yord-yside/2)),(int)floor(maxdia*(xord-xside/2)),(int)floor(maxdia*(yord+yside/2)),(int)floor(maxdia*(xord+xside/2)),(int)floor(maxdia*(yord+yside/2)));
      fprintf(fpr,"2 2 0 4 0 7 900 0 -1 0.000 0 0 -1 0 0 5\n");
      fprintf(fpr,"\t %d %d %d %d %d %d %d %d %d %d\n",(int)floor(maxdia*(xord+xside/2)),(int)floor(maxdia*(yord+yside/2)),(int)floor(maxdia*(xord+xside/2)),(int)floor(maxdia*(yord-yside/2)),(int)floor(maxdia*(xord-xside/2)),(int)floor(maxdia*(yord-yside/2)),(int)floor(maxdia*(xord-xside/2)),(int)floor(maxdia*(yord+yside/2)),(int)floor(maxdia*(xord+xside/2)),(int)floor(maxdia*(yord+yside/2)));}}
  else /* if (l0!=NULL) */{ /* fill out ptree */
    if (parent!=NULL && l0->parent==NULL){ /* first descent into parent->childllitem */ 
      assert(parent->childllitem==l0); firstentry=1;}
    if (firstentry){ litemadd(L,parent->region);}
    if (firstentry){ 
      switch(L->length){
      case 0: break;
      case 1: 
	r1 = (struct region *) L->last->item;
	colorscale(0,parent->relevance,rmax,rmin,&rcolor,&gcolor,&bcolor);
	colorcode = crop((int)floor(512*rcolor),0,511);
	xord = 0.5 + 0.1*posorpre;
	yord = 1.25*(double)nregions/2.0 + 1*(double)r1->label*rord;
	side = rord;
	fprintf(fpr,"2 2 0 0 -1 %d %d 0 %d 0.000 0 0 -1 0 0 5\n",/*colorcode*/colorcode+32,/*depth*/1+r1->label,/*fill*/20);
	fprintf(fpr,"\t %d %d %d %d %d %d %d %d %d %d\n",(int)floor(maxdia*(xord+side/2)),(int)floor(maxdia*(yord+side/2)),(int)floor(maxdia*(xord+side/2)),(int)floor(maxdia*(yord-side/2)),(int)floor(maxdia*(xord-side/2)),(int)floor(maxdia*(yord-side/2)),(int)floor(maxdia*(xord-side/2)),(int)floor(maxdia*(yord+side/2)),(int)floor(maxdia*(xord+side/2)),(int)floor(maxdia*(yord+side/2)));
	colorscale(0,parent->weight,wmax,wmin,&rcolor,&gcolor,&bcolor);
	colorcode = crop((int)floor(512*rcolor),0,511);
	fprintf(fpw,"2 2 0 0 -1 %d %d 0 %d 0.000 0 0 -1 0 0 5\n",/*colorcode*/colorcode+32,/*depth*/1+r1->label,/*fill*/20);
	fprintf(fpw,"\t %d %d %d %d %d %d %d %d %d %d\n",(int)floor(maxdia*(xord+side/2)),(int)floor(maxdia*(yord+side/2)),(int)floor(maxdia*(xord+side/2)),(int)floor(maxdia*(yord-side/2)),(int)floor(maxdia*(xord-side/2)),(int)floor(maxdia*(yord-side/2)),(int)floor(maxdia*(xord-side/2)),(int)floor(maxdia*(yord+side/2)),(int)floor(maxdia*(xord+side/2)),(int)floor(maxdia*(yord+side/2)));
	break;
      case 2:
	r1 = (struct region *) L->last->parent->item;
	r2 = (struct region *) L->last->item;
	colorscale(0,parent->relevance,rmax,rmin,&rcolor,&gcolor,&bcolor);
	colorcode = crop((int)floor(512*rcolor),0,511);
	xord = 0.5 + 0.75*posorpre;
	yord = 1.25*(double)r1->label + 1*(double)r2->label*rord;
	side = rord;
	fprintf(fpr,"2 2 0 0 -1 %d %d 0 %d 0.000 0 0 -1 0 0 5\n",/*colorcode*/colorcode+32,/*depth*/1+r2->label,/*fill*/20);
	fprintf(fpr,"\t %d %d %d %d %d %d %d %d %d %d\n",(int)floor(maxdia*(xord+side/2)),(int)floor(maxdia*(yord+side/2)),(int)floor(maxdia*(xord+side/2)),(int)floor(maxdia*(yord-side/2)),(int)floor(maxdia*(xord-side/2)),(int)floor(maxdia*(yord-side/2)),(int)floor(maxdia*(xord-side/2)),(int)floor(maxdia*(yord+side/2)),(int)floor(maxdia*(xord+side/2)),(int)floor(maxdia*(yord+side/2)));
	colorscale(0,parent->weight,wmax,wmin,&rcolor,&gcolor,&bcolor);
	colorcode = crop((int)floor(512*rcolor),0,511);
	fprintf(fpw,"2 2 0 0 -1 %d %d 0 %d 0.000 0 0 -1 0 0 5\n",/*colorcode*/colorcode+32,/*depth*/1+r2->label,/*fill*/20);
	fprintf(fpw,"\t %d %d %d %d %d %d %d %d %d %d\n",(int)floor(maxdia*(xord+side/2)),(int)floor(maxdia*(yord+side/2)),(int)floor(maxdia*(xord+side/2)),(int)floor(maxdia*(yord-side/2)),(int)floor(maxdia*(xord-side/2)),(int)floor(maxdia*(yord-side/2)),(int)floor(maxdia*(xord-side/2)),(int)floor(maxdia*(yord+side/2)),(int)floor(maxdia*(xord+side/2)),(int)floor(maxdia*(yord+side/2)));
	break;
      case 3:
	r1 = (struct region *) L->last->parent->parent->item;
	r2 = (struct region *) L->last->parent->item;
	r3 = (struct region *) L->last->item;
	colorscale(0,parent->relevance,rmax,rmin,&rcolor,&gcolor,&bcolor);
	colorcode = crop((int)floor(512*rcolor),0,511);
	xord = 1.75*posorpre + (double)r2->label*rord;
	yord = 1.25*(double)r1->label + 1*(double)r3->label*rord;
	side = rord;
	fprintf(fpr,"2 2 0 0 -1 %d %d 0 %d 0.000 0 0 -1 0 0 5\n",/*colorcode*/colorcode+32,/*depth*/1+r3->label,/*fill*/20);
	fprintf(fpr,"\t %d %d %d %d %d %d %d %d %d %d\n",(int)floor(maxdia*(xord+side/2)),(int)floor(maxdia*(yord+side/2)),(int)floor(maxdia*(xord+side/2)),(int)floor(maxdia*(yord-side/2)),(int)floor(maxdia*(xord-side/2)),(int)floor(maxdia*(yord-side/2)),(int)floor(maxdia*(xord-side/2)),(int)floor(maxdia*(yord+side/2)),(int)floor(maxdia*(xord+side/2)),(int)floor(maxdia*(yord+side/2)));
	colorscale(0,parent->weight,wmax,wmin,&rcolor,&gcolor,&bcolor);
	colorcode = crop((int)floor(512*rcolor),0,511);
	fprintf(fpw,"2 2 0 0 -1 %d %d 0 %d 0.000 0 0 -1 0 0 5\n",/*colorcode*/colorcode+32,/*depth*/1+r3->label,/*fill*/20);
	fprintf(fpw,"\t %d %d %d %d %d %d %d %d %d %d\n",(int)floor(maxdia*(xord+side/2)),(int)floor(maxdia*(yord+side/2)),(int)floor(maxdia*(xord+side/2)),(int)floor(maxdia*(yord-side/2)),(int)floor(maxdia*(xord-side/2)),(int)floor(maxdia*(yord-side/2)),(int)floor(maxdia*(xord-side/2)),(int)floor(maxdia*(yord+side/2)),(int)floor(maxdia*(xord+side/2)),(int)floor(maxdia*(yord+side/2)));
	break;
      case 4:
	r1 = (struct region *) L->last->parent->parent->parent->item;
	r2 = (struct region *) L->last->parent->parent->item;
	r3 = (struct region *) L->last->parent->item;
	r4 = (struct region *) L->last->item;
	colorscale(0,parent->relevance,rmax,rmin,&rcolor,&gcolor,&bcolor);
	colorcode = crop((int)floor(512*rcolor),0,511);
	xord = 3.0*posorpre + (double)r2->label*rord;
	yord = 1.25*(double)r1->label + 1*(double)r3->label*rord;
	side = (double)(1 + r4->label)*rord*rord;
	fprintf(fpr,"2 2 0 0 -1 %d %d 0 %d 0.000 0 0 -1 0 0 5\n",/*colorcode*/colorcode+32,/*depth*/1+r4->label,/*fill*/20);
	fprintf(fpr,"\t %d %d %d %d %d %d %d %d %d %d\n",(int)floor(maxdia*(xord+side/2)),(int)floor(maxdia*(yord+side/2)),(int)floor(maxdia*(xord+side/2)),(int)floor(maxdia*(yord-side/2)),(int)floor(maxdia*(xord-side/2)),(int)floor(maxdia*(yord-side/2)),(int)floor(maxdia*(xord-side/2)),(int)floor(maxdia*(yord+side/2)),(int)floor(maxdia*(xord+side/2)),(int)floor(maxdia*(yord+side/2)));
	colorscale(0,parent->weight,wmax,wmin,&rcolor,&gcolor,&bcolor);
	colorcode = crop((int)floor(512*rcolor),0,511);
	fprintf(fpw,"2 2 0 0 -1 %d %d 0 %d 0.000 0 0 -1 0 0 5\n",/*colorcode*/colorcode+32,/*depth*/1+r4->label,/*fill*/20);
	fprintf(fpw,"\t %d %d %d %d %d %d %d %d %d %d\n",(int)floor(maxdia*(xord+side/2)),(int)floor(maxdia*(yord+side/2)),(int)floor(maxdia*(xord+side/2)),(int)floor(maxdia*(yord-side/2)),(int)floor(maxdia*(xord-side/2)),(int)floor(maxdia*(yord-side/2)),(int)floor(maxdia*(xord-side/2)),(int)floor(maxdia*(yord+side/2)),(int)floor(maxdia*(xord+side/2)),(int)floor(maxdia*(yord+side/2)));
	break;      
      default: /* do nothing */ break;}}
    if (l0->kidl!=NULL){ ptree2jpg_helper(posorpre,nregions,wmax,wmin,rmax,rmin,parent,l0->kidl,fpw,fpr,L);}
    if (l0->item!=NULL){ pn=(struct pnode *)l0->item; ptree2jpg_helper(posorpre,nregions,wmax,wmin,rmax,rmin,pn,pn->childllitem,fpw,fpr,L);}
    if (l0->kidr!=NULL){ ptree2jpg_helper(posorpre,nregions,wmax,wmin,rmax,rmin,parent,l0->kidr,fpw,fpr,L);}
    if (firstentry){ llistkilllast(L);}}
}

void ptree2star_helper(int posorpre,int nlegs,int nregions,int maxlevel,double wmax,double wmin,double rmax,double rmin,struct pnode *parent,struct llitem *l0,FILE *fpws,FILE *fprs,struct llist *L)
{
  /* runs through llitem *l0 of pnode items, the relevance is dumped up to depth 3
     any pnode order works, but we do reachfirst.
  */
  int firstentry=0;
  int nr=0,maxdia=10000;
  struct pnode *pn=NULL;
  struct litem *l=NULL;
  double rcolor=0,gcolor=0,bcolor=0;
  int colorcode=0;
  struct region *rg=NULL;
  int nl=0;
  double xord=0,yord=0,xord2=0,yord2=0,rord=0;
  double vx=0,vy=0,vrad=0,epsilon=0;
  if (parent!=NULL && l0->parent==NULL){ /* first descent into parent->childllitem */ 
    assert(parent->childllitem==l0); firstentry=1;}
  if (firstentry){ litemadd(L,parent->region);}
  if (firstentry){ 
    l=L->last;nl=0;xord=0,yord=0,vx=0;vy=0;vrad=PI/(double)nregions;vrad=0.5*sin(vrad)/(1+sin(vrad));epsilon=1.75;
    while (l!=NULL){
      rg=(struct region *)l->item;
      vx = cos(2*PI*(rg->label+0.5)/(double)nregions)*(0.5-vrad);
      vy = sin(2*PI*(rg->label+0.5)/(double)nregions)*(0.5-vrad);
      xord += pow(3+epsilon,nl)*vx; yord += pow(3+epsilon,nl)*vy;
      nl+=1;
      l=l->parent;}
    xord /= pow(3+epsilon,maxlevel>0?maxlevel:nlegs);
    yord /= pow(3+epsilon,maxlevel>0?maxlevel:nlegs);
    rord = vrad*pow(3+epsilon,-(maxlevel>0?maxlevel:nlegs));
    xord2 = xord + 0.5*posorpre; yord2 = yord;
    colorscale(0,parent->relevance,rmax,rmin,&rcolor,&gcolor,&bcolor);colorcode=crop((int)floor(512*rcolor),0,511);
    fprintf(fprs,"2 1 0 0 0 %d %d 0 %d 0.000 0 0 -1 0 0 %d\n",/*colorcode*/colorcode+32,/*depth*/1+L->length,/*fill*/20,/*npoints*/nregions+1); fprintf(fprs,"\t"); for (nr=0;nr<=nregions;nr++){ fprintf(fprs,"%d %d ",(int)floor(maxdia*(xord2+rord*cos(2*PI*((double)nr+0.5)/(double)nregions))),(int)floor(maxdia*(yord2+rord*sin(2*PI*((double)nr+0.5)/(double)nregions))));} fprintf(fprs,"\n");
    colorscale(0,log(maximum(1,parent->weight)),log(wmax),log(wmin),&rcolor,&gcolor,&bcolor);colorcode=crop((int)floor(512*rcolor),0,511);
    fprintf(fpws,"2 1 0 0 0 %d %d 0 %d 0.000 0 0 -1 0 0 %d\n",/*colorcode*/colorcode+32,/*depth*/1+L->length,/*fill*/20,/*npoints*/nregions+1); fprintf(fpws,"\t"); for (nr=0;nr<=nregions;nr++){ fprintf(fpws,"%d %d ",(int)floor(maxdia*(xord2+rord*cos(2*PI*((double)nr+0.5)/(double)nregions))),(int)floor(maxdia*(yord2+rord*sin(2*PI*((double)nr+0.5)/(double)nregions))));} fprintf(fpws,"\n");}
  if (l0->kidl!=NULL){ ptree2star_helper(posorpre,nlegs,nregions,maxlevel,wmax,wmin,rmax,rmin,parent,l0->kidl,fpws,fprs,L);}
  if (l0->item!=NULL){ 
    if (maxlevel==-1 || L->length<maxlevel){
      pn=(struct pnode *)l0->item; 
      ptree2star_helper(posorpre,nlegs,nregions,maxlevel,wmax,wmin,rmax,rmin,pn,pn->childllitem,fpws,fprs,L);}}
  if (l0->kidr!=NULL){ ptree2star_helper(posorpre,nlegs,nregions,maxlevel,wmax,wmin,rmax,rmin,parent,l0->kidr,fpws,fprs,L);}
  if (firstentry){ llistkilllast(L);}
}

void ptreedump_starter(struct ptree *p,char *fgvn,int dump_type,int fullname_flag,double wmax,double wmin,double rmax,double rmin)
{
  /* dumps ptree *p as a sparse array, first dumping
     p->nregions
     p->nlegs
     p->legtime
     p->length 
     p->update_every
     p->total_time
     assumes fgvn does NOT start with "./"
     unless fullname_flag==1
  */
  char filename[256];
  FILE *fp=NULL;
  struct llist *L=NULL;
  int bitbybit=0,continue_flag=0;
  if (dump_type==2){
    if (fullname_flag==1){
      sprintf(filename,"%s",fgvn);}
    else{
      if (GLOBAL_PTREE_BITBYBIT==0){
	if (fgvn==NULL){ sprintf(filename,"./ptree_%srecord_time%d",GLOBAL_STRING_2,(int)floor(GLOBAL_time));}
	else /* if (fgvn!=NULL) */{ sprintf(filename,"./ptree_%s_%srecord_time%d",fgvn,GLOBAL_STRING_2,(int)floor(GLOBAL_time));}}
      else /* if (GLOBAL_PTREE_BITBYBIT==1) */{
	bitbybit=0;continue_flag=0;
	do{
	  bitbybit+=1;
	  if (fgvn==NULL){ sprintf(filename,"./ptree_%srecord_%d",GLOBAL_STRING_2,bitbybit);}
	  else /* if (fgvn!=NULL) */{ sprintf(filename,"./ptree_%s_%srecord_%d",fgvn,GLOBAL_STRING_2,bitbybit);}
	  if ((fp=fopen(filename,"r"))==NULL){ continue_flag=0;} else{ fclose(fp); continue_flag=1;}}
	while (continue_flag);}}
    if ((fp=fopen(filename,"w"))==NULL){ printf("warning, cannot create %s in ptreedump_starter",filename); fp=stdout;}
    fwrite(&(p->nregions),sizeof(int),1,fp);
    fwrite(&(p->nlegs),sizeof(int),1,fp);
    fwrite(&(p->legtime),sizeof(int),1,fp);
    fwrite(&(p->length),sizeof(int),1,fp);
    fwrite(&(p->update_every),sizeof(double),1,fp);
    fwrite(&(p->total_time),sizeof(double),1,fp);
    L=llistmake();
    ptreedump_helper(1,NULL,p->postree,fp,L);
    llisttfree(L); L=NULL;
    L=llistmake();
    ptreedump_helper(0,NULL,p->pretree,fp,L);
    llisttfree(L); L=NULL;
    if (fp!=stdout){ fclose(fp);}
    if (fullname_flag!=0 || GLOBAL_PTREE_BITBYBIT==0){ ptree2jpg_starter(p,filename,wmax,wmin,rmax,rmin);}}
}

void ptreedump_helper(int posorpre,struct pnode *parent,struct llitem *l0,FILE *fp,struct llist *L)
{
  /* runs through llitem *l0 of pnode items, the format is:
     parent->weight,parent->relevance,label_0,label_1,...,label_n,posorpre ? -1 : -2
     any pnode order works, but we do reachfirst.
  */
  int firstentry=0;
  int minusone=-1,minustwo=-2;
  struct pnode *pn=NULL;
  struct litem *l=NULL;
  struct region *r=NULL;
  if (parent!=NULL && l0->parent==NULL){ /* first descent into parent->childllitem */ 
    assert(parent->childllitem==l0); firstentry=1;}
  if (firstentry){ litemadd(L,parent->region);}
  if (firstentry){ 
    fwrite(&(parent->weight),sizeof(double),1,fp);
    fwrite(&(parent->relevance),sizeof(double),1,fp);
    l = L->first;
    while (l!=NULL){
      r = (struct region *) l->item;
      fwrite(&(r->label),sizeof(int),1,fp);
      l=l->child;}
    if (posorpre){ fwrite(&minusone,sizeof(int),1,fp);}
    else{ fwrite(&minustwo,sizeof(int),1,fp);}}
  if (l0->kidl!=NULL){ ptreedump_helper(posorpre,parent,l0->kidl,fp,L);}
  if (l0->item!=NULL){ pn=(struct pnode *)l0->item; ptreedump_helper(posorpre,pn,pn->childllitem,fp,L);}
  if (l0->kidr!=NULL){ ptreedump_helper(posorpre,parent,l0->kidr,fp,L);}
  if (firstentry){ llistkilllast(L);}
}

struct ptree * ptreadback(char *filename)
{
  int verbose=0;
  char filename2[256];
  FILE *fp=NULL;
  struct ptree *p=NULL;
  int nr=0;
  struct pnode *pn=NULL;
  struct llitem *l0=NULL,*l1=NULL;
  int readout=0;
  double weight=0,relevance=0;
  struct llist *L=NULL;
  struct litem *l=NULL;
  int *label=NULL;
  int event_threshold = GLOBAL_PTREE_EVENT_THRESHOLD;
  double event_within = GLOBAL_PTREE_EVENT_WITHIN;
  int region_type = GLOBAL_PTREE_REGION_TYPE;
  int exit_flag=0,depth=0;
  if (verbose){ printf(" %% \n");}
  if (verbose){ printf(" %% [entering ptreadback] with filename %s\n",filename==NULL ? "<null>" : filename);}
  if (filename==NULL){ sprintf(filename2,"./ptree_%srecord",GLOBAL_STRING_2);} else{ sprintf(filename2,"%s",filename);}
  if (verbose){ printf(" %% decided on filename %s\n",filename2);}
  if ((fp=fopen(filename2,"r"))==NULL){ printf("warning, cannot read %s in ptreadback\n",filename2);}
  else{
    p = (struct ptree *) tmalloc(sizeof(struct ptree));
    fread(&(p->nregions),sizeof(int),1,fp);
    fread(&(p->nlegs),sizeof(int),1,fp);
    fread(&(p->legtime),sizeof(int),1,fp);
    fread(&(p->length),sizeof(int),1,fp);
    fread(&(p->update_every),sizeof(double),1,fp);
    fread(&(p->total_time),sizeof(double),1,fp);
    if (verbose){ printf(" %% nregions %d, nlegs %d, legtime %d, length %d, update_every %f, total_time %f\n",p->nregions,p->nlegs,p->legtime,p->length,p->update_every,p->total_time);}
    p->regionra = (struct region **) tcalloc(p->nregions,sizeof(struct region*));
    if (verbose){ printf(" %% now making regions\n");}
    regionramake(p,event_within,event_threshold,region_type);
    assert(p->length>p->nlegs*p->legtime);
    p->eventra = (struct llitem **) tmalloc(p->length*sizeof(struct llitem *)); for (nr=0;nr<p->length;nr++){ p->eventra[nr]=llitemmake();}
    p->tab=0;
    p->update_last=GLOBAL_TI;
    if (verbose){ printf(" %% now making pretree and postree\n");}
    p->pretree=llitemmake();
    p->postree=llitemmake();
    p->wh=histmake(1,+1,-1);
    p->rh=histmake(1,+1,-1);
    exit_flag=0;
    while(exit_flag!=3){
      if ((readout = fread(&weight,sizeof(double),1,fp))==1 && (readout = fread(&relevance,sizeof(double),1,fp))==1){
	if (verbose){ printf(" %% read weight %f, relevance %f\n",weight,relevance);}
	L=llistmake();exit_flag=0;
	while (!exit_flag){
	  label = (int *) tmalloc(sizeof(int));
	  readout = fread(label,sizeof(int),1,fp);
	  if (readout!=1){ tfree(label); label=NULL; exit_flag=3;}
	  else{ 
	    if (*label>=0){ litemadd(L,label); exit_flag=0;}
	    else if (*label==-1){ tfree(label); label=NULL; exit_flag=1;}
	    else if (*label==-2){ tfree(label); label=NULL; exit_flag=2;}
	    else{ printf("warning, funny terminator %d in ptreadback\n",*label);}}}
	if (L->length>0 && (exit_flag==1 || exit_flag==2)){ /* postree or pretree */ 
	  l=L->first;depth=1; if (exit_flag==1){ l0=p->postree;pn=NULL;} else{ l0=p->pretree;pn=NULL;}
	  while(l!=NULL && depth<L->length){ /* traverse through nonlast elements of llist *L */
	    label = (int *) l->item;
	    if ((l1=llitemaddorfind(0,l0,p->regionra[*label],&region2pnode_compare_label))==NULL){
	      l1 = llitemaddorfind(1,l0,pnodemake(pn,p->regionra[*label],0,0),&pnode2pnode_compare_label);}
	    pn = (struct pnode *) l1->item;
	    l0 = pn->childllitem;
	    l=l->child;
	    depth+=1;}
	  if (depth==L->length){ /* at end of llist *L, must add pnode */
	    label = (int *) l->item;
	    if ((l1=llitemaddorfind(0,l0,p->regionra[*label],&region2pnode_compare_label))==NULL){
	      l1 = llitemaddorfind(1,l0,pnodemake(pn,p->regionra[*label],weight,relevance),&pnode2pnode_compare_label);}
	    else{ pn = (struct pnode *)l1->item; pn->weight += weight; pn->relevance=relevance;}}}
	llisttfree3(L);L=NULL;}
      else{ if (verbose){ printf(" %% couldn't read weight and relevance\n");} exit_flag=3;}}
    if (fp!=stdout){ fclose(fp);}
    if (verbose){ printf(" %% now balancing pretree and postree\n");}
    llitembalance(p->pretree); p->pretree = llitemclimb(p->pretree);
    llitembalance(p->postree); p->postree = llitemclimb(p->postree);}
  return p;
}

void ptreeplusequals_starter(struct ptree *p0,struct ptree *p1)
{
  ptreeplusequals_helper(p0,NULL,p0->postree,p1,NULL,p1->postree);
  ptreeplusequals_helper(p0,NULL,p0->pretree,p1,NULL,p1->pretree);
  pnodebalance_starter(NULL,p0->postree);
  pnodebalance_starter(NULL,p0->pretree);
  p0->total_time += p1->total_time;
}

void ptreeplusequals_helper(struct ptree *p0,struct pnode *parent0,struct llitem *l0,struct ptree *p1,struct pnode *parent1,struct llitem *l1)
{
  /* *p0 += *p1 */
  int firstentry=0;
  struct pnode *pn0=NULL,*pn1=NULL;
  struct llitem *l2=NULL;
  if (parent0!=NULL){ assert(parent0->childllitem==l0);}
  if (parent1!=NULL && l1->parent==NULL){ /* first descent into parent1->childllitem */ 
    assert(parent1->childllitem==l1); firstentry=1;}
  if (l1->kidl!=NULL){ ptreeplusequals_helper(p0,parent0,l0,p1,parent1,l1->kidl);}
  if (l1->kidr!=NULL){ ptreeplusequals_helper(p0,parent0,l0,p1,parent1,l1->kidr);}
  if (l1->item!=NULL){ 
    pn1 = (struct pnode *)l1->item; 
    if ((l2 = llitemaddorfind(0,l0,p0->regionra[pn1->region->label],&region2pnode_compare_label))==NULL){
      l2 = llitemaddorfind(1,l0,pnodemake(parent0,p0->regionra[pn1->region->label],0,0),&pnode2pnode_compare_label);}
    pn0 = (struct pnode *)l2->item;
    ptreeplusequals_helper(p0,pn0,pn0->childllitem,p1,pn1,pn1->childllitem);}
  if (firstentry){ 
    parent0->weight += parent1->weight;
    parent0->relevance = parent1->relevance;}
}

struct ptree * ptreesetequalto_starter(struct ptree *p1)
{
  struct ptree *p=NULL;
  int nr=0;
  p = (struct ptree *) tmalloc(sizeof(struct ptree));
  p->nregions = p1->nregions;
  p->regionra = (struct region **) tcalloc(p->nregions,sizeof(struct region*));
  for (nr=0;nr<p->nregions;nr++){ 
    p->regionra[nr] = (struct region *) tmalloc(sizeof(struct region)); 
    p->regionra[nr]->label = nr; /* critical to have label nr match location in regionra */
    p->regionra[nr]->last_event = p1->regionra[nr]->last_event;
    p->regionra[nr]->event_within = p1->regionra[nr]->event_within;
    p->regionra[nr]->event_threshold = p1->regionra[nr]->event_threshold;
    p->regionra[nr]->neuronllist = llistcopy(p1->regionra[nr]->neuronllist);
    p->regionra[nr]->pn = p1->regionra[nr]->pn;}
  p->nlegs=p1->nlegs;
  p->legtime=p1->legtime;
  p->length = p1->length;
  p->eventra = (struct llitem **) tmalloc(p->length*sizeof(struct llitem *));
  for (nr=0;nr<p->length;nr++){ 
    p->eventra[nr]=llitemmake(); 
    llitemgrowllitem(p->eventra[nr],p1->eventra[nr],&region2region_compare_label); 
    llitembalance(p->eventra[nr]);
    p->eventra[nr]=llitemclimb(p->eventra[nr]);}
  p->tab=p1->tab;
  p->update_every=p1->update_every;
  p->update_last=p1->update_last;
  p->total_time=p1->total_time;
  p->gate_flag=p1->gate_flag;
  ptreemakepospre(p);
  p->wh = histmake(1,+1,-1);
  p->rh = histmake(1,+1,-1);
  ptreesetequalto_helper(p,NULL,p->postree,p1,NULL,p1->postree);
  ptreesetequalto_helper(p,NULL,p->pretree,p1,NULL,p1->pretree);
  pnodebalance_starter(NULL,p->postree);
  pnodebalance_starter(NULL,p->pretree);
  return p;
}

void ptreesetequalto_helper(struct ptree *p0,struct pnode *parent0,struct llitem *l0,struct ptree *p1,struct pnode *parent1,struct llitem *l1)
{
  /* *p0 = *p1 */
  int firstentry=0;
  struct pnode *pn0=NULL,*pn1=NULL;
  struct llitem *l2=NULL;
  if (parent0!=NULL){ assert(parent0->childllitem==l0);}
  if (parent1!=NULL && l1->parent==NULL){ /* first descent into parent1->childllitem */ 
    assert(parent1->childllitem==l1); firstentry=1;}
  if (l1->kidl!=NULL){ ptreesetequalto_helper(p0,parent0,l0,p1,parent1,l1->kidl);}
  if (l1->kidr!=NULL){ ptreesetequalto_helper(p0,parent0,l0,p1,parent1,l1->kidr);}
  if (l1->item!=NULL){ 
    pn1 = (struct pnode *)l1->item; 
    if ((l2 = llitemaddorfind(0,l0,p0->regionra[pn1->region->label],&region2pnode_compare_label))==NULL){
      l2 = llitemaddorfind(1,l0,pnodemake(parent0,p0->regionra[pn1->region->label],0,0),&pnode2pnode_compare_label);}
    pn0 = (struct pnode *)l2->item;
    ptreesetequalto_helper(p0,pn0,pn0->childllitem,p1,pn1,pn1->childllitem);}
  if (firstentry){ 
    parent0->weight = parent1->weight;
    parent0->relevance = parent1->relevance;}
}

struct ptree * ptreesubtptree_starter(struct ptree *px,struct ptree *py,int ratematch,int reldif_flag)
{
  /* ratematch matches l2 norm of leg-0, reldif_flag computes relative difference per node */
  struct llist *pnLx=NULL,*pnLy=NULL;
  struct litem *lx=NULL,*ly=NULL;
  struct pnode *pnx=NULL,*pny=NULL;
  double xy=0,yy=0,ynormalizer=0;
  struct ptree *pz=NULL;
  if (ratematch){
    pnLx=llistmake(); pnLy=llistmake();
    pnode2llist_starter(NULL,px->postree,1,0,pnLx); pnode2llist_starter(NULL,py->postree,1,0,pnLy);
    llistsort(pnLx->first,pnLx->last,pnLx->length,&pnode2pnode_compare_label);
    llistsort(pnLy->first,pnLy->last,pnLy->length,&pnode2pnode_compare_label);
    lx=pnLx->first; ly=pnLy->first; xy=0; yy=0;
    while (lx!=NULL && ly!=NULL){
      pnx=(struct pnode *)lx->item; pny=(struct pnode *)ly->item;
      xy += pnx->weight*pny->weight; yy += pny->weight*pny->weight;      
      lx=lx->child; ly=ly->child;}
    llisttfree(pnLx); pnLx=NULL; llisttfree(pnLy); pnLy=NULL;
    ynormalizer=xy/yy; if (!finite(ynormalizer)){ ynormalizer=1;}}
  else /* if (!ratematch) */{ ynormalizer=1;}
  pz = ptreesetequalto_starter(px);
  ptreesubtptree_helper(pz,NULL,pz->postree,py,NULL,py->postree,ynormalizer,reldif_flag);
  ptreesubtptree_helper(pz,NULL,pz->pretree,py,NULL,py->pretree,ynormalizer,reldif_flag);
  if (reldif_flag){
    ptreecomplement_helper(py,NULL,py->postree,pz,NULL,pz->postree,1,0);
    ptreecomplement_helper(py,NULL,py->pretree,pz,NULL,pz->pretree,1,0);}
  return pz;
}

void ptreecomplement_helper(struct ptree *px,struct pnode *parentx,struct llitem *lx,struct ptree *py,struct pnode *parenty,struct llitem *ly,double set_weight,double set_relevance)
{
  /* every node in *py which doesn't have a corresponding node in *py gets set to set_weight,set_relevance */
  int firstentry=0;
  struct pnode *pnx=NULL,*pny=NULL;
  struct llitem *l2=NULL;
  if (parentx!=NULL){ assert(parentx->childllitem==lx);}
  if (parenty!=NULL && ly->parent==NULL){ /* first descent into parenty->childllitem */ 
    assert(parenty->childllitem==ly); firstentry=1;}
  if (ly->kidl!=NULL){ ptreecomplement_helper(px,parentx,lx,py,parenty,ly->kidl,set_weight,set_relevance);}
  if (ly->kidr!=NULL){ ptreecomplement_helper(px,parentx,lx,py,parenty,ly->kidr,set_weight,set_relevance);}
  if (ly->item!=NULL){ 
    pny = (struct pnode *)ly->item; 
    if (lx!=NULL){
      if ((l2 = llitemaddorfind(0,lx,px->regionra[pny->region->label],&region2pnode_compare_label))!=NULL){
	pnx = (struct pnode *)l2->item;
	ptreecomplement_helper(px,pnx,pnx->childllitem,py,pny,pny->childllitem,set_weight,set_relevance);}
      else{ ptreecomplement_helper(px,NULL,NULL,py,pny,pny->childllitem,set_weight,set_relevance);}}
    else /* if (lx==NULL) */{ ptreecomplement_helper(px,NULL,NULL,py,pny,pny->childllitem,set_weight,set_relevance);}}
  if (firstentry){ 
    if (parentx==NULL && lx==NULL){ parenty->weight = set_weight; parenty->relevance = set_relevance;}}
}

void ptreesubtptree_helper(struct ptree *px,struct pnode *parentx,struct llitem *lx,struct ptree *py,struct pnode *parenty,struct llitem *ly,double ynormalizer,int reldif_flag)
{
  /* *px -= *py*ynormalizer */
  int firstentry=0;
  struct pnode *pnx=NULL,*pny=NULL;
  struct llitem *l2=NULL;
  if (parentx!=NULL){ assert(parentx->childllitem==lx);}
  if (parenty!=NULL && ly->parent==NULL){ /* first descent into parenty->childllitem */ 
    assert(parenty->childllitem==ly); firstentry=1;}
  if (ly->kidl!=NULL){ ptreesubtptree_helper(px,parentx,lx,py,parenty,ly->kidl,ynormalizer,reldif_flag);}
  if (ly->kidr!=NULL){ ptreesubtptree_helper(px,parentx,lx,py,parenty,ly->kidr,ynormalizer,reldif_flag);}
  if (ly->item!=NULL){ 
    pny = (struct pnode *)ly->item; 
    if ((l2 = llitemaddorfind(0,lx,px->regionra[pny->region->label],&region2pnode_compare_label))==NULL){
      l2 = llitemaddorfind(1,lx,pnodemake(parentx,px->regionra[pny->region->label],0,0),&pnode2pnode_compare_label);}
    pnx = (struct pnode *)l2->item;
    ptreesubtptree_helper(px,pnx,pnx->childllitem,py,pny,pny->childllitem,ynormalizer,reldif_flag);}
  if (firstentry){ 
    if (reldif_flag){ parentx->weight = fabs(parentx->weight-parenty->weight*ynormalizer)/maximum(1,maximum(parenty->weight*ynormalizer,parentx->weight));}
    else /* if (!reldif_flag) */{ parentx->weight -= parenty->weight*ynormalizer;}
    parentx->relevance = 0;}
}

void pnodefun_starter(struct pnode *parent,struct llitem *l0,int maxlevel,int minlevel,int level,void (*fun)(double *,double *))
{
  /* performs (*fun)(&(parent->weight),&(parent->relevance)) for each node */
  int firstentry=0;
  struct pnode *pn0=NULL;
  if (parent!=NULL && l0->parent==NULL){ /* first descent into parent->childllitem */ 
    assert(parent->childllitem==l0); firstentry=1;}
  if (l0->kidl!=NULL){ pnodefun_starter(parent,l0->kidl,maxlevel,minlevel,level,fun);}
  if (l0->kidr!=NULL){ pnodefun_starter(parent,l0->kidr,maxlevel,minlevel,level,fun);}
  if (l0->item!=NULL){ 
    if (maxlevel==-1 || level<maxlevel){ 
      pn0 = (struct pnode *)l0->item; 
      pnodefun_starter(pn0,pn0->childllitem,maxlevel,minlevel,level+1,fun);}}
  if (firstentry){ if (minlevel==-1 || level>=minlevel){ if (fun!=NULL){ (*fun)(&(parent->weight),&(parent->relevance));}}}
}

void pnodefrob_starter(struct pnode *parent,struct llitem *l0,int maxlevel,int level,double *wnorm,double *rnorm)
{
  /* returns frobenius norm */
  int firstentry=0;
  struct pnode *pn0=NULL;
  if (parent!=NULL && l0->parent==NULL){ /* first descent into parent->childllitem */ 
    assert(parent->childllitem==l0); firstentry=1;}
  if (l0->kidl!=NULL){ pnodefrob_starter(parent,l0->kidl,maxlevel,level,wnorm,rnorm);}
  if (l0->kidr!=NULL){ pnodefrob_starter(parent,l0->kidr,maxlevel,level,wnorm,rnorm);}
  if (l0->item!=NULL){ 
    if (maxlevel==-1 || level<maxlevel){ 
      pn0 = (struct pnode *)l0->item; 
      pnodefrob_starter(pn0,pn0->childllitem,maxlevel,level+1,wnorm,rnorm);}}
  if (firstentry){ 
    if (wnorm!=NULL){ *wnorm += parent->weight*parent->weight;}
    if (rnorm!=NULL){ *rnorm += parent->relevance*parent->relevance;}}
}

void ptreex2(struct ptree *p0,struct pnode *parent0,struct llitem *l0,struct ptree *p1,struct pnode *parent1,struct llitem *l1,int retain_self,int maxlevel,int level,double *output_x2,double *output_x2p,double *output_fabs,double *output_fabs2)
{
  /* computes x2 estimate with p0 expected and p1 sample 
     returns estimators for chi_squared, forced_finite_chi_squared, fabsenius, squared_fabsenius */
  int verbose=0,nr=0;
  char lvlchar[64];
  int firstentry=0;
  struct pnode *pn0=NULL,*pn1=NULL;
  struct llitem *l2=NULL;
  double weight_expected=0,weight_sample=0,weight_minimum=0;
  if (verbose){ sprintf(lvlchar," %%"); for (nr=0;nr<level;nr++){sprintf(lvlchar,"%s %%",lvlchar);}}
  if (verbose>1){ 
    printf(" %s [entering ptreex2]\n",lvlchar);
    printf(" %s p0 %d parent0 %d->%d l0 %d \n",lvlchar,(int)p0,(int)parent0,parent0==NULL ? 0 : (int)(parent0->childllitem),(int)l0);
    printf(" %s p1 %d parent1 %d->%d l1 %d\n",lvlchar,(int)p1,(int)parent1,parent1==NULL ? 0 : (int)(parent1->childllitem),(int)l1);
    printf(" %s retain_self %d maxlevel %d level %d\n",lvlchar,retain_self,maxlevel,level);
    if (output_x2!=NULL){ printf(" %s output_x2 %f\n",lvlchar,*output_x2);}
    if (output_x2p!=NULL){ printf(" %s output_x2p %f\n",lvlchar,*output_x2p);}
    if (output_fabs!=NULL){ printf(" %s output_fabs %f\n",lvlchar,*output_fabs);}
    if (output_fabs2!=NULL){ printf(" %s output_fabs2 %f\n",lvlchar,*output_fabs2);}}
  if (parent1!=NULL){ assert(parent1->childllitem==l1);}
  if (parent0!=NULL && l0->parent==NULL){ /* first descent into parent0->childllitem */ 
    if (verbose>1){ printf(" %s first descent into parent0->childllitem, setting firstentry to 1\n",lvlchar);}
    assert(parent0->childllitem==l0); firstentry=1;}
  if (l0->kidl!=NULL){ 
    if (verbose>1){ printf(" %s l0->kidl found, moving left to new llitem l0->kidl=%d\n",lvlchar,(int)(l0->kidl));}
    ptreex2(p0,parent0,l0->kidl,p1,parent1,l1,retain_self,maxlevel,level,output_x2,output_x2p,output_fabs,output_fabs2);}
  if (l0->kidr!=NULL){ 
    if (verbose>1){ printf(" %s l0->kidr found, moving right to new llitem l0->kidr=%d\n",lvlchar,(int)(l0->kidr));}
    ptreex2(p0,parent0,l0->kidr,p1,parent1,l1,retain_self,maxlevel,level,output_x2,output_x2p,output_fabs,output_fabs2);}
  if (l0->item!=NULL){ 
    if (verbose>1){ printf(" %s l0->item %d exists\n",lvlchar,(int)(struct pnode *)(l0->item));}
    pn0 = (struct pnode *)l0->item; 
    if (maxlevel==-1 || level<maxlevel){
      if (verbose>1){ printf(" %s level %d less than maxlevel %d, continuing...\n",lvlchar,level,maxlevel);}
      if (l1!=NULL && (l2=llitemaddorfind(0,l1,p1->regionra[pn0->region->label],&region2pnode_compare_label))!=NULL){
	pn1 = (struct pnode *)l2->item;
	if (verbose>1){ printf(" %s pn0 label %d matches pn1 label %d, descending...\n",lvlchar,pn0->region->label,pn1->region->label);}
	ptreex2(p0,pn0,pn0->childllitem,p1,pn1,pn1->childllitem,retain_self,maxlevel,level+1,output_x2,output_x2p,output_fabs,output_fabs2);}
      else{
	pn1 = NULL;
	if (verbose>1){ printf(" %s pn0 label %d not found within l1, descending...\n",lvlchar,pn0->region->label);}
	ptreex2(p0,pn0,pn0->childllitem,p1,NULL,NULL,retain_self,maxlevel,level+1,output_x2,output_x2p,output_fabs,output_fabs2);}}}
  if (firstentry){ 
    if (!retain_self){ 
      weight_expected = (parent0->weight - (parent1==NULL ? 0 : parent1->weight))*p1->total_time/(p0->total_time-p1->total_time); 
      weight_minimum = 1.0*p1->total_time/(p0->total_time-p1->total_time);}
    else /* if (retain_self) */{
      weight_expected = parent0->weight*p1->total_time/p0->total_time; weight_minimum = 1.0*p1->total_time/p0->total_time;}
    if (verbose>0){ printf(" %s first entry, and weight_expected = %f...",lvlchar,weight_expected);}
    if (parent1!=NULL){ 
      if (verbose>0){ printf(" parent1 exists, so using sample weight %f\n",parent1->weight);}
      weight_sample = parent1->weight;}
    else{       
      if (verbose>0){ printf(" parent1 does not exist, so using sample weight 0\n");}
      weight_sample = 0;}
    if (output_x2!=NULL){ *output_x2 += pow(weight_sample - weight_expected,2)/pow((weight_expected > 0 ? weight_expected : 1),2);}
    if (output_x2p!=NULL){ *output_x2p += pow(weight_sample - weight_expected,2)/maximum(weight_minimum,weight_expected);}
    if (output_fabs!=NULL){ *output_fabs += fabs(weight_sample - weight_expected);}
    if (output_fabs2!=NULL){ *output_fabs2 += pow(weight_sample - weight_expected,2);}}
}

void ptreerelent_breadth(struct ptree *p0,struct pnode *parent0,struct llitem *l0,struct ptree *p1,struct pnode *parent1,struct llitem *l1,int retain_self,int maxlevel,int minlevel,int level,double *output,double *rnormalizer,double *onormalizer)
{
  /* computes "breadthwise" relative entropy with p0 expected and p1 observed 
     this assumes that rnormalizer and onormalizer are both maxlevel long, 
     and contain the mean weights at each level */
  int verbose=0,nr=0;
  char lvlchar[64];
  int firstentry=0;
  struct pnode *pn0=NULL,*pn1=NULL;
  struct llitem *l2=NULL;
  double weight_expected=0,weight_sample=0;
  if (verbose){ sprintf(lvlchar," %%"); for (nr=0;nr<level;nr++){sprintf(lvlchar,"%s %%",lvlchar);}}
  if (verbose>1){ 
    printf(" %s [entering ptreerelent_breadth]\n",lvlchar);
    printf(" %s p0 %d parent0 %d->%d l0 %d \n",lvlchar,(int)p0,(int)parent0,parent0==NULL ? 0 : (int)(parent0->childllitem),(int)l0);
    printf(" %s p1 %d parent1 %d->%d l1 %d\n",lvlchar,(int)p1,(int)parent1,parent1==NULL ? 0 : (int)(parent1->childllitem),(int)l1);
    printf(" %s retain_self %d maxlevel %d minlevel %d level %d\n",lvlchar,retain_self,maxlevel,minlevel,level);
    if (output!=NULL){ printf(" %s output %f\n",lvlchar,*output);}
    printf(" rnormalizer:\n"); raprintf(rnormalizer,"double",1,maxlevel,lvlchar);
    printf(" onormalizer:\n"); raprintf(onormalizer,"double",1,maxlevel,lvlchar);}
  if (parent1!=NULL){ assert(parent1->childllitem==l1);}
  if (parent0!=NULL && l0->parent==NULL){ /* first descent into parent0->childllitem */ 
    if (verbose>1){ printf(" %s first descent into parent0->childllitem, setting firstentry to 1\n",lvlchar);}
    assert(parent0->childllitem==l0); firstentry=1;}
  if (l0->kidl!=NULL){ 
    if (verbose>1){ printf(" %s l0->kidl found, moving left to new llitem l0->kidl=%d\n",lvlchar,(int)(l0->kidl));}
    ptreerelent_breadth(p0,parent0,l0->kidl,p1,parent1,l1,retain_self,maxlevel,minlevel,level,output,rnormalizer,onormalizer);}
  if (l0->kidr!=NULL){ 
    if (verbose>1){ printf(" %s l0->kidr found, moving right to new llitem l0->kidr=%d\n",lvlchar,(int)(l0->kidr));}
    ptreerelent_breadth(p0,parent0,l0->kidr,p1,parent1,l1,retain_self,maxlevel,minlevel,level,output,rnormalizer,onormalizer);}
  if (l0->item!=NULL){ 
    if (verbose>1){ printf(" %s l0->item %d exists\n",lvlchar,(int)(struct pnode *)(l0->item));}
    pn0 = (struct pnode *)l0->item; 
    if (maxlevel==-1 || level<maxlevel){
      if (verbose>1){ printf(" %s level %d less than maxlevel %d, continuing...\n",lvlchar,level,maxlevel);}
      if (l1!=NULL && (l2=llitemaddorfind(0,l1,p1->regionra[pn0->region->label],&region2pnode_compare_label))!=NULL){
	pn1 = (struct pnode *)l2->item;
	if (verbose>1){ printf(" %s pn0 label %d matches pn1 label %d, descending...\n",lvlchar,pn0->region->label,pn1->region->label);}
	ptreerelent_breadth(p0,pn0,pn0->childllitem,p1,pn1,pn1->childllitem,retain_self,maxlevel,minlevel,level+1,output,rnormalizer,onormalizer);}
      else{
	pn1 = NULL;
	if (verbose>1){ printf(" %s pn0 label %d not found within l1, descending...\n",lvlchar,pn0->region->label);}
	ptreerelent_breadth(p0,pn0,pn0->childllitem,p1,NULL,NULL,retain_self,maxlevel,minlevel,level+1,output,rnormalizer,onormalizer);}}}
  if (firstentry){ if (minlevel==-1 || level>=minlevel){
    if (!retain_self){ 
      weight_expected = maximum(1,(parent0->weight - (parent1==NULL ? 0 : parent1->weight)))/maximum(1,rnormalizer[level-1]-onormalizer[level-1]);}
    else /* if (retain_self) */{
      weight_expected = maximum(1,parent0->weight)/maximum(1,rnormalizer[level-1]);}
    if (verbose>0){ printf(" %s first entry, and weight_expected = %f...",lvlchar,weight_expected);}
    if (parent1!=NULL){ 
      weight_sample = parent1->weight/maximum(1,onormalizer[level-1]);
      if (verbose>0){ printf(" parent1 exists, so using sample weight %f\n",weight_sample);}}
    else{       
      if (verbose>0){ printf(" parent1 does not exist, so using sample weight 0\n");}
      weight_sample = 0;}
    if (output!=NULL){ *output += weight_sample>0 ? weight_sample*log(weight_sample/weight_expected) : 0;}}}
}

void ptreerelent_depth(struct ptree *p0,struct pnode *parent0,struct llitem *l0,struct ptree *p1,struct pnode *parent1,struct llitem *l1,int retain_self,int maxlevel,int level,double *output)
{
  /* computes "depthwise" relative entropy with p0 expected and p1 observed */
  int verbose=0,nr=0;
  char lvlchar[64];
  int firstentry=0;
  struct pnode *pn0=NULL,*pn1=NULL;
  struct llitem *l2=NULL;
  double weight_expected_prior=0,weight_expected=0,weight_observed_prior=0,weight_observed=0;
  double p_expected_prior=0,p_expected=0,p_observed_prior=0,p_observed=0;
  int maxdepth = (maxlevel==-1 ? p0->nlegs : minimum(p0->nlegs,maxlevel));
  if (verbose){ sprintf(lvlchar," %%"); for (nr=0;nr<level;nr++){sprintf(lvlchar,"%s %%",lvlchar);}}
  if (verbose>1){ 
    printf(" %s [entering ptreerelent_depth]\n",lvlchar);
    printf(" %s p0 %d parent0 %d->%d l0 %d \n",lvlchar,(int)p0,(int)parent0,parent0==NULL ? 0 : (int)(parent0->childllitem),(int)l0);
    printf(" %s p1 %d parent1 %d->%d l1 %d\n",lvlchar,(int)p1,(int)parent1,parent1==NULL ? 0 : (int)(parent1->childllitem),(int)l1);
    printf(" %s retain_self %d maxlevel %d level %d (maxdepth %d)\n",lvlchar,retain_self,maxlevel,level,maxdepth);
    if (output!=NULL){ printf(" %s output %f\n",lvlchar,*output);}}
  if (parent1!=NULL){ assert(parent1->childllitem==l1);}
  if (parent0!=NULL && l0->parent==NULL){ /* first descent into parent0->childllitem */ 
    if (verbose>1){ printf(" %s first descent into parent0->childllitem, setting firstentry to 1\n",lvlchar);}
    assert(parent0->childllitem==l0); firstentry=1;}
  if (l0->kidl!=NULL){ 
    if (verbose>1){ printf(" %s l0->kidl found, moving left to new llitem l0->kidl=%d\n",lvlchar,(int)(l0->kidl));}
    ptreerelent_depth(p0,parent0,l0->kidl,p1,parent1,l1,retain_self,maxlevel,level,output);}
  if (l0->kidr!=NULL){ 
    if (verbose>1){ printf(" %s l0->kidr found, moving right to new llitem l0->kidr=%d\n",lvlchar,(int)(l0->kidr));}
    ptreerelent_depth(p0,parent0,l0->kidr,p1,parent1,l1,retain_self,maxlevel,level,output);}
  if (l0->item!=NULL){ 
    if (verbose>1){ printf(" %s l0->item %d exists\n",lvlchar,(int)(struct pnode *)(l0->item));}
    pn0 = (struct pnode *)l0->item; 
    if (maxlevel==-1 || level<maxlevel){
      if (verbose>1){ printf(" %s level %d less than maxlevel %d, continuing...\n",lvlchar,level,maxlevel);}
      if (l1!=NULL && (l2=llitemaddorfind(0,l1,p1->regionra[pn0->region->label],&region2pnode_compare_label))!=NULL){
	pn1 = (struct pnode *)l2->item;
	if (verbose>1){ printf(" %s pn0 label %d matches pn1 label %d, descending...\n",lvlchar,pn0->region->label,pn1->region->label);}
	ptreerelent_depth(p0,pn0,pn0->childllitem,p1,pn1,pn1->childllitem,retain_self,maxlevel,level+1,output);}
      else{
	pn1 = NULL;
	if (verbose>1){ printf(" %s pn0 label %d not found within l1, descending...\n",lvlchar,pn0->region->label);}
	ptreerelent_depth(p0,pn0,pn0->childllitem,p1,NULL,NULL,retain_self,maxlevel,level+1,output);}}}
  if (firstentry){ 
    if (parent0->parent==NULL){ 
      assert(level==1);
      weight_expected_prior = p0->total_time;
      if (verbose>1){ printf(" %s %d==1 level event, no parent, must be a firing rate, weight_expected_prior=%f\n",lvlchar,level,weight_expected_prior);}}
    else{
      weight_expected_prior = parent0->parent->weight;
      if (verbose>1){ printf(" %s %d>1 level event, parent exists, weight_expected_prior=%f\n",lvlchar,level,weight_expected_prior);}}
    if (level==1){
      weight_observed_prior=p1->total_time;
      if (verbose>1){ printf(" %s level 1, must use weight_observed_prior=%f\n",lvlchar,weight_observed_prior);}}
    else{
      if (parent1!=NULL){
	assert(parent1->parent!=NULL);
	weight_observed_prior=parent1->parent->weight;
	if (verbose>1){ printf(" %s level %d, parent exists, weight_observed_prior=%f\n",lvlchar,level,weight_observed_prior);}}
      else{ 
	weight_observed_prior=0;
	if (verbose>1){ printf(" %s level %d, parent doesn't exist, weight_observed_prior=%f\n",lvlchar,level,weight_observed_prior);}}}
    weight_expected = parent0->weight;
    if (verbose>1){ printf(" %s first entry, weight_expected=%f...",lvlchar,weight_expected);}
    if (parent1!=NULL){ 
      weight_observed = parent1->weight;
      if (verbose>1){ printf(" %s parent1 exists, so weight_observed %f\n",lvlchar,weight_observed);}}
    else{       
      weight_observed = 0;
      if (verbose>1){ printf(" %s parent1 does not exist, so weight_observed %f\n",lvlchar,weight_observed);}}
    if (!retain_self){
      p_expected = minimum(1,maximum(0,(weight_expected-weight_observed)/(p0->total_time-p1->total_time)));
      p_expected_prior = minimum(1,maximum(0,((weight_expected_prior-weight_observed_prior)-(weight_expected-weight_observed))/(p0->total_time-p1->total_time)));
      p_observed = minimum(1,maximum(0,weight_observed/p1->total_time));
      p_observed_prior = minimum(1,maximum(0,(weight_observed_prior-weight_observed)/p1->total_time));}
    else /* if (retain_self) */{
      p_expected = minimum(1,maximum(0,weight_expected/p0->total_time));
      p_expected_prior = minimum(1,maximum(0,(weight_expected_prior-weight_expected)/p0->total_time));
      p_observed = minimum(1,maximum(0,weight_observed/p1->total_time));
      p_observed_prior = minimum(1,maximum(0,(weight_observed_prior-weight_observed)/p1->total_time));}
    if (verbose>1){ printf(" %s probabilities e %f ep %f o %f op %f,\n",lvlchar,p_expected,p_expected_prior,p_observed,p_observed_prior);}
    if (p_observed_prior>0 && p_expected_prior>0){ if (output!=NULL){ 
      *output += pow(p0->nregions,maxdepth-level)*p_observed_prior*log(p_observed_prior/p_expected_prior);}}
    if (level==maxdepth){ if (p_observed>0 && p_expected>0){ if (output!=NULL){ 
      *output += p_observed*log(p_observed/p_expected);}}}
    else if (level<maxdepth){ if (p_observed>0 && p_expected>0){ if (output!=NULL){ 
      *output += pow(p0->nregions,maxdepth-(level+1))*(p0->nregions-llitemlength(l1))*p_observed*log(p_observed/p_expected);}}}}
}

void pnodeclear_starter(struct pnode *parent,struct llitem *l0)
{
  int firstentry=0;
  struct pnode *pn=NULL;
  if (parent!=NULL && l0->parent==NULL){ /* first descent into parent->childllitem */ 
    assert(parent->childllitem==l0); firstentry=1;}
  if (l0->kidl!=NULL){ pnodeclear_starter(parent,l0->kidl);}
  if (l0->kidr!=NULL){ pnodeclear_starter(parent,l0->kidr);}
  if (l0->item!=NULL){ pn=(struct pnode *)l0->item; pnodeclear_starter(pn,pn->childllitem);}
  if (firstentry){ if (parent!=NULL){ parent->weight=0; parent->relevance=0;}}
}

void pnodehist_starter(struct pnode *parent,struct llitem *l0,int maxlevel,int level,struct hist *histw,struct hist *histr)
{
  /* puts weights and relevances into histw and histr */
  int firstentry=0;
  struct pnode *pn=NULL;
  if (parent!=NULL && l0->parent==NULL){ /* first descent into parent->childllitem */ 
    assert(parent->childllitem==l0); firstentry=1;}
  if (l0->kidl!=NULL){ pnodehist_starter(parent,l0->kidl,maxlevel,level,histw,histr);}
  if (l0->kidr!=NULL){ pnodehist_starter(parent,l0->kidr,maxlevel,level,histw,histr);}
  if (l0->item!=NULL){
    if (maxlevel==-1 || level<maxlevel){ 
      pn=(struct pnode *)l0->item; 
      pnodehist_starter(pn,pn->childllitem,maxlevel,level+1,histw,histr);}}
  if (firstentry){ if (parent!=NULL){ 
    if (histw!=NULL){ histadd(histw,parent->weight,1);} if (histr!=NULL){ histadd(histr,parent->relevance,1);}}}
}

void pnodestats_starter(struct pnode *parent,struct llitem *l0,int maxlevel,int minlevel,int level,int *nelements,double *wmax,double *wmin,double *wmean,double *wstdev,double *rmax,double *rmin,double *rmean,double *rstdev)
{
  /* start out with nelements==NULL, stdev only accumulates if mean!=NULL */
  int firstentry=0,veryfirstentry=0,madenelements_flag=0;
  struct pnode *pn=NULL;
  if (parent==NULL && l0->parent==NULL){ /* very first descent into ptree */
    if (nelements==NULL){ madenelements_flag=1; nelements=(int *)tmalloc(sizeof(int));} else{ madenelements_flag=0;} *nelements=0;
    if (wmin!=NULL){ *wmin=0;} if (wmax!=NULL){ *wmax=0;} if (wmean!=NULL){ *wmean=0; if (wstdev!=NULL){ *wstdev=0;}}
    if (rmin!=NULL){ *rmin=0;} if (rmax!=NULL){ *rmax=0;} if (rmean!=NULL){ *rmean=0; if (rstdev!=NULL){ *rstdev=0;}}
    veryfirstentry=1;}
  if (parent!=NULL && l0->parent==NULL){ /* first descent into parent->childllitem */ 
    assert(parent->childllitem==l0); assert(nelements!=NULL); firstentry=1;}
  if (l0->kidl!=NULL){ pnodestats_starter(parent,l0->kidl,maxlevel,minlevel,level,nelements,wmax,wmin,wmean,wstdev,rmax,rmin,rmean,rstdev);}
  if (l0->kidr!=NULL){ pnodestats_starter(parent,l0->kidr,maxlevel,minlevel,level,nelements,wmax,wmin,wmean,wstdev,rmax,rmin,rmean,rstdev);}
  if (l0->item!=NULL){
    if (maxlevel==-1 || level<maxlevel){ 
      pn=(struct pnode *)l0->item; 
      pnodestats_starter(pn,pn->childllitem,maxlevel,minlevel,level+1,nelements,wmax,wmin,wmean,wstdev,rmax,rmin,rmean,rstdev);}}
  if (firstentry){ if (parent!=NULL){ if (nelements!=NULL){ if (minlevel==-1 || level>=minlevel){
    if (*nelements==0){ if (wmax!=NULL){ *wmax=parent->weight;} if (wmin!=NULL){ *wmin=parent->weight;}}
    else if (*nelements>0){ if (wmax!=NULL){ *wmax=maximum(*wmax,parent->weight);} if (wmin!=NULL){ *wmin=minimum(*wmin,parent->weight);}}
    if (wmean!=NULL){ *wmean += parent->weight; if (wstdev!=NULL){ *wstdev += pow(parent->weight,2);}}
    if (*nelements==0){ if (rmax!=NULL){ *rmax=parent->relevance;} if (rmin!=NULL){ *rmin=parent->relevance;}}
    else if (*nelements>0){ if (rmax!=NULL){ *rmax=maximum(*rmax,parent->relevance);} if (rmin!=NULL){ *rmin=minimum(*rmin,parent->relevance);}}
    if (rmean!=NULL){ *rmean += parent->relevance; if (rstdev!=NULL){ *rstdev += pow(parent->relevance,2);}}
    *nelements += 1;}}}}
  if (veryfirstentry){ if (nelements!=NULL){ 
    if (wmean!=NULL){ *wmean /= maximum(1.0,(double)(*nelements)); if (wstdev!=NULL){ *wstdev /= maximum(1.0,(double)(*nelements)); *wstdev -= pow(*wmean,2);}}
    if (rmean!=NULL){ *rmean /= maximum(1.0,(double)(*nelements)); if (rstdev!=NULL){ *rstdev /= maximum(1.0,(double)(*nelements)); *rstdev -= pow(*rmean,2); *rstdev = sqrt(*rstdev);}}
    if (madenelements_flag==1){ tfree(nelements);nelements=NULL;}}}
}

void pnode2llist_starter(struct pnode *parent,struct llitem *l0,int maxlevel,int level,struct llist *L)
{
  /* start out with empty llist *L */
  int firstentry=0;
  struct pnode *pn=NULL;
  if (parent!=NULL && l0->parent==NULL){ /* first descent into parent->childllitem */ 
    assert(parent->childllitem==l0); firstentry=1;}
  if (l0->kidl!=NULL){ pnode2llist_starter(parent,l0->kidl,maxlevel,level,L);}
  if (l0->kidr!=NULL){ pnode2llist_starter(parent,l0->kidr,maxlevel,level,L);}
  if (l0->item!=NULL){
    if (maxlevel==-1 || level<maxlevel){ 
      pn=(struct pnode *)l0->item; 
      pnode2llist_starter(pn,pn->childllitem,maxlevel,level+1,L);}}
  if (firstentry){ if (parent!=NULL){ litemadd(L,parent);}}
}

void pnode2hist_starter(struct pnode *parent,struct llitem *l0,int maxlevel,int level,struct hist *wh,struct hist *rh)
{
  /* adds weights to *wh and relevances to *rh */
  int firstentry=0;
  struct pnode *pn=NULL;
  if (parent!=NULL && l0->parent==NULL){ /* first descent into parent->childllitem */ 
    assert(parent->childllitem==l0); firstentry=1;}
  if (l0->kidl!=NULL){ pnode2hist_starter(parent,l0->kidl,maxlevel,level,wh,rh);}
  if (l0->kidr!=NULL){ pnode2hist_starter(parent,l0->kidr,maxlevel,level,wh,rh);}
  if (l0->item!=NULL){
    if (maxlevel==-1 || level<maxlevel){ 
      pn=(struct pnode *)l0->item; 
      pnode2hist_starter(pn,pn->childllitem,maxlevel,level+1,wh,rh);}}
  if (firstentry){ if (parent!=NULL){ 
    if (wh!=NULL){ histadd(wh,parent->weight,1);} 
    if (rh!=NULL && parent->relevance!=0){ histadd(rh,parent->relevance,1);}}}
}

void pnodetimesd_starter(struct pnode *parent,struct llitem *l0,double weightx,double relevancex)
{
  int firstentry=0;
  struct pnode *pn=NULL;
  if (parent!=NULL && l0->parent==NULL){ /* first descent into parent->childllitem */ 
    assert(parent->childllitem==l0); firstentry=1;}
  if (l0->kidl!=NULL){ pnodetimesd_starter(parent,l0->kidl,weightx,relevancex);}
  if (l0->kidr!=NULL){ pnodetimesd_starter(parent,l0->kidr,weightx,relevancex);}
  if (l0->item!=NULL){ pn=(struct pnode *)l0->item; pnodetimesd_starter(pn,pn->childllitem,weightx,relevancex);}
  if (firstentry){ if (parent!=NULL){ parent->weight*=weightx; parent->relevance*=relevancex;}}
}

double ptreex2_match(int nE,int nS,double *x1ra,int *minindex)
{
  /* checks to see if x1ra[nS] is the minimum value
     if x1ra[nS] matches nmin of the values, we pass back 1.0/nmin; */
  int ne=0,nmin=0;
  double min=0,wtol=0.000000001,wr=0;
  stats("double",x1ra,nE,NULL,&min,NULL,NULL);
  nmin=0; for (ne=0;ne<nE;ne++){ if (fabs(x1ra[ne]-min)<wtol){ nmin += 1; if (minindex!=NULL){ *minindex=ne;}}}
  if (fabs(x1ra[nS]-min)<wtol){ wr = 1.0/(double)nmin; /* (nmin==1); */} else{ wr=0;}
  if (minindex!=NULL && nmin>1){ *minindex=-1;}
  return wr;
}

double ptreex1_smatch_old(int nE,int nS,double *we,double wa,int *minindex)
{
  /* additive matching criterion for nE expected weights *we and actual weight wa, given that the actual index should be nS 
     if wa matches nmin of the maxima of *we, we pass back (nmin==1) */
  int ne=0,nmin=0;
  double *wf=NULL,min=0,wtol=0.000000001,wr=0;
  wf = (double *) tcalloc(nE,sizeof(double));
  for (ne=0;ne<nE;ne++){ wf[ne] = fabs(we[ne]-wa);}
  stats("double",wf,nE,NULL,&min,NULL,NULL);
  nmin=0; for (ne=0;ne<nE;ne++){ if (fabs(wf[ne]-min)<wtol){ nmin += 1; if (minindex!=NULL){ *minindex=ne;}}}
  if (fabs(wf[nS]-min)<wtol){ wr=(nmin==1);} else{ wr=0;}
  tfree(wf); wf=NULL;
  if (minindex!=NULL && nmin>1){ *minindex=-1;}
  return wr;
}

double ptreex1_smatch(int nE,int nS,double *we,double wa,int *minindex)
{
  /* additive matching criterion for nE expected weights *we and actual weight wa, given that the actual index should be nS 
     if wa matches nmin of the maxima of *we, we pass back 1.0/nmin */
  int ne=0,nmin=0;
  double *wf=NULL,min=0,wtol=0.0000001,wr=0;
  wf = (double *) tcalloc(nE,sizeof(double));
  for (ne=0;ne<nE;ne++){ wf[ne] = fabs(we[ne]-wa);}
  stats("double",wf,nE,NULL,&min,NULL,NULL);
  nmin=0; for (ne=0;ne<nE;ne++){ if (fabs(wf[ne]-min)<wtol){ nmin += 1; if (minindex!=NULL){ *minindex=ne;}}}
  if (fabs(wf[nS]-min)<wtol){ wr=1.0/(double)nmin;} else{ wr=0;}
  tfree(wf); wf=NULL;
  if (minindex!=NULL && nmin>1){ *minindex=-1;}
  return wr;
}

double ptreex1_xmatch(int nE,int nS,double *we,double wa,int *maxindex)
{
  /* multiplicative matching criterion for nE expected weights *we and actual weight wa, given that the actual index should be nS */
  int ne=0,nmax=0;
  double *wf=NULL,max=0,wtol=0.0000001,wr=0;
  wf = (double *) tcalloc(nE,sizeof(double));
  for (ne=0;ne<nE;ne++){ wf[ne] = ((we[ne]!=0 && wa!=0) ? minimum(we[ne]/wa,wa/we[ne]) : 0);}
  stats("double",wf,nE,&max,NULL,NULL,NULL);
  nmax=0; for (ne=0;ne<nE;ne++){ if (fabs(wf[ne]-max)<wtol){ nmax += 1; if (maxindex!=NULL){ *maxindex=ne;}}}
  if (fabs(wf[nS]-max)<wtol){ wr=1.0/(double)nmax;} else{ wr=0;}
  tfree(wf); wf=NULL;
  if (maxindex!=NULL && nmax>1){ *maxindex=-1;}
  return wr;
}

void ptreex1(int nE,struct ptree **pE,struct pnode **parentE,struct llitem **lE,int nS,struct ptree *pS,struct pnode *parentS,struct llitem *lS,struct ptree *pR,struct pnode *parentR,struct llitem *lR,int maxlevel,int level)
{
  /* computes individual x1 estimate of sample ptree (*pS) against array of nE expected ptrees (**pE),
     output stored in record ptree (*pR), 
     for which it is assumed that the structure of *pR is the union of all the **pE (or more precisely, all possible *pS) */
  int verbose=0;
  char lvlchar[64];
  int firstentry=0,ne=0;
  struct pnode *pnR=NULL,*pnS=NULL,**pnE=NULL;
  struct llitem *lS2=NULL,**lE2=NULL;
  double *weight_expected=NULL,weight_actual=0;
  if (verbose){ sprintf(lvlchar," %%"); for (ne=0;ne<level;ne++){ sprintf(lvlchar,"%s %%",lvlchar);}}
  if (verbose>1){ printf(" %s \n",lvlchar);}
  if (verbose>1){ 
    printf(" %s [entering ptreex1] with nE %d and nS %d\n",lvlchar,nE,nS);
    for (ne=0;ne<nE;ne++){ printf(" %s pE[%d] %d parentE[%d] %d->%d lE[%d] %d \n",lvlchar,ne,(int)pE[ne],ne,(int)parentE[ne],parentE[ne]==NULL ? 0 : (int)(parentE[ne]->childllitem),ne,(int)lE[ne]);}
    printf(" %s pS %d parentS %d->%d lS %d\n",lvlchar,(int)pS,(int)parentS,parentS==NULL ? 0 : (int)(parentS->childllitem),(int)lS);
    printf(" %s pR %d parentR %d->%d lR %d\n",lvlchar,(int)pR,(int)parentR,parentR==NULL ? 0 : (int)(parentR->childllitem),(int)lR);
    printf(" %s maxlevel %d level %d\n",lvlchar,maxlevel,level);}
  if (parentS!=NULL){ assert(parentS->childllitem==lS);}
  for (ne=0;ne<nE;ne++){ if (parentE[ne]!=NULL){ assert(parentE[ne]->childllitem==lE[ne]);}}
  if (parentR!=NULL && lR->parent==NULL){ /* first descent into parentR->childllitem */
    assert(parentR->childllitem==lR);
    if (verbose>1){ printf(" %s first descent into parentR->childllitem, setting firstentry to 1\n",lvlchar);}
    firstentry=1;}
  if (lR->kidl!=NULL){ 
    if (verbose>1){ printf(" %s stepping left to lR->kidl %d\n",lvlchar,(int)(lR->kidl));}
    ptreex1(nE,pE,parentE,lE,nS,pS,parentS,lS,pR,parentR,lR->kidl,maxlevel,level);}
  if (lR->kidr!=NULL){
    if (verbose>1){ printf(" %s stepping right to lR->kidr %d\n",lvlchar,(int)(lR->kidr));}
    ptreex1(nE,pE,parentE,lE,nS,pS,parentS,lS,pR,parentR,lR->kidr,maxlevel,level);}
  if (lR->item!=NULL){
    if (verbose>1){ printf(" %s lR->item %d exists, delving deeper\n",lvlchar,(int)(lR->item));}
    pnR = (struct pnode *)lR->item;
    if (level < maxlevel){
      if (verbose>1){ printf(" %s level %d less than maxlevel %d, continuing...\n",lvlchar,level,maxlevel);}
      if (lS!=NULL && (lS2 = llitemaddorfind(0,lS,pS->regionra[pnR->region->label],&region2pnode_compare_label))!=NULL){
	pnS = (struct pnode *)lS2->item;}
      else{ pnS=NULL;}
      lS2 = (pnS==NULL ? NULL : pnS->childllitem);
      if (verbose>1){ 
	if (pnS!=NULL){ printf(" %s pnS label %d matches pnR label %d, ...\n",lvlchar,pnS->region->label,pnR->region->label);}
	else{ printf(" %s pnS label doesn't match pnR label %d, ...\n",lvlchar,pnR->region->label);}}
      pnE = (struct pnode **) tcalloc(nE,sizeof(struct pnode *));
      lE2 = (struct llitem **) tcalloc(nE,sizeof(struct llitem *));
      for (ne=0;ne<nE;ne++){
	if (lE[ne]!=NULL && (lE2[ne] = llitemaddorfind(0,lE[ne],pE[ne]->regionra[pnR->region->label],&region2pnode_compare_label))!=NULL){
	  pnE[ne] = (struct pnode *)lE2[ne]->item;}
	else{ pnE[ne] = NULL;}
	lE2[ne] = (pnE[ne]==NULL ? NULL : pnE[ne]->childllitem);
	if (verbose>1){ 
	  if (pnE[ne]!=NULL){ 
	    printf(" %s pnE[%d] label %d matches pnR label %d, ...\n",lvlchar,ne,pnE[ne]->region->label,pnR->region->label);}
	  else{ printf(" %s pnE[%d] label doesn't match pnR label %d, ...\n",lvlchar,ne,pnR->region->label);}}}
      if (verbose>1){ printf(" %s descending\n",lvlchar);}
      ptreex1(nE,pE,pnE,lE2,nS,pS,pnS,lS2,pR,pnR,pnR->childllitem,maxlevel,level+1);
      tfree(pnE); pnE=NULL;
      tfree(lE2); lE2=NULL;}}
  if (firstentry){
    if (verbose>1){ printf(" %s first entry for actual index %d,\n",lvlchar,nS);}
    weight_expected = (double *) tcalloc(nE,sizeof(double));
    for (ne=0;ne<nE;ne++){ if (parentE[ne]!=NULL){ weight_expected[ne] = parentE[ne]->weight*pS->total_time/pE[ne]->total_time;}}
    if (verbose>1){ printf(" %s index expcted weights:",lvlchar); for (ne=0;ne<nE;ne++){ printf(" %f",weight_expected[ne]);} printf("\n");}
    if (parentS!=NULL){ weight_actual = parentS->weight;} else{ weight_actual=0;}
    if (verbose>1){ printf(" %s index actual weight:  %f\n",lvlchar,weight_actual);}
    if (verbose>1){ printf(" %s computing best index\n",lvlchar);}
    parentR->weight += ptreex1_smatch(nE,nS,weight_expected,weight_actual,NULL);
    parentR->relevance += ptreex1_xmatch(nE,nS,weight_expected,weight_actual,NULL);
    if (verbose>1){ printf(" %s indexing... weight set to %f and relevance set to %f\n",lvlchar,parentR->weight,parentR->relevance);}
    tfree(weight_expected);weight_expected=NULL;}
}

void ptreex1sbias(int nE,struct ptree **pE,struct pnode **parentE,struct llitem **lE,int nS,struct ptree *pS,struct pnode *parentS,struct llitem *lS,struct ptree *pR,struct pnode *parentR,struct llitem *lR,int maxlevel,int level,int nX,double *percentages,double *sbias)
{
  /* calculates the biased sum-based x1 index estimator for sample *pS with actual index nS
     in comparison to nE expected ptrees **pE
     using reference ptree *pR which has been normalized to have entries between 0 and 1 (probabilities).
     Assumes that *percentages is an ascending array of nX elements, 
     and *sbias is some pregenerated array of nbases*nX*maxelevel elements
  */
  int verbose=0;
  char lvlchar[64];
  int firstentry=0,ne=0,nx=0,nl=0;
  struct pnode *pnR=NULL,*pnS=NULL,**pnE=NULL;
  struct llitem *lS2=NULL,**lE2=NULL;
  double *weight_expected=NULL,weight_actual=0;
  int minindex=0;
  double nadd=0;
  if (verbose){ sprintf(lvlchar," %%"); for (ne=0;ne<level;ne++){ sprintf(lvlchar,"%s %%",lvlchar);}}
  if (verbose>1){ printf(" %s \n",lvlchar);}
  if (verbose>1){ 
    printf(" %s [entering ptreex1sbias] with nE %d and nS %d and nX %d\n",lvlchar,nE,nS,nX);
    for (ne=0;ne<nE;ne++){ printf(" %s pE[%d] %d parentE[%d] %d->%d lE[%d] %d \n",lvlchar,ne,(int)pE[ne],ne,(int)parentE[ne],parentE[ne]==NULL ? 0 : (int)(parentE[ne]->childllitem),ne,(int)lE[ne]);}
    printf(" %s pS %d parentS %d->%d lS %d\n",lvlchar,(int)pS,(int)parentS,parentS==NULL ? 0 : (int)(parentS->childllitem),(int)lS);
    printf(" %s pR %d parentR %d->%d lR %d\n",lvlchar,(int)pR,(int)parentR,parentR==NULL ? 0 : (int)(parentR->childllitem),(int)lR);
    printf(" %s maxlevel %d level %d\n",lvlchar,maxlevel,level);
    printf(" %s percentages:",lvlchar); for (nx=0;nx<nX;nx++){ printf(" %0.3f",percentages[nx]);} printf("\n");}
  if (parentS!=NULL){ assert(parentS->childllitem==lS);}
  for (ne=0;ne<nE;ne++){ if (parentE[ne]!=NULL){ assert(parentE[ne]->childllitem==lE[ne]);}}
  if (parentR!=NULL && lR->parent==NULL){ /* first descent into parentR->childllitem */
    assert(parentR->childllitem==lR);
    if (verbose>1){ printf(" %s first descent into parentR->childllitem, setting firstentry to 1\n",lvlchar);}
    firstentry=1;}
  if (lR->kidl!=NULL){ 
    if (verbose>1){ printf(" %s stepping left to lR->kidl %d\n",lvlchar,(int)(lR->kidl));}
    ptreex1sbias(nE,pE,parentE,lE,nS,pS,parentS,lS,pR,parentR,lR->kidl,maxlevel,level,nX,percentages,sbias);}
  if (lR->kidr!=NULL){
    if (verbose>1){ printf(" %s stepping right to lR->kidr %d\n",lvlchar,(int)(lR->kidr));}
    ptreex1sbias(nE,pE,parentE,lE,nS,pS,parentS,lS,pR,parentR,lR->kidr,maxlevel,level,nX,percentages,sbias);}
  if (lR->item!=NULL){
    if (verbose>1){ printf(" %s lR->item %d exists, delving deeper\n",lvlchar,(int)(lR->item));}
    pnR = (struct pnode *)lR->item;
    if (level < maxlevel){
      if (verbose>1){ printf(" %s level %d less than maxlevel %d, continuing...\n",lvlchar,level,maxlevel);}
      if (lS!=NULL && (lS2 = llitemaddorfind(0,lS,pS->regionra[pnR->region->label],&region2pnode_compare_label))!=NULL){
	pnS = (struct pnode *)lS2->item;}
      else{ pnS=NULL;}
      lS2 = (pnS==NULL ? NULL : pnS->childllitem);
      if (verbose>1){ 
	if (pnS!=NULL){ printf(" %s pnS label %d matches pnR label %d, ...\n",lvlchar,pnS->region->label,pnR->region->label);}
	else{ printf(" %s pnS label doesn't match pnR label %d, ...\n",lvlchar,pnR->region->label);}}
      pnE = (struct pnode **) tcalloc(nE,sizeof(struct pnode *));
      lE2 = (struct llitem **) tcalloc(nE,sizeof(struct llitem *));
      for (ne=0;ne<nE;ne++){
	if (lE[ne]!=NULL && (lE2[ne] = llitemaddorfind(0,lE[ne],pE[ne]->regionra[pnR->region->label],&region2pnode_compare_label))!=NULL){
	  pnE[ne] = (struct pnode *)lE2[ne]->item;}
	else{ pnE[ne] = NULL;}
	lE2[ne] = (pnE[ne]==NULL ? NULL : pnE[ne]->childllitem);
	if (verbose>1){ 
	  if (pnE[ne]!=NULL){ 
	    printf(" %s pnE[%d] label %d matches pnR label %d, ...\n",lvlchar,ne,pnE[ne]->region->label,pnR->region->label);}
	  else{ printf(" %s pnE[%d] label doesn't match pnR label %d, ...\n",lvlchar,ne,pnR->region->label);}}}
      if (verbose>1){ printf(" %s descending\n",lvlchar);}
      ptreex1sbias(nE,pE,pnE,lE2,nS,pS,pnS,lS2,pR,pnR,pnR->childllitem,maxlevel,level+1,nX,percentages,sbias);
      tfree(pnE);
      tfree(lE2);}}
  if (firstentry){
    if (verbose>1){ printf(" %s first entry at level %d, actual index %d, element's value is %f\n",lvlchar,level,nS,parentR->weight);}
    if (parentR->weight>=percentages[0]){ /* not bad */
      weight_expected = (double *) tcalloc(nE,sizeof(double));
      for (ne=0;ne<nE;ne++){ if (parentE[ne]!=NULL){ weight_expected[ne] = parentE[ne]->weight*pS->total_time/pE[ne]->total_time;}}
      if (parentS!=NULL){ weight_actual = parentS->weight;} else{ weight_actual=0;}
      nadd = ptreex1_smatch(nE,nS,weight_expected,weight_actual,&minindex);
      if (verbose>1){ printf(" %s expected weights:",lvlchar); for (ne=0;ne<nE;ne++){ printf(" %0.5f",weight_expected[ne]);} printf("\n"); printf(" %s actual weight %0.2f\n",lvlchar,weight_actual); printf(" %s determined minindex %d\n",lvlchar,minindex);}
      /* if we want to only consider unique minimums */
      if (minindex!=-1){ 
        if (verbose>1){ printf(" %s unique minimum found at index %d, with reference %0.4f, sbias starting at:\n",lvlchar,minindex,parentR->weight); printf(" %s\n",lvlchar); for (nl=0;nl<maxlevel;nl++){ printf(" %s level %d\n",lvlchar,nl); raprintf(&(sbias[0+0*nE+nl*nX*nE]),"double",nE,nX,lvlchar);}}
	for (nx=0;nx<nX;nx++){ if (parentR->weight>=percentages[nx]){ for (nl=level-1;nl<maxlevel;nl++){ sbias[minindex+nx*nE+nl*nE*nX] += 1;}}}
        if (verbose>1){ printf(" %s now sbias:\n",lvlchar); printf(" %s\n",lvlchar); for (nl=0;nl<maxlevel;nl++){ printf(" %s level %d\n",lvlchar,nl); raprintf(&(sbias[0+0*nE+nl*nX*nE]),"double",nE,nX,lvlchar);}}}
      /* on the other hand, if nonstrict minimums won't affect the answer if nbases==2 */
      tfree(weight_expected);}}
}

void pnodeprintf_label(void *v)
{
  /* prints pn->region->label */
  struct pnode *pn=NULL;
  if (v!=NULL){ pn=(struct pnode *)v; printf("pnode %d label %d\n",(int)pn,pn->region->label);}
  else{ printf("NULL\n");}
}

void pnodeshizuffle_starter(struct pnode *parent,struct llitem *l0,struct ptree *p)
{
  int verbose=0;
  int firstentry=0;
  struct pnode *pn=NULL;
  struct llist *L=NULL;
  struct litem *l=NULL;
  struct pnode **permutation=NULL,*pntemp=NULL;
  int npasses=5;
  int nr=0,pr=0,index1=0,index2=0;
  if (verbose){ printf(" %% \n");}
  if (verbose){ printf(" %% [entering pnodeshizuffle_starter]\n");}
  if (parent!=NULL && l0->parent==NULL){ /* first descent into parent->childllitem */ 
    assert(parent->childllitem==l0); firstentry=1;}
  if (l0->kidl!=NULL){ 
    if (verbose){ printf(" %% passing left\n");} pnodeshizuffle_starter(parent,l0->kidl,p);}
  if (l0->kidr!=NULL){ 
    if (verbose){ printf(" %% passing right\n");} pnodeshizuffle_starter(parent,l0->kidr,p);}
  if (l0->item!=NULL){ 
    if (verbose){ printf(" %% passing down\n");} pn=(struct pnode *)l0->item; pnodeshizuffle_starter(pn,pn->childllitem,p);}
  if (firstentry){ 
    if (parent!=NULL && l0!=NULL && l0->item!=NULL){ 
      pr = parent->region->label;
      if (verbose){ printf(" %% firstentry, starting with parent region %d, and childllitem:\n",pr); llitemprintf(l0,&pnodeprintf_label);}
      L=llistmake();
      llitem2llist(l0,L); 
      llistsort(L->first,L->last,L->length,&pnode2pnode_compare_label); /* should be unnecessary */
      if (verbose){ printf(" %% constructed llist:\n"); l=L->first;while(l!=NULL){ pn=(struct pnode *)l->item;printf("%d\n",pn->region->label);l=l->child;}}
      permutation = (struct pnode **) tcalloc(p->nregions,sizeof(struct pnode *));
      nr=0;l=L->first;
      while(l!=NULL && nr<p->nregions){ 
	pn=(struct pnode *)l->item; 
	if (pn->region->label==nr){ permutation[nr]=pn; nr+=1; l=l->child;}
	else{ permutation[nr]=NULL; nr+=1;}}
      llisttfree(L);L=NULL;
      if (verbose){ printf(" %% constructed permutation:\n"); for (nr=0;nr<p->nregions;nr++){ if (permutation[nr]!=NULL){ printf("pnode %d label %d\n",(int)permutation[nr],permutation[nr]->region->label);} else{ printf("NULL\n");}}}
      for (nr=0;nr<npasses*p->nregions;nr++){
	index1=periodize((int)floor(p->nregions*rand01),0,p->nregions);
	index2=periodize((int)floor(p->nregions*rand01),0,p->nregions);
	if (index1!=pr && index2!=pr){
	  pntemp = permutation[index1];
	  permutation[index1] = permutation[index2];
	  permutation[index2] = pntemp;}}
      for (nr=0;nr<p->nregions;nr++){ if (permutation[nr]!=NULL){ permutation[nr]->region = p->regionra[nr];}}
      if (verbose){ printf(" %% transposed into permutation:\n"); for (nr=0;nr<p->nregions;nr++){ if (permutation[nr]!=NULL){ printf("pnode %d label %d\n",(int)permutation[nr],permutation[nr]->region->label);} else{ printf("NULL\n");}}}
      L=llistmake();
      for (nr=0;nr<p->nregions;nr++){ if (permutation[nr]!=NULL){ litemadd(L,permutation[nr]);}}
      if (verbose){ printf(" %% transposed into llist:\n"); l=L->first;while(l!=NULL){ pn=(struct pnode *)l->item;printf("%d\n",pn->region->label);l=l->child;}}
      llitemtfree(l0,NULL);l0=llitemmake();
      llist2llitem(L,l0,&pnode2pnode_compare_label);
      llisttfree(L);L=NULL;
      tfree(permutation);permutation=NULL;
      parent->childllitem=l0;
      llitembalance(parent->childllitem);parent->childllitem=llitemclimb(parent->childllitem);l0=parent->childllitem;
      if (verbose){ printf(" %% firstentry, ending with parent region %d, and childllitem:\n",pr); llitemprintf(parent->childllitem,&pnodeprintf_label);}}}
}

void ptreex2_disthist(double *distra,double *x2ra,int nbases)
{
  /* adds (1.0/(double)nmin) to the indices of distra which correspond to the minima of x2ra */
  double max=0,min=0,wtol=0.0000001;
  int nr=0,nmin=0;
  stats("double",x2ra,nbases,&max,&min,NULL,NULL);
  for (nr=0;nr<nbases;nr++){ if (fabs(x2ra[nr]-min)<wtol){ nmin += 1;}} 
  nmin=maximum(1,nmin);
  for (nr=0;nr<nbases;nr++){ if (fabs(x2ra[nr]-min)<wtol){ distra[nr] += (1.0/(double)nmin);}}
}

double logfactorial(double n)
{
/*   return 0.5*(log(2.0*n + 1.0/3.0)+log(PI)) + n*log(n>0 ? n : 1) - n; */
  return lgamma(n+1);
}

double binomial(double n,double N,double p)
{
  /* Probability that given N choices of 2 options (option 1 has probability p), option 1 occurs n times */
  /* approximate binomial distribution is used */
  int verbose=0;
  double n2=0;
  double wtol=0.0000001;
  double b=0;
  if (verbose){ printf(" %% binomial(%f,%f,%f)",n,N,p);}
  N=maximum(0,N);n=maximum(0,n);p=minimum(1,maximum(0,p));
  n2=N-n;
  if (n>=0 && n2>=0){
    if (fabs(N-n)<wtol){ b = exp(n*log(p));}
    else if (fabs(N-n2)<wtol){ b = exp(n2*log(1-p));}
    else /* if (n<N && n2<N) */{ b = exp(logfactorial(N) - logfactorial(n) - logfactorial(n2) + n*log(p) + n2*log(1-p));}}
  else /* if (n<0 || n2<0) */{ b = 0;}
  if (verbose){ printf("=%f\n",b);}
  return b;
}

double Pmode(double n,double N,int m,double *p)
{
  /* Probability that given N choices of m options (option i has probability p[i]), the mode occurs more than n times */
  /* GLOBAL_QUADRATURE_N, GLOBAL_QUADRATURE_Z, GLOBAL_QUADRATURE_W are used as quadrature nodes */
  /* for integrals spanning more than quadthresh integers we approximate
     \int_{0}^{n}d\theta binomial(\theta,N,p[0])*Pmode(n,N-\theta,m-1,&(p[1]))
     + \int_{n}^{N}d\theta binomial(\theta,N,p[0]) */
  /* my intuition says that p should be sorted in descending order */
  int verbose=0;
  int quadthresh=100;
  int nl=0;
  double sumprob=0;
  double icenter=0,iradius=0,ix=0;
  double integral1=0,integral2=0;
  if (verbose){ printf(" %% [entering Pmode]\n %% n=%f,N=%f,m=%d",n,N,m); raprintf(p,"double",1,m," p = ");}
  if (n>=N){ if (verbose){ printf(" %% n>=N, exiting\n");} return 0;}
  else /* if (n<N) */{
    if (m<1){ printf(" %% error, m=%d in Pmode\n",m);}
    else if (m==1){ if (verbose){ printf(" %% only one choice\n");} return 1;}
    else if (m>1){ 
      if (n<(N/(double)m)){ 
	if (verbose){ printf(" %% too few choices, mode at least N/m=%f\n",N/(double)m);}
	return 1;}
      else /* if n>=(N/(double)m) */{
	if (verbose){ printf(" %% many choices, mode may not be %f or more\n",n+1);}
	stats("double",p,m,NULL,NULL,&sumprob,NULL); sumprob *= m; if (sumprob<=0){ sumprob=1;}
	if (verbose){ printf(" %% probability sums to %f\n",sumprob);}
	if (n>quadthresh){
	  integral1=0;icenter=(n+0.5+0.0)/2.0;iradius=(n+0.5-0.0)/2.0;
	  for (nl=0;nl<GLOBAL_QUADRATURE_N;nl++){
	    ix = GLOBAL_QUADRATURE_Z[nl]*iradius + icenter;
	    integral1 += GLOBAL_QUADRATURE_W[nl]*binomial(ix,N,p[0]/sumprob)*Pmode(n,N-ix,m-1,&(p[1]));}
	  integral1 *= iradius;}
	else /* if domain is small */{
	  integral1=0;
	  for (nl=0;nl<=n;nl++){
	    integral1 += binomial(nl,N,p[0]/sumprob)*Pmode(n,N-nl,m-1,&(p[1]));}}
	if ((N-n-1)>quadthresh){
	  integral2=0;icenter=(N+n+0.5)/2.0;iradius=(N-n-0.5)/2.0;
	  for (nl=0;nl<GLOBAL_QUADRATURE_N;nl++){
	    ix = GLOBAL_QUADRATURE_Z[nl]*iradius + icenter;
	    integral2 += GLOBAL_QUADRATURE_W[nl]*binomial(ix,N,p[0]/sumprob);}
	  integral2 *= iradius;}
	else /* if domain is small */{ 
	  integral2=0;
	  for (nl=n+1;nl<=N;nl++){
	    integral2 += binomial(nl,N,p[0]/sumprob);}}
	if (verbose){ printf(" %% returning %f + %f\n",integral1,integral2);}
	return integral1+integral2;}}}
  return 0;
}

double multithresh_single(double *pra,int nbases,double threshold)
{
  /* returns minimum number of independent observations which guarantees (up to threshold) 
     that the most probable observation (with probability p[0]) is the mode */  
  /* GLOBAL_QUADRATURE_N, GLOBAL_QUADRATURE_Z, GLOBAL_QUADRATURE_W are used as quadrature nodes */
  /* calculate 
     \int_{N/m}^{N}d\theta binomial(\theta,N,p[0])*(1-Pmode(\theta,N-\theta,m-1,&(p[1]))) */
  int quadthresh=100;
  int nb=0,nl=0;
  struct llist *L=NULL;
  struct litem *l=NULL;
  double *pra2=NULL,sumprob=0;
  double icenter=0,iradius=0,ix=0;
  double integral=0;
  int nobs=0,nobsmax=100;
  L=llistmake();
  for (nb=0;nb<nbases;nb++){ litemadd(L,&(pra[nb]));}
  llistsort(L->first,L->last,L->length,&double_compare);
  pra2 = (double *) tcalloc(nbases,sizeof(double));
  l=L->last; nb=0; while (l!=NULL){ pra2[nb]=*(double *)l->item; l=l->parent; nb++;}
  stats("double",pra2,nbases,NULL,NULL,&sumprob,NULL); sumprob *= nbases; if (sumprob<=0){ sumprob=1;}
  llisttfree(L);L=NULL;
  nobs=0;
  do{
    nobs += 1;
    if (nobs>quadthresh){
      integral=0;icenter=((double)nobs/(double)nbases + (double)nobs)/2.0;iradius=((double)nobs-(double)nobs/(double)nbases)/2.0;
      for (nl=0;nl<GLOBAL_QUADRATURE_N;nl++){
	ix = GLOBAL_QUADRATURE_Z[nl]*iradius + icenter;
	integral += GLOBAL_QUADRATURE_W[nl]*binomial(ix,(double)nobs,pra2[0]/sumprob)*(1-Pmode(ix,(double)nobs-ix,nbases-1,&(pra2[1])));}
      integral *= iradius;}
    else /* if domain is small */{
      integral=0;
      for (nl=(int)floor(nobs/(double)nbases);nl<=nobs;nl++){
	integral += binomial(nl,nobs,pra2[0]/sumprob)*(1-Pmode(nl,nobs-nl,nbases-1,&(pra2[1])));}}}
  while (nobs<nobsmax && integral<threshold);
  tfree(pra2);pra2=NULL;
  return nobs;
}

double multithresh(double *distra,int nbases,double threshold)
{
  /* probabilities are normalized versions of distra[observed_base + actual_base*nbases]
     returns minimum number of independent observations which guarantees (up to threshold) observing actual_base as the mode 
     or -1 if not possible */
  int verbose=0;
  double mean=0,*pra=NULL,max=0,wtol=0.0000001,temp=0;
  int nr=0,nr2=0;
  double nobs=0,nobsmax=100;
  if (verbose){ printf(" %% [entering multithresh]\n");}
  pra = (double *) tcalloc(nbases*nbases,sizeof(double));
  for (nr=0;nr<nbases;nr++){ 
    stats("double",&(distra[0+nr*nbases]),nbases,NULL,NULL,&mean,NULL); 
    for (nr2=0;nr2<nbases;nr2++){ pra[nr2+nr*nbases] = distra[nr2+nr*nbases]/(mean*(double)nbases);}}
  if (verbose){ raprintf(distra,"double",nbases,nbases," %% distra "); raprintf(pra,"double",nbases,nbases," %% pra ");}
  nobs=0;
  for (nr=0;nr<nbases;nr++){
    stats("double",&(pra[0+nr*nbases]),nbases,&max,NULL,NULL,NULL);
    if (fabs(pra[nr+nr*nbases]-max)>wtol){ 
      if (verbose){ printf(" %% actual_base %d not possible\n",nr);} return nobsmax;}
    else{ 
      temp=pra[0+nr*nbases]; pra[0+nr*nbases]=pra[nr+nr*nbases]; pra[nr+nr*nbases]=temp; 
      if (verbose){ printf(" %% shuffled pra for %d\n",nr); raprintf(&(pra[0+nr*nbases]),"double",1,nbases," %% pra ");}
      nobs += (multithresh_single(&(pra[0+nr*nbases]),nbases,threshold))/(double)nbases;}}
  tfree(pra);pra=NULL;
  if (verbose){ printf(" %% returning %f\n",nobs);}
  return nobs;
}

void llistmedian(struct llist *L1,struct llist *L2,int (*compare)(void *,void *),int sort_flag,void **iout,int *L1less,int *L2more)
{
  /* assumes that L1 and L2 are (disjoint and distinct) llists of comparable objects
     Here we return the address iout of the item (either within L1 or L2) which is the comparative median, 
     that is, which maximizes the number of elements of L1 that are < iout, and of L2 which are >= iout 
     these numbers are stored as *L1less and *L2more */
  int verbose=0;
  struct llist *L1t=NULL,*L2t=NULL,*L3t=NULL;
  struct litem *l1=NULL,*l2=NULL,*l3=NULL;
  int L3less=0,L3more=0,L3tot=0;
  if (verbose){ if (compare==&double_compare){ 
    printf(" %% double comparison of\n");llistprintf2(L1);printf(" %% and\n");llistprintf2(L2);}}
  if (sort_flag){ L1t=L1;L2t=L2;} else{ L1t=llistcopy(L1);L2t=llistcopy(L2);} L3t=llistmake();
  llistsort(L1t->first,L1t->last,L1t->length,compare);
  llistsort(L2t->first,L2t->last,L2t->length,compare);
  l1=L1t->first;l2=L2t->first;
  while (l1!=NULL || l2!=NULL){
    if (0){}
    else if (l1==NULL){ litemadd(L3t,l2->item); l2=l2->child;}
    else if (l2==NULL){ litemadd(L3t,l1->item); l1=l1->child;}
    else if (l1!=NULL && l2!=NULL){
      if ((*compare)(l1->item,l2->item)>0){ litemadd(L3t,l2->item); l2=l2->child;}
      else{ litemadd(L3t,l1->item); l1=l1->child;}}}
  if (verbose){ if (compare==&double_compare){ 
    printf(" %% sorted...\n");llistprintf2(L1t);llistprintf2(L2t);llistprintf2(L3t);}}
  if (L1t->length==0 && L2t->length==0){ /* do nothing */ *L1less=0;*L2more=0;*iout=NULL;}
  else if (L1t->length==0 && L2t->length>0){ /* easy */ *L1less=0;*L2more=L2t->length;*iout=(void *)L2t->first->item;}
  else if (L1t->length>0 && L2t->length==0){ /* easy */ *L1less=L1t->length;*L2more=0;*iout=(void *)L1t->last->item;}
  else if (L1t->length>0 && L2t->length>0){ 
    l1=L1t->first;l2=L2t->first;l3=L3t->first; *iout=l3->item;
    if ((*compare)((void *)l3->item,(void *)l1->item)==0){ *L1less=0;*L2more=L2->length;}
    else if ((*compare)((void *)l3->item,(void *)l2->item)==0){ *L1less=0;*L2more=L2->length;}
    L3less=*L1less;L3more=*L2more;L3tot=L3less+L3more;
    if (verbose){ if (compare==&double_compare){
      printf(" %% starting with l1->item %f l2->item %f l3->item %f\n",*(double *)l1->item,*(double *)l2->item,*(double *)l3->item);
      printf(" %% *iout %f L3less %d L3more %d L3tot %d\n",*(double *)*iout,L3less,L3more,L3tot);}}
    while (l3!=NULL){
      if (l1!=NULL && l3->item==l1->item){ 
	if (verbose){ printf(" %%  belongs to first llist \n");}
	l1=l1->child; l3=l3->child; L3less += 1;
	if (l3!=NULL && L3less+L3more > L3tot){ *iout=l3->item; *L1less=L3less;*L2more=L3more; L3tot=L3less+L3more;}}
      else if (l2!=NULL && l3->item==l2->item){ 
	if (verbose){ printf(" %%  belongs to second llist \n");}
	l2=l2->child; l3=l3->child; L3more -= 1;
	if (l3!=NULL && L3less+L3more > L3tot){ *iout=l3->item; *L1less=L3less;*L2more=L3more; L3tot=L3less+L3more;}}
      if (verbose){ if (compare==&double_compare){ 
	printf(" %% %% now l1->item %f l2->item %f l3->item %f\n",l1==NULL?0:*(double *)l1->item,l2==NULL?0:*(double *)l2->item,l3==NULL?0:*(double *)l3->item);
	printf(" %% %% *iout %f L3less %d L3more %d L3tot %d\n",*(double *)*iout,L3less,L3more,L3tot);}}}
    if (verbose){ if (compare==&double_compare){
      printf(" %% ending with *iout %f *L1less %d *L2more %d L3tot %d\n",*(double *)*iout,*L1less,*L2more,L3tot);}}}
  if (!sort_flag){ llisttfree(L1t);L1t=NULL;llisttfree(L2t);L2t=NULL;} llisttfree(L3t);L3t=NULL;
}

void discdata_threshold(double *discdata,int nbases,int maxlevel,int nrecords,char *outputname)
{
  /* given the indicators of discdata (sorted so that discdata[0 + level*nbases + nr*nbases*maxlevel] should be lowest) */
  int verbose=0;
  int nl=0,nb=0,nr=0;
  double val=0;
  struct llist *L=NULL,*Lsml1=NULL,*Lsml2=NULL,*Lbig1=NULL,*Lbig2=NULL;
  double wtol=0.0000001;
  int Lsmlless=0,Lsmlmore=0,Lbigless=0,Lbigmore=0;
  void *ismlout=NULL,*ibigout=NULL;
  double min1=0,min2=0,max1=0,max2=0;
  double *smlmarg=NULL,*smlpass=NULL,*bigmarg=NULL,*bigpass=NULL,*classifysml=NULL,*classifybig=NULL,*threshsml=NULL,*threshbig=NULL;
  double *outprint=NULL;
  if (maxlevel>1 && nbases>1){
    smlmarg = (double *) tcalloc(maxlevel*nrecords,sizeof(double));
    smlpass = (double *) tcalloc(maxlevel*nrecords,sizeof(double));
    bigmarg = (double *) tcalloc(maxlevel*nrecords,sizeof(double));
    bigpass = (double *) tcalloc(maxlevel*nrecords,sizeof(double));
    threshsml = (double *) tcalloc(maxlevel,sizeof(double));
    threshbig = (double *) tcalloc(maxlevel,sizeof(double));
    outprint = (double *) tcalloc(maxlevel,sizeof(double));
    for (nr=0;nr<nrecords;nr++){
      if (verbose){ 
	printf(" %% examining record %d\n",nr);
	raprintf(&(discdata[0 + 0*nbases + nr*nbases*maxlevel]),"double",nbases,maxlevel," %% ");}
      for (nl=0;nl<maxlevel;nl++){
	val = discdata[0+nl*nbases+nr*nbases*maxlevel];
	L = llistmake();
	for (nb=0;nb<nbases;nb++){ litemadd(L,&(discdata[nb+nl*nbases+nr*nbases*maxlevel]));}
	if (verbose){ printf(" %% level %d ",nl); llistprintf2(L);}
	llistsort(L->first,L->last,L->length,&double_compare);
	if (verbose){ printf(" %% sorted... "); llistprintf2(L);}
	min1=*(double *)L->first->item; min2=*(double *)L->first->child->item;
	max1=*(double *)L->last->item; max2=*(double *)L->last->parent->item;
	smlpass[nl+nr*maxlevel] = (fabs(val-min1)<wtol?+1:-1); smlmarg[nl+nr*maxlevel] = fabs((min2-min1)/min1);
	bigpass[nl+nr*maxlevel] = (fabs(val-max1)<wtol?+1:-1); bigmarg[nl+nr*maxlevel] = fabs((max2-max1)/max2);
	if (verbose){ printf(" %% min1 %f min2 %f max2 %f max1 %f, smlmarg %f smlpass %f, bigmarg %f bigpass %f\n",min1,min2,max2,max2,smlmarg[nl+nr*maxlevel],smlpass[nl+nr*maxlevel],bigmarg[nl+nr*maxlevel],bigpass[nl+nr*maxlevel]);}
	llisttfree(L);L=NULL;}}
    classifysml = (double *) tcalloc(nrecords*maxlevel,sizeof(double)); classifybig = (double *) tcalloc(nrecords*maxlevel,sizeof(double));
    for (nr=0;nr<nrecords;nr++){ 
      classifysml[nr+0*nrecords]=(smlpass[0+nr*maxlevel]>=0);classifybig[nr+0*nrecords]=(bigpass[0+nr*maxlevel]>=0);}
    for (nl=1;nl<maxlevel;nl++){
      if (verbose){ printf(" %% performing check for levels (%d,%d)\n",(nl-1),nl);}
      Lsml1=llistmake();Lsml2=llistmake(); Lbig1=llistmake();Lbig2=llistmake();
      for (nr=0;nr<nrecords;nr++){
	if (0){}
	else if (!classifysml[nr+(nl-1)*nrecords] && smlpass[nl+nr*maxlevel]>=0){ /* population 1 */ 
	  litemadd(Lsml1,&(smlmarg[(nl-1)+nr*maxlevel]));}
	else if (classifysml[nr+(nl-1)*nrecords] && smlpass[nl+nr*maxlevel]<0){ /* population 2 */ 
	  litemadd(Lsml2,&(smlmarg[(nl-1)+nr*maxlevel]));}
	if (0){}
	else if (!classifybig[nr+(nl-1)*nrecords] && bigpass[nl+nr*maxlevel]>=0){ /* population 1 */ 
	  litemadd(Lbig1,&(bigmarg[(nl-1)+nr*maxlevel]));}
	else if (classifybig[nr+(nl-1)*nrecords] && bigpass[nl+nr*maxlevel]<0){ /* population 2 */ 
	  litemadd(Lbig2,&(bigmarg[(nl-1)+nr*maxlevel]));}}
      if (verbose){ 
	printf(" %% sml bad to good "); llistprintf2(Lsml1);
	printf(" %% sml good to bad "); llistprintf2(Lsml2);
	printf(" %% big bad to good "); llistprintf2(Lbig1);
	printf(" %% big good to bad "); llistprintf2(Lbig2);}
      llistmedian(Lsml1,Lsml2,&double_compare,1,&ismlout,&Lsmlless,&Lsmlmore); threshsml[(nl-1)]=(ismlout==NULL?0:*(double *)ismlout);
      if (verbose){ printf(" %% sml divider %f, left %d, right %d\n",threshsml[(nl-1)],Lsmlless,Lsmlmore);}
      llistmedian(Lbig1,Lbig2,&double_compare,1,&ibigout,&Lbigless,&Lbigmore); threshbig[(nl-1)]=(ibigout==NULL?0:*(double *)ibigout);
      if (verbose){ printf(" %% big divider %f, left %d, right %d\n",threshsml[(nl-1)],Lbigless,Lbigmore);}
      for (nr=0;nr<nrecords;nr++){
	if (verbose){ printf(" %% classifysml[%d]=%f,smlpass[%d]=%f,smlmarg[%d]=%f,smlpass[%d]=%f,smlmarg[%d]=%f\n",nr+(nl-1)*nrecords,classifysml[nr+(nl-1)*nrecords],(nl-1)+nr*maxlevel,smlpass[(nl-1)+nr*maxlevel],(nl-1)+nr*maxlevel,smlmarg[(nl-1)+nr*maxlevel],nl+nr*maxlevel,smlpass[nl+nr*maxlevel],nl+nr*maxlevel,smlmarg[nl+nr*maxlevel]);}
	if (classifysml[nr+(nl-1)*nrecords]){
	  if (smlmarg[(nl-1)+nr*maxlevel]>=threshsml[(nl-1)]){ classifysml[nr+nl*nrecords] = 1;}
	  else /* if (smlmarg[(nl-1)+nr*maxlevel]<threshsml[(nl-1)]) */{ 
	    if (smlpass[nl+nr*maxlevel]>=0){ classifysml[nr+nl*nrecords] = 1;}
	    else /* if (smlpass[nl+nr*maxlevel]<0) */{ classifysml[nr+nl*nrecords] = 0;}}}
	else /* if (!classifysml[nr+(nl-1)*nrecords]) */{
	  if (smlmarg[(nl-1)+nr*maxlevel]>=threshsml[(nl-1)]){ classifysml[nr+nl*nrecords] = 0;}
	  else /* if (smlmarg[(nl-1)+nr*maxlevel]<threshsml[(nl-1)]) */{ 
	    if (smlpass[nl+nr*maxlevel]>=0){ classifysml[nr+nl*nrecords] = 1;}
	    else /* if (smlpass[nl+nr*maxlevel]<0) */{ classifysml[nr+nl*nrecords] = 0;}}}
	if (verbose){ printf(" %% classifysml[%d]=%f\n",nr+nl*nrecords,classifysml[nr+nl*nrecords]);}}
      for (nr=0;nr<nrecords;nr++){
	if (verbose){ printf(" %% classifybig[%d]=%f,bigpass[%d]=%f,bigmarg[%d]=%f,bigpass[%d]=%f,bigmarg[%d]=%f\n",nr+(nl-1)*nrecords,classifybig[nr+(nl-1)*nrecords],(nl-1)+nr*maxlevel,bigpass[(nl-1)+nr*maxlevel],(nl-1)+nr*maxlevel,bigmarg[(nl-1)+nr*maxlevel],nl+nr*maxlevel,bigpass[nl+nr*maxlevel],nl+nr*maxlevel,bigmarg[nl+nr*maxlevel]);}
	if (classifybig[nr+(nl-1)*nrecords]){
	  if (bigmarg[(nl-1)+nr*maxlevel]>=threshbig[(nl-1)]){ classifybig[nr+nl*nrecords] = 1;}
	  else /* if (bigmarg[(nl-1)+nr*maxlevel]<threshbig[(nl-1)]) */{ 
	    if (bigpass[nl+nr*maxlevel]>=0){ classifybig[nr+nl*nrecords] = 1;}
	    else /* if (bigpass[nl+nr*maxlevel]<0) */{ classifybig[nr+nl*nrecords] = 0;}}}
	else /* if (!classifybig[nr+(nl-1)*nrecords]) */{
	  if (bigmarg[(nl-1)+nr*maxlevel]>=threshbig[(nl-1)]){ classifybig[nr+nl*nrecords] = 0;}
	  else /* if (bigmarg[(nl-1)+nr*maxlevel]<threshbig[(nl-1)]) */{ 
	    if (bigpass[nl+nr*maxlevel]>=0){ classifybig[nr+nl*nrecords] = 1;}
	    else /* if (bigpass[nl+nr*maxlevel]<0) */{ classifybig[nr+nl*nrecords] = 0;}}}
	if (verbose){ printf(" %% classifybig[%d]=%f\n",nr+nl*nrecords,classifybig[nr+nl*nrecords]);}}
      llisttfree(Lsml1);Lsml1=NULL; llisttfree(Lsml2);Lsml2=NULL;
      llisttfree(Lbig1);Lbig1=NULL; llisttfree(Lbig2);Lbig2=NULL;}
    for (nl=0;nl<maxlevel;nl++){
      stats("double",&(classifysml[0+nl*nrecords]),nrecords,NULL,NULL,&max1,NULL);
      stats("double",&(classifybig[0+nl*nrecords]),nrecords,NULL,NULL,&max2,NULL);
      if (verbose){ printf(" %% level %d, smldisc %f, bigdisc %f\n",nl,max1,max2);}
      outprint[nl]=maximum(max1,max2);}
    ra2jpg2(outprint,"double",1,maxlevel,0,1.0,0.0,outputname);
    tfree(smlmarg);tfree(smlpass);tfree(classifysml);tfree(threshsml);
    tfree(bigmarg);tfree(bigpass);tfree(classifybig);tfree(threshbig);
    tfree(outprint);}
}

void ptreeobsdist(int nbases,struct ptree *pR,struct pnode *parentR,struct llitem *lR,int nb,struct ptree *p1,struct pnode *parent1,struct llitem *l1,int maxlevel,int level)
{
  /* creates llistra (nbases elements) at each node, corresponding to observation distribution of each event-chain */
  int verbose=1,nr=0,nb2=0;
  char lvlchar[64];
  int firstentry=0;
  struct pnode *pnR=NULL,*pn1=NULL;
  struct llitem *l2=NULL;
  struct llist **Lra=NULL;
  double wtol=0.0000001;
  double *temp=NULL;
  if (verbose){ sprintf(lvlchar," %%"); for (nr=0;nr<level;nr++){sprintf(lvlchar,"%s %%",lvlchar);}}
  if (verbose>1){ 
    printf(" %s [entering ptreeobsdist]\n",lvlchar);
    printf(" %s pR %d parentR %d->%d lR %d \n",lvlchar,(int)pR,(int)parentR,parentR==NULL ? 0 : (int)(parentR->childllitem),(int)lR);
    printf(" %s p1 %d parent1 %d->%d l1 %d\n",lvlchar,(int)p1,(int)parent1,parent1==NULL ? 0 : (int)(parent1->childllitem),(int)l1);
    printf(" %s nbases %d maxlevel %d level %d\n",lvlchar,nbases,maxlevel,level);}
  if (parentR!=NULL){ assert(parentR->childllitem==lR);}
  if (parent1!=NULL && l1->parent==NULL){ /* first descent into parent1->childllitem */ 
    if (verbose>1){ printf(" %s first descent into parent1->childllitem, setting firstentry to 1\n",lvlchar);}
    assert(parent1->childllitem==l1); firstentry=1;}
  if (l1->kidl!=NULL){ 
    if (verbose>1){ printf(" %s l1->kidl found, moving left to new llitem l1->kidl=%d\n",lvlchar,(int)(l1->kidl));}
    ptreeobsdist(nbases,pR,parentR,lR,nb,p1,parent1,l1->kidl,maxlevel,level);}
  if (l1->kidr!=NULL){ 
    if (verbose>1){ printf(" %s l1->kidr found, moving right to new llitem l1->kidr=%d\n",lvlchar,(int)(l1->kidr));}
    ptreeobsdist(nbases,pR,parentR,lR,nb,p1,parent1,l1->kidr,maxlevel,level);}
  if (l1->item!=NULL){ 
    if (verbose>1){ printf(" %s l1->item %d exists\n",lvlchar,(int)(struct pnode *)(l1->item));}
    pn1 = (struct pnode *)l1->item; 
    if (maxlevel==-1 || level<maxlevel){
      if (verbose>1){ printf(" %s level %d less than maxlevel %d, continuing...\n",lvlchar,level,maxlevel);}
      if (lR!=NULL && (l2=llitemaddorfind(0,lR,pR->regionra[pn1->region->label],&region2pnode_compare_label))!=NULL){
	pnR = (struct pnode *)l2->item;
	if (verbose>1){ printf(" %s pn1 label %d matches pnR label %d, descending...\n",lvlchar,pn1->region->label,pnR->region->label);}
	ptreeobsdist(nbases,pR,pnR,pnR->childllitem,nb,p1,pn1,pn1->childllitem,maxlevel,level+1);}
      else{
	pnR = NULL; 
	if (verbose>1){ printf(" %s pn1 label %d not found,terminating\n",lvlchar,pn1->region->label);}}}}
  if (firstentry){
    if (parent1!=NULL){ 
      assert(parentR!=NULL);
      if (parentR->temp==NULL){ 
	Lra = (struct llist **) tcalloc(nbases,sizeof(struct llist *)); parentR->temp = (void *) Lra;
	for (nb2=0;nb2<nbases;nb2++){ Lra[nb2] = llistmake();}}
      Lra = (struct llist **) parentR->temp;
      if (fabs(parent1->weight)>wtol){ 
	if (verbose>2){ printf(" %s starting with llist:",lvlchar); llistprintf2(Lra[nb]);}
	temp=(double *)tcalloc(1,sizeof(double));*temp=parent1->weight;litemadd(Lra[nb],temp);
	if (verbose>2){ printf(" %s adding %f to llist %d at pR %d:",lvlchar,parent1->weight,nb,(int)parentR); llistprintf2(Lra[nb]);}}}}
}

void obsdistdisc(int nbases,int nrecords,struct llist **Lra,double *discriminability,double *reliability)
{
  /* given llistra (nbases elements) of observation distributions, calculates discriminability (assuming good statistics for now) */
  int verbose=0;
  int nb=0,nr=0,nv=0;
  double *maxra=NULL,max=0,max2=0,mean=0;
  struct hist **histra=NULL;
  struct litem *l=NULL;
  double nright=0,nwrong=0,nweak=0;
  double wtol=0.0000001;
  int nmax=0;
  maxra = (double *) tcalloc(nbases,sizeof(double));
  for (nb=0;nb<nbases;nb++){ lliststats(Lra[nb],&(maxra[nb]),NULL,NULL,NULL);}
  stats("double",maxra,nbases,&max,NULL,NULL,NULL); tfree(maxra);maxra=NULL;
  histra = (struct hist **) tcalloc(nbases,sizeof(struct hist *));
  for (nb=0;nb<nbases;nb++){ 
    if (verbose){ printf(" %% dealing with nb %d llist:\n",nb); llistprintf2(Lra[nb]);}
    histra[nb]=histmake((int)max,max,0); 
    l=Lra[nb]->first; nr=0; while(l!=NULL){ histadd(histra[nb],*(double *)l->item,1); nr+=1; l=l->child;}
    histadd(histra[nb],0,nrecords-nr);
    if (verbose){ printf(" %% now have hist:\n"); histprintf(histra[nb]," %% ");}}
  nright=0;nwrong=0;nweak=0;
  for (nv=0;nv<(int)max;nv++){
    maxra=(double *)tcalloc(nbases,sizeof(double));
    for (nb=0;nb<nbases;nb++){ if (nv>=histra[nb]->nbins){ maxra[nb]=0;} else{ maxra[nb]=histra[nb]->data[nv];}}
    if (verbose){ printf(" %% comparing bin %d: \n",nv); raprintf(maxra,"double",1,nbases," %% ");}
    stats("double",maxra,nbases,&max2,NULL,&mean,NULL);
    if ((nbases*mean)>1.5){
      nmax=0; for (nb=0;nb<nbases;nb++){ if (fabs(maxra[nb]-max2)<wtol){ nmax +=1;}}
      for (nb=0;nb<nbases;nb++){ 
	if (fabs(maxra[nb]-max2)<wtol){ nright += 1.0/(double)nmax*maxra[nb]; nwrong += (1-1.0/(double)nmax)*maxra[nb];}
	else{ nwrong += maxra[nb];}}
      if (verbose){ printf(" %% non-singleton, nmax %d nright %f nwrong %f\n",nmax,nright,nwrong);}}
    else if ((nbases*mean)<=1.5 && fabs(nbases*mean)>wtol){ 
      nweak += 1; if (verbose){ printf(" %% singleton, ignoring, nweak %f\n",nweak);}}
    else if (fabs(nbases*mean)<wtol){ if (verbose){ printf(" %% both empty, ignoring\n");}}
    tfree(maxra);maxra=NULL;}
  if (nright+nwrong>0){ *discriminability = nright/(nright+nwrong);} else{ *discriminability = 1.0/(double)nbases;}
  *reliability = 1-nweak/maximum(1,nright+nwrong+nweak);
  if (verbose){ printf(" %% now discriminability %f reliability %f\n",*discriminability,*reliability);}
  for (nb=0;nb<nbases;nb++){ histtfree(histra[nb]);histra[nb]=NULL;} tfree(histra);histra=NULL;
}

void ptreeobsdistend(int nbases,int nrecords,struct ptree *pR,struct pnode *parentR,struct llitem *lR,int level)
{
  /* checks discriminability of each event-chain, and then tfrees llistra (nbases elements) at each node */
  int verbose=0,nr=0,nb=0;
  char lvlchar[64];
  int firstentry=0;
  struct pnode *pnR=NULL;
  struct llist **Lra=NULL;
  if (verbose){ sprintf(lvlchar," %%"); for (nr=0;nr<level;nr++){sprintf(lvlchar,"%s %%",lvlchar);}}
  if (verbose>1){ 
    printf(" %s [entering ptreeobsdistend]\n",lvlchar);
    printf(" %s pR %d parentR %d->%d lR %d \n",lvlchar,(int)pR,(int)parentR,parentR==NULL ? 0 : (int)(parentR->childllitem),(int)lR);
    printf(" %s nbases %d nrecords %d level %d\n",lvlchar,nbases,nrecords,level);}
  if (parentR!=NULL && lR->parent==NULL){ /* first descent into parentR->childllitem */ 
    if (verbose>1){ printf(" %s first descent into parentR->childllitem, setting firstentry to 1\n",lvlchar);}
    assert(parentR->childllitem==lR); firstentry=1;}
  if (lR->kidl!=NULL){ 
    if (verbose>1){ printf(" %s lR->kidl found, moving left to new llitem lR->kidl=%d\n",lvlchar,(int)(lR->kidl));}
    ptreeobsdistend(nbases,nrecords,pR,parentR,lR->kidl,level);}
  if (lR->kidr!=NULL){ 
    if (verbose>1){ printf(" %s lR->kidr found, moving right to new llitem lR->kidr=%d\n",lvlchar,(int)(lR->kidr));}
    ptreeobsdistend(nbases,nrecords,pR,parentR,lR->kidr,level);}
  if (lR->item!=NULL){ 
    if (verbose>1){ printf(" %s lR->item %d exists\n",lvlchar,(int)(struct pnode *)(lR->item));}
    pnR = (struct pnode *)lR->item; 
    ptreeobsdistend(nbases,nrecords,pR,pnR,pnR->childllitem,level+1);}
  if (firstentry){
    if (parentR!=NULL){ 
      if (parentR->temp!=NULL){ 
	Lra = (struct llist **) parentR->temp;
	obsdistdisc(nbases,nrecords,Lra,&(parentR->weight),&(parentR->relevance));
	for (nb=0;nb<nbases;nb++){ llisttfree2(Lra[nb]); Lra[nb]=NULL;} tfree(Lra);
	if (verbose>1){ printf(" %s freeing Lra at pR\n",lvlchar);}}}
    else{ if (verbose>1){ printf(" %s Lra not found, weight must be 0\n",lvlchar);} parentR->weight=0;}}
}

void ptree_trialaverage_helper(int nbases,int maxlevel,int tclump,int verbose,int X_flag,int x_flag,int s_flag,int d_flag,char **filename_base,char *outputname)
{
  /* assumes that outputname does NOT start with "./" */
  int old_region_type;
  int nb=0,nb2=0,nmin=0; 
  double ntemp=0;
  char **filename=NULL,allstring[16];
  char *gs2=GLOBAL_STRING_2;
  char outputfilename[256],outputfilename2[256];
  char outputfilename_disc_x2p[256];
  char outputfilename_disc_depthent[256];
  char outputfilename_disc_breadthent[256];
  char outputfilename_obsdist[256];
  char plotfilename2[256],plotfilename_x2p[256],plotfilename_x2pper[256];
  char plotfilename_depthent[256],plotfilename_depthentper[256];
  char plotfilename_breadthent[256],plotfilename_breadthentper[256];
  FILE **fp=NULL,*op=NULL,*op2=NULL;
  struct ptree **p0=NULL,**p1=NULL,*pR=NULL,*ptmp=NULL;
  struct pnode **pn0=NULL;
  struct llitem **l0pos=NULL,**l0pre=NULL;
  double *nright_x2=NULL,*nwrong_x2=NULL;
  double *nright_x2p=NULL,*nwrong_x2p=NULL;
  double *nright_fabs=NULL,*nwrong_fabs=NULL;
  double *nright_fabs2=NULL,*nwrong_fabs2=NULL;
  double *nright_depthent=NULL,*nwrong_depthent=NULL;
  double *nright_breadthent=NULL,*nwrong_breadthent=NULL;
  double *temp_x2=NULL,*temp_x2p=NULL,*temp_fabs=NULL,*temp_fabs2=NULL,*temp_depthent=NULL,*temp_breadthent=NULL;
  double *discra=NULL;
  double *dist_x2p=NULL;
  double *nright_x2p_perstimulus=NULL,*nwrong_x2p_perstimulus=NULL;
  double *nright_depthent_perstimulus=NULL,*nwrong_depthent_perstimulus=NULL;
  double *nright_breadthent_perstimulus=NULL,*nwrong_breadthent_perstimulus=NULL;
  int bitbybit=0,continue_flag=0,level=0,tcounter=0;
  double wmax=0,wmin=0,rmax=0,rmin=0;
  double *percentages=NULL,*sbias=NULL,maxbias=0,biastol=0.00001;
  int nx,nX=0,nR=0;
  double *nbias=NULL;
  struct hist *temphist=NULL;
  int nelements=0;
  double *referencenormalizers=NULL;
  double *observednormalizers=NULL;
  double *discdata_x2p=NULL;
  double *discdata_depthent=NULL;
  double *discdata_breadthent=NULL;
  if (verbose>1){ 
    printf(" %% [entering ptree_trialaverage_helper]\n");
    printf(" %% nbases %d maxlevel %d tclump %d verbose %d X_flag %d x_flag %d s_flag %d d_flag %d\n",nbases,maxlevel,tclump,verbose,X_flag,x_flag,s_flag,d_flag);
    for (nb=0;nb<nbases;nb++){ printf(" %% filename_base[%d]=%s\n",nb,filename_base[nb]);}
    printf(" %% outputname=%s\n",outputname);}
  if (outputname==NULL){ 
    sprintf(outputfilename,"./outprint_T%d_s%d_d%d_%srecord.m",tclump,s_flag,d_flag,gs2); 
    sprintf(outputfilename2,"./outdist_T%d_s%d_d%d_%srecord.m",tclump,s_flag,d_flag,gs2); 
    sprintf(outputfilename_disc_x2p,"./outdiscdata_x2p_T%d_s%d_d%d_%srecord",tclump,s_flag,d_flag,gs2); 
    sprintf(outputfilename_disc_depthent,"./outdiscdata_depthent_T%d_s%d_d%d_%srecord",tclump,s_flag,d_flag,gs2); 
    sprintf(outputfilename_disc_breadthent,"./outdiscdata_breadthent_T%d_s%d_d%d_%srecord",tclump,s_flag,d_flag,gs2); 
    sprintf(outputfilename_obsdist,"./outobsdist_s%d_d%d_%srecord",s_flag,d_flag,gs2);
    sprintf(plotfilename2,"./outdist_s%d_d%d_%srecord",s_flag,d_flag,gs2);
    sprintf(plotfilename_x2p,"./outprint_x2p_s%d_d%d_%srecord",s_flag,d_flag,gs2);
    sprintf(plotfilename_x2pper,"./outprintper_x2p_s%d_d%d_%srecord",s_flag,d_flag,gs2);
    sprintf(plotfilename_depthent,"./outprint_depthent_s%d_d%d_%srecord",s_flag,d_flag,gs2);
    sprintf(plotfilename_depthentper,"./outprintper_depthent_s%d_d%d_%srecord",s_flag,d_flag,gs2);
    sprintf(plotfilename_breadthent,"./outprint_breadthent_s%d_d%d_%srecord",s_flag,d_flag,gs2);
    sprintf(plotfilename_breadthentper,"./outprintper_breadthent_s%d_d%d_%srecord",s_flag,d_flag,gs2);}
  else{ 
    sprintf(outputfilename,"./outprint%s_T%d_s%d_d%d_%srecord.m",outputname,tclump,s_flag,d_flag,gs2); 
    sprintf(outputfilename2,"./outdist%s_T%d_s%d_d%d_%srecord.m",outputname,tclump,s_flag,d_flag,gs2); 
    sprintf(outputfilename_disc_x2p,"./outdiscdata_x2p_%s_T%d_s%d_d%d_%srecord",outputname,tclump,s_flag,d_flag,gs2); 
    sprintf(outputfilename_disc_depthent,"./outdiscdata_depthent_%s_T%d_s%d_d%d_%srecord",outputname,tclump,s_flag,d_flag,gs2); 
    sprintf(outputfilename_disc_breadthent,"./outdiscdata_breadthent_%s_T%d_s%d_d%d_%srecord",outputname,tclump,s_flag,d_flag,gs2); 
    sprintf(outputfilename_obsdist,"./outobsdist_%s_s%d_d%d_%srecord",outputname,s_flag,d_flag,gs2);
    sprintf(plotfilename2,"./outdist%s_s%d_d%d_%srecord",outputname,s_flag,d_flag,gs2);
    sprintf(plotfilename_x2p,"./outprint_x2p_%s_s%d_d%d_%srecord",outputname,s_flag,d_flag,gs2);
    sprintf(plotfilename_x2pper,"./outprintper_x2p_%s_s%d_d%d_%srecord",outputname,s_flag,d_flag,gs2);
    sprintf(plotfilename_depthent,"./outprint_depthent_%s_s%d_d%d_%srecord",outputname,s_flag,d_flag,gs2);
    sprintf(plotfilename_depthentper,"./outprintper_depthent_%s_s%d_d%d_%srecord",outputname,s_flag,d_flag,gs2);
    sprintf(plotfilename_breadthent,"./outprint_breadthent_%s_s%d_d%d_%srecord",outputname,s_flag,d_flag,gs2);
    sprintf(plotfilename_breadthentper,"./outprintper_breadthent_%s_s%d_d%d_%srecord",outputname,s_flag,d_flag,gs2);}
  if ((op=fopen(outputfilename,"w"))==NULL){ if (verbose>1){ printf(" %% writing to stdout\n");} op=stdout;}
  if ((op2=fopen(outputfilename2,"w"))==NULL){ if (verbose>1){ printf(" %% writing to stdout\n");} op2=stdout;}
  sprintf(allstring,"_s%d_all",s_flag);
  nright_x2 = (double *) tcalloc(maxlevel,sizeof(double));
  nwrong_x2 = (double *) tcalloc(maxlevel,sizeof(double));
  nright_x2p = (double *) tcalloc(maxlevel,sizeof(double));
  nwrong_x2p = (double *) tcalloc(maxlevel,sizeof(double));
  nright_fabs = (double *) tcalloc(maxlevel,sizeof(double));
  nwrong_fabs = (double *) tcalloc(maxlevel,sizeof(double));
  nright_fabs2 = (double *) tcalloc(maxlevel,sizeof(double));
  nwrong_fabs2 = (double *) tcalloc(maxlevel,sizeof(double));
  nright_depthent = (double *) tcalloc(maxlevel,sizeof(double));
  nwrong_depthent = (double *) tcalloc(maxlevel,sizeof(double));
  nright_breadthent = (double *) tcalloc(maxlevel,sizeof(double));
  nwrong_breadthent = (double *) tcalloc(maxlevel,sizeof(double));
  discra = (double *) tcalloc(maxlevel*nbases,sizeof(double));
  filename = (char **) tcalloc(nbases,sizeof(char *));
  for (nb=0;nb<nbases;nb++){ 
    filename[nb] = (char *) tcalloc(256,sizeof(char));}
  fp = (FILE **) tcalloc(nbases,sizeof(FILE *));
  p0 = (struct ptree **) tcalloc(nbases,sizeof(struct ptree *));
  p1 = (struct ptree **) tcalloc(nbases,sizeof(struct ptree *));
  temp_x2 = (double *) tcalloc(nbases,sizeof(double));
  temp_x2p = (double *) tcalloc(nbases,sizeof(double));
  temp_fabs = (double *) tcalloc(nbases,sizeof(double));
  temp_fabs2 = (double *) tcalloc(nbases,sizeof(double));
  temp_depthent = (double *) tcalloc(nbases,sizeof(double));
  temp_breadthent = (double *) tcalloc(nbases,sizeof(double));
  dist_x2p = (double *) tcalloc(nbases*nbases*maxlevel,sizeof(double)); /* deduced_base + actual_base*nbases + level*nbases*nbases */
  nright_x2p_perstimulus = (double *) tcalloc(maxlevel*nbases,sizeof(double));
  nwrong_x2p_perstimulus = (double *) tcalloc(maxlevel*nbases,sizeof(double));
  nright_depthent_perstimulus = (double *) tcalloc(maxlevel*nbases,sizeof(double));
  nwrong_depthent_perstimulus = (double *) tcalloc(maxlevel*nbases,sizeof(double));
  nright_breadthent_perstimulus = (double *) tcalloc(maxlevel*nbases,sizeof(double));
  nwrong_breadthent_perstimulus = (double *) tcalloc(maxlevel*nbases,sizeof(double));
  referencenormalizers = (double *) tcalloc(maxlevel*nbases,sizeof(double));
  observednormalizers = (double *) tcalloc(maxlevel,sizeof(double));
  old_region_type=GLOBAL_PTREE_REGION_TYPE; GLOBAL_PTREE_REGION_TYPE=0;
  if (verbose>1){ printf(" %% setting GLOBAL_PTREE_REGION_TYPE to %d\n",GLOBAL_PTREE_REGION_TYPE);}
  for (nb=0;nb<nbases;nb++){
    if (verbose){ printf(" %% checking for %s%s... ",filename_base[nb],allstring);}
    sprintf(filename[nb],"./%s%s",filename_base[nb],allstring);
    if (checktofind(filename[nb])){ if (verbose){ printf("found\n");} p0[nb] = ptreadback(filename[nb]);}
    else /* if not found */{
      if (verbose){ printf("not found, making...\n");}
      bitbybit=1;
      sprintf(filename[nb],"%s_%d",filename_base[nb],bitbybit);
      p0[nb] = ptreadback(filename[nb]); 
      if (s_flag==1){ if (verbose>1){ printf("shizuffling %d\n",s_flag);} pnodeshizuffle_starter(NULL,p0[nb]->postree,p0[nb]); pnodeshizuffle_starter(NULL,p0[nb]->pretree,p0[nb]);}
      bitbybit=2; continue_flag=1;
      do{
	sprintf(filename[nb],"./%s_%d",filename_base[nb],bitbybit);
	if (verbose>1){ printf(" %% trying %s... ",filename[nb]);}
	if (checktofind(filename[nb])){ 
	  if (verbose>1){ printf("found\n");}
	  continue_flag=1; 
	  p1[nb] = ptreadback(filename[nb]);
	  if (s_flag==1){ if (verbose>1){ printf("shizuffling %d\n",s_flag);} pnodeshizuffle_starter(NULL,p1[nb]->postree,p1[nb]); pnodeshizuffle_starter(NULL,p1[nb]->pretree,p1[nb]);}
	  ptreeplusequals_starter(p0[nb],p1[nb]);
	  ptreetfree(p1[nb]); p1[nb]=NULL;}
	else /* if not found */{ if (verbose>1){ printf("not found, stopping\n");} continue_flag=0;}
	bitbybit+=1;}
      while (continue_flag);
      sprintf(filename[nb],"./%s%s",filename_base[nb],allstring);
      if (verbose){ printf(" %% dumping %s\n",filename[nb]);}
      ptreerate(p0[nb]); ptreedump_starter(p0[nb],filename[nb],2,1,0,0,+1,-1);}
    if (verbose){ printf(" %% p0[%d]->total_time %f\n",nb,p0[nb]->total_time);}
    for (level=0;level<maxlevel;level++){ nelements=0; pnodestats_starter(NULL,p0[nb]->postree,level+1,level+1,0,&nelements,NULL,NULL,&(referencenormalizers[level+nb*maxlevel]),NULL,NULL,NULL,NULL,NULL); referencenormalizers[level+nb*maxlevel]*=nelements;}
    if (verbose>1){ for (level=0;level<maxlevel;level++){ printf(" %% p0[%d] level %d normalizer = %f\n",nb,level,referencenormalizers[level+nb*maxlevel]);}}}  
  if (X_flag>=1){
    sprintf(filename[0],"./%s_%d",filename_base[0],1);
    if (verbose>1){ printf(" %% reading time from %s... ",filename[0]);}
    if (checktofind(filename[0])){ p1[0] = ptreadback(filename[0]);} else{ printf(" %% warning, %s not found\n",filename[0]);}
    nR=0; for (nb=0;nb<nbases;nb++){ nR+=ceil(p0[nb]->total_time);} nR=ceil(2*nR/(double)(tclump*p1[0]->total_time));
    ptreetfree(p1[0]);p1[0]=NULL;
    if (verbose>1){ printf(" %% %% %% allocating discdata size %d with %d records\n",nbases*maxlevel*nR,nR);}
    discdata_x2p = (double *) tcalloc((nbases*maxlevel)*nR,sizeof(double));
    discdata_depthent = (double *) tcalloc((nbases*maxlevel)*nR,sizeof(double));
    discdata_breadthent = (double *) tcalloc((nbases*maxlevel)*nR,sizeof(double));
    nR=0;
    for (nb=0;nb<nbases;nb++){
      bitbybit=1; continue_flag=1;
      do{
	sprintf(filename[nb],"./%s_%d",filename_base[nb],bitbybit);
	if (verbose>1){ printf(" %% checking %s... ",filename[nb]);}
	if (checktofind(filename[nb])){ 
	  if (verbose>1){ printf("found, initiating tcounter\n");}
	  p1[nb] = ptreadback(filename[nb]);
	  if (s_flag==1){ if (verbose>1){ printf("shizuffling %d\n",s_flag);} pnodeshizuffle_starter(NULL,p1[nb]->postree,p1[nb]); pnodeshizuffle_starter(NULL,p1[nb]->pretree,p1[nb]);}
	  bitbybit += 1; continue_flag=1;}
	else /* if not found */{ if (verbose>1){ printf("not found, stopping\n");} continue_flag=0;}
	tcounter=1;
	while (tcounter<tclump && continue_flag){
	  sprintf(filename[nb],"./%s_%d",filename_base[nb],bitbybit);
	  if (verbose>1){ printf(" %% checking %s... ",filename[nb]);}
	  if (checktofind(filename[nb])){ 
	    if (verbose>1){ printf("found with tcounter %d\n",tcounter);}
	    bitbybit += 1; continue_flag=1;
	    ptmp = ptreadback(filename[nb]);
	    if (s_flag==1){ if (verbose>1){ printf("shizuffling %d\n",s_flag);} pnodeshizuffle_starter(NULL,ptmp->postree,ptmp); pnodeshizuffle_starter(NULL,ptmp->pretree,ptmp);}
	    ptreeplusequals_starter(p1[nb],ptmp);
	    ptreetfree(ptmp);ptmp=NULL;}
	  else /* if not found */{ if (verbose>1){ printf("not found, stopping\n");} continue_flag=0;}
	  tcounter += 1;}
	if (continue_flag){
	  if (verbose>1){ printf(" %% accumulated %d-%d total time %0.2f\n",bitbybit-1-tcounter+1,bitbybit-1,p1[nb]->total_time);}
	  nR += 1;
	  for (level=0;level<maxlevel;level++){ nelements=0; pnodestats_starter(NULL,p1[nb]->postree,level+1,level+1,0,&nelements,NULL,NULL,&(observednormalizers[level]),NULL,NULL,NULL,NULL,NULL); observednormalizers[level]*=nelements;}
	  if (verbose>2){ raprintf(observednormalizers,"double",1,maxlevel," %% %% observednormalizers ");}
	  for (level=1;level<=maxlevel;level++){
	    for (nb2=0;nb2<nbases;nb2++){
	      if (verbose>2){ printf(" %% %% %% nb %d level %d nb2 %d\n",nb,level,nb2);}
	      if (verbose>2){ printf(" %% %% %% calling ptreex2\n");}
	      temp_x2[nb2]=0; temp_x2p[nb2]=0; temp_fabs[nb2]=0; temp_fabs2[nb2]=0; temp_depthent[nb2]=0; temp_breadthent[nb2]=0;
	      ptreex2(p0[nb2],NULL,p0[nb2]->postree,p1[nb],NULL,p1[nb]->postree,nb==nb2 ? d_flag : 1,level,0,&(temp_x2[nb2]),&(temp_x2p[nb2]),&(temp_fabs[nb2]),&(temp_fabs2[nb2]));
	      if (verbose>2){ printf(" %% %% %% calling ptreerelent_depth\n");}
	      ptreerelent_depth(p0[nb2],NULL,p0[nb2]->postree,p1[nb],NULL,p1[nb]->postree,nb==nb2 ? d_flag : 1,level,0,&(temp_depthent[nb2]));
	      if (verbose>2){ printf(" %% %% %% calling ptreerelent_breadth\n");}
	      ptreerelent_breadth(p0[nb2],NULL,p0[nb2]->postree,p1[nb],NULL,p1[nb]->postree,nb==nb2 ? d_flag : 1,level,level,0,&(temp_breadthent[nb2]),&(referencenormalizers[0+nb2*maxlevel]),observednormalizers);}
	    if (verbose>2){ printf(" %% %% discdata tab %d+%d*nbases+%d*nbases*maxlevel\n",periodize(nb2-nb,0,nbases),level-1,nR-1);}
	    for (nb2=0;nb2<nbases;nb2++){
	      discdata_x2p[periodize(nb2-nb,0,nbases) + (level-1)*nbases + (nR-1)*nbases*maxlevel] = temp_x2p[nb2];
	      discdata_depthent[periodize(nb2-nb,0,nbases) + (level-1)*nbases + (nR-1)*nbases*maxlevel] = temp_depthent[nb2];
	      discdata_breadthent[periodize(nb2-nb,0,nbases) + (level-1)*nbases + (nR-1)*nbases*maxlevel] = temp_breadthent[nb2];}
	    ptreex2_disthist(&(dist_x2p[0 + nb*nbases + (level-1)*nbases*nbases]),temp_x2p,nbases);
	    if (verbose>1){
	      printf(" %% %% level %d, nb %d,temp_x2p[0]=%f, temp_x2p[1]=%f\n",level,nb,temp_x2p[0],temp_x2p[1]);
	      printf(" %% %% level %d, nb %d,temp_breadthent[0]=%f, temp_breadthent[1]=%f\n",level,nb,temp_breadthent[0],temp_breadthent[1]);
	      printf(" %% %% level %d, nb %d, temp_depthent[0]=%f, temp_depthent[1]=%f\n",level,nb,temp_depthent[0],temp_depthent[1]);}
	    ntemp = ptreex2_match(nbases,nb,temp_x2,NULL); nright_x2[level-1] += ntemp; nwrong_x2[level-1] += 1-ntemp;
	    ntemp = ptreex2_match(nbases,nb,temp_x2p,NULL); nright_x2p[level-1] += ntemp; nwrong_x2p[level-1] += 1-ntemp;
	    nright_x2p_perstimulus[(level-1)+nb*maxlevel] += ntemp; nwrong_x2p_perstimulus[(level-1)+nb*maxlevel] += 1-ntemp;
	    ntemp = ptreex2_match(nbases,nb,temp_fabs,NULL); nright_fabs[level-1] += ntemp; nwrong_fabs[level-1] += 1-ntemp;
	    ntemp = ptreex2_match(nbases,nb,temp_fabs2,NULL); nright_fabs2[level-1] += ntemp; nwrong_fabs2[level-1] += 1-ntemp;
	    ntemp = ptreex2_match(nbases,nb,temp_depthent,NULL); nright_depthent[level-1] += ntemp; nwrong_depthent[level-1] += 1-ntemp;
	    nright_depthent_perstimulus[(level-1)+nb*maxlevel] += ntemp; nwrong_depthent_perstimulus[(level-1)+nb*maxlevel] += 1-ntemp;
	    ntemp = ptreex2_match(nbases,nb,temp_breadthent,NULL); nright_breadthent[level-1] += ntemp; nwrong_breadthent[level-1] += 1-ntemp;
	    nright_breadthent_perstimulus[(level-1)+nb*maxlevel] += ntemp; nwrong_breadthent_perstimulus[(level-1)+nb*maxlevel] += 1-ntemp;
	    if (verbose>2){
	      printf(" %% ptreex2 for %s->postree nl %d got ",filename[nb],level);
	      for (nb2=0;nb2<nbases;nb2++){
		printf("x2 %f x2p %f fabs %f fabs2 %f ",temp_x2[nb2],temp_x2p[nb2],temp_fabs[nb2],temp_fabs2[nb2]);}
	      printf("correct index %d",nb);
	      printf("\n");}}}
	if (verbose>2){ printf(" %% %% freeing p1\n");} 
	if (p1[nb]!=NULL){ ptreetfree(p1[nb]);p1[nb]=NULL;} 
	if (verbose>2){ printf(" %% %% freed p1\n");}}
      while (continue_flag);}
    for (level=0;level<maxlevel;level++){ discra[level] = nright_x2p[level]/(nright_x2p[level]+nwrong_x2p[level]);}
    ra2jpg2(discra,"double",1,maxlevel,0,1.0,0.0,plotfilename_x2p);
    for (level=0;level<maxlevel;level++){ discra[level] = nright_depthent[level]/(nright_depthent[level]+nwrong_depthent[level]);}
    ra2jpg2(discra,"double",1,maxlevel,0,1.0,0.0,plotfilename_depthent);
    for (level=0;level<maxlevel;level++){ discra[level] = nright_breadthent[level]/(nright_breadthent[level]+nwrong_breadthent[level]);}
    ra2jpg2(discra,"double",1,maxlevel,0,1.0,0.0,plotfilename_breadthent);
    for (level=0;level<maxlevel;level++){ discra[level] = multithresh(&(dist_x2p[0 + 0*nbases + level*nbases*nbases]),nbases,0.95);}
    ra2jpg2(discra,"double",1,maxlevel,0,100,0.0,plotfilename2);
    for (nb=0;nb<nbases;nb++){ for (level=0;level<maxlevel;level++){ discra[level+nb*maxlevel]=nright_x2p_perstimulus[level+nb*maxlevel]/(nright_x2p_perstimulus[level+nb*maxlevel]+nwrong_x2p_perstimulus[level+nb*maxlevel]);}} ra2jpg2(discra,"double",maxlevel,nbases,1,1.0,0.0,plotfilename_x2pper);
    for (nb=0;nb<nbases;nb++){ for (level=0;level<maxlevel;level++){ discra[level+nb*maxlevel]=nright_depthent_perstimulus[level+nb*maxlevel]/(nright_depthent_perstimulus[level+nb*maxlevel]+nwrong_depthent_perstimulus[level+nb*maxlevel]);}} ra2jpg2(discra,"double",maxlevel,nbases,1,1.0,0.0,plotfilename_depthentper);
    for (nb=0;nb<nbases;nb++){ for (level=0;level<maxlevel;level++){ discra[level+nb*maxlevel]=nright_breadthent_perstimulus[level+nb*maxlevel]/(nright_breadthent_perstimulus[level+nb*maxlevel]+nwrong_breadthent_perstimulus[level+nb*maxlevel]);}} ra2jpg2(discra,"double",maxlevel,nbases,1,1.0,0.0,plotfilename_breadthentper);
    if (verbose){ 
      for (level=0;level<maxlevel;level++){ 
	printf(" %% level %d (%0.1f,%0.1f) or (%0.1f,%0.1f) or (%0.1f,%0.1f) or (%0.1f,%0.1f) or (%0.1f,%0.1f) or (%0.1f,%0.1f)\n",level+1,nright_x2[level],nwrong_x2[level],nright_x2p[level],nwrong_x2p[level],nright_fabs[level],nwrong_fabs[level],nright_fabs2[level],nwrong_fabs2[level],nright_depthent[level],nwrong_depthent[level],nright_breadthent[level],nwrong_breadthent[level]);}
      for (level=0;level<maxlevel;level++){ 
	fprintf(op," x2_%s_STLR(%d,%d,%d)=%0.1f; x2_%s_STLW(%d,%d,%d)=%0.1f;\n",gs2,s_flag+1,tclump,level+1,nright_x2p[level],gs2,s_flag+1,tclump,level+1,nwrong_x2p[level]);}
      fprintf(op," x2_%s_STL=x2_%s_STLR./(x2_%s_STLR + x2_%s_STLW);\n",gs2,gs2,gs2,gs2);
      fprintf(op," plot(permute(x2_%s_STL(s_to_plot,:,:),[3 2 1]),'-');\n",gs2);
      for (level=0;level<maxlevel;level++){ for (nb=0;nb<nbases;nb++){ for (nb2=0;nb2<nbases;nb2++){
	fprintf(op2,"dist_x2p_%s(%d,%d,%d)=%f;\n",gs2,nb2+1,nb+1,level+1,dist_x2p[nb2+nb*nbases+level*nbases*nbases]);}}}}
    if (verbose>1){ printf(" %% testing discdata_x2p\n");}
    discdata_threshold(discdata_x2p,nbases,maxlevel,nR,outputfilename_disc_x2p);
    if (verbose>1){ printf(" %% testing discdata_depthent\n");}
    discdata_threshold(discdata_depthent,nbases,maxlevel,nR,outputfilename_disc_depthent);
    if (verbose>1){ printf(" %% testing discdata_breadthent\n");}
    discdata_threshold(discdata_breadthent,nbases,maxlevel,nR,outputfilename_disc_breadthent);
    radump(discdata_x2p,"double",nbases*maxlevel,nR,outputfilename_disc_x2p);
    radump(discdata_depthent,"double",nbases*maxlevel,nR,outputfilename_disc_depthent);
    radump(discdata_breadthent,"double",nbases*maxlevel,nR,outputfilename_disc_breadthent);
    tfree(discdata_x2p);discdata_x2p=NULL;
    tfree(discdata_depthent);discdata_depthent=NULL;
    tfree(discdata_breadthent);discdata_breadthent=NULL;}
  if (x_flag==-1){ /* approximation of observation distribution for each event-chain */
    if (verbose){ printf(" %% rereading %s%s... ",filename_base[0],allstring);}
    sprintf(filename[0],"./%s%s",filename_base[0],allstring);
    if (checktofind(filename[0])){ if (verbose){ printf("found\n");} pR = ptreadback(filename[0]);}
    else /* if not found */{ if (verbose){ printf("warning! couldn't find %s, aborting\n",filename[0]);} exit(EXIT_SUCCESS);}
    for (nb=1;nb<nbases;nb++){
      if (verbose){ printf(" %% adding on %s%s...\n",filename_base[nb],allstring);}
      ptreeplusequals_starter(pR,p0[nb]);}
    pnodeclear_starter(NULL,pR->postree); pnodeclear_starter(NULL,pR->pretree);
    nR=0;
    for (nb=0;nb<nbases;nb++){
      bitbybit=1; continue_flag=1;
      do{
	sprintf(filename[nb],"./%s_%d",filename_base[nb],bitbybit);
	if (verbose>1){ printf(" %% checking %s... ",filename[nb]);}
	if (checktofind(filename[nb])){ 
	  if (verbose>1){ printf("found, initiating tcounter\n");}
	  p1[nb] = ptreadback(filename[nb]);
	  if (s_flag==1){ if (verbose>1){ printf("shizuffling %d\n",s_flag);} pnodeshizuffle_starter(NULL,p1[nb]->postree,p1[nb]); pnodeshizuffle_starter(NULL,p1[nb]->pretree,p1[nb]);}
	  bitbybit += 1; continue_flag=1;}
	else /* if not found */{ if (verbose>1){ printf("not found, stopping\n");} continue_flag=0;}
	tcounter=1;
	while (tcounter<tclump && continue_flag){
	  sprintf(filename[nb],"./%s_%d",filename_base[nb],bitbybit);
	  if (verbose>1){ printf(" %% checking %s... ",filename[nb]);}
	  if (checktofind(filename[nb])){ 
	    if (verbose>1){ printf("found with tcounter %d\n",tcounter);}
	    bitbybit += 1; continue_flag=1;
	    ptmp = ptreadback(filename[nb]);
	    if (s_flag==1){ if (verbose>1){ printf("shizuffling %d\n",s_flag);} pnodeshizuffle_starter(NULL,ptmp->postree,ptmp); pnodeshizuffle_starter(NULL,ptmp->pretree,ptmp);}
	    ptreeplusequals_starter(p1[nb],ptmp);
	    ptreetfree(ptmp);ptmp=NULL;}
	  else /* if not found */{ if (verbose>1){ printf("not found, stopping\n");} continue_flag=0;}
	  tcounter += 1;}
	if (continue_flag){
	  if (verbose>1){ printf(" %% accumulated %d-%d total time %0.2f\n",bitbybit-1-tcounter+1,bitbybit-1,p1[nb]->total_time);}
	  nR += 1;
	  if (verbose>1){ printf(" %% calling ptreeobsdist for %s->postree\n",filename[nb]);}
	  ptreeobsdist(nbases,pR,NULL,pR->postree,nb,p1[nb],NULL,p1[nb]->postree,maxlevel,0);}
	if (p1[nb]!=NULL){ ptreetfree(p1[nb]);p1[nb]=NULL;}}
      while (continue_flag);}
    ptreeobsdistend(nbases,nR,pR,NULL,pR->postree,0);
    temphist=histmake(128,1,0);
    pnode2hist_starter(NULL,pR->postree,-1,0,temphist,NULL);
    if (verbose>1){ histprintf(temphist,"discriminability: ");}
    histdump(temphist,0,outputfilename_obsdist," ",0);
    histtfree(temphist);temphist=NULL;
    ptreetfree(pR);pR=NULL;}
  if (x_flag==1){
    if (verbose){ printf(" %% rereading %s%s... ",filename_base[0],allstring);}
    sprintf(filename[0],"./%s%s",filename_base[0],allstring);
    if (checktofind(filename[0])){ if (verbose){ printf("found\n");} pR = ptreadback(filename[0]);}
    else /* if not found */{ if (verbose){ printf("warning! couldn't find %s, aborting\n",filename[0]);} exit(EXIT_SUCCESS);}
    for (nb=1;nb<nbases;nb++){
      if (verbose){ printf(" %% adding on %s%s...\n",filename_base[nb],allstring);}
      ptreeplusequals_starter(pR,p0[nb]);}
    pnodeclear_starter(NULL,pR->postree); pnodeclear_starter(NULL,pR->pretree);
    if (verbose>1){ printf(" %% making pn0 array and l0pos,l0pre arrays\n");}
    pn0 = (struct pnode **) tcalloc(nbases,sizeof(struct pnode *));
    l0pos = (struct llitem **) tcalloc(nbases,sizeof(struct llitem *));
    l0pre = (struct llitem **) tcalloc(nbases,sizeof(struct llitem *));
    for (nb=0;nb<nbases;nb++){ pn0[nb] = NULL; l0pos[nb] = p0[nb]->postree; l0pre[nb] = p0[nb]->pretree;}
    nR=0;
    for (nb=0;nb<nbases;nb++){
      bitbybit=1; continue_flag=1;
      do{
	sprintf(filename[nb],"./%s_%d",filename_base[nb],bitbybit);
	if (verbose>1){ printf(" %% checking %s... ",filename[nb]);}
	if (checktofind(filename[nb])){ 
	  if (verbose>1){ printf("found, initiating tcounter\n");}
	  p1[nb] = ptreadback(filename[nb]);
	  if (s_flag==1){ if (verbose>1){ printf("shizuffling %d\n",s_flag);} pnodeshizuffle_starter(NULL,p1[nb]->postree,p1[nb]); pnodeshizuffle_starter(NULL,p1[nb]->pretree,p1[nb]);}
	  bitbybit += 1; continue_flag=1;}
	else /* if not found */{ if (verbose>1){ printf("not found, stopping\n");} continue_flag=0;}
	tcounter=1;
	while (tcounter<tclump && continue_flag){
	  sprintf(filename[nb],"./%s_%d",filename_base[nb],bitbybit);
	  if (verbose>1){ printf(" %% checking %s... ",filename[nb]);}
	  if (checktofind(filename[nb])){ 
	    if (verbose>1){ printf("found with tcounter %d\n",tcounter);}
	    bitbybit += 1; continue_flag=1;
	    ptmp = ptreadback(filename[nb]);
	    if (s_flag==1){ if (verbose>1){ printf("shizuffling %d\n",s_flag);} pnodeshizuffle_starter(NULL,ptmp->postree,ptmp); pnodeshizuffle_starter(NULL,ptmp->pretree,ptmp);}
	    ptreeplusequals_starter(p1[nb],ptmp);
	    ptreetfree(ptmp);ptmp=NULL;}
	  else /* if not found */{ if (verbose>1){ printf("not found, stopping\n");} continue_flag=0;}
	  tcounter += 1;}
	if (continue_flag){
	  if (verbose>1){ printf(" %% accumulated %d-%d total time %0.2f\n",bitbybit-1-tcounter+1,bitbybit-1,p1[nb]->total_time);}
	  nR += 1;
	  if (verbose>1){ printf(" %% calling ptreex1 for %s->postree\n",filename[nb]);}
	  ptreex1(nbases,p0,pn0,l0pos,nb,p1[nb],NULL,p1[nb]->postree,pR,NULL,pR->postree,maxlevel,0);
	  if (verbose>1){ printf(" %% calling ptreex1 for %s->pretree\n",filename[nb]);}
	  ptreex1(nbases,p0,pn0,l0pre,nb,p1[nb],NULL,p1[nb]->pretree,pR,NULL,pR->pretree,maxlevel,0);}
	if (p1[nb]!=NULL){ ptreetfree(p1[nb]);p1[nb]=NULL;}}
      while (continue_flag);}
    pnodetimesd_starter(NULL,pR->postree,1.0/(double)nR,1.0/(double)nR);
    pnodetimesd_starter(NULL,pR->pretree,1.0/(double)nR,1.0/(double)nR);
    tfree(pn0); pn0=NULL; tfree(l0pos); l0pos=NULL; tfree(l0pre); l0pre=NULL;
    if (verbose>2){ pnodeprintf(NULL,pR->postree,maxlevel,0);}
    if (outputname==NULL){sprintf(filename[0],"./ptree_x1_t%d_%srecord%s",tclump,gs2,allstring);}
    else{ sprintf(filename[0],"./ptree_%s_x1_t%d_%srecord%s",outputname,tclump,gs2,allstring);}
    if (verbose){ printf(" %% dumping %s, with %d records found\n",filename[0],nR);}
    pnodestats_starter(NULL,pR->postree,-1,-1,0,NULL,&wmax,&wmin,NULL,NULL,&rmax,&rmin,NULL,NULL);
    ptreedump_starter(pR,filename[0],2,1,wmax,wmin,rmax,rmin);
    ptreetfree(pR);pR=NULL;}
  if (x_flag>=1){
    if (outputname==NULL){sprintf(filename[0],"./ptree_x1_t%d_%srecord%s",tclump,gs2,allstring);}
    else{ sprintf(filename[0],"./ptree_%s_x1_t%d_%srecord%s",outputname,tclump,gs2,allstring);}
    if (verbose){ printf(" %% checking %s...",filename[0]);}
    if (checktofind(filename[0])){
      if (verbose){ printf("found\n");}
      pR = ptreadback(filename[0]);}
    else /* if not found */{ if (verbose){ printf("not found\n");} exit(EXIT_SUCCESS);}
    if (verbose>1){ printf(" %% making pn0 array and l0pos,l0pre arrays and percentages and nbias\n");}
    pn0 = (struct pnode **) tcalloc(nbases,sizeof(struct pnode *));
    l0pos = (struct llitem **) tcalloc(nbases,sizeof(struct llitem *));
    l0pre = (struct llitem **) tcalloc(nbases,sizeof(struct llitem *));
    for (nb=0;nb<nbases;nb++){ pn0[nb] = NULL; l0pos[nb] = p0[nb]->postree; l0pre[nb] = p0[nb]->pretree;}
    nX=11;
    percentages = (double *) tcalloc(nX,sizeof(double)); 
    pnodestats_starter(NULL,pR->postree,-1,-1,0,NULL,&(percentages[nX-1]),&(percentages[0]),NULL,NULL,NULL,NULL,NULL,NULL);
    percentages[0]=0.50;
    for (nx=1;nx<nX-1;nx++){ percentages[nx] = percentages[0] + (percentages[nX-1]-percentages[0])*(double)nx/(double)(nX-1);}
    temphist=histmake(nX,percentages[nX-1],percentages[0]);
    pnodehist_starter(NULL,pR->postree,-1,0,temphist,NULL);
    if (verbose){
      for (nx=0;nx<nX;nx++){ fprintf(op," percentages(%d)=%0.16lf;\n",nx+1,percentages[nx]);}
      for (nx=0;nx<temphist->nbins;nx++){
	fprintf(op," x1_%s_hist_T%dS%d(%d) = %0.16lf;\n",gs2,tclump,s_flag,nx+1,temphist->data[nx]);}}
    histtfree(temphist);temphist=NULL;
    nbias = (double *) tcalloc(nX*maxlevel,sizeof(double));
    nR=0;
    for (nb=0;nb<nbases;nb++){
      bitbybit=1; continue_flag=1;
      do{
	sprintf(filename[nb],"./%s_%d",filename_base[nb],bitbybit);
	if (verbose>1){ printf(" %% checking %s... ",filename[nb]);}
	if (checktofind(filename[nb])){ 
	  if (verbose>1){ printf("found, initiating tcounter\n");}
	  p1[nb] = ptreadback(filename[nb]);
	  if (s_flag==1){ if (verbose>1){ printf("shizuffling %d\n",s_flag);} pnodeshizuffle_starter(NULL,p1[nb]->postree,p1[nb]); pnodeshizuffle_starter(NULL,p1[nb]->pretree,p1[nb]);}
	  bitbybit += 1; continue_flag=1;}
	else /* if not found */{ if (verbose>1){ printf("not found, stopping\n");} continue_flag=0;}
	tcounter=1;
	while (tcounter<tclump && continue_flag){
	  sprintf(filename[nb],"./%s_%d",filename_base[nb],bitbybit);
	  if (verbose>1){ printf(" %% checking %s... ",filename[nb]);}
	  if (checktofind(filename[nb])){ 
	    if (verbose>1){ printf("found with tcounter %d\n",tcounter);}
	    bitbybit += 1; continue_flag=1;
	    ptmp = ptreadback(filename[nb]);
	    if (s_flag==1){ if (verbose>1){ printf("shizuffling %d\n",s_flag);} pnodeshizuffle_starter(NULL,ptmp->postree,ptmp); pnodeshizuffle_starter(NULL,ptmp->pretree,ptmp);}
	    ptreeplusequals_starter(p1[nb],ptmp);
	    ptreetfree(ptmp);ptmp=NULL;}
	  else /* if not found */{ if (verbose>1){ printf("not found, stopping\n");} continue_flag=0;}
	  tcounter += 1;}
	if (continue_flag){
	  if (verbose>1){ printf(" %% accumulated %d-%d total time %0.2f\n",bitbybit-1-tcounter+1,bitbybit-1,p1[nb]->total_time);}
	  nR += 1;
	  if (verbose>1){ printf(" %% calling ptreex1sbias for %s->postree (base number %d)\n",filename[nb],nb);}
	  if (verbose>1){ printf(" %% creating sbias\n");}
	  sbias = (double *) tcalloc(nbases*nX*maxlevel,sizeof(double));
	  ptreex1sbias(nbases,p0,pn0,l0pos,nb,p1[nb],NULL,p1[nb]->postree,pR,NULL,pR->postree,maxlevel,0,nX,percentages,sbias);
	  if (verbose>2){ for (level=0;level<maxlevel;level++){ printf(" %% sbias nl %d given by:\n",level); raprintf(&(sbias[0 + 0*nbases + level*nbases*nX]),"double",nbases,nX," %% ");}}
	  if (verbose>2){ printf(" %% now turning nbias from:\n"); raprintf(nbias,"double",nX,maxlevel," %% ");}
	  for (level=0;level<maxlevel;level++){ for (nx=0;nx<nX;nx++){
	    stats("double",&(sbias[0 + nx*nbases + level*nbases*nX]),nbases,&maxbias,NULL,NULL,NULL);
	    if (verbose>1){ printf(" %% nl %d nx %d found maxbias %f...",level,nx,maxbias);}
	    nmin=0;
	    for (nb2=0;nb2<nbases;nb2++){
	      if (maxbias>0 && fabs(sbias[nb2 + nx*nbases + level*nbases*nX]-maxbias)<biastol){ nmin+=1;}}
	    if (maxbias>0 && fabs(sbias[nb + nx*nbases + level*nbases*nX]-maxbias)<biastol){
	      if (verbose>1){ printf("at correct index %d\n",nb);}
	      nbias[nx + level*nX] += 1.0/(double)nmin;}
	    else /* if (fabs(sbias[nb+nx*nbases+level*nbases*nX]-maxbias)>=biastol) */{
	      if (verbose>1){ printf("at wrong index\n");}}}}
	  if (verbose>2){ printf(" %% into:\n"); raprintf(nbias,"double",nX,maxlevel," %% ");}
	  tfree(sbias); sbias=NULL;}
	if (p1[nb]!=NULL){ ptreetfree(p1[nb]);p1[nb]=NULL;}}
      while (continue_flag);}
    if (verbose){ 
      printf(" %% number of records = %d\n",nR);
      fprintf(op," nR_%s_T(%d) = %d;\n",gs2,tclump,nR);
      for (level=0;level<maxlevel;level++){ for (nx=0;nx<nX;nx++){ 
	printf(" %% level %d nx %d percentage %0.2f nbias %0.1f\n",level+1,nx,percentages[nx],nbias[nx+level*nX]);}}
      for (level=0;level<maxlevel;level++){ for (nx=0;nx<nX;nx++){
	fprintf(op," x1_%s_STLPB(%d,%d,%d,%d)=%0.1f;\n",gs2,s_flag+1,tclump,level+1,nx+1,nbias[nx+level*nX]);}}
      for (level=0;level<maxlevel;level++){
	fprintf(op," temp=x1_%s_STLPB(%d,%d,1:%d,:)/nR_%s_T(%d);\n",gs2,s_flag+1,tclump,level+1,gs2,tclump);
	fprintf(op," x1_%s_STL(%d,%d,%d)=max(temp(:));\n",gs2,s_flag+1,tclump,level+1);}
      fprintf(op," plot(permute(x1_%s_STL(s_to_plot,:,:),[3 2 1]),':');\n",gs2);}
    tfree(nbias);nbias=NULL;
    tfree(percentages);percentages=NULL;
    tfree(pn0); pn0=NULL; tfree(l0pos); l0pos=NULL; tfree(l0pre); l0pre=NULL;
    ptreetfree(pR);pR=NULL;}
  for (nb=0;nb<nbases;nb++){
    ptreetfree(p0[nb]);p0[nb]=NULL;
    tfree(filename[nb]);}
  tfree(p0);tfree(p1);
  tfree(filename);
  tfree(temp_x2);
  tfree(temp_x2p);
  tfree(temp_fabs);
  tfree(temp_fabs2);
  tfree(temp_depthent);
  tfree(temp_breadthent);
  tfree(nright_x2);tfree(nwrong_x2);
  tfree(nright_x2p);tfree(nwrong_x2p);
  tfree(nright_fabs);tfree(nwrong_fabs);
  tfree(nright_fabs2);tfree(nwrong_fabs2);
  tfree(nright_depthent);tfree(nwrong_depthent);
  tfree(nright_breadthent);tfree(nwrong_breadthent);
  tfree(discra);
  tfree(dist_x2p);
  tfree(nright_x2p_perstimulus);tfree(nwrong_x2p_perstimulus);
  tfree(nright_depthent_perstimulus);tfree(nwrong_depthent_perstimulus);
  tfree(nright_breadthent_perstimulus);tfree(nwrong_breadthent_perstimulus);
  tfree(referencenormalizers);
  tfree(observednormalizers);
  tfree(fp);
  GLOBAL_PTREE_REGION_TYPE=old_region_type;
  if (op!=stdout){ fclose(op); op=NULL;}
  if (op2!=stdout){ fclose(op2); op2=NULL;}
}

void ptree_trialaverage(int argc,char **argv)
{
  int verbose=2;
  int still_have_options = argc-2;
  int nbases=0,nb=0;
  char helpfile[1024];
  char **filename_base=NULL;
  int maxlevel=1,tclump=1;
  int X_flag=1,x_flag=1,s_flag=0,d_flag=0;
  char *outputname=NULL;
  if (verbose>1){ printf(" %% [entering ptree_trialaverage]\n");}
  sprintf(helpfile,"%% [-Vverboselevel] -Nnumber_of_bases [-Mmaxlevel] [-Ttclump] [-XX_flag] [-xx_flag] [-ss_flag] [-dd_flag] -Fbase_name1 -Fbase_name2 ... [-Ooutputname]\n");
  if (still_have_options==0){ printf(helpfile); exit(EXIT_SUCCESS);}
  while (still_have_options){
    switch(getopt(argc,argv,"N:n:F:f:M:m:T:t:V:v:X:x:s:d:O:o:")){
    case 'M': case 'm':
      maxlevel = maximum(1,atoi(optarg));
      if (verbose){ printf(" %% maxlevel set to %d\n",maxlevel);}
      break;
    case 'T': case 't':
      tclump = maximum(1,atoi(optarg));
      if (verbose){ printf(" %% tclump set to %d\n",tclump);}
      break;
    case 'V': case 'v':
      verbose = atoi(optarg);
      if (verbose){ printf(" %% verboselevel set to %d\n",verbose);}
      break;
    case 'N': case 'n': 
      nbases = maximum(0,atoi(optarg));
      filename_base = (char **) tcalloc(nbases,sizeof(char *));
      for (nb=0;nb<nbases;nb++){ filename_base[nb] = (char *) tcalloc(256,sizeof(char));}
      nb=0;
      if (verbose){ printf(" %% nbases read as %d\n",nbases);}
      break;
    case 'F': case 'f': 
      sprintf(filename_base[nb],"%s",optarg); 
      if (verbose){ printf(" %% filename_base %d given by %s...",nb,filename_base[nb]);}
      nb = minimum(nbases-1,nb+1);
      if (verbose){ printf(" now nb %d\n",nb);}
      break;
    case 'X':
      X_flag = atoi(optarg);
      if (verbose){ printf(" %% now X_flag=%d\n",X_flag);}
      break;
    case 'x':
      x_flag = atoi(optarg);
      if (verbose){ printf(" %% now x_flag=%d\n",x_flag);}
      break;
    case 's':
      s_flag = atoi(optarg);
      if (verbose){ printf(" %% now s_flag=%d\n",s_flag);}
      break;
    case 'd':
      d_flag = atoi(optarg);
      if (verbose){ printf(" %% now d_flag=%d\n",d_flag);}
      break;
    case 'O': case 'o': 
      if (outputname==NULL){ outputname = (char *) tcalloc(128,sizeof(char));}
      sprintf(outputname,"%s",optarg);
      if (verbose){ printf(" %% now outputname=%s\n",outputname);}
      break;
    default: printf(helpfile); exit(EXIT_SUCCESS); break;}
    still_have_options--;}
  if (verbose>1){ printf(" %% calling ptree_trialaverage_helper\n");}
  ptree_trialaverage_helper(nbases,maxlevel,tclump,verbose,X_flag,x_flag,s_flag,d_flag,filename_base,outputname);
  for (nb=0;nb<nbases;nb++){ tfree(filename_base[nb]);} tfree(filename_base);
}

void pnode2distribution_starter(struct ptree *p0,struct pnode *parent0,struct llitem *l0,struct ptree *p1,struct pnode *parent1,struct llitem *l1,int maxlevel,int level)
{
  /* adds weights of tree *p1 to llist stored using *temp in reference tree *p0 */
  int firstentry=0;
  struct pnode *pn0=NULL,*pn1=NULL;
  struct llitem *l2=NULL;
  double *temp=NULL;
  if (parent1!=NULL && l1->parent==NULL){ /* first descent into parent1->childllitem */ 
    assert(parent1->childllitem==l1); firstentry=1;}
  if (l1->kidl!=NULL){ pnode2distribution_starter(p0,parent0,l0,p1,parent1,l1->kidl,maxlevel,level);}
  if (l1->kidr!=NULL){ pnode2distribution_starter(p0,parent0,l0,p1,parent1,l1->kidr,maxlevel,level);}
  if (l1->item!=NULL){
    if (maxlevel==-1 || level<maxlevel){ 
      pn1=(struct pnode *)l1->item; 
      if ((l2=llitemaddorfind(0,l0,pn1->region,&region2pnode_compare_label))!=NULL){
	pn0=(struct pnode *)l2->item;
	pnode2distribution_starter(p0,pn0,pn0->childllitem,p1,pn1,pn1->childllitem,maxlevel,level+1);}}}
  if (firstentry){ if (parent0!=NULL && parent1!=NULL){ 
    if (parent0->temp==NULL){ parent0->temp=(void *)llistmake();}
    temp = (double *) tcalloc(1,sizeof(double)); *temp = parent1->weight; litemadd((struct llist *)parent0->temp,temp);}}
}

void pnode2distribution_helper(struct ptree *p0,struct pnode *parent0,struct llitem *l0,int maxlevel,int level,struct llist **Lra,int nm,int nR,struct hist **Phistra)
{
  /* constructs nm llists from the weight llists stored in *temp of reference tree *p0 
     Lra[moment + level*nm] stores moments obtained by various observations at level 
     in addition constructs the "average" distribution Phistra, 
     assuming that (struct llist *)(parent0->temp)->length<=nR and Phistra is maxlevel elements long */
  int firstentry=0;
  struct pnode *pn0=NULL;
  int nr=0,nl=0;
  double *mra=NULL,*temp=NULL;
  struct llist *L=NULL;
  struct litem *l=NULL;
  if (parent0!=NULL && l0->parent==NULL){ /* first descent into parent0->childllitem */ assert(parent0->childllitem==l0); firstentry=1;}
  if (l0->kidl!=NULL){ pnode2distribution_helper(p0,parent0,l0->kidl,maxlevel,level,Lra,nm,nR,Phistra);}
  if (l0->kidr!=NULL){ pnode2distribution_helper(p0,parent0,l0->kidr,maxlevel,level,Lra,nm,nR,Phistra);}
  if (l0->item!=NULL){ if (maxlevel==-1||level<maxlevel){ pn0=(struct pnode *)l0->item; pnode2distribution_helper(p0,pn0,pn0->childllitem,maxlevel,level+1,Lra,nm,nR,Phistra);}}
  if (firstentry){ if (parent0!=NULL){ if (parent0->temp!=NULL){ 
    L = (struct llist *)parent0->temp;
    mra = (double *) tcalloc(maximum(3,nm),sizeof(double));
    lliststats2(L,mra,maximum(3,nm),nR-L->length);
    for (nr=0;nr<nm;nr++){ temp = (double *) tcalloc(1,sizeof(double)); *temp = mra[nr]; litemadd(Lra[nr + level*nm],temp);}
    l=L->first;nr=0;
    while (l!=NULL && nr<nR){
      for (nl=level-1;nl<maxlevel;nl++){ histadd(Phistra[nl],((*(double *)l->item)-mra[1])/sqrt(mra[2]),1);}
      nr+=1;
      l=l->child;}
    for (nl=level-1;nl<maxlevel;nl++){ histadd(Phistra[nl],(0-mra[1])/sqrt(mra[2]),nR-nr);}
    tfree(mra);mra=NULL;}}}
}

void pnode2distribution_ender(struct ptree *p0,struct pnode *parent0,struct llitem *l0,int maxlevel,int level)
{
  /* frees the llists stored in *temp of reference tree *p0 */
  int firstentry=0;
  struct pnode *pn0=NULL;
  if (parent0!=NULL && l0->parent==NULL){ /* first descent into parent0->childllitem */  assert(parent0->childllitem==l0); firstentry=1;}
  if (l0->kidl!=NULL){ pnode2distribution_ender(p0,parent0,l0->kidl,maxlevel,level);}
  if (l0->kidr!=NULL){ pnode2distribution_ender(p0,parent0,l0->kidr,maxlevel,level);}
  if (l0->item!=NULL){
    if (maxlevel==-1||level<maxlevel){ pn0=(struct pnode *)l0->item; pnode2distribution_ender(p0,pn0,pn0->childllitem,maxlevel,level+1);}}
  if (firstentry){ if (parent0!=NULL){ if (parent0->temp!=NULL){ llisttfree2((struct llist *)parent0->temp);parent0->temp=NULL;}}}
}

void ptree_trialdistribution_helper(int nbases,int maxlevel,int tclump,int verbose,char **filename_base,char *outputname)
{
  /* assumes that outputname does NOT start with "./" */
  int old_region_type;
  int s_flag=0;
  int nmoments=4;
  int nb=0,nm=0,nl=0;
  char *filename=NULL,allstring[16];
  struct ptree **p0=NULL,**p1=NULL,*ptmp=NULL;
  int bitbybit=0,continue_flag=0,tcounter=0;
  struct llist **Lra=NULL;
  double max=0,min=0,mean=0,stdev=0;
  struct hist *temphist=NULL;
  int nR=0;
  struct hist **Phistra=NULL;
  if (verbose>1){ 
    printf(" %% [entering ptree_trialdistribution_helper]\n");
    printf(" %% nbases %d maxlevel %d tclump %d verbose %d s_flag %d\n",nbases,maxlevel,tclump,verbose,s_flag);
    for (nb=0;nb<nbases;nb++){ printf(" %% filename_base[%d]=%s\n",nb,filename_base[nb]);}
    printf(" %% outputname=%s\n",outputname);}
  sprintf(allstring,"_s%d_all",s_flag);
  filename = (char *) tcalloc(256,sizeof(char));
  p0 = (struct ptree **) tcalloc(nbases,sizeof(struct ptree *));
  p1 = (struct ptree **) tcalloc(nbases,sizeof(struct ptree *));
  old_region_type=GLOBAL_PTREE_REGION_TYPE; GLOBAL_PTREE_REGION_TYPE=0;
  if (verbose>1){ printf(" %% setting GLOBAL_PTREE_REGION_TYPE to %d\n",GLOBAL_PTREE_REGION_TYPE);}
  for (nb=0;nb<nbases;nb++){
    if (verbose){ printf(" %% checking for %s%s... ",filename_base[nb],allstring);}
    sprintf(filename,"./%s%s",filename_base[nb],allstring);
    if (checktofind(filename)){ if (verbose){ printf("found\n");} p0[nb] = ptreadback(filename);}
    else /* if not found */{ printf("error, not found, exiting\n"); exit(EXIT_FAILURE);}}
  for (nb=0;nb<nbases;nb++){
    bitbybit=1; continue_flag=1;nR=0;
    do{
      sprintf(filename,"./%s_%d",filename_base[nb],bitbybit); if (verbose>1){ printf(" %% checking %s... ",filename);}
      if (checktofind(filename)){ 
	if (verbose>1){ printf("found, initiating tcounter\n");} 
	p1[nb] = ptreadback(filename); if (s_flag==1){ if (verbose>1){ printf("shizuffling %d\n",bitbybit);} pnodeshizuffle_starter(NULL,p1[nb]->postree,p1[nb]); pnodeshizuffle_starter(NULL,p1[nb]->pretree,p1[nb]);}
	bitbybit += 1; continue_flag=1;}
      else /* if not found */{ continue_flag=0;}
      tcounter=1;
      while (tcounter<tclump && continue_flag){
	sprintf(filename,"./%s_%d",filename_base[nb],bitbybit); if (verbose>1){ printf(" %% checking %s... ",filename);}
	if (checktofind(filename)){
	  if (verbose>1){ printf("found, with tcounter %d\n",tcounter);} 
	  ptmp = ptreadback(filename); if (s_flag==1){ if (verbose>1){ printf("shizuffling %d\n",bitbybit);} pnodeshizuffle_starter(NULL,ptmp->postree,ptmp); pnodeshizuffle_starter(NULL,ptmp->pretree,ptmp);}
	  ptreeplusequals_starter(p1[nb],ptmp);  
	  ptreetfree(ptmp);ptmp=NULL;
	  bitbybit += 1; continue_flag=1;}
	else /* if not found */{ if (verbose>1){ printf("not found, stopping\n");} continue_flag=0;}
	tcounter += 1;}
      if (continue_flag){
	nR += 1;
	if (verbose>1){ printf(" %% accumulated %d-%d total time %0.2f\n",bitbybit-1-tcounter+1,bitbybit-1,p1[nb]->total_time);}
	pnode2distribution_starter(p0[nb],NULL,p0[nb]->postree,p1[nb],NULL,p1[nb]->postree,maxlevel,0);}
      if (p1[nb]!=NULL){ ptreetfree(p1[nb]);p1[nb]=NULL;}}
    while (continue_flag);
    if (verbose){ printf(" %% found %d total records\n",nR);}
    Lra = (struct llist **) tcalloc((maxlevel+1)*nmoments,sizeof(struct llist *));
    for (nm=0;nm<nmoments*(maxlevel+1);nm++){ Lra[nm]=llistmake();}
    Phistra=(struct hist **) tcalloc(maxlevel,sizeof(struct hist *));
    for (nl=0;nl<maxlevel;nl++){ Phistra[nl]=histmake(128,0+16*STD_VIEW,0-1*STD_VIEW);}
    pnode2distribution_helper(p0[nb],NULL,p0[nb]->postree,maxlevel,0,Lra,nmoments,nR,Phistra);
    for (nl=0;nl<maxlevel;nl++){
      sprintf(filename,"./%s_dist_Phist_level%d",filename_base[nb],nl);
      histdump(Phistra[nl],0,filename," ",0); histdump(Phistra[nl],1,filename," ",0);}
    for (nl=0;nl<maxlevel;nl++){ histtfree(Phistra[nl]);Phistra[nl]=NULL;} tfree(Phistra);Phistra=NULL;
    for (nm=0;nm<nmoments;nm++){ for (nl=0;nl<maxlevel;nl++){ 
      lliststats(Lra[nm + nl*nmoments],&max,&min,&mean,&stdev);
      temphist = histmake(128,mean+STD_VIEW*stdev,mean-STD_VIEW*stdev);
      llist2hist(Lra[nm + nl*nmoments],temphist);
      sprintf(filename,"./%s_dist_moment%d_level%d",filename_base[nb],nm,nl); 
      histdump(temphist,0,filename," ",0); histdump(temphist,1,filename," ",0);
      histtfree(temphist);temphist=NULL;}}
    pnode2distribution_ender(p0[nb],NULL,p0[nb]->postree,maxlevel,0);
    for (nm=0;nm<nmoments*(maxlevel+1);nm++){ llisttfree2(Lra[nm]);Lra[nm]=NULL;}
    tfree(Lra);Lra=NULL;}
  for (nb=0;nb<nbases;nb++){
    ptreetfree(p0[nb]);p0[nb]=NULL;}
  tfree(p0);tfree(p1);
  tfree(filename);
  GLOBAL_PTREE_REGION_TYPE=old_region_type;
}

void ptree_trialdistribution(int argc,char **argv)
{
  int verbose=2;
  int still_have_options = argc-2;
  int nbases=0,nb=0;
  char helpfile[1024];
  char **filename_base=NULL;
  int maxlevel=1,tclump=1;
  char *outputname=NULL;
  if (verbose>1){ printf(" %% [entering ptree_trialdistribution]\n");}
  sprintf(helpfile,"%% [-Vverboselevel] -Nnumber_of_bases [-Mmaxlevel] [-Ttclump] -Fbase_name1 -Fbase_name2 ... [-Ooutputname]\n");
  if (still_have_options==0){ printf(helpfile); exit(EXIT_SUCCESS);}
  while (still_have_options){
    switch(getopt(argc,argv,"N:n:F:f:M:m:T:t:V:v:X:x:s:O:o:")){
    case 'M': case 'm':
      maxlevel = maximum(1,atoi(optarg));
      if (verbose){ printf(" %% maxlevel set to %d\n",maxlevel);}
      break;
    case 'T': case 't':
      tclump = maximum(1,atoi(optarg));
      if (verbose){ printf(" %% tclump set to %d\n",tclump);}
      break;
    case 'V': case 'v':
      verbose = atoi(optarg);
      if (verbose){ printf(" %% verboselevel set to %d\n",verbose);}
      break;
    case 'N': case 'n': 
      nbases = maximum(0,atoi(optarg));
      filename_base = (char **) tcalloc(nbases,sizeof(char *));
      for (nb=0;nb<nbases;nb++){ filename_base[nb] = (char *) tcalloc(256,sizeof(char));}
      nb=0;
      if (verbose){ printf(" %% nbases read as %d\n",nbases);}
      break;
    case 'F': case 'f': 
      sprintf(filename_base[nb],"%s",optarg); 
      if (verbose){ printf(" %% filename_base %d given by %s...",nb,filename_base[nb]);}
      nb = minimum(nbases-1,nb+1);
      if (verbose){ printf(" now nb %d\n",nb);}
      break;
    case 'O': case 'o': 
      if (outputname==NULL){ outputname = (char *) tcalloc(128,sizeof(char));}
      sprintf(outputname,"%s",optarg);
      if (verbose){ printf(" %% now outputname=%s\n",outputname);}
      break;
    default: printf(helpfile); exit(EXIT_SUCCESS); break;}
    still_have_options--;}
  if (verbose>1){ printf(" %% calling ptree_trialdistribution_helper\n");}
  ptree_trialdistribution_helper(nbases,maxlevel,tclump,verbose,filename_base,outputname);
  for (nb=0;nb<nbases;nb++){ tfree(filename_base[nb]);} tfree(filename_base);
}

/* Here are the closet functions */

struct closet * closetmake(int rmin,int rmax,int cmin,int cmax,int legtime,double minworth)
{
  int verbose=0;
  struct closet *c=NULL;
  int nr=0,nc=0;
  struct neuronarray *Nra=GLOBAL_Nra;
  struct neuron *n=NULL;
  if (verbose){ printf(" %% [entering closetmake]\n");}
  c = (struct closet *) tcalloc(1,sizeof(struct closet));
  c->neuronllist = llistmake();
  for (nr=rmin;nr<rmax;nr++){ for (nc=cmin;nc<cmax;nc++){
    n = nget(Nra,nr,nc);
    if (verbose){ printf(" %% adding neuron (%d,%d) to c->neuronllist\n",n->row,n->col);}
    litemadd(c->neuronllist,n);}}
  c->nlegs=0;
  c->legtime=legtime;
  c->minworth=minworth;
  c->p=NULL;
  if (verbose){ printf(" %% calling closetremake\n");}
  closetremake(c,0);
  return c;
}

void closettfree(struct closet *c)
{
  llisttfree(c->neuronllist); c->neuronllist=NULL;
  ptreetfree(c->p); c->p=NULL; 
  tfree(c); c=NULL;
}

void closetremake(struct closet *c,int legadd)
{
  int verbose=0;
  int baseregions=1000000;
  struct llist *Ll=NULL,*Ln=NULL;
  struct litem *l=NULL,*l2=NULL;
  struct neuron *n=NULL;
  int nregions=0,depth=0;
  if (verbose){ printf(" %% [entering closetremake] with legadd %d\n",legadd);}
  legadd=maximum(0,legadd);
  if (c->p==NULL){ /* do nothing */ c->nlegs=0;}
  else /* if (c->p!=NULL) */{
    nregions = maximum(2,minimum(c->neuronllist->length,(int)floor(pow(baseregions,1.0/(double)(c->nlegs+legadd+1)))));
    pnodebalance_starter(NULL,c->p->postree); pnodebalance_starter(NULL,c->p->pretree); ptreerate(c->p);
    Ll=ptreextract_starter(c->p,c->minworth,nregions); /* llistsort(Ll->first,Ll->last,Ll->length,&double_compare); */
    if (verbose){ l=Ll->first; while (l!=NULL){ printf(" %% detected label %d\n",(int)*(double *)l->item); l=l->child;}}
    Ln = llistmake();
    l=Ll->first; l2=c->neuronllist->first; depth=0;
    while (l!=NULL && l2!=NULL){
      if (depth==(int)*(double *)l->item){ litemadd(Ln,l2->item); l=l->child;}
      l2=l2->child; depth+=1;}
    llisttfree2(Ll); Ll=NULL;
    llisttfree(c->neuronllist); c->neuronllist=NULL;
    ptreetfree(c->p); c->p=NULL; 
    c->neuronllist=Ln; c->nlegs += legadd;}
  if (verbose){ printf(" %% c->neuronllist of length %d:\n",c->neuronllist->length);}
  c->p=ptreemake(c->neuronllist->length,1,1,0,c->nlegs,c->legtime);
  l=c->neuronllist->first;depth=0;
  while (l!=NULL && depth<c->p->nregions){
    n = (struct neuron *) l->item;
    if (verbose){ printf(" %% %% adding neuron (%d,%d) to region with label %d\n",n->row,n->col,c->p->regionra[depth]->label);}
    litemadd(c->p->regionra[depth]->neuronllist,n);
    l=l->child; depth+=1;}
}

/* Here are the yggdrasil functions */

struct yggdrasil * yggdrasilmake(int nregions,double event_within,int event_threshold,int region_type,int nlegs,int legtime,int ppnregions,int ppnlegs,int pplegtime,int weight_minimum)
{
  struct yggdrasil *Y=NULL;
  Y = (struct yggdrasil *) tmalloc(sizeof(struct yggdrasil));
  Y->p = ptreemake(nregions,event_within,event_threshold,region_type,nlegs,legtime);
  Y->ppnregions = ppnregions;
  Y->ppnlegs = ppnlegs;
  Y->pplegtime = pplegtime;
  Y->pp = NULL;
  Y->still_observing = 1;
  Y->weight_minimum = maximum(weight_minimum,1);
  return Y;
}

void yggdrasiltfree(struct yggdrasil *Y)
{
  ptreetfree(Y->p); Y->p=NULL;
  if (Y->pp!=NULL){ ptreetfree(Y->pp); Y->pp=NULL;}
  tfree(Y); Y=NULL;
}

int pnode2pnode_labels_in_common(struct pnode *pn0,struct pnode *pn1)
{
  /* two populated pnodes may have labels in common */
  struct llist *L0=llistmake();
  struct llist *L1=llistmake();
  struct litem *l0=NULL,*l1=NULL;
  double *temp=NULL;
  int in_common=0,compare=0;
  while (pn0!=NULL){ temp=(double *)tmalloc(sizeof(double));*temp=pn0->region->label;litemadd(L0,temp);pn0=pn0->parent;}
  while (pn1!=NULL){ temp=(double *)tmalloc(sizeof(double));*temp=pn1->region->label;litemadd(L1,temp);pn1=pn1->parent;}
  llistsort(L0->first,L0->last,L0->length,&double_compare);
  llistsort(L1->first,L1->last,L1->length,&double_compare);
  l0=L0->first; l1=L1->first; in_common=0;
  while (l0!=NULL && l1!=NULL){
    compare = double_compare(l0->item,l1->item);
    if (compare==1){ /* l0->item bigger */ l1=l1->child;}
    else if (compare==-1){ /* l1->item bigger */ l0=l0->child;}
    else if (compare==0){ /* l0->item same size as l1->item */ l0=l0->child; l1=l1->child; in_common += 1;}}
  llisttfree2(L0);L0=NULL;
  llisttfree2(L1);L1=NULL;
  return in_common;
}

void yggdrasilupdate(struct yggdrasil *Y,double t,double DT,double weight)
{
  /* updates ptree Y->p and yggdrasil *Y */
  int verbose=0;
  double wmax=0;
  double weight_maximum=0;
  struct llist *pnL=NULL;
  struct litem *l=NULL,*l2=NULL;
  struct ptree *p=NULL;
  struct pnode *pn=NULL,*pn2=NULL;
  int nr=0,in_common=0,in_common_threshold=2;
  if (verbose>2){ printf(" %% [starting yggdrasilupdate] with Y=%d, t=%0.3f, DT=%0.3f, weight=%0.3f\n",(int)Y,t,DT,weight);}
  if (verbose>2){ printf(" %% updating Y->p\n");}
  ptreeupdate(Y->p,t,DT,weight);
  if (verbose>2){ printf(" %% still_observing %d\n",Y->still_observing);}
  if (Y->still_observing){
    pnodestats_starter(NULL,Y->p->postree,1,-1,0,NULL,&wmax,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
    /* fix later */
    weight_maximum = pow(Y->weight_minimum*pow(Y->p->total_time/Y->p->legtime,Y->p->nlegs),1.0/(double)(Y->p->nlegs+1));
    weight_maximum = 1000;
    if (verbose>2){ printf(" %% max weight %0.1f, weight_maximum %0.1f\n",wmax,weight_maximum);}
    if (wmax>weight_maximum){
      pnodestats_starter(NULL,Y->p->postree,-1,-1,0,NULL,&wmax,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
      if (verbose>2){ printf(" %% more accurate max weight %0.1f, rating Y->p\n",wmax);}
      ptreerate(Y->p);
      if (verbose>2){ printf(" %% making Y->pp\n");}
      pnL = llistmake();
      pnode2llist_starter(NULL,Y->p->postree,-1,0,pnL);
      if (verbose>1){
	printf(" %% found relevant events, llist of length %d: \n",pnL->length);
	l=pnL->first;
	while (l!=NULL){
	  pn = (struct pnode *)l->item;
	  printf(" %% %%");
	  if (pn!=NULL){ printf(" weight %f, relevance %f ... %d",pn->weight,pn->relevance,pn->region->label); pn=pn->parent;}
	  while (pn!=NULL){ printf("<--%d",pn->region->label); pn=pn->parent;}
	  printf("\n");
	  l=l->child;}}
      llistsort(pnL->first,pnL->last,pnL->length,&pnode2pnode_compare_relevance);
      pn=NULL; if (pnL->first!=NULL){ pn = (struct pnode *)pnL->first->item;}
      while (pn!=NULL && pn->relevance<=0){ 
	litemminus(pnL,pnL->first->item); pn=NULL; if (pnL->first!=NULL){ pn = (struct pnode *)pnL->first->item;}}
      l=pnL->first;
      while (l!=NULL){
	pn=(struct pnode *)l->item;
	l2=l->child; in_common=0;
	while (l2!=NULL){
	  pn2 = (struct pnode *)l2->item;
	  in_common = maximum(in_common,pnode2pnode_labels_in_common(pn,pn2));
	  l2=l2->child;}
	if (in_common<in_common_threshold){ l=l->child;}
	else /* if (in_common>=in_common_threshold) */{ l=l->child; litemminus(pnL,l->parent->item);}}
      if (verbose>1){
	printf(" %% sorted and pruned llist of length %d: \n",pnL->length);
	l=pnL->first;
	while (l!=NULL){
	  pn = (struct pnode *)l->item;
	  printf(" %% %%");
	  if (pn!=NULL){ printf(" weight %f, relevance %f ... %d",pn->weight,pn->relevance,pn->region->label); pn=pn->parent;}
	  while (pn!=NULL){ printf("<--%d",pn->region->label); pn=pn->parent;}
	  printf("\n");
	  l=l->child;}}
      Y->pp = (struct ptree *) tmalloc(sizeof(struct ptree));
      p = Y->pp;
      l=pnL->last;nr=0;
      while(l!=NULL && nr<Y->ppnregions){ 
	pn=(struct pnode *)l->item; 
	if (pn->relevance<=0){ l=NULL;} 
	else{ l=l->parent;nr++;}}
      Y->ppnregions = maximum(1,nr);
      p->nregions = Y->ppnregions;
      p->regionra = (struct region **) tcalloc(p->nregions,sizeof(struct region*));
      l=pnL->last;nr=0;
      while (l!=NULL && nr<p->nregions){ /* default setup */
	p->regionra[nr] = (struct region *) tmalloc(sizeof(struct region));
	p->regionra[nr]->label = nr; /* critical to have label nr match location in regionra */
	p->regionra[nr]->last_event = t-rand01;
	p->regionra[nr]->event_within = 0;
	p->regionra[nr]->event_threshold = 0;
	p->regionra[nr]->neuronllist = llistmake(); /* not used */
	p->regionra[nr]->pn = (struct pnode *)l->item;	
	l=l->parent; nr++;}
      l=pnL->last;nr=0;
      while (l!=NULL && nr<p->nregions){ /* more careful setup */
	pn = (struct pnode *)l->item;
	if (pn->weight>Y->weight_minimum){
	  p->regionra[nr]->pn = (struct pnode *)l->item;	
	  l=l->parent; nr++;}
	else{ l=l->parent;}}
      if (verbose>1){
	printf(" %% found relevant events, regionra: \n");
	for (nr=0;nr<Y->pp->nregions;nr++){
	  pn = Y->pp->regionra[nr]->pn;
	  printf(" %% %%");
	  while (pn!=NULL){ printf(" <--%d",pn->region->label); pn=pn->parent;}
	  printf("\n %% %% weight %f, relevance %f\n",Y->pp->regionra[nr]->pn->weight,Y->pp->regionra[nr]->pn->relevance);}}
      llisttfree(pnL); pnL=NULL;
      p->nlegs=Y->ppnlegs;
      p->legtime=Y->pplegtime;
      p->length = p->nlegs*p->legtime+1;
      p->eventra = (struct llitem **) tmalloc(p->length*sizeof(struct llitem *));
      for (nr=0;nr<p->length;nr++){ p->eventra[nr]=llitemmake();}
      p->tab=0;
      p->update_every=1;
      p->update_last=t;
      p->total_time=0;
      p->gate_flag=1;
      ptreemakepospre(p);
      p->wh = histmake(1,+1,-1);
      p->rh = histmake(1,+1,-1);
      Y->still_observing=0;}
    else{ /* do nothing */ }}
  if (!Y->still_observing){
    /* update Y->pp */
    if (verbose>2){ printf(" %% updating Y->pp\n");}
    p = Y->pp;
    for (nr=0;nr<p->nregions;nr++){
      if (p->gate_flag && yggdrasil_region_has_event(p->regionra[nr],Y->p,t,DT) /* determine activity and set last_event_flag */){ 
	if (llitemaddorfind(0,p->eventra[p->tab],p->regionra[nr],&region2region_compare_label)!=NULL){ 
	  if (verbose>2){ printf(" %% region %d active, but already in llitem\n",nr);}}
	else{
	  if (verbose>2){ printf(" %% region %d active, adding to llitem\n",nr);}
	  llitemaddorfind(1,p->eventra[p->tab],p->regionra[nr],&region2region_compare_label);}}}
    p->total_time += DT;
    if (t >= p->update_last + p->update_every){ ptreeupdate_helper(p,t,DT,weight,verbose);}}
}

/* Here are the bonsai functions */

struct bonsai * bonsaimake(struct neuronarray *Nra,int nregions,double event_within,int event_threshold,int region_type,int nlegs,int legtime)
{
  struct bonsai *b=NULL;
  char fgvn[512];
  b = (struct bonsai *) tmalloc(sizeof(struct bonsai));
  b->Nra = Nra;
  b->p = (struct ptree **) tmalloc(2*sizeof(struct ptree *));
  b->p[0] = ptreemake(nregions,event_within,event_threshold,region_type,nlegs,legtime);
  b->p[1] = ptreemake(nregions,event_within,event_threshold,region_type,nlegs,legtime);
  sprintf(fgvn,"nra_%srecord_perturbed_save",GLOBAL_STRING_2); nrasave(b->Nra,fgvn);
  b->nvsf = 0;
  b->switch_last = GLOBAL_TI;
  b->switch_every = nlegs*legtime+1;
  return b;
}

void bonsaitfree(struct bonsai *b)
{
  char command[1024];
  ptreetfree(b->p[0]);b->p[0]=NULL;
  ptreetfree(b->p[1]);b->p[1]=NULL;
  tfree(b->p);b->p=NULL;
  sprintf(command,"""rm"" ./nra_%srecord_perturbed_save",GLOBAL_STRING_2); system(command);
  sprintf(command,"""rm"" ./nra_%srecord_original_save",GLOBAL_STRING_2); system(command);
  tfree(b);b=NULL;
}

void bonsaiupdate(struct bonsai *b,double t,double DT,double weight)
{
  int verbose=0;
  char fgvn[512],command[1024];
  int nr=0;
  struct llist *pnL=NULL;
  struct litem *l=NULL;
  struct pnode *pn=NULL;
  double wtotal=0,wrand=0;
  struct neuron *n=NULL;
  if (verbose){ printf(" %% [entering bonsaiupdate] with t=%0.1f, DT=%0.1f, weight %0.1f\n",t,DT,weight);}
  if (t >= b->switch_last + b->switch_every){ /* rate, then switch mode */
    if (verbose){ printf(" %% switching mode from %d to %d, rating b->p[%d]\n",b->nvsf,1-b->nvsf,b->nvsf);}
    ptreerate(b->p[b->nvsf]);
    if (b->nvsf==0){
      if (verbose){ printf(" %% just finished normal run, saving original and loading perturbed\n");}
      sprintf(fgvn,"nra_%srecord_original_save",GLOBAL_STRING_2); nrasave(b->Nra,fgvn);
      sprintf(fgvn,"nra_%srecord_perturbed_save",GLOBAL_STRING_2); nraload(b->Nra,fgvn);}
    else if (b->nvsf==1){
    if (verbose){ printf(" %% just finished perturbed run, loading original and overwriting perturbed\n");}
      sprintf(fgvn,"nra_%srecord_original_save",GLOBAL_STRING_2); nraload(b->Nra,fgvn);
      if (verbose){ printf(" %% shuffling save files\n");}
      sprintf(command,"mv ./nra_%srecord_original_save ./nra_%srecord_perturbed_save;",GLOBAL_STRING_2,GLOBAL_STRING_2); system(command);}
    b->nvsf = 1 - b->nvsf;
    if (verbose){ printf(" %% destroying b->p[%d]->eventra\n",b->nvsf);}
    for (nr=0;nr<b->p[b->nvsf]->length;nr++){ 
      llitemtfree(b->p[b->nvsf]->eventra[nr],NULL); b->p[b->nvsf]->eventra[nr]=NULL; b->p[b->nvsf]->eventra[nr]=llitemmake();} 
    b->p[b->nvsf]->tab=0;
    if (b->nvsf==1){
      if (verbose){ printf(" %% forcing an event...\n");}
      pnL=llistmake();
      pnode2llist_starter(NULL,b->p[0]->postree,1,0,pnL);
      if (verbose>1){
	printf(" %% found %d baseline events: \n",pnL->length);
	l=pnL->first;
	while (l!=NULL){
	  pn = (struct pnode *)l->item;
	  printf(" %% %%");
	  if (pn!=NULL){ printf(" weight %f, relevance %f ... %d",pn->weight,pn->relevance,pn->region->label); pn=pn->parent;}
	  while (pn!=NULL){ printf("<--%d",pn->region->label); pn=pn->parent;}
	  printf("\n");
	  l=l->child;}}
      wtotal=0;
      l=pnL->first; while (l!=NULL){ pn = (struct pnode *)l->item; wtotal += pn->weight; l=l->child;}
      if (wtotal>0){
	wrand = (int)floor(rand01*wtotal);
	if (verbose){ printf(" %% wtotal %0.1f, wrand %0.1f\n",wtotal,wrand);}
	l=pnL->first; while (l!=NULL){ pn = (struct pnode *)l->item; wrand -= pn->weight; if (wrand<0){ l=NULL;} else{ l=l->child;}}
	if (pn!=NULL){ 
	  if (verbose){ printf(" %% settled on region %d\n",pn->region->label);}
	  l=pn->region->neuronllist->first;
	  while (l!=NULL){
	    n = (struct neuron *)l->item;
	    *(n->V) = EIFvsIF ? 0.99*VOLTAGE_THRESHOLD_EIF : 0.99*VOLTAGE_THRESHOLD;
	    if (verbose){ printf(" %% setting n(%d,%d) to have voltage %0.1f\n",n->row,n->col,*(n->V));}
	    l=l->child;}}}
      else /* if (wtotal==0) */{ if (verbose){ printf(" %% no events yet, not forcing\n");}}
      llisttfree(pnL);pnL=NULL;}
    b->switch_last = maximum(t,b->switch_last + b->switch_every);
    if (verbose){ printf(" %% last switch set to %0.1f\n",b->switch_last);}}
  if (verbose){ printf(" %% updating ptree b->p[%d]\n",b->nvsf);}
  ptreeupdate(b->p[b->nvsf],t,DT,weight);
}

void bonsaidump(struct bonsai *b,char *fgvn,int dump_type)
{
  /* assumes fgvn does NOT start with "./" */
  char filename[512];
  char *gs2=GLOBAL_STRING_2;
  struct ptree *p2=NULL;
  int lmin=0,nl=0;
  double *ra0=NULL,*ra1=NULL,*ra2=NULL;
  FILE *fp=NULL;
  if (dump_type==2){
    ptreerate(b->p[0]);
    sprintf(filename,"./ptree_%s_0_%srecord",fgvn,gs2);
    ptreedump_starter(b->p[0],filename,dump_type,1,0,0,+1,-1);
    ptreerate(b->p[1]);
    sprintf(filename,"./ptree_%s_1_%srecord",fgvn,gs2);
    ptreedump_starter(b->p[1],filename,dump_type,1,0,0,+1,-1);}
  if (dump_type==0 || dump_type==1){
    p2 = ptreesubtptree_starter(b->p[0],b->p[1],1,0);
    ra0 = (double *) tcalloc(b->p[0]->nlegs+2,sizeof(double));
    ra1 = (double *) tcalloc(b->p[1]->nlegs+2,sizeof(double));
    ra2 = (double *) tcalloc(p2->nlegs+2,sizeof(double));
    lmin=minimum(p2->nlegs+2,minimum(b->p[0]->nlegs+2,b->p[1]->nlegs+2));
    for (nl=0;nl<lmin;nl++){
      pnodefrob_starter(NULL,b->p[0]->postree,nl-1,0,&(ra0[nl]),NULL);
      pnodefrob_starter(NULL,b->p[1]->postree,nl-1,0,&(ra1[nl]),NULL);
      pnodefrob_starter(NULL,p2->postree,nl-1,0,&(ra2[nl]),NULL);
      ra2[nl] = ra2[nl]/(ra0[nl]>0 ? ra0[nl] : 1);}
    if (dump_type==0){
      sprintf(filename,"./ptree_%s_norms_%srecord.m",fgvn,gs2);
      if ((fp=fopen(filename,"w"))!=NULL){
	for (nl=0;nl<lmin;nl++){ fprintf(fp,"ptree_%s_norms_%srecord(%d,:)= [%f %f %f];\n",fgvn,gs2,nl+1,ra0[nl],ra1[nl],ra2[nl]);}
	fprintf(fp,"figure;clf;hold on;;\n");
	fprintf(fp,"subplot(3,1,1);hold on;plot(ptree_%s_norms_%srecord(:,1);title('%s_0');\n",fgvn,gs2,fgvn);
	fprintf(fp,"subplot(3,1,2);hold on;plot(ptree_%s_norms_%srecord(:,2);title('%s_1');\n",fgvn,gs2,fgvn);
	fprintf(fp,"subplot(3,1,3);hold on;plot(ptree_%s_norms_%srecord(:,3);title('relative difference');\n",fgvn,gs2);
	fclose(fp);fp=NULL;}}
    else if (dump_type==1){
      sprintf(filename,"./ptree_%s_norm0_%srecord",fgvn,gs2);
      ra2jpg2(ra0,"double",1,lmin,0,0,0,filename);
      sprintf(filename,"./ptree_%s_norm1_%srecord",fgvn,gs2);
      ra2jpg2(ra1,"double",1,lmin,0,0,0,filename);
      sprintf(filename,"./ptree_%s_reldif_%srecord",fgvn,gs2);
      ra2jpg2(ra2,"double",1,lmin,0,1,0,filename);}
    tfree(ra0);ra0=NULL; tfree(ra1);ra1=NULL; tfree(ra2);ra2=NULL;
    ptreetfree(p2);p2=NULL;}
}

/* Here are the hydra functions */

void ptree_inputswitchon(int flag){ 
  if (LGN_BOTHER){ INPUT_SPACEANGLE+=pow(-1,flag)*INPUT_SPACEANGLE_BACON;} else{ LGN_BACKRATE*=1+flag*INPUT_SPACEANGLE_BACON;}}
void ptree_inputswitchoff(int flag){ 
  if (LGN_BOTHER){ INPUT_SPACEANGLE-=pow(-1,flag)*INPUT_SPACEANGLE_BACON;} else{ LGN_BACKRATE/=1+flag*INPUT_SPACEANGLE_BACON;}}

struct hydra * hydramake(double justontime,double stayontime,double dump_every,void (*switchon)(int),void (*switchoff)(int))
{
  struct hydra *p=NULL;
  p = (struct hydra *) tmalloc(sizeof(struct hydra));
  p->justontime = justontime;
  p->stayontime = stayontime;
  p->switch_every_two = 4*(p->justontime+p->stayontime);
  p->last_switch = GLOBAL_TI;
  p->dump_every = dump_every;
  p->pjuston_1 = ptreemake(GLOBAL_PTREE_NREGIONS,GLOBAL_PTREE_EVENT_WITHIN,GLOBAL_PTREE_EVENT_THRESHOLD,GLOBAL_PTREE_REGION_TYPE,GLOBAL_PTREE_NLEGS,GLOBAL_PTREE_LEGTIME);
  p->pstayon_1 = ptreemake(GLOBAL_PTREE_NREGIONS,GLOBAL_PTREE_EVENT_WITHIN,GLOBAL_PTREE_EVENT_THRESHOLD,GLOBAL_PTREE_REGION_TYPE,GLOBAL_PTREE_NLEGS,GLOBAL_PTREE_LEGTIME);
  p->pjustoff_1 = ptreemake(GLOBAL_PTREE_NREGIONS,GLOBAL_PTREE_EVENT_WITHIN,GLOBAL_PTREE_EVENT_THRESHOLD,GLOBAL_PTREE_REGION_TYPE,GLOBAL_PTREE_NLEGS,GLOBAL_PTREE_LEGTIME);
  p->pjuston_2 = ptreemake(GLOBAL_PTREE_NREGIONS,GLOBAL_PTREE_EVENT_WITHIN,GLOBAL_PTREE_EVENT_THRESHOLD,GLOBAL_PTREE_REGION_TYPE,GLOBAL_PTREE_NLEGS,GLOBAL_PTREE_LEGTIME);
  p->pstayon_2 = ptreemake(GLOBAL_PTREE_NREGIONS,GLOBAL_PTREE_EVENT_WITHIN,GLOBAL_PTREE_EVENT_THRESHOLD,GLOBAL_PTREE_REGION_TYPE,GLOBAL_PTREE_NLEGS,GLOBAL_PTREE_LEGTIME);
  p->pjustoff_2 = ptreemake(GLOBAL_PTREE_NREGIONS,GLOBAL_PTREE_EVENT_WITHIN,GLOBAL_PTREE_EVENT_THRESHOLD,GLOBAL_PTREE_REGION_TYPE,GLOBAL_PTREE_NLEGS,GLOBAL_PTREE_LEGTIME);
  p->pstayoff = ptreemake(GLOBAL_PTREE_NREGIONS,GLOBAL_PTREE_EVENT_WITHIN,GLOBAL_PTREE_EVENT_THRESHOLD,GLOBAL_PTREE_REGION_TYPE,GLOBAL_PTREE_NLEGS,GLOBAL_PTREE_LEGTIME);
  p->switchon = switchon;
  p->switchoff = switchoff;
  p->pjuston_1->gate_flag=0; p->pstayon_1->gate_flag=0; p->pjustoff_1->gate_flag=0;
  p->pjuston_2->gate_flag=0; p->pstayon_2->gate_flag=0; p->pjustoff_2->gate_flag=0;
  p->pstayoff->gate_flag=0;
  return p;
}

void hydratfree(struct hydra *p)
{
  ptreetfree(p->pjuston_1); p->pjuston_1=NULL;
  ptreetfree(p->pstayon_1); p->pstayon_1=NULL;
  ptreetfree(p->pjustoff_1); p->pjustoff_1=NULL;
  ptreetfree(p->pjuston_2); p->pjuston_2=NULL;
  ptreetfree(p->pstayon_2); p->pstayon_2=NULL;
  ptreetfree(p->pjustoff_2); p->pjustoff_2=NULL;
  ptreetfree(p->pstayoff); p->pstayoff=NULL;
  tfree(p);p=NULL;
}

void hydraupdate(struct hydra *p,double t,double DT,double weight)
{
  int verbose=0;
  double t2=0,t3=0;
  if (verbose){ printf(" %% [entering hydraupdate] t=%f,DT=%f,weight=%f\n",t,DT,weight);}
  ptreeupdate(p->pjustoff_2,t,DT,weight);
  ptreeupdate(p->pstayoff,t,DT,weight);
  ptreeupdate(p->pjuston_1,t,DT,weight);
  ptreeupdate(p->pstayon_1,t,DT,weight);
  ptreeupdate(p->pjustoff_1,t,DT,weight);
  ptreeupdate(p->pjuston_2,t,DT,weight);
  ptreeupdate(p->pstayon_2,t,DT,weight);
  t2 = t-p->last_switch;
  p->pjustoff_2->gate_flag=0;
  p->pstayoff->gate_flag=0;
  p->pjuston_1->gate_flag=0;
  p->pstayon_1->gate_flag=0;
  p->pjustoff_1->gate_flag=0;
  p->pjuston_2->gate_flag=0;
  p->pstayon_2->gate_flag=0;
  if (0){/* do nothing */}
  else if (t2>0 && t2<=(0+p->justontime)){ 
    if (verbose){ printf(" %% justoff_2 \n");}
    t3=t2-0;
    p->pjustoff_2->gate_flag=1;
    if (p->dump_every>0 && t3>0 && (int)floor(t3/p->dump_every)<(int)floor((t3+DT)/p->dump_every)){ if (p->pjustoff_2!=NULL){ if (GLOBAL_PTREE_BITBYBIT){ if (verbose){ printf(" %% %% dumping p->pjustoff_2\n");} ptreerate(p->pjustoff_2); ptreedump_starter(p->pjustoff_2,"justoff_2",2,0,0,0,+1,-1); if (verbose){ printf(" %% %% resetting p->pjustoff_2\n");} ptreereset(p->pjustoff_2);}}}}
  else if (t2>(0+p->justontime) && t2<=(0+p->justontime+p->stayontime)){ 
    if (verbose){ printf(" %% stayoff \n");}
    t3=t2-(0+p->justontime);
    p->pstayoff->gate_flag=1;
    if (p->dump_every>0 && t3>0 && (int)floor(t3/p->dump_every)<(int)floor((t3+DT)/p->dump_every)){ if (p->pstayoff!=NULL){ if (GLOBAL_PTREE_BITBYBIT){ if (verbose){ printf(" %% %% dumping p->pstayoff\n");} ptreerate(p->pstayoff); ptreedump_starter(p->pstayoff,"stayoff",2,0,0,0,+1,-1); if (verbose){ printf(" %% %% resetting p->pstayoff\n");} ptreereset(p->pstayoff);}}}}
  else if (t2>(0+p->justontime+p->stayontime) && t2<=(0+2*p->justontime+p->stayontime)){ 
    if (verbose){ printf(" %% juston_1 \n");}
    t3=t2-(0+p->justontime+p->stayontime);
    p->pjuston_1->gate_flag=1;
    if (p->dump_every>0 && t3>0 && (int)floor(t3/p->dump_every)<(int)floor((t3+DT)/p->dump_every)){ if (p->pjuston_1!=NULL){ if (GLOBAL_PTREE_BITBYBIT){ if (verbose){ printf(" %% %% dumping p->pjuston_1\n");} ptreerate(p->pjuston_1); ptreedump_starter(p->pjuston_1,"juston_1",2,0,0,0,+1,-1); if (verbose){ printf(" %% %% resetting p->pjuston_1\n");} ptreereset(p->pjuston_1);}}}}
  else if (t2>(0+2*p->justontime+p->stayontime) && t2<=(0+2*p->justontime+2*p->stayontime)){ 
    if (verbose){ printf(" %% stayon_1 \n");}
    t3=t2-(0+2*p->justontime+p->stayontime);
    p->pstayon_1->gate_flag=1;
    if (p->dump_every>0 && t3>0 && (int)floor(t3/p->dump_every)<(int)floor((t3+DT)/p->dump_every)){ if (p->pstayon_1!=NULL){ if (GLOBAL_PTREE_BITBYBIT){ if (verbose){ printf(" %% %% dumping p->pstayon_1\n");} ptreerate(p->pstayon_1); ptreedump_starter(p->pstayon_1,"stayon_1",2,0,0,0,+1,-1); if (verbose){ printf(" %% %% resetting p->pstayon_1\n");} ptreereset(p->pstayon_1);}}}}
  else if (t2>(0+2*p->justontime+2*p->stayontime) && t2<=(0+3*p->justontime+2*p->stayontime)){ 
    if (verbose){ printf(" %% justoff_1 \n");}
    t3=t2-(0+2*p->justontime+2*p->stayontime);
    p->pjustoff_1->gate_flag=1;
    if (p->dump_every>0 && t3>0 && (int)floor(t3/p->dump_every)<(int)floor((t3+DT)/p->dump_every)){ if (p->pjustoff_1!=NULL){ if (GLOBAL_PTREE_BITBYBIT){ if (verbose){ printf(" %% %% dumping p->pjustoff_1\n");} ptreerate(p->pjustoff_1); ptreedump_starter(p->pjustoff_1,"justoff_1",2,0,0,0,+1,-1); if (verbose){ printf(" %% %% resetting p->pjustoff_1\n");} ptreereset(p->pjustoff_1);}}}}
  else if (t2>(0+3*p->justontime+2*p->stayontime) && t2<=(0+3*p->justontime+3*p->stayontime)){ 
    if (verbose){ printf(" %% stayoff \n");}
    t3=t2-(0+3*p->justontime+2*p->stayontime);
    p->pstayoff->gate_flag=1;
    if (p->dump_every>0 && t3>0 && (int)floor(t3/p->dump_every)<(int)floor((t3+DT)/p->dump_every)){ if (p->pstayoff!=NULL){ if (GLOBAL_PTREE_BITBYBIT){ if (verbose){ printf(" %% %% dumping p->pstayoff\n");} ptreerate(p->pstayoff); ptreedump_starter(p->pstayoff,"stayoff",2,0,0,0,+1,-1); if (verbose){ printf(" %% %% resetting p->pstayoff\n");} ptreereset(p->pstayoff);}}}}
  else if (t2>(0+3*p->justontime+3*p->stayontime) && t2<=(0+4*p->justontime+3*p->stayontime)){ 
    if (verbose){ printf(" %% juston_2 \n");}
    t3=t2-(0+3*p->justontime+3*p->stayontime);
    p->pjuston_2->gate_flag=1;
    if (p->dump_every>0 && t3>0 && (int)floor(t3/p->dump_every)<(int)floor((t3+DT)/p->dump_every)){ if (p->pjuston_2!=NULL){ if (GLOBAL_PTREE_BITBYBIT){ if (verbose){ printf(" %% %% dumping p->pjuston_2\n");} ptreerate(p->pjuston_2); ptreedump_starter(p->pjuston_2,"juston_2",2,0,0,0,+1,-1); if (verbose){ printf(" %% %% resetting p->pjuston_2\n");} ptreereset(p->pjuston_2);}}}}
  else if (t2>(0+4*p->justontime+3*p->stayontime) && t2<=(0+4*p->justontime+4*p->stayontime)){ 
    if (verbose){ printf(" %% stayon_2 \n");}
    t3=t2-(0+4*p->justontime+3*p->stayontime);
    p->pstayon_2->gate_flag=1;
    if (p->dump_every>0 && t3>0 && (int)floor(t3/p->dump_every)<(int)floor((t3+DT)/p->dump_every)){ if (p->pstayon_2!=NULL){ if (GLOBAL_PTREE_BITBYBIT){ if (verbose){ printf(" %% %% dumping p->pstayon_2\n");} ptreerate(p->pstayon_2); ptreedump_starter(p->pstayon_2,"stayon_2",2,0,0,0,+1,-1); if (verbose){ printf(" %% %% resetting p->pstayon_2\n");} ptreereset(p->pstayon_2);}}}}
  if (t2<(0+1*p->justontime+1*p->stayontime) && t2+DT>=(0+1*p->justontime+1*p->stayontime)){ (p->switchon)(1);}
  if (t2<(0+2*p->justontime+2*p->stayontime) && t2+DT>=(0+2*p->justontime+2*p->stayontime)){ (p->switchoff)(1);}
  if (t2<(0+3*p->justontime+3*p->stayontime) && t2+DT>=(0+3*p->justontime+3*p->stayontime)){ (p->switchon)(2);}
  if (t2<(0+4*p->justontime+4*p->stayontime) && t2+DT>=(0+4*p->justontime+4*p->stayontime)){ (p->switchoff)(2);}
  if (t2<(0+4*p->justontime+4*p->stayontime) && t2+DT>=(0+4*p->justontime+4*p->stayontime)){ p->last_switch += p->switch_every_two;}
}

void hydradump(struct hydra *p)
{
  /* assumes that GLOBAL_PTREE_BITBYBIT was on during evolution, also SUITE_TINDEXMAX should equal 1 */
  char **fnamebase=NULL;
  char *gs2=GLOBAL_STRING_2;
  int tindex=0,sindex=0,dindex=0;
  char outputname[256];
  int bitbybit=0,continue_flag=0;  
  if (GLOBAL_PTREE_BITBYBIT){
    fnamebase = (char **) tcalloc(7,sizeof(char *));
    fnamebase[0] = (char *) tcalloc(256,sizeof(char));
    fnamebase[1] = (char *) tcalloc(256,sizeof(char));
    fnamebase[2] = (char *) tcalloc(256,sizeof(char));
    fnamebase[3] = (char *) tcalloc(256,sizeof(char));
    fnamebase[4] = (char *) tcalloc(256,sizeof(char));
    fnamebase[5] = (char *) tcalloc(256,sizeof(char));
    fnamebase[6] = (char *) tcalloc(256,sizeof(char));
    sprintf(fnamebase[0],"ptree_juston_1_%srecord",gs2); sprintf(fnamebase[1],"ptree_stayon_1_%srecord",gs2);
    for (tindex=1;tindex<=SUITE_TINDEXMAX;tindex++){ for (sindex=0;sindex<=SUITE_SINDEXMAX;sindex++){ for (dindex=0;dindex<=SUITE_DINDEXMAX;dindex++){
      sprintf(outputname,"j1vss1"); ptree_trialaverage_helper(2,GLOBAL_PTREE_NLEGS+1,tindex,1,1,0,sindex,dindex,&fnamebase[0],outputname);}}}
    sprintf(fnamebase[0],"ptree_juston_1_%srecord",gs2); sprintf(fnamebase[1],"ptree_justoff_1_%srecord",gs2);
    for (tindex=1;tindex<=SUITE_TINDEXMAX;tindex++){ for (sindex=0;sindex<=SUITE_SINDEXMAX;sindex++){ for (dindex=0;dindex<=SUITE_DINDEXMAX;dindex++){
      sprintf(outputname,"j1vso1"); ptree_trialaverage_helper(2,GLOBAL_PTREE_NLEGS+1,tindex,1,1,0,sindex,dindex,&fnamebase[0],outputname);}}}
    sprintf(fnamebase[0],"ptree_juston_1_%srecord",gs2); sprintf(fnamebase[1],"ptree_stayoff_%srecord",gs2);
    for (tindex=1;tindex<=SUITE_TINDEXMAX;tindex++){ for (sindex=0;sindex<=SUITE_SINDEXMAX;sindex++){ for (dindex=0;dindex<=SUITE_DINDEXMAX;dindex++){
      sprintf(outputname,"j1vsoo"); ptree_trialaverage_helper(2,GLOBAL_PTREE_NLEGS+1,tindex,1,1,0,sindex,dindex,&fnamebase[0],outputname);}}}
    sprintf(fnamebase[0],"ptree_juston_1_%srecord",gs2); sprintf(fnamebase[1],"ptree_juston_2_%srecord",gs2);
    for (tindex=1;tindex<=SUITE_TINDEXMAX;tindex++){ for (sindex=0;sindex<=SUITE_SINDEXMAX;sindex++){ for (dindex=0;dindex<=SUITE_DINDEXMAX;dindex++){
      sprintf(outputname,"j1vsj2"); ptree_trialaverage_helper(2,GLOBAL_PTREE_NLEGS+1,tindex,1,1,0,sindex,dindex,&fnamebase[0],outputname);}}}
    sprintf(fnamebase[0],"ptree_juston_1_%srecord",gs2); sprintf(fnamebase[1],"ptree_stayon_2_%srecord",gs2);
    for (tindex=1;tindex<=SUITE_TINDEXMAX;tindex++){ for (sindex=0;sindex<=SUITE_SINDEXMAX;sindex++){ for (dindex=0;dindex<=SUITE_DINDEXMAX;dindex++){
      sprintf(outputname,"j1vss2"); ptree_trialaverage_helper(2,GLOBAL_PTREE_NLEGS+1,tindex,1,1,0,sindex,dindex,&fnamebase[0],outputname);}}}
    sprintf(fnamebase[0],"ptree_juston_1_%srecord",gs2); sprintf(fnamebase[1],"ptree_justoff_2_%srecord",gs2);
    for (tindex=1;tindex<=SUITE_TINDEXMAX;tindex++){ for (sindex=0;sindex<=SUITE_SINDEXMAX;sindex++){ for (dindex=0;dindex<=SUITE_DINDEXMAX;dindex++){
      sprintf(outputname,"j1vso2"); ptree_trialaverage_helper(2,GLOBAL_PTREE_NLEGS+1,tindex,1,1,0,sindex,dindex,&fnamebase[0],outputname);}}}
    sprintf(fnamebase[0],"ptree_stayon_1_%srecord",gs2); sprintf(fnamebase[1],"ptree_justoff_1_%srecord",gs2);
    for (tindex=1;tindex<=SUITE_TINDEXMAX;tindex++){ for (sindex=0;sindex<=SUITE_SINDEXMAX;sindex++){ for (dindex=0;dindex<=SUITE_DINDEXMAX;dindex++){
      sprintf(outputname,"s1vso1"); ptree_trialaverage_helper(2,GLOBAL_PTREE_NLEGS+1,tindex,1,1,0,sindex,dindex,&fnamebase[0],outputname);}}}
    sprintf(fnamebase[0],"ptree_stayon_1_%srecord",gs2); sprintf(fnamebase[1],"ptree_stayoff_%srecord",gs2);
    for (tindex=1;tindex<=SUITE_TINDEXMAX;tindex++){ for (sindex=0;sindex<=SUITE_SINDEXMAX;sindex++){ for (dindex=0;dindex<=SUITE_DINDEXMAX;dindex++){
      sprintf(outputname,"s1vsoo"); ptree_trialaverage_helper(2,GLOBAL_PTREE_NLEGS+1,tindex,1,1,0,sindex,dindex,&fnamebase[0],outputname);}}}
    sprintf(fnamebase[0],"ptree_stayon_1_%srecord",gs2); sprintf(fnamebase[1],"ptree_juston_2_%srecord",gs2);
    for (tindex=1;tindex<=SUITE_TINDEXMAX;tindex++){ for (sindex=0;sindex<=SUITE_SINDEXMAX;sindex++){ for (dindex=0;dindex<=SUITE_DINDEXMAX;dindex++){
      sprintf(outputname,"s1vsj2"); ptree_trialaverage_helper(2,GLOBAL_PTREE_NLEGS+1,tindex,1,1,0,sindex,dindex,&fnamebase[0],outputname);}}}
    sprintf(fnamebase[0],"ptree_stayon_1_%srecord",gs2); sprintf(fnamebase[1],"ptree_stayon_2_%srecord",gs2);
    for (tindex=1;tindex<=SUITE_TINDEXMAX;tindex++){ for (sindex=0;sindex<=SUITE_SINDEXMAX;sindex++){ for (dindex=0;dindex<=SUITE_DINDEXMAX;dindex++){
      sprintf(outputname,"s1vss2"); ptree_trialaverage_helper(2,GLOBAL_PTREE_NLEGS+1,tindex,1,1,0,sindex,dindex,&fnamebase[0],outputname);}}}
    sprintf(fnamebase[0],"ptree_stayon_1_%srecord",gs2); sprintf(fnamebase[1],"ptree_justoff_2_%srecord",gs2);
    for (tindex=1;tindex<=SUITE_TINDEXMAX;tindex++){ for (sindex=0;sindex<=SUITE_SINDEXMAX;sindex++){ for (dindex=0;dindex<=SUITE_DINDEXMAX;dindex++){
      sprintf(outputname,"s1vso2"); ptree_trialaverage_helper(2,GLOBAL_PTREE_NLEGS+1,tindex,1,1,0,sindex,dindex,&fnamebase[0],outputname);}}}
    sprintf(fnamebase[0],"ptree_justoff_1_%srecord",gs2); sprintf(fnamebase[1],"ptree_stayoff_%srecord",gs2);
    for (tindex=1;tindex<=SUITE_TINDEXMAX;tindex++){ for (sindex=0;sindex<=SUITE_SINDEXMAX;sindex++){ for (dindex=0;dindex<=SUITE_DINDEXMAX;dindex++){
      sprintf(outputname,"o1vsoo"); ptree_trialaverage_helper(2,GLOBAL_PTREE_NLEGS+1,tindex,1,1,0,sindex,dindex,&fnamebase[0],outputname);}}}
    sprintf(fnamebase[0],"ptree_justoff_1_%srecord",gs2); sprintf(fnamebase[1],"ptree_juston_2_%srecord",gs2);
    for (tindex=1;tindex<=SUITE_TINDEXMAX;tindex++){ for (sindex=0;sindex<=SUITE_SINDEXMAX;sindex++){ for (dindex=0;dindex<=SUITE_DINDEXMAX;dindex++){
      sprintf(outputname,"o1vsj2"); ptree_trialaverage_helper(2,GLOBAL_PTREE_NLEGS+1,tindex,1,1,0,sindex,dindex,&fnamebase[0],outputname);}}}
    sprintf(fnamebase[0],"ptree_justoff_1_%srecord",gs2); sprintf(fnamebase[1],"ptree_stayon_2_%srecord",gs2);
    for (tindex=1;tindex<=SUITE_TINDEXMAX;tindex++){ for (sindex=0;sindex<=SUITE_SINDEXMAX;sindex++){ for (dindex=0;dindex<=SUITE_DINDEXMAX;dindex++){
      sprintf(outputname,"o1vss2"); ptree_trialaverage_helper(2,GLOBAL_PTREE_NLEGS+1,tindex,1,1,0,sindex,dindex,&fnamebase[0],outputname);}}}
    sprintf(fnamebase[0],"ptree_justoff_1_%srecord",gs2); sprintf(fnamebase[1],"ptree_justoff_2_%srecord",gs2);
    for (tindex=1;tindex<=SUITE_TINDEXMAX;tindex++){ for (sindex=0;sindex<=SUITE_SINDEXMAX;sindex++){ for (dindex=0;dindex<=SUITE_DINDEXMAX;dindex++){
      sprintf(outputname,"o1vso2"); ptree_trialaverage_helper(2,GLOBAL_PTREE_NLEGS+1,tindex,1,1,0,sindex,dindex,&fnamebase[0],outputname);}}}
    sprintf(fnamebase[0],"ptree_stayoff_%srecord",gs2); sprintf(fnamebase[1],"ptree_juston_2_%srecord",gs2);
    for (tindex=1;tindex<=SUITE_TINDEXMAX;tindex++){ for (sindex=0;sindex<=SUITE_SINDEXMAX;sindex++){ for (dindex=0;dindex<=SUITE_DINDEXMAX;dindex++){
      sprintf(outputname,"oovsj2"); ptree_trialaverage_helper(2,GLOBAL_PTREE_NLEGS+1,tindex,1,1,0,sindex,dindex,&fnamebase[0],outputname);}}}
    sprintf(fnamebase[0],"ptree_stayoff_%srecord",gs2); sprintf(fnamebase[1],"ptree_stayon_2_%srecord",gs2);
    for (tindex=1;tindex<=SUITE_TINDEXMAX;tindex++){ for (sindex=0;sindex<=SUITE_SINDEXMAX;sindex++){ for (dindex=0;dindex<=SUITE_DINDEXMAX;dindex++){
      sprintf(outputname,"oovss2"); ptree_trialaverage_helper(2,GLOBAL_PTREE_NLEGS+1,tindex,1,1,0,sindex,dindex,&fnamebase[0],outputname);}}}
    sprintf(fnamebase[0],"ptree_stayoff_%srecord",gs2); sprintf(fnamebase[1],"ptree_justoff_2_%srecord",gs2);
    for (tindex=1;tindex<=SUITE_TINDEXMAX;tindex++){ for (sindex=0;sindex<=SUITE_SINDEXMAX;sindex++){ for (dindex=0;dindex<=SUITE_DINDEXMAX;dindex++){
      sprintf(outputname,"oovso2"); ptree_trialaverage_helper(2,GLOBAL_PTREE_NLEGS+1,tindex,1,1,0,sindex,dindex,&fnamebase[0],outputname);}}}
    sprintf(fnamebase[0],"ptree_juston_2_%srecord",gs2); sprintf(fnamebase[1],"ptree_stayon_2_%srecord",gs2);
    for (tindex=1;tindex<=SUITE_TINDEXMAX;tindex++){ for (sindex=0;sindex<=SUITE_SINDEXMAX;sindex++){ for (dindex=0;dindex<=SUITE_DINDEXMAX;dindex++){
      sprintf(outputname,"j2vss2"); ptree_trialaverage_helper(2,GLOBAL_PTREE_NLEGS+1,tindex,1,1,0,sindex,dindex,&fnamebase[0],outputname);}}}
    sprintf(fnamebase[0],"ptree_juston_2_%srecord",gs2); sprintf(fnamebase[1],"ptree_justoff_2_%srecord",gs2);
    for (tindex=1;tindex<=SUITE_TINDEXMAX;tindex++){ for (sindex=0;sindex<=SUITE_SINDEXMAX;sindex++){ for (dindex=0;dindex<=SUITE_DINDEXMAX;dindex++){
      sprintf(outputname,"j2vso2"); ptree_trialaverage_helper(2,GLOBAL_PTREE_NLEGS+1,tindex,1,1,0,sindex,dindex,&fnamebase[0],outputname);}}}
    sprintf(fnamebase[0],"ptree_stayon_2_%srecord",gs2); sprintf(fnamebase[1],"ptree_justoff_2_%srecord",gs2);
    for (tindex=1;tindex<=SUITE_TINDEXMAX;tindex++){ for (sindex=0;sindex<=SUITE_SINDEXMAX;sindex++){ for (dindex=0;dindex<=SUITE_DINDEXMAX;dindex++){
      sprintf(outputname,"s2vso2"); ptree_trialaverage_helper(2,GLOBAL_PTREE_NLEGS+1,tindex,1,1,0,sindex,dindex,&fnamebase[0],outputname);}}}
    sprintf(fnamebase[0],"ptree_juston_1_%srecord",gs2);
    sprintf(fnamebase[1],"ptree_stayon_1_%srecord",gs2);
    sprintf(fnamebase[2],"ptree_justoff_1_%srecord",gs2);
    sprintf(fnamebase[3],"ptree_juston_2_%srecord",gs2);
    sprintf(fnamebase[4],"ptree_stayon_2_%srecord",gs2);
    sprintf(fnamebase[5],"ptree_justoff_2_%srecord",gs2);
    sprintf(fnamebase[6],"ptree_stayoff_%srecord",gs2);
    for (tindex=1;tindex<=SUITE_TINDEXMAX;tindex++){ for (sindex=0;sindex<=SUITE_SINDEXMAX;sindex++){ for (dindex=0;dindex<=SUITE_DINDEXMAX;dindex++){
      sprintf(outputname,"j1s1o1j2s2o2oo");
      ptree_trialaverage_helper(7,GLOBAL_PTREE_NLEGS+1,tindex,1,1,0,sindex,dindex,&fnamebase[0],outputname);}}}
    tfree(fnamebase[0]);fnamebase[0]=NULL;
    tfree(fnamebase[1]);fnamebase[1]=NULL;
    tfree(fnamebase[2]);fnamebase[3]=NULL;
    tfree(fnamebase[3]);fnamebase[3]=NULL;
    tfree(fnamebase[4]);fnamebase[4]=NULL;
    tfree(fnamebase[5]);fnamebase[5]=NULL;
    tfree(fnamebase[6]);fnamebase[6]=NULL;
    tfree(fnamebase);fnamebase=NULL;
    if (SUITE_BITBYBIT_REMOVE){ 
      bitbybit=1;continue_flag=1;
      do{ 
	sprintf(outputname,"./ptree_juston_1_%srecord_%d",gs2,bitbybit); continue_flag=checktofind(outputname);
	if (continue_flag){ sprintf(outputname,"""rm"" ./ptree_juston_1_%srecord_%d",gs2,bitbybit); system(outputname); bitbybit += 1;}}
      while (continue_flag);
      bitbybit=1;continue_flag=1;
      do{ 
	sprintf(outputname,"./ptree_stayon_1_%srecord_%d",gs2,bitbybit); continue_flag=checktofind(outputname);
	if (continue_flag){ sprintf(outputname,"""rm"" ./ptree_stayon_1_%srecord_%d",gs2,bitbybit); system(outputname); bitbybit += 1;}}
      while (continue_flag);
      bitbybit=1;continue_flag=1;
      do{ 
	sprintf(outputname,"./ptree_justoff_1_%srecord_%d",gs2,bitbybit); continue_flag=checktofind(outputname);
	if (continue_flag){ sprintf(outputname,"""rm"" ./ptree_justoff_1_%srecord_%d",gs2,bitbybit); system(outputname); bitbybit += 1;}}
      while (continue_flag);
      bitbybit=1;continue_flag=1;
      do{ 
	sprintf(outputname,"./ptree_juston_2_%srecord_%d",gs2,bitbybit); continue_flag=checktofind(outputname);
	if (continue_flag){ sprintf(outputname,"""rm"" ./ptree_juston_2_%srecord_%d",gs2,bitbybit); system(outputname); bitbybit += 1;}}
      while (continue_flag);
      bitbybit=1;continue_flag=1;
      do{ 
	sprintf(outputname,"./ptree_stayon_2_%srecord_%d",gs2,bitbybit); continue_flag=checktofind(outputname);
	if (continue_flag){ sprintf(outputname,"""rm"" ./ptree_stayon_2_%srecord_%d",gs2,bitbybit); system(outputname); bitbybit += 1;}}
      while (continue_flag);
      bitbybit=1;continue_flag=1;
      do{ 
	sprintf(outputname,"./ptree_justoff_2_%srecord_%d",gs2,bitbybit); continue_flag=checktofind(outputname);
	if (continue_flag){ sprintf(outputname,"""rm"" ./ptree_justoff_2_%srecord_%d",gs2,bitbybit); system(outputname); bitbybit += 1;}}
      while (continue_flag);
      bitbybit=1;continue_flag=1;
      do{ 
	sprintf(outputname,"./ptree_stayoff_%srecord_%d",gs2,bitbybit); continue_flag=checktofind(outputname);
	if (continue_flag){ sprintf(outputname,"""rm"" ./ptree_stayoff_%srecord_%d",gs2,bitbybit); system(outputname); bitbybit += 1;}}
      while (continue_flag);}}
  else /* if (!GLOBAL_PTREE_BITBYBIT) */{
    ptreerate(p->pjuston_1); ptreedump_starter(p->pjuston_1,"j1",2,0,0,0,+1,-1);
    ptreerate(p->pstayon_1); ptreedump_starter(p->pstayon_1,"s1",2,0,0,0,+1,-1);
    ptreerate(p->pjustoff_1); ptreedump_starter(p->pjustoff_1,"o1",2,0,0,0,+1,-1);
    ptreerate(p->pjuston_2); ptreedump_starter(p->pjuston_2,"j2",2,0,0,0,+1,-1);
    ptreerate(p->pstayon_2); ptreedump_starter(p->pstayon_2,"s2",2,0,0,0,+1,-1);
    ptreerate(p->pjustoff_2); ptreedump_starter(p->pjustoff_2,"o2",2,0,0,0,+1,-1);
    ptreerate(p->pstayoff); ptreedump_starter(p->pstayoff,"oo",2,0,0,0,+1,-1);}
}

/* Here are the lyapunov functions */

void nrasave(struct neuronarray *Nra,char *fgvn)
{
  /* does a bare bones save state of neuronarray,
     assuming that architecture doesn't change.
     intended to be reloaded with nraload 
     assumes fgvn does NOT start with "./" */
  char filename[256];
  FILE *fp=NULL;
  int nr=0,nc=0;
  struct neuron *n=NULL;
  sprintf(filename,"./%s",fgvn);
  if ((fp=fopen(filename,"w"))==NULL){ printf(" %% warning, couldn't open %s in nrasave\n",filename); fp=stdout;}
  fwrite(&(GLOBAL_time),sizeof(double),1,fp);
  fwrite(Nra->Vra,sizeof(double),Nra->rows*Nra->cols,fp);
  fwrite(Nra->sAra,sizeof(double),Nra->rows*Nra->cols,fp);
  fwrite(Nra->sNra,sizeof(double),Nra->rows*Nra->cols,fp);
  fwrite(Nra->sGra,sizeof(double),Nra->rows*Nra->cols,fp);
  fwrite(Nra->VSra,sizeof(double),Nra->rows*Nra->cols,fp);
  for (nr=0;nr<Nra->rows;nr++){ for (nc=0;nc<Nra->cols;nc++){
    n=nget(Nra,nr,nc);
    fwrite(&(n->spikeinput_flag),sizeof(int),1,fp);
    fwrite(&(n->spikeinput_multiplicity),sizeof(int),1,fp);
    fwrite(&(n->spikeinput_rseed),sizeof(long long int),1,fp);
    fwrite(&(n->spikeinput_time),sizeof(double),1,fp);
    fwrite(&(n->spikelast),sizeof(double),1,fp);
    fwrite(&(n->spiketime),sizeof(double),1,fp);
    fwrite(&(n->spikenext),sizeof(double),1,fp);}}
  if (fp!=stdout){ fclose(fp);}
};

void nraload(struct neuronarray *Nra,char *fgvn)
{
  /* this does a bare bones load state 
     assuming state was written with nrasave,
     assumes fgvn does NOT start with "./" */
  double told=0;
  char filename[256];
  FILE *fp=NULL;
  int nr=0,nc=0;
  struct neuron *n=NULL;
  double temp=0;
  sprintf(filename,"./%s",fgvn);
  if ((fp=fopen(filename,"r"))==NULL){ printf(" %% warning, couldn't open %s in nrasave\n",filename); fp=stdout;}
  fread(&(told),sizeof(double),1,fp); told = GLOBAL_time-told;
  fread(Nra->Vra,sizeof(double),Nra->rows*Nra->cols,fp);
  fread(Nra->sAra,sizeof(double),Nra->rows*Nra->cols,fp);
  fread(Nra->sNra,sizeof(double),Nra->rows*Nra->cols,fp);
  fread(Nra->sGra,sizeof(double),Nra->rows*Nra->cols,fp);
  fread(Nra->VSra,sizeof(double),Nra->rows*Nra->cols,fp);
  for (nr=0;nr<Nra->rows;nr++){ for (nc=0;nc<Nra->cols;nc++){
    n=nget(Nra,nr,nc);
    fread(&(n->spikeinput_flag),sizeof(int),1,fp);
    fread(&(n->spikeinput_multiplicity),sizeof(int),1,fp);
    fread(&(n->spikeinput_rseed),sizeof(long long int),1,fp);
    fread(&(temp),sizeof(double),1,fp); n->spikeinput_time = temp+told;
    fread(&(temp),sizeof(double),1,fp); n->spikelast = temp+told;
    fread(&(temp),sizeof(double),1,fp); n->spiketime = temp+told;
    fread(&(temp),sizeof(double),1,fp); n->spikenext = temp+told;}}
  if (fp!=stdout){ fclose(fp);}
}

double nraperturb(struct neuronarray *Nra,double jiggle)
{
  /* returns the norm of the perturbation */
  int voltage_jiggle=0;
  int nr=0,nc=0,tab=0;
  double temp=0,temp2=0;
  double vt=EIFvsIF ? VOLTAGE_THRESHOLD_EIF : VOLTAGE_THRESHOLD;
  for (nr=0;nr<Nra->rows;nr++){ for (nc=0;nc<Nra->cols;nc++){
    tab = nr + nc*Nra->rows;
    /* don't perturb voltage */
    if (voltage_jiggle){ temp=jiggle*2*(rand01-0.5)*(vt-Nra->Vra[tab]); temp2 += pow(temp,2); Nra->Vra[tab] += temp;}
    temp=exp(jiggle*2*(rand01-0.5)); temp2 += pow(temp*Nra->sAra[tab],2); Nra->sAra[tab] *= temp;
    temp=exp(jiggle*2*(rand01-0.5)); temp2 += pow(temp*Nra->sNra[tab],2); Nra->sNra[tab] *= temp;
    temp=exp(jiggle*2*(rand01-0.5)); temp2 += pow(temp*Nra->sGra[tab],2); Nra->sGra[tab] *= temp;
    slavecalc(Nra->sAra[tab],Nra->sNra[tab],Nra->sGra[tab],NULL,&(Nra->VSra[tab]),NULL);}}
  return sqrt(temp2)/(double)(Nra->rows*Nra->cols);
}

struct lyapunov * lyapunovmake(struct neuronarray *Nra,double update_every,double jiggle)
{
  /* jiggle should be between 0 and 1 */
  struct lyapunov *y=NULL;
  char fgvn[128];
  y = (struct lyapunov *) tmalloc(sizeof(struct lyapunov));
  y->Nra = Nra;
  y->total_time = 0;
  y->update_last = GLOBAL_TI;
  y->update_every = maximum(1,update_every);
  y->original_vs_perturbed = +1;
  y->perturbnorm = 0;
  y->jiggle = jiggle;
  sprintf(fgvn,"nra_%srecord_perturbed_save",GLOBAL_STRING_2); nrasave(y->Nra,fgvn);
  y->L = llistmake(); /* llist of lyapunov exponents determined by measuring conductance profiles */
  y->pL = llistmake(); /* llist of lyapunov exponents determined by measuring best pcs */
  return y;
}

void lyapunovtfree(struct lyapunov *y)
{
  char command[512];
  llisttfree2(y->L); y->L=NULL;
  llisttfree2(y->pL); y->pL=NULL;
  sprintf(command,"""rm"" ./nra_%srecord_perturbed_save",GLOBAL_STRING_2); system(command);
  sprintf(command,"""rm"" ./nra_%srecord_original_save",GLOBAL_STRING_2); system(command);
  tfree(y);
}

void lyapunovupdate(struct lyapunov *y,double t,double DT)
{
  /* runs alternate universe, renormalizing alternate trajectory to maximize unstable direction */
  int verbose=0;
  char fgvn[128];
  struct neuronarray *Nra=y->Nra;
  int nr=0,nc=0,tab=0,tab2=0;
  double *sveco=NULL,*svecp=NULL,*temp=NULL,tempp=0;
  if (verbose){ printf(" %% [entering lyapunovupdate] with t=%0.2f DT=%0.2f\n",t,DT);}
  y->total_time += DT;
  if (t >= y->update_last + y->update_every){
    if (verbose){ printf(" %% time to act\n");}
    if (y->original_vs_perturbed == +1){ /* original run, save reload and perturb*/
      if (verbose){ printf(" %% original run, saving...\n");}
      sprintf(fgvn,"nra_%srecord_original_save",GLOBAL_STRING_2); nrasave(Nra,fgvn);
      if (verbose){ printf(" %% noting state vector...\n");}
      sveco = (double *) tcalloc(Nra->rows*Nra->cols*3,sizeof(double));
      for (nr=0;nr<Nra->rows;nr++){ for (nc=0;nc<Nra->cols;nc++){
	tab = nr + nc*Nra->rows; 
	sveco[tab + 0*Nra->rows*Nra->cols] = Nra->sAra[tab];
	sveco[tab + 1*Nra->rows*Nra->cols] = Nra->sNra[tab];
	sveco[tab + 2*Nra->rows*Nra->cols] = Nra->sGra[tab];}}
      if (verbose){ printf(" %% and reload...\n");}
      sprintf(fgvn,"nra_%srecord_perturbed_save",GLOBAL_STRING_2); nraload(Nra,fgvn); 
      if (verbose){ printf(" %% noting state vector...\n");}
      svecp = (double *) tcalloc(Nra->rows*Nra->cols*3,sizeof(double));
      for (nr=0;nr<Nra->rows;nr++){ for (nc=0;nc<Nra->cols;nc++){
	tab = nr + nc*Nra->rows; 
	svecp[tab + 0*Nra->rows*Nra->cols] = Nra->sAra[tab];
	svecp[tab + 1*Nra->rows*Nra->cols] = Nra->sNra[tab];
	svecp[tab + 2*Nra->rows*Nra->cols] = Nra->sGra[tab];}}
      tempp=0; for (nr=0;nr<Nra->rows*Nra->cols*3;nr++){ tempp+=pow(svecp[nr]-sveco[nr],2);} 
      tempp = sqrt(tempp/(double)(Nra->rows*Nra->cols));
      temp = (double *) tmalloc(sizeof(double));
      *temp = log(tempp/(y->perturbnorm<=0 ? 1 : y->perturbnorm))/(double)y->update_every; if (!finite(*temp)){ *temp=0;}
      if (verbose){ printf(" %% found difference in norms %f, perturbnorm %f, log(d/p)/t %f\n",tempp,y->perturbnorm,*temp);}
      if (verbose){ printf(" %% adding to y->L\n");}
      litemadd(y->L,temp);
      if (tempp==0){ if (verbose){ printf(" %% no difference, must rejiggle\n");} y->perturbnorm = nraperturb(Nra,y->jiggle);}
      else /* if (tempp>0) */{ 
	if (verbose){ printf(" %% renormalize...\n");}
	y->perturbnorm=0;
	for (nr=0;nr<Nra->rows;nr++){ for (nc=0;nc<Nra->cols;nc++){
	  tab = nr + nc*Nra->rows; 
	  tab2 = tab + 0*Nra->rows*Nra->cols; Nra->sAra[tab] = maximum(0,sveco[tab2] + (svecp[tab2]-sveco[tab2])*y->jiggle/tempp);
	  y->perturbnorm += pow(Nra->sAra[tab]-sveco[tab2],2);
	  tab2 = tab + 1*Nra->rows*Nra->cols; Nra->sNra[tab] = maximum(0,sveco[tab2] + (svecp[tab2]-sveco[tab2])*y->jiggle/tempp);
	  y->perturbnorm += pow(Nra->sNra[tab]-sveco[tab2],2);
	  tab2 = tab + 2*Nra->rows*Nra->cols; Nra->sGra[tab] = maximum(0,sveco[tab2] + (svecp[tab2]-sveco[tab2])*y->jiggle/tempp);
	  y->perturbnorm += pow(Nra->sGra[tab]-sveco[tab2],2);
	  slavecalc(Nra->sAra[tab],Nra->sNra[tab],Nra->sGra[tab],NULL,&(Nra->VSra[tab]),NULL);}}
	y->perturbnorm = sqrt(y->perturbnorm/(double)(Nra->rows*Nra->cols));}
      tfree(sveco);sveco=NULL;tfree(svecp);svecp=NULL;
      y->original_vs_perturbed *= -1;}
    else if (y->original_vs_perturbed == -1){ /* perturbed run, save and reload */
      if (verbose){ printf(" %% perturbed run, saving perturbed state\n");}
      sprintf(fgvn,"nra_%srecord_perturbed_save",GLOBAL_STRING_2); nrasave(Nra,fgvn);      
      if (verbose){ printf(" %% perturbed run, open original save\n");}
      sprintf(fgvn,"nra_%srecord_original_save",GLOBAL_STRING_2); nraload(Nra,fgvn);
      y->original_vs_perturbed *= -1;}
    y->update_last = maximum(t,y->update_last + y->update_every);}
}

void lyapunovupdate_old(struct lyapunov *y,double t,double DT)
{
  /* runs an alternate universe of perturbed states which are randomly jiggled to shadow original trajectory */
  int verbose=0;
  char fgvn[128],command[512];
  FILE *fp=NULL;
  struct neuronarray *Nra=y->Nra;
  struct neuron *n=NULL;
  double *ra=NULL;
  int nr=0,nc=0,tab=0;
  double *temp=NULL,temp2=0;
  struct strobetrace *st=NULL;
  int rows=0,cols=0,s=0,nr2=0,nc2=0,na=0;
  double *tempVSnow=NULL,*tempVSold=NULL;
  double *tempVSnow2=NULL,*tempVSold2=NULL;
  double *pcsranow=NULL,*pcsraold=NULL;
  double tempp=0;
  int *maxindexnow=NULL,nmaxnow=0,*minindexnow=NULL,nminnow=0;
  int *maxindexold=NULL,nmaxold=0,*minindexold=NULL,nminold=0;
  if (verbose){ printf(" %% [entering lyapunovupdate] with t=%0.2f DT=%0.2f\n",t,DT);}
  y->total_time += DT;
  if (t >= y->update_last + y->update_every){
    if (verbose){ printf(" %% time to act\n");}
    if (y->original_vs_perturbed == +1){ /* original run, simply save */
      if (verbose){ printf(" %% original run, simply save...\n");}
      sprintf(fgvn,"nra_%srecord_original_save",GLOBAL_STRING_2); nrasave(Nra,fgvn);
      if (verbose){ printf(" %% and reload...\n");}
      sprintf(fgvn,"nra_%srecord_perturbed_save",GLOBAL_STRING_2); nraload(Nra,fgvn); 
      if (verbose){ printf(" %% and perturb.\n");}
      y->perturbnorm = nraperturb(Nra,y->jiggle);
      y->original_vs_perturbed *= -1;}
    else if (y->original_vs_perturbed == -1){ /* perturbed run, compare norms and reload */
      if (verbose){ printf(" %% perturbed run, open original save\n");}
      sprintf(fgvn,"nra_%srecord_original_save",GLOBAL_STRING_2);
      temp2=0;
      if ((fp=fopen(fgvn,"r"))!=NULL){ 
	ra = (double *) tmalloc(Nra->rows*Nra->cols*5*sizeof(double));
	fread(ra,sizeof(double),1,fp);
	fread(ra,sizeof(double),Nra->rows*Nra->cols*5,fp);
	for (nr=0;nr<Nra->rows;nr++){ for (nc=0;nc<Nra->cols;nc++){
	  /* don't measure voltage */
	  tab = nr + nc*Nra->rows + 1*Nra->rows*Nra->cols;
	  temp2 += pow(Nra->sAra[nr + nc*Nra->rows]-ra[tab],2);
	  tab = nr + nc*Nra->rows + 2*Nra->rows*Nra->cols;
	  temp2 += pow(Nra->sNra[nr + nc*Nra->rows]-ra[tab],2);
	  tab = nr + nc*Nra->rows + 3*Nra->rows*Nra->cols;
	  temp2 += pow(Nra->sGra[nr + nc*Nra->rows]-ra[tab],2);
	  tab = nr + nc*Nra->rows + 4*Nra->rows*Nra->cols;
	  temp2 += pow(Nra->VSra[nr + nc*Nra->rows]-ra[tab],2);}}
	if (STROBETRACE_BOTHER){
	  if (verbose){ printf(" %% \n");}
	  if (verbose){ printf(" %% also computing most similar orientation\n");}
	  st = GLOBAL_STROBETRACE;
	  rows = st->pcspiesacross*PIE_ROW_DIA;
	  cols = st->pcspiestall*PIE_COL_DIA;
	  s = maximum(0,GLOBAL_SPACE_SMOOTHER);
	  tempVSnow = (double *) tcalloc((rows+2*s)*(cols+2*s),sizeof(double));  
	  tempVSold = (double *) tcalloc((rows+2*s)*(cols+2*s),sizeof(double));  
	  for (nr=0;nr<rows+2*s;nr++){ for (nc=0;nc<cols+2*s;nc++){
	    nr2 = periodize(st->n->row+nr-rows/2-s,0,Nra->rows);
	    nc2 = periodize(st->n->col+nc-cols/2-s,0,Nra->cols);
	    n = nget(Nra,nr2,nc2);
	    tempVSnow[nr+nc*(rows+2*s)] = *(n->VS);
	    tempVSold[nr+nc*(rows+2*s)] = ra[nr2 + nc2*Nra->rows + 4*Nra->rows*Nra->cols];}}
	  tempVSnow2 = spacesmear(tempVSnow,rows+2*s,cols+2*s,s);
	  tempVSold2 = spacesmear(tempVSold,rows+2*s,cols+2*s,s);
	  tfree(tempVSnow);tempVSnow=NULL;tfree(tempVSold);tempVSold=NULL;
	  tempVSnow = (double *) tcalloc(rows*cols,sizeof(double));
	  tempVSold = (double *) tcalloc(rows*cols,sizeof(double));
	  for (nr=0;nr<rows;nr++){ for (nc=0;nc<cols;nc++){
	    tempVSnow[nr + nc*rows] = tempVSnow2[(nr+s) + (nc+s)*(rows+2*s)];
	    tempVSold[nr + nc*rows] = tempVSold2[(nr+s) + (nc+s)*(rows+2*s)];}}
	  tfree(tempVSnow2);tempVSnow2=NULL;tfree(tempVSold2);tempVSold2=NULL;
	  pcsranow = (double *) tcalloc(st->nangles,sizeof(double));
	  pcsraold = (double *) tcalloc(st->nangles,sizeof(double));
	  for (na=0;na<st->nangles;na++){
	    pcsranow[na] = correlation(tempVSnow,&(st->pcsra[0+0*rows+na*rows*cols]),rows*cols);
	    pcsraold[na] = correlation(tempVSold,&(st->pcsra[0+0*rows+na*rows*cols]),rows*cols);}
	  tfree(tempVSnow);tempVSnow=NULL; tfree(tempVSold);tempVSold=NULL;
	  if (verbose){ printf(" %% obtained pcsranow:\n"); raprintf(pcsranow,"double",1,st->nangles,"pcsranow");}
	  if (verbose){ printf(" %% obtained pcsraold:\n"); raprintf(pcsraold,"double",1,st->nangles,"pcsranow");}
	  maxmindex("double",pcsranow,st->nangles,&maxindexnow,&nmaxnow,&minindexnow,&nminnow,0.000001);
	  maxmindex("double",pcsraold,st->nangles,&maxindexold,&nmaxold,&minindexold,&nminold,0.000001);
	  tempp = 0;
	  tempp += pow(maxindexnow[periodize((int)rand01*nmaxnow,0,nmaxnow)] - maxindexold[periodize((int)rand01*nmaxold,0,nmaxold)],2);
	  tempp += pow(minindexnow[periodize((int)rand01*nminnow,0,nminnow)] - minindexold[periodize((int)rand01*nminold,0,nminold)],2);
	  tempp /= pow(st->nangles,2);
	  tfree(maxindexnow);maxindexnow=NULL;tfree(minindexnow);minindexnow=NULL;
	  tfree(maxindexold);maxindexold=NULL;tfree(minindexold);minindexold=NULL;
	  tfree(pcsranow);pcsranow=NULL; tfree(pcsraold);pcsraold=NULL;}
	tfree(ra); ra=NULL;
	fclose(fp);}
      temp = (double *) tmalloc(sizeof(double));
      *temp = log(temp2/(y->perturbnorm<=0 ? 1 : y->perturbnorm))/(double)y->update_every; if (!finite(*temp)){ *temp=0;}
      if (verbose){ printf(" %% found difference in norms %f, perturbnorm %f, log(d/p)/t %f\n",temp2,y->perturbnorm,*temp);}
      if (verbose){ printf(" %% adding to y->L\n");}
      litemadd(y->L,temp);
      if (STROBETRACE_BOTHER){
	temp = (double *) tmalloc(sizeof(double));
	*temp = log(tempp/(y->perturbnorm<=0 ? 1 : y->perturbnorm))/(double)y->update_every; if (!finite(*temp)){ *temp=0;}
	if (verbose){ printf(" %% found difference in pnorms %f, perturbnorm %f, log(d/p)/t %f\n",tempp,y->perturbnorm,*temp);}
	if (verbose){ printf(" %% adding to y->pL\n");}
	litemadd(y->pL,temp);}
      if (verbose){ printf(" %% loading original save\n");}
      sprintf(fgvn,"nra_%srecord_original_save",GLOBAL_STRING_2); nraload(Nra,fgvn);
      if (verbose){ printf(" %% shuffling save files\n");}
      sprintf(command,"mv ./nra_%srecord_original_save ./nra_%srecord_perturbed_save;",GLOBAL_STRING_2,GLOBAL_STRING_2); system(command);
      y->original_vs_perturbed *= -1;}
    y->update_last = maximum(t,y->update_last + y->update_every);}
}

void lyapunovdump(struct lyapunov *y,int dump_type)
{
  char filename[256],gs2[512];
  double *ra=NULL;
  int nb=0;
  FILE *fp=NULL;
  struct litem *l=NULL;
  sprintf(gs2,"%srecord",GLOBAL_STRING_2);
  if (dump_type==0){
    llistsort(y->L->first,y->L->last,y->L->length,&double_compare);
    sprintf(filename,"./lyapunov_%s.m",gs2);
    if ((fp=fopen(filename,"w"))!=NULL){
      l=y->L->first; nb=0;
      while (l!=NULL && nb<y->L->length){ 
	fprintf(fp,"ly_%s_L(%d) = %f;\n",gs2,nb+1,*(double *)l->item);
	nb++; l=l->child;}
      fprintf(fp,"hist(ly_%s_L,32);\n",gs2);
      if (STROBETRACE_BOTHER){
	llistsort(y->pL->first,y->pL->last,y->pL->length,&double_compare);
	l=y->pL->first; nb=0;
	while (l!=NULL && nb<y->pL->length){ 
	  fprintf(fp,"ly_%s_pL(%d) = %f;\n",gs2,nb+1,*(double *)l->item);
	  nb++; l=l->child;}
	fprintf(fp,"hist(ly_%s_pL,32);\n",gs2);}
      fclose(fp);fp=NULL;}}
  else if (dump_type==1){
    llistsort(y->L->first,y->L->last,y->L->length,&double_compare);
    ra = (double *) tcalloc(y->L->length,sizeof(double));
    l=y->L->first; nb=0;
    while (l!=NULL && nb<y->L->length){ ra[nb]=*(double *)l->item; nb++; l=l->child;}
    sprintf(filename,"./lyapunov_%s_L",gs2);
    ra2jpg(ra,"double",1,y->L->length,0,filename,0);
    tfree(ra);
    if (STROBETRACE_BOTHER){
      llistsort(y->pL->first,y->pL->last,y->pL->length,&double_compare);
      ra = (double *) tcalloc(y->pL->length,sizeof(double));
      l=y->pL->first; nb=0;
      while (l!=NULL && nb<y->pL->length){ ra[nb]=*(double *)l->item; nb++; l=l->child;}
      sprintf(filename,"./lyapunov_%s_pL",gs2);
      ra2jpg(ra,"double",1,y->pL->length,0,filename,0);
      tfree(ra);}}
}

/* Here are the avalanche functions */

struct avalanche * avalanchemake(struct neuronarray *Nra,double update_every,int logdmin,int logdmax)
{
  int nd=0,d=0,nt2s=0;
  struct avalanche *a=NULL;
  a = (struct avalanche *) tmalloc(sizeof(struct avalanche));
  a->Nra = Nra;
  a->update_last = GLOBAL_TI;
  a->update_every = update_every;
  a->logdmin = logdmin;
  a->logdmax = logdmax;
  a->mraprev = (double **) tcalloc(a->logdmax*5,sizeof(double *));
  for (nd=a->logdmin;nd<a->logdmax;nd++){ d=(int)pow(2,nd); for (nt2s=0;nt2s<5;nt2s++){
    a->mraprev[nd+nt2s*a->logdmax] = (double *) tcalloc((Nra->rows/d+1)*(Nra->cols/d+1),sizeof(double));}}
  a->mranext = (double **) tcalloc(a->logdmax*5,sizeof(double *));
  for (nd=a->logdmin;nd<a->logdmax;nd++){ d=(int)pow(2,nd); for (nt2s=0;nt2s<5;nt2s++){
    a->mranext[nd+nt2s*a->logdmax] = (double *) tcalloc((Nra->rows/d+1)*(Nra->cols/d+1),sizeof(double));}}
  a->mhist = (double **) tcalloc(a->logdmax*5,sizeof(double *));
  for (nd=a->logdmin;nd<a->logdmax;nd++){ d=(int)pow(2,nd); for (nt2s=0;nt2s<5;nt2s++){
    a->mhist[nd+nt2s*a->logdmax] = (double *) tcalloc((Nra->rows/d+1)*(Nra->cols/d+1),sizeof(double));}}
  return a;
}

void avalanchetfree(struct avalanche *a)
{
  int nd=0,d=-0,nt2s=0;
  for (nd=a->logdmin;nd<a->logdmax;nd++){ d=(int)pow(2,nd); for (nt2s=0;nt2s<5;nt2s++){
    tfree(a->mraprev[nd+nt2s*a->logdmax]); a->mraprev[nd+nt2s*a->logdmax]=NULL;}}
  tfree(a->mraprev); a->mraprev=NULL;
  for (nd=a->logdmin;nd<a->logdmax;nd++){ d=(int)pow(2,nd); for (nt2s=0;nt2s<5;nt2s++){
    tfree(a->mranext[nd+nt2s*a->logdmax]); a->mranext[nd+nt2s*a->logdmax]=NULL;}}
  tfree(a->mranext); a->mranext=NULL;
  for (nd=a->logdmin;nd<a->logdmax;nd++){ d=(int)pow(2,nd); for (nt2s=0;nt2s<5;nt2s++){
    tfree(a->mhist[nd+nt2s*a->logdmax]); a->mhist[nd+nt2s*a->logdmax]=NULL;}}
  tfree(a->mhist); a->mhist=NULL;
  tfree(a); a=NULL;
}

void avalancheupdate(struct avalanche *a,double t,double DT)
{
  int verbose=0;
  int nr=0,nr2=0,nc=0,nc2=0,nd=0,d=0,nt2s=0,tab=0;
  struct neuronarray *Nra=a->Nra;
  struct neuron *n=NULL;
  if (verbose){ printf(" %% [entering avalancheupdate] with t=%f DT=%f update_last=%f update_next=%f\n",t,DT,a->update_last,a->update_last+a->update_every);}
  for (nr=0;nr<Nra->rows;nr++){ for (nc=0;nc<Nra->cols;nc++){
    n = nget(Nra,nr,nc);
    if (n->spikelast==n->spiketime && n->spiketime>=t && n->spiketime <= t+DT){
      if (verbose){ printf(" %% neuron (%d,%d) of type %d spiked at time %f\n",n->row,n->col,(int)*(n->t2s),n->spiketime);}
      if (n->spiketime < a->update_last + a->update_every){
	if (verbose){ printf(" %% current stage\n");}
	for (nd=a->logdmin;nd<a->logdmax;nd++){ 
	  d=(int)pow(2,nd); nr2=nr/d; nc2=nc/d;
	  if (verbose){ printf(" %% (%d,%d) level %d d %d --> (%d,%d)\n",nr,nc,nd,d,nr2,nc2);}
	  a->mraprev[nd + (int)*(n->t2s)*a->logdmax][nr2 + nc2*(Nra->rows/d+1)] += 1;
	  a->mraprev[nd + 4*a->logdmax][nr2 + nc2*(Nra->rows/d+1)] += 1;}}
      else if (n->spiketime >=a->update_last + a->update_every){
	if (verbose){ printf(" %% next stage\n");}
	for (nd=a->logdmin;nd<a->logdmax;nd++){ 
	  d=(int)pow(2,nd); nr2=nr/d; nc2=nc/d; 
	  if (verbose){ printf(" %% (%d,%d) level %d d %d --> (%d,%d)\n",nr,nc,nd,d,nr2,nc2);}
	  a->mranext[nd + (int)*(n->t2s)*a->logdmax][nr2 + nc2*(Nra->rows/d+1)] += 1;
	  a->mranext[nd + 4*a->logdmax][nr2 + nc2*(Nra->rows/d+1)] += 1;}}}}}
  if (verbose>1){ 
    printf(" %% mraprev\n");
    for (nd=a->logdmin;nd<a->logdmax;nd++){ 
      d=(int)pow(2,nd); printf(" %% level %d, ",nd);
      for (nt2s=0;nt2s<5;nt2s++){ 
	printf(" nt2s %d\n",nt2s);
	for (nr=0;nr<(Nra->rows/d+1);nr++){ 
	  for (nc=0;nc<(Nra->cols/d+1);nc++){ 
	    printf("%d ",(int)a->mraprev[nd+nt2s*a->logdmax][nr+nc*(Nra->rows/d+1)]);}
	  printf("\n");}}}
    printf(" %% mranext\n");
    for (nd=a->logdmin;nd<a->logdmax;nd++){ 
      d=(int)pow(2,nd); printf(" %% level %d, ",nd);
      for (nt2s=0;nt2s<5;nt2s++){ 
	printf(" nt2s %d\n",nt2s);
	for (nr=0;nr<(Nra->rows/d+1);nr++){ 
	  for (nc=0;nc<(Nra->cols/d+1);nc++){ 
	    printf("%d ",(int)a->mranext[nd+nt2s*a->logdmax][nr+nc*(Nra->rows/d+1)]);}
	  printf("\n");}}}}
  if (t+DT > a->update_last + a->update_every){
    if (verbose){ printf(" %% need to update, adding to mhist\n");}
    for (nd=a->logdmin;nd<a->logdmax;nd++){ d=(int)pow(2,nd); for (nt2s=0;nt2s<5;nt2s++){
      tab=0;
      for (nr=0;nr<(Nra->rows/d+1);nr++){ for (nc=0;nc<(Nra->cols/d+1);nc++){
	tab += (a->mraprev[nd+nt2s*a->logdmax][nr + nc*(Nra->rows/d+1)] > 0);}}
      a->mhist[nd+nt2s*a->logdmax][tab] += 1;}}
    if (verbose>1){
      printf(" %% mhist\n");
      for (nd=a->logdmin;nd<a->logdmax;nd++){ 
	d=(int)pow(2,nd); printf(" %% level %d, ",nd);
	for (nt2s=0;nt2s<5;nt2s++){ 
	  printf(" nt2s %d\n",nt2s);
	  for (nr=0;nr<(Nra->rows/d+1)*(Nra->cols/d+1);nr++){ 
	    printf("%d ",(int)a->mhist[nd+nt2s*a->logdmax][nr]);}
	  printf("\n");}}}
    if (verbose){ printf(" %% freeing mraprev\n");}
    for (nd=a->logdmin;nd<a->logdmax;nd++){ d=(int)pow(2,nd); for (nt2s=0;nt2s<5;nt2s++){
      tfree(a->mraprev[nd+nt2s*a->logdmax]); a->mraprev[nd+nt2s*a->logdmax]=NULL;}}
    tfree(a->mraprev); a->mraprev=NULL;
    if (verbose){ printf(" %% swapping mraprev with mranext\n");}
    a->mraprev=a->mranext;
    if (verbose){ printf(" %% initializing mranext\n");}
    a->mranext = (double **) tcalloc(a->logdmax*5,sizeof(double *));
    for (nd=a->logdmin;nd<a->logdmax;nd++){ d=(int)pow(2,nd); for (nt2s=0;nt2s<5;nt2s++){
      a->mranext[nd+nt2s*a->logdmax] = (double *) tcalloc((Nra->rows/d+1)*(Nra->cols/d+1),sizeof(double));}}
    a->update_last = maximum(t+DT,a->update_last + a->update_every);}
}

void avalanchedump(struct avalanche *a,int dump_type,char *fgvn)
{
  /* assumes fgvn does NOT start with "./" */
  char filename[512];
  int nd=0,d=0,nt2s=0,nr=0;
  FILE *fp=NULL;
  if (dump_type==0){
    sprintf(filename,"./%s.m",fgvn);
    if ((fp=fopen(filename,"w"))!=NULL){
      fprintf(fp,";\n");
      fprintf(fp,"%% %s{a}(b,c) --> level a, histbin b, t2s c;\n",fgvn);
      for (nd=a->logdmin;nd<a->logdmax;nd++){ d=(int)pow(2,nd); for (nt2s=0;nt2s<5;nt2s++){
	for (nr=0;nr<(a->Nra->rows/d+1)*(a->Nra->cols/d+1);nr++){
	  fprintf(fp,"%s{%d}(%d,%d) = %0.16f;\n",fgvn,nd+1,nr+1,nt2s+1,a->mhist[nd+nt2s*a->logdmax][nr]);}}}
      for (nd=a->logdmin;nd<a->logdmax;nd++){ 
	d=(int)pow(2,nd);
	fprintf(fp,"figure;clf;hold on;\n");
	for (nt2s=0;nt2s<5;nt2s++){
	  fprintf(fp,"subplot(1,5,%d);loglog(%s{%d}(:,%d);title('nt2s%d');\n",nt2s+1,fgvn,nd+1,nt2s+1,nt2s);}}
      fclose(fp);fp=NULL;}}
  else if (dump_type==1){
    for (nd=a->logdmin;nd<a->logdmax;nd++){ d=(int)pow(2,nd); for (nt2s=0;nt2s<5;nt2s++){
      sprintf(filename,"./%s_level%d_t2s%d",fgvn,nd,nt2s);
      ra2jpg(a->mhist[nd+nt2s*a->logdmax],"double",1,(a->Nra->rows/d+1)*(a->Nra->cols/d+1),0,filename,3);}}}
}

/* Here are the power functions */

struct power * powermake(struct neuronarray *Nra,int length,int howmany)
{
  int nr=0;
  struct power *p=NULL;
  p = (struct power *) tmalloc(sizeof(struct power));
  p->update_last = GLOBAL_TI;
  p->update_every=1;
  howmany=minimum(Nra->rows*Nra->cols,howmany);
  p->N=howmany;
  p->length=length;
  p->n = (struct neuron **) tcalloc(p->N,sizeof(struct neuron *));
  p->tst = strobemake(length,1,0);
  p->raster = (struct strobe **) tcalloc(p->N,sizeof(struct strobe *));
  p->Vstra = (struct strobe **) tcalloc(p->N,sizeof(struct strobe *));
  p->sAstra = (struct strobe **) tcalloc(p->N,sizeof(struct strobe *));
  p->sNstra = (struct strobe **) tcalloc(p->N,sizeof(struct strobe *));
  p->sGstra = (struct strobe **) tcalloc(p->N,sizeof(struct strobe *));
  p->VSstra = (struct strobe **) tcalloc(p->N,sizeof(struct strobe *));
  for (nr=0;nr<p->N;nr++){
    p->n[nr] = Nra->N[nr];
    p->raster[nr] = strobemake(p->length,p->update_every,0);
    p->Vstra[nr] = strobemake(p->length,p->update_every,0);
    p->sAstra[nr] = strobemake(p->length,p->update_every,0);
    p->sNstra[nr] = strobemake(p->length,p->update_every,0);
    p->sGstra[nr] = strobemake(p->length,p->update_every,0);
    p->VSstra[nr] = strobemake(p->length,p->update_every,0);}
  p->Vpower = (double *) tcalloc(p->length,sizeof(double));
  p->sApower = (double *) tcalloc(p->length,sizeof(double));
  p->sNpower = (double *) tcalloc(p->length,sizeof(double));
  p->sGpower = (double *) tcalloc(p->length,sizeof(double));
  p->VSpower = (double *) tcalloc(p->length,sizeof(double));
  p->VSbarpower = (double *) tcalloc(p->length,sizeof(double));
  p->update_number=0;
  p->autocorrelation = (double *) tcalloc(p->length,sizeof(double));
  p->crosscorrelation = (double *) tcalloc(p->length,sizeof(double));
  return p;
}

void powertfree(struct power *p)
{
  int nr=0;
  for (nr=0;nr<p->N;nr++){
    strobetfree(p->raster[nr]);p->raster[nr]=NULL;
    strobetfree(p->Vstra[nr]);p->Vstra[nr]=NULL;
    strobetfree(p->sAstra[nr]);p->sAstra[nr]=NULL;
    strobetfree(p->sNstra[nr]);p->sNstra[nr]=NULL;
    strobetfree(p->sGstra[nr]);p->sGstra[nr]=NULL;
    strobetfree(p->VSstra[nr]);p->VSstra[nr]=NULL;}
  tfree(p->n);p->n=NULL;
  strobetfree(p->tst);p->tst=NULL;
  tfree(p->raster);p->raster=NULL;
  tfree(p->Vstra);p->Vstra=NULL;
  tfree(p->sAstra);p->sAstra=NULL;
  tfree(p->sNstra);p->sNstra=NULL;
  tfree(p->sGstra);p->sGstra=NULL;
  tfree(p->VSstra);p->VSstra=NULL;
  tfree(p->Vpower);
  tfree(p->sApower);
  tfree(p->sNpower);
  tfree(p->sGpower);
  tfree(p->VSpower);
  tfree(p->VSbarpower);
  tfree(p->autocorrelation);
  tfree(p->crosscorrelation);
}

void powerupdate(struct power *p,double t,double DT)
{
  int nr=0,tab_old=0,tab_new=0,nl=0;
  struct neuron *n=NULL;
  double *temp=NULL;
  double *ra=NULL;
  double **rara=NULL;
  tab_old = p->tst->tab;
  strobeupdate_old(p->tst,t,t);
  tab_new = p->tst->tab;
  for (nr=0;nr<p->N;nr++){ 
    n=p->n[nr];
    strobeupdate(p->raster[nr],t,DT,n->spikelast==n->spiketime && n->spiketime>=t && n->spiketime<=t+DT);
    strobeupdate(p->Vstra[nr],t,DT,*(n->V));
    strobeupdate(p->sAstra[nr],t,DT,*(n->sA));
    strobeupdate(p->sNstra[nr],t,DT,*(n->sN));
    strobeupdate(p->sGstra[nr],t,DT,*(n->sG));
    strobeupdate(p->VSstra[nr],t,DT,*(n->VS));}
  if (tab_old==p->length-1 && tab_new==0){ /* time to update */
    ra = (double *) tcalloc(p->VSstra[0]->length,sizeof(double));
    rara = (double **) tcalloc(p->N,sizeof(double *));
    for (nr=0;nr<p->N;nr++){
      n=p->n[nr];
      temp = ra2power(NULL,p->Vstra[nr]->data,p->Vstra[nr]->length,NULL,0,0);
      raplusequals(p->Vpower,p->length,temp);
      tfree(temp);temp=NULL;
      temp = ra2power(NULL,p->sAstra[nr]->data,p->sAstra[nr]->length,NULL,0,0);
      raplusequals(p->sApower,p->length,temp);
      tfree(temp);temp=NULL;
      temp = ra2power(NULL,p->sNstra[nr]->data,p->sNstra[nr]->length,NULL,0,0);
      raplusequals(p->sNpower,p->length,temp);
      tfree(temp);temp=NULL;
      temp = ra2power(NULL,p->sGstra[nr]->data,p->sGstra[nr]->length,NULL,0,0);
      raplusequals(p->sGpower,p->length,temp);
      tfree(temp);temp=NULL;
      temp = ra2power(NULL,p->VSstra[nr]->data,p->VSstra[nr]->length,NULL,0,0);
      raplusequals(p->VSpower,p->length,temp);
      tfree(temp);temp=NULL;
      rara[nr] = (double *) tcalloc(p->length,sizeof(double));
      for (nl=0;nl<p->length;nl++){ rara[nr][nl] = p->Vstra[nr]->data[nl]; ra[nl] += p->VSstra[nr]->data[nl];}
      p->update_number += 1;}
    rara2corr(NULL,NULL,rara,p->N,p->length,p->autocorrelation,p->crosscorrelation);
    temp = ra2power(NULL,ra,p->Vstra[0]->length,NULL,0,0); raplusequals(p->VSbarpower,p->length,temp); tfree(temp); temp=NULL;
    tfree(ra);ra=NULL;
    for (nr=0;nr<p->N;nr++){ if (rara[nr]!=NULL){ tfree(rara[nr]); rara[nr]=NULL;}}
    tfree(rara);rara=NULL;}
}

void powerdump(struct power *p,char *fgvn,int dump_type)
{
  char filename[1024];
  char gs2[512];
  FILE *fp=NULL;
  int nr=0;
  if (fgvn==NULL){ sprintf(gs2,"%srecord",GLOBAL_STRING_2);}
  else{ sprintf(gs2,"%srecord_%s",GLOBAL_STRING_2,fgvn);}
  sprintf(filename,"power_%s_raster",gs2); stradump(p->raster,p->N,dump_type,filename);
  sprintf(filename,"power_%s_Vstra",gs2); stradump(p->Vstra,p->N,dump_type,filename);
  sprintf(filename,"power_%s_sAstra",gs2); stradump(p->sAstra,p->N,dump_type,filename);
  sprintf(filename,"power_%s_sNstra",gs2); stradump(p->sNstra,p->N,dump_type,filename);
  sprintf(filename,"power_%s_sGstra",gs2); stradump(p->sGstra,p->N,dump_type,filename);
  sprintf(filename,"power_%s_VSstra",gs2); stradump(p->VSstra,p->N,dump_type,filename);
  if (p->update_number>0){
    if (dump_type==0){
      sprintf(filename,"./power_%s.m",gs2);
      if ((fp=fopen(filename,"w"))!=NULL){
	for (nr=0;nr<p->length;nr++){
	  fprintf(fp,"power_%s_V(%d) = %0.16f;\n",gs2,nr+1,p->Vpower[nr]/(double)p->update_number);
	  fprintf(fp,"power_%s_sA(%d) = %0.16f;\n",gs2,nr+1,p->sApower[nr]/(double)p->update_number);
	  fprintf(fp,"power_%s_sN(%d) = %0.16f;\n",gs2,nr+1,p->sNpower[nr]/(double)p->update_number);
	  fprintf(fp,"power_%s_sG(%d) = %0.16f;\n",gs2,nr+1,p->sGpower[nr]/(double)p->update_number);
	  fprintf(fp,"power_%s_VS(%d) = %0.16f;\n",gs2,nr+1,p->VSpower[nr]/(double)p->update_number);
	  fprintf(fp,"power_%s_VSbar(%d) = %0.16f;\n",gs2,nr+1,p->VSbarpower[nr]/(double)p->update_number);
	  fprintf(fp,"power_%s_autocorrelation(%d) = %0.16f;\n",gs2,nr+1,p->autocorrelation[nr]*p->N/(double)p->update_number);
	  fprintf(fp,"power_%s_crosscorrelation(%d) = %0.16f;\n",gs2,nr+1,p->crosscorrelation[nr]*p->N/(double)p->update_number);}
	fprintf(fp,"figure;clf;hold on;\n");
	fprintf(fp,"subplot(1,8,1);loglog(power_%s_V);title('V');\n",gs2);
	fprintf(fp,"subplot(1,8,2);loglog(power_%s_sA);title('sA');\n",gs2);
	fprintf(fp,"subplot(1,8,3);loglog(power_%s_sN);title('sN');\n",gs2);
	fprintf(fp,"subplot(1,8,4);loglog(power_%s_sG);title('sG');\n",gs2);
	fprintf(fp,"subplot(1,8,5);loglog(power_%s_VS);title('VS');\n",gs2);
	fprintf(fp,"subplot(1,8,6);loglog(power_%s_VSbar);title('VS');\n",gs2);
	fprintf(fp,"subplot(1,8,7);loglog(power_%s_autocorrelation);title('ac');\n",gs2);
	fprintf(fp,"subplot(1,8,8);loglog(power_%s_crosscorrelation);title('xc');\n",gs2);
	fclose(fp);fp=NULL;}}
    else if (dump_type==1){
      sprintf(filename,"./power_%s_V",gs2);
      ra2jpg(p->Vpower,"double",1,p->length,0,filename,3);
      sprintf(filename,"./power_%s_sA",gs2);
      ra2jpg(p->sApower,"double",1,p->length,0,filename,3);
      sprintf(filename,"./power_%s_sN",gs2);
      ra2jpg(p->sNpower,"double",1,p->length,0,filename,3);
      sprintf(filename,"./power_%s_sG",gs2);
      ra2jpg(p->sGpower,"double",1,p->length,0,filename,3);
      sprintf(filename,"./power_%s_VS",gs2);
      ra2jpg(p->VSpower,"double",1,p->length,0,filename,3);
      sprintf(filename,"./power_%s_VSbar",gs2);
      ra2jpg(p->VSbarpower,"double",1,p->length,0,filename,3);
      sprintf(filename,"./power_%s_autocorrelation",gs2);
      ra2jpg(p->autocorrelation,"double",1,p->length,0,filename,0);
      sprintf(filename,"./power_%s_crosscorrelation",gs2);
      ra2jpg(p->crosscorrelation,"double",1,p->length,0,filename,0);}}
}

/* Here are the taof function */

struct taof * taofmake(struct neuronarray *Nra,double length,double step_every)
{
  struct taof *tf=NULL;
  int nr=0,nc=0;
  tf = (struct taof *) tcalloc(1,sizeof(struct taof));
  tf->Nra = Nra;
  tf->length = length;
  tf->step_every = step_every;
  tf->input_contrast = strobemake(tf->length,tf->step_every,0);
  tf->input_spaceangle = strobemake(tf->length,tf->step_every,0);
  tf->f1r = (struct strobe **) tcalloc(tf->Nra->rows*tf->Nra->cols,sizeof(struct strobe *));
  tf->f1i = (struct strobe **) tcalloc(tf->Nra->rows*tf->Nra->cols,sizeof(struct strobe *));
  tf->f0 = (struct strobe **) tcalloc(tf->Nra->rows*tf->Nra->cols,sizeof(struct strobe *));
  for (nr=0;nr<tf->Nra->rows;nr++){ for (nc=0;nc<tf->Nra->cols;nc++){
    tf->f1r[nr+nc*tf->Nra->rows] = strobemake(tf->length,tf->step_every,0);
    tf->f1i[nr+nc*tf->Nra->rows] = strobemake(tf->length,tf->step_every,0);
    tf->f0[nr+nc*tf->Nra->rows] = strobemake(tf->length,tf->step_every,0);}}
  return tf;
}

void taoftfree(struct taof *tf)
{
  int nr=0,nc=0;
  strobetfree(tf->input_contrast);tf->input_contrast=NULL;
  strobetfree(tf->input_spaceangle);tf->input_spaceangle=NULL;
  for (nr=0;nr<tf->Nra->rows;nr++){ for (nc=0;nc<tf->Nra->cols;nc++){
    strobetfree(tf->f1r[nr+nc*tf->Nra->rows]); tf->f1r[nr+nc*tf->Nra->rows] = NULL;
    strobetfree(tf->f1i[nr+nc*tf->Nra->rows]); tf->f1i[nr+nc*tf->Nra->rows] = NULL;
    strobetfree(tf->f0[nr+nc*tf->Nra->rows]); tf->f0[nr+nc*tf->Nra->rows] = NULL;}}
  tfree(tf->f1r);tf->f1r=NULL; 
  tfree(tf->f1i);tf->f1i=NULL; 
  tfree(tf->f0);tf->f0=NULL;
  tfree(tf);tf=NULL;
}

void taofupdate(struct taof *tf,double t,double DT)
{
  int nr=0,nc=0;
  struct neuron *n=NULL;
  double spikeflag=0,temp=0,tempc=0,temps=0;
  strobeupdate(tf->input_contrast,t,DT,INPUT_CONTRAST);
  strobeupdate(tf->input_spaceangle,t,DT,INPUT_SPACEANGLE);
  temp = 2*PI*(t-tf->input_contrast->last_time)/(GLOBAL_SPF*GRATING_DRIFT);
  tempc = cos(temp); temps = sin(temp);
  for (nr=0;nr<tf->Nra->rows;nr++){ for (nc=0;nc<tf->Nra->cols;nc++){
    n=nget(tf->Nra,nr,nc);
    spikeflag = (n->spikelast==n->spiketime && n->spiketime>=t && n->spiketime<=t+DT);
    strobeupdate_sum(tf->f1r[nr+nc*tf->Nra->rows],t,DT,spikeflag*tempc/tf->step_every);
    strobeupdate_sum(tf->f1i[nr+nc*tf->Nra->rows],t,DT,spikeflag*temps/tf->step_every);
    strobeupdate_sum(tf->f0[nr+nc*tf->Nra->rows],t,DT,spikeflag/tf->step_every);}}
}

double * taofdump_helper(double *ra,int length,int nheight,int leftmost_vs_rightmost,double *max_height_out,double *min_height_out)
{
  /* given ra[index]=height, returns ra_transpose[height]=index 
     the flag leftmost_vs_rightmost determines which multiply valued index is chosen */
  int verbose=0;
  int nl=0,nr=0;
  double tolerance=0.000001;
  double *ra_transpose=NULL;
  double max_height=0,min_height=0,temp_height=0;
  if (verbose){ printf(" %% [entering taofdump_helper]\n");}
  nheight = maximum(2,nheight);
  ra_transpose = (double *) tcalloc(nheight,sizeof(double));
  stats("double",ra,length,&max_height,&min_height,NULL,NULL); 
  if (verbose){ printf(" %% height within range (%f,%f)\n",min_height,max_height);}
  for (nl=0;nl<nheight;nl++){
    temp_height = min_height + nl*(max_height-min_height)/(nheight-1);
    if (verbose){ printf(" %% \n");}
    if (verbose){ printf(" %% nl %d temp_height %f, ",nl,temp_height);}
    if (leftmost_vs_rightmost){
      nr=0; ra_transpose[nl]=-1;
      while (nr+1<length && ra_transpose[nl]<0){
	if ((temp_height-minimum(ra[nr],ra[nr+1]))>=-tolerance && (maximum(ra[nr],ra[nr+1])-temp_height)>=-tolerance){ 
	  ra_transpose[nl] = nr+0.5;
	  if ((ra[nr]-temp_height>=-tolerance && temp_height-ra[nr+1]>=-tolerance) || (ra[nr+1]-temp_height>=-tolerance && temp_height-ra[nr]>=-tolerance)){ ra_transpose[nl]= nr + (temp_height-ra[nr])/(ra[nr+1]-ra[nr]);}}
	nr+=1;}
      if (ra_transpose[nl]<0 && nr+1>=length){ ra_transpose[nl]=length-1;}}
    else /* if (!leftmost_vs_rightmost) */{
      nr=length-1; ra_transpose[nl]=-1;
      while (nr>0 && ra_transpose[nl]<0){
	if ((temp_height-minimum(ra[nr],ra[nr-1]))>=-tolerance && (maximum(ra[nr],ra[nr-1])-temp_height)>=-tolerance){ 
	  ra_transpose[nl] = nr+0.5;
	  if ((ra[nr]-temp_height>=-tolerance && temp_height-ra[nr-1]>=-tolerance) || (ra[nr-1]-temp_height>=-tolerance && temp_height-ra[nr]>=-tolerance)){ ra_transpose[nl]= nr-1 + (temp_height-ra[nr-1])/(ra[nr]-ra[nr-1]);}}
	nr-=1;}
      if (ra_transpose[nl]<0 && nr-1<0){ ra_transpose[nl]=0;}}
    if (verbose){ printf(" %% ra_transpose[%d] set to %2f\n",nl,ra_transpose[nl]);}}
  if (max_height_out!=NULL){ *max_height_out = max_height;}
  if (min_height_out!=NULL){ *min_height_out = min_height;}
  return ra_transpose;
}

void taofdump(struct taof *tf,char *fgvn,int dump_type)
{
  /* dump_flag:
     3: histogram of totalspikes_contrast_increase vs totalspikes_contrast_decrease 
     4: histogram of half width (at max height) (in terms of contrast) of hysteresis loop
  */
  char filename[1024],command[2048],tempchar[512];
  char gs2[512];
  FILE **fpra=NULL,*fptmp=NULL;
  int nr=0,nc=0,nl=0,nl2=0,nl3=0,tab=0,nt2s=0;
  double maxdia = 10000;
  double maxc=0,minc=0;
  double maxme=0,minme=0,maxmi=0,minmi=0,maxtmp=0,mintmp=0;
  double f0leg1=0,f0leg2=0;
  double rcolor=0,gcolor=0,bcolor=0;
  struct neuron *n=NULL;
  int colorcode=0,index=0,length=0;
  double xpos=0,ypos=0;
  int drawtype_flag=1;
  int jpeg_flag=(ON_MY_COMPUTER ? 1 : 0),remove_flag=0;
  double *temp=NULL;  
  struct llist **Lra=NULL;
  struct hist **hra=NULL;
  double *ra1=NULL,*ra2=NULL,*ra1r=NULL,*ra1l=NULL,*ra2r=NULL,*ra2l=NULL;
  int halflength=0,within_flag_1=0,within_flag_2=0;
  double max1=0,max2=0,min1=0,min2=0,halfheight=0,index_positive=0,index_negative=0,index1l=0,index1r=0,index2l=0,index2r=0;
  if (fgvn==NULL){ sprintf(gs2,"%srecord",GLOBAL_STRING_2);}
  else{ sprintf(gs2,"%srecord_%s",GLOBAL_STRING_2,fgvn);}
  if (dump_type==0 || dump_type==1){
    sprintf(filename,"taof_%s_input_contrast",gs2); stradump(&(tf->input_contrast),1,dump_type,filename);
    sprintf(filename,"taof_%s_input_spaceangle",gs2); stradump(&(tf->input_spaceangle),1,dump_type,filename);
    sprintf(filename,"taof_%s_f1r",gs2); stradump(tf->f1r,tf->Nra->rows*tf->Nra->cols,dump_type,filename);
    sprintf(filename,"taof_%s_f1i",gs2); stradump(tf->f1i,tf->Nra->rows*tf->Nra->cols,dump_type,filename);
    sprintf(filename,"taof_%s_f0",gs2); stradump(tf->f0,tf->Nra->rows*tf->Nra->cols,dump_type,filename);
    fpra = (FILE **) tcalloc(4,sizeof(FILE *));
    for (nt2s=0;nt2s<4;nt2s++){
      sprintf(filename,"./taof_%s_f0_t2s%d.fig",gs2,nt2s);
      if ((fpra[nt2s]=fopen(filename,"w"))==NULL){ printf(" warning, couldn't open %s in taofdump\n",filename); fpra[nt2s]=stdout;}
      fprintf(fpra[nt2s],"%s",FIG_PREAMBLE); fprintf(fpra[nt2s],"%s",FIG_PREAMBLE_COLOR_7);}
    maxme=(double)30/(double)1024;minme=0;maxmi=(double)300/(double)1024;minmi=0;
    stats("double",tf->input_contrast->data,tf->input_contrast->length,&maxc,&minc,NULL,NULL);
    for (nr=0;nr<tf->Nra->rows;nr++){ for (nc=0;nc<tf->Nra->cols;nc++){
      n=nget(tf->Nra,nr,nc); index = n->row+n->col*tf->Nra->rows; length = minimum(tf->input_contrast->length,tf->f0[index]->length);
      fptmp = fpra[(int)*(n->t2s)];maxtmp=(n->type==+1?maxme:maxmi);mintmp=(n->type==+1?minme:minmi);
      f0leg1=0;f0leg2=0;
      for (nl=0;nl<length;nl++){
	tab = periodize(tf->f0[index]->tab+nl,0,tf->f0[index]->length);
	if (nl<=(tf->length-1)/2){ f0leg1+=tf->f0[index]->data[tab];}
	if (nl>=(tf->length-1)/2){ f0leg2+=tf->f0[index]->data[tab];}}
      colorscale(0,(f0leg1-f0leg2)/maximum(f0leg1,f0leg2),+1,-1,&rcolor,&gcolor,&bcolor);
      colorcode = crop((int)floor(512*rcolor),0,511);
      switch (drawtype_flag){
      case 0: /* firing rate versus contrast */
	fprintf(fptmp,"2 1 0 5 %d 7 %d 0 -1 0.000 0 0 -1 0 0 %d\n\t",/*color*/colorcode+32,/*depth*/index%999+1,/*npoints*/length);
	for (nl=0;nl<length;nl++){
	  tab = periodize(tf->input_contrast->tab+nl,0,tf->input_contrast->length);
	  xpos = (tf->input_contrast->data[tab]-minc)/(maxc-minc);
	  tab = periodize(tf->f0[index]->tab+nl,0,tf->f0[index]->length);
	  ypos = (tf->f0[index]->data[tab]-mintmp)/(maxtmp-mintmp);
	  fprintf(fptmp," %d %d",(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos));}
	fprintf(fptmp,"\n");
	break;
      case 1: /* firing rate on upstroke versus firing rate on downstroke */
	nl=0; tab = periodize(tf->input_contrast->tab+nl,0,tf->input_contrast->length);
	while (fabs(tf->input_contrast->data[tab]-maxc)>0.0000001 && nl<length){ 
	  nl+=1; tab = periodize(tf->input_contrast->tab+nl,0,tf->input_contrast->length);}
	if (nl>=length){ printf(" %% warning! odd tf->input_contrast in taofdump\n"); nl=0;}
	nl2=nl; nl3=minimum(nl2,length-1-nl2);
	fprintf(fptmp,"2 1 0 5 %d 7 %d 0 -1 0.000 0 0 -1 0 0 %d\n\t",/*color*/colorcode+32,/*depth*/index%999+1,/*npoints*/nl3+1);
	for (nl=0;nl<=nl3;nl++){
	  tab = periodize(tf->f0[index]->tab+nl2-nl,0,tf->f0[index]->length);
	  xpos = (tf->f0[index]->data[tab]-mintmp)/(maxtmp-mintmp);
	  tab = periodize(tf->f0[index]->tab+nl2+nl,0,tf->f0[index]->length);
	  ypos = (tf->f0[index]->data[tab]-mintmp)/(maxtmp-mintmp);
	  fprintf(fptmp," %d %d",(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos));}
	fprintf(fptmp,"\n");
	break;
      default: break;}}}
    fprintf(fpra[0],"4 0 0 %d 0 %d %d 0.0000 4 %d %d %d %d %0.16f\\001\n",/*depth*/1,/*font*/12,/*point*/24,/*textheight*/(int)(24*13.5),/*textwidth*/(int)(24*9*18),/*xpos*/(int)(0.1*maxdia),/*ypos*/(int)maxdia-(int)(1.0*maxdia),maxmi);
    fprintf(fpra[0],"4 0 0 %d 0 %d %d 0.0000 4 %d %d %d %d %0.16f\\001\n",/*depth*/1,/*font*/12,/*point*/24,/*textheight*/(int)(24*13.5),/*textwidth*/(int)(24*9*18),/*xpos*/(int)(0.1*maxdia),/*ypos*/(int)maxdia-(int)(0.0*maxdia),minmi);
    fprintf(fpra[1],"4 0 0 %d 0 %d %d 0.0000 4 %d %d %d %d %0.16f\\001\n",/*depth*/1,/*font*/12,/*point*/24,/*textheight*/(int)(24*13.5),/*textwidth*/(int)(24*9*18),/*xpos*/(int)(0.1*maxdia),/*ypos*/(int)maxdia-(int)(1.0*maxdia),maxmi);
    fprintf(fpra[1],"4 0 0 %d 0 %d %d 0.0000 4 %d %d %d %d %0.16f\\001\n",/*depth*/1,/*font*/12,/*point*/24,/*textheight*/(int)(24*13.5),/*textwidth*/(int)(24*9*18),/*xpos*/(int)(0.1*maxdia),/*ypos*/(int)maxdia-(int)(0.0*maxdia),minmi);
    fprintf(fpra[2],"4 0 0 %d 0 %d %d 0.0000 4 %d %d %d %d %0.16f\\001\n",/*depth*/1,/*font*/12,/*point*/24,/*textheight*/(int)(24*13.5),/*textwidth*/(int)(24*9*18),/*xpos*/(int)(0.1*maxdia),/*ypos*/(int)maxdia-(int)(1.0*maxdia),maxme);
    fprintf(fpra[2],"4 0 0 %d 0 %d %d 0.0000 4 %d %d %d %d %0.16f\\001\n",/*depth*/1,/*font*/12,/*point*/24,/*textheight*/(int)(24*13.5),/*textwidth*/(int)(24*9*18),/*xpos*/(int)(0.1*maxdia),/*ypos*/(int)maxdia-(int)(0.0*maxdia),minme);
    fprintf(fpra[3],"4 0 0 %d 0 %d %d 0.0000 4 %d %d %d %d %0.16f\\001\n",/*depth*/1,/*font*/12,/*point*/24,/*textheight*/(int)(24*13.5),/*textwidth*/(int)(24*9*18),/*xpos*/(int)(0.1*maxdia),/*ypos*/(int)maxdia-(int)(1.0*maxdia),maxme);
    fprintf(fpra[3],"4 0 0 %d 0 %d %d 0.0000 4 %d %d %d %d %0.16f\\001\n",/*depth*/1,/*font*/12,/*point*/24,/*textheight*/(int)(24*13.5),/*textwidth*/(int)(24*9*18),/*xpos*/(int)(0.1*maxdia),/*ypos*/(int)maxdia-(int)(0.0*maxdia),minme);
    for (nt2s=0;nt2s<4;nt2s++){ 
      if (fpra[nt2s]!=stdout){ fclose(fpra[nt2s]); fpra[nt2s]=NULL;}
      sprintf(filename,"./taof_%s_f0_t2s%d.fig",gs2,nt2s);
      if (tf->Nra->t2stotal[nt2s]>0){ if (jpeg_flag){ sprintf(command,"fig2dev -Ljpeg -q %d %s ./taof_%s_f0_t2s%d.jpg;",/*quality*/5,filename,gs2,nt2s); system(command);}}
      if (remove_flag){ sprintf(command,"rm %s;",filename); system(command);}}
    tfree(fpra); fpra=NULL;}
  else if (dump_type==3){
    stats("double",tf->f0[0]->data,tf->f0[0]->length,&max2,&min2,NULL,NULL);
    for (index=0;index<tf->Nra->rows*tf->Nra->cols;index++){ stats("double",tf->f0[index]->data,tf->f0[index]->length,&max1,&min1,NULL,NULL); max2 = maximum(max2,max1); min2 = minimum(min2,min1);}
    fpra = (FILE **) tcalloc(4,sizeof(FILE *));
    for (nt2s=0;nt2s<4;nt2s++){
      sprintf(filename,"./taof_%s_curves_t2s%d.fig",gs2,nt2s);
      if ((fpra[nt2s]=fopen(filename,"w"))==NULL){ printf(" warning, couldn't open %s in taofdump\n",filename); fpra[nt2s]=stdout;}
      fprintf(fpra[nt2s],"%s",FIG_PREAMBLE); fprintf(fpra[nt2s],"%s",FIG_PREAMBLE_COLOR_7);}
    Lra = (struct llist **) tcalloc(4,sizeof(struct llist *)); for (nt2s=0;nt2s<4;nt2s++){ Lra[nt2s] = llistmake();}
    for (nr=0;nr<tf->Nra->rows;nr++){ for (nc=0;nc<tf->Nra->cols;nc++){
      n=nget(tf->Nra,nr,nc); index = n->row+n->col*tf->Nra->rows; length = minimum(tf->input_contrast->length,tf->f0[index]->length);
      f0leg1=0;f0leg2=0;
      for (nl=0;nl<length;nl++){
	tab = periodize(tf->f0[index]->tab+nl,0,tf->f0[index]->length);
	if (nl<=(tf->length-1)/2){ f0leg1+=tf->f0[index]->data[tab];}
	if (nl>=(tf->length-1)/2){ f0leg2+=tf->f0[index]->data[tab];}}
      temp=(double *)tcalloc(1,sizeof(double));*temp=f0leg1-f0leg2;litemadd(Lra[(int)*(n->t2s)],temp);}}
    for (nt2s=0;nt2s<4;nt2s++){ if (Lra[nt2s]->length>0){ llistsort(Lra[nt2s]->first,Lra[nt2s]->last,Lra[nt2s]->length,&double_compare);}}
    for (nr=0;nr<tf->Nra->rows;nr++){ for (nc=0;nc<tf->Nra->cols;nc++){
      n=nget(tf->Nra,nr,nc); index = n->row+n->col*tf->Nra->rows; length = minimum(tf->input_contrast->length,tf->f0[index]->length);
      f0leg1=0;f0leg2=0;
      for (nl=0;nl<length;nl++){
	tab = periodize(tf->f0[index]->tab+nl,0,tf->f0[index]->length);
	if (nl<=(tf->length-1)/2){ f0leg1+=tf->f0[index]->data[tab];}
	if (nl>=(tf->length-1)/2){ f0leg2+=tf->f0[index]->data[tab];}}
      min1 = *(double *)(Lra[(int)*(n->t2s)]->first->item);
      max1 = *(double *)(Lra[(int)*(n->t2s)]->last->item);
      if (fabs(f0leg1-f0leg2-max1)<(max1-min1)/128.0 || fabs(f0leg1-f0leg2-min1)<(max1-min1)/128.0){
	halflength = (length-1)/2;
	fprintf(fpra[(int)*(n->t2s)],"2 1 0 5 %d 7 %d 0 -1 0.000 0 0 -1 0 0 %d\n\t",/*color*/0,/*depth*/0,/*npoints*/halflength);
	for (nl=0;nl<=halflength;nl++){
	  tab = periodize(tf->f0[index]->tab+nl,0,tf->f0[index]->length);
	  xpos = ((double)nl+0.5)/(double)halflength;
	  ypos = (tf->f0[index]->data[tab] - min2)/(max2-min2);
	  fprintf(fpra[(int)*(n->t2s)]," %d %d",(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos));}
	fprintf(fpra[(int)*(n->t2s)],"\n");
	fprintf(fpra[(int)*(n->t2s)],"2 1 0 5 %d 7 %d 0 -1 0.000 0 0 -1 0 0 %d\n\t",/*color*/1,/*depth*/0,/*npoints*/halflength);
	for (nl=halflength;nl<length;nl++){
	  tab = periodize(tf->f0[index]->tab+nl,0,tf->f0[index]->length);
	  xpos = ((double)length-1-nl+0.5)/(double)halflength;
	  ypos = (tf->f0[index]->data[tab]-min2)/(max2-min2);
	  fprintf(fpra[(int)*(n->t2s)]," %d %d",(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos));}
	fprintf(fpra[(int)*(n->t2s)],"\n");}}}
    for (nt2s=0;nt2s<4;nt2s++){
      fprintf(fpra[nt2s],"4 0 0 %d 0 %d %d 0.0000 4 %d %d %d %d %0.16f\\001\n",/*depth*/1,/*font*/12,/*point*/24,/*textheight*/(int)(24*13.5),/*textwidth*/(int)(24*9*18),/*xpos*/(int)(0.1*maxdia),/*ypos*/(int)maxdia-(int)(1.0*maxdia),max2);
      fprintf(fpra[nt2s],"4 0 0 %d 0 %d %d 0.0000 4 %d %d %d %d %0.16f\\001\n",/*depth*/1,/*font*/12,/*point*/24,/*textheight*/(int)(24*13.5),/*textwidth*/(int)(24*9*18),/*xpos*/(int)(0.1*maxdia),/*ypos*/(int)maxdia-(int)(0.0*maxdia),min2);}
    hra = (struct hist **) tcalloc(4,sizeof(struct hist *)); 
    for (nt2s=0;nt2s<4;nt2s++){ 
      if (Lra[nt2s]->length>0){
	lliststats(Lra[nt2s],&maxtmp,&mintmp,NULL,NULL);
	hra[nt2s] = histmake(64,maxtmp,mintmp);
	llist2hist(Lra[nt2s],hra[nt2s]);
	sprintf(filename,"./taof_hist_totalspikes_t2s%d_%s",nt2s,gs2); sprintf(tempchar,"taof_hist_totalspikes_t2s%d_%s",nt2s,gs2);
	histdump(hra[nt2s],0,filename,tempchar,0);
	histtfree(hra[nt2s]);hra[nt2s]=NULL;}}
    for (nt2s=0;nt2s<4;nt2s++){ 
      if (fpra[nt2s]!=stdout){ fclose(fpra[nt2s]); fpra[nt2s]=NULL;}
      if (Lra[nt2s]->length>0){
	sprintf(filename,"./taof_%s_curves_t2s%d.fig",gs2,nt2s);
	if (jpeg_flag){ sprintf(command,"fig2dev -Ljpeg -q %d %s ./taof_%s_curves_t2s%d.jpg;",/*quality*/5,filename,gs2,nt2s); system(command);}
	if (remove_flag){ sprintf(command,"rm %s;",filename); system(command);}}
      llisttfree2(Lra[nt2s]);Lra[nt2s]=NULL;}
    tfree(fpra); fpra=NULL;
    tfree(Lra);Lra=NULL;
    tfree(hra);hra=NULL;}
  else if (dump_type==4){
    Lra = (struct llist **) tcalloc(4,sizeof(struct llist *)); for (nt2s=0;nt2s<4;nt2s++){ Lra[nt2s] = llistmake();}
    for (nr=0;nr<tf->Nra->rows;nr++){ for (nc=0;nc<tf->Nra->cols;nc++){
      n=nget(tf->Nra,nr,nc); index = n->row+n->col*tf->Nra->rows; length = minimum(tf->input_contrast->length,tf->f0[index]->length);
      halflength = (length-1)/2;
      ra1 = (double *) tcalloc(halflength+1,sizeof(double));
      ra2 = (double *) tcalloc(halflength+1,sizeof(double));
      for (nl=0;nl<length;nl++){
	tab = periodize(tf->f0[index]->tab+nl,0,tf->f0[index]->length);
	if (nl<=halflength){ ra1[nl] = tf->f0[index]->data[tab];}
	if (nl>=halflength){ ra2[length-1-nl] = tf->f0[index]->data[tab];}}
      ra1l = taofdump_helper(ra1,halflength+1,2*halflength,1,&max1,&min1);
      ra1r = taofdump_helper(ra1,halflength+1,2*halflength,0,&max1,&min1);
      ra2l = taofdump_helper(ra2,halflength+1,2*halflength,1,&max2,&min2);
      ra2r = taofdump_helper(ra2,halflength+1,2*halflength,0,&max2,&min2);
      halfheight = (minimum(min1,min2) + maximum(max1,max2))/2.0;
      tab = (int)floor(2*halflength*(halfheight - min1)/(max1 - min1));
      if (tab>=0 && tab<2*halflength){ within_flag_1 = 0; index1l = ra1l[tab]; index1r = ra1r[tab];}
      else if (tab<0){ within_flag_1 = -1;} else if (tab>=2*halflength){ within_flag_1 = +1;}
      tab = (int)floor(2*halflength*(halfheight - min2)/(max2 - min2));
      if (tab>=0 && tab<2*halflength){ within_flag_2 = 0; index2l = ra2l[tab]; index2r = ra2r[tab];}
      else if (tab<0){ within_flag_2 = -1;} else if (tab>=2*halflength){ within_flag_2 = +1;}
      index_positive=0,index_negative=0;
      if (within_flag_1==0 && within_flag_2==0){
	index_positive = index1r - index2l;
	index_negative = index2r - index1l;}
      else if (within_flag_1==0){
	switch (within_flag_2){
	case  1: index_positive = index1r - index1l; index_negative = 0; break;
	case -1: index_negative = index1r - index1l; index_positive = 0; break;
	default: break;}}
      else if (within_flag_2==0){
	switch (within_flag_1){
	case  1: index_negative = index2r - index2l; index_positive = 0; break;
	case -1: index_positive = index2r - index2l; index_negative = 0; break;
	default: break;}}
      temp = (double *)tcalloc(1,sizeof(double));
      if (fabs(index_positive) > fabs(index_negative)){ *temp=index_positive;} 
      else if (fabs(index_positive) < fabs(index_negative)){ *temp=-index_negative;} 
      else if (fabs(index_positive) == fabs(index_negative)){ *temp=(index_positive + index_negative)/2.0;} 
      litemadd(Lra[(int)*(n->t2s)],temp);
      tfree(ra1);ra1=NULL; tfree(ra1r);ra1r=NULL; tfree(ra1l);ra1l=NULL;
      tfree(ra2);ra2=NULL; tfree(ra2r);ra2r=NULL; tfree(ra2l);ra2l=NULL;}}
    hra = (struct hist **) tcalloc(4,sizeof(struct hist *)); 
    for (nt2s=0;nt2s<4;nt2s++){ 
      if (Lra[nt2s]->length>0){
	lliststats(Lra[nt2s],&maxtmp,&mintmp,NULL,NULL);
	hra[nt2s] = histmake(64,maxtmp,mintmp);
	llist2hist(Lra[nt2s],hra[nt2s]);
	sprintf(filename,"./taof_hist_halfwidth_t2s%d_%s",nt2s,gs2); sprintf(tempchar,"taof_hist_halfwidth_t2s%d_%s",nt2s,gs2);
	histdump(hra[nt2s],0,filename,tempchar,0);
	histtfree(hra[nt2s]);hra[nt2s]=NULL;}
      llisttfree2(Lra[nt2s]);Lra[nt2s]=NULL;}
    tfree(Lra);Lra=NULL;
    tfree(hra);hra=NULL;}
}

/* Here are the seidcorr functions */

struct seidcorr * seidcorrmake(struct neuronarray *Nra,int space_bin_size,int time_bin_size,int length,double time_start)
{
  int nr=0,nc=0;
  struct seidcorr *s=NULL;
  s = (struct seidcorr *) tcalloc(1,sizeof(struct seidcorr));
  s->Nra = Nra;
  s->space_bin_size = space_bin_size;
  s->time_bin_size = time_bin_size;
  s->length = length;
  s->time_start = time_start;
  s->time_end = s->time_start + s->length*s->time_bin_size;
  s->rows = (int)floor(s->Nra->rows/s->space_bin_size);
  s->cols = (int)floor(s->Nra->cols/s->space_bin_size);
  s->Vstra = (struct strobe **) tcalloc(s->rows*s->cols,sizeof(struct strobe *));
  s->sAstra = (struct strobe **) tcalloc(s->rows*s->cols,sizeof(struct strobe *));
  s->sNstra = (struct strobe **) tcalloc(s->rows*s->cols,sizeof(struct strobe *));
  s->sGstra = (struct strobe **) tcalloc(s->rows*s->cols,sizeof(struct strobe *));
  s->VSstra = (struct strobe **) tcalloc(s->rows*s->cols,sizeof(struct strobe *));
  for (nr=0;nr<s->rows;nr++){ for (nc=0;nc<s->cols;nc++){ 
    s->Vstra[nr+nc*s->rows] = strobemake(s->length,s->time_bin_size,0); s->Vstra[nr+nc*s->rows]->last_time = s->time_start;
    s->sAstra[nr+nc*s->rows] = strobemake(s->length,s->time_bin_size,0); s->sAstra[nr+nc*s->rows]->last_time = s->time_start;
    s->sNstra[nr+nc*s->rows] = strobemake(s->length,s->time_bin_size,0); s->sNstra[nr+nc*s->rows]->last_time = s->time_start;
    s->sGstra[nr+nc*s->rows] = strobemake(s->length,s->time_bin_size,0); s->sGstra[nr+nc*s->rows]->last_time = s->time_start;
    s->VSstra[nr+nc*s->rows] = strobemake(s->length,s->time_bin_size,0); s->VSstra[nr+nc*s->rows]->last_time = s->time_start;}}
  return s;
}

void seidcorrtfree(struct seidcorr *s)
{
  int nr=0,nc=0;
  for (nr=0;nr<s->rows;nr++){ for (nc=0;nc<s->cols;nc++){ 
    strobetfree(s->Vstra[nr+nc*s->rows]); s->Vstra[nr+nc*s->rows]=NULL;
    strobetfree(s->sAstra[nr+nc*s->rows]); s->sAstra[nr+nc*s->rows]=NULL;
    strobetfree(s->sNstra[nr+nc*s->rows]); s->sNstra[nr+nc*s->rows]=NULL;
    strobetfree(s->sGstra[nr+nc*s->rows]); s->sGstra[nr+nc*s->rows]=NULL;
    strobetfree(s->VSstra[nr+nc*s->rows]); s->VSstra[nr+nc*s->rows]=NULL;}}
  tfree(s->Vstra);s->Vstra=NULL;
  tfree(s->sAstra);s->sAstra=NULL;
  tfree(s->sNstra);s->sNstra=NULL;
  tfree(s->sGstra);s->sGstra=NULL;
  tfree(s->VSstra);s->VSstra=NULL;
  tfree(s);s=NULL;
}

void seidcorrupdate(struct seidcorr *s,double t,double DT)
{
  int nr=0,nc=0,nr2=0,nc2=0;
  double *ra=NULL;
  if (t>=s->time_start && t<=s->time_end){
    ra = spacesmear(s->Nra->Vra,s->Nra->rows,s->Nra->cols,s->space_bin_size/2);
    for (nr=0;nr<s->rows;nr++){ for (nc=0;nc<s->cols;nc++){
      nr2 = nr*s->space_bin_size; nc2 = nc*s->space_bin_size;
      strobeupdate(s->Vstra[nr+nc*s->rows],t,DT,ra[nr2+nc2*s->Nra->rows]);}}
    tfree(ra);ra=NULL;
    ra = spacesmear(s->Nra->sAra,s->Nra->rows,s->Nra->cols,s->space_bin_size/2);
    for (nr=0;nr<s->rows;nr++){ for (nc=0;nc<s->cols;nc++){
      nr2 = nr*s->space_bin_size; nc2 = nc*s->space_bin_size;
      strobeupdate(s->sAstra[nr+nc*s->rows],t,DT,ra[nr2+nc2*s->Nra->rows]);}}
    tfree(ra);ra=NULL;
    ra = spacesmear(s->Nra->sNra,s->Nra->rows,s->Nra->cols,s->space_bin_size/2);
    for (nr=0;nr<s->rows;nr++){ for (nc=0;nc<s->cols;nc++){
      nr2 = nr*s->space_bin_size; nc2 = nc*s->space_bin_size;
      strobeupdate(s->sNstra[nr+nc*s->rows],t,DT,ra[nr2+nc2*s->Nra->rows]);}}
    tfree(ra);ra=NULL;
    ra = spacesmear(s->Nra->sGra,s->Nra->rows,s->Nra->cols,s->space_bin_size/2);
    for (nr=0;nr<s->rows;nr++){ for (nc=0;nc<s->cols;nc++){
      nr2 = nr*s->space_bin_size; nc2 = nc*s->space_bin_size;
      strobeupdate(s->sGstra[nr+nc*s->rows],t,DT,ra[nr2+nc2*s->Nra->rows]);}}
    tfree(ra);ra=NULL;
    ra = spacesmear(s->Nra->VSra,s->Nra->rows,s->Nra->cols,s->space_bin_size/2);
    for (nr=0;nr<s->rows;nr++){ for (nc=0;nc<s->cols;nc++){
      nr2 = nr*s->space_bin_size; nc2 = nc*s->space_bin_size;
      strobeupdate(s->VSstra[nr+nc*s->rows],t,DT,ra[nr2+nc2*s->Nra->rows]);}}
    tfree(ra);ra=NULL;}
  if (t<s->time_end && t+DT>=s->time_end){ seidcorrdump(s);}
}

void seidcorrdump(struct seidcorr *s)
{
  char filename[1024];
  int nr=0,nc=0,nl=0,tab=0;
  int *ira=NULL;
  double *ra=NULL;
  sprintf(filename,"./seidcorr_%s_%d_Vra",GLOBAL_STRING,GLOBAL_RECORD_NUMBER);
  if (checktofind(filename)){ printf(" %% warning, filename %s already exists in seidcorrdump\n",filename);}
  else /* if (!checktofind(filename)) */{ 
    ra=(double *) tcalloc(s->rows*s->cols*s->length,sizeof(double));
    for (nr=0;nr<s->rows;nr++){ for (nc=0;nc<s->cols;nc++){ for (nl=0;nl<s->length;nl++){
      tab = periodize(s->Vstra[nr+nc*s->rows]->tab+nl,0,s->Vstra[nr+nc*s->rows]->length);
      ra[nr+nc*s->rows+nl*s->rows*s->cols] = s->Vstra[nr+nc*s->rows]->data[tab];}}}
    ira = (int *) tcalloc(3,sizeof(int)); ira[0]=s->rows; ira[1]=s->cols; ira[2]=s->length;
    multidradump(ra,3,ira,filename);
    tfree(ira);ira=NULL;
    tfree(ra);ra=NULL;}
  sprintf(filename,"./seidcorr_%s_%d_sAra",GLOBAL_STRING,GLOBAL_RECORD_NUMBER);
  if (checktofind(filename)){ printf(" %% warning, filename %s already exists in seidcorrdump\n",filename);}
  else /* if (!checktofind(filename)) */{ 
    ra=(double *) tcalloc(s->rows*s->cols*s->length,sizeof(double));
    for (nr=0;nr<s->rows;nr++){ for (nc=0;nc<s->cols;nc++){ for (nl=0;nl<s->length;nl++){
      tab = periodize(s->sAstra[nr+nc*s->rows]->tab+nl,0,s->sAstra[nr+nc*s->rows]->length);
      ra[nr+nc*s->rows+nl*s->rows*s->cols] = s->sAstra[nr+nc*s->rows]->data[tab];}}}
    ira = (int *) tcalloc(3,sizeof(int)); ira[0]=s->rows; ira[1]=s->cols; ira[2]=s->length;
    multidradump(ra,3,ira,filename);
    tfree(ira);ira=NULL;
    tfree(ra);ra=NULL;}
  sprintf(filename,"./seidcorr_%s_%d_sNra",GLOBAL_STRING,GLOBAL_RECORD_NUMBER);
  if (checktofind(filename)){ printf(" %% warning, filename %s already exists in seidcorrdump\n",filename);}
  else /* if (!checktofind(filename)) */{ 
    ra=(double *) tcalloc(s->rows*s->cols*s->length,sizeof(double));
    for (nr=0;nr<s->rows;nr++){ for (nc=0;nc<s->cols;nc++){ for (nl=0;nl<s->length;nl++){
      tab = periodize(s->sNstra[nr+nc*s->rows]->tab+nl,0,s->sNstra[nr+nc*s->rows]->length);
      ra[nr+nc*s->rows+nl*s->rows*s->cols] = s->sNstra[nr+nc*s->rows]->data[tab];}}}
    ira = (int *) tcalloc(3,sizeof(int)); ira[0]=s->rows; ira[1]=s->cols; ira[2]=s->length;
    multidradump(ra,3,ira,filename);
    tfree(ira);ira=NULL;
    tfree(ra);ra=NULL;}
  sprintf(filename,"./seidcorr_%s_%d_sGra",GLOBAL_STRING,GLOBAL_RECORD_NUMBER);
  if (checktofind(filename)){ printf(" %% warning, filename %s already exists in seidcorrdump\n",filename);}
  else /* if (!checktofind(filename)) */{ 
    ra=(double *) tcalloc(s->rows*s->cols*s->length,sizeof(double));
    for (nr=0;nr<s->rows;nr++){ for (nc=0;nc<s->cols;nc++){ for (nl=0;nl<s->length;nl++){
      tab = periodize(s->sGstra[nr+nc*s->rows]->tab+nl,0,s->sGstra[nr+nc*s->rows]->length);
      ra[nr+nc*s->rows+nl*s->rows*s->cols] = s->sGstra[nr+nc*s->rows]->data[tab];}}}
    ira = (int *) tcalloc(3,sizeof(int)); ira[0]=s->rows; ira[1]=s->cols; ira[2]=s->length;
    multidradump(ra,3,ira,filename);
    tfree(ira);ira=NULL;
    tfree(ra);ra=NULL;}
  sprintf(filename,"./seidcorr_%s_%d_VSra",GLOBAL_STRING,GLOBAL_RECORD_NUMBER);
  if (checktofind(filename)){ printf(" %% warning, filename %s already exists in seidcorrdump\n",filename);}
  else /* if (!checktofind(filename)) */{ 
    ra=(double *) tcalloc(s->rows*s->cols*s->length,sizeof(double));
    for (nr=0;nr<s->rows;nr++){ for (nc=0;nc<s->cols;nc++){ for (nl=0;nl<s->length;nl++){
      tab = periodize(s->VSstra[nr+nc*s->rows]->tab+nl,0,s->VSstra[nr+nc*s->rows]->length);
      ra[nr+nc*s->rows+nl*s->rows*s->cols] = s->VSstra[nr+nc*s->rows]->data[tab];}}}
    ira = (int *) tcalloc(3,sizeof(int)); ira[0]=s->rows; ira[1]=s->cols; ira[2]=s->length;
    multidradump(ra,3,ira,filename);
    tfree(ira);ira=NULL;
    tfree(ra);ra=NULL;}
}

void seidcorr_compile_helper(int verbose,int maxindex,int example_flag,int movie_flag)
{
  char *gs=GLOBAL_STRING,filename[1024],filename2[1024],command[2048];
  int *ira=NULL,ndim=0,nrows=0,ncols=0,nlength=0,nr=0,nc=0,nl=0,nd=0,tab=0;
  int found_flag=0,total_found=0;
  double *ra_movie=NULL,*ra=NULL,*ra_sum=NULL,*ra_cor=NULL,*ra_cor_sum=NULL,mean=0,stdev=0,max=0,min=0;
  double *ra_cor_rad_v=NULL,*ra_cor_rad_n=NULL,dist=0;
  int var_index=0,nr2=0,nc2=0,nr3=0,nc3=0,movie_rowoffset=0,movie_coloffset=0;
  double movie_mean=0,movie_std=0;
  double maxdia = 10000,xpos=0,ypos=0;
  char var_type[4];
  FILE *fp=NULL;
  if (verbose>1){ printf(" %% [entering seidcorr_compile_helper_helper]\n");}
  if (example_flag){
    if (verbose){ printf(" %% \n");}
    if (verbose){ printf(" %% starting example for seidcorr_compile_helper\n");}
    total_found = 4; nrows = 64, ncols = 64; nlength = 8;
    if (verbose){ printf(" %% first generating %d data files, each of dimension (%d,%d,%d)\n",total_found,nrows,ncols,nlength);}
    for (nd=0;nd<total_found;nd++){
      ra = (double *) tcalloc(nrows*ncols*nlength,sizeof(double));
      for (nr=0;nr<nrows;nr++){ for (nc=0;nc<ncols;nc++){ for (nl=0;nl<nlength;nl++){
	tab = nr+nc*nrows+nl*nrows*ncols;
	ra[tab] = sin(2*PI*(double)nl/(double)nlength + 2*PI*(double)nr/(double)nrows + 2*PI*(double)nd/(double)total_found)*cos(2*PI*(double)nl/(double)nlength + 2*PI*(double)nc/(double)ncols + 2*PI*(double)nd/(double)total_found);}}}
      ira = (int *) tcalloc(3,sizeof(int)); ira[0]=nrows; ira[1]=ncols; ira[2]=nlength;
      sprintf(filename,"./seidcorr_%s_%d_%sra","example",nd,"X");
      multidradump(ra,3,ira,filename);
      sprintf(filename,"./seidcorr_%s_%d_%sra.pnm","example",nd,"X");
      WritePNMfile_color(ra,nrows,ncols*nlength,0,0,filename,7);
      tfree(ira);ira=NULL;
      tfree(ra);ra=NULL;}
    if (verbose){ printf(" %% now attempting to read the data files and calculate the average over %d data files\n",total_found);}
    nd=0;found_flag=1;total_found=0;ira=NULL,ra=NULL;ra_sum=NULL;
    do{
      sprintf(filename,"./seidcorr_%s_%d_%sra","example",nd,"X");
      if (checktofind(filename)){ 
	if (verbose>1){ printf(" %% found filename %s\n",filename);}
	found_flag=1; total_found+=1;
	ndim=multidraread(filename,&ira,&ra); assert(ndim==3); nrows=ira[0];ncols=ira[1];nlength=ira[2];
	if (ra_sum==NULL){ ra_sum = (double *) tcalloc(nrows*ncols*nlength,sizeof(double));}
	raplusequals(ra_sum,nrows*ncols*nlength,ra);
	if (ira!=NULL){ tfree(ira);ira=NULL;} if (ra!=NULL){ tfree(ra);ra=NULL;}}
      else /* if (!checktofind(filename)) */{ if (verbose>1){ printf(" %% didn't find filename %s\n",filename);} found_flag=0;}
      nd+=1;}
    while ((maxindex>-1 && nd<=maxindex) || (maxindex<=-1 && found_flag));
    if (verbose){ printf(" %% found %d total files, final nrows %d ncols %d nlength %d\n",total_found,nrows,ncols,nlength);}
    sprintf(filename,"./seidcorr_%s_sum_%sra.pnm","example","X");
    WritePNMfile_color(ra_sum,nrows,ncols*nlength,0,0,filename,7);
    if (verbose){ printf(" %% now attempting to reread the data files, and compute correlation while subtracting off the average\n");}
    ratimesequals(ra_sum,nrows*ncols*nlength,-1.0/(double)total_found);
    if (verbose){ printf(" %% first: doing things the easy way\n");}
    ra_cor = tcalloc(nrows*ncols*nlength,sizeof(double));
    ra_cor_sum = tcalloc(nrows*ncols,sizeof(double));
    nd=0;found_flag=1;
    do{
      sprintf(filename,"./seidcorr_%s_%d_%sra","example",nd,"X");
      if (checktofind(filename)){ 
	if (verbose>1){ printf(" %% found filename %s\n",filename);}
	found_flag=1; 
	ndim=multidraread(filename,&ira,&ra); assert(ndim==3); nrows=ira[0];ncols=ira[1];nlength=ira[2];
	raplusequals(ra,nrows*ncols*nlength,ra_sum);
	for (nr=0;nr<nrows;nr++){ for (nc=0;nc<ncols;nc++){ 
	  tab=nr+nc*nrows; mean=0;stdev=0;
	  for (nl=0;nl<nlength;nl++){ mean+=ra[tab+nl*nrows*ncols];}
	  for (nl=0;nl<nlength;nl++){ stdev+=pow(ra[tab+nl*nrows*ncols],2);}
	  mean/=nlength; stdev/=nlength;
	  stdev = sqrt(maximum(0,stdev-pow(mean,2)));
	  for (nl=0;nl<nlength;nl++){ ra[tab+nl*nrows*ncols] = (ra[tab+nl*nrows*ncols]-mean)/(stdev<=0?1:stdev);}}}
	for (nl=0;nl<nlength;nl++){ 
	  tab=0+0*nrows+nl*nrows*ncols;
	  fftwconvolve(NULL,NULL,&(ra[tab]),&(ra[tab]),&(ra_cor[tab]),nrows,ncols);}
	sprintf(filename2,"./seidcorr_%s_%d_conv_%sra.pnm","example",nd,"X");
	WritePNMfile_color(ra_cor,nrows,ncols*nlength,1,-1,filename2,7);
	for (nl=0;nl<nlength;nl++){ raplusequals(ra_cor_sum,nrows*ncols,&(ra_cor[0+0*nrows+nl*nrows*ncols]));}
	if (ira!=NULL){ tfree(ira);ira=NULL;} if (ra!=NULL){ tfree(ra);ra=NULL;}}
      else /* if (!checktofind(filename)) */{ if (verbose>1){ printf(" %% didn't find filename %s\n",filename);} found_flag=0;}
      nd+=1;}
    while ((maxindex>-1 && nd<=maxindex) || (maxindex<=-1 && found_flag));
    ratimesequals(ra_cor_sum,nrows*ncols,1.0/(double)(total_found*nlength*nrows*ncols));
    stats("double",ra_cor_sum,nrows*ncols,&max,&min,NULL,NULL);
    if (verbose){ printf(" %% ra_cor_sum ranges from %f to %f\n",max,min);}
    sprintf(filename,"./seidcorr_%s_cor_%sra","example","X"); radump(ra_cor_sum,"double",nrows,ncols,filename);
    sprintf(filename,"./seidcorr_%s_cor_%sra.pnm","example","X"); WritePNMfile_color(ra_cor_sum,nrows,ncols,1,-1,filename,7);
    tfree(ra_cor);ra_cor=NULL;  
    tfree(ra_cor_sum);ra_cor_sum=NULL;
    if (verbose){ printf(" %% second: doing things the slow way\n");}
    ra_cor = tcalloc(nrows*ncols,sizeof(double));
    ra_cor_sum = tcalloc(nrows*ncols,sizeof(double));
    nd=0;found_flag=1;
    do{
      sprintf(filename,"./seidcorr_%s_%d_%sra","example",nd,"X");
      if (checktofind(filename)){ 
	if (verbose>1){ printf(" %% found filename %s\n",filename);}
	found_flag=1; 
	ndim=multidraread(filename,&ira,&ra); assert(ndim==3); nrows=ira[0];ncols=ira[1];nlength=ira[2];
	raplusequals(ra,nrows*ncols*nlength,ra_sum);
	sprintf(filename2,"./seidcorr_%s_%d_minus_avg_%sra.pnm","example",nd,"X");
	WritePNMfile_color(ra,nrows,ncols*nlength,0,0,filename2,7);
	for (nr=0;nr<nrows;nr++){ for (nc=0;nc<ncols;nc++){ 
	  tab=nr+nc*nrows; mean=0;stdev=0;
	  for (nl=0;nl<nlength;nl++){ mean+=ra[tab+nl*nrows*ncols];}
	  for (nl=0;nl<nlength;nl++){ stdev+=pow(ra[tab+nl*nrows*ncols],2);}
	  mean/=nlength; stdev/=nlength;
	  stdev = sqrt(maximum(0,stdev-pow(mean,2)));
	  for (nl=0;nl<nlength;nl++){ ra[tab+nl*nrows*ncols] = (ra[tab+nl*nrows*ncols]-mean)/(stdev<=0?1:stdev);}}}
	sprintf(filename2,"./seidcorr_%s_%d_minus_avg_mean0var1_%sra.pnm","example",nd,"X");
	WritePNMfile_color(ra,nrows,ncols*nlength,0,0,filename2,7);
	for (nr=0;nr<nrows;nr++){ for (nc=0;nc<ncols;nc++){
	  ra_cor[nr+nc*nrows]=0;
	  for (nr2=0;nr2<nrows;nr2++){ for (nc2=0;nc2<ncols;nc2++){
	    nr3 = periodize(nr2-nr,0,nrows); nc3 = periodize(nc2-nc,0,ncols);
	    for (nl=0;nl<nlength;nl++){
	      ra_cor[nr+nc*nrows] += ra[nr2+nc2*nrows+nl*nrows*ncols]*ra[nr3+nc3*nrows+nl*nrows*ncols];}}}}}
	ratimesequals(ra_cor,nrows*ncols,1.0/(double)(nrows*ncols*nlength));
	stats("double",ra_cor,nrows*ncols,&max,&min,NULL,NULL);
	if (verbose){ printf(" %% ra_cor %d ranges from %f to %f\n",nd,max,min);}
	sprintf(filename2,"./seidcorr_%s_%d_cor_slow_%sra.pnm","example",nd,"X");
	WritePNMfile_color(ra_cor,nrows,ncols,1,-1,filename2,7);
	raplusequals(ra_cor_sum,nrows*ncols,ra_cor);
	if (ira!=NULL){ tfree(ira);ira=NULL;} if (ra!=NULL){ tfree(ra);ra=NULL;}}
      else /* if (!checktofind(filename)) */{ if (verbose>1){ printf(" %% didn't find filename %s\n",filename);} found_flag=0;}
      nd+=1;}
    while ((maxindex>-1 && nd<=maxindex) || (maxindex<=-1 && found_flag));
    ratimesequals(ra_cor_sum,nrows*ncols,1.0/(double)(total_found));
    sprintf(filename,"./seidcorr_%s_cor_sum_slow_%sra","example","X"); radump(ra_cor_sum,"double",nrows,ncols,filename);
    sprintf(filename,"./seidcorr_%s_cor_sum_slow_%sra.pnm","example","X"); WritePNMfile_color(ra_cor_sum,nrows,ncols,1,-1,filename,7);
    tfree(ra_cor);ra_cor=NULL;  
    tfree(ra_cor_sum);ra_cor_sum=NULL;
    tfree(ra_sum);ra_sum=NULL;}
  else /* if (!example_flag) */{
    for (var_index=0;var_index<5;var_index++){
      switch (var_index){
      case 0: sprintf(var_type,"V"); movie_rowoffset=0; movie_coloffset=0; break;
      case 1: sprintf(var_type,"VS"); movie_rowoffset=0; movie_coloffset=1; break;
      case 2: sprintf(var_type,"sA"); movie_rowoffset=1; movie_coloffset=0; break;
      case 3: sprintf(var_type,"sN"); movie_rowoffset=1; movie_coloffset=1; break;
      case 4: sprintf(var_type,"sG"); movie_rowoffset=1; movie_coloffset=2; break;
      default: sprintf(var_type,"VS"); movie_rowoffset=0; movie_coloffset=1; break;}
      nd=0;found_flag=1;total_found=0;ira=NULL,ra=NULL;ra_sum=NULL;
      do{
	sprintf(filename,"./seidcorr_%s_%d_%sra",gs,nd,var_type);
	if (checktofind(filename)){ 
	  if (verbose>1){ printf(" %% found filename %s\n",filename);}
	  found_flag=1; total_found+=1;
	  ndim=multidraread(filename,&ira,&ra); assert(ndim==3); nrows=ira[0];ncols=ira[1];nlength=ira[2];
	  if (movie_flag){ if (ra_movie==NULL){ ra_movie = (double *) tcalloc(nrows*ncols*6*nlength,sizeof(double));}}
	  if (ra_sum==NULL){ ra_sum = (double *) tcalloc(nrows*ncols*nlength,sizeof(double));}
	  raplusequals(ra_sum,nrows*ncols*nlength,ra);
	  if (ira!=NULL){ tfree(ira);ira=NULL;} if (ra!=NULL){ tfree(ra);ra=NULL;}}
	else /* if (!checktofind(filename)) */{ if (verbose>1){ printf(" %% didn't find filename %s\n",filename);} found_flag=0;}
	nd+=1;}
      while ((maxindex>-1 && nd<=maxindex) || (maxindex<=-1 && found_flag));
      if (verbose){ printf(" %% found %d total files, final nrows %d ncols %d nlength %d\n",total_found,nrows,ncols,nlength);}
      if (movie_flag){ 
	stats("double",ra_sum,nrows*ncols*nlength,NULL,NULL,&movie_mean,&movie_std);
	for (nl=0;nl<nlength;nl++){ 
	  ra = (double *) tcalloc(nrows*ncols,sizeof(double));
	  raplugin(ra,nrows,ncols,&(ra_sum[0+0*nrows+nl*nrows*ncols]),nrows,ncols,0,0);
	  raaddequals(ra,nrows*ncols,-movie_mean); ratimesequals(ra,nrows*ncols,1.0/STD_VIEW/movie_std);
	  for (nr=0;nr<nrows*ncols;nr++){ ra[nr]=crop(ra[nr],-1,1);}
	  raplugin(&(ra_movie[0+0+nl*nrows*ncols*6]),nrows*2,ncols*3,ra,nrows,ncols,movie_rowoffset*nrows,movie_coloffset*ncols);
	  if (ra!=NULL){ tfree(ra);ra=NULL;}}}
      ra_cor = tcalloc(nrows*ncols*nlength,sizeof(double));
      ra_cor_sum = tcalloc(nrows*ncols,sizeof(double));
      nd=0;found_flag=1;
      do{
	sprintf(filename,"./seidcorr_%s_%d_%sra",gs,nd,var_type);
	if (checktofind(filename)){ 
	  if (verbose>1){ printf(" %% found filename %s\n",filename);}
	  found_flag=1; 
	  ndim=multidraread(filename,&ira,&ra); assert(ndim==3); nrows=ira[0];ncols=ira[1];nlength=ira[2];
	  if (verbose>2){
	    sprintf(filename,"./seidcorr_%s_%d_%sra.pnm",gs,nd,var_type);
	    WritePNMfile_color(ra,nrows,ncols*nlength,0,0,filename,7);}
	  raplusequals(ra,nrows*ncols*nlength,ra_sum);
	  for (nr=0;nr<nrows;nr++){ for (nc=0;nc<ncols;nc++){ 
	    tab=nr+nc*nrows; mean=0;stdev=0;
	    for (nl=0;nl<nlength;nl++){ mean+=ra[tab+nl*nrows*ncols];}
	    for (nl=0;nl<nlength;nl++){ stdev+=pow(ra[tab+nl*nrows*ncols],2);}
	    mean/=nlength; stdev/=nlength;
	    stdev = sqrt(maximum(0,stdev-pow(mean,2)));
	    for (nl=0;nl<nlength;nl++){ ra[tab+nl*nrows*ncols] = (ra[tab+nl*nrows*ncols]-mean)/(stdev<=0?1:stdev);}}}
	  for (nl=0;nl<nlength;nl++){ 
	    tab=0+0*nrows+nl*nrows*ncols;
	    fftwconvolve(NULL,NULL,&(ra[tab]),&(ra[tab]),&(ra_cor[tab]),nrows,ncols);}
	  for (nl=0;nl<nlength;nl++){ raplusequals(ra_cor_sum,nrows*ncols,&(ra_cor[0+0*nrows+nl*nrows*ncols]));}
	  if (ira!=NULL){ tfree(ira);ira=NULL;} if (ra!=NULL){ tfree(ra);ra=NULL;}}
	else /* if (!checktofind(filename)) */{ if (verbose>1){ printf(" %% didn't find filename %s\n",filename);} found_flag=0;}
	nd+=1;}
      while ((maxindex>-1 && nd<=maxindex) || (maxindex<=-1 && found_flag));
      ratimesequals(ra_cor_sum,nrows*ncols,1.0/(double)(total_found*nlength*nrows*ncols));
      sprintf(filename,"./seidcorr_%s_cor_%sra",gs,var_type); radump(ra_cor_sum,"double",nrows,ncols,filename);
      sprintf(filename,"./seidcorr_%s_cor_%sra.pnm",gs,var_type); WritePNMfile_color(ra_cor_sum,nrows,ncols,1,-1,filename,7);
      nd = (int)floor(sqrt(pow(nrows,2)+pow(ncols,2))/2)+2;
      ra_cor_rad_v = (double *) tcalloc(nd,sizeof(double)); ra_cor_rad_n = (double *) tcalloc(nd,sizeof(double));
      for (nr=0;nr<nrows;nr++){ for (nc=0;nc<ncols;nc++){
	nr2=periodize(nr,-nrows/2,nrows/2); nc2=periodize(nc,-ncols/2,ncols/2); dist = sqrt(pow(nr2,2)+pow(nc2,2));
	tab = maximum(0,minimum(nd-1,(int)floor(dist)));
	ra_cor_rad_v[tab] += ra_cor_sum[nr+nc*nrows]; ra_cor_rad_n[tab] += 1;}}
      sprintf(filename2,"seidcorr_%s_cor_rad_%s",gs,var_type);
      sprintf(filename,"%s.fig",filename2);
      if ((fp=fopen(filename,"w"))==NULL){ printf(" %% Warning, cannot create %s in seidcorr_compile_helper\n",filename); fp=stdout;}
      fprintf(fp,"%s",FIG_PREAMBLE);
      fprintf(fp,"2 1 0 5 %d 7 %d 0 -1 0.000 0 0 -1 0 0 %d\n\t",/*color*/0,/*depth*/50,/*npoints*/nd);
      for (nr=0;nr<nd;nr++){
	xpos = (double)(nr+0.5)/(double)nd; ypos = (double)(ra_cor_rad_v[nr]/=maximum(1,ra_cor_rad_n[nr]));
	printf("%f\n",ypos);
	fprintf(fp," %d %d",(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos));}
      fprintf(fp,"\n");
      fprintf(fp,"2 1 0 5 %d 7 %d 0 -1 0.000 0 0 -1 0 0 %d\n\t",/*color*/0,/*depth*/50,/*npoints*/4+1);
      fprintf(fp," %d %d",(int)floor(maxdia*+1),(int)maxdia-(int)floor(maxdia*+1));
      fprintf(fp," %d %d",(int)floor(maxdia*+1),(int)maxdia-(int)floor(maxdia*-1));
      fprintf(fp," %d %d",(int)floor(maxdia*0),(int)maxdia-(int)floor(maxdia*-1));
      fprintf(fp," %d %d",(int)floor(maxdia*0),(int)maxdia-(int)floor(maxdia*+1));
      fprintf(fp," %d %d",(int)floor(maxdia*+1),(int)maxdia-(int)floor(maxdia*+1));
      fprintf(fp,"\n");
      if (fp!=stdout){ fclose(fp);}
      sprintf(command,"fig2dev -Ljpeg -q %d %s.fig %s.jpg;",/*quality*/5,filename2,filename2); system(command);
      tfree(ra_cor);ra_cor=NULL;  
      tfree(ra_sum);ra_sum=NULL;
      tfree(ra_cor_sum);ra_cor_sum=NULL;
      tfree(ra_cor_rad_v);ra_cor_rad_v=NULL;
      tfree(ra_cor_rad_n);ra_cor_rad_n=NULL;}
    if (movie_flag && ra_movie!=NULL){ 
      for (nl=0;nl<nlength;nl++){ 
	sprintf(filename,"./seidcorr_%s_movie_frame.%d.pnm",gs,nl);
	WritePNMfile_color(&(ra_movie[0+0+nl*nrows*ncols*6]),nrows*2,ncols*3,1,-1,filename,7);}
      if (ra_movie!=NULL){ tfree(ra_movie);ra_movie=NULL;}
      if (verbose){ printf(" %% attempting to make movie out of seidcorr_%s_movie_frame.*.pnm --- should work on sphinx. requires /usr/bin/ppmtompeg/\n",gs);}
      sprintf(filename,"seidcorr_%s_movie.param",gs);
      if ((fp=fopen(filename,"w"))==NULL){ printf("cannot open %s in seidcorr_compile_helper, writing to stdout\n",filename); fp=stdout;}
      fprintf(fp,"PATTERN \t IBBPBBPBBPBBPBBP\n");
      fprintf(fp,"OUTPUT \t seidcorr_%s_movie.mpg\n",gs);
      fprintf(fp,"BASE_FILE_FORMAT \t PNM\n");
      fprintf(fp,"INPUT_CONVERT *\n");
      fprintf(fp,"GOP_SIZE \t 16\n");
      fprintf(fp,"SLICES_PER_FRAME \t 1\n");
      fprintf(fp,"INPUT_DIR \t ./\n");
      fprintf(fp,"INPUT\n");
      fprintf(fp,"seidcorr_%s_movie_frame.*.pnm [%d-%d]\n",gs,0,nlength-1);
      fprintf(fp,"END_INPUT\n");
      fprintf(fp,"PIXEL \t HALF\n");
      fprintf(fp,"RANGE \t 10\n");
      fprintf(fp,"PSEARCH_ALG \t LOGARITHMIC\n");
      fprintf(fp,"BSEARCH_ALG \t CROSS2\n");
      fprintf(fp,"IQSCALE \t 8\n");
      fprintf(fp,"PQSCALE \t 10\n");
      fprintf(fp,"BQSCALE \t 25\n");
      fprintf(fp,"REFERENCE_FRAME ORIGINAL\n");
      fprintf(fp,"BIT_RATE 1000000\n");
      fprintf(fp,"BUFFER_SIZE 327680\n");
      fprintf(fp,"FRAME_RATE 30\n");
      if (fp!=stdout){ fclose(fp);}
      sprintf(filename,"seidcorr_%s_movie_maker.sh",gs);
      if ((fp=fopen(filename,"w"))==NULL){ printf("cannot open %s in seidcorr_compile_helper, writing to stdout\n",filename); fp=stdout;}
      fprintf(fp,"/usr/bin/ppmtompeg ./seidcorr_%s_movie.param;\n",gs);
      for (nl=0;nl<nlength;nl++){ fprintf(fp,"rm ./seidcorr_%s_movie_frame.%d.pnm;\n",gs,nl);}
      if (fp!=stdout){ fclose(fp); sprintf(command,"nohup nice sh ./seidcorr_%s_movie_maker.sh & \n",gs); system(command);}}}
}

void seidcorr_compile(int argc,char **argv)
{
  int verbose=2;
  int still_have_options = argc-2;
  char helpfile[1024];
  int maxindex=-1;
  int example_flag=0;
  int movie_flag=0;
  if (verbose>1){ printf(" %% [entering seidcorr_compile]\n");}
  sprintf(helpfile,"%% [-Vverboselevel] [-Mmaxindex] [-X (example_flag=1)]\n");
  if (still_have_options==0){ printf(helpfile); exit(EXIT_SUCCESS);}
  while (still_have_options){
    switch(getopt(argc,argv,"V:M:X:m")){
    case 'V': 
      verbose = atoi(optarg);
      if (verbose){ printf(" %% verboselevel set to %d\n",verbose);}
      break;
    case 'M': 
      maxindex = maximum(-1,atoi(optarg));
      if (verbose){ printf(" %% maxindex set to %d\n",maxindex);}
      break;
    case 'X':
      example_flag = 1;
      if (verbose){ printf(" %% example_flag set to %d\n",example_flag);}
      break;
    case 'm':
      movie_flag = 1;
      if (verbose){ printf(" %% movie_flag set to %d\n",movie_flag);}
    default: printf(helpfile); exit(EXIT_SUCCESS); break;}
    still_have_options--;}
  seidcorr_compile_helper(verbose,maxindex,example_flag,movie_flag);
}
  
/* Here are the suite functions */

void cortex_data_tfree()
{
  int verbose=0;
  if (verbose){ printf(" %% [entering cortex_data_clear] will tfree all data structure\n");}
  if (CORTEX_BOTHER){
    if (RTC_BOTHER!=0 || GLOBAL_RTC!=NULL){ rtctfree(GLOBAL_RTC); GLOBAL_RTC=NULL; RTC_BOTHER=0;}
    if (STROBETRACE_BOTHER!=0 || GLOBAL_STROBETRACE!=NULL){ strobetracetfree(GLOBAL_STROBETRACE); GLOBAL_STROBETRACE=NULL; STROBETRACE_BOTHER=0;}
    if (TUNINGCURVE_BOTHER!=0 || GLOBAL_TUNINGCURVE!=NULL){ tuningcurvetfree(GLOBAL_TUNINGCURVE); GLOBAL_TUNINGCURVE=NULL; TUNINGCURVE_BOTHER=0;}
    if (LMITRI_BOTHER!=0 || GLOBAL_LMITRI!=NULL){ lmitritfree(GLOBAL_LMITRI); GLOBAL_LMITRI=NULL; LMITRI_BOTHER=0;}
    if (PTREE_BOTHER!=0 || GLOBAL_PTREE!=NULL){ ptreetfree(GLOBAL_PTREE); GLOBAL_PTREE=NULL; PTREE_BOTHER=0;}
    if (CLOSET_BOTHER!=0 || GLOBAL_CLOSET!=NULL){ closettfree(GLOBAL_CLOSET); GLOBAL_CLOSET=NULL; CLOSET_BOTHER=0;}
    if (POWER_BOTHER!=0 || GLOBAL_POWER!=NULL){ powertfree(GLOBAL_POWER); GLOBAL_POWER=NULL; POWER_BOTHER=0;}
    if (TAOF_BOTHER!=0 || GLOBAL_TAOF!=NULL){ taoftfree(GLOBAL_TAOF); GLOBAL_TAOF=NULL; TAOF_BOTHER=0;}
    if (SEIDCORR_BOTHER!=0 || GLOBAL_SEIDCORR!=NULL){ seidcorrtfree(GLOBAL_SEIDCORR); GLOBAL_SEIDCORR=NULL; SEIDCORR_BOTHER=0;}
    if (YGGDRASIL_BOTHER!=0 || GLOBAL_YGGDRASIL!=NULL){ yggdrasiltfree(GLOBAL_YGGDRASIL); GLOBAL_YGGDRASIL=NULL; YGGDRASIL_BOTHER=0;}
    if (BONSAI_BOTHER!=0 || GLOBAL_BONSAI!=NULL){ bonsaitfree(GLOBAL_BONSAI); GLOBAL_BONSAI=NULL; BONSAI_BOTHER=0;}
    if (HYDRA_BOTHER!=0 || GLOBAL_HYDRA!=NULL){ hydratfree(GLOBAL_HYDRA); GLOBAL_HYDRA=NULL; HYDRA_BOTHER=0;}
    if (LYAPUNOV_BOTHER!=0 || GLOBAL_LYAPUNOV!=NULL){ lyapunovtfree(GLOBAL_LYAPUNOV); GLOBAL_LYAPUNOV=NULL; LYAPUNOV_BOTHER=0;}}
}

void lgnswitch(struct lgn *p,double t,double DT)
{
  /* govern input and data-structures */
  int verbose=1;
  double time_0=0,time_1=0,time_2=0,time_delay=0,time_3=0,time_4=0,time_5=0,time_6=0,time_cycle=0,time_7=0,time_dumpevery=0,time_8=0,time_9=0,time_10=0;
  char **fnamebase=NULL;
  char *gs2=GLOBAL_STRING_2;
  char outputname[512];
  char command[1024];
  int nr=0,tindex=0,sindex=0,dindex=0;
  int nr1=0,nr2=0;
  double *timera=NULL;
  int nt=0,ntimes=0,nowtime=0;
  char tempchar[1024];
  int taof_flag=0;
  int suite_10_idslmi_flag=0;
  int suite_11_lmitri_flag=0;
  if (verbose>1){ printf(" %% [enter lgnswitch] with t %0.2f and DT %0.2f\n",t,DT);}
  switch (SUITE_BOTHER){
  case 1:
    time_0=0;
    time_1=time_0+1024*4;
    time_2=time_1+512;
    time_delay=16; /* fix for lmi */
    time_3=time_2+512;
    time_4=time_3+512;
    time_5=time_4+512;
    time_6=time_5+512;
    time_cycle=1024;
    time_7=time_6+time_cycle*SUITE_NSECONDS;
    time_dumpevery = SUITE_DUMPEVERY;
    time_8=time_7+1024*SUITE_NSECONDS;
    time_9=time_8+1024*SUITE_NSECONDS;
    time_10=time_9+1024*SUITE_NSECONDS;
    if (verbose>2){ printf(" times are: %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f \n",time_0,time_1,time_2,time_3,time_4,time_5,time_6,time_7,time_8,time_9,time_10);}
    if (t==0){ /* clear everything */
      if (verbose){ printf(" %% \n");}
      if (verbose){ printf(" %% entering first phase of suite one -- background\n");}
      GLOBAL_TF=0;
      if (GRATING_VS_LMI!=-10){ GRATING_VS_LMI = -10;}
      if (CORTEX_BOTHER){
	if (RTC_BOTHER!=0){ rtctfree(GLOBAL_RTC); GLOBAL_RTC=NULL; RTC_BOTHER=0;}
	if (STROBETRACE_BOTHER!=0){ strobetracetfree(GLOBAL_STROBETRACE); GLOBAL_STROBETRACE=NULL; STROBETRACE_BOTHER=0;}
	if (TUNINGCURVE_BOTHER!=0){ tuningcurvetfree(GLOBAL_TUNINGCURVE); GLOBAL_TUNINGCURVE=NULL; TUNINGCURVE_BOTHER=0;}
	if (LMITRI_BOTHER!=0){ lmitritfree(GLOBAL_LMITRI); GLOBAL_LMITRI=NULL; LMITRI_BOTHER=0;}
	if (PTREE_BOTHER!=0){ ptreetfree(GLOBAL_PTREE); GLOBAL_PTREE=NULL; PTREE_BOTHER=0;}
	if (CLOSET_BOTHER!=0){ closettfree(GLOBAL_CLOSET); GLOBAL_CLOSET=NULL; CLOSET_BOTHER=0;}
	if (POWER_BOTHER!=0){ powertfree(GLOBAL_POWER); GLOBAL_POWER=NULL; POWER_BOTHER=0;}
	if (YGGDRASIL_BOTHER!=0){ yggdrasiltfree(GLOBAL_YGGDRASIL); GLOBAL_YGGDRASIL=NULL; YGGDRASIL_BOTHER=0;}
	if (BONSAI_BOTHER!=0){ bonsaitfree(GLOBAL_BONSAI); GLOBAL_BONSAI=NULL; BONSAI_BOTHER=0;}
	if (HYDRA_BOTHER!=0){ hydratfree(GLOBAL_HYDRA); GLOBAL_HYDRA=NULL; HYDRA_BOTHER=0;}
	if (LYAPUNOV_BOTHER!=0){ lyapunovtfree(GLOBAL_LYAPUNOV); GLOBAL_LYAPUNOV=NULL; LYAPUNOV_BOTHER=0;}}}
    else if (t>=0 && t<time_1){ /* background */
      if (GRATING_VS_LMI!=-10){ GRATING_VS_LMI = -10;}
      if (CORTEX_BOTHER){
	if (PTREE_BOTHER==0){
	  if (verbose){ printf(" %% %% making ptree\n");}
	  PTREE_BOTHER=1; GLOBAL_PTREE_BITBYBIT=0; GLOBAL_PTREE_REGION_TYPE=3;
	  GLOBAL_PTREE = ptreemake(GLOBAL_PTREE_NREGIONS,GLOBAL_PTREE_EVENT_WITHIN,GLOBAL_PTREE_EVENT_THRESHOLD,GLOBAL_PTREE_REGION_TYPE,GLOBAL_PTREE_NLEGS,GLOBAL_PTREE_LEGTIME);}
	if (STROBETRACE_BOTHER==0){
	  if (verbose){ printf(" %% %% making strobetrace\n");}
	  STROBETRACE_BOTHER=1;
	  GLOBAL_STROBETRACE = strobetracemake(nget(GLOBAL_Nra,0,0),GLOBAL_STROBETRACE_NANGLES,time_1-time_0,1.0,0,GLOBAL_STROBETRACE_AVALANCHE_UPDATE_EVERY,GLOBAL_STROBETRACE_AVALANCHE_LOGDMIN,GLOBAL_STROBETRACE_AVALANCHE_LOGDMAX);}}}
    else if (t>=time_1 && t<time_2){ /* line motion illusion */
      if (GRATING_VS_LMI!=1){ 
	if (verbose){ printf(" %% entering second phase of suite one -- line motion illusion\n");}
	GRATING_VS_LMI = 1; STIMULUS_ONSET_TIME = time_1+time_delay;}
      if (CORTEX_BOTHER){
	if (PTREE_BOTHER!=0){
	  if (verbose){ printf(" %% %% dumping ptree\n");}
	  ptreerate(GLOBAL_PTREE); 
	  ptreedump_starter(GLOBAL_PTREE,"background",2,0,0,0,+1,-1); 
	  ptreetfree(GLOBAL_PTREE); GLOBAL_PTREE=NULL; PTREE_BOTHER=0;}
	if (STROBETRACE_BOTHER!=0){ 
	  if (verbose){ printf(" %% %% dumping strobetrace\n");}
	  strobetracedump(GLOBAL_STROBETRACE,0); 
	  strobetracedump(GLOBAL_STROBETRACE,1); 
	  strobetracetfree(GLOBAL_STROBETRACE); GLOBAL_STROBETRACE=NULL; STROBETRACE_BOTHER=0;}
	if (LMITRI_BOTHER==0){ 
	  if (verbose){ printf(" %% %% making lmitri\n");}
	  LMITRI_BOTHER=1;
	  GLOBAL_LMITRI = lmitrimake(STIMULUS_ONSET_TIME,GLOBAL_LMITRI_TIMELENGTH,GLOBAL_LMITRI_ROW_MAX,GLOBAL_LMITRI_ROW_MIN);}}}
    else if (t>=time_2 && t<time_3){ /* another line motion illusion */
      if (STIMULUS_ONSET_TIME!=time_2+time_delay){
	if (verbose){ printf(" %% entering third phase of suite one -- line motion illusion\n");}
	STIMULUS_ONSET_TIME = time_2+time_delay;}
      if (CORTEX_BOTHER){ 
	if (GLOBAL_LMITRI->time_start!=STIMULUS_ONSET_TIME){
	  GLOBAL_LMITRI->time_start = STIMULUS_ONSET_TIME;}}}
    else if (t>=time_3 && t<time_4){ /* yet another line motion illusion */
      if (STIMULUS_ONSET_TIME!=time_3+time_delay){
	if (verbose){ printf(" %% entering fourth phase of suite one -- line motion illusion\n");}
	STIMULUS_ONSET_TIME = time_3+time_delay;}
      if (CORTEX_BOTHER){ 
	if (GLOBAL_LMITRI->time_start!=STIMULUS_ONSET_TIME){
	  GLOBAL_LMITRI->time_start = STIMULUS_ONSET_TIME;}}}
    else if (t>=time_4 && t<time_5){ /* and another line motion illusion */
      if (STIMULUS_ONSET_TIME!=time_4+time_delay){
	if (verbose){ printf(" %% entering fifth phase of suite one -- line motion illusion\n");}
	STIMULUS_ONSET_TIME = time_4+time_delay;}
      if (CORTEX_BOTHER){ 
	if (GLOBAL_LMITRI->time_start!=STIMULUS_ONSET_TIME){
	  GLOBAL_LMITRI->time_start = STIMULUS_ONSET_TIME;}}}
    else if (t>=time_5 && t<time_6){ /* and a final line motion illusion */
      if (STIMULUS_ONSET_TIME!=time_5+time_delay){
	if (verbose){ printf(" %% entering sixth phase of suite one -- line motion illusion\n");}
	STIMULUS_ONSET_TIME = time_5+time_delay;}
      if (CORTEX_BOTHER){ 
	if (GLOBAL_LMITRI->time_start!=STIMULUS_ONSET_TIME){
	  GLOBAL_LMITRI->time_start = STIMULUS_ONSET_TIME;}}}
    else if (t>=time_6 && t<time_7){ /* dahump pulsed grating */
      STIMULUS_ONSET_TIME = time_6;
      if (GRATING_VS_LMI!=-1){ 
	if (verbose){ printf(" %% entering seventh phase of suite one -- pulsed grating\n");}
	GRATING_VS_LMI = -1; INPUT_SPACEANGLE= 0 /* dahump angle (vs rtc) */; CYCLE_LENGTH = time_cycle;}
      if (CORTEX_BOTHER){
	if (LMITRI_BOTHER!=0){
	  if (verbose){ printf(" %% %% dumping lmitri\n");}
	  lmitridump(GLOBAL_LMITRI,1); lmitritfree(GLOBAL_LMITRI); GLOBAL_LMITRI=NULL; LMITRI_BOTHER=0;}
	if (PTREE_BOTHER==0){
	  if (verbose){ printf(" %% %% making ptree\n");}
	  PTREE_BOTHER=1; GLOBAL_PTREE_BITBYBIT=0; GLOBAL_PTREE_REGION_TYPE=3;
	  GLOBAL_PTREE = ptreemake(GLOBAL_PTREE_NREGIONS,GLOBAL_PTREE_EVENT_WITHIN,GLOBAL_PTREE_EVENT_THRESHOLD,GLOBAL_PTREE_REGION_TYPE,GLOBAL_PTREE_NLEGS,GLOBAL_PTREE_LEGTIME);}
	if (STROBETRACE_BOTHER==0){
	  if (verbose){ printf(" %% %% making strobetrace\n");}	  
	  STROBETRACE_BOTHER=1;
	  GLOBAL_STROBETRACE = strobetracemake(nget(GLOBAL_Nra,0,0),GLOBAL_STROBETRACE_NANGLES,time_cycle,1.0,1,GLOBAL_STROBETRACE_AVALANCHE_UPDATE_EVERY,GLOBAL_STROBETRACE_AVALANCHE_LOGDMIN,GLOBAL_STROBETRACE_AVALANCHE_LOGDMAX);}}}
    else if (t>=time_7 && t<time_8){ /* tuningcurve1 */
      STIMULUS_ONSET_TIME = time_7; OUTPUT_DUMP_EVERY = time_dumpevery; /* start dumping data */
      if (GRATING_VS_LMI!=0){ 
	if (verbose){ printf(" %% entering eigth phase of suite one -- tuningcurve1\n");}
	GRATING_VS_LMI = 0; INPUT_SPACEANGLE = 0;}
      if (CORTEX_BOTHER){
	if (GLOBAL_PTREE_BITBYBIT==0 && PTREE_BOTHER!=0){ 
	  if (verbose){ printf(" %% %% dumping ptree\n");}
	  ptreerate(GLOBAL_PTREE); 
	  ptreedump_starter(GLOBAL_PTREE,"dahump",2,0,0,0,+1,-1); 
	  ptreetfree(GLOBAL_PTREE); GLOBAL_PTREE=NULL; PTREE_BOTHER=0;
	  if (verbose){ printf(" %% %% making ptree\n");}
	  PTREE_BOTHER=1;GLOBAL_PTREE_BITBYBIT=1; GLOBAL_PTREE_REGION_TYPE=2;
	  GLOBAL_PTREE = ptreemake(GLOBAL_PTREE_NREGIONS,GLOBAL_PTREE_EVENT_WITHIN,GLOBAL_PTREE_EVENT_THRESHOLD,GLOBAL_PTREE_REGION_TYPE,GLOBAL_PTREE_NLEGS,GLOBAL_PTREE_LEGTIME);}
	if (STROBETRACE_BOTHER!=0){ 
	  if (verbose){ printf(" %% %% dumping strobetrace\n");}
	  strobetracedump(GLOBAL_STROBETRACE,0); 
	  strobetracedump(GLOBAL_STROBETRACE,1); 
	  strobetracetfree(GLOBAL_STROBETRACE); GLOBAL_STROBETRACE=NULL; STROBETRACE_BOTHER=0;}
	if (TUNINGCURVE_BOTHER==0){
	  if (verbose){ printf(" %% %% making tuningcurve\n");}
	  TUNINGCURVE_BOTHER=1;
	  GLOBAL_TUNINGCURVE = tuningcurvemake(GLOBAL_TUNINGCURVE_NANGLES,GLOBAL_TUNINGCURVE_NRADIUS);}	
	if (OUTPUT_DUMP_EVERY>0 && (int)floor((GLOBAL_time-time_7)/OUTPUT_DUMP_EVERY)<(int)floor((GLOBAL_time-time_7+GLOBAL_DT)/OUTPUT_DUMP_EVERY)){
	  if (PTREE_BOTHER){
	    if (verbose){ printf(" %% %% dumping ptree\n");}
	    ptreerate(GLOBAL_PTREE);
	    ptreedump_starter(GLOBAL_PTREE,"tuningcurve1",2,0,0,0,+1,-1);
	    if (GLOBAL_PTREE_BITBYBIT){
	      if (verbose){ printf(" %% %% resetting ptree\n");}
	      ptreereset(GLOBAL_PTREE);}}}}}
    else if (t>=time_8 && t<time_9){ /* rtc experiment */
      STIMULUS_ONSET_TIME = time_8;
      if (GRATING_VS_LMI!=-2){ 
	if (verbose){ printf(" %% entering ninth phase of suite one -- rtc\n");}
	GRATING_VS_LMI = -2; OUTPUT_DUMP_EVERY = 0;}
      if (CORTEX_BOTHER){
	if (GLOBAL_PTREE_BITBYBIT==1 && PTREE_BOTHER!=0){ 
	  if (verbose){ printf(" %% %% dumping ptree\n");}
	  ptreerate(GLOBAL_PTREE); 
	  ptreedump_starter(GLOBAL_PTREE,"tuningcurve1",2,0,0,0,+1,-1); 
	  ptreetfree(GLOBAL_PTREE); GLOBAL_PTREE=NULL; PTREE_BOTHER=0;
	  if (verbose){ printf(" %% %% making ptree\n");}
	  PTREE_BOTHER=1;GLOBAL_PTREE_BITBYBIT=0; GLOBAL_PTREE_REGION_TYPE=3;
	  GLOBAL_PTREE = ptreemake(GLOBAL_PTREE_NREGIONS,GLOBAL_PTREE_EVENT_WITHIN,GLOBAL_PTREE_EVENT_THRESHOLD,GLOBAL_PTREE_REGION_TYPE,GLOBAL_PTREE_NLEGS,GLOBAL_PTREE_LEGTIME);}
	if (TUNINGCURVE_BOTHER!=0){ 
	  if (verbose){ printf(" %% %% dumping tuningcurve\n");}
	  tuningcurvedump(GLOBAL_TUNINGCURVE,1); tuningcurvetfree(GLOBAL_TUNINGCURVE); GLOBAL_TUNINGCURVE=NULL; TUNINGCURVE_BOTHER=0;}
	if (RTC_BOTHER==0){
	  if (verbose){ printf(" %% %% making rtc\n");}
	  RTC_BOTHER=1;
	  GLOBAL_RTC = rtcmake(GLOBAL_RTC_LENGTH,GLOBAL_RTC_FRAMELENGTH,GLOBAL_RTC_NANGLES,GLOBAL_RTC_NPHASES);}}}
    else if (t>=time_9 && t<time_10){ /* tuningcurve2 */
      STIMULUS_ONSET_TIME = time_9; OUTPUT_DUMP_EVERY = time_dumpevery; /* start dumping data */
      if (GRATING_VS_LMI!=0){ 
	if (verbose){ printf(" %% %% entering tenth phase of suite one -- tuningcurve2\n");}
	GRATING_VS_LMI = 0; INPUT_SPACEANGLE = 0 + INPUT_SPACEANGLE_BACON; /* six degrees of separation */}
      if (CORTEX_BOTHER){
	if (GLOBAL_PTREE_BITBYBIT==0 && PTREE_BOTHER!=0){
	  if (verbose){ printf(" %% %% dumping ptree\n");}
	  ptreerate(GLOBAL_PTREE);
	  ptreedump_starter(GLOBAL_PTREE,"rtc",2,0,0,0,+1,-1);
	  ptreetfree(GLOBAL_PTREE); GLOBAL_PTREE=NULL; PTREE_BOTHER=0;
	  if (verbose){ printf(" %% %% making ptree\n");}
	  PTREE_BOTHER=1;GLOBAL_PTREE_BITBYBIT=1; GLOBAL_PTREE_REGION_TYPE=2;
	  GLOBAL_PTREE = ptreemake(GLOBAL_PTREE_NREGIONS,GLOBAL_PTREE_EVENT_WITHIN,GLOBAL_PTREE_EVENT_THRESHOLD,GLOBAL_PTREE_REGION_TYPE,GLOBAL_PTREE_NLEGS,GLOBAL_PTREE_LEGTIME);}
	if (RTC_BOTHER!=0){ 
	  if (verbose){ printf(" %% %% dumping rtc\n");}
	  rtcdump(GLOBAL_RTC,0); 
	  rtcdump(GLOBAL_RTC,1); 
	  rtctfree(GLOBAL_RTC); GLOBAL_RTC=NULL; RTC_BOTHER=0;}
	if (OUTPUT_DUMP_EVERY>0 && (int)floor((GLOBAL_time-time_9)/OUTPUT_DUMP_EVERY)<(int)floor((GLOBAL_time-time_9+GLOBAL_DT)/OUTPUT_DUMP_EVERY)){
	  if (PTREE_BOTHER){
	    if (verbose){ printf(" %% dumping ptree\n");}
	    ptreerate(GLOBAL_PTREE);
	    ptreedump_starter(GLOBAL_PTREE,"tuningcurve2",2,0,0,0,+1,-1);
	    if (GLOBAL_PTREE_BITBYBIT){
	      if (verbose){ printf(" %% %% resetting ptree\n");}
	      ptreereset(GLOBAL_PTREE);}}}}}
    else if (t>=time_10){ /* finish */
      if (STIMULUS_ONSET_TIME!=time_10){
	if (verbose){ printf(" %% entering eleventh phase of suite one -- finishing\n");}
	STIMULUS_ONSET_TIME = time_10;
	GRATING_VS_LMI = -10; GLOBAL_TF = time_10+1;
	if (CORTEX_BOTHER){
	  if (GLOBAL_PTREE_BITBYBIT==1 && PTREE_BOTHER!=0){ 
	    if (verbose){ printf(" %% %% dumping ptree\n");}
	    ptreerate(GLOBAL_PTREE);
	    ptreedump_starter(GLOBAL_PTREE,"tuningcurve2",2,0,0,0,+1,-1);
	    ptreetfree(GLOBAL_PTREE); GLOBAL_PTREE=NULL; PTREE_BOTHER=0;}}
	fnamebase = (char **) tcalloc(2,sizeof(char *));
	fnamebase[0] = (char *) tcalloc(256,sizeof(char)); sprintf(fnamebase[0],"ptree_tuningcurve1_%srecord",gs2);
	fnamebase[1] = (char *) tcalloc(256,sizeof(char)); sprintf(fnamebase[1],"ptree_tuningcurve2_%srecord",gs2);
	for (tindex=1;tindex<=SUITE_TINDEXMAX;tindex++){ for (sindex=0;sindex<=SUITE_SINDEXMAX;sindex++){ for (dindex=0;dindex<=SUITE_DINDEXMAX;dindex++){
	  sprintf(outputname,"a12"); ptree_trialaverage_helper(2,GLOBAL_PTREE_NLEGS+1,tindex,1,1,0,sindex,dindex,fnamebase,outputname);}}}
	tfree(fnamebase[0]);fnamebase[0]=NULL;
	tfree(fnamebase[1]);fnamebase[1]=NULL;
	tfree(fnamebase);fnamebase=NULL;
	if (SUITE_BITBYBIT_REMOVE){ 
	  for (nr=0;nr<(int)floor((double)SUITE_NSECONDS*1024/(double)SUITE_DUMPEVERY)+1;nr++){
	    sprintf(command,"""rm"" ./ptree_tuningcurve1_%srecord_%d ;",gs2,nr); system(command);
	    sprintf(command,"""rm"" ./ptree_tuningcurve2_%srecord_%d ;",gs2,nr); system(command);}}}}
    break;
  case 2: 
    time_0=0;
    time_1=time_0+1024*4;
    time_2=time_1+1024*5;
    time_3=time_2+1024*6;
    time_4=time_3+1024*7;
    time_5=time_4+1024*8;
    time_6=time_5+1024*8;
    time_7=time_6+1024*64;
    time_8=time_6+1024*64;
    time_dumpevery=SUITE_DUMPEVERY;
    if (t==0){ /* clear everything */
      if (verbose){ printf(" %% \n");}
      if (verbose){ printf(" %% entering first phase of suite two -- background\n");}
      GLOBAL_TF=0;
      if (GRATING_VS_LMI!=-10){ GRATING_VS_LMI = -10;}
      if (CORTEX_BOTHER){
	if (RTC_BOTHER!=0){ rtctfree(GLOBAL_RTC); GLOBAL_RTC=NULL; RTC_BOTHER=0;}
	if (STROBETRACE_BOTHER!=0){ strobetracetfree(GLOBAL_STROBETRACE); GLOBAL_STROBETRACE=NULL; STROBETRACE_BOTHER=0;}
	if (TUNINGCURVE_BOTHER!=0){ tuningcurvetfree(GLOBAL_TUNINGCURVE); GLOBAL_TUNINGCURVE=NULL; TUNINGCURVE_BOTHER=0;}
	if (LMITRI_BOTHER!=0){ lmitritfree(GLOBAL_LMITRI); GLOBAL_LMITRI=NULL; LMITRI_BOTHER=0;}
	if (PTREE_BOTHER!=0){ ptreetfree(GLOBAL_PTREE); GLOBAL_PTREE=NULL; PTREE_BOTHER=0;}
	if (CLOSET_BOTHER!=0){ closettfree(GLOBAL_CLOSET); GLOBAL_CLOSET=NULL; CLOSET_BOTHER=0;}
	if (POWER_BOTHER!=0){ powertfree(GLOBAL_POWER); GLOBAL_POWER=NULL; POWER_BOTHER=0;}
	if (YGGDRASIL_BOTHER!=0){ yggdrasiltfree(GLOBAL_YGGDRASIL); GLOBAL_YGGDRASIL=NULL; YGGDRASIL_BOTHER=0;}
	if (BONSAI_BOTHER!=0){ bonsaitfree(GLOBAL_BONSAI); GLOBAL_BONSAI=NULL; BONSAI_BOTHER=0;}
	if (HYDRA_BOTHER!=0){ hydratfree(GLOBAL_HYDRA); GLOBAL_HYDRA=NULL; HYDRA_BOTHER=0;}
	if (LYAPUNOV_BOTHER!=0){ lyapunovtfree(GLOBAL_LYAPUNOV); GLOBAL_LYAPUNOV=NULL; LYAPUNOV_BOTHER=0;}}}
    else if (t>=0 && t<time_1){ /* background */
      if (GRATING_VS_LMI!=-10){ GRATING_VS_LMI = -10;}
      if (CORTEX_BOTHER){
	if (CLOSET_BOTHER==0){
	  if (verbose){ printf(" %% %% making closet\n");}
	  CLOSET_BOTHER=1; GLOBAL_CLOSET=closetmake(0,PIE_ROW_DIA,0,PIE_COL_DIA,GLOBAL_PTREE_LEGTIME,0.1);}}}
    else if (t>=time_1 && t<time_2){ /* background */
      if (GRATING_VS_LMI!=-10){ GRATING_VS_LMI = -10;}
      if (CORTEX_BOTHER){
	if (CLOSET_BOTHER==1){ if (GLOBAL_CLOSET->nlegs<1){       
	  if (verbose){ printf(" %% entering second phase of suite two -- background\n");}
	  ptreerate(GLOBAL_CLOSET->p);
	  closetremake(GLOBAL_CLOSET,1-GLOBAL_CLOSET->nlegs);}}}}
    else if (t>=time_2 && t<time_3){ /* background */
      if (GRATING_VS_LMI!=-10){ GRATING_VS_LMI = -10;}
      if (CORTEX_BOTHER){
	if (CLOSET_BOTHER==1){ if (GLOBAL_CLOSET->nlegs<2){       
	  if (verbose){ printf(" %% entering third phase of suite two -- background\n");}
	  ptreerate(GLOBAL_CLOSET->p);
	  closetremake(GLOBAL_CLOSET,2-GLOBAL_CLOSET->nlegs);}}}}
    else if (t>=time_3 && t<time_4){ /* background */
      if (GRATING_VS_LMI!=-10){ GRATING_VS_LMI = -10;}
      if (CORTEX_BOTHER){
	if (CLOSET_BOTHER==1){ if (GLOBAL_CLOSET->nlegs<3){ 
	  if (verbose){ printf(" %% entering fourth phase of suite two -- background\n");}
	  ptreerate(GLOBAL_CLOSET->p);
	  closetremake(GLOBAL_CLOSET,3-GLOBAL_CLOSET->nlegs);}}}}
    else if (t>=time_4 && t<time_5){ /* background */
      if (GRATING_VS_LMI!=-10){ GRATING_VS_LMI = -10;}
      if (CORTEX_BOTHER){
	if (CLOSET_BOTHER==1){ if (GLOBAL_CLOSET->nlegs<4){ 
	  if (verbose){ printf(" %% entering fifth phase of suite two -- background\n");}
	  ptreerate(GLOBAL_CLOSET->p);
	  closetremake(GLOBAL_CLOSET,4-GLOBAL_CLOSET->nlegs);}}}}
    else if (t>=time_5 && t<time_6){ /* background */
      if (GRATING_VS_LMI!=-10){ GRATING_VS_LMI = -10;}
      if (CORTEX_BOTHER){
	if (CLOSET_BOTHER==1){ if (GLOBAL_CLOSET->nlegs<5){ 
	  if (verbose){ printf(" %% entering sixth phase of suite two -- background\n");}
	  ptreerate(GLOBAL_CLOSET->p);
	  closetremake(GLOBAL_CLOSET,5-GLOBAL_CLOSET->nlegs);}}}}
    else if (t>=time_6 && t<time_7){ /* tuningcurve1 */
      STIMULUS_ONSET_TIME = time_7; OUTPUT_DUMP_EVERY = time_dumpevery; /* start dumping data */
      if (GRATING_VS_LMI!=0){ 
	if (verbose){ printf(" %% entering seventh phase of suite two -- tuningcurve1\n");}
	GRATING_VS_LMI = 0; INPUT_SPACEANGLE = 0;}
      if (CORTEX_BOTHER){
	if (CLOSET_BOTHER==1){ if (GLOBAL_CLOSET->nlegs<6){ 
	  ptreerate(GLOBAL_CLOSET->p);
	  closetremake(GLOBAL_CLOSET,6-GLOBAL_CLOSET->nlegs);}}
	if (GLOBAL_PTREE_BITBYBIT==0){ GLOBAL_PTREE_BITBYBIT=1;}
	if (OUTPUT_DUMP_EVERY>0 && (int)floor((GLOBAL_time-time_6)/OUTPUT_DUMP_EVERY)<(int)floor((GLOBAL_time-time_6+GLOBAL_DT)/OUTPUT_DUMP_EVERY)){
	  if (CLOSET_BOTHER){
	    if (verbose){ printf(" %% %% dumping closet\n");}
	    ptreerate(GLOBAL_CLOSET->p);
	    ptreedump_starter(GLOBAL_CLOSET->p,"tuningcurve1",2,0,0,0,+1,-1);
	    if (GLOBAL_PTREE_BITBYBIT){
	      if (verbose){ printf(" %% %% remaking closet\n");}
	      closetremake(GLOBAL_CLOSET,0);}}}}}
    else if (t>=time_7 && t<time_8){ /* tuningcurve2 */
      if (STIMULUS_ONSET_TIME!=time_8){ 
	if (verbose){ printf(" %% entering eigth phase of suite two -- tuningcurve2\n");}
	STIMULUS_ONSET_TIME = time_8;
	if (verbose){ printf(" %% remaking closet\n");}
	if (CLOSET_BOTHER==1){ closetremake(GLOBAL_CLOSET,0);}}
      OUTPUT_DUMP_EVERY = time_dumpevery; /* start dumping data */
      if (GRATING_VS_LMI!=0){ 
	GRATING_VS_LMI = 0; INPUT_SPACEANGLE = 0 + INPUT_SPACEANGLE_BACON; /* six degrees of separation */}
      if (CORTEX_BOTHER){
	if (CLOSET_BOTHER==1){ if (GLOBAL_CLOSET->nlegs<6){ 
	  ptreerate(GLOBAL_CLOSET->p);
	  closetremake(GLOBAL_CLOSET,6-GLOBAL_CLOSET->nlegs);}}
	if (GLOBAL_PTREE_BITBYBIT==0){ GLOBAL_PTREE_BITBYBIT=1;}
	if (OUTPUT_DUMP_EVERY>0 && (int)floor((GLOBAL_time-time_7)/OUTPUT_DUMP_EVERY)<(int)floor((GLOBAL_time-time_7+GLOBAL_DT)/OUTPUT_DUMP_EVERY)){
	  if (CLOSET_BOTHER){
	    if (verbose){ printf(" %% %% dumping closet\n");}
	    ptreerate(GLOBAL_CLOSET->p);
	    ptreedump_starter(GLOBAL_CLOSET->p,"tuningcurve2",2,0,0,0,+1,-1);
	    if (GLOBAL_PTREE_BITBYBIT){
	      if (verbose){ printf(" %% %% remaking closet\n");}
	      closetremake(GLOBAL_CLOSET,0);}}}}}
    else if (t>=time_8){ /* finish */
      GLOBAL_TF = time_8+1;
      if (CORTEX_BOTHER){
	if (CLOSET_BOTHER!=0){
	  if (verbose){ printf(" %% entering ninth phase of suite two -- finishing\n");}
	  closettfree(GLOBAL_CLOSET); GLOBAL_CLOSET=NULL; CLOSET_BOTHER=0;}}}
    break;
  case 3:
    time_0=0;
    time_1=time_0+1024*4;
    time_2=time_1+512;
    time_delay=16; /* fix for lmi */
    time_3=time_2+512;
    time_4=time_3+512;
    time_5=time_4+512;
    time_6=time_5+512;
    if (verbose>2){ printf(" times are: %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f\n",time_0,time_1,time_2,time_3,time_4,time_5,time_6);}
    if (t==0){ /* clear everything */
      if (verbose){ printf(" %% \n");}
      if (verbose){ printf(" %% entering first phase of suite three -- background\n");}
      GLOBAL_TF=0;
      if (GRATING_VS_LMI!=-10){ GRATING_VS_LMI = -10;}
      if (CORTEX_BOTHER){
	if (RTC_BOTHER!=0){ rtctfree(GLOBAL_RTC); GLOBAL_RTC=NULL; RTC_BOTHER=0;}
	if (STROBETRACE_BOTHER!=0){ strobetracetfree(GLOBAL_STROBETRACE); GLOBAL_STROBETRACE=NULL; STROBETRACE_BOTHER=0;}
	if (TUNINGCURVE_BOTHER!=0){ tuningcurvetfree(GLOBAL_TUNINGCURVE); GLOBAL_TUNINGCURVE=NULL; TUNINGCURVE_BOTHER=0;}
	if (LMITRI_BOTHER!=0){ lmitritfree(GLOBAL_LMITRI); GLOBAL_LMITRI=NULL; LMITRI_BOTHER=0;}
	if (PTREE_BOTHER!=0){ ptreetfree(GLOBAL_PTREE); GLOBAL_PTREE=NULL; PTREE_BOTHER=0;}
	if (CLOSET_BOTHER!=0){ closettfree(GLOBAL_CLOSET); GLOBAL_CLOSET=NULL; CLOSET_BOTHER=0;}
	if (POWER_BOTHER!=0){ powertfree(GLOBAL_POWER); GLOBAL_POWER=NULL; POWER_BOTHER=0;}
	if (YGGDRASIL_BOTHER!=0){ yggdrasiltfree(GLOBAL_YGGDRASIL); GLOBAL_YGGDRASIL=NULL; YGGDRASIL_BOTHER=0;}
	if (BONSAI_BOTHER!=0){ bonsaitfree(GLOBAL_BONSAI); GLOBAL_BONSAI=NULL; BONSAI_BOTHER=0;}
	if (HYDRA_BOTHER!=0){ hydratfree(GLOBAL_HYDRA); GLOBAL_HYDRA=NULL; HYDRA_BOTHER=0;}
	if (LYAPUNOV_BOTHER!=0){ lyapunovtfree(GLOBAL_LYAPUNOV); GLOBAL_LYAPUNOV=NULL; LYAPUNOV_BOTHER=0;}}}
    else if (t>=0 && t<time_1){ /* background */
      if (GRATING_VS_LMI!=-10){ GRATING_VS_LMI = -10;}
      if (CORTEX_BOTHER){
	if (PTREE_BOTHER==0){
	  if (verbose){ printf(" %% %% making ptree\n");}
	  PTREE_BOTHER=1; GLOBAL_PTREE_BITBYBIT=0; GLOBAL_PTREE_REGION_TYPE=3;
	  GLOBAL_PTREE = ptreemake(GLOBAL_PTREE_NREGIONS,GLOBAL_PTREE_EVENT_WITHIN,GLOBAL_PTREE_EVENT_THRESHOLD,GLOBAL_PTREE_REGION_TYPE,GLOBAL_PTREE_NLEGS,GLOBAL_PTREE_LEGTIME);}
	if (STROBETRACE_BOTHER==0){
	  if (verbose){ printf(" %% %% making strobetrace\n");}
	  STROBETRACE_BOTHER=1;
	  GLOBAL_STROBETRACE = strobetracemake(nget(GLOBAL_Nra,0,0),GLOBAL_STROBETRACE_NANGLES,time_1-time_0,1.0,0,GLOBAL_STROBETRACE_AVALANCHE_UPDATE_EVERY,GLOBAL_STROBETRACE_AVALANCHE_LOGDMIN,GLOBAL_STROBETRACE_AVALANCHE_LOGDMAX);}}}
    else if (t>=time_1 && t<time_2){ /* line motion illusion */
      if (GRATING_VS_LMI!=1){ 
	if (verbose){ printf(" %% entering second phase of suite three -- line motion illusion\n");}
	GRATING_VS_LMI = 1; STIMULUS_ONSET_TIME = time_1+time_delay;}
      if (CORTEX_BOTHER){
	if (PTREE_BOTHER!=0){
	  if (verbose){ printf(" %% %% dumping ptree\n");}
	  ptreerate(GLOBAL_PTREE); 
	  ptreedump_starter(GLOBAL_PTREE,"background",2,0,0,0,+1,-1); 
	  ptreetfree(GLOBAL_PTREE); GLOBAL_PTREE=NULL; PTREE_BOTHER=0;}
	if (STROBETRACE_BOTHER!=0){ 
	  if (verbose){ printf(" %% %% dumping strobetrace\n");}
	  strobetracedump(GLOBAL_STROBETRACE,0); 
	  strobetracedump(GLOBAL_STROBETRACE,1); 
	  strobetracetfree(GLOBAL_STROBETRACE); GLOBAL_STROBETRACE=NULL; STROBETRACE_BOTHER=0;}
	if (LMITRI_BOTHER==0){ 
	  if (verbose){ printf(" %% %% making lmitri\n");}
	  LMITRI_BOTHER=1;
	  GLOBAL_LMITRI = lmitrimake(STIMULUS_ONSET_TIME,GLOBAL_LMITRI_TIMELENGTH,GLOBAL_LMITRI_ROW_MAX,GLOBAL_LMITRI_ROW_MIN);}}}
    else if (t>=time_2 && t<time_3){ /* another line motion illusion */
      if (STIMULUS_ONSET_TIME!=time_2+time_delay){
	if (verbose){ printf(" %% entering third phase of suite three -- line motion illusion\n");}
	STIMULUS_ONSET_TIME = time_2+time_delay;}
      if (CORTEX_BOTHER){ 
	if (GLOBAL_LMITRI->time_start!=STIMULUS_ONSET_TIME){
	  GLOBAL_LMITRI->time_start = STIMULUS_ONSET_TIME;}}}
    else if (t>=time_3 && t<time_4){ /* yet another line motion illusion */
      if (STIMULUS_ONSET_TIME!=time_3+time_delay){
	if (verbose){ printf(" %% entering fourth phase of suite three -- line motion illusion\n");}
	STIMULUS_ONSET_TIME = time_3+time_delay;}
      if (CORTEX_BOTHER){ 
	if (GLOBAL_LMITRI->time_start!=STIMULUS_ONSET_TIME){
	  GLOBAL_LMITRI->time_start = STIMULUS_ONSET_TIME;}}}
    else if (t>=time_4 && t<time_5){ /* and another line motion illusion */
      if (STIMULUS_ONSET_TIME!=time_4+time_delay){
	if (verbose){ printf(" %% entering fifth phase of suite three -- line motion illusion\n");}
	STIMULUS_ONSET_TIME = time_4+time_delay;}
      if (CORTEX_BOTHER){ 
	if (GLOBAL_LMITRI->time_start!=STIMULUS_ONSET_TIME){
	  GLOBAL_LMITRI->time_start = STIMULUS_ONSET_TIME;}}}
    else if (t>=time_5 && t<time_6){ /* and a final line motion illusion */
      if (STIMULUS_ONSET_TIME!=time_5+time_delay){
	if (verbose){ printf(" %% entering sixth phase of suite three -- line motion illusion\n");}
	STIMULUS_ONSET_TIME = time_5+time_delay;}
      if (CORTEX_BOTHER){ 
	if (GLOBAL_LMITRI->time_start!=STIMULUS_ONSET_TIME){
	  GLOBAL_LMITRI->time_start = STIMULUS_ONSET_TIME;}}}
    else if (t>=time_6){ /* finish */
      STIMULUS_ONSET_TIME = time_6;
      if (GRATING_VS_LMI!=-10);{
	if (verbose){ printf(" %% entering seventh phase of suite three -- finishing\n");}
	GRATING_VS_LMI = -10; GLOBAL_TF = time_6+1;}
      if (CORTEX_BOTHER){
	if (LMITRI_BOTHER){
	  if (verbose){ printf(" %% %% dumping lmitri\n");}
	  lmitridump(GLOBAL_LMITRI,1); lmitritfree(GLOBAL_LMITRI); GLOBAL_LMITRI=NULL; LMITRI_BOTHER=0;}}}
    break;
  case 4:
    time_0=0;
    time_1=time_0+1024*SUITE_NSECONDS;
    time_2=time_1+1024*SUITE_NSECONDS;
    time_dumpevery=SUITE_DUMPEVERY;
    if (verbose>2){ printf(" times are: %0.2f, %0.2f, %0.2f \n",time_0,time_1,time_2);}
    if (t==0){ /* clear everything */
      if (verbose){ printf(" %% \n");}
      if (verbose){ printf(" %% initializing suite four\n");}
      GLOBAL_TF=0;
      if (GRATING_VS_LMI!=-10){ GRATING_VS_LMI = -10;}
      if (CORTEX_BOTHER){
	if (RTC_BOTHER!=0){ rtctfree(GLOBAL_RTC); GLOBAL_RTC=NULL; RTC_BOTHER=0;}
	if (STROBETRACE_BOTHER!=0){ strobetracetfree(GLOBAL_STROBETRACE); GLOBAL_STROBETRACE=NULL; STROBETRACE_BOTHER=0;}
	if (TUNINGCURVE_BOTHER!=0){ tuningcurvetfree(GLOBAL_TUNINGCURVE); GLOBAL_TUNINGCURVE=NULL; TUNINGCURVE_BOTHER=0;}
	if (LMITRI_BOTHER!=0){ lmitritfree(GLOBAL_LMITRI); GLOBAL_LMITRI=NULL; LMITRI_BOTHER=0;}
	if (PTREE_BOTHER!=0){ ptreetfree(GLOBAL_PTREE); GLOBAL_PTREE=NULL; GLOBAL_PTREE_BITBYBIT=0; PTREE_BOTHER=0;}
	if (CLOSET_BOTHER!=0){ closettfree(GLOBAL_CLOSET); GLOBAL_CLOSET=NULL; CLOSET_BOTHER=0;}
	if (POWER_BOTHER!=0){ powertfree(GLOBAL_POWER); GLOBAL_POWER=NULL; POWER_BOTHER=0;}
	if (YGGDRASIL_BOTHER!=0){ yggdrasiltfree(GLOBAL_YGGDRASIL); GLOBAL_YGGDRASIL=NULL; YGGDRASIL_BOTHER=0;}
	if (BONSAI_BOTHER!=0){ bonsaitfree(GLOBAL_BONSAI); GLOBAL_BONSAI=NULL; BONSAI_BOTHER=0;}
	if (HYDRA_BOTHER!=0){ hydratfree(GLOBAL_HYDRA); GLOBAL_HYDRA=NULL; HYDRA_BOTHER=0;}
	if (LYAPUNOV_BOTHER!=0){ lyapunovtfree(GLOBAL_LYAPUNOV); GLOBAL_LYAPUNOV=NULL; LYAPUNOV_BOTHER=0;}}}
    else if (t>=0 && t<time_1){ /* tuningcurve1 */
      STIMULUS_ONSET_TIME = time_0; OUTPUT_DUMP_EVERY = time_dumpevery; /* start dumping data */
      if (GRATING_VS_LMI!=0){ 
	if (verbose){ printf(" %% entering first phase of suite one -- tuningcurve1\n");}
	GRATING_VS_LMI = 0; INPUT_SPACEANGLE = 0;}
      if (CORTEX_BOTHER){
	if (GLOBAL_PTREE_BITBYBIT==0 && PTREE_BOTHER==0){ 
	  if (verbose){ printf(" %% %% making ptree\n");}
	  PTREE_BOTHER=1;GLOBAL_PTREE_BITBYBIT=1; GLOBAL_PTREE_REGION_TYPE=4;
	  GLOBAL_PTREE = ptreemake(GLOBAL_PTREE_NREGIONS,GLOBAL_PTREE_EVENT_WITHIN,GLOBAL_PTREE_EVENT_THRESHOLD,GLOBAL_PTREE_REGION_TYPE,GLOBAL_PTREE_NLEGS,GLOBAL_PTREE_LEGTIME);}
	if (OUTPUT_DUMP_EVERY>0 && (int)floor((GLOBAL_time-time_0)/OUTPUT_DUMP_EVERY)<(int)floor((GLOBAL_time-time_0+GLOBAL_DT)/OUTPUT_DUMP_EVERY)){
	  if (PTREE_BOTHER){
	    if (verbose){ printf(" %% %% dumping ptree\n");}
	    ptreerate(GLOBAL_PTREE);
	    ptreedump_starter(GLOBAL_PTREE,"tuningcurve1",2,0,0,0,+1,-1);
	    if (GLOBAL_PTREE_BITBYBIT){
	      if (verbose){ printf(" %% %% resetting ptree\n");}
	      ptreereset(GLOBAL_PTREE);}}}}}
    else if (t>=time_1 && t<time_2){ /* tuningcurve2 */
      STIMULUS_ONSET_TIME = time_1; OUTPUT_DUMP_EVERY = time_dumpevery; /* start dumping data */
      if (INPUT_SPACEANGLE==0){
	if (verbose){ printf(" %% entering second phase of suite four -- tuningcurve2\n");}
	if (CORTEX_BOTHER){
	  ptreetfree(GLOBAL_PTREE); GLOBAL_PTREE=NULL; GLOBAL_PTREE_BITBYBIT=0; PTREE_BOTHER=0;}
	GRATING_VS_LMI=0; INPUT_SPACEANGLE = 0 + INPUT_SPACEANGLE_BACON; /* six degrees of separation */}
      if (CORTEX_BOTHER){
	if (GLOBAL_PTREE_BITBYBIT==0 && PTREE_BOTHER==0){ 
	  if (verbose){ printf(" %% %% making ptree\n");}
	  PTREE_BOTHER=1;GLOBAL_PTREE_BITBYBIT=1; GLOBAL_PTREE_REGION_TYPE=4;
	  GLOBAL_PTREE = ptreemake(GLOBAL_PTREE_NREGIONS,GLOBAL_PTREE_EVENT_WITHIN,GLOBAL_PTREE_EVENT_THRESHOLD,GLOBAL_PTREE_REGION_TYPE,GLOBAL_PTREE_NLEGS,GLOBAL_PTREE_LEGTIME);}
	if (OUTPUT_DUMP_EVERY>0 && (int)floor((GLOBAL_time-time_1)/OUTPUT_DUMP_EVERY)<(int)floor((GLOBAL_time-time_1+GLOBAL_DT)/OUTPUT_DUMP_EVERY)){
	  if (PTREE_BOTHER){
	    if (verbose){ printf(" %% %% dumping ptree\n");}
	    ptreerate(GLOBAL_PTREE);
	    ptreedump_starter(GLOBAL_PTREE,"tuningcurve2",2,0,0,0,+1,-1);
	    if (GLOBAL_PTREE_BITBYBIT){
	      if (verbose){ printf(" %% %% resetting ptree\n");}
	      ptreereset(GLOBAL_PTREE);}}}}}
    else if (t>=time_2){ /* finish */
      STIMULUS_ONSET_TIME = time_2;
      if (GRATING_VS_LMI!=-10){ 
	if (verbose){ printf(" %% entering third phase of suite four -- finishing\n");}
	if (CORTEX_BOTHER){
	  ptreetfree(GLOBAL_PTREE); GLOBAL_PTREE=NULL; GLOBAL_PTREE_BITBYBIT=0; PTREE_BOTHER=0;}
	GRATING_VS_LMI = -10; GLOBAL_TF = time_2+1;}}
    break;
  case 5:
    time_0=0;
    time_1=time_0+1024*SUITE_NSECONDS;
    time_2=time_1+1024*SUITE_NSECONDS;
    time_3=time_2+1024*SUITE_NSECONDS;
    time_4=time_3+1024*SUITE_NSECONDS;
    time_dumpevery=SUITE_DUMPEVERY;
    if (verbose>2){ printf(" times are: %0.2f, %0.2f, %0.2f, %0.2f, %0.2f \n",time_0,time_1,time_2,time_3,time_4);}
    if (t==0){ /* clear everything */
      if (verbose){ printf(" %% \n");}
      if (verbose){ printf(" %% initializing suite five\n");}
      GLOBAL_TF=0;
      STIMULUS_ONSET_TIME=-1;
      if (GRATING_VS_LMI!=-10){ GRATING_VS_LMI = -10;}
      if (CORTEX_BOTHER){
	if (RTC_BOTHER!=0){ rtctfree(GLOBAL_RTC); GLOBAL_RTC=NULL; RTC_BOTHER=0;}
	if (STROBETRACE_BOTHER!=0){ strobetracetfree(GLOBAL_STROBETRACE); GLOBAL_STROBETRACE=NULL; STROBETRACE_BOTHER=0;}
	if (TUNINGCURVE_BOTHER!=0){ tuningcurvetfree(GLOBAL_TUNINGCURVE); GLOBAL_TUNINGCURVE=NULL; TUNINGCURVE_BOTHER=0;}
	if (LMITRI_BOTHER!=0){ lmitritfree(GLOBAL_LMITRI); GLOBAL_LMITRI=NULL; LMITRI_BOTHER=0;}
	if (PTREE_BOTHER!=0){ ptreetfree(GLOBAL_PTREE); GLOBAL_PTREE=NULL; GLOBAL_PTREE_BITBYBIT=0; PTREE_BOTHER=0;}
	if (CLOSET_BOTHER!=0){ closettfree(GLOBAL_CLOSET); GLOBAL_CLOSET=NULL; CLOSET_BOTHER=0;}
	if (POWER_BOTHER!=0){ powertfree(GLOBAL_POWER); GLOBAL_POWER=NULL; POWER_BOTHER=0;}
	if (YGGDRASIL_BOTHER!=0){ yggdrasiltfree(GLOBAL_YGGDRASIL); GLOBAL_YGGDRASIL=NULL; YGGDRASIL_BOTHER=0;}
	if (BONSAI_BOTHER!=0){ bonsaitfree(GLOBAL_BONSAI); GLOBAL_BONSAI=NULL; BONSAI_BOTHER=0;}
	if (HYDRA_BOTHER!=0){ hydratfree(GLOBAL_HYDRA); GLOBAL_HYDRA=NULL; HYDRA_BOTHER=0;}
	if (LYAPUNOV_BOTHER!=0){ lyapunovtfree(GLOBAL_LYAPUNOV); GLOBAL_LYAPUNOV=NULL; LYAPUNOV_BOTHER=0;}}}
    else if (t>=0 && t<time_1){ /* rate1 and yggdrasil*/
      if (STIMULUS_ONSET_TIME != time_0){
	if (verbose){ printf(" %% entering first phase of suite five -- rate1\n");}
	GRATING_VS_LMI = -10;
	STIMULUS_ONSET_TIME = time_0; OUTPUT_DUMP_EVERY = time_dumpevery; /* start dumping data */
	if (CORTEX_BOTHER){
	  PTREE_BOTHER=1; GLOBAL_PTREE_BITBYBIT=1; GLOBAL_PTREE_REGION_TYPE=1;
	  GLOBAL_PTREE = ptreemake(GLOBAL_PTREE_NREGIONS,GLOBAL_PTREE_EVENT_WITHIN,GLOBAL_PTREE_EVENT_THRESHOLD,GLOBAL_PTREE_REGION_TYPE,GLOBAL_PTREE_NLEGS,GLOBAL_PTREE_LEGTIME);
	  POWER_BOTHER=1; GLOBAL_POWER = powermake(GLOBAL_Nra,GLOBAL_POWER_LENGTH,GLOBAL_POWER_HOWMANY);
	  YGGDRASIL_BOTHER=1; GLOBAL_YGGDRASIL = yggdrasilmake(GLOBAL_PTREE_NREGIONS,GLOBAL_PTREE_EVENT_WITHIN,GLOBAL_PTREE_EVENT_THRESHOLD,GLOBAL_PTREE_REGION_TYPE,GLOBAL_PTREE_NLEGS,GLOBAL_PTREE_LEGTIME,GLOBAL_YGGDRASIL_PPNREGIONS,GLOBAL_YGGDRASIL_PPNLEGS,GLOBAL_YGGDRASIL_PPLEGTIME,GLOBAL_YGGDRASIL_WEIGHT_MINIMUM);}}
      if (CORTEX_BOTHER){
	if (OUTPUT_DUMP_EVERY>0 && (int)floor((GLOBAL_time-time_0)/OUTPUT_DUMP_EVERY)<(int)floor((GLOBAL_time-time_0+GLOBAL_DT)/OUTPUT_DUMP_EVERY)){
	  if (PTREE_BOTHER){
	    if (verbose){ printf(" %% %% dumping ptree\n");}
	    ptreerate(GLOBAL_PTREE);
	    ptreedump_starter(GLOBAL_PTREE,"rate1",2,0,0,0,+1,-1);
	    if (GLOBAL_PTREE_BITBYBIT){
	      if (verbose){ printf(" %% %% resetting ptree\n");}
	      ptreereset(GLOBAL_PTREE);}}}}}
    else if (t>=time_1 && t<time_2){ /* rate2 */
      if (STIMULUS_ONSET_TIME != time_1){
	if (verbose){ printf(" %% entering second phase of suite five -- rate2\n");}
	GRATING_VS_LMI = -10; LGN_BACKRATE *= 1 + INPUT_SPACEANGLE_BACON;
	STIMULUS_ONSET_TIME = time_1; OUTPUT_DUMP_EVERY = time_dumpevery; /* start dumping data */
	if (CORTEX_BOTHER){
	  if (GLOBAL_PTREE!=NULL){ ptreetfree(GLOBAL_PTREE); GLOBAL_PTREE=NULL;} GLOBAL_PTREE_BITBYBIT=0; PTREE_BOTHER=0;
	  PTREE_BOTHER=1; GLOBAL_PTREE_BITBYBIT=1; GLOBAL_PTREE_REGION_TYPE=1;
	  GLOBAL_PTREE = ptreemake(GLOBAL_PTREE_NREGIONS,GLOBAL_PTREE_EVENT_WITHIN,GLOBAL_PTREE_EVENT_THRESHOLD,GLOBAL_PTREE_REGION_TYPE,GLOBAL_PTREE_NLEGS,GLOBAL_PTREE_LEGTIME);
	  if (POWER_BOTHER){ 
	    if (verbose){ printf(" %% %% dumping power\n");} 
	    powerdump(GLOBAL_POWER,"rate1",0); powerdump(GLOBAL_POWER,"rate1",1);}
	  if (GLOBAL_POWER!=NULL){ powertfree(GLOBAL_POWER); GLOBAL_POWER=NULL;}
	  POWER_BOTHER=1; GLOBAL_POWER = powermake(GLOBAL_Nra,GLOBAL_POWER_LENGTH,GLOBAL_POWER_HOWMANY);
	  if (YGGDRASIL_BOTHER){
	    if (verbose){ printf(" %% %% dumping yggdrasil\n");}
	    ptreerate(GLOBAL_YGGDRASIL->p);
	    ptreedump_starter(GLOBAL_YGGDRASIL->p,"yggdrasil_p",2,-1,0,0,+1,-1);
	    if (GLOBAL_YGGDRASIL->pp!=NULL){ 
	      ptreerate(GLOBAL_YGGDRASIL->pp);
	      ptreedump_starter(GLOBAL_YGGDRASIL->pp,"yggdrasil_pp",2,-1,0,0,+1,-1);}
	    yggdrasiltfree(GLOBAL_YGGDRASIL); GLOBAL_YGGDRASIL=NULL; YGGDRASIL_BOTHER=0;}}}
      if (CORTEX_BOTHER){
	if (OUTPUT_DUMP_EVERY>0 && (int)floor((GLOBAL_time-time_1)/OUTPUT_DUMP_EVERY)<(int)floor((GLOBAL_time-time_1+GLOBAL_DT)/OUTPUT_DUMP_EVERY)){
	  if (PTREE_BOTHER){
	    if (verbose){ printf(" %% %% dumping ptree\n");}
	    ptreerate(GLOBAL_PTREE);
	    ptreedump_starter(GLOBAL_PTREE,"rate2",2,0,0,0,+1,-1);
	    if (GLOBAL_PTREE_BITBYBIT){
	      if (verbose){ printf(" %% %% resetting ptree\n");}
	      ptreereset(GLOBAL_PTREE);}}}}}
    else if (t>=time_2 && t<time_3){ /* lyapunov */
      if (STIMULUS_ONSET_TIME != time_2){
	if (verbose){ printf(" %% entering third phase of suite five -- lyapunov\n");}
	GRATING_VS_LMI = -10; LGN_BACKRATE /= 1 + INPUT_SPACEANGLE_BACON;
	STIMULUS_ONSET_TIME = time_2; OUTPUT_DUMP_EVERY = 0; /* stop dumping data */
	if (CORTEX_BOTHER){
	  if (GLOBAL_PTREE!=NULL){ ptreetfree(GLOBAL_PTREE); GLOBAL_PTREE=NULL;} GLOBAL_PTREE_BITBYBIT=0; PTREE_BOTHER=0;
	  if (POWER_BOTHER){ 
	    if (verbose){ printf(" %% %% dumping power\n");} 
	    powerdump(GLOBAL_POWER,"rate2",0); powerdump(GLOBAL_POWER,"rate2",1);}
	  if (GLOBAL_POWER!=NULL){ powertfree(GLOBAL_POWER); GLOBAL_POWER=NULL;} POWER_BOTHER=0;
	  if (GLOBAL_LYAPUNOV!=NULL){ lyapunovtfree(GLOBAL_LYAPUNOV); GLOBAL_LYAPUNOV=NULL;} LYAPUNOV_BOTHER=0;
	  LYAPUNOV_BOTHER=1;
	  GLOBAL_LYAPUNOV = lyapunovmake(GLOBAL_Nra,GLOBAL_LYAPUNOV_UPDATE_EVERY,GLOBAL_LYAPUNOV_JIGGLE);}}}
    else if (t>=time_3 && t<time_4){ /* bonsai */
      if (STIMULUS_ONSET_TIME != time_3){
	if (verbose){ printf(" %% entering fourth phase of suite five -- bonsai\n");}
	GRATING_VS_LMI=-10;
	STIMULUS_ONSET_TIME=time_3;
	if (CORTEX_BOTHER){
	  if (LYAPUNOV_BOTHER){ lyapunovdump(GLOBAL_LYAPUNOV,0); lyapunovdump(GLOBAL_LYAPUNOV,1);}
	  if (GLOBAL_LYAPUNOV!=NULL){ lyapunovtfree(GLOBAL_LYAPUNOV); GLOBAL_LYAPUNOV=NULL;} LYAPUNOV_BOTHER=0;
	  BONSAI_BOTHER=1; GLOBAL_BONSAI = bonsaimake(GLOBAL_Nra,GLOBAL_PTREE_NREGIONS,GLOBAL_PTREE_EVENT_WITHIN,GLOBAL_PTREE_EVENT_THRESHOLD,GLOBAL_PTREE_REGION_TYPE,GLOBAL_PTREE_NLEGS,GLOBAL_PTREE_LEGTIME);}}}
    else if (t>=time_4){ /* finish */
      if (STIMULUS_ONSET_TIME != time_4){
	if (verbose){ printf(" %% entering fifth phase of suite five -- finishing\n");}
	GRATING_VS_LMI=-10;
	STIMULUS_ONSET_TIME=time_4;
	GLOBAL_TF = time_4+1;
	if (CORTEX_BOTHER){
	  connectionsdump(GLOBAL_Nra,1,"connections");
	  if (GLOBAL_BONSAI!=NULL){ 
	    bonsaidump(GLOBAL_BONSAI,"bonsai",0); bonsaidump(GLOBAL_BONSAI,"bonsai",1); bonsaidump(GLOBAL_BONSAI,"bonsai",2);
	    bonsaitfree(GLOBAL_BONSAI); GLOBAL_BONSAI=NULL; BONSAI_BOTHER=0;}}
	fnamebase = (char **) tcalloc(2,sizeof(char *));
	fnamebase[0] = (char *) tcalloc(256,sizeof(char)); sprintf(fnamebase[0],"ptree_rate1_%srecord",gs2);
	fnamebase[1] = (char *) tcalloc(256,sizeof(char)); sprintf(fnamebase[1],"ptree_rate2_%srecord",gs2);
	for (tindex=1;tindex<=SUITE_TINDEXMAX;tindex++){ for (sindex=0;sindex<=SUITE_SINDEXMAX;sindex++){ for (dindex=0;dindex<=SUITE_DINDEXMAX;dindex++){
	  sprintf(outputname,"r12"); ptree_trialaverage_helper(2,GLOBAL_PTREE_NLEGS+1,tindex,1,1,0,sindex,dindex,fnamebase,outputname);}}}
	tfree(fnamebase[0]);fnamebase[0]=NULL;
	tfree(fnamebase[1]);fnamebase[1]=NULL;
	tfree(fnamebase);fnamebase=NULL;
	if (SUITE_BITBYBIT_REMOVE){ 
	  for (nr=0;nr<(int)floor((double)SUITE_NSECONDS*1024/(double)SUITE_DUMPEVERY)+1;nr++){
	    sprintf(command,"""rm"" ./ptree_rate1_%srecord_%d ;",gs2,nr); system(command);
	    sprintf(command,"""rm"" ./ptree_rate2_%srecord_%d ;",gs2,nr); system(command);}}}}
    break;
  case 6:
    time_0=0;
    time_1=time_0+1024*SUITE_NSECONDS;
    time_2=time_1+1024*SUITE_NSECONDS;
    time_3=time_2+1024*SUITE_NSECONDS;
    time_dumpevery=SUITE_DUMPEVERY;
    if (verbose>2){ printf(" times are: %0.2f, %0.2f, %0.2f, %0.2f \n",time_0,time_1,time_2,time_3);}
    if (t==0){ /* clear everything */
      if (verbose){ printf(" %% \n");}
      if (verbose){ printf(" %% initializing suite six\n");}
      GLOBAL_TF=0;
      STIMULUS_ONSET_TIME=-1;
      if (GRATING_VS_LMI!=-10){ GRATING_VS_LMI = -10;}
      if (CORTEX_BOTHER){
	if (RTC_BOTHER!=0){ rtctfree(GLOBAL_RTC); GLOBAL_RTC=NULL; RTC_BOTHER=0;}
	if (STROBETRACE_BOTHER!=0){ strobetracetfree(GLOBAL_STROBETRACE); GLOBAL_STROBETRACE=NULL; STROBETRACE_BOTHER=0;}
	if (TUNINGCURVE_BOTHER!=0){ tuningcurvetfree(GLOBAL_TUNINGCURVE); GLOBAL_TUNINGCURVE=NULL; TUNINGCURVE_BOTHER=0;}
	if (LMITRI_BOTHER!=0){ lmitritfree(GLOBAL_LMITRI); GLOBAL_LMITRI=NULL; LMITRI_BOTHER=0;}
	if (PTREE_BOTHER!=0){ ptreetfree(GLOBAL_PTREE); GLOBAL_PTREE=NULL; GLOBAL_PTREE_BITBYBIT=0; PTREE_BOTHER=0;}
	if (CLOSET_BOTHER!=0){ closettfree(GLOBAL_CLOSET); GLOBAL_CLOSET=NULL; CLOSET_BOTHER=0;}
	if (POWER_BOTHER!=0){ powertfree(GLOBAL_POWER); GLOBAL_POWER=NULL; POWER_BOTHER=0;}
	if (YGGDRASIL_BOTHER!=0){ yggdrasiltfree(GLOBAL_YGGDRASIL); GLOBAL_YGGDRASIL=NULL; YGGDRASIL_BOTHER=0;}
	if (BONSAI_BOTHER!=0){ bonsaitfree(GLOBAL_BONSAI); GLOBAL_BONSAI=NULL; BONSAI_BOTHER=0;}
	if (HYDRA_BOTHER!=0){ hydratfree(GLOBAL_HYDRA); GLOBAL_HYDRA=NULL; HYDRA_BOTHER=0;}
	if (LYAPUNOV_BOTHER!=0){ lyapunovtfree(GLOBAL_LYAPUNOV); GLOBAL_LYAPUNOV=NULL; LYAPUNOV_BOTHER=0;}}}
    else if (t>=0 && t<time_1){ /* rate1 */
      if (STIMULUS_ONSET_TIME != time_0){
	if (verbose){ printf(" %% entering first phase of suite six -- rate1\n");}
	GRATING_VS_LMI = -10;
	STIMULUS_ONSET_TIME = time_0; OUTPUT_DUMP_EVERY = time_dumpevery; /* start dumping data */
	if (CORTEX_BOTHER){
	  PTREE_BOTHER=1; GLOBAL_PTREE_BITBYBIT=1; GLOBAL_PTREE_REGION_TYPE=1;
	  GLOBAL_PTREE = ptreemake(GLOBAL_PTREE_NREGIONS,GLOBAL_PTREE_EVENT_WITHIN,GLOBAL_PTREE_EVENT_THRESHOLD,GLOBAL_PTREE_REGION_TYPE,GLOBAL_PTREE_NLEGS,GLOBAL_PTREE_LEGTIME);
	  POWER_BOTHER=1; GLOBAL_POWER = powermake(GLOBAL_Nra,GLOBAL_POWER_LENGTH,GLOBAL_POWER_HOWMANY);}}
      if (CORTEX_BOTHER){
	if (OUTPUT_DUMP_EVERY>0 && (int)floor((GLOBAL_time-time_0)/OUTPUT_DUMP_EVERY)<(int)floor((GLOBAL_time-time_0+GLOBAL_DT)/OUTPUT_DUMP_EVERY)){
	  if (PTREE_BOTHER){
	    if (verbose){ printf(" %% %% dumping ptree\n");}
	    ptreerate(GLOBAL_PTREE);
	    ptreedump_starter(GLOBAL_PTREE,"rate1",2,0,0,0,+1,-1);
	    if (GLOBAL_PTREE_BITBYBIT){
	      if (verbose){ printf(" %% %% resetting ptree\n");}
	      ptreereset(GLOBAL_PTREE);}}}}}
    else if (t>=time_1 && t<time_2){ /* rate2 */
      if (STIMULUS_ONSET_TIME != time_1){
	if (verbose){ printf(" %% entering second phase of suite six -- rate2\n");}
	GRATING_VS_LMI = -10; LGN_BACKRATE *= 1 + INPUT_SPACEANGLE_BACON;
	STIMULUS_ONSET_TIME = time_1; OUTPUT_DUMP_EVERY = time_dumpevery; /* start dumping data */
	if (CORTEX_BOTHER){
	  if (GLOBAL_PTREE!=NULL){ ptreetfree(GLOBAL_PTREE); GLOBAL_PTREE=NULL;} GLOBAL_PTREE_BITBYBIT=0; PTREE_BOTHER=0;
	  PTREE_BOTHER=1; GLOBAL_PTREE_BITBYBIT=1; GLOBAL_PTREE_REGION_TYPE=1;
	  GLOBAL_PTREE = ptreemake(GLOBAL_PTREE_NREGIONS,GLOBAL_PTREE_EVENT_WITHIN,GLOBAL_PTREE_EVENT_THRESHOLD,GLOBAL_PTREE_REGION_TYPE,GLOBAL_PTREE_NLEGS,GLOBAL_PTREE_LEGTIME);
	  if (POWER_BOTHER){ 
	    if (verbose){ printf(" %% %% dumping power\n");} 
	    powerdump(GLOBAL_POWER,"rate1",0); powerdump(GLOBAL_POWER,"rate1",1);}
	  if (GLOBAL_POWER!=NULL){ powertfree(GLOBAL_POWER); GLOBAL_POWER=NULL;}
	  POWER_BOTHER=1; GLOBAL_POWER = powermake(GLOBAL_Nra,GLOBAL_POWER_LENGTH,GLOBAL_POWER_HOWMANY);}}
      if (CORTEX_BOTHER){
	if (OUTPUT_DUMP_EVERY>0 && (int)floor((GLOBAL_time-time_1)/OUTPUT_DUMP_EVERY)<(int)floor((GLOBAL_time-time_1+GLOBAL_DT)/OUTPUT_DUMP_EVERY)){
	  if (PTREE_BOTHER){
	    if (verbose){ printf(" %% %% dumping ptree\n");}
	    ptreerate(GLOBAL_PTREE);
	    ptreedump_starter(GLOBAL_PTREE,"rate2",2,0,0,0,+1,-1);
	    if (GLOBAL_PTREE_BITBYBIT){
	      if (verbose){ printf(" %% %% resetting ptree\n");}
	      ptreereset(GLOBAL_PTREE);}}}}}
    else if (t>=time_2 && t<time_3){ /* rate3 */
      if (STIMULUS_ONSET_TIME != time_2){
	if (verbose){ printf(" %% entering third phase of suite six -- rate3\n");}
	GRATING_VS_LMI = -10; LGN_BACKRATE /= 1 + INPUT_SPACEANGLE_BACON; LGN_STRENGTH *= 1 + INPUT_SPACEANGLE_BACON; 
	STIMULUS_ONSET_TIME = time_2; OUTPUT_DUMP_EVERY = time_dumpevery; /* start dumping data */
	if (CORTEX_BOTHER){
	  if (GLOBAL_PTREE!=NULL){ ptreetfree(GLOBAL_PTREE); GLOBAL_PTREE=NULL;} GLOBAL_PTREE_BITBYBIT=0; PTREE_BOTHER=0;
	  PTREE_BOTHER=1; GLOBAL_PTREE_BITBYBIT=1; GLOBAL_PTREE_REGION_TYPE=1;
	  GLOBAL_PTREE = ptreemake(GLOBAL_PTREE_NREGIONS,GLOBAL_PTREE_EVENT_WITHIN,GLOBAL_PTREE_EVENT_THRESHOLD,GLOBAL_PTREE_REGION_TYPE,GLOBAL_PTREE_NLEGS,GLOBAL_PTREE_LEGTIME);
	  if (POWER_BOTHER){ 
	    if (verbose){ printf(" %% %% dumping power\n");} 
	    powerdump(GLOBAL_POWER,"rate2",0); powerdump(GLOBAL_POWER,"rate2",1);}
	  if (GLOBAL_POWER!=NULL){ powertfree(GLOBAL_POWER); GLOBAL_POWER=NULL;}
	  POWER_BOTHER=1; GLOBAL_POWER = powermake(GLOBAL_Nra,GLOBAL_POWER_LENGTH,GLOBAL_POWER_HOWMANY);}}
      if (CORTEX_BOTHER){
	if (OUTPUT_DUMP_EVERY>0 && (int)floor((GLOBAL_time-time_2)/OUTPUT_DUMP_EVERY)<(int)floor((GLOBAL_time-time_2+GLOBAL_DT)/OUTPUT_DUMP_EVERY)){
	  if (PTREE_BOTHER){
	    if (verbose){ printf(" %% %% dumping ptree\n");}
	    ptreerate(GLOBAL_PTREE);
	    ptreedump_starter(GLOBAL_PTREE,"rate3",2,0,0,0,+1,-1);
	    if (GLOBAL_PTREE_BITBYBIT){
	      if (verbose){ printf(" %% %% resetting ptree\n");}
	      ptreereset(GLOBAL_PTREE);}}}}}
    else if (t>=time_3){ /* finish */
      if (STIMULUS_ONSET_TIME != time_3){
	if (verbose){ printf(" %% entering fourth phase of suite six -- finishing\n");}
	GRATING_VS_LMI=-10;
	STIMULUS_ONSET_TIME=time_3;
	GLOBAL_TF = time_3+1;
	if (CORTEX_BOTHER){
	  connectionsdump(GLOBAL_Nra,1,"connections");
	  if (GLOBAL_PTREE!=NULL){ ptreetfree(GLOBAL_PTREE); GLOBAL_PTREE=NULL;} GLOBAL_PTREE_BITBYBIT=0; PTREE_BOTHER=0;
	  if (POWER_BOTHER){ 
	    if (verbose){ printf(" %% %% dumping power\n");} 
	    powerdump(GLOBAL_POWER,"rate3",0); powerdump(GLOBAL_POWER,"rate3",1);}
	  if (GLOBAL_POWER!=NULL){ powertfree(GLOBAL_POWER); GLOBAL_POWER=NULL;} POWER_BOTHER=0;}
	fnamebase = (char **) tcalloc(3,sizeof(char *));
	fnamebase[0] = (char *) tcalloc(256,sizeof(char)); sprintf(fnamebase[0],"ptree_rate1_%srecord",gs2);
	fnamebase[1] = (char *) tcalloc(256,sizeof(char)); sprintf(fnamebase[1],"ptree_rate2_%srecord",gs2);
	fnamebase[2] = (char *) tcalloc(256,sizeof(char)); sprintf(fnamebase[2],"ptree_rate3_%srecord",gs2);
	for (tindex=1;tindex<=SUITE_TINDEXMAX;tindex++){ for (sindex=0;sindex<=SUITE_SINDEXMAX;sindex++){ for (dindex=0;dindex<=SUITE_DINDEXMAX;dindex++){
	  sprintf(outputname,"r12");ptree_trialaverage_helper(2,GLOBAL_PTREE_NLEGS+1,tindex,1,1,-1,sindex,dindex,&fnamebase[0],outputname);
	  sprintf(outputname,"r23");ptree_trialaverage_helper(2,GLOBAL_PTREE_NLEGS+1,tindex,1,1,-1,sindex,dindex,&fnamebase[1],outputname);
	  sprintf(outputname,"r123");ptree_trialaverage_helper(3,GLOBAL_PTREE_NLEGS+1,tindex,1,1,-1,sindex,dindex,&fnamebase[0],outputname);}}}
	sprintf(fnamebase[1],"ptree_rate3_%srecord",gs2);
	for (tindex=1;tindex<=SUITE_TINDEXMAX;tindex++){ for (sindex=0;sindex<=SUITE_SINDEXMAX;sindex++){ for (dindex=0;dindex<=SUITE_DINDEXMAX;dindex++){
	  sprintf(outputname,"r13");ptree_trialaverage_helper(2,GLOBAL_PTREE_NLEGS+1,tindex,1,1,-1,sindex,dindex,&fnamebase[0],outputname);}}}
	tfree(fnamebase[0]);fnamebase[0]=NULL;
	tfree(fnamebase[1]);fnamebase[1]=NULL;
	tfree(fnamebase[2]);fnamebase[2]=NULL;
	tfree(fnamebase);fnamebase=NULL;
	if (SUITE_BITBYBIT_REMOVE){ 
	  for (nr=1;nr<(int)floor((double)SUITE_NSECONDS*1024/(double)SUITE_DUMPEVERY)+1;nr++){
	    sprintf(command,"""rm"" ./ptree_rate1_%srecord_%d ;",gs2,nr); system(command);
	    sprintf(command,"""rm"" ./ptree_rate2_%srecord_%d ;",gs2,nr); system(command);
	    sprintf(command,"""rm"" ./ptree_rate3_%srecord_%d ;",gs2,nr); system(command);}}}}
    break;
  case 7:
    time_0=0;
    time_1=time_0+1024*SUITE_NSECONDS;
    time_2=time_1+1024*SUITE_NSECONDS;
    time_3=time_2+1024*SUITE_NSECONDS;
    time_4=time_3+1024*SUITE_NSECONDS;
    time_5=time_4+1024*SUITE_NSECONDS;
    time_6=time_5+1024*SUITE_NSECONDS;
    time_7=time_6+1024*SUITE_NSECONDS;
    time_dumpevery=SUITE_DUMPEVERY;
    if (verbose>2){ printf(" times are: %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f \n",time_0,time_1,time_2,time_3,time_4,time_5,time_6,time_7);}
    if (t==0){ /* clear everything */
      if (verbose){ printf(" %% \n");}
      if (verbose){ printf(" %% initializing suite seven\n");}
      GLOBAL_TF=0;
      STIMULUS_ONSET_TIME=-1;
      if (GRATING_VS_LMI!=-10){ GRATING_VS_LMI = -10;}
      if (CORTEX_BOTHER){
	if (RTC_BOTHER!=0){ rtctfree(GLOBAL_RTC); GLOBAL_RTC=NULL; RTC_BOTHER=0;}
	if (STROBETRACE_BOTHER!=0){ strobetracetfree(GLOBAL_STROBETRACE); GLOBAL_STROBETRACE=NULL; STROBETRACE_BOTHER=0;}
	if (TUNINGCURVE_BOTHER!=0){ tuningcurvetfree(GLOBAL_TUNINGCURVE); GLOBAL_TUNINGCURVE=NULL; TUNINGCURVE_BOTHER=0;}
	if (LMITRI_BOTHER!=0){ lmitritfree(GLOBAL_LMITRI); GLOBAL_LMITRI=NULL; LMITRI_BOTHER=0;}
	if (PTREE_BOTHER!=0){ ptreetfree(GLOBAL_PTREE); GLOBAL_PTREE=NULL; GLOBAL_PTREE_BITBYBIT=0; PTREE_BOTHER=0;}
	if (CLOSET_BOTHER!=0){ closettfree(GLOBAL_CLOSET); GLOBAL_CLOSET=NULL; CLOSET_BOTHER=0;}
	if (POWER_BOTHER!=0){ powertfree(GLOBAL_POWER); GLOBAL_POWER=NULL; POWER_BOTHER=0;}
	if (YGGDRASIL_BOTHER!=0){ yggdrasiltfree(GLOBAL_YGGDRASIL); GLOBAL_YGGDRASIL=NULL; YGGDRASIL_BOTHER=0;}
	if (BONSAI_BOTHER!=0){ bonsaitfree(GLOBAL_BONSAI); GLOBAL_BONSAI=NULL; BONSAI_BOTHER=0;}
	if (HYDRA_BOTHER!=0){ hydratfree(GLOBAL_HYDRA); GLOBAL_HYDRA=NULL; HYDRA_BOTHER=0;}
	if (LYAPUNOV_BOTHER!=0){ lyapunovtfree(GLOBAL_LYAPUNOV); GLOBAL_LYAPUNOV=NULL; LYAPUNOV_BOTHER=0;}}}
    else if (t>=0 && t<time_1){ /* rate1 */
      if (STIMULUS_ONSET_TIME != time_0){
	if (verbose){ printf(" %% entering first phase of suite seven -- rate1\n");}
	GRATING_VS_LMI = -10; LGN_DUMB=0;
	STIMULUS_ONSET_TIME = time_0; OUTPUT_DUMP_EVERY = time_dumpevery; /* start dumping data */
	if (CORTEX_BOTHER){
	  PTREE_BOTHER=1; GLOBAL_PTREE_BITBYBIT=1; GLOBAL_PTREE_REGION_TYPE=1;
	  GLOBAL_PTREE = ptreemake(GLOBAL_PTREE_NREGIONS,GLOBAL_PTREE_EVENT_WITHIN,GLOBAL_PTREE_EVENT_THRESHOLD,GLOBAL_PTREE_REGION_TYPE,GLOBAL_PTREE_NLEGS,GLOBAL_PTREE_LEGTIME);
	  POWER_BOTHER=1; GLOBAL_POWER = powermake(GLOBAL_Nra,GLOBAL_POWER_LENGTH,GLOBAL_POWER_HOWMANY);}}
      if (CORTEX_BOTHER){
	if (OUTPUT_DUMP_EVERY>0 && (int)floor((GLOBAL_time-time_0)/OUTPUT_DUMP_EVERY)<(int)floor((GLOBAL_time-time_0+GLOBAL_DT)/OUTPUT_DUMP_EVERY)){
	  if (PTREE_BOTHER){
	    if (verbose){ printf(" %% %% dumping ptree\n");}
	    ptreerate(GLOBAL_PTREE);
	    ptreedump_starter(GLOBAL_PTREE,"rate1",2,0,0,0,+1,-1);
	    if (GLOBAL_PTREE_BITBYBIT){
	      if (verbose){ printf(" %% %% resetting ptree\n");}
	      ptreereset(GLOBAL_PTREE);}}}}}
    else if (t>=time_1 && t<time_2){ /* rate2 */
      if (STIMULUS_ONSET_TIME != time_1){
	if (verbose){ printf(" %% entering second phase of suite seven -- rate2\n");}
	GRATING_VS_LMI = -10; LGN_BACKRATE *= 1 + INPUT_SPACEANGLE_BACON; LGN_DUMB=0;
	STIMULUS_ONSET_TIME = time_1; OUTPUT_DUMP_EVERY = time_dumpevery; /* start dumping data */
	if (CORTEX_BOTHER){
	  if (GLOBAL_PTREE!=NULL){ ptreetfree(GLOBAL_PTREE); GLOBAL_PTREE=NULL;} GLOBAL_PTREE_BITBYBIT=0; PTREE_BOTHER=0;
	  PTREE_BOTHER=1; GLOBAL_PTREE_BITBYBIT=1; GLOBAL_PTREE_REGION_TYPE=1;
	  GLOBAL_PTREE = ptreemake(GLOBAL_PTREE_NREGIONS,GLOBAL_PTREE_EVENT_WITHIN,GLOBAL_PTREE_EVENT_THRESHOLD,GLOBAL_PTREE_REGION_TYPE,GLOBAL_PTREE_NLEGS,GLOBAL_PTREE_LEGTIME);
	  if (POWER_BOTHER){ 
	    if (verbose){ printf(" %% %% dumping power\n");} 
	    powerdump(GLOBAL_POWER,"rate1",0); powerdump(GLOBAL_POWER,"rate1",1);}
	  if (GLOBAL_POWER!=NULL){ powertfree(GLOBAL_POWER); GLOBAL_POWER=NULL;}
	  POWER_BOTHER=1; GLOBAL_POWER = powermake(GLOBAL_Nra,GLOBAL_POWER_LENGTH,GLOBAL_POWER_HOWMANY);}}
      if (CORTEX_BOTHER){
	if (OUTPUT_DUMP_EVERY>0 && (int)floor((GLOBAL_time-time_1)/OUTPUT_DUMP_EVERY)<(int)floor((GLOBAL_time-time_1+GLOBAL_DT)/OUTPUT_DUMP_EVERY)){
	  if (PTREE_BOTHER){
	    if (verbose){ printf(" %% %% dumping ptree\n");}
	    ptreerate(GLOBAL_PTREE);
	    ptreedump_starter(GLOBAL_PTREE,"rate2",2,0,0,0,+1,-1);
	    if (GLOBAL_PTREE_BITBYBIT){
	      if (verbose){ printf(" %% %% resetting ptree\n");}
	      ptreereset(GLOBAL_PTREE);}}}}}
    else if (t>=time_2 && t<time_3){ /* rate3 */
      if (STIMULUS_ONSET_TIME != time_2){
	if (verbose){ printf(" %% entering third phase of suite seven -- rate3\n");}
	GRATING_VS_LMI = -10; LGN_BACKRATE /= 1 + INPUT_SPACEANGLE_BACON; LGN_STRENGTH *= 1 + INPUT_SPACEANGLE_BACON; LGN_DUMB=0;
	STIMULUS_ONSET_TIME = time_2; OUTPUT_DUMP_EVERY = time_dumpevery; /* start dumping data */
	if (CORTEX_BOTHER){
	  if (GLOBAL_PTREE!=NULL){ ptreetfree(GLOBAL_PTREE); GLOBAL_PTREE=NULL;} GLOBAL_PTREE_BITBYBIT=0; PTREE_BOTHER=0;
	  PTREE_BOTHER=1; GLOBAL_PTREE_BITBYBIT=1; GLOBAL_PTREE_REGION_TYPE=1;
	  GLOBAL_PTREE = ptreemake(GLOBAL_PTREE_NREGIONS,GLOBAL_PTREE_EVENT_WITHIN,GLOBAL_PTREE_EVENT_THRESHOLD,GLOBAL_PTREE_REGION_TYPE,GLOBAL_PTREE_NLEGS,GLOBAL_PTREE_LEGTIME);
	  if (POWER_BOTHER){ 
	    if (verbose){ printf(" %% %% dumping power\n");} 
	    powerdump(GLOBAL_POWER,"rate2",0); powerdump(GLOBAL_POWER,"rate2",1);}
	  if (GLOBAL_POWER!=NULL){ powertfree(GLOBAL_POWER); GLOBAL_POWER=NULL;}
	  POWER_BOTHER=1; GLOBAL_POWER = powermake(GLOBAL_Nra,GLOBAL_POWER_LENGTH,GLOBAL_POWER_HOWMANY);}}
      if (CORTEX_BOTHER){
	if (OUTPUT_DUMP_EVERY>0 && (int)floor((GLOBAL_time-time_2)/OUTPUT_DUMP_EVERY)<(int)floor((GLOBAL_time-time_2+GLOBAL_DT)/OUTPUT_DUMP_EVERY)){
	  if (PTREE_BOTHER){
	    if (verbose){ printf(" %% %% dumping ptree\n");}
	    ptreerate(GLOBAL_PTREE);
	    ptreedump_starter(GLOBAL_PTREE,"rate3",2,0,0,0,+1,-1);
	    if (GLOBAL_PTREE_BITBYBIT){
	      if (verbose){ printf(" %% %% resetting ptree\n");}
	      ptreereset(GLOBAL_PTREE);}}}}}
    else if (t>=time_3 && t<time_4){ /* rate4 */
      if (STIMULUS_ONSET_TIME != time_3){
	if (verbose){ printf(" %% entering fourth phase of suite seven -- rate4\n");}
	GRATING_VS_LMI = -10; LGN_STRENGTH /= 1 + INPUT_SPACEANGLE_BACON; LGN_DUMB=1;
	STIMULUS_ONSET_TIME = time_3; OUTPUT_DUMP_EVERY = time_dumpevery; /* start dumping data */
	if (CORTEX_BOTHER){
	  if (GLOBAL_PTREE!=NULL){ ptreetfree(GLOBAL_PTREE); GLOBAL_PTREE=NULL;} GLOBAL_PTREE_BITBYBIT=0; PTREE_BOTHER=0;
	  PTREE_BOTHER=1; GLOBAL_PTREE_BITBYBIT=1; GLOBAL_PTREE_REGION_TYPE=1;
	  GLOBAL_PTREE = ptreemake(GLOBAL_PTREE_NREGIONS,GLOBAL_PTREE_EVENT_WITHIN,GLOBAL_PTREE_EVENT_THRESHOLD,GLOBAL_PTREE_REGION_TYPE,GLOBAL_PTREE_NLEGS,GLOBAL_PTREE_LEGTIME);
	  if (POWER_BOTHER){ 
	    if (verbose){ printf(" %% %% dumping power\n");} 
	    powerdump(GLOBAL_POWER,"rate3",0); powerdump(GLOBAL_POWER,"rate3",1);}
	  if (GLOBAL_POWER!=NULL){ powertfree(GLOBAL_POWER); GLOBAL_POWER=NULL;}
	  POWER_BOTHER=1; GLOBAL_POWER = powermake(GLOBAL_Nra,GLOBAL_POWER_LENGTH,GLOBAL_POWER_HOWMANY);}}
      if (CORTEX_BOTHER){
	if (OUTPUT_DUMP_EVERY>0 && (int)floor((GLOBAL_time-time_3)/OUTPUT_DUMP_EVERY)<(int)floor((GLOBAL_time-time_3+GLOBAL_DT)/OUTPUT_DUMP_EVERY)){
	  if (PTREE_BOTHER){
	    if (verbose){ printf(" %% %% dumping ptree\n");}
	    ptreerate(GLOBAL_PTREE);
	    ptreedump_starter(GLOBAL_PTREE,"rate4",2,0,0,0,+1,-1);
	    if (GLOBAL_PTREE_BITBYBIT){
	      if (verbose){ printf(" %% %% resetting ptree\n");}
	      ptreereset(GLOBAL_PTREE);}}}}}
    else if (t>=time_4 && t<time_5){ /* rate5 */
      if (STIMULUS_ONSET_TIME != time_4){
	if (verbose){ printf(" %% entering fifth phase of suite seven -- rate5\n");}
	GRATING_VS_LMI = -10; LGN_DUMB=2;
	STIMULUS_ONSET_TIME = time_4; OUTPUT_DUMP_EVERY = time_dumpevery; /* start dumping data */
	if (CORTEX_BOTHER){
	  if (GLOBAL_PTREE!=NULL){ ptreetfree(GLOBAL_PTREE); GLOBAL_PTREE=NULL;} GLOBAL_PTREE_BITBYBIT=0; PTREE_BOTHER=0;
	  PTREE_BOTHER=1; GLOBAL_PTREE_BITBYBIT=1; GLOBAL_PTREE_REGION_TYPE=1;
	  GLOBAL_PTREE = ptreemake(GLOBAL_PTREE_NREGIONS,GLOBAL_PTREE_EVENT_WITHIN,GLOBAL_PTREE_EVENT_THRESHOLD,GLOBAL_PTREE_REGION_TYPE,GLOBAL_PTREE_NLEGS,GLOBAL_PTREE_LEGTIME);
	  if (POWER_BOTHER){ 
	    if (verbose){ printf(" %% %% dumping power\n");} 
	    powerdump(GLOBAL_POWER,"rate4",0); powerdump(GLOBAL_POWER,"rate4",1);}
	  if (GLOBAL_POWER!=NULL){ powertfree(GLOBAL_POWER); GLOBAL_POWER=NULL;}
	  POWER_BOTHER=1; GLOBAL_POWER = powermake(GLOBAL_Nra,GLOBAL_POWER_LENGTH,GLOBAL_POWER_HOWMANY);}}
      if (CORTEX_BOTHER){
	if (OUTPUT_DUMP_EVERY>0 && (int)floor((GLOBAL_time-time_4)/OUTPUT_DUMP_EVERY)<(int)floor((GLOBAL_time-time_4+GLOBAL_DT)/OUTPUT_DUMP_EVERY)){
	  if (PTREE_BOTHER){
	    if (verbose){ printf(" %% %% dumping ptree\n");}
	    ptreerate(GLOBAL_PTREE);
	    ptreedump_starter(GLOBAL_PTREE,"rate5",2,0,0,0,+1,-1);
	    if (GLOBAL_PTREE_BITBYBIT){
	      if (verbose){ printf(" %% %% resetting ptree\n");}
	      ptreereset(GLOBAL_PTREE);}}}}}
    else if (t>=time_5 && t<time_6){ /* rate6 */
      if (STIMULUS_ONSET_TIME != time_5){
	if (verbose){ printf(" %% entering sixth phase of suite seven -- rate6\n");}
	GRATING_VS_LMI = -10; LGN_DUMB=3;
	STIMULUS_ONSET_TIME = time_5; OUTPUT_DUMP_EVERY = time_dumpevery; /* start dumping data */
	if (CORTEX_BOTHER){
	  if (GLOBAL_PTREE!=NULL){ ptreetfree(GLOBAL_PTREE); GLOBAL_PTREE=NULL;} GLOBAL_PTREE_BITBYBIT=0; PTREE_BOTHER=0;
	  PTREE_BOTHER=1; GLOBAL_PTREE_BITBYBIT=1; GLOBAL_PTREE_REGION_TYPE=1;
	  GLOBAL_PTREE = ptreemake(GLOBAL_PTREE_NREGIONS,GLOBAL_PTREE_EVENT_WITHIN,GLOBAL_PTREE_EVENT_THRESHOLD,GLOBAL_PTREE_REGION_TYPE,GLOBAL_PTREE_NLEGS,GLOBAL_PTREE_LEGTIME);
	  if (POWER_BOTHER){ 
	    if (verbose){ printf(" %% %% dumping power\n");} 
	    powerdump(GLOBAL_POWER,"rate5",0); powerdump(GLOBAL_POWER,"rate5",1);}
	  if (GLOBAL_POWER!=NULL){ powertfree(GLOBAL_POWER); GLOBAL_POWER=NULL;}
	  POWER_BOTHER=1; GLOBAL_POWER = powermake(GLOBAL_Nra,GLOBAL_POWER_LENGTH,GLOBAL_POWER_HOWMANY);}}
      if (CORTEX_BOTHER){
	if (OUTPUT_DUMP_EVERY>0 && (int)floor((GLOBAL_time-time_5)/OUTPUT_DUMP_EVERY)<(int)floor((GLOBAL_time-time_5+GLOBAL_DT)/OUTPUT_DUMP_EVERY)){
	  if (PTREE_BOTHER){
	    if (verbose){ printf(" %% %% dumping ptree\n");}
	    ptreerate(GLOBAL_PTREE);
	    ptreedump_starter(GLOBAL_PTREE,"rate6",2,0,0,0,+1,-1);
	    if (GLOBAL_PTREE_BITBYBIT){
	      if (verbose){ printf(" %% %% resetting ptree\n");}
	      ptreereset(GLOBAL_PTREE);}}}}}
    else if (t>=time_6 && t<time_7){ /* rate7 */
      if (STIMULUS_ONSET_TIME != time_6){
	if (verbose){ printf(" %% entering seventh phase of suite seven -- rate7\n");}
	GRATING_VS_LMI = -10; LGN_DUMB=4;
	STIMULUS_ONSET_TIME = time_6; OUTPUT_DUMP_EVERY = time_dumpevery; /* start dumping data */
	if (CORTEX_BOTHER){
	  if (GLOBAL_PTREE!=NULL){ ptreetfree(GLOBAL_PTREE); GLOBAL_PTREE=NULL;} GLOBAL_PTREE_BITBYBIT=0; PTREE_BOTHER=0;
	  PTREE_BOTHER=1; GLOBAL_PTREE_BITBYBIT=1; GLOBAL_PTREE_REGION_TYPE=1;
	  GLOBAL_PTREE = ptreemake(GLOBAL_PTREE_NREGIONS,GLOBAL_PTREE_EVENT_WITHIN,GLOBAL_PTREE_EVENT_THRESHOLD,GLOBAL_PTREE_REGION_TYPE,GLOBAL_PTREE_NLEGS,GLOBAL_PTREE_LEGTIME);
	  if (POWER_BOTHER){ 
	    if (verbose){ printf(" %% %% dumping power\n");} 
	    powerdump(GLOBAL_POWER,"rate6",0); powerdump(GLOBAL_POWER,"rate6",1);}
	  if (GLOBAL_POWER!=NULL){ powertfree(GLOBAL_POWER); GLOBAL_POWER=NULL;}
	  POWER_BOTHER=1; GLOBAL_POWER = powermake(GLOBAL_Nra,GLOBAL_POWER_LENGTH,GLOBAL_POWER_HOWMANY);}}
      if (CORTEX_BOTHER){
	if (OUTPUT_DUMP_EVERY>0 && (int)floor((GLOBAL_time-time_6)/OUTPUT_DUMP_EVERY)<(int)floor((GLOBAL_time-time_6+GLOBAL_DT)/OUTPUT_DUMP_EVERY)){
	  if (PTREE_BOTHER){
	    if (verbose){ printf(" %% %% dumping ptree\n");}
	    ptreerate(GLOBAL_PTREE);
	    ptreedump_starter(GLOBAL_PTREE,"rate7",2,0,0,0,+1,-1);
	    if (GLOBAL_PTREE_BITBYBIT){
	      if (verbose){ printf(" %% %% resetting ptree\n");}
	      ptreereset(GLOBAL_PTREE);}}}}}
    else if (t>=time_7){ /* finish */
      if (STIMULUS_ONSET_TIME != time_7){
	if (verbose){ printf(" %% entering eight phase of suite seven -- finishing\n");}
	GRATING_VS_LMI=-10;
	STIMULUS_ONSET_TIME=time_7;
	GLOBAL_TF = time_7+1;
	if (CORTEX_BOTHER){
	  connectionsdump(GLOBAL_Nra,1,"connections");
	  if (GLOBAL_PTREE!=NULL){ ptreetfree(GLOBAL_PTREE); GLOBAL_PTREE=NULL;} GLOBAL_PTREE_BITBYBIT=0; PTREE_BOTHER=0;
	  if (POWER_BOTHER){ 
	    if (verbose){ printf(" %% %% dumping power\n");} 
	    powerdump(GLOBAL_POWER,"rate7",0); powerdump(GLOBAL_POWER,"rate7",1);}
	  if (GLOBAL_POWER!=NULL){ powertfree(GLOBAL_POWER); GLOBAL_POWER=NULL;} POWER_BOTHER=0;}
	fnamebase = (char **) tcalloc(7,sizeof(char *));
	fnamebase[0] = (char *) tcalloc(256,sizeof(char));
	fnamebase[1] = (char *) tcalloc(256,sizeof(char)); 
	for (nr1=1;nr1<=7;nr1++){ for (nr2=nr1+1;nr2<=7;nr2++){
	  sprintf(fnamebase[0],"ptree_rate%d_%srecord",nr1,gs2);
	  sprintf(fnamebase[1],"ptree_rate%d_%srecord",nr2,gs2);
	  for (tindex=1;tindex<=SUITE_TINDEXMAX;tindex++){ for (sindex=0;sindex<=SUITE_SINDEXMAX;sindex++){ for (dindex=0;dindex<=SUITE_DINDEXMAX;dindex++){
	    sprintf(outputname,"r%d%d",nr1,nr2);
	    ptree_trialaverage_helper(2,GLOBAL_PTREE_NLEGS+1,tindex,1,1,-1,sindex,dindex,&fnamebase[0],outputname);}}}}}
	sprintf(fnamebase[0],"ptree_rate1_%srecord",gs2);
	sprintf(fnamebase[1],"ptree_rate2_%srecord",gs2);
	fnamebase[2] = (char *) tcalloc(256,sizeof(char)); sprintf(fnamebase[2],"ptree_rate3_%srecord",gs2);
	fnamebase[3] = (char *) tcalloc(256,sizeof(char)); sprintf(fnamebase[3],"ptree_rate4_%srecord",gs2);
	fnamebase[4] = (char *) tcalloc(256,sizeof(char)); sprintf(fnamebase[4],"ptree_rate5_%srecord",gs2);
	fnamebase[5] = (char *) tcalloc(256,sizeof(char)); sprintf(fnamebase[5],"ptree_rate6_%srecord",gs2);
	fnamebase[6] = (char *) tcalloc(256,sizeof(char)); sprintf(fnamebase[6],"ptree_rate7_%srecord",gs2);
	for (tindex=1;tindex<=SUITE_TINDEXMAX;tindex++){ for (sindex=0;sindex<=SUITE_SINDEXMAX;sindex++){ for (dindex=0;dindex<=SUITE_DINDEXMAX;dindex++){
	  sprintf(outputname,"rALL"); ptree_trialaverage_helper(7,GLOBAL_PTREE_NLEGS+1,tindex,1,1,-1,sindex,dindex,&fnamebase[0],outputname);}}}
	tfree(fnamebase[0]);fnamebase[0]=NULL;
	tfree(fnamebase[1]);fnamebase[1]=NULL;
	tfree(fnamebase[2]);fnamebase[2]=NULL;
	tfree(fnamebase[3]);fnamebase[3]=NULL;
	tfree(fnamebase[4]);fnamebase[4]=NULL;
	tfree(fnamebase[5]);fnamebase[5]=NULL;
	tfree(fnamebase[6]);fnamebase[6]=NULL;
	tfree(fnamebase);fnamebase=NULL;
	if (SUITE_BITBYBIT_REMOVE){ 
	  for (nr=1;nr<(int)floor((double)SUITE_NSECONDS*1024/(double)SUITE_DUMPEVERY)+1;nr++){
	    sprintf(command,"""rm"" ./ptree_rate1_%srecord_%d ;",gs2,nr); system(command);
	    sprintf(command,"""rm"" ./ptree_rate2_%srecord_%d ;",gs2,nr); system(command);
	    sprintf(command,"""rm"" ./ptree_rate3_%srecord_%d ;",gs2,nr); system(command);
	    sprintf(command,"""rm"" ./ptree_rate4_%srecord_%d ;",gs2,nr); system(command);
	    sprintf(command,"""rm"" ./ptree_rate5_%srecord_%d ;",gs2,nr); system(command);
	    sprintf(command,"""rm"" ./ptree_rate6_%srecord_%d ;",gs2,nr); system(command);
	    sprintf(command,"""rm"" ./ptree_rate7_%srecord_%d ;",gs2,nr); system(command);}}}}
    break;
  case 8:
    time_0=0;
    time_1=time_0+1024*SUITE_NSECONDS;
    if (verbose>2){ printf(" times are: %0.2f, %0.2f, \n",time_0,time_1);}
    if (t==0){ /* clear everything */
      if (verbose){ printf(" %% \n");}
      if (verbose){ printf(" %% initializing suite eight\n");}
      GLOBAL_TF=0;
      STIMULUS_ONSET_TIME=-1;
      if (GRATING_VS_LMI!=-10){ GRATING_VS_LMI = -10;}
      if (CORTEX_BOTHER){
	if (RTC_BOTHER!=0){ rtctfree(GLOBAL_RTC); GLOBAL_RTC=NULL; RTC_BOTHER=0;}
	if (STROBETRACE_BOTHER!=0){ strobetracetfree(GLOBAL_STROBETRACE); GLOBAL_STROBETRACE=NULL; STROBETRACE_BOTHER=0;}
	if (TUNINGCURVE_BOTHER!=0){ tuningcurvetfree(GLOBAL_TUNINGCURVE); GLOBAL_TUNINGCURVE=NULL; TUNINGCURVE_BOTHER=0;}
	if (LMITRI_BOTHER!=0){ lmitritfree(GLOBAL_LMITRI); GLOBAL_LMITRI=NULL; LMITRI_BOTHER=0;}
	if (PTREE_BOTHER!=0){ ptreetfree(GLOBAL_PTREE); GLOBAL_PTREE=NULL; GLOBAL_PTREE_BITBYBIT=0; PTREE_BOTHER=0;}
	if (CLOSET_BOTHER!=0){ closettfree(GLOBAL_CLOSET); GLOBAL_CLOSET=NULL; CLOSET_BOTHER=0;}
	if (POWER_BOTHER!=0){ powertfree(GLOBAL_POWER); GLOBAL_POWER=NULL; POWER_BOTHER=0;}
	if (YGGDRASIL_BOTHER!=0){ yggdrasiltfree(GLOBAL_YGGDRASIL); GLOBAL_YGGDRASIL=NULL; YGGDRASIL_BOTHER=0;}
	if (BONSAI_BOTHER!=0){ bonsaitfree(GLOBAL_BONSAI); GLOBAL_BONSAI=NULL; BONSAI_BOTHER=0;}
	if (HYDRA_BOTHER!=0){ hydratfree(GLOBAL_HYDRA); GLOBAL_HYDRA=NULL; HYDRA_BOTHER=0;}
	if (LYAPUNOV_BOTHER!=0){ lyapunovtfree(GLOBAL_LYAPUNOV); GLOBAL_LYAPUNOV=NULL; LYAPUNOV_BOTHER=0;}}}
    else if (t>=0 && t<time_1){ /* rate1vsrate2step  */
      if (STIMULUS_ONSET_TIME != time_0){
	if (verbose){ printf(" %% entering first phase of suite eight -- rate1vsrate2step\n");}
	GRATING_VS_LMI=0;
	STIMULUS_ONSET_TIME = time_0; /* start dumping data within hydra */
	if (CORTEX_BOTHER){
	  HYDRA_BOTHER=1; GLOBAL_PTREE_BITBYBIT=1; /* GLOBAL_PTREE_REGION_TYPE=1; */
	  GLOBAL_HYDRA = hydramake(GLOBAL_HYDRA_JUSTONTIME,GLOBAL_HYDRA_STAYONTIME,GLOBAL_HYDRA_DUMP_EVERY,&ptree_inputswitchon,&ptree_inputswitchoff);}}}
    else if (t>=time_1){ /* finish */
      if (STIMULUS_ONSET_TIME != time_1){
	if (verbose){ printf(" %% entering second phase of suite eight -- finishing\n");}
	GRATING_VS_LMI=-10;
	STIMULUS_ONSET_TIME=time_1;
	GLOBAL_TF = time_1+1;
	if (CORTEX_BOTHER){
	  connectionsdump(GLOBAL_Nra,1,"connections");
	  if (GLOBAL_HYDRA!=NULL){ hydradump(GLOBAL_HYDRA); hydratfree(GLOBAL_HYDRA); GLOBAL_HYDRA=NULL; HYDRA_BOTHER=0;}}}}
    break;
  case 9:
    ntimes=4+1;
    timera = (double *) tcalloc(ntimes,sizeof(double));
    for (nt=0;nt<ntimes;nt++){ timera[nt] = nt*1024*SUITE_NSECONDS;}
    time_dumpevery = SUITE_DUMPEVERY;
    if (verbose>2){ printf(" times are: ");raprintf(timera,"double",1,ntimes,"timera ");}
    if (t==0){ /* clear everything */
      if (verbose){ printf(" %% \n");}
      if (verbose){ printf(" %% initializing suite nine\n");}
      GLOBAL_TF=0;
      STIMULUS_ONSET_TIME=-1;
      if (GRATING_VS_LMI!=-10){ GRATING_VS_LMI = -10;}
      if (CORTEX_BOTHER){
	if (RTC_BOTHER!=0){ rtctfree(GLOBAL_RTC); GLOBAL_RTC=NULL; RTC_BOTHER=0;}
	if (STROBETRACE_BOTHER!=0){ strobetracetfree(GLOBAL_STROBETRACE); GLOBAL_STROBETRACE=NULL; STROBETRACE_BOTHER=0;}
	if (TUNINGCURVE_BOTHER!=0){ tuningcurvetfree(GLOBAL_TUNINGCURVE); GLOBAL_TUNINGCURVE=NULL; TUNINGCURVE_BOTHER=0;}
	if (LMITRI_BOTHER!=0){ lmitritfree(GLOBAL_LMITRI); GLOBAL_LMITRI=NULL; LMITRI_BOTHER=0;}
	if (PTREE_BOTHER!=0){ ptreetfree(GLOBAL_PTREE); GLOBAL_PTREE=NULL; GLOBAL_PTREE_BITBYBIT=0; PTREE_BOTHER=0;}
	if (CLOSET_BOTHER!=0){ closettfree(GLOBAL_CLOSET); GLOBAL_CLOSET=NULL; CLOSET_BOTHER=0;}
	if (POWER_BOTHER!=0){ powertfree(GLOBAL_POWER); GLOBAL_POWER=NULL; POWER_BOTHER=0;}
	if (YGGDRASIL_BOTHER!=0){ yggdrasiltfree(GLOBAL_YGGDRASIL); GLOBAL_YGGDRASIL=NULL; YGGDRASIL_BOTHER=0;}
	if (BONSAI_BOTHER!=0){ bonsaitfree(GLOBAL_BONSAI); GLOBAL_BONSAI=NULL; BONSAI_BOTHER=0;}
	if (HYDRA_BOTHER!=0){ hydratfree(GLOBAL_HYDRA); GLOBAL_HYDRA=NULL; HYDRA_BOTHER=0;}
	if (LYAPUNOV_BOTHER!=0){ lyapunovtfree(GLOBAL_LYAPUNOV); GLOBAL_LYAPUNOV=NULL; LYAPUNOV_BOTHER=0;}}}
    else /* if (t>0) */{
      nowtime = 0; for (nt=0;nt<ntimes-1;nt++){ if (t>=timera[nt] && t<timera[nt+1]){ nowtime=nt;}} 
      if (t>=timera[ntimes-1]){ nowtime=ntimes-1;}
      if (nowtime<ntimes-1){
	if (STIMULUS_ONSET_TIME != timera[nowtime]){
	  if (verbose){ printf(" %% entering %d_th phase of suite nine... ",nowtime);}
	  GRATING_VS_LMI = -10; LGN_DUMB = -1; GLOBAL_ODOR += (nowtime>0)*GLOBAL_ODOR_BACON;
	  if (verbose){ printf("setting LGN_DUMB=%d,GLOBAL_ODOR=%d\n",LGN_DUMB,GLOBAL_ODOR);}
	  STIMULUS_ONSET_TIME = timera[nowtime]; OUTPUT_DUMP_EVERY = time_dumpevery; /* start dumping data */
	  if (CORTEX_BOTHER){
	    if (GLOBAL_PTREE!=NULL){ ptreetfree(GLOBAL_PTREE); GLOBAL_PTREE=NULL;} GLOBAL_PTREE_BITBYBIT=0; PTREE_BOTHER=0;
	    PTREE_BOTHER=1; GLOBAL_PTREE_BITBYBIT=1; GLOBAL_PTREE_REGION_TYPE=1;
	    GLOBAL_PTREE = ptreemake(GLOBAL_PTREE_NREGIONS,GLOBAL_PTREE_EVENT_WITHIN,GLOBAL_PTREE_EVENT_THRESHOLD,GLOBAL_PTREE_REGION_TYPE,GLOBAL_PTREE_NLEGS,GLOBAL_PTREE_LEGTIME);}}
	if (CORTEX_BOTHER){
	  if (OUTPUT_DUMP_EVERY>0 && (int)floor((GLOBAL_time-timera[nowtime])/OUTPUT_DUMP_EVERY)<(int)floor((GLOBAL_time-timera[nowtime]+GLOBAL_DT)/OUTPUT_DUMP_EVERY)){
	    if (PTREE_BOTHER){
	      if (verbose){ printf(" %% %% dumping ptree\n");}
	      ptreerate(GLOBAL_PTREE); 
	      sprintf(tempchar,"rate%d",nowtime); ptreedump_starter(GLOBAL_PTREE,tempchar,2,0,0,0,+1,-1);
	      if (GLOBAL_PTREE_BITBYBIT){
		if (verbose){ printf(" %% %% resetting ptree\n");}
		ptreereset(GLOBAL_PTREE);}}}}}
      else /* if (nowtime==ntimes-1) */{
	if (STIMULUS_ONSET_TIME != timera[nowtime]){
	  if (verbose){ printf(" %% entering %d_th phase of suite nine -- finishing\n",nowtime);}
	  GRATING_VS_LMI=-10;
	  STIMULUS_ONSET_TIME=timera[nowtime];
	  GLOBAL_TF = timera[nowtime]+1;
	  if (CORTEX_BOTHER){
	    connectionsdump(GLOBAL_Nra,1,"connections");
	    if (GLOBAL_PTREE!=NULL){ ptreetfree(GLOBAL_PTREE); GLOBAL_PTREE=NULL;} GLOBAL_PTREE_BITBYBIT=0; PTREE_BOTHER=0;}
	  fnamebase = (char **) tcalloc(2,sizeof(char *));
	  fnamebase[0] = (char *) tcalloc(256,sizeof(char));
	  fnamebase[1] = (char *) tcalloc(256,sizeof(char)); 
	  for (nr1=0;nr1<=ntimes-2;nr1++){ for (nr2=nr1+1;nr2<=ntimes-2;nr2++){
	    sprintf(fnamebase[0],"ptree_rate%d_%srecord",nr1,gs2);
	    sprintf(fnamebase[1],"ptree_rate%d_%srecord",nr2,gs2);
	    for (tindex=1;tindex<=SUITE_TINDEXMAX;tindex++){ for (sindex=0;sindex<=SUITE_SINDEXMAX;sindex++){ for (dindex=0;dindex<=SUITE_DINDEXMAX;dindex++){
	    sprintf(outputname,"r%d%d",nr1,nr2);
	    ptree_trialaverage_helper(2,GLOBAL_PTREE_NLEGS+1,tindex,1,1,-1,sindex,dindex,&fnamebase[0],outputname);}}}}}
	  tfree(fnamebase[0]);fnamebase[0]=NULL;
	  tfree(fnamebase[1]);fnamebase[1]=NULL;
	  tfree(fnamebase);fnamebase=NULL;
	  if (SUITE_BITBYBIT_REMOVE){ 
	    for (nr=1;nr<(int)floor((double)SUITE_NSECONDS*1024/(double)SUITE_DUMPEVERY)+1;nr++){
	      for (nt=0;nt<ntimes-1;nt++){
		sprintf(command,"""rm"" ./ptree_rate%d_%srecord_%d ;",nt,gs2,nr); system(command);}}}}}}
    tfree(timera);timera=NULL;
    break;
  case 10:
    time_0=0;
    time_1=time_0+2048;
    time_2=time_1+512;
    time_delay=16; /* fix for lmi */
    time_3=time_2+512;
    time_4=time_3+512;
    time_5=time_4+512;
    time_6=time_5+512;
    taof_flag=1;
    time_7=time_6+512;
    suite_10_idslmi_flag=0;
    if (!suite_10_idslmi_flag){ time_7=time_0;}
    ntimes=GLOBAL_TAOF_LENGTH+2;
    timera = (double *) tcalloc(ntimes,sizeof(double));
    for (nt=0;nt<ntimes;nt++){ timera[nt] = time_7 + nt*GLOBAL_TAOF_STEP_EVERY;}
    if (verbose && t==0){ printf(" times are: %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f\n",time_0,time_1,time_2,time_3,time_4,time_5,time_6); printf(" times are: ");raprintf(timera,"double",1,ntimes,"timera ");}
    if (t==0){ /* clear everything */
      if (verbose){ printf(" %% \n");}
      if (verbose){ printf(" %% entering first phase of suite ten -- background\n");}
      GLOBAL_TF=0;
      if (GRATING_VS_LMI!=-10){ GRATING_VS_LMI = -10;}
      if (CORTEX_BOTHER){
	if (RTC_BOTHER!=0){ rtctfree(GLOBAL_RTC); GLOBAL_RTC=NULL; RTC_BOTHER=0;}
	if (STROBETRACE_BOTHER!=0){ strobetracetfree(GLOBAL_STROBETRACE); GLOBAL_STROBETRACE=NULL; STROBETRACE_BOTHER=0;}
	if (TUNINGCURVE_BOTHER!=0){ tuningcurvetfree(GLOBAL_TUNINGCURVE); GLOBAL_TUNINGCURVE=NULL; TUNINGCURVE_BOTHER=0;}
	if (LMITRI_BOTHER!=0){ lmitritfree(GLOBAL_LMITRI); GLOBAL_LMITRI=NULL; LMITRI_BOTHER=0;}
	if (PTREE_BOTHER!=0){ ptreetfree(GLOBAL_PTREE); GLOBAL_PTREE=NULL; PTREE_BOTHER=0;}
	if (CLOSET_BOTHER!=0){ closettfree(GLOBAL_CLOSET); GLOBAL_CLOSET=NULL; CLOSET_BOTHER=0;}
	if (POWER_BOTHER!=0){ powertfree(GLOBAL_POWER); GLOBAL_POWER=NULL; POWER_BOTHER=0;}
	if (TAOF_BOTHER!=0){ taoftfree(GLOBAL_TAOF); GLOBAL_TAOF=NULL; TAOF_BOTHER=0;}
	if (YGGDRASIL_BOTHER!=0){ yggdrasiltfree(GLOBAL_YGGDRASIL); GLOBAL_YGGDRASIL=NULL; YGGDRASIL_BOTHER=0;}
	if (BONSAI_BOTHER!=0){ bonsaitfree(GLOBAL_BONSAI); GLOBAL_BONSAI=NULL; BONSAI_BOTHER=0;}
	if (HYDRA_BOTHER!=0){ hydratfree(GLOBAL_HYDRA); GLOBAL_HYDRA=NULL; HYDRA_BOTHER=0;}
	if (LYAPUNOV_BOTHER!=0){ lyapunovtfree(GLOBAL_LYAPUNOV); GLOBAL_LYAPUNOV=NULL; LYAPUNOV_BOTHER=0;}}}
/*     else if (t>=0 && t<time_1){ /\* background *\/ */
/*       if (GRATING_VS_LMI!=-10){ GRATING_VS_LMI = -10;} */
/*       if (CORTEX_BOTHER && GLOBAL_STROBETRACE==NULL){ */
/* 	if (verbose){ printf(" %% %% making strobetrace\n");} */
/* 	STROBETRACE_BOTHER=1; */
/* 	GLOBAL_STROBETRACE = strobetracemake(nget(GLOBAL_Nra,0,0),GLOBAL_STROBETRACE_NANGLES,time_1-time_0,1.0,0,GLOBAL_STROBETRACE_AVALANCHE_UPDATE_EVERY,GLOBAL_STROBETRACE_AVALANCHE_LOGDMIN,GLOBAL_STROBETRACE_AVALANCHE_LOGDMAX);}} */
/*     else if (t>=time_1 && t<time_2){ /\* line motion illusion *\/ */
/*       if (GRATING_VS_LMI!=1){  */
/* 	if (verbose){ printf(" %% entering second phase of suite ten -- line motion illusion\n");} */
/* 	if (CORTEX_BOTHER && STROBETRACE_BOTHER){ */
/* 	  if (verbose){ printf(" %% %% dumping strobetrace\n");} */
/* 	  strobetracedump(GLOBAL_STROBETRACE,0);  */
/* 	  strobetracedump(GLOBAL_STROBETRACE,1);  */
/* 	  strobetracetfree(GLOBAL_STROBETRACE); GLOBAL_STROBETRACE=NULL; STROBETRACE_BOTHER=0;} */
/* 	GRATING_VS_LMI = 1; STIMULUS_ONSET_TIME = time_1+time_delay;} */
/*       if (CORTEX_BOTHER){ */
/* 	if (LMITRI_BOTHER==0){  */
/* 	  if (verbose){ printf(" %% %% making lmitri\n");} */
/* 	  LMITRI_BOTHER=1; */
/* 	  GLOBAL_LMITRI = lmitrimake(STIMULUS_ONSET_TIME,GLOBAL_LMITRI_TIMELENGTH,GLOBAL_LMITRI_ROW_MAX,GLOBAL_LMITRI_ROW_MIN);}}} */
/*     else if (t>=time_2 && t<time_3){ /\* another line motion illusion *\/ */
/*       if (STIMULUS_ONSET_TIME!=time_2+time_delay){ */
/* 	if (verbose){ printf(" %% entering third phase of suite ten -- line motion illusion\n");} */
/* 	STIMULUS_ONSET_TIME = time_2+time_delay;} */
/*       if (CORTEX_BOTHER){  */
/* 	if (GLOBAL_LMITRI->time_start!=STIMULUS_ONSET_TIME){ */
/* 	  GLOBAL_LMITRI->time_start = STIMULUS_ONSET_TIME;}}} */
/*     else if (t>=time_3 && t<time_4){ /\* yet another line motion illusion *\/ */
/*       if (STIMULUS_ONSET_TIME!=time_3+time_delay){ */
/* 	if (verbose){ printf(" %% entering fourth phase of suite ten -- line motion illusion\n");} */
/* 	STIMULUS_ONSET_TIME = time_3+time_delay;} */
/*       if (CORTEX_BOTHER){  */
/* 	if (GLOBAL_LMITRI->time_start!=STIMULUS_ONSET_TIME){ */
/* 	  GLOBAL_LMITRI->time_start = STIMULUS_ONSET_TIME;}}} */
/*     else if (t>=time_4 && t<time_5){ /\* and another line motion illusion *\/ */
/*       if (STIMULUS_ONSET_TIME!=time_4+time_delay){ */
/* 	if (verbose){ printf(" %% entering fifth phase of suite ten -- line motion illusion\n");} */
/* 	STIMULUS_ONSET_TIME = time_4+time_delay;} */
/*       if (CORTEX_BOTHER){  */
/* 	if (GLOBAL_LMITRI->time_start!=STIMULUS_ONSET_TIME){ */
/* 	  GLOBAL_LMITRI->time_start = STIMULUS_ONSET_TIME;}}} */
/*     else if (t>=time_5 && t<time_6){ /\* and a final line motion illusion *\/ */
/*       if (STIMULUS_ONSET_TIME!=time_5+time_delay){ */
/* 	if (verbose){ printf(" %% entering sixth phase of suite ten -- line motion illusion\n");} */
/* 	STIMULUS_ONSET_TIME = time_5+time_delay;} */
/*       if (CORTEX_BOTHER){  */
/* 	if (GLOBAL_LMITRI->time_start!=STIMULUS_ONSET_TIME){ */
/* 	  GLOBAL_LMITRI->time_start = STIMULUS_ONSET_TIME;}}} */
/*     else if (t>=time_6 && t<time_7){ /\* wait *\/ */
/*       if (fabs(STIMULUS_ONSET_TIME-time_6)>0.5){ */
/* 	if (!taof_flag){ GLOBAL_TF = time_6+1;} */
/* 	STIMULUS_ONSET_TIME = time_6; */
/* 	if (verbose){ printf(" %% entering seventh phase of suite ten -- waiting\n");} */
/* 	GRATING_VS_LMI = -10; */
/* 	if (CORTEX_BOTHER){ */
/* 	  if (LMITRI_BOTHER){ */
/* 	    if (verbose){ printf(" %% %% dumping lmitri\n");} */
/* 	    lmitridump(GLOBAL_LMITRI,1); lmitritfree(GLOBAL_LMITRI); GLOBAL_LMITRI=NULL; LMITRI_BOTHER=0;}} */
/* 	time_dumpevery = GLOBAL_TAOF_STEP_EVERY; */
/* 	if (CORTEX_BOTHER){ */
/* 	  if (RTC_BOTHER!=0){ rtctfree(GLOBAL_RTC); GLOBAL_RTC=NULL; RTC_BOTHER=0;} */
/* 	  if (STROBETRACE_BOTHER!=0){ strobetracetfree(GLOBAL_STROBETRACE); GLOBAL_STROBETRACE=NULL; STROBETRACE_BOTHER=0;} */
/* 	  if (TUNINGCURVE_BOTHER!=0){ tuningcurvetfree(GLOBAL_TUNINGCURVE); GLOBAL_TUNINGCURVE=NULL; TUNINGCURVE_BOTHER=0;} */
/* 	  if (LMITRI_BOTHER!=0){ lmitritfree(GLOBAL_LMITRI); GLOBAL_LMITRI=NULL; LMITRI_BOTHER=0;} */
/* 	  if (PTREE_BOTHER!=0){ ptreetfree(GLOBAL_PTREE); GLOBAL_PTREE=NULL; GLOBAL_PTREE_BITBYBIT=0; PTREE_BOTHER=0;} */
/* 	  if (CLOSET_BOTHER!=0){ closettfree(GLOBAL_CLOSET); GLOBAL_CLOSET=NULL; CLOSET_BOTHER=0;} */
/* 	  if (POWER_BOTHER!=0){ powertfree(GLOBAL_POWER); GLOBAL_POWER=NULL; POWER_BOTHER=0;} */
/* 	  if (TAOF_BOTHER!=0){ taoftfree(GLOBAL_TAOF); GLOBAL_TAOF=NULL; TAOF_BOTHER=0;} */
/* 	  if (YGGDRASIL_BOTHER!=0){ yggdrasiltfree(GLOBAL_YGGDRASIL); GLOBAL_YGGDRASIL=NULL; YGGDRASIL_BOTHER=0;} */
/* 	  if (BONSAI_BOTHER!=0){ bonsaitfree(GLOBAL_BONSAI); GLOBAL_BONSAI=NULL; BONSAI_BOTHER=0;} */
/* 	  if (HYDRA_BOTHER!=0){ hydratfree(GLOBAL_HYDRA); GLOBAL_HYDRA=NULL; HYDRA_BOTHER=0;} */
    /* 	  if (LYAPUNOV_BOTHER!=0){ lyapunovtfree(GLOBAL_LYAPUNOV); GLOBAL_LYAPUNOV=NULL; LYAPUNOV_BOTHER=0;}}}} */
    else if (t>=time_7){ /* taof */
      nowtime = 0; for (nt=0;nt<ntimes-1;nt++){ if (t>=timera[nt] && t<timera[nt+1]){ nowtime=nt;}} 
      if (t>=timera[ntimes-1]){ nowtime=ntimes-1;}
      if (nowtime<ntimes-1){
	if (fabs(STIMULUS_ONSET_TIME-timera[nowtime])>0.5){
	  if (verbose){ printf(" %% entering %d_th phase of suite ten... ",nowtime+8*suite_10_idslmi_flag);}
	  GRATING_VS_LMI = 0; INPUT_CONTRAST = 0.75*(nowtime>0)*(1-2*abs((ntimes-1)/2-nowtime)/(double)(ntimes-1)); /* 0.75 max */
	  if (verbose){ printf("setting INPUT_CONTRAST=%0.4f\n",INPUT_CONTRAST);}
	  STIMULUS_ONSET_TIME = timera[nowtime]; OUTPUT_DUMP_EVERY = time_dumpevery; /* start dumping data */
	  if (CORTEX_BOTHER){
	    if (GLOBAL_TAOF==NULL){ GLOBAL_TAOF = taofmake(GLOBAL_Nra,GLOBAL_TAOF_LENGTH,GLOBAL_TAOF_STEP_EVERY); TAOF_BOTHER=1;}}}
	if (CORTEX_BOTHER){
	  if (OUTPUT_DUMP_EVERY>0 && (int)floor((GLOBAL_time-timera[nowtime])/OUTPUT_DUMP_EVERY)<(int)floor((GLOBAL_time-timera[nowtime]+GLOBAL_DT)/OUTPUT_DUMP_EVERY)){
	    if (GLOBAL_TAOF){ taofdump(GLOBAL_TAOF,NULL,3);taofdump(GLOBAL_TAOF,NULL,4);}}}}
      else /* if (nowtime==ntimes-1) */{
	if (STIMULUS_ONSET_TIME != timera[nowtime]){
	  if (verbose){ printf(" %% entering %d_th phase of suite ten -- finishing\n",nowtime);}
	  GRATING_VS_LMI=-10;
	  STIMULUS_ONSET_TIME=timera[nowtime];
	  GLOBAL_TF = timera[nowtime]+1;
	  if (CORTEX_BOTHER){ if (GLOBAL_TAOF){ taofdump(GLOBAL_TAOF,NULL,3);taofdump(GLOBAL_TAOF,NULL,4);}}}}}
    tfree(timera);timera=NULL;
    break;
  case 11:
    suite_11_lmitri_flag=0;
    time_0=1;
    ntimes=16; /* fix for seidcorr */
    time_delay=512;
    time_dumpevery=1024;
    timera = (double *) tcalloc(ntimes,sizeof(double));
    for (nt=0;nt<ntimes;nt++){ timera[nt] = time_0 + nt*time_dumpevery;}
    if (verbose>2){ printf(" times are: ");raprintf(timera,"double",1,ntimes,"timera ");}
    if (t==0){ /* clear everything */
      if (verbose){ printf(" %% \n");}
      if (verbose){ printf(" %% entering first phase of suite eleven -- background\n");}
      GLOBAL_TF=0; STIMULUS_ONSET_TIME=0; if (GRATING_VS_LMI!=-10){ GRATING_VS_LMI = -10;}
      if (CORTEX_BOTHER){ cortex_data_tfree();}}
    else if (t>=time_0){ /* testing seidcorr */
      nowtime = 0; for (nt=0;nt<ntimes-1;nt++){ if (t>=timera[nt] && t<timera[nt+1]){ nowtime=nt;}} 
      if (t>=timera[ntimes-1]){ nowtime=ntimes-1;}
      if (nowtime<ntimes-1){
	if (fabs(STIMULUS_ONSET_TIME-timera[nowtime]-time_delay)>0.5){
	  if (verbose){ printf(" %% entering %d_th phase of suite eleven... ",nowtime);}
	  GLOBAL_RECORD_NUMBER=nowtime;
	  if (verbose){ printf("setting GLOBAL_RECORD_NUMBER to %d\n",GLOBAL_RECORD_NUMBER);}
	  STIMULUS_ONSET_TIME = timera[nowtime]+time_delay; OUTPUT_DUMP_EVERY = 0; /* don't dump data */
	  if (verbose){ printf("setting STIMULUS_ONSET_TIME to %f\n",STIMULUS_ONSET_TIME);}
	  if (CORTEX_BOTHER){
	    if (SEIDCORR_BOTHER || GLOBAL_SEIDCORR!=NULL){ seidcorrtfree(GLOBAL_SEIDCORR); GLOBAL_SEIDCORR=NULL;}
	    if (verbose){ printf(" %% %% making seidcorr to measure");}
	    SEIDCORR_BOTHER=1;
	    GLOBAL_SEIDCORR = seidcorrmake(GLOBAL_Nra,GLOBAL_SEIDCORR_SPACE_BIN_SIZE,GLOBAL_SEIDCORR_TIME_BIN_SIZE,GLOBAL_SEIDCORR_LENGTH,STIMULUS_ONSET_TIME);
	    if (verbose){ printf(" time %f to %f\n",GLOBAL_SEIDCORR->time_start,GLOBAL_SEIDCORR->time_end);}
	    if (suite_11_lmitri_flag){
	      if (LMITRI_BOTHER==0){ 
		if (verbose){ printf(" %% %% making lmitri\n");}
		LMITRI_BOTHER=1;
		GLOBAL_LMITRI = lmitrimake(timera[nowtime],GLOBAL_LMITRI_TIMELENGTH,GLOBAL_LMITRI_ROW_MAX,GLOBAL_LMITRI_ROW_MIN);}
	      else /* if (LMITRI_BOTHER!=0) */{
		if (GLOBAL_LMITRI->time_start!=timera[nowtime]){ GLOBAL_LMITRI->time_start = timera[nowtime];}}
	      if (verbose){ printf(" %% setting GLOBAL_LMITRI->time_start %d\n",GLOBAL_LMITRI->time_start);}}}}
	if (CORTEX_BOTHER){ if (OUTPUT_DUMP_EVERY>0 && (int)floor((GLOBAL_time-timera[nowtime])/OUTPUT_DUMP_EVERY)<(int)floor((GLOBAL_time-timera[nowtime]+GLOBAL_DT)/OUTPUT_DUMP_EVERY)){
	  if (GLOBAL_LMITRI){ lmitridump(GLOBAL_LMITRI,1);}}}}
      else /* if (nowtime==ntimes-1) */{
	if (STIMULUS_ONSET_TIME != timera[nowtime]){
	  if (verbose){ printf(" %% entering %d_th phase of suite eleven -- finishing\n",nowtime);}
	  GRATING_VS_LMI=-10;
	  STIMULUS_ONSET_TIME=timera[nowtime];
	  GLOBAL_TF = timera[nowtime]+1;
	  if (CORTEX_BOTHER){ 
	    if (GLOBAL_LMITRI){ lmitridump(GLOBAL_LMITRI,1);} 
	    seidcorr_compile_helper(0,-1,0,1);
	    cortex_data_tfree();}}}}
    tfree(timera);timera=NULL;
    break;
  default: break;}
}
