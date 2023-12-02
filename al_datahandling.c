/* generically, dump_type satisfies:
   0: ascii
   1: colorscaled pnm or jpeg
   2: rawbits data dump
*/

/* These are the odor functions */

struct odor * odormake(int ntypes,int *lengthra,int base,char **strra,struct odor *o2)
{
  /* assumes that strra[nr] is of form %d%d...%d (with basedigits*lengthra[nr] components) */
  /* this sort of strra can be obtained via odorfprintf */
  /* this can also be used with strra==NULL, in which case o1 is a copy of o2 */
  int verbose=0;
  int nr1=0,nr2=0,nb=0;
  char *temp=NULL;
  struct odor *o1=NULL;
  if (verbose){ printf(" %% [entering odormake] ntypes %d base %d\n",ntypes,base); raprintf(lengthra,"int",1,ntypes,"%% lengthra:"); if (strra!=NULL){ for (nr1=0;nr1<ntypes;nr1++){ printf("%s\n",strra[nr1]);}}}
  o1 = (struct odor *) tcalloc(1,sizeof(struct odor));
  if (strra!=NULL){
    o1->ntypes = ntypes;
    o1->lengthra = (int *) tcalloc(o1->ntypes,sizeof(int));
    o1->base=base; o1->basedigits=maximum(1,(int)ceil(log10(maximum(1,o1->base))));
    temp = (char *) tcalloc(o1->basedigits,sizeof(char));
    for (nr1=0;nr1<o1->ntypes;nr1++){ o1->lengthra[nr1]=lengthra[nr1];}
    if (verbose){ 
      raprintf(o1->lengthra,"int",1,ntypes,"%% o1->lengthra:"); 
      printf("%% o1->base %d o1->basedigits %d\n",o1->base,o1->basedigits);}
    o1->rara = (double **) tcalloc(o1->ntypes,sizeof(double *)); 
    for (nr1=0;nr1<o1->ntypes;nr1++){ 
      o1->rara[nr1] = (double *) tcalloc(o1->lengthra[nr1],sizeof(double));
      for (nr2=0;nr2<minimum(strlen(strra[nr1])/o1->basedigits,o1->lengthra[nr1]);nr2++){ 
	for (nb=0;nb<o1->basedigits;nb++){ temp[nb]=strra[nr1][nb+nr2*o1->basedigits];}
	o1->rara[nr1][nr2] = (double)atoi(temp);}}
    tfree(temp);temp=NULL;}
  else if (strra==NULL && o2!=NULL){
    o1->ntypes = o2->ntypes;
    o1->lengthra = (int *) tcalloc(o1->ntypes,sizeof(int));
    for (nr1=0;nr1<o1->ntypes;nr1++){ o1->lengthra[nr1]=o2->lengthra[nr1];}
    o1->base=o2->base; o1->basedigits=o2->basedigits;
    o1->rara = (double **) tcalloc(o1->ntypes,sizeof(double *)); 
    for (nr1=0;nr1<o1->ntypes;nr1++){ 
      o1->rara[nr1] = (double *) tcalloc(o1->lengthra[nr1],sizeof(double));
      for (nr2=0;nr2<o1->lengthra[nr1];nr2++){ 
	o1->rara[nr1][nr2] = o2->rara[nr1][nr2];}}}
  else /* if (strra==NULL && o2==NULL) */{ printf(" %% Warning! o1 not defined in odormake\n");}
  return o1;
}

void odortfree(struct odor *o)
{
  int nr1=0;
  for (nr1=0;nr1<o->ntypes;nr1++){ tfree(o->rara[nr1]);o->rara[nr1]=NULL;} tfree(o->rara);o->rara=NULL;tfree(o->lengthra);o->lengthra=NULL;
  tfree(o);o=NULL;
}

void odorplusequals(struct odor *o1,double ntimes,struct odor *o2)
{
  int nr1=0,nr2=0;
  if (o1->ntypes!=o2->ntypes){ printf(" %% warning, incompatible ntypes %d, %d in odorplusequals\n",o1->ntypes,o2->ntypes);}
  for (nr1=0;nr1<o1->ntypes;nr1++){ 
    if (o1->lengthra[nr1]!=o2->lengthra[nr1]){ printf(" %% warning, o1->lra[%d]=%d,o2->lra[%d]=%d in odorplusequals\n",nr1,o1->lengthra[nr1],nr1,o2->lengthra[nr1]);}
    for (nr2=0;nr2<o1->lengthra[nr1];nr2++){ 
      o1->rara[nr1][nr2] += (ntimes*o2->rara[nr1][nr2]);}}
}

void odorfprintf(FILE *fp,struct odor *o)
{
  int nr1=0,nr2=0,nb=0; 
  for (nr1=0;nr1<o->ntypes;nr1++){ 
    for (nr2=0;nr2<o->lengthra[nr1];nr2++){ 
      for (nb=0;nb<o->basedigits-maximum(1,(int)ceil(log10(maximum(1,o->rara[nr1][nr2]))));nb++){ fprintf(fp,"0");}
      fprintf(fp,"%d",(int)o->rara[nr1][nr2]);} 
    fprintf(fp,"%s",nr1<o->ntypes-1 ? ",":";");}
  fprintf(fp,"\n");
}

void odorfprintf_full(FILE *fp,struct odor *o)
{
  int nr1=0,nr2=0; 
  for (nr1=0;nr1<o->ntypes;nr1++){ 
    for (nr2=0;nr2<o->lengthra[nr1];nr2++){ fprintf(fp," %0.16f",(double)o->rara[nr1][nr2]);} 
    fprintf(fp,"%s",nr1<o->ntypes-1 ? ",":";");}
  fprintf(fp,"\n");
}

void odor_full_fwrite(char *filename,struct odor *o)
{
  int verbose=0;
  int nt=0,nr=0;
  FILE *fp=NULL;
  if ((fp=fopen(filename,"w"))==NULL){ printf(" %% Warning, couldn't open %s in odor_full_fwrite, writing to stdout\n",filename); fp=stdout;}
  fprintf(fp,"ntypes= %d;\n",o->ntypes); if (verbose){ printf(" %% ntypes set\n");}
  fprintf(fp,"lengthra= "); for (nt=0;nt<o->ntypes;nt++){ fprintf(fp,"%d%s",o->lengthra[nt],(nt<o->ntypes-1)?",":";");} fprintf(fp,"\n");
  if (verbose){ printf(" %% lengthra set\n");}
  fprintf(fp,"base= %d;\n",o->base); if (verbose){ printf(" %% base set\n");}
  fprintf(fp,"basedigits= %d;\n",o->basedigits); if (verbose){ printf(" %% basedigits set\n");}
  for (nt=0;nt<o->ntypes;nt++){ 
    fprintf(fp,"rara[%d]= ",nt);
    for (nr=0;nr<o->lengthra[nt];nr++){
      fprintf(fp,"%0.16f%s",o->rara[nt][nr],(nr<o->lengthra[nt]-1)?",":";");}
    fprintf(fp,"\n");
    if (verbose){ printf(" %% rara[%d] set\n",nt);}}
  if (fp!=stdout){ fclose(fp); fp=NULL;}
}

struct odor * odor_full_fread(char *filename)
{
  /* here we assume that filename is a file dumped with odor_full_fwrite */
  int verbose=0;
  int nt=0,nr=0,length=0;
  FILE *fp=NULL;
  char space[128],vname[128];
  struct odor *o=NULL;
  if (verbose){ printf(" %% [entering odor_full_fread]\n");}
  if ((fp=fopen(filename,"r"))==NULL){ printf(" %% Warning, couldn't open %s in odor_full_fread\n",filename);}
  o = (struct odor *) tcalloc(1,sizeof(struct odor));
  fscanf(fp,"%[^=]",vname); fscanf(fp,"%c",space); fscanf(fp,"%c",space);
  fscanf(fp,"%d",&(o->ntypes)); if (verbose){ printf(" %% o->%s read to be %d\n",vname,o->ntypes);}
  fscanf(fp,"%c",space); fscanf(fp,"%c",space);
  fscanf(fp,"%[^=]",vname);fscanf(fp,"%c",space);fscanf(fp,"%c",space); 
  length=o->ntypes; if (o->lengthra!=NULL){ tfree(o->lengthra); o->lengthra=NULL;} o->lengthra = (int *) tcalloc(length,sizeof(int)); for (nr=0;nr<length;nr++){ fscanf(fp,"%d",&(o->lengthra[nr])); if (nr<length-1){ fscanf(fp,"%c",space);} if (verbose){ printf(" %% o->%s[%d] read to be %d\n",vname,nr,o->lengthra[nr]);}}
  fscanf(fp,"%c",space); fscanf(fp,"%c",space);
  fscanf(fp,"%[^=]",vname);fscanf(fp,"%c",space);fscanf(fp,"%c",space);
  fscanf(fp,"%d",&(o->base)); if (verbose){ printf(" %% o->%s read to be %d\n",vname,o->base);}
  fscanf(fp,"%c",space); fscanf(fp,"%c",space);
  fscanf(fp,"%[^=]",vname);fscanf(fp,"%c",space);fscanf(fp,"%c",space);
  fscanf(fp,"%d",&(o->basedigits)); if (verbose){ printf(" %% o->%s read to be %d\n",vname,o->basedigits);}
  fscanf(fp,"%c",space); fscanf(fp,"%c",space);
  o->rara = (double **) tcalloc(o->ntypes,sizeof(double *));
  for (nt=0;nt<o->ntypes;nt++){
    fscanf(fp,"%[^=]",vname);fscanf(fp,"%c",space);fscanf(fp,"%c",space);
    length=o->lengthra[nt]; if (o->rara[nt]!=NULL){ tfree(o->rara[nt]); o->rara[nt]=NULL;} o->rara[nt] = (double *) tcalloc(length,sizeof(double)); for (nr=0;nr<length;nr++){ fscanf(fp,"%lf",&(o->rara[nt][nr])); if (nr<length-1){ fscanf(fp,"%c",space);} if (verbose){ printf(" %% o->%s[%d] read to be %0.16f\n",vname,nr,o->rara[nt][nr]);}}
    fscanf(fp,"%c",space); fscanf(fp,"%c",space);}
  if (verbose){ odorfprintf_full(stdout,o);}
  if (fp!=stdout){ fclose(fp); fp=NULL;}
  return o;
}

void granule_local_findint(char *filename,int *input1,int *input2)
{
  /* returns the first nonextant _i?j?_i?j?_ within *filename. does NOT correct for changes in input1max and input2max */
  int verbose=0;
  int length = strlen(filename);
  int nl=0,nl1=0,nl2=0,nltemp=0;
  char tempchar[1],fgvn[512],filename2[1024];
  int tempint1=0,tempint2=0,tempint3=0,tempint4=0;
  FILE *fp=NULL;
  if (verbose){ printf(" %% [entering granule_local_findint] filename %s of length %d\n",filename,length);}
  sprintf(tempchar," ");
  while (nl<length-1 && !(filename[nl]=='_' && filename[nl+1]=='i')){ nl+=1;}
  if (nl<length-1 && filename[nl]=='_'){
    if (verbose){ printf(" %% found ""_"" in filename %s at position %d\n",filename,nl);}
    nl+=1; *tempchar = filename[nl]; if (verbose){ printf(" %% assuming that %s is an ""i""\n",tempchar);} 
    nl+=1; *tempchar = filename[nl]; nl1=0;
    while (nl+nl1<length && *tempchar!='j'){ nl1+=1; *tempchar = filename[nl+nl1];}
    if (nl+nl1<length && *tempchar=='j'){
      if (verbose){ printf(" %% found nl1 %d in filename %s at position %d+%d\n",nl1,filename,nl,nl1);}
      tempint1=0; 
      for (nltemp=0;nltemp<nl1;nltemp++){ *tempchar = filename[nl+nltemp]; tempint1 += (int)pow(10,nl1-nltemp-1)*atoi(tempchar);}
      if (verbose){ printf(" %% read tempint1 %d\n",tempint1);}
      *tempchar = filename[nl+nl1];
      if (verbose){ printf(" %% assuming that %s is a ""j""\n",tempchar);}
      nl+=1; *tempchar = filename[nl+nl1]; nl2=0;
      while (nl+nl1+nl2<length && *tempchar!='_'){ nl2+=1; *tempchar = filename[nl+nl1+nl2];}
      if (nl+nl1+nl2<length && *tempchar=='_'){
	if (verbose){ printf(" %% found nl2 %d in filename %s at position %d+%d+%d\n",nl2,filename,nl,nl1,nl2);}
	tempint2=0; 
	for (nltemp=0;nltemp<nl2;nltemp++){ *tempchar = filename[nl+nl1+nltemp]; tempint2 += (int)pow(10,nl2-nltemp-1)*atoi(tempchar);}
	if (verbose){ printf(" %% read tempint2 %d\n",tempint2);}}
      else /**/{ printf(" %% warning! aborted at second step nl2 %d in granule_local_findint\n",nl2);}
      if (nl+nl1+nl2<length-1 && filename[nl+nl1+nl2]=='_'){
	nl=nl+nl1+nl2;
	if (verbose){ printf(" %% found ""_"" in filename %s at position %d\n",filename,nl);}
	nl+=1; *tempchar = filename[nl]; if (verbose){ printf(" %% assuming that %s is an ""i""\n",tempchar);} 
	nl+=1; *tempchar = filename[nl]; nl1=0;
	while (nl+nl1<length && *tempchar!='j'){ nl1+=1; *tempchar = filename[nl+nl1];}
	if (nl+nl1<length && *tempchar=='j'){
	  if (verbose){ printf(" %% found nl1 %d in filename %s at position %d+%d\n",nl1,filename,nl,nl1);}
	  tempint3=0; 
	  for (nltemp=0;nltemp<nl1;nltemp++){ *tempchar = filename[nl+nltemp]; tempint3 += (int)pow(10,nl1-nltemp-1)*atoi(tempchar);}
	  if (verbose){ printf(" %% read tempint3 %d\n",tempint3);}
	  *tempchar = filename[nl+nl1];
	  if (verbose){ printf(" %% assuming that %s is a ""j""\n",tempchar);}
	  nl+=1; *tempchar = filename[nl+nl1]; nl2=0;
	  while (nl+nl1+nl2<length && *tempchar!='_'){ nl2+=1; *tempchar = filename[nl+nl1+nl2];}
	  if (nl+nl1+nl2<length && *tempchar=='_'){
	    if (verbose){ printf(" %% found nl2 %d in filename %s at position %d+%d+%d\n",nl2,filename,nl,nl1,nl2);}
	    tempint4=0; 
	    for (nltemp=0;nltemp<nl2;nltemp++){ *tempchar = filename[nl+nl1+nltemp]; tempint4 += (int)pow(10,nl2-nltemp-1)*atoi(tempchar);}
	    if (verbose){ printf(" %% read tempint4 %d\n",tempint4);}}
	  else /**/{ printf(" %% warning! aborted at third step nl2 %d in granule_local_findint\n",nl2);}}}}}
  else /**/{ printf(" %% warning! aborted at first step nl %d in granule_local_findint\n",nl);}
  sprintf(fgvn,"i%dj%d",tempint1,tempint2);
  sprintf(filename2,"./ptree_%s_%srecord",fgvn,GLOBAL_STRING_2);
  if (!checktofind(filename2)){ 
    if (input1!=NULL){ *input1=tempint1;}
    if (input2!=NULL){ *input2=tempint2;}
    if (verbose){ printf(" %% input (%d,%d) not yet made in granule_local_findint\n",tempint1,tempint2);}
    if ((fp=fopen(filename2,"w"))==NULL){ printf(" %% Warning, couldn't open %s in granule_findint, writing to stdout\n",filename2); fp=stdout;} if (fp!=stdout){ fclose(fp);}}
  else /* if found */{
    if (verbose>1){ printf(" %% input (%d,%d) already made in granule_findint\n",tempint1,tempint2);}
    sprintf(fgvn,"i%dj%d",tempint3,tempint4);
    sprintf(filename2,"./ptree_%s_%srecord",fgvn,GLOBAL_STRING_2);
    if (!checktofind(filename2)){ 
      if (input1!=NULL){ *input1=tempint3;}
      if (input2!=NULL){ *input2=tempint4;}
      if (verbose){ printf(" %% input (%d,%d) not yet made in granule_local_findint\n",tempint3,tempint4);}
      if ((fp=fopen(filename2,"w"))==NULL){ printf(" %% Warning, couldn't open %s in granule_findint, writing to stdout\n",filename2); fp=stdout;} if (fp!=stdout){ fclose(fp);}}
    else /* if found */{
      if (verbose>1){ printf(" %% input (%d,%d) already made in granule_findint\n",tempint3,tempint4);}
      printf(" %% Warning, no more work to do in granule_local_findint, exiting\n"); exit(EXIT_SUCCESS);}}
}

void granule_findint(int input1max,int input2max,int *input1,int *input2)
{
  /* searches array */
  int verbose=0;
  char fgvn[512],filename[1024];
  int nr1=0,nr2=0,notfound_flag=0;
  int input2max_scale = (int)(input2max*sqrt(3)/1.5);
  FILE *fp=NULL;
  for (nr1=0;nr1<input1max;nr1++){ for (nr2=0;nr2<input2max_scale;nr2++){
    sprintf(fgvn,"i%dj%d",nr1,nr2);
    sprintf(filename,"./ptree_%s_%srecord",fgvn,GLOBAL_STRING_2);
    if (!checktofind(filename)){ 
      notfound_flag=1;
      if (input1!=NULL){ *input1=nr1;}
      if (input2!=NULL){ *input2=nr2;}
      if (verbose){ printf(" %% input (%d,%d) not yet made in granule_findint\n",nr1,nr2);}
      if ((fp=fopen(filename,"w"))==NULL){ printf(" %% Warning, couldn't open %s in granule_findint, writing to stdout\n",filename); fp=stdout;}
      if (fp!=stdout){ fclose(fp);}
      nr1=input1max+1;nr2=input2max_scale+1; /* break */}
    else /* if (checktofind(filename)) */{
      if (verbose>1){ printf(" %% input (%d,%d) already made in granule_findint\n",nr1,nr2);}}}}
  if (!notfound_flag){ printf(" %% Warning, no more work to do in granule_findint, exiting\n"); exit(EXIT_SUCCESS);}
}

void granule_int2input(struct odor *o1,int input1max,int input1,int input2max,int input2)
{
  int verbose=0;
  double tempx = (double)(input2*1.5/sqrt(3))/(double)input2max;
  double tempy = (double)(input1+0.5-(double)(input2%2)/2.0)/(double)input1max;
  int nr1=0,nr2=0;
  for (nr1=0;nr1<o1->ntypes;nr1++){ for (nr2=0;nr2<o1->lengthra[nr1];nr2++){
    o1->rara[nr1][nr2] = o1->base*(nr2%2==0 ? tempx : tempy);}}
  if (verbose){ printf("odor set to ");odorfprintf_full(stdout,o1);}
  if (verbose){ printf("odor will be written as ");odorfprintf(stdout,o1);}
  if (verbose){ granule_input2int(o1,input1max,&nr1,input2max,&nr2); printf("odor will be recorded as i%dj%d\n",nr1,nr2);}
}

void granule_input2int(struct odor *o1,int input1max,int *input1,int input2max,int *input2)
{
  int nr1=0;
  if (input2!=NULL){ *input2 = (int)(o1->rara[0][0]/(double)o1->base*input2max*sqrt(3)/1.5+0.25);}
  if (input1!=NULL){ 
    nr1=0; while (nr1<o1->ntypes && o1->lengthra[nr1]<2){ nr1+=1;} 
    if (nr1<o1->ntypes && o1->lengthra[nr1]>=2){ 
      *input1 = (int)(o1->rara[nr1][1]/(double)o1->base*input1max - 0.5 + (double)(input2!=NULL?(*input2)%2:0)/2.0 + 0.25);}
    else if (nr1>=o1->ntypes || o1->lengthra[nr1]<2){ printf(" %% warning! input2 not set in granule_input2int\n");}}
}

void granule_int2input_or_input2int(int direction_flag,struct odor *o1,int input1max,int *input1,int input2max,int *input2)
{
  int verbose=0,fig_verbose=1;
  int remove_flag=0,jpg_flag=0,eps_flag=(ON_MY_COMPUTER ? 1 : 0);
  char tempchar[512];
  int iorigin=0,jorigin=0;
  double inputx=0,inputy=0,xorigin=0,yorigin=0,rad_input_max=0;
  int input2max_scale = (int)(input2max*sqrt(3)/1.5);
  double ***Ira=NULL,*Inorm=NULL,*Icomp=NULL,*Ipair=NULL;
  int *Imaxindex=NULL,Inmax=0,*Iminindex=NULL,Inmin=0;
  int nv=0,nvec=0,nt=0,ntypes=0,nr=0,*lengthra=0,input_k=0,input_kp1=0;
  double theta_input=0,theta_k=0,theta_kp1=0,theta_internal=0,temp=0,temp2=0,rad_input=0,rad_min=0,rad_max=0;
  double tolerance=0.000001;
  int nr1=0,nr2=0,nr3=0,nr4=0,ns=0;
  char filename1a[1024],filename1b[1024];
  double maxdia = 20000;
  FILE *fp=NULL;
  char command[1024];
  double xord1=0,yord1=0,rcolor=0,gcolor=0,bcolor=0;
  int colorcode=0;
  if (verbose){ 
    printf(" %% [entering granule_int2input_or_input2int] direction_flag %d\n",direction_flag); 
    if (verbose>1){ printf(" %% odor given as "); odorfprintf_full(stdout,o1);}}
  iorigin = (double)input1max/2.0; jorigin = (double)input2max_scale/2.0;
  xorigin = (sqrt(3)/2.0*jorigin); yorigin = iorigin - 0.5*(jorigin%2);
  rad_input_max = minimum(yorigin,xorigin);
  if (verbose>1){ printf(" %% origin (%d,%d)-->(%0.2f,%0.2f), rad_input_max %0.2f\n",iorigin,jorigin,xorigin,yorigin,rad_input_max);}
  nvec = 3; ntypes = o1->ntypes; lengthra = o1->lengthra; rad_min=0.1; rad_max=1.5;
  Ira = (double ***) tcalloc(nvec,sizeof(double **));
  for (nv=0;nv<nvec;nv++){
    Ira[nv] = (double **) tcalloc(ntypes,sizeof(double *));
    for (nt=0;nt<ntypes;nt++){ 
      Ira[nv][nt] = (double *) tcalloc(lengthra[nt],sizeof(double));
      for (nr=0;nr<lengthra[nt];nr++){ Ira[nv][nt][nr] = (nr%nvec == nv);}}}
  Inorm = (double *) tcalloc(nvec,sizeof(double));
  for (nv=0;nv<nvec;nv++){ for (nt=0;nt<ntypes;nt++){ Inorm[nv] += ra2ra_dot(Ira[nv][nt],Ira[nv][nt],lengthra[nt]);}}
  for (nv=0;nv<nvec;nv++){ for (nt=0;nt<ntypes;nt++){ ratimesequals(Ira[nv][nt],lengthra[nt],1.0/(double)sqrt(Inorm[nv]));}}
  if (verbose>1){
    for (nv=0;nv<nvec;nv++){ for (nt=0;nt<ntypes;nt++){ 
      sprintf(tempchar," %% Ira nv %d nt %d",nv,nt);
      raprintf(Ira[nv][nt],"double",1,lengthra[nt],tempchar);}}}
  if (direction_flag==1){
    inputx = sqrt(3)/2.0*(*input2); inputy = (*input1) - 0.5*((*input2)%2);
    if (verbose){ printf(" %% input given as i%dj%d (x=%0.2f,y=%0.2f)\n",*input1,*input2,inputx,inputy);}
    theta_input = periodize(atan2(inputy - yorigin,inputx - xorigin),0,2*PI);
    rad_input = sqrt(pow(inputx - xorigin,2)+pow(inputy - yorigin,2));
    input_k = periodize((int)(theta_input*nvec/(2*PI)),0,nvec);
    input_kp1 = periodize(input_k+1,0,nvec);
    theta_k = 2*PI/(double)nvec*input_k; theta_kp1 = 2*PI/(double)nvec*(input_k+1);
    theta_internal = (theta_input-theta_k)/(theta_kp1-theta_k)*PI/2;
    if (verbose){ printf(" %% theta_input %0.2f rad_input %0.2f, input_k %d theta_k %0.2f, theta_kp1 %0.2f, theta_internal %0.2f\n",theta_input,rad_input,input_k,theta_k,theta_kp1,theta_internal);}
    for (nt=0;nt<ntypes;nt++){ for (nr=0;nr<lengthra[nt];nr++){
      temp = cos(theta_internal)*Ira[input_k][nt][nr] + sin(theta_internal)*Ira[input_kp1][nt][nr];
      temp2 = rad_min + rad_input/rad_input_max*(rad_max-rad_min);
      o1->rara[nt][nr] = o1->base*temp*temp2;}}
    if (verbose){ printf(" %% odor set to ");odorfprintf_full(stdout,o1);}}
  else if (direction_flag==0){
    if (verbose){ printf(" %% odor given as ");odorfprintf_full(stdout,o1);}
    Icomp = (double *) tcalloc(nvec,sizeof(double));
    for (nv=0;nv<nvec;nv++){ for (nt=0;nt<ntypes;nt++){ 
      Icomp[nv] += ra2ra_dot(Ira[nv][nt],o1->rara[nt],lengthra[nt])/(double)o1->base;}}
    Ipair = (double *) tcalloc(nvec,sizeof(double));
    for (nv=0;nv<nvec;nv++){
      input_k = nv; input_kp1 = periodize(input_k+1,0,nvec);
      Ipair[nv] = sqrt(pow(Icomp[input_k],2)+pow(Icomp[input_kp1],2));}
    if (verbose>1){ raprintf(Ipair,"double",1,nvec,"%% Ipair");}
    maxmindex("double",Ipair,nvec,&Imaxindex,&Inmax,&Iminindex,&Inmin,tolerance);
    input_k = Imaxindex[0]; input_kp1 = periodize(input_k+1,0,nvec);
    tfree(Imaxindex);Imaxindex=NULL;
    tfree(Iminindex);Iminindex=NULL;
    if (fabs(Icomp[input_k]-0)>tolerance){ theta_internal = atan(Icomp[input_kp1]/Icomp[input_k]);}
    else /* if (fabs(Icomp[input_k]-0)<=tolerance) */{ theta_internal = PI/2;}
    theta_k = 2*PI/(double)nvec*input_k; theta_kp1 = 2*PI/(double)nvec*(input_k+1);
    theta_input = theta_internal*(2/PI)*(theta_kp1-theta_k) + theta_k;
    rad_input = (sqrt(pow(Icomp[input_kp1],2)+pow(Icomp[input_k],2))-rad_min)/(rad_max-rad_min)*rad_input_max;
    inputx = rad_input*cos(theta_input)+xorigin; inputy = rad_input*sin(theta_input)+yorigin;
    if (verbose){ printf(" %% theta_input %0.2f rad_input %0.2f, input_k %d theta_k %0.2f, theta_kp1 %0.2f, theta_internal %0.2f\n",theta_input,rad_input,input_k,theta_k,theta_kp1,theta_internal);}
    if (input2!=NULL){ *input2 = adi_round(2.0/sqrt(3)*inputx); if (input1!=NULL){ *input1 = adi_round(inputy + 0.5*((*input2)%2));}}
    if (verbose){ printf(" %% input set to i%dj%d (x=%0.2f,y=%0.2f)\n",*input1,*input2,inputx,inputy);}}
  else /* if (direction_flag!=0 && direction_flag!=1) */{
    if (verbose){ printf(" %% testing granule_int2input_or_input2int\n");}
    sprintf(filename1a,"./granule_int2input_or_input2int_test");
    sprintf(filename1b,"%s.fig",filename1a);
    if ((fp=fopen(filename1b,"w"))==NULL){ printf(" %% Can't create %s in granule_int2input_or_input2int\n",filename1b); fp=stdout;}
    fprintf(fp,"%s",FIG_PREAMBLE); fprintf(fp,"%s",FIG_PREAMBLE_COLOR_7);
    for (nr1=0;nr1<input1max;nr1++){ for (nr2=0;nr2<input2max_scale;nr2++){
      granule_int2input_or_input2int(1,o1,input1max,&nr1,input2max,&nr2);
      Icomp = (double *) tcalloc(nvec,sizeof(double));
      for (nv=0;nv<nvec;nv++){ for (nt=0;nt<ntypes;nt++){ 
      Icomp[nv] += ra2ra_dot(Ira[nv][nt],o1->rara[nt],lengthra[nt])/(double)o1->base;}}
      xord1 = (double)nr2*1.5/sqrt(3); yord1 = nr1 - (double)(nr2%2)/2.0; temp=0;
      for (nv=0;nv<nvec;nv++){
	temp += pow(Icomp[nv],2);
	colorscale(0,Icomp[nv],rad_max,0,&rcolor,&gcolor,&bcolor);colorcode=crop((int)floor(512*rcolor),0,511);
	fprintf(fp,"2 1 0 0 0 %d %d 0 %d 0.000 0 0 -1 0 0 %d\n",/*colorcode*/colorcode+32,/*depth*/nv,/*fill*/20,/*npoints*/6+1); fprintf(fp,"\t"); for (ns=0;ns<6+1;ns++){ fprintf(fp,"%d %d ",(int)floor(maxdia*(xord1+cos(2*PI*ns/6.0)/sqrt(3))/(double)input2max),(int)floor(maxdia*(yord1+sin(2*PI*ns/6.0)/sqrt(3))/(double)input1max));} fprintf(fp,"\n");}
      colorscale(0,sqrt(temp),rad_max,0,&rcolor,&gcolor,&bcolor);colorcode=crop((int)floor(512*rcolor),0,511);
      fprintf(fp,"2 1 0 0 0 %d %d 0 %d 0.000 0 0 -1 0 0 %d\n",/*colorcode*/colorcode+32,/*depth*/nvec,/*fill*/20,/*npoints*/6+1); fprintf(fp,"\t"); for (ns=0;ns<6+1;ns++){ fprintf(fp,"%d %d ",(int)floor(maxdia*(xord1+cos(2*PI*ns/6.0)/sqrt(3))/(double)input2max),(int)floor(maxdia*(yord1+sin(2*PI*ns/6.0)/sqrt(3))/(double)input1max));} fprintf(fp,"\n");
      if (Icomp!=NULL){ tfree(Icomp);Icomp=NULL;}
      if (fig_verbose){ fprintf(fp,"4 0 0 %d 0 %d %d 0.0000 4 %d %d %d %d i%dj%d\\001\n",/*depth*/nvec+1,/*font*/5,/*point*/5,/*textheight*/(int)(5*13.5),/*textwidth*/(int)(5*9*5),/*xpos*/(int)floor(maxdia*xord1/(double)input2max),/*ypos*/(int)floor(maxdia*yord1/(double)input1max),nr1,nr2);}
      granule_int2input_or_input2int(0,o1,input1max,&nr3,input2max,&nr4);
      if (nr1!=nr3 || nr2!=nr4){
	printf(" warning! input (%d,%d) --> odor:\n",nr1,nr2); odorfprintf_full(stdout,o1);
	printf(" yet odor gives back input (%d,%d)\n",nr3,nr4);}}}
    if (jpg_flag){ sprintf(command,"fig2dev -Ljpeg -q %d %s.fig %s.jpg;",/*quality*/5,filename1a,filename1a); system(command);}
    if (eps_flag){ sprintf(command,"fig2dev -Leps %s.fig %s.eps;",filename1a,filename1a); system(command);}
    if (remove_flag){ sprintf(command,"rm %s.fig;",filename1a); system(command);}
    if (fp!=stdout){ fclose(fp);}}
  if (Ira!=NULL){ for(nv=0;nv<nvec;nv++){ if (Ira[nv]!=NULL){ for (nt=0;nt<ntypes;nt++){ if (Ira[nv][nt]!=NULL){ tfree(Ira[nv][nt]);Ira[nv][nt]=NULL;}} tfree(Ira[nv]); Ira[nv]=NULL;}} tfree(Ira); Ira=NULL;}
  if (Inorm!=NULL){ tfree(Inorm); Inorm=NULL;} if (Icomp!=NULL){ tfree(Icomp); Icomp=NULL;} if (Ipair!=NULL){ tfree(Ipair); Ipair=NULL;}
}

void granule_local_plot()
{
  /* fix later -- must be edited on a case by case basis */
  int verbose=0;
  int length=0,maxlevel=0,nr=0,nl=0;
  char filename1[1024],filename2[1024];
  char filename_base1[1024],filename_base2[1024];
  char tempchar[512];
  double *temp=NULL,*ra=NULL;
  struct ptree *p1=NULL,*p2=NULL;
  int nrecords = 1024*SUITE_NSECONDS/SUITE_DUMPEVERY;
  int no_error1=0,no_error2=0;
  if (verbose){ printf(" %% [entering granule_local_plot]\n");}
  sprintf(filename_base1,"ptree_i23j15_granuleC_");
  sprintf(filename_base2,"ptree_i22j16_granuleC_");
  length=19;
  maxlevel=3;
  ra = (double *) tcalloc(length*maxlevel,sizeof(double));
  for (nr=0;nr<length;nr++){
    sprintf(filename1,"./%s%d0record",filename_base1,nr);
    sprintf(filename2,"./%s%d0record",filename_base2,nr);
    if (checktofind(filename1)){
      if (verbose){ printf(" %% base file %s\n",filename1);}
      p1=ptree_readbacktemp(filename1,&file2obsdisthist,&no_error1);
      if (no_error1){
	if (checktofind(filename2)){
	  if (verbose){ printf(" %% base file %s\n",filename2);}
	  p2=ptree_readbacktemp(filename2,&file2obsdisthist,&no_error2);
	  if (no_error2){
	    temp = (double *) tcalloc(2,sizeof(double)); 
	    for (nl=0;nl<maxlevel;nl++){ 
	      temp[0]=nrecords; temp[1]=0;
	      pnode2pnode_void_operate(p1,NULL,p1->postree,p2,NULL,p2->postree,nl+1,0,&pnode_obsdisthist2obsdisthist_appraise,temp,0);
	      printf("filenames %s vs %s, yields %f/%f\n",filename1,filename2,temp[1],temp[0]);
	      ra[nr+nl*length] = maximum(nrecords,temp[1]);}
	    tfree(temp);temp=NULL;
	    pnode_obs2dist_starter(0,0,0,NULL,p2->postree,-1,0,2);
	    ptreetfree(p2);p2=NULL;}}
	pnode_obs2dist_starter(0,0,0,NULL,p1->postree,-1,0,2);
	ptreetfree(p1);p1=NULL;}}}
  for (nl=0;nl<maxlevel;nl++){
    if (verbose){ raprintf(&(ra[0+nl*length]),"double",1,length,"%% ra");}
    sprintf(tempchar,"./granule_local_plot_level%d",nl);
    ra2jpg(ra,"double",1,length,0,tempchar,0);}
  tfree(ra);ra=NULL;
}

void granule_plot(int input1max,int input2max)
{
  int verbose=0,fig_verbose=0;
  int remove_flag=0,jpg_flag=0,eps_flag=(ON_MY_COMPUTER ? 1 : 0);
  int input2max_scale = (int)(input2max*sqrt(3)/1.5);
  int nr1=0,nr2=0,nl=0,maxlevel=0,ni=0,nr1a=0,nr2a=0,nr3=0,tab=0,depthmin=45,depthmax=55,depth=0;
  char filename1a[1024],filename1b[1024],*gs2=GLOBAL_STRING_2,tempchar[16];
  double *temp=NULL,*ra=NULL,*ra2=NULL,*ra_angle=NULL,*ra_reweight=NULL;
  struct ptree *p1=NULL,*p2=NULL,*pW=NULL;
  double *bipoo_ra=NULL;
  double time_total=SUITE_DUMPEVERY; int nrecords = 1024*SUITE_NSECONDS/SUITE_DUMPEVERY;
  double maxdia = 20000;
  FILE *fp=NULL;
  char command[1024];
  double xord1=0,yord1=0,xord2=0,yord2=0,xord3=0,yord3=0,xord4=0,yord4=0,xord5=0,yord5=0;
  double rcolor=0,gcolor=0,bcolor=0;
  double max=0,min=0,mean=0,stdev=0;
  int colorcode=0;
  int no_error1=0,no_error2=0;
  if (verbose){ printf(" %% [entering granule_plot] input1max %d input2max %d\n",input1max,input2max);}
  for (ni=0;ni<2;ni++){
    switch(ni){ 
    case 0: sprintf(tempchar,"_"); maxlevel=GLOBAL_PTREE_NLEGS+1; break; 
    case 1: sprintf(tempchar,"_input_"); maxlevel=1; break; 
    default: break;}
    ra = (double *) tcalloc(3*input1max*input2max_scale*maxlevel,sizeof(double));
    ra_angle = (double *) tcalloc(3*input1max*input2max_scale,sizeof(double));
    ra_reweight = (double *) tcalloc(3*input1max*input2max_scale,sizeof(double));
    for (nr1=0;nr1<input1max;nr1++){ for (nr2=0;nr2<input2max_scale;nr2++){
      sprintf(filename1a,"./ptree%si%dj%d_%srecord",tempchar,nr1,nr2,gs2);
      if (checktofind(filename1a)){
	if (verbose){ printf(" %% base file %s\n",filename1a);}
	p1=ptree_readbacktemp(filename1a,&file2obsdisthist,&no_error1);
	if (no_error1){
	  bipoo_ra = granule_obtain_bipoo_matrix(filename1a,time_total,(double)nrecords);
	  if (nr1+1<input1max){ 
	    sprintf(filename1b,"./ptree%si%dj%d_%srecord",tempchar,nr1+1,nr2,gs2);
	    if (checktofind(filename1b)){
	      if (verbose){ printf(" %% compared to file %s\n",filename1b);}
	      p2=ptree_readbacktemp(filename1b,&file2obsdisthist,&no_error2);
	      if (no_error2){
		temp = (double *) tcalloc(2,sizeof(double)); 
		for (nl=0;nl<maxlevel;nl++){ 
		  temp[0]=nrecords; temp[1]=0;
		  pnode2pnode_void_operate(p1,NULL,p1->postree,p2,NULL,p2->postree,nl+1,0,&pnode_obsdisthist2obsdisthist_appraise,temp,0);
		  ra[0+nr1*3+nr2*3*input1max+nl*3*input1max*input2max_scale] = maximum(nrecords,temp[1]);}
		pW = granule_obtain_reweight(filename1a,filename1b,time_total,(double)nrecords); ra_angle[0+nr1*3+nr2*3*input1max] = granule_obtain_reweight_vs_bipoo_angle(pW,bipoo_ra); pnodefabswr_starter(NULL,pW->postree,0,1); pnodestats_starter(NULL,pW->postree,2,2,0,NULL,NULL,NULL,NULL,NULL,&(ra_reweight[0+nr1*3+nr2*3*input1max]),NULL,NULL,NULL); ptreetfree(pW); pW=NULL;
		tfree(temp);temp=NULL;
		pnode_obs2dist_starter(0,0,0,NULL,p2->postree,-1,0,2);
		ptreetfree(p2);p2=NULL;}
	      else /* if (!no_error2) */{ printf(" %% warning! ptree %s improperly read in granule_plot\n",filename1b);}}}
	  if (nr2+1<input2max_scale){
	    sprintf(filename1b,"./ptree%si%dj%d_%srecord",tempchar,nr1,nr2+1,gs2);
	    if (checktofind(filename1b)){
	      if (verbose){ printf(" %% compared to file %s\n",filename1b);}
	      p2=ptree_readbacktemp(filename1b,&file2obsdisthist,&no_error2);
	      if (no_error2){
		temp = (double *) tcalloc(2,sizeof(double)); 
		for (nl=0;nl<maxlevel;nl++){
		  temp[0]=nrecords; temp[1]=0;
		  pnode2pnode_void_operate(p1,NULL,p1->postree,p2,NULL,p2->postree,nl+1,0,&pnode_obsdisthist2obsdisthist_appraise,temp,0);
		  ra[1+nr1*3+nr2*3*input1max+nl*3*input1max*input2max_scale] = maximum(nrecords,temp[1]);}
		pW = granule_obtain_reweight(filename1a,filename1b,time_total,(double)nrecords); ra_angle[1+nr1*3+nr2*3*input1max] = granule_obtain_reweight_vs_bipoo_angle(pW,bipoo_ra); pnodefabswr_starter(NULL,pW->postree,0,1); pnodestats_starter(NULL,pW->postree,2,2,0,NULL,NULL,NULL,NULL,NULL,&(ra_reweight[1+nr1*3+nr2*3*input1max]),NULL,NULL,NULL); ptreetfree(pW); pW=NULL;
		tfree(temp);temp=NULL;
		pnode_obs2dist_starter(0,0,0,NULL,p2->postree,-1,0,2);
		ptreetfree(p2);p2=NULL;}
	      else /* if (!no_error2) */{ printf(" %% warning! ptree %s improperly read in granule_plot\n",filename1b);}}}
	  if (nr2+1<input2max_scale && nr1+(int)pow(-1,nr2%2)>=0 && nr1+(int)pow(-1,nr2%2)<input1max){
	    sprintf(filename1b,"./ptree%si%dj%d_%srecord",tempchar,nr1+(int)pow(-1,nr2%2),nr2+1,gs2);
	    if (checktofind(filename1b)){
	      if (verbose){ printf(" %% compared to file %s\n",filename1b);}
	      p2=ptree_readbacktemp(filename1b,&file2obsdisthist,&no_error2);
	      if (no_error2){
		temp = (double *) tcalloc(2,sizeof(double)); 
		for (nl=0;nl<maxlevel;nl++){
		  temp[0]=nrecords; temp[1]=0;
		  pnode2pnode_void_operate(p1,NULL,p1->postree,p2,NULL,p2->postree,nl+1,0,&pnode_obsdisthist2obsdisthist_appraise,temp,0);
		  ra[2+nr1*3+nr2*3*input1max+nl*3*input1max*input2max_scale] = maximum(nrecords,temp[1]);}
		pW = granule_obtain_reweight(filename1a,filename1b,time_total,(double)nrecords); ra_angle[2+nr1*3+nr2*3*input1max] = granule_obtain_reweight_vs_bipoo_angle(pW,bipoo_ra); pnodefabswr_starter(NULL,pW->postree,0,1); pnodestats_starter(NULL,pW->postree,2,2,0,NULL,NULL,NULL,NULL,NULL,&(ra_reweight[2+nr1*3+nr2*3*input1max]),NULL,NULL,NULL); ptreetfree(pW); pW=NULL;
		tfree(temp);temp=NULL;
		pnode_obs2dist_starter(0,0,0,NULL,p2->postree,-1,0,2);
		ptreetfree(p2);p2=NULL;}
	      else /* if (!no_error2) */{ printf(" %% warning! ptree %s improperly read in granule_plot\n",filename1b);}}}
	  tfree(bipoo_ra);bipoo_ra=NULL;
	  pnode_obs2dist_starter(0,0,0,NULL,p1->postree,-1,0,2);
	  ptreetfree(p1);p1=NULL;}
	else /* if (!no_error1) */{ printf(" %% warning! ptree %s improperly read in granule_plot\n",filename1a);}}}}
    stats("double",ra,3*input1max*input2max_scale*maxlevel,&max,&min,&mean,&stdev); max = mean+STD_VIEW*stdev; min = mean-STD_VIEW*stdev;
    max = nrecords*1.75;
    min = nrecords*1.0;
    for (nl=0;nl<maxlevel;nl++){
      ra2 = &(ra[0 + 0*3 + 0*3*input1max + nl*3*input1max*input2max_scale]);
      if (verbose){ raprintf(ra2,"double",3*input1max,input2max_scale,"%% ra2:");}
      sprintf(filename1a,"./ptree%sadjacencyplot_maxlevel%d_%srecord",tempchar,nl+1,gs2);
      sprintf(filename1b,"%s.fig",filename1a);
      if ((fp=fopen(filename1b,"w"))==NULL){ printf(" %% Warning, cannot create %s in granule_plot\n",filename1b); fp=stdout;}
      fprintf(fp,"%s",FIG_PREAMBLE); fprintf(fp,"%s",FIG_PREAMBLE_COLOR_7);
      for (nr1=0;nr1<input1max;nr1++){ for (nr2=0;nr2<input2max_scale;nr2++){
	xord1 = (double)nr2*1.5/sqrt(3); yord1 = nr1 - (double)(nr2%2)/2.0;
	for (nr3=0;nr3<3;nr3++){
	  switch(nr3){
	  case 0: nr1a = nr1+1; nr2a = nr2+0; break;
	  case 1: nr1a = nr1+0; nr2a = nr2+1; break;
	  case 2: nr1a = nr1+(int)pow(-1,nr2%2); nr2a = nr2+1; break;
	  default: break;}
	  if (nr1a>=0 && nr1a<input1max && nr2a>=0 && nr2a<input2max_scale){
	    xord2 = (double)nr2a*1.5/sqrt(3); yord2 = nr1a - (double)(nr2a%2)/2.0;
	    xord3 = (xord1+xord2)/2.0 + sqrt(3)/6.0*(yord2-yord1); yord3 = (yord1+yord2)/2.0 + sqrt(3)/6.0*(xord1-xord2);
	    xord4 = (xord1+xord2)/2.0 - sqrt(3)/6.0*(yord2-yord1); yord4 = (yord1+yord2)/2.0 - sqrt(3)/6.0*(xord1-xord2);
	    xord5 = (xord1+xord2+xord3+xord4)/4.0; yord5 = (yord1+yord2+yord3+yord4)/4.0;
	    tab = nr3+nr1*3+nr2*3*input1max; periodify("int",&tab,&depthmin,&depthmax,&depth);
	    colorscale(0,ra2[tab],max,min,&rcolor,&gcolor,&bcolor);colorcode=crop((int)floor(512*rcolor),0,511);
	    fprintf(fp,"2 1 0 0 0 %d %d 0 %d 0.000 0 0 -1 0 0 %d\n",/*colorcode*/colorcode+32,/*depth*/depth,/*fill*/20,/*npoints*/4+1); fprintf(fp,"\t"); fprintf(fp,"%d %d ",(int)floor(maxdia*xord1/(double)input2max),(int)floor(maxdia*yord1/(double)input1max)); fprintf(fp,"%d %d ",(int)floor(maxdia*xord3/(double)input2max),(int)floor(maxdia*yord3/(double)input1max)); fprintf(fp,"%d %d ",(int)floor(maxdia*xord2/(double)input2max),(int)floor(maxdia*yord2/(double)input1max)); fprintf(fp,"%d %d ",(int)floor(maxdia*xord4/(double)input2max),(int)floor(maxdia*yord4/(double)input1max)); fprintf(fp,"%d %d ",(int)floor(maxdia*xord1/(double)input2max),(int)floor(maxdia*yord1/(double)input1max)); fprintf(fp,"\n"); if (fig_verbose){ fprintf(fp,"4 0 0 %d 0 %d %d 0.0000 4 %d %d %d %d %d%d_%d%d\\001\n",/*depth*/depthmax+1,/*font*/5,/*point*/5,/*textheight*/(int)(5*13.5),/*textwidth*/(int)(5*9*5),/*xpos*/(int)floor(maxdia*xord5/(double)input2max),/*ypos*/(int)floor(maxdia*yord5/(double)input1max),nr1,nr2,nr1a,nr2a);}}}}}
      if (jpg_flag){ sprintf(command,"fig2dev -Ljpeg -q %d %s.fig %s.jpg;",/*quality*/5,filename1a,filename1a); system(command);}
      if (eps_flag){ sprintf(command,"fig2dev -Leps %s.fig %s.eps;",filename1a,filename1a); system(command);}
      if (remove_flag){ sprintf(command,"rm %s.fig;",filename1a); system(command);}
      if (fp!=stdout){ fclose(fp);}}
    stats("double",ra_reweight,3*input1max*input2max_scale,&max,&min,&mean,&stdev); max = mean+STD_VIEW*stdev; min = mean-STD_VIEW*stdev;
    /*     max = 1.0; */
    /*     min = 0.0; */
    sprintf(filename1a,"./ptree%sreweight_adjacencyplot_%srecord",tempchar,gs2);
    sprintf(filename1b,"%s.fig",filename1a);
    if ((fp=fopen(filename1b,"w"))==NULL){ printf(" %% Warning, cannot create %s in granule_plot\n",filename1b); fp=stdout;}
    fprintf(fp,"%s",FIG_PREAMBLE); fprintf(fp,"%s",FIG_PREAMBLE_COLOR_7);
    for (nr1=0;nr1<input1max;nr1++){ for (nr2=0;nr2<input2max_scale;nr2++){
      xord1 = (double)nr2*1.5/sqrt(3); yord1 = nr1 - (double)(nr2%2)/2.0;
      for (nr3=0;nr3<3;nr3++){
	switch(nr3){
	case 0: nr1a = nr1+1; nr2a = nr2+0; break;
	case 1: nr1a = nr1+0; nr2a = nr2+1; break;
	case 2: nr1a = nr1+(int)pow(-1,nr2%2); nr2a = nr2+1; break;
	default: break;}
	if (nr1a>=0 && nr1a<input1max && nr2a>=0 && nr2a<input2max_scale){
	  xord2 = (double)nr2a*1.5/sqrt(3); yord2 = nr1a - (double)(nr2a%2)/2.0;
	  xord3 = (xord1+xord2)/2.0 + sqrt(3)/6.0*(yord2-yord1); yord3 = (yord1+yord2)/2.0 + sqrt(3)/6.0*(xord1-xord2);
	  xord4 = (xord1+xord2)/2.0 - sqrt(3)/6.0*(yord2-yord1); yord4 = (yord1+yord2)/2.0 - sqrt(3)/6.0*(xord1-xord2);
	  xord5 = (xord1+xord2+xord3+xord4)/4.0; yord5 = (yord1+yord2+yord3+yord4)/4.0;
	  tab = nr3+nr1*3+nr2*3*input1max; periodify("int",&tab,&depthmin,&depthmax,&depth);
	  colorscale(0,ra_reweight[tab],max,min,&rcolor,&gcolor,&bcolor);colorcode=crop((int)floor(512*rcolor),0,511);
	  fprintf(fp,"2 1 0 0 0 %d %d 0 %d 0.000 0 0 -1 0 0 %d\n",/*colorcode*/colorcode+32,/*depth*/depth,/*fill*/20,/*npoints*/4+1); fprintf(fp,"\t"); fprintf(fp,"%d %d ",(int)floor(maxdia*xord1/(double)input2max),(int)floor(maxdia*yord1/(double)input1max)); fprintf(fp,"%d %d ",(int)floor(maxdia*xord3/(double)input2max),(int)floor(maxdia*yord3/(double)input1max)); fprintf(fp,"%d %d ",(int)floor(maxdia*xord2/(double)input2max),(int)floor(maxdia*yord2/(double)input1max)); fprintf(fp,"%d %d ",(int)floor(maxdia*xord4/(double)input2max),(int)floor(maxdia*yord4/(double)input1max)); fprintf(fp,"%d %d ",(int)floor(maxdia*xord1/(double)input2max),(int)floor(maxdia*yord1/(double)input1max)); fprintf(fp,"\n"); if (fig_verbose){ fprintf(fp,"4 0 0 %d 0 %d %d 0.0000 4 %d %d %d %d %d%d_%d%d\\001\n",/*depth*/depthmax+1,/*font*/5,/*point*/5,/*textheight*/(int)(5*13.5),/*textwidth*/(int)(5*9*5),/*xpos*/(int)floor(maxdia*xord5/(double)input2max),/*ypos*/(int)floor(maxdia*yord5/(double)input1max),nr1,nr2,nr1a,nr2a);}}}}}    
    if (jpg_flag){ sprintf(command,"fig2dev -Ljpeg -q %d %s.fig %s.jpg;",/*quality*/5,filename1a,filename1a); system(command);}
    if (eps_flag){ sprintf(command,"fig2dev -Leps %s.fig %s.eps;",filename1a,filename1a); system(command);}
    if (remove_flag){ sprintf(command,"rm %s.fig;",filename1a); system(command);}
    if (fp!=stdout){ fclose(fp);}
    stats("double",ra_angle,3*input1max*input2max_scale,&max,&min,&mean,&stdev); max = mean+STD_VIEW*stdev; min = mean-STD_VIEW*stdev;
    max = PI;
    min = 0.0;
    sprintf(filename1a,"./ptree%sangle_adjacencyplot_%srecord",tempchar,gs2);
    sprintf(filename1b,"%s.fig",filename1a);
    if ((fp=fopen(filename1b,"w"))==NULL){ printf(" %% Warning, cannot create %s in granule_plot\n",filename1b); fp=stdout;}
    fprintf(fp,"%s",FIG_PREAMBLE); fprintf(fp,"%s",FIG_PREAMBLE_COLOR_7);
    for (nr1=0;nr1<input1max;nr1++){ for (nr2=0;nr2<input2max_scale;nr2++){
      xord1 = (double)nr2*1.5/sqrt(3); yord1 = nr1 - (double)(nr2%2)/2.0;
      for (nr3=0;nr3<3;nr3++){
	switch(nr3){
	case 0: nr1a = nr1+1; nr2a = nr2+0; break;
	case 1: nr1a = nr1+0; nr2a = nr2+1; break;
	case 2: nr1a = nr1+(int)pow(-1,nr2%2); nr2a = nr2+1; break;
	default: break;}
	if (nr1a>=0 && nr1a<input1max && nr2a>=0 && nr2a<input2max_scale){
	  xord2 = (double)nr2a*1.5/sqrt(3); yord2 = nr1a - (double)(nr2a%2)/2.0;
	  xord3 = (xord1+xord2)/2.0 + sqrt(3)/6.0*(yord2-yord1); yord3 = (yord1+yord2)/2.0 + sqrt(3)/6.0*(xord1-xord2);
	  xord4 = (xord1+xord2)/2.0 - sqrt(3)/6.0*(yord2-yord1); yord4 = (yord1+yord2)/2.0 - sqrt(3)/6.0*(xord1-xord2);
	  xord5 = (xord1+xord2+xord3+xord4)/4.0; yord5 = (yord1+yord2+yord3+yord4)/4.0;
	  tab = nr3+nr1*3+nr2*3*input1max; periodify("int",&tab,&depthmin,&depthmax,&depth);
	  colorscale(0,ra_angle[tab],max,min,&rcolor,&gcolor,&bcolor);colorcode=crop((int)floor(512*rcolor),0,511);
	  fprintf(fp,"2 1 0 0 0 %d %d 0 %d 0.000 0 0 -1 0 0 %d\n",/*colorcode*/colorcode+32,/*depth*/depth,/*fill*/20,/*npoints*/4+1); fprintf(fp,"\t"); fprintf(fp,"%d %d ",(int)floor(maxdia*xord1/(double)input2max),(int)floor(maxdia*yord1/(double)input1max)); fprintf(fp,"%d %d ",(int)floor(maxdia*xord3/(double)input2max),(int)floor(maxdia*yord3/(double)input1max)); fprintf(fp,"%d %d ",(int)floor(maxdia*xord2/(double)input2max),(int)floor(maxdia*yord2/(double)input1max)); fprintf(fp,"%d %d ",(int)floor(maxdia*xord4/(double)input2max),(int)floor(maxdia*yord4/(double)input1max)); fprintf(fp,"%d %d ",(int)floor(maxdia*xord1/(double)input2max),(int)floor(maxdia*yord1/(double)input1max)); fprintf(fp,"\n"); if (fig_verbose){ fprintf(fp,"4 0 0 %d 0 %d %d 0.0000 4 %d %d %d %d %d%d_%d%d\\001\n",/*depth*/depthmax+1,/*font*/5,/*point*/5,/*textheight*/(int)(5*13.5),/*textwidth*/(int)(5*9*5),/*xpos*/(int)floor(maxdia*xord5/(double)input2max),/*ypos*/(int)floor(maxdia*yord5/(double)input1max),nr1,nr2,nr1a,nr2a);}}}}}    
    if (jpg_flag){ sprintf(command,"fig2dev -Ljpeg -q %d %s.fig %s.jpg;",/*quality*/5,filename1a,filename1a); system(command);}
    if (eps_flag){ sprintf(command,"fig2dev -Leps %s.fig %s.eps;",filename1a,filename1a); system(command);}
    if (remove_flag){ sprintf(command,"rm %s.fig;",filename1a); system(command);}
    if (fp!=stdout){ fclose(fp);}
    tfree(ra);ra=NULL;
    tfree(ra_angle); ra_angle=NULL;
    tfree(ra_reweight); ra_reweight=NULL;}
}

double granule_obtain_reweight_vs_bipoo_angle(struct ptree *pW,double *ra_bipoo)
{
  /* returns angle between second level of pW->postree (relevance) and ra, as vectors */
  int verbose=0;
  int nr=0,nc=0,tab=0;
  double *ra_reweight,temp_d=0,temp1=0,temp2=0,angle=0;
  struct llist *LW=NULL;
  if (verbose){ printf(" %% [entering granule_obtain_reweight_vs_bipoo_angle]\n");}
  if (verbose){ raprintf(ra_bipoo,"double",pW->nregions,pW->nregions," %% granule_obtain_reweight_vs_bipoo_angle --> bipoo:");}
  ra_reweight = (double *) tcalloc(pW->nregions*pW->nregions,sizeof(double));
  for (nr=0;nr<pW->nregions;nr++){ for (nc=0;nc<pW->nregions;nc++){
    LW=llistmake();
    litemadd(LW,pW->regionra[nr]); litemadd(LW,pW->regionra[nc]);
    pnode_dig_and_evaluate(1,LW->first,pW,NULL,pW->postree,&pnode_evaluate_relevance,&(ra_reweight[nr+nc*pW->nregions]));
    llisttfree(LW);LW=NULL;}}
  if (verbose){ raprintf(ra_reweight,"double",pW->nregions,pW->nregions," %% granule_obtain_reweight_vs_bipoo_angle --> reweight:");}
  temp_d=0;temp1=0;temp2=0;
  for (nr=0;nr<pW->nregions;nr++){ for (nc=0;nc<pW->nregions;nc++){
    tab = nr+nc*pW->nregions;
    temp_d += ra_reweight[tab]*ra_bipoo[tab];
    temp1 += pow(ra_reweight[tab],2);
    temp2 += pow(ra_bipoo[tab],2);}}
  angle = acos(temp_d/(sqrt(temp1)*sqrt(temp2)));
  if (!finite(angle)){ if (verbose){ printf(" %% error! infinite angle, temp1 %f, temp2 %f, temp_d%f\n",temp1,temp2,temp_d);} angle=PI/2;}
  if (verbose){ printf(" %% angle %f\n",angle);}
  tfree(ra_reweight);ra_reweight=NULL;
  return angle;
}

double * granule_obtain_bipoo_matrix(char *filename1,double time_total,double nrecords)
{
  /* assume that filename? starts with a "./" */
  int verbose=0;
  struct ptree *p1=NULL;
  double *ra=NULL;
  void **vrara=NULL;
  int level = 2;
  int no_error1=0;
  int nr=0,nc=0;
  if (verbose){ printf(" %% [entering granule_obtain_bipoo_matrix]\n");}
  if (checktofind(filename1)){
    if (verbose){ printf(" %% base file %s\n",filename1);}
    p1=ptree_readbacktemp(filename1,&file2obsdisthist,&no_error1);}
  if (no_error1){
    if (verbose){ pnodeprintf(NULL,p1->postree,-1,0);}
    ra = (double *) tcalloc(p1->nregions*p1->nregions,sizeof(double));
    vrara = (void **) tcalloc(5,sizeof(void *));
    vrara[0] = p1; vrara[1] = ra; vrara[2] = &level; vrara[3] = &time_total; vrara[4] = &nrecords;
    pnode2pnode_void_operate(p1,NULL,p1->postree,NULL,NULL,NULL,2,0,&pnode_obsdisthist2mtrx,vrara,0);
    if (verbose){ printf(" %% lvl2_ra read to be:\n"); raprintf(ra,"double",p1->nregions,p1->nregions,"%% lvl2:");}
    for (nr=0;nr<p1->nregions;nr++){ for (nc=nr;nc<p1->nregions;nc++){ 
      ra[nr+nc*p1->nregions] = ra[nr+nc*p1->nregions]-ra[nc+nr*p1->nregions];
      ra[nc+nr*p1->nregions] = -ra[nr+nc*p1->nregions];}}
    if (verbose){ printf(" %% bipoo_ra read to be:\n"); raprintf(ra,"double",p1->nregions,p1->nregions,"%% bipoo:");}
    pnode_obs2dist_starter(0,0,0,NULL,p1->postree,-1,0,2); ptreetfree(p1);p1=NULL;
    tfree(vrara);vrara=NULL;}
  else /* if (!no_error1 || !no_error2) */{ printf(" %% warning! improperly read p1 in granule_obtain_bipoo_matrix\n");}
  return ra;
}

struct ptree * ptree_obtain_bipoo_reweight(char *filename1)
{
  /* This uses the ptree stored in filename1 to construct a full bipoo_reweight ptree pW
     assumes that filename1 starts with a "./" */
  int verbose=0;
  struct ptree *p1=NULL,*pW=NULL;
  int nr1=0,nr2=0;
  struct llist *L12=NULL,*L21=NULL;
  double weight12=0,weight21=0;
  struct llitem *l1=NULL;
  struct pnode *pn=NULL;
  if (verbose){ printf(" %% [entering ptree_obtain_bipoo_reweight]\n");}
  if (checktofind(filename1)){
    if (verbose){ printf(" %% base file %s\n",filename1);}
    p1=ptreadback(filename1,0);
    if (verbose){ pnodeprintf(NULL,p1->postree,-1,0);}
    pW = ptreemake(p1->nregions,1,1,0,1,p1->legtime);
    for (nr1=0;nr1<p1->nregions;nr1++){ for (nr2=0;nr2<p1->nregions;nr2++){
      L12=llistmake(); litemadd(L12,p1->regionra[nr1]); litemadd(L12,p1->regionra[nr2]);
      L21=llistmake(); litemadd(L21,p1->regionra[nr2]); litemadd(L21,p1->regionra[nr1]);
      pnode_dig_and_evaluate(1,L12->first,p1,NULL,p1->postree,&pnode_evaluate_weight,&weight12);
      pnode_dig_and_evaluate(1,L21->first,p1,NULL,p1->postree,&pnode_evaluate_weight,&weight21);
      if (verbose){ printf(" %% evaluating region1 %d region2 %d ... weight12 %0.1f weight21 %0.1f\n",nr1,nr2,weight12,weight21);}
      llisttfree(L12);L12=NULL;
      llisttfree(L21);L21=NULL;      
      if ((l1=llitemaddorfind(0,pW->postree,pW->regionra[nr1],&region2pnode_compare_label))==NULL){
	l1 = llitemaddorfind(1,pW->postree,pnodemake(NULL,pW->regionra[nr1],0,0),&pnode2pnode_compare_label);}
      pn = (struct pnode *)l1->item;
      if ((l1=llitemaddorfind(0,pn->childllitem,pW->regionra[nr2],&region2pnode_compare_label))==NULL){
	l1 = llitemaddorfind(1,pn->childllitem,pnodemake(pn,pW->regionra[nr2],0,0),&pnode2pnode_compare_label);}
      pn = (struct pnode *)l1->item; pn->relevance = (weight12-weight21)/p1->total_time;
      if ((l1=llitemaddorfind(0,pW->postree,pW->regionra[nr2],&region2pnode_compare_label))==NULL){
	l1 = llitemaddorfind(1,pW->postree,pnodemake(NULL,pW->regionra[nr2],0,0),&pnode2pnode_compare_label);}
      pn = (struct pnode *)l1->item;
      if ((l1=llitemaddorfind(0,pn->childllitem,pW->regionra[nr1],&region2pnode_compare_label))==NULL){
	l1 = llitemaddorfind(1,pn->childllitem,pnodemake(pn,pW->regionra[nr1],0,0),&pnode2pnode_compare_label);}
      pn = (struct pnode *)l1->item; pn->relevance = (weight21-weight12)/p1->total_time;}}
    if (verbose){ pnodeprintf(NULL,pW->postree,-1,0);}
    ptreetfree(p1);p1=NULL;}
  else /* if (!checktofind(filename1)) */{ printf(" %% warning! improperly read p1 in ptree_obtain_bipoo_reweight\n");}
  return pW;
}

void ptree_test_mcpit()
{
  int verbose=GLOBAL_verbose;
  int nregions=0,tab=0,nr1=0,nr2=0,nr3=0,nv=0,nl=0,nl1=0,nl2=0,nt1=0,nt2=0,nb=0,nodor=0,nreweight=0;
  double tolerance=0.000001,temp1=0;
  char filename_reweight[1024],filename_basic_reweight0[1024],filename_basic_reweight1[1024],filename[1024],tempchar[1024];
  struct ptree *p0=NULL,*p1=NULL,*pW=NULL;
  double *connectivity_full=NULL,*connectivity_mcpit_effective=NULL,*connectivity_full_skeleton=NULL,*event_tree_1=NULL,*event_tree_2=NULL,*event_tree_3=NULL,*event_tree_1_shift_basic_reweight=NULL,background_rate_estimated=0,*connectivity_mcpit_estimated=NULL,*asym_connectivity_mcpit_estimated=NULL,*connectivity_reweight=NULL,*connectivity_mcpit_reweight_effective=NULL;
  double *event_tree_1_basic_reweight=NULL,*event_tree_2_basic_reweight=NULL,*event_tree_3_basic_reweight=NULL,background_rate_estimated_basic_reweight=0,*connectivity_mcpit_basic_reweight_estimated=NULL;
  struct neuronarray *Nra=GLOBAL_Nra;
  struct neuron *n1=NULL,*n2=NULL;
  struct llist *L=NULL;
  int *ira=NULL;
  double *ra=NULL,*ra2=NULL,*ra3=NULL,*ra4=NULL,*ra5=NULL,*input_shift=NULL;
  int reweight_strength_old=0;
  struct ptree *reweight_old=0;
  struct odor *o1=NULL,*o2=NULL;
  int found_flag=0;
  sprintf(filename_basic_reweight0,"./ptree_%srecord_basic_reweight0",GLOBAL_STRING_2);
  p0 = ptreadback(filename_basic_reweight0,0); if (verbose>2){ pnodeprintf(NULL,p0->postree,-1,0);} 
  nregions = p0->nregions;
  connectivity_full=(double *) tcalloc(Nra->nsval*nregions*nregions,sizeof(double));
  connectivity_full_skeleton=(double *) tcalloc(nregions*nregions,sizeof(double));
  for (nr1=0;nr1<nregions;nr1++){ for (nr2=0;nr2<nregions;nr2++){
    nl1=0;nt1=0;nl=nr1;
    while (nl>0){ if (nl>=Nra->lengthra[nt1]){ nl-=Nra->lengthra[nt1]; nt1+=1;} else if (nl<Nra->lengthra[nt1]){ nl1=nl;nl=0;}}
    n1=nget(Nra,nt1,nl1);
    nl2=0;nt2=0;nl=nr2;
    while (nl>0){ if (nl>=Nra->lengthra[nt2]){ nl-=Nra->lengthra[nt2]; nt2+=1;} else if (nl<Nra->lengthra[nt2]){ nl2=nl;nl=0;}}
    n2=nget(Nra,nt2,nl2);
    reweight_strength_old = GLOBAL_REWEIGHT_STRENGTH; GLOBAL_REWEIGHT_STRENGTH=0;
    connectivity_full_skeleton[nr1+nr2*nregions]=slink(n1,n2,NULL,&(connectivity_full[0+nr1*Nra->nsval+nr2*nregions*Nra->nsval]));
    GLOBAL_REWEIGHT_STRENGTH = reweight_strength_old;}}
  if (verbose>1){ raprintf(connectivity_full,"double",Nra->nsval*nregions,nregions,"connectivity_full: ");}
  if (verbose>1){raprintf(connectivity_full_skeleton,"double",nregions,nregions,"connectivity_full_skeleton: ");}
  connectivity_mcpit_effective=(double *) tcalloc(nregions*nregions,sizeof(double)); 
  for (nr1=0;nr1<nregions;nr1++){ for (nr2=0;nr2<nregions;nr2++){
    connectivity_mcpit_effective[nr1+nr2*nregions] = 0;
    for (nv=0;nv<Nra->nsval;nv++){
      connectivity_mcpit_effective[nr1+nr2*nregions] += connectivity_full[nv+nr1*Nra->nsval+nr2*nregions*Nra->nsval]*(VOLTAGE_[GLOBAL_INDEXING_CHECKOUT_sra[nv]]-VOLTAGE_THRESHOLD_S);}}}
  if (verbose>1){ raprintf(connectivity_mcpit_effective,"double",nregions,nregions,"connectivity_mcpit_effective: ");}
  nb = (int)(p0->total_time/p0->legtime);
  if (verbose>1){ printf(" %% ptree_basic_reweight0 total time %f, nbins %d\n",p0->total_time,nb);}
  event_tree_1=(double *)tcalloc(nregions,sizeof(double));
  for (nr1=0;nr1<nregions;nr1++){ 
    L=llistmake();litemadd(L,p0->regionra[nr1]);
    pnode_dig_and_evaluate(1,L->first,p0,NULL,p0->postree,&pnode_evaluate_weight,&(event_tree_1[nr1])); event_tree_1[nr1] /= nb;
    llisttfree(L);L=NULL;}
  if (verbose>1){ raprintf(event_tree_1,"double",nregions,1,"event_tree_1: ");}
  event_tree_2=(double *)tcalloc(nregions*nregions,sizeof(double));
  for (nr1=0;nr1<nregions;nr1++){ for (nr2=0;nr2<nregions;nr2++){
    L=llistmake();litemadd(L,p0->regionra[nr1]);litemadd(L,p0->regionra[nr2]);
    pnode_dig_and_evaluate(1,L->first,p0,NULL,p0->postree,&pnode_evaluate_weight,&(event_tree_2[nr1+nr2*nregions])); 
    event_tree_2[nr1+nr2*nregions] /= nb;
    llisttfree(L);L=NULL;}}
  if (verbose>1){ raprintf(event_tree_2,"double",nregions,nregions,"event_tree_2: ");}
  event_tree_3=(double *)tcalloc(nregions*nregions*nregions,sizeof(double));
  for (nr1=0;nr1<nregions;nr1++){ for (nr2=0;nr2<nregions;nr2++){ for (nr3=0;nr3<nregions;nr3++){
    L=llistmake();litemadd(L,p0->regionra[nr1]);litemadd(L,p0->regionra[nr2]);litemadd(L,p0->regionra[nr3]);
    pnode_dig_and_evaluate(1,L->first,p0,NULL,p0->postree,&pnode_evaluate_weight,&(event_tree_3[nr1+nr2*nregions+nr3*nregions*nregions])); 
    event_tree_3[nr1+nr2*nregions+nr3*nregions*nregions] /= nb;
    llisttfree(L);L=NULL;}}}
  if (verbose>1){ for (nr1=0;nr1<nregions;nr1++){ raprintf(&(event_tree_3[nr1*nregions*nregions]),"double",nregions,nregions,"event_tree_3[nr1]: ");}}
  tab=0;background_rate_estimated=0;
  for (nr1=0;nr1<nregions;nr1++){ for (nr2=0;nr2<nregions;nr2++){ for (nr3=0;nr3<nregions;nr3++){
    if (nr1!=nr2 && nr1!=nr3 && nr2!=nr3){
      temp1 = event_tree_2[nr1+nr2*nregions] - event_tree_1[nr1]*event_tree_1[nr2] - event_tree_2[nr2+nr3*nregions] + event_tree_1[nr2]*event_tree_1[nr3];
      if (fabs(temp1)>tolerance){
	tab += 1;
	background_rate_estimated += (event_tree_3[nr1+nr2*nregions+nr3*nregions*nregions]-event_tree_1[nr1]*event_tree_1[nr2]*event_tree_1[nr3])/temp1;}}}}}
  background_rate_estimated /= maximum(1,tab);
  if (verbose>=0){ printf(" background_rate_estimated = %0.16f\n",background_rate_estimated);}
  stats("double",event_tree_1,nregions,NULL,NULL,&background_rate_estimated,NULL);
  if (verbose>=0){ printf(" background_rate_estimated = %0.16f\n",background_rate_estimated);}
  connectivity_mcpit_estimated=(double *)tcalloc(nregions*nregions,sizeof(double));
  for (nr1=0;nr1<nregions;nr1++){ for (nr2=0;nr2<nregions;nr2++){
    if (nr1!=nr2){ 
      connectivity_mcpit_estimated[nr1+nr2*nregions] = event_tree_2[nr1+nr2*nregions]-event_tree_1[nr1]*event_tree_1[nr2];
      connectivity_mcpit_estimated[nr1+nr2*nregions] /= background_rate_estimated*(1-background_rate_estimated);}}}
  if (verbose>1){ raprintf(connectivity_mcpit_estimated,"double",nregions,nregions,"connectivity_mcpit_estimated: ");}
  sprintf(filename,"./connectivity_full_%srecord",GLOBAL_STRING_2);
  ira=(int *)tcalloc(3,sizeof(int));ira[0]=Nra->nsval;ira[1]=nregions;ira[2]=nregions;
  multidradump(connectivity_full,3,ira,filename);
  tfree(ira);ira=NULL;
  if (verbose>1){ 
    multidraread(filename,&ira,&ra); raprintf(ra,"double",ira[0]*ira[1],ira[2],"connectivity_full: "); 
    tfree(ira);ira=NULL;tfree(ra);ra=NULL;}
  sprintf(filename,"./connectivity_mcpit_effective_%srecord",GLOBAL_STRING_2);
  ira=(int *)tcalloc(2,sizeof(int));ira[0]=nregions;ira[1]=nregions;
  multidradump(connectivity_mcpit_effective,2,ira,filename);
  tfree(ira);ira=NULL;
  if (verbose>1){ 
    multidraread(filename,&ira,&ra); raprintf(ra,"double",ira[0],ira[1],"connectivity_mcpit_effective: "); 
    tfree(ira);ira=NULL;tfree(ra);ra=NULL;}
  sprintf(filename,"./connectivity_mcpit_estimated_%srecord",GLOBAL_STRING_2);
  ira=(int *)tcalloc(2,sizeof(int));ira[0]=nregions;ira[1]=nregions;
  multidradump(connectivity_mcpit_estimated,2,ira,filename);
  tfree(ira);ira=NULL;
  if (verbose>1){ 
    multidraread(filename,&ira,&ra); raprintf(ra,"double",ira[0],ira[1],"connectivity_mcpit_estimated: "); 
    tfree(ira);ira=NULL;tfree(ra);ra=NULL;}
  if (verbose>=0){ printf("correlation between connectivity_mcpit_estimated and connectivity_mcpit_effective = %f\n",correlation(connectivity_mcpit_estimated,connectivity_mcpit_effective,nregions*nregions));}
  if (p0!=NULL){ ptreetfree(p0);p0=NULL;}
  sprintf(filename_reweight,"./ptree_%srecord_reweight",GLOBAL_STRING_2);
  pW = ptreadback(filename_reweight,0); if (verbose>2){ pnodeprintf(NULL,pW->postree,-1,0);}
  connectivity_reweight=(double *) tcalloc(Nra->nsval*nregions*nregions,sizeof(double));
  for (nr1=0;nr1<nregions;nr1++){ for (nr2=0;nr2<nregions;nr2++){
    nl1=0;nt1=0;nl=nr1;
    while (nl>0){ if (nl>=Nra->lengthra[nt1]){ nl-=Nra->lengthra[nt1]; nt1+=1;} else if (nl<Nra->lengthra[nt1]){ nl1=nl;nl=0;}}
    n1=nget(Nra,nt1,nl1);
    nl2=0;nt2=0;nl=nr2;
    while (nl>0){ if (nl>=Nra->lengthra[nt2]){ nl-=Nra->lengthra[nt2]; nt2+=1;} else if (nl<Nra->lengthra[nt2]){ nl2=nl;nl=0;}}
    n2=nget(Nra,nt2,nl2);
    reweight_old = Nra->reweight; Nra->reweight=pW;
    reweight_strength_old = GLOBAL_REWEIGHT_STRENGTH; GLOBAL_REWEIGHT_STRENGTH=SUITE_4_REWEIGHT_STRENGTH;
    slink(n1,n2,NULL,&(connectivity_reweight[0+nr1*Nra->nsval+nr2*nregions*Nra->nsval]));
    Nra->reweight = reweight_old;
    GLOBAL_REWEIGHT_STRENGTH = reweight_strength_old;}}
  if (verbose>1){ raprintf(connectivity_reweight,"double",Nra->nsval*nregions,nregions,"connectivity_reweight: ");}
  connectivity_mcpit_reweight_effective=(double *) tcalloc(nregions*nregions,sizeof(double)); 
  for (nr1=0;nr1<nregions;nr1++){ for (nr2=0;nr2<nregions;nr2++){
    tab=nr1+nr2*nregions; connectivity_mcpit_reweight_effective[tab] = 0;
    for (nv=0;nv<Nra->nsval;nv++){
      connectivity_mcpit_reweight_effective[tab] += connectivity_reweight[nv+nr1*Nra->nsval+nr2*nregions*Nra->nsval]*(VOLTAGE_[GLOBAL_INDEXING_CHECKOUT_sra[nv]]-VOLTAGE_THRESHOLD_S);}
    connectivity_mcpit_reweight_effective[tab] -= connectivity_mcpit_effective[tab];}}
  if (verbose>1){ raprintf(connectivity_mcpit_reweight_effective,"double",nregions,nregions,"connectivity_mcpit_reweight_effective: ");}
  asym_connectivity_mcpit_estimated = (double *) tcalloc(nregions*nregions,sizeof(double));
  for (nr1=0;nr1<nregions;nr1++){ for (nr2=0;nr2<nregions;nr2++){ asym_connectivity_mcpit_estimated[nr1+nr2*nregions]=connectivity_mcpit_estimated[nr1+nr2*nregions]-connectivity_mcpit_estimated[nr2+nr1*nregions];}}
  if (verbose>1){ raprintf(asym_connectivity_mcpit_estimated,"double",nregions,nregions,"asym_connectivity_mcpit_estimated: ");}
  if (verbose>=0){ printf("correlation between connectivity_mcpit_reweight_effective and asym_connectivity_mcpit_estimated = %f\n",correlation(connectivity_mcpit_reweight_effective,asym_connectivity_mcpit_estimated,nregions*nregions));}
  if (pW!=NULL){ ptreetfree(pW);pW=NULL;}
  sprintf(filename_basic_reweight1,"./ptree_%srecord_basic_reweight1",GLOBAL_STRING_2);
  p1 = ptreadback(filename_basic_reweight1,0); if (verbose>2){ pnodeprintf(NULL,p1->postree,-1,0);}
  nb = (int)(p1->total_time/p1->legtime);
  if (verbose>1){ printf(" %% ptree_basic_reweight1 total time %f, nbins %d\n",p1->total_time,nb);}
  event_tree_1_basic_reweight=(double *)tcalloc(nregions,sizeof(double));
  for (nr1=0;nr1<nregions;nr1++){ 
    L=llistmake();litemadd(L,p1->regionra[nr1]);
    pnode_dig_and_evaluate(1,L->first,p1,NULL,p1->postree,&pnode_evaluate_weight,&(event_tree_1_basic_reweight[nr1])); event_tree_1_basic_reweight[nr1] /= nb;
    llisttfree(L);L=NULL;}
  if (verbose>1){ raprintf(event_tree_1_basic_reweight,"double",nregions,1,"event_tree_1_basic_reweight: ");}
  event_tree_2_basic_reweight=(double *)tcalloc(nregions*nregions,sizeof(double));
  for (nr1=0;nr1<nregions;nr1++){ for (nr2=0;nr2<nregions;nr2++){
    L=llistmake();litemadd(L,p1->regionra[nr1]);litemadd(L,p1->regionra[nr2]);
    pnode_dig_and_evaluate(1,L->first,p1,NULL,p1->postree,&pnode_evaluate_weight,&(event_tree_2_basic_reweight[nr1+nr2*nregions])); 
    event_tree_2_basic_reweight[nr1+nr2*nregions] /= nb;
    llisttfree(L);L=NULL;}}
  if (verbose>1){ raprintf(event_tree_2_basic_reweight,"double",nregions,nregions,"event_tree_2_basic_reweight: ");}
  event_tree_3_basic_reweight=(double *)tcalloc(nregions*nregions*nregions,sizeof(double));
  for (nr1=0;nr1<nregions;nr1++){ for (nr2=0;nr2<nregions;nr2++){ for (nr3=0;nr3<nregions;nr3++){
    L=llistmake();litemadd(L,p1->regionra[nr1]);litemadd(L,p1->regionra[nr2]);litemadd(L,p1->regionra[nr3]);
    pnode_dig_and_evaluate(1,L->first,p1,NULL,p1->postree,&pnode_evaluate_weight,&(event_tree_3_basic_reweight[nr1+nr2*nregions+nr3*nregions*nregions])); 
    event_tree_3_basic_reweight[nr1+nr2*nregions+nr3*nregions*nregions] /= nb;
    llisttfree(L);L=NULL;}}}
  if (verbose>1){ for (nr1=0;nr1<nregions;nr1++){ raprintf(&(event_tree_3_basic_reweight[nr1*nregions*nregions]),"double",nregions,nregions,"event_tree_3_basic_reweight[nr1]: ");}}
  tab=0;background_rate_estimated_basic_reweight=0;
  for (nr1=0;nr1<nregions;nr1++){ for (nr2=0;nr2<nregions;nr2++){ for (nr3=0;nr3<nregions;nr3++){
    if (nr1!=nr2 && nr1!=nr3 && nr2!=nr3){
      temp1 = event_tree_2_basic_reweight[nr1+nr2*nregions] - event_tree_1_basic_reweight[nr1]*event_tree_1_basic_reweight[nr2] - event_tree_2_basic_reweight[nr2+nr3*nregions] + event_tree_1_basic_reweight[nr2]*event_tree_1_basic_reweight[nr3];
      if (fabs(temp1)>tolerance){
	tab += 1;
	background_rate_estimated_basic_reweight += (event_tree_3_basic_reweight[nr1+nr2*nregions+nr3*nregions*nregions]-event_tree_1_basic_reweight[nr1]*event_tree_1_basic_reweight[nr2]*event_tree_1_basic_reweight[nr3])/temp1;}}}}}
  background_rate_estimated_basic_reweight /= maximum(1,tab);
  if (verbose>=0){ printf(" background_rate_estimated_basic_reweight = %0.16f\n",background_rate_estimated_basic_reweight);}
  stats("double",event_tree_1_basic_reweight,nregions,NULL,NULL,&background_rate_estimated_basic_reweight,NULL);
  if (verbose>=0){ printf(" background_rate_estimated_basic_reweight = %0.16f\n",background_rate_estimated_basic_reweight);}
  connectivity_mcpit_basic_reweight_estimated=(double *)tcalloc(nregions*nregions,sizeof(double));
  for (nr1=0;nr1<nregions;nr1++){ for (nr2=0;nr2<nregions;nr2++){
    if (nr1!=nr2){ 
      connectivity_mcpit_basic_reweight_estimated[nr1+nr2*nregions] = event_tree_2_basic_reweight[nr1+nr2*nregions]-event_tree_1_basic_reweight[nr1]*event_tree_1_basic_reweight[nr2];
      connectivity_mcpit_basic_reweight_estimated[nr1+nr2*nregions] /= background_rate_estimated_basic_reweight*(1-background_rate_estimated_basic_reweight);}}}
  if (verbose>1){ raprintf(connectivity_mcpit_basic_reweight_estimated,"double",nregions,nregions,"connectivity_mcpit_basic_reweight_estimated: ");}
  ra = ra2ra_minus(connectivity_mcpit_basic_reweight_estimated,connectivity_mcpit_estimated,nregions*nregions);
  if (verbose>1){ raprintf(ra,"double",nregions,nregions,"connectivity_mcpit_estimated shift: ");}
  if (verbose>=0){ printf("correlation between connectivity_mcpit_reweight_effective and asym_connectivity_mcpit_estimated = %f\n",correlation(connectivity_mcpit_reweight_effective,asym_connectivity_mcpit_estimated,nregions*nregions));}
  if (verbose>=0){ printf("correlation between connectivity_mcpit_reweight_effective and connectivity_mcpit_estimated shift = %f\n",correlation(connectivity_mcpit_reweight_effective,ra,nregions*nregions));}
  if (verbose>=0){ printf("correlation between connectivity_mcpit_estimated shift and asym_connectivity_mcpit_estimated = %f\n",correlation(ra,asym_connectivity_mcpit_estimated,nregions*nregions));}
  tfree(ra);ra=NULL;
  event_tree_1_shift_basic_reweight = (double *) tcalloc(nregions,sizeof(double));
  for (nr1=0;nr1<nregions;nr1++){ 
    L=llistmake();litemadd(L,p1->regionra[nr1]);
    pnode_dig_and_evaluate(1,L->first,p1,NULL,p1->postree,&pnode_evaluate_weight,&(event_tree_1_shift_basic_reweight[nr1])); 
    event_tree_1_shift_basic_reweight[nr1] /= nb;
    event_tree_1_shift_basic_reweight[nr1] -= event_tree_1[nr1];
    llisttfree(L);L=NULL;}
  if (verbose>1){ raprintf(event_tree_1_shift_basic_reweight,"double",nregions,1,"event_tree_1 for ptree_basic_reweight1 - event_tree_1 for ptree_basic_reweight0: ");}
  ra2 = (double *) tcalloc(nregions,sizeof(double));
  for (nr1=0;nr1<nregions;nr1++){ 
    ra2[nr1]=0;
    for (nr2=0;nr2<nregions;nr2++){ ra2[nr1] += event_tree_2[nr2+nr1*nregions]-event_tree_2[nr1+nr2*nregions];}}
  if (verbose>=0){ printf("correlation between event_tree_1_shift_basic_reweight and predicted shift = %f\n",correlation(event_tree_1_shift_basic_reweight,ra2,nregions));}
  tfree(ra2);ra2=NULL;
  if (p1!=NULL){ ptreetfree(p1);p1=NULL;}
  nodor=0;
  sprintf(tempchar,"odor_%srecord_odor%d_reweight%d",GLOBAL_STRING_2,nodor,0); found_flag = checktofind(tempchar);
  sprintf(tempchar,"odor_%srecord_odor%d_reweight%d",GLOBAL_STRING_2,nodor,1); found_flag = (found_flag && checktofind(tempchar));
  while (found_flag){
    nreweight=0;
    if (verbose>1){ printf(" %% analysing nodor %d reweight %d\n",nodor,nreweight);}
    sprintf(filename,"./odor_%srecord_basic_reweight0",GLOBAL_STRING_2);
    o1 = odor_full_fread(filename);
    sprintf(filename,"./odor_%srecord_odor%d_reweight%d",GLOBAL_STRING_2,nodor,nreweight);
    o2 = odor_full_fread(filename);
    input_shift = (double *) tcalloc(nregions,sizeof(double));
    for (nr1=0;nr1<nregions;nr1++){ 
      nl1=0;nt1=0;nl=nr1;
      while (nl>0){ if (nl>=Nra->lengthra[nt1]){ nl-=Nra->lengthra[nt1]; nt1+=1;} else if (nl<Nra->lengthra[nt1]){ nl1=nl;nl=0;}}
      input_shift[nr1] = o2->rara[nt1][nl1] - o1->rara[nt1][nl1];}
    if (verbose>1){ raprintf(input_shift,"double",nregions,1,"input_shift: ");}  
    ra2 = (double *) tcalloc(nregions,sizeof(double));
    for (nr1=0;nr1<nregions;nr1++){ 
      ra2[nr1]=input_shift[nr1];
      for (nr2=0;nr2<nregions;nr2++){ ra2[nr1] += input_shift[nr2]*connectivity_mcpit_estimated[nr2+nr1*nregions];}}
    if (verbose>1){
      sprintf(tempchar,"estimated shift in event_tree_1 for odor%d_reweight%d: ",nodor,nreweight);
      raprintf(ra2,"double",nregions,1,tempchar);}
    sprintf(filename,"./ptree_%srecord_odor%d_reweight%d",GLOBAL_STRING_2,nodor,nreweight);
    p1 = ptreadback(filename,0); if (verbose>2){ pnodeprintf(NULL,p1->postree,-1,0);}
    nb = (int)(p1->total_time/p1->legtime);
    if (verbose>1){ printf(" %% ptree_odor%d_reweight%d total time %f, nbins %d\n",nodor,nreweight,p1->total_time,nb);}
    ra3 = (double *)tcalloc(nregions,sizeof(double)); /* will hold the event_tree_1 for odor%d_reweight0 */
    ra=(double *)tcalloc(nregions,sizeof(double));
    for (nr1=0;nr1<nregions;nr1++){ 
      L=llistmake();litemadd(L,p1->regionra[nr1]);
      pnode_dig_and_evaluate(1,L->first,p1,NULL,p1->postree,&pnode_evaluate_weight,&(ra[nr1])); ra[nr1] /= nb;
      ra3[nr1] = ra[nr1];
      ra[nr1] -= event_tree_1[nr1];
      llisttfree(L);L=NULL;}
    if (verbose>1){ 
      sprintf(tempchar,"event_tree_1 for ptree_odor%d_reweight%d: ",nodor,nreweight);
      raprintf(ra3,"double",nregions,1,tempchar);}
    if (verbose>1){ 
      sprintf(tempchar,"event_tree_1 for ptree_odor%d_reweight%d - event_tree_1 for ptree_basic_reweight0: ",nodor,nreweight);
      raprintf(ra,"double",nregions,1,tempchar);}
    if (verbose>=0){ printf("correlation between shift in ptree_odor%d_reweight%d event_tree_1 and predicted shift = %f\n",nodor,nreweight,correlation(ra,ra2,nregions));}
    tfree(ra2);ra2=NULL;
    tfree(ra);ra=NULL;
    ra2 = (double *) tcalloc(nregions*nregions,sizeof(double));
    for (nr1=0;nr1<nregions;nr1++){ for (nr2=0;nr2<nregions;nr2++){
      ra2[nr1+nr2*nregions]=0;
      for (nr3=0;nr3<nregions;nr3++){ 
	ra2[nr1+nr2*nregions] += input_shift[nr3]*(background_rate_estimated*(connectivity_mcpit_estimated[nr3+nr1*nregions] + connectivity_mcpit_estimated[nr3+nr2*nregions]) + (nr3==nr1)*(event_tree_1[nr2] + (1-2*background_rate_estimated)*connectivity_mcpit_estimated[nr3+nr2*nregions]) + (nr3==nr2)*event_tree_1[nr1]);}}}
    if (verbose>1){
      sprintf(tempchar,"estimated shift in event_tree_2 for odor%d_reweight%d: ",nodor,nreweight);
      raprintf(ra2,"double",nregions,nregions,tempchar);}
    ra=(double *)tcalloc(nregions*nregions,sizeof(double));
    for (nr1=0;nr1<nregions;nr1++){ for (nr2=0;nr2<nregions;nr2++){
      L=llistmake();litemadd(L,p1->regionra[nr1]);litemadd(L,p1->regionra[nr2]);
      pnode_dig_and_evaluate(1,L->first,p1,NULL,p1->postree,&pnode_evaluate_weight,&(ra[nr1+nr2*nregions])); ra[nr1+nr2*nregions] /= nb;
      ra[nr1+nr2*nregions] -= event_tree_2[nr1+nr2*nregions];
      llisttfree(L);L=NULL;}}
    if (verbose>1){ 
      sprintf(tempchar,"event_tree_2 for ptree_odor%d_reweight%d - event_tree_1 for ptree_basic_reweight0: ",nodor,nreweight);
      raprintf(ra,"double",nregions,nregions,tempchar);}
    if (verbose>=0){ printf("correlation between shift in ptree_odor%d_reweight%d event_tree_2 and predicted shift = %f\n",nodor,nreweight,correlation(ra,ra2,nregions*nregions));}
    tfree(ra2);ra2=NULL;
    tfree(ra);ra=NULL;
    if (p1!=NULL){ ptreetfree(p1);p1=NULL;}
    nreweight=1;
    if (verbose>1){ printf(" %% analysing nodor %d reweight %d\n",nodor,nreweight);}
    sprintf(filename,"./ptree_%srecord_odor%d_reweight%d",GLOBAL_STRING_2,nodor,nreweight);
    p1 = ptreadback(filename,0); if (verbose>2){ pnodeprintf(NULL,p1->postree,-1,0);}
    nb = (int)(p1->total_time/p1->legtime);
    if (verbose>1){ printf(" %% ptree_odor%d_reweight%d total time %f, nbins %d\n",nodor,nreweight,p1->total_time,nb);}
    ra=(double *)tcalloc(nregions,sizeof(double));
    for (nr1=0;nr1<nregions;nr1++){ 
      L=llistmake();litemadd(L,p1->regionra[nr1]);
      pnode_dig_and_evaluate(1,L->first,p1,NULL,p1->postree,&pnode_evaluate_weight,&(ra[nr1])); ra[nr1] /= nb;
      ra[nr1] -= ra3[nr1];
      llisttfree(L);L=NULL;}
    if (verbose>1){ 
      sprintf(tempchar,"event_tree_1 for ptree_odor%d_reweight%d - event_tree_1 for ptree_odor%d_reweight0: ",nodor,nreweight,nodor);
      raprintf(ra,"double",nregions,1,tempchar);}
    ra2 = (double *) tcalloc(nregions,sizeof(double));
    for (nr1=0;nr1<nregions;nr1++){
      ra2[nr1] = ra[nr1] - event_tree_1_shift_basic_reweight[nr1];}
    ra5 = (double *) tcalloc(nregions,sizeof(double));
    for (nr1=0;nr1<nregions;nr1++){ 
      ra5[nr1]=0;
      for (nr2=0;nr2<nregions;nr2++){
	ra5[nr1] += connectivity_mcpit_estimated[nr1+nr2*nregions];}}
    ra4 = (double *) tcalloc(nregions,sizeof(double));
    for (nr1=0;nr1<nregions;nr1++){
      ra4[nr1]=0;
      for (nr2=0;nr2<nregions;nr2++){
	ra4[nr1] += input_shift[nr2]*(connectivity_mcpit_estimated[nr1+nr2*nregions] - 2*connectivity_mcpit_estimated[nr2+nr1*nregions] + (nr1==nr2)*ra5[nr1]);}}
    if (verbose>=0){ printf("correlation between (e_1_odor_reweight-e_1_odor)-(e_1_basic_reweight-e_1_basic) and predicted shift = %f\n",correlation(ra2,ra4,nregions));}
    tfree(ra2);ra2=NULL;
    tfree(ra);ra=NULL;
    tfree(ra3);ra3=NULL;
    tfree(ra5);ra5=NULL;
    tfree(ra4);ra4=NULL;
    if (p1!=NULL){ ptreetfree(p1);p1=NULL;}
    tfree(input_shift);input_shift=NULL;
    if (o1!=NULL){ odortfree(o1);o1=NULL;} if (o2!=NULL){ odortfree(o2);o2=NULL;}
    nodor += 1;
    sprintf(tempchar,"odor_%srecord_odor%d_reweight%d",GLOBAL_STRING_2,nodor,0); found_flag = checktofind(tempchar);
    sprintf(tempchar,"odor_%srecord_odor%d_reweight%d",GLOBAL_STRING_2,nodor,1); found_flag = (found_flag && checktofind(tempchar));}
  tfree(connectivity_full);tfree(connectivity_mcpit_effective);tfree(connectivity_full_skeleton);tfree(event_tree_1);tfree(event_tree_2);tfree(event_tree_3);tfree(event_tree_1_shift_basic_reweight);tfree(connectivity_mcpit_estimated);tfree(asym_connectivity_mcpit_estimated);tfree(connectivity_reweight);tfree(connectivity_mcpit_reweight_effective);tfree(event_tree_1_basic_reweight);tfree(event_tree_2_basic_reweight);tfree(event_tree_3_basic_reweight);tfree(connectivity_mcpit_basic_reweight_estimated);
}

struct ptree * granule_obtain_reweight(char *filename1,char *filename2,double time_total,double nrecords)
{
  /* assume that filename? starts with a "./" */
  int verbose=0;
  struct ptree *p1=NULL,*p2=NULL,*pW=NULL;
  void **vrara=NULL;
  int nr=0;
  int no_error1=0,no_error2=0;
  if (verbose){ printf(" %% [entering granule_obtain_reweight]\n");}
  if (checktofind(filename1)){
    if (verbose){ printf(" %% base file %s\n",filename1);}
    p1=ptree_readbacktemp(filename1,&file2obsdisthist,&no_error1);}
  if (checktofind(filename2)){
    if (verbose){ printf(" %% base file %s\n",filename2);}
    p2=ptree_readbacktemp(filename2,&file2obsdisthist,&no_error2);}
  if (no_error1 && no_error2){
    pW = ptreemake(p1->nregions,1,1,0,p1->nlegs,p1->legtime);
    vrara = (void **) tcalloc(6,sizeof(void *));
    vrara[0] = p1; vrara[1] = p2; vrara[2] = pW; vrara[3] = &time_total; vrara[4] = &nrecords; nr=1; vrara[5] = &nr;
    pnode2pnode_void_operate(p1,NULL,p1->postree,p2,NULL,p2->postree,2,0,&pnode_obsdisthist2obsdisthist_reweight,vrara,0);
    tfree(vrara);vrara=NULL;}
  else /* if (!no_error1 || !no_error2) */{ printf(" %% warning! improperly read p1 and p2 in granule_obtain_reweight\n");}
  if (no_error1){ pnode_obs2dist_starter(0,0,0,NULL,p1->postree,-1,0,2); ptreetfree(p1);p1=NULL;}
  if (no_error2){ pnode_obs2dist_starter(0,0,0,NULL,p2->postree,-1,0,2); ptreetfree(p2);p2=NULL;}
  return pW;
}

int regionra_neuron2label(struct ptree *p,struct neuronarray *Nra,struct neuron *n)
{
  /* does not require p->regionra[0]->neuronllist to be populated, assumes region_type=1 */
  int verbose=0;
  int nl=0,label=-1;
  if (verbose){ printf(" %% [entering regionra_neuron2label] with neuron (%d,%d)\n",n->type,n->index);}
  label=0; for (nl=0;nl<n->type;nl++){ label += Nra->lengthra[nl];} label += n->index;
  if (label<p->nregions){ /* do nothing */ if (verbose){ printf(" %% label %d out of %d\n",label,p->nregions);}} 
  else /* if (label>=p->nregions) */{ if (verbose){ printf(" %% warning! label out of bounds\n");} label = -1;}
  return label;
}

double reweight2adder(struct ptree *p,struct neuronarray *Nra,struct neuron *s,struct neuron *n,double strength)
{
  int verbose=0;
  double adder=0;
  int labels=0,labeln=0;
  struct llitem *l0=NULL;
  struct pnode *pn=NULL;
  if (verbose){ printf(" %% [entering reweight2adder], neuron link (%d,%d)-->(%d,%d)\n",s->type,s->index,n->type,n->index);}
  labels = regionra_neuron2label(p,Nra,s); labeln = regionra_neuron2label(p,Nra,n);
  if (verbose){ printf(" %% s->label %d, n->label %d\n",labels,labeln);}
  if (labels>=0 && labeln>=0){
    if ((l0=llitemaddorfind(0,p->postree,p->regionra[labels],&region2pnode_compare_label))!=NULL){
      pn = (struct pnode *) l0->item;
      if (verbose){ printf(" (%d,%d)->%d found...",s->type,s->index,labels);}
      if ((l0=llitemaddorfind(0,pn->childllitem,p->regionra[labeln],&region2pnode_compare_label))!=NULL){
	pn = (struct pnode *) l0->item;
	if (verbose){ printf(" (%d,%d)->%d found, weight reweight %f\n",n->type,n->index,labeln,pn->relevance);}
	adder += pn->relevance*strength;}
      else /* if not found */{ if (verbose){ printf(" (%d,%d)->%d not found\n",n->type,n->index,labeln);}}}
    else /* if not found */{ if (verbose){ printf(" (%d,%d)->%d not found\n",s->type,s->index,labels);}}}
  if (verbose){ printf(" %% [finishing reweight2adder], adder %f\n",adder);}
  return adder;
}

double reweight2multiplier(struct ptree *p,struct neuronarray *Nra,struct neuron *s,struct neuron *n,double strength)
{
  int verbose=0;
  double multiplier=1;
  int labels=0,labeln=0;
  struct llitem *l0=NULL;
  struct pnode *pn=NULL;
  if (verbose){ printf(" %% [entering reweight2multiplier], neuron link (%d,%d)-->(%d,%d)\n",s->type,s->index,n->type,n->index);}
  labels = regionra_neuron2label(p,Nra,s); labeln = regionra_neuron2label(p,Nra,n);
  if (verbose){ printf(" %% s->label %d, n->label %d\n",labels,labeln);}
  if (labels>=0 && labeln>=0){
    if ((l0=llitemaddorfind(0,p->postree,p->regionra[labels],&region2pnode_compare_label))!=NULL){
      pn = (struct pnode *) l0->item;
      if (verbose){ printf(" (%d,%d)->%d found...",s->type,s->index,labels);}
      if ((l0=llitemaddorfind(0,pn->childllitem,p->regionra[labeln],&region2pnode_compare_label))!=NULL){
	pn = (struct pnode *) l0->item;
	if (verbose){ printf(" (%d,%d)->%d found, weight reweight %f\n",n->type,n->index,labeln,pn->relevance);}
	multiplier *= exp(pn->relevance*strength);}
      else /* if not found */{ if (verbose){ printf(" (%d,%d)->%d not found\n",n->type,n->index,labeln);}}}
    else /* if not found */{ if (verbose){ printf(" (%d,%d)->%d not found\n",s->type,s->index,labels);}}}
  if (verbose){ printf(" %% [finishing reweight2multiplier], multiplier %f\n",multiplier);}
  return multiplier;
}

double * nsphere_index_read(int ndim,int npoints,int np)
{
  /* will generate diagnostic pictures if np<0 */
  int verbose=0;
  FILE *fp=NULL;
  char filename[1024];
  int nr=0,nptemp=0,nd=0;
  double *ra=NULL,*ra2=NULL,*ra3=NULL;
  if (verbose){ printf(" %% [entering nsphere_index_read], ndim %d npoints %d np %d\n",ndim,npoints,nd);}
  if (np>=0 && np<npoints){
    ra = (double *) tcalloc(ndim,sizeof(double));
    sprintf(filename,"nsphere_index_ndim%d_npoints%d",ndim,npoints);
    if (checktofind(filename)){
      if ((fp=fopen(filename,"r"))==NULL){ printf(" %% Can't open %s in nsphere_index_read\n",filename); fp=stdout;}
      for (nr=0;nr<np;nr++){ fread(&nptemp,sizeof(int),1,fp); fread(ra,sizeof(double),ndim,fp);}
      fread(&nptemp,sizeof(int),1,fp); fread(ra,sizeof(double),ndim,fp);
      if (nptemp!=np){ printf(" warning! fishy matrix read in nsphere_index_read\n");}}
    if (fp!=stdout){ fclose(fp);fp=NULL;}}
  else /* if (np<0 || np>=npoints) */{
    if(verbose){ printf(" %% testing nsphere_index...\n");}
    ra2 = (double *) tcalloc(ndim*npoints,sizeof(double));
    ra3 = (double *) tcalloc(ndim,sizeof(double));
    sprintf(filename,"nsphere_index_ndim%d_npoints%d",ndim,npoints);
    if (checktofind(filename)){
      if ((fp=fopen(filename,"r"))==NULL){ printf(" %% Can't open %s in nsphere_index_read\n",filename); fp=stdout;}
      for (nr=0;nr<npoints;nr++){ 
	fread(&nptemp,sizeof(int),1,fp); fread(ra3,sizeof(double),ndim,fp);
	if (verbose){ printf(" %% at point %d: ",nptemp); raprintf(ra3,"double",1,ndim," ");}
	for (nd=0;nd<ndim;nd++){ ra2[nr+nd*npoints]=ra3[nd];}}
      for (nd=0;nd<ndim;nd++){
	sprintf(filename,"./nsphere_index_ndim%d_npoints%d_dim%d.pnm",ndim,npoints,nd);
	WritePNMfile_color(&(ra2[0+nd*npoints]),(int)sqrt(npoints),(int)sqrt(npoints),+1,-1,filename,7);}}
    tfree(ra3);ra3=NULL;
    tfree(ra2);ra2=NULL;
    if (fp!=stdout){ fclose(fp);fp=NULL;}}
  return ra;
}

void nsphere_index_generate(int ndim,int npoints)
{
  /* generate roughly uniformly distributed points on the nsphere, and store them in a lookup table 
     the first point is the origin (center of the nsphere) */
  int verbose=0;
  int nr=0,nd=0,error_flag=0;
  size_t doublesize=sizeof(double);
  double *ra=NULL,temp=0;
  struct llitem *l0=NULL;
  void **void_parameters=NULL;
  struct llist *L=NULL;
  struct litem *l=NULL;
  char filename[1024];
  FILE *fp=NULL;
  if (verbose){ printf(" %% [entering nsphere_index_generate], ndim %d npoints %d\n",ndim,npoints);}
  ra = (double *) tcalloc(ndim*npoints+ndim,sizeof(double));
  l0=llitemmake();
  void_parameters = (void **) tcalloc(3,sizeof(void *));
  for (nr=0;nr<npoints;nr++){
    do{
      temp=0; for (nd=0;nd<ndim;nd++){ ra[nd+nr*ndim]=randn(); temp += pow(ra[nd+nr*ndim],2);}
      if (temp>0){ for (nd=0;nd<ndim;nd++){ ra[nd+nr*ndim]/=sqrt(temp);} error_flag=0;} else /* if (temp<=0) */{ error_flag=1;}}
    while (error_flag);
    if (verbose){ printf(" at step %d found array:",nr); raprintf(&(ra[0+nr*ndim]),"double",1,ndim," ");}
    void_parameters[0] = &ndim;
    void_parameters[1] = &doublesize;
    void_parameters[2] = &double_compare;
    llitemaddorfind_generic(1,l0,&(ra[0+nr*ndim]),&ra2ra_generic_compare_dictionary,void_parameters);}
  sprintf(filename,"nsphere_index_ndim%d_npoints%d",ndim,npoints);
  if ((fp=fopen(filename,"w"))==NULL){ printf(" %% Can't create %s in nsphere_index_generate\n",filename); fp=stdout;}
  L=llistmake(); llitem2llist(l0,L); nr=0;
  l=L->first; while (l!=NULL){ 
    fwrite(&nr,sizeof(int),1,fp);
    if (l==L->first){ fwrite(&(ra[0+ndim*npoints]),sizeof(double),ndim,fp);}
    else /* if (l!=L->first) */{ fwrite(l->item,sizeof(double),ndim,fp);}
    nr+=1;
    l=l->child;}
  if (fp!=stdout){fclose(fp);fp=NULL;}
  llisttfree(L);L=NULL;
  llitemtfree(l0,NULL);l0=NULL;
  tfree(ra);ra=NULL;
  tfree(void_parameters);void_parameters=NULL;  
}

void synaptic_sphere_findint2input_or_input2int(int direction_flag,struct ptree *pW,int input1max,int *input1)
{
  /* assumes that pW is made, with pW->nregions given */
  int verbose=0;
  int ndim = pW->nregions*pW->nregions-pW->nregions;
  char fgvn[512],filename[1024];
  int nr1=0,nr2=0,nr3=0,notfound_flag=0,input1_temp=0;
  FILE *fp=NULL;
  struct llitem *l0=NULL;
  double *ra=NULL;
  struct pnode *pnW=NULL;
  struct llist *LW=NULL;
  if (verbose){ printf(" %% [entering synaptic_sphere_findint2input_or_input2int] direction_flag %d\n",direction_flag);}
  sprintf(filename,"nsphere_index_ndim%d_npoints%d",ndim,input1max);
  if (!checktofind(filename)){ nsphere_index_generate(ndim,input1max);}
  if (direction_flag==1){
    if (verbose){ printf(" %% searching for index...\n");}
    notfound_flag=0;
    for (nr1=0;nr1<input1max;nr1++){
      sprintf(fgvn,"s%d",nr1);
      sprintf(filename,"./ptree_%s_%srecord",fgvn,GLOBAL_STRING_2);
      if (!checktofind(filename)){ 
	notfound_flag=1; input1_temp=nr1;
	if (input1!=NULL){ *input1=nr1;}
	if (verbose){ printf(" %% input (%d) not yet made in synaptic_sphere_findint2input_or_input2int\n",nr1);}
	if ((fp=fopen(filename,"w"))==NULL){ printf(" %% Warning, couldn't open %s in synaptic_sphere_findint2input_or_input2findint, writing to stdout\n",filename); fp=stdout;}
	if (fp!=stdout){ fclose(fp);}
	nr1=input1max+1; /* break */}
      else /* if (checktofind(filename)) */{
	if (verbose>1){ printf(" %% input (%d) already made in synaptic_sphere_findint2input_or_input2int\n",nr1);}}}
    if (!notfound_flag){ printf(" %% Warning, no work for synaptic_sphere_findint2input_or_input2int, exiting\n"); exit(EXIT_SUCCESS);}
    else /* if (notfound_flag) */{
      ra = nsphere_index_read(ndim,input1max,input1_temp); nr3=0;
      if (verbose){ raprintf(ra,"double",1,ndim," %% ra_reweight read to be: ");}
      for (nr2=0;nr2<pW->nregions;nr2++){ for (nr1=0;nr1<pW->nregions;nr1++){ 
	if (nr2!=nr1){
	  if (verbose){ printf(" %% attempting to search for synapse %d-->%d\n",nr1,nr2);}
	  if ((l0=llitemaddorfind(0,pW->postree,pW->regionra[nr1],&region2pnode_compare_label))!=NULL){ 
	    pnW = (struct pnode *) l0->item;}
	  else /* if not found */{ 
	    l0 = llitemaddorfind(1,pW->postree,pnodemake(NULL,pW->regionra[nr1],input1_temp,0),pnode2pnode_compare_label);
	    pnW = (struct pnode *) l0->item;}
	  pnW->weight = input1_temp; pnW->relevance = 0;
	  if ((l0=llitemaddorfind(0,pnW->childllitem,pW->regionra[nr2],&region2pnode_compare_label))!=NULL){ 
	    pnW = (struct pnode *) l0->item;}
	  else /* if not found */{ 
	    l0 = llitemaddorfind(1,pnW->childllitem,pnodemake(NULL,pW->regionra[nr2],input1_temp,0),pnode2pnode_compare_label);
	    pnW = (struct pnode *) l0->item;}
	  pnW->weight = input1_temp; pnW->relevance = ra[nr3];	  
	  if (verbose>2){ printf(" %% currently pW is:\n"); pnodeprintf(NULL,pW->postree,-1,0);}
	  nr3+=1;}}}
      if (verbose){ printf(" just set input to %d, and reweight to:\n",input1_temp); pnodeprintf(NULL,pW->postree,-1,0);}
      tfree(ra);ra=NULL;
      if (verbose>2){
	ra = (double *) tcalloc(pW->nregions*pW->nregions,sizeof(double));
	for (nr1=0;nr1<pW->nregions;nr1++){ for (nr2=0;nr2<pW->nregions;nr2++){
	  LW=llistmake();
	  litemadd(LW,pW->regionra[nr1]); litemadd(LW,pW->regionra[nr2]);
	  pnode_dig_and_evaluate(1,LW->first,pW,NULL,pW->postree,&pnode_evaluate_relevance,&(ra[nr1+nr2*pW->nregions]));
	  llisttfree(LW);LW=NULL;}}
	raprintf(ra,"double",pW->nregions,pW->nregions," %% reweight:");
	tfree(ra);ra=NULL;}}}
  else if (direction_flag==0){
    if (pW!=NULL && pW->postree!=NULL && pW->postree->item!=NULL){
      pnW = (struct pnode *) pW->postree->item;
      if (verbose){ printf(" %% searching pW, it seems as though input==%d\n",(int)pnW->weight);}
      if (input1!=NULL){ *input1 = (int)pnW->weight;}}
    else /* if NULL */{ printf(" %% warning! pW incomplete in synaptic_sphere_findint2input_or_input2int\n");}}
}

void synaptic_sphere_local_findint(char *filename,int *input1)
{
  /* returns the first _s?_ within *filename */
  int verbose=0;
  int length = strlen(filename);
  int nl=0,nl1=0,nltemp=0;
  char tempchar[1];
  int tempint1=0;
  if (verbose){ printf(" %% [entering synaptic_sphere_local_findint] filename %s of length %d\n",filename,length);}
  sprintf(tempchar," ");
  while (nl<length-1 && !(filename[nl]=='_' && filename[nl+1]=='s')){ nl+=1;}
  if (nl<length-1 && filename[nl]=='_'){
    if (verbose){ printf(" %% found ""_"" in filename %s at position %d\n",filename,nl);}
    nl+=1; *tempchar = filename[nl]; if (verbose){ printf(" %% assuming that %s is an ""s""\n",tempchar);} 
    nl+=1; *tempchar = filename[nl]; nl1=0;
    while (nl+nl1<length && *tempchar!='_'){ nl1+=1; *tempchar = filename[nl+nl1];}
    if (nl+nl1<length && *tempchar=='_'){
      if (verbose){ printf(" %% found nl1 %d in filename %s at position %d+%d\n",nl1,filename,nl,nl1);}
      tempint1=0; 
      for (nltemp=0;nltemp<nl1;nltemp++){ *tempchar = filename[nl+nltemp]; tempint1 += (int)pow(10,nl1-nltemp-1)*atoi(tempchar);}
      if (verbose){ printf(" %% read tempint1 %d\n",tempint1);}
      *tempchar = filename[nl+nl1];
      if (verbose){ printf(" %% assuming that %s is a ""_""\n",tempchar);}}}
  else /**/{ printf(" %% warning! aborted at first step nl %d in synaptic_sphere_local_findint\n",nl);}
  if (input1!=NULL){ *input1=tempint1;}
}

void synaptic_sphere_prediction_check(int ndim,int input1max,char *filename1,char *filename2,double time_total,double nrecords,double *forwards_cosangle,double *backwards_cosangle)
{
  /* assume that filename? starts with a "./"
     vector points from 1 to 2 (i.e., 2-1)
     assumes that nsphere_index_read produces synaptic_reweight which perturbes 1 (and not 2) 
     errors are given in terms of cos(theta), ranging from -1 (antiparallel) to +1 (parallel) */
  int verbose=0;
  struct ptree *p1=NULL,*p2=NULL,*pW=NULL,*pD=NULL;
  int nr1=0,nr2=0,nr3=0,secondorfirstpass=0;
  double *ra1=NULL,*ra2=NULL,*raW_real=NULL,*raW_predict=NULL;
  void **vrara=NULL;
  int no_error1=0,no_error2=0;
  struct llist *LW=NULL;
  double normW=0,normD=0,dotDW=0;
  double cosangle_forward=0,cosangle_backward=0;
  if (verbose){ printf(" %% [entering synaptic_sphere_prediction_check]\n");}
  if (checktofind(filename1)){
    if (verbose){ printf(" %% base file %s\n",filename1);}
    p1=ptree_readbacktemp(filename1,&file2obsdisthist,&no_error1);}
  if (checktofind(filename2)){
    if (verbose){ printf(" %% base file %s\n",filename2);}
    p2=ptree_readbacktemp(filename2,&file2obsdisthist,&no_error2);}
  if (no_error1 && no_error2){
    if (verbose){ printf(" %% both files successfully found\n");}
    synaptic_sphere_local_findint(filename1,&nr1); synaptic_sphere_local_findint(filename2,&nr2);
    if (verbose){ printf(" %% corresponding to points %d and %d respectively\n",nr1,nr2);}
    ra1 = nsphere_index_read(ndim,input1max,nr1); ra2 = nsphere_index_read(ndim,input1max,nr2);
    if (verbose){ raprintf(ra1,"double",1,ndim,"%% ra1:"); raprintf(ra2,"double",1,ndim,"%% ra2:");}
    raW_real = (double *) tcalloc(p1->nregions*p1->nregions,sizeof(double));
    nr3=0; 
    for (nr2=0;nr2<p1->nregions;nr2++){ for (nr1=0;nr1<p1->nregions;nr1++){ 
      if (nr2!=nr1){ raW_real[nr1+nr2*p1->nregions] = ra2[nr3]-ra1[nr3]; nr3+=1;}}}
    pW = ptreemake(p1->nregions,1,1,0,p1->nlegs,p1->legtime);
    vrara = (void **) tcalloc(6,sizeof(void *));
    vrara[0] = p1; vrara[1] = raW_real; vrara[2] = pW; vrara[3] = &time_total; vrara[4] = &nrecords; 
    secondorfirstpass=1; vrara[5] = &secondorfirstpass;
    pnode2pnode_void_operate(p1,NULL,p1->postree,NULL,NULL,NULL,2,0,&pnode_obsdisthist2reweight_obsdisthist_mean_shift,vrara,0);
    secondorfirstpass=0; vrara[5] = &secondorfirstpass;
    pnode2pnode_void_operate(p1,NULL,p1->postree,NULL,NULL,NULL,2,0,&pnode_obsdisthist2reweight_obsdisthist_mean_shift,vrara,0);
    if (verbose){ printf(" %% forwards prediction pW:\n"); pnodeprintf(NULL,pW->postree,-1,0);}
    tfree(vrara);vrara=NULL;
    pD = ptreemake(p1->nregions,1,1,0,p1->nlegs,p1->legtime);
    vrara = (void **) tcalloc(5,sizeof(void *));
    vrara[0] = p2; vrara[1] = p1; vrara[2] = pD; vrara[3] = &time_total; vrara[4] = &nrecords; 
    pnode2pnode_void_operate(p2,NULL,p2->postree,p1,NULL,p1->postree,2,0,&pnode_obsdisthist2obsdisthist_meansubtmean,vrara,0);
    if (verbose){ printf(" %% forwards actual pD:\n"); pnodeprintf(NULL,pD->postree,-1,0);}
    tfree(vrara);vrara=NULL;
    vrara = (void **) tcalloc(6,sizeof(void *));
    normD=0;normW=0;dotDW=0; vrara[0]=NULL;vrara[1]=NULL;vrara[2]=NULL; vrara[3] = &normD; vrara[4] = &normW; vrara[5] = &dotDW;
    pnode2pnode_void_operate(pD,NULL,pD->postree,pW,NULL,pW->postree,2,0,&pnode_normanddots,vrara,0);
    tfree(vrara);vrara=NULL;
    cosangle_forward = dotDW/sqrt(normD*normW);
    if (verbose){ printf(" %% forwards error %f\n",cosangle_forward);}
    if (forwards_cosangle!=NULL){ *forwards_cosangle = cosangle_forward;}
    ptreetfree(pW);pW=NULL;
    ptreetfree(pD);pD=NULL;    
    pW = ptreemake(p1->nregions,1,1,0,p1->nlegs,p1->legtime);
    vrara = (void **) tcalloc(6,sizeof(void *));
    vrara[0] = p1; vrara[1] = p2; vrara[2] = pW; vrara[3] = &time_total; vrara[4] = &nrecords; nr1=0; vrara[5] = &nr1;
    pnode2pnode_void_operate(p1,NULL,p1->postree,p2,NULL,p2->postree,2,0,&pnode_obsdisthist2obsdisthist_reweight,vrara,0);
    if (verbose>2){ printf(" %% backwards prediction pW:\n"); pnodeprintf(NULL,pW->postree,-1,0);}
    raW_predict = (double *) tcalloc(pW->nregions*pW->nregions,sizeof(double));
    for (nr1=0;nr1<pW->nregions;nr1++){ for (nr2=0;nr2<pW->nregions;nr2++){
      LW=llistmake();
      litemadd(LW,pW->regionra[nr1]); litemadd(LW,pW->regionra[nr2]);
      pnode_dig_and_evaluate(1,LW->first,pW,NULL,pW->postree,&pnode_evaluate_relevance,&(raW_predict[nr1+nr2*pW->nregions]));
      llisttfree(LW);LW=NULL;}}
    if (verbose){ raprintf(raW_real,"double",p1->nregions,p1->nregions,"%% raW_real:");}
    for (nr1=0;nr1<p1->nregions;nr1++){ if (raW_real[nr1+nr1*p1->nregions]==0){ raW_predict[nr1+nr1*p1->nregions]=0;}}
    if (verbose){ raprintf(raW_predict,"double",pW->nregions,pW->nregions," %% backwards reweight matrix:");}
    cosangle_backward = radotmean(raW_real,raW_predict,p1->nregions*p1->nregions)/sqrt(radotmean(raW_real,raW_real,p1->nregions*p1->nregions)*radotmean(raW_predict,raW_predict,p1->nregions*p1->nregions));
    if (verbose){ printf(" %% backwards error %f\n",cosangle_backward);}
    if (backwards_cosangle!=NULL){ *backwards_cosangle = cosangle_backward;}
    ptreetfree(pW);pW=NULL;
    tfree(raW_real);raW_real=NULL;
    tfree(raW_predict);raW_predict=NULL;
    tfree(ra1);ra1=NULL; tfree(ra2);ra2=NULL;
    tfree(vrara);vrara=NULL;}
  else /* if (!no_error1 || !no_error2) */{ printf(" %% warning! improperly read p1 and p2 in synaptic_sphere_prediction_check\n");}
  if (no_error1){ pnode_obs2dist_starter(0,0,0,NULL,p1->postree,-1,0,2); ptreetfree(p1);p1=NULL;}
  if (no_error2){ pnode_obs2dist_starter(0,0,0,NULL,p2->postree,-1,0,2); ptreetfree(p2);p2=NULL;}
}

void synaptic_sphere_prediction_plot(int npoints)
{
  /* slow and stupid for now, just calls synaptic_sphere_prediction_check for each pair */
  int nr1=0;
  char filename1[1024],filename2[1024];
  char *gs2=GLOBAL_STRING_2;
  double *ra_fcos=NULL,*ra_bcos=NULL;
  int ndim = GLOBAL_PTREE_NREGIONS*(GLOBAL_PTREE_NREGIONS-1);
  double time_total=SUITE_DUMPEVERY; int nrecords = 1024*SUITE_NSECONDS/SUITE_DUMPEVERY;
  ra_fcos = (double *) tcalloc(npoints,sizeof(double));
  ra_bcos = (double *) tcalloc(npoints,sizeof(double));
  sprintf(filename1,"./ptree_s%d_%srecord",0,gs2);
  for (nr1=0;nr1<npoints;nr1++){
    sprintf(filename2,"./ptree_s%d_%srecord",nr1,gs2);
    synaptic_sphere_prediction_check(ndim,npoints,filename1,filename2,time_total,nrecords,&(ra_fcos[nr1]),&(ra_bcos[nr1]));}
  sprintf(filename1,"./synaptic_sphere_prediction_fcos_ndim%d_npoints%d_%s.pnm",ndim,npoints,gs2);
  WritePNMfile_color(ra_fcos,(int)sqrt(npoints),(int)sqrt(npoints),+1,-1,filename1,7);
  sprintf(filename1,"./synaptic_sphere_prediction_bcos_ndim%d_npoints%d_%s.pnm",ndim,npoints,gs2);
  WritePNMfile_color(ra_bcos,(int)sqrt(npoints),(int)sqrt(npoints),+1,-1,filename1,7);
  tfree(ra_fcos);ra_fcos=NULL;
  tfree(ra_bcos);ra_bcos=NULL;
}

/* Here are system functions */

int cleanupoutput()
{
  /* call this after calling everything else to clean up the output */
  int output_flag=1;
  char command[256];
  char *s=GLOBAL_STRING_2;
  printf(" %% trying to bundle everything into %s.tar.gz\n",s);
  sprintf(command,"nice -19 tar cvf ./%sbundle.tar *%s*;",s,s); output_flag *= system(command);
  sprintf(command,"nice -19 rm *%srecord*;",s); output_flag *= system(command);
  sprintf(command,"nice -19 gzip ./%sbundle.tar;",s); output_flag *= system(command);
  return output_flag;
}

/* Here are the strobe and strobetrace functions */

struct strobe * strobemake(int length,double update_timestep,int cycle_bother,int lpower_bother,int lpower_window_length,int lpower_update_every)
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
  st->lpower_bother = lpower_bother;
  st->lpower_window_length = lpower_window_length;
  st->lpower_window_length = maximum(1,lpower_window_length); while (st->length%st->lpower_window_length>0){ st->lpower_window_length-=1;}
  st->lpower_update_every = maximum(1,lpower_update_every);
  st->lpower_length = st->lpower_window_length*st->length/st->lpower_update_every;
  st->lpowerdata = NULL;
  if (st->lpower_bother){ st->lpowerdata = (double *) tcalloc(st->lpower_length,sizeof(double));}
  return st;
}

void strobeupdate(struct strobe *st,double t,double DT,double val)
{
  /* this averages, rather than strobes, which is better for spike statistics
     in addition, time-steps are resolved accurately */
  int verbose=0;
  int i=0;
  double oldstep=0,newstep=0,nexttime=st->last_time+st->update_timestep;
  int oldtab=0,newtab=0,tab2=0,tab3=0,nt=0;
  double *temp=NULL,*temp2=NULL,tempmean=0;
  if (t+DT < nexttime){ 
    st->data[st->tab] += val*DT/st->update_timestep;}
  else /* if (t+DT >= nexttime) */{ 
    oldstep = minimum(st->update_timestep,maximum(0,nexttime-t)); newstep = minimum(st->update_timestep,maximum(0,t+DT-nexttime));
    oldtab = st->tab; 
    st->data[oldtab] += val*oldstep/st->update_timestep;
    newtab = st->tab+1; 
    if (newtab==st->length){ 
      newtab=0; st->cyclenum += 1;
      if (st->cycle_bother){ for (i=0;i<st->length;i++){ st->cycledata[i] += st->data[i];}}}
    st->data[newtab] = val*newstep/st->update_timestep;
    if (st->lpower_bother && (newtab%st->lpower_update_every==0)){
      tab3 = oldtab/st->lpower_update_every * st->lpower_window_length;
      if (verbose){ printf(" %% updating lpowerdata[%d]",tab3);}
      if (oldtab-(st->lpower_window_length-1)>=0){
	tab2 = oldtab-(st->lpower_window_length-1);
	if (verbose){ printf(" contiguous array, oldtab %d,tab2 %d",oldtab,tab2);}
	if (st->cycle_bother){
	  temp = ra2power(NULL,&(st->data[tab2]),st->lpower_window_length,NULL,0,1);
	  raplusequals(&(st->lpowerdata[tab3]),st->lpower_window_length,temp);
	  tfree(temp);temp=NULL;}
	else /* if (!st->cycle_bother) */{ 
	  ra2power(NULL,&(st->data[tab2]),st->lpower_window_length,&(st->lpowerdata[tab3]),0,1);
	  if (verbose){ raprintf(&(st->lpowerdata[tab3]),"double",1,st->lpower_window_length,"lpowerdata");}}}
      else if (oldtab-(st->lpower_window_length-1)<0){
	tab2 = periodize(oldtab-(st->lpower_window_length-1),0,st->length);
	if (verbose){ printf(" discontiguous array, oldtab %d,tab2 %d",oldtab,tab2);}
	temp2 = (double *) tcalloc(st->lpower_window_length,sizeof(double));
	if (st->cyclenum>0){
	  for (nt=0;nt<st->lpower_window_length;nt++){ temp2[nt] = st->data[periodize(tab2+nt,0,st->length)];}}
	else /* if (st->cyclenum==0) */{
	  stats("double",st->data,oldtab,NULL,NULL,&tempmean,NULL);
	  for (nt=0;nt<=st->lpower_window_length-1-oldtab;nt++){ temp2[nt] = tempmean;}
	  for (nt=st->lpower_window_length-oldtab;nt<st->lpower_window_length;nt++){ 
	    temp2[nt] = st->data[periodize(tab2+nt,0,st->length)];}}
	if (verbose){ raprintf(temp2,"double",1,st->lpower_window_length,"temp");}
	if (st->cycle_bother){
	  temp = ra2power(NULL,temp2,st->lpower_window_length,NULL,0,1);
	  raplusequals(&(st->lpowerdata[tab3]),st->lpower_window_length,temp);
	  tfree(temp);temp=NULL;}
	else /* if (!st->cycle_bother) */{ 
	  ra2power(NULL,temp2,st->lpower_window_length,&(st->lpowerdata[tab3]),0,1);}
	tfree(temp2);temp2=NULL;}
      if (verbose){ raprintf(&(st->lpowerdata[tab3]),"double",1,st->lpower_window_length,"lpowerdata");}}
    st->tab = newtab;
    st->last_time = maximum(t,st->last_time + st->update_timestep);}
}

void strobeupdate_bkp(struct strobe *st,double t,double DT,double val)
{
  /* this averages, rather than strobes, which is better for spike statistics */
  int i=0;
  st->data[st->tab] += val*DT/st->update_timestep;
  if (t+DT >= st->last_time+st->update_timestep){
    st->tab++;
    if (st->tab >= st->length){ 
      st->tab=0; st->cyclenum += 1;
      if (st->cycle_bother){ for (i=0;i<st->length;i++){ st->cycledata[i] += st->data[i];}}}
    st->data[st->tab] = 0;
    st->last_time = maximum(t,st->last_time + st->update_timestep);}
}

void strobeupdate_old(struct strobe *st,double t,double val)
{
  /* this strobes, as opposed to averaging, which is not good for discontinuous data like spike statistics */
  int i=0;
  if (t >= st->last_time+st->update_timestep){
    st->data[st->tab] = val; st->tab++;
    if (st->tab >= st->length){ 
      st->tab=0; st->cyclenum += 1;
      if (st->cycle_bother){ for (i=0;i<st->length;i++){ st->cycledata[i] += st->data[i];}}}
    st->last_time = maximum(t,st->last_time + st->update_timestep);}
}

void stradump(struct strobe **stra,int ralength,int dump_type,char *filename_base)
{
  /* assumes filename_base does NOT start with "./" */
  int wrap_at=1024,rows=0,cols=0;
  char filename_base2[256],filename2[512];
  int sttab=stra[0]->tab,stlength=stra[0]->length,stcycle_bother=stra[0]->cycle_bother,stcyclenum=stra[0]->cyclenum,stpower_bother=stra[0]->lpower_bother;
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
      if (fp!=stdout){ fclose(fp);}}
    if (stpower_bother){
      sprintf(filename2,"%s_power.m",filename_base);
      if ((fp=fopen(filename2,"w"))==NULL){ printf(" %% Warning, couldn't open %s in stradump, writing to stdout\n",filename2); fp=stdout;}
      fprintf(fp,"clear %s_power_ra;\n",filename_base);
      for (na=0;na<ralength;na++){
	fprintf(fp,"%% array %d of %d;\n",na+1,ralength);
	for (nt=0;nt<stra[na]->lpower_window_length;nt++){
	  for (nt2=0;nt2<stra[na]->lpower_length/stra[na]->lpower_window_length;nt2++){
	    fprintf(fp,"%s_power_ra(%d,%d,%d)=%0.16f;\n",filename_base,na+1,nt+1,nt2+1,stra[na]->lpowerdata[nt+nt2*stra[na]->lpower_window_length]);}}}
      fprintf(fp,"%% now plotting;\n");
      for (nt=0;nt<ralength;nt++){
	fprintf(fp,"figure;clf;hold on;imagesc(%s_power_ra(%d,:,:));title('%s_power_ra_%d');hold off;\n",filename_base,nt+1,filename_base,nt+1);}
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
      sprintf(filename_base2,"./%s_cycle",filename_base);
      if (wrap_at>0 && stlength>wrap_at){
	rows = (ralength+1)*((int)ceil((double)stlength/(double)wrap_at));
	cols = wrap_at;
	ra = (double *) tcalloc(rows*cols,sizeof(double));
	for (na=0;na<ralength;na++){ for (nt=0;nt<stlength;nt++){
	  ra[na + ((nt)/wrap_at)*(ralength+1) + ((nt)%wrap_at)*rows] = stra[na]->cycledata[nt]/stcyclenum;}}
	if (ralength>1){ sprintf(filename2,"%s.pnm",filename_base2); WritePNMfile_color(ra,rows,cols,0,0,filename2,7);}
	ra2jpg(ra,"double",rows,cols,0,filename_base2,0);
	tfree(ra);}
      else{ 
	ra = (double *) tcalloc(ralength*stlength,sizeof(double));
	for (na=0;na<ralength;na++){ for (nt=0;nt<stlength;nt++){
	  ra[na + (nt)*ralength] = stra[na]->cycledata[nt]/stcyclenum;}}
	if (ralength>1){ sprintf(filename2,"%s.pnm",filename_base2); WritePNMfile_color(ra,ralength,stlength,0,0,filename2,7);}
	ra2jpg(ra,"double",ralength,stlength,0,filename_base2,0);
	tfree(ra);}}
    if (stpower_bother){
      for (na=0;na<ralength;na++){
	sprintf(filename_base2,"./%s_lpower%d",filename_base,na);
	sprintf(filename2,"%s.pnm",filename_base2);
	WritePNMfile_color(stra[na]->lpowerdata,stra[na]->lpower_window_length,stra[na]->lpower_length/stra[na]->lpower_window_length,0,0,filename2,7); 
	ra2jpg(stra[na]->lpowerdata,"double",stra[na]->lpower_window_length,stra[na]->lpower_length/stra[na]->lpower_window_length,0,filename_base2,3);}}}
}

void strobetfree(struct strobe *st){ 
  if (st->data!=NULL){ tfree(st->data); st->data=NULL;} 
  if (st->cycledata!=NULL){ tfree(st->cycledata); st->cycledata=NULL;} 
  if (st->lpowerdata!=NULL){ tfree(st->lpowerdata); st->lpowerdata=NULL;}
  tfree(st); st=NULL;}

void strobereset(struct strobe *st)
{
  st->tab = 0;
  st->last_time = GLOBAL_TI;
  rareset(st->data,"double",st->length,NULL);
  if (st->cycle_bother){ rareset(st->cycledata,"double",st->length,NULL);}
  if (st->lpower_bother){ rareset(st->lpowerdata,"double",st->lpower_length,NULL);}
}

double * stra2ra(struct strobe **stra,int ralength,int window_length,int cycle_flag,double *output)
{
  int verbose=0;
  int length=stra[0]->length/window_length,nr=0,nt=0,tab=0;
  double *datara=NULL;
  if (verbose){ printf(" %% [entering stra2ra] ralength %d window_length %d\n",ralength,window_length);}
  if (verbose){ printf(" %% simply sums over window_length %f\n",window_length*stra[0]->update_timestep);}
  if (output!=NULL){ datara=output; rareset(datara,"double",ralength*length,NULL);}
  else /* if (output==NULL) */{ datara = (double *) tcalloc(ralength*length,sizeof(double));}
  for (nr=0;nr<ralength;nr++){ for (nt=0;nt<stra[nr]->length;nt++){ 
    if (cycle_flag){ datara[nt/window_length+nr*length] += stra[nr]->cycledata[nt];}
    else /* if (!cycle_flag) */{ 
      tab=periodize(stra[nr]->tab+nt,0,stra[nr]->length);
      datara[nt/window_length+nr*length] += stra[nr]->data[tab];}}}
  if (verbose){ raprintf(datara,"double",length,ralength,"datara");}
  return datara;
}

double * ra2pca(double *ra,int length,int ntracks,double *output)
{
  int verbose=0;
  int nr=0,nt=0,nr2=0,tab=0;
  double *meanra=NULL,*datara=NULL,*dataproject=NULL;
  doublereal *cormat=NULL,*eigval=NULL,*work=NULL;
  char uplo='U',jobz='V';
  integer info=0,lcormat=ntracks,lwork=3*ntracks;
  int ndimensions=3,nproject=minimum(ndimensions,ntracks);
  if (verbose){ printf(" %% [entering ra2pca] length %d ntracks %d\n",length,ntracks);}
  if (verbose){ printf(" %% ra stored as ra[nt + nr*length]\n");}
  if (verbose){ raprintf(ra,"double",length,ntracks,"ra");}
  meanra = (double *) tcalloc(ntracks,sizeof(double));
  for (nr=0;nr<ntracks;nr++){ stats("double",&(ra[0+nr*length]),length,NULL,NULL,&(meanra[nr]),NULL);}
  if (verbose){ raprintf(meanra,"double",1,ntracks,"meanra");}
  datara = (double *) tcalloc(length*ntracks,sizeof(double));
  for (nr=0;nr<ntracks;nr++){ for (nt=0;nt<length;nt++){ 
    datara[nt+nr*length] += ra[nt+nr*length]-meanra[nr];}}
  if (verbose){ raprintf(datara,"double",length,ntracks,"datara");}
  cormat = (doublereal *) tcalloc(ntracks*ntracks,sizeof(doublereal));
  for (nr=0;nr<ntracks;nr++){ for (nr2=nr;nr2<ntracks;nr2++){
    tab = nr+nr2*ntracks;
    cormat[tab]=0; for (nt=0;nt<length;nt++){ cormat[tab] += datara[nt+nr*length]*datara[nt+nr2*length];} 
    cormat[tab] /= maximum(1,length-1);
    cormat[nr2+nr*ntracks]=cormat[tab];
    if (verbose){ printf(" covariance %d,%d = %f\n",nr,nr2,cormat[tab]);}}}
  eigval = (doublereal *) tcalloc(ntracks,sizeof(doublereal));
  work = (doublereal *) tcalloc(ntracks*3,sizeof(doublereal));
  dsyev_(&jobz,&uplo,&lcormat,cormat,&lcormat,eigval,work,&lwork,&info);
  if (verbose){ raprintf(eigval,"double",1,ntracks,"eigval"); raprintf(cormat,"double",1,ntracks*ntracks,"eigvec");}
  if (info!=0){ printf("warning! info=%d in stra2pca\n",(int)info);}
  if (output!=NULL){ dataproject=output; rareset(dataproject,"double",length*ndimensions,NULL);}
  else /* if (output==NULL) */{ dataproject = (double *) tcalloc(length*ndimensions,sizeof(double));}
  for (nr=0;nr<nproject;nr++){ for (nt=0;nt<length;nt++){ for (nr2=0;nr2<ntracks;nr2++){
    dataproject[nt + nr*length] += datara[nt+nr2*length]*(double)cormat[(ntracks-(nr+1))*ntracks+nr2];}}}
  if (verbose){ raprintf(dataproject,"double",length,ndimensions,"dataproject");}
  tfree(meanra);meanra=NULL;
  tfree(datara);datara=NULL;
  tfree(cormat);cormat=NULL;
  tfree(eigval);eigval=NULL;
  tfree(work);work=NULL;
  return dataproject;
}

void pca2jpg(double *pcara,int length,char *filename_base)
{
  /* plots ndimensions=3 columns of pcara (each length elements) */
  /* filename_base should start with "./" */
  int ndimensions=3,nsides=8;
  int backarrow_flag=1,dot_flag=0,use_stdev=0,remove_flag=0,jpeg_flag=(ON_MY_COMPUTER ? 1 : 0);
  double maxdia = 10000;
  FILE *fp=NULL;
  char filename[512];
  double *maxra=NULL,*minra=NULL,*meanra=NULL,*stdevra=NULL;
  int nr=0,nt=0,nl;
  int colorcode=0;
  double rcolor=0,gcolor=0,bcolor=0;
  char command[1024];
  double xpos=0,ypos=0,zpos=0,rord=0;
  maxra = (double *) tcalloc(ndimensions,sizeof(double));
  minra = (double *) tcalloc(ndimensions,sizeof(double));
  meanra = (double *) tcalloc(ndimensions,sizeof(double));
  stdevra = (double *) tcalloc(ndimensions,sizeof(double));
  for (nr=0;nr<ndimensions;nr++){ 
    stats("double",&(pcara[0+nr*length]),length,&(maxra[nr]),&(minra[nr]),&(meanra[nr]),&(stdevra[nr]));
    if (use_stdev){ maxra[nr]=meanra[nr]+STD_VIEW*stdevra[nr]; minra[nr]=meanra[nr]-STD_VIEW*stdevra[nr];}}
  sprintf(filename,"%s.fig",filename_base);
  if ((fp=fopen(filename,"w"))==NULL){ printf(" %% Warning, cannot create %s in pca2jpg\n",filename); fp=stdout;}
  fprintf(fp,"%s",FIG_PREAMBLE); fprintf(fp,"%s",FIG_PREAMBLE_COLOR_7);
  for (nt=0;nt<length;nt++){ 
    xpos = (pcara[nt+0*length]-minra[0])/(maxra[0]-minra[0]);
    ypos = (pcara[nt+1*length]-minra[1])/(maxra[1]-minra[1]);
    zpos = (pcara[nt+2*length]-minra[2])/(maxra[2]-minra[2]);
    rord = 0.01+0.04*zpos;
    colorscale(0,zpos,1,0,&rcolor,&gcolor,&bcolor);colorcode=crop((int)floor(512*rcolor),0,511);
    if (dot_flag){
      fprintf(fp,"2 1 0 0 0 %d %d 0 %d 0.000 0 0 -1 0 0 %d\n",/*colorcode*/colorcode+32,/*depth*/(int)(1+998*(1-zpos)),/*fill*/20,/*npoints*/nsides+1); fprintf(fp,"\t"); for (nl=0;nl<=nsides;nl++){ fprintf(fp,"%d %d ",(int)floor(maxdia*(xpos+rord*cos(2*PI*((double)nl+0.5)/(double)nsides))),(int)maxdia-(int)floor(maxdia*(ypos+rord*sin(2*PI*((double)nl+0.5)/(double)nsides))));} fprintf(fp,"\n");}
    if (nt>0){
      fprintf(fp,"2 1 0 5 %d 7 %d 0 -1 0.000 0 0 -1 0 %d %d\n",/*color*/colorcode+32,/*depth*/(int)(1+998*(1-zpos)),/*backarrow*/backarrow_flag,/*npoints*/2);
      if (backarrow_flag){ fprintf(fp,"\t0 0 5.00 300.00 600.00\n\t");}
      fprintf(fp," %d %d",(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos));
      xpos = (pcara[(nt-1)+0*length]-minra[0])/(maxra[0]-minra[0]); ypos = (pcara[(nt-1)+1*length]-minra[1])/(maxra[1]-minra[1]);
      fprintf(fp," %d %d",(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos));
      fprintf(fp,"\n");}}
  fprintf(fp,"4 0 0 %d 0 %d %d 0.0000 4 %d %d %d %d x=%0.2e\\001\n",/*depth*/1,/*font*/12,/*point*/24,/*textheight*/(int)(24*13.5),/*textwidth*/(int)(24*9*18),/*xpos*/(int)(0.1*maxdia),/*ypos*/(int)maxdia-(int)(1.0*maxdia),maxra[0]);
  fprintf(fp,"4 0 0 %d 0 %d %d 0.0000 4 %d %d %d %d y=%0.2e\\001\n",/*depth*/1,/*font*/12,/*point*/24,/*textheight*/(int)(24*13.5),/*textwidth*/(int)(24*9*18),/*xpos*/(int)(0.1*maxdia),/*ypos*/(int)maxdia-(int)(0.9*maxdia),maxra[1]);
  fprintf(fp,"4 0 0 %d 0 %d %d 0.0000 4 %d %d %d %d z=%0.2e\\001\n",/*depth*/1,/*font*/12,/*point*/24,/*textheight*/(int)(24*13.5),/*textwidth*/(int)(24*9*18),/*xpos*/(int)(0.1*maxdia),/*ypos*/(int)maxdia-(int)(0.8*maxdia),maxra[2]);
  fprintf(fp,"4 0 0 %d 0 %d %d 0.0000 4 %d %d %d %d x=%0.2e\\001\n",/*depth*/1,/*font*/12,/*point*/24,/*textheight*/(int)(24*13.5),/*textwidth*/(int)(24*9*18),/*xpos*/(int)(0.1*maxdia),/*ypos*/(int)maxdia-(int)(0.0*maxdia),minra[0]);
  fprintf(fp,"4 0 0 %d 0 %d %d 0.0000 4 %d %d %d %d y=%0.2e\\001\n",/*depth*/1,/*font*/12,/*point*/24,/*textheight*/(int)(24*13.5),/*textwidth*/(int)(24*9*18),/*xpos*/(int)(0.1*maxdia),/*ypos*/(int)maxdia-(int)(0.1*maxdia),minra[1]);
  fprintf(fp,"4 0 0 %d 0 %d %d 0.0000 4 %d %d %d %d z=%0.2e\\001\n",/*depth*/1,/*font*/12,/*point*/24,/*textheight*/(int)(24*13.5),/*textwidth*/(int)(24*9*18),/*xpos*/(int)(0.1*maxdia),/*ypos*/(int)maxdia-(int)(0.2*maxdia),minra[2]);
  if (fp!=stdout){ fclose(fp);}
  if (jpeg_flag){ sprintf(command,"fig2dev -Ljpeg -q %d %s.fig %s.jpg;",/*quality*/5,filename_base,filename_base); system(command);}
  if (remove_flag){ sprintf(command,"rm %s.fig;",filename_base); system(command);}
}

/* Here are the power functions */

struct power * powermake(struct neuronarray *Nra,int length,int indexing_ntype_length,int *indexing_ntype_checkout,int *indexing_ntype_refile,int indexing_nvar_length,int *indexing_nvar_checkout,int *indexing_nvar_refile,double *maxra,double *minra,int correlation_bother)
{
  int verbose=0;
  int nr=0,nt=0,nv=0;
  struct power *p=NULL;
  if (verbose){ 
    printf(" %% [entering powermake]\n");
    raprintf(indexing_ntype_checkout,"int",1,indexing_ntype_length,"ntype_checkout");
    raprintf(indexing_ntype_refile,"int",1,Nra->ntypes,"ntype_refile");
    raprintf(indexing_nvar_checkout,"int",1,indexing_nvar_length,"nvar_checkout");
    raprintf(indexing_nvar_refile,"int",1,Nra->nvars,"nvar_refile");
    raprintf(maxra,"double",1,Nra->nvars,"maxra");
    raprintf(minra,"double",1,Nra->nvars,"minra");}
  p = (struct power *) tmalloc(sizeof(struct power));
  p->update_every=1;
  p->length=length;
  p->Nra=Nra;
  p->indexing_nvar_length = indexing_nvar_length;
  p->indexing_nvar_checkout = (int *) tcalloc(p->indexing_nvar_length,sizeof(int));
  for (nv=0;nv<p->indexing_nvar_length;nv++){ p->indexing_nvar_checkout[nv] = indexing_nvar_checkout[nv];}
  p->indexing_nvar_refile = (int *) tcalloc(p->Nra->nvars,sizeof(int));
  for (nv=0;nv<p->Nra->nvars;nv++){ p->indexing_nvar_refile[nv] = indexing_nvar_refile[nv];}
  p->indexing_ntype_length = indexing_ntype_length;
  p->indexing_ntype_checkout = (int *) tcalloc(p->indexing_ntype_length,sizeof(int));
  for (nt=0;nt<p->indexing_ntype_length;nt++){ p->indexing_ntype_checkout[nt] = indexing_ntype_checkout[nt];}
  p->indexing_ntype_refile = (int *) tcalloc(p->Nra->ntypes,sizeof(int));
  for (nt=0;nt<p->Nra->ntypes;nt++){ p->indexing_ntype_refile[nt] = indexing_ntype_refile[nt];}
  p->tst = strobemake(p->length,p->update_every,0,0,0,0);
  p->strarara = (struct strobe ****) tcalloc(p->indexing_ntype_length,sizeof(struct strobe ***));
  p->vpowrarara = (double ***) tcalloc(p->indexing_ntype_length,sizeof(double **));
  for (nt=0;nt<p->indexing_ntype_length;nt++){
    p->strarara[nt] = (struct strobe ***) tcalloc(p->indexing_nvar_length,sizeof(struct strobe **));
    p->vpowrarara[nt] = (double **) tcalloc(p->indexing_nvar_length,sizeof(double *));
    for (nv=0;nv<p->indexing_nvar_length;nv++){
      if (p->Nra->lengthra[p->indexing_ntype_checkout[nt]]!=0){
	p->strarara[nt][nv] = (struct strobe **) tcalloc(p->Nra->lengthra[p->indexing_ntype_checkout[nt]],sizeof(struct strobe *));
	for (nr=0;nr<p->Nra->lengthra[p->indexing_ntype_checkout[nt]];nr++){
	  p->strarara[nt][nv][nr] = strobemake(p->length,p->update_every,0,0,0,0);}
	p->vpowrarara[nt][nv] = (double *) tcalloc(p->length,sizeof(double));}
      else /* if (p->Nra->lengthra[p->indexing_ntype_checkout[nt]]==0) */{ p->strarara[nt][nv] = NULL; p->vpowrarara[nt][nv] = NULL;}}}
  p->firstvarpower = (double *) tcalloc(p->length,sizeof(double));
  p->correlation_bother = correlation_bother;
  if (p->correlation_bother){ 
    p->autocorrelation = (double *) tcalloc(p->length,sizeof(double));
    p->crosscorrelation = (double *) tcalloc(p->length,sizeof(double));}
  p->maxra= (double *) tcalloc(p->Nra->nvars,sizeof(double));
  p->minra= (double *) tcalloc(p->Nra->nvars,sizeof(double));
  for (nv=0;nv<Nra->nvars;nv++){ p->maxra[nv] = maxra[nv]; p->minra[nv] = minra[nv];}
  p->update_number=0;
  return p;
}

void powertfree(struct power *p)
{
  int verbose=0;
  int nr=0,nt=0,nv=0;
  if (verbose){ printf(" %% [entering powertfree]\n");}
  strobetfree(p->tst);p->tst=NULL;
  for (nt=0;nt<p->indexing_ntype_length;nt++){
    for (nv=0;nv<p->indexing_nvar_length;nv++){
      if (p->strarara[nt][nv]!=NULL){
	for (nr=0;nr<p->Nra->lengthra[p->indexing_ntype_checkout[nt]];nr++){
	  strobetfree(p->strarara[nt][nv][nr]); p->strarara[nt][nv][nr]=NULL;}
	tfree(p->strarara[nt][nv]);p->strarara[nt][nv]=NULL;}
      if (p->vpowrarara[nt][nv]!=NULL){
	tfree(p->vpowrarara[nt][nv]);p->vpowrarara[nt][nv]=NULL;}}
    tfree(p->strarara[nt]);p->strarara[nt]=NULL;
    tfree(p->vpowrarara[nt]);p->vpowrarara[nt]=NULL;}
  tfree(p->indexing_nvar_checkout); p->indexing_nvar_checkout=NULL;
  tfree(p->indexing_nvar_refile); p->indexing_nvar_refile=NULL;
  tfree(p->indexing_ntype_checkout); p->indexing_ntype_checkout=NULL;
  tfree(p->indexing_ntype_refile); p->indexing_ntype_refile=NULL;
  tfree(p->strarara);p->strarara=NULL;
  tfree(p->vpowrarara);p->vpowrarara=NULL;
  tfree(p->firstvarpower);p->firstvarpower=NULL;
  if (p->correlation_bother){
    if (p->autocorrelation!=NULL){ tfree(p->autocorrelation);p->autocorrelation=NULL;}
    if (p->crosscorrelation!=NULL){ tfree(p->crosscorrelation);p->crosscorrelation=NULL;}}
  tfree(p->maxra);p->maxra=NULL;
  tfree(p->minra);p->minra=NULL;
  tfree(p);p=NULL;
  if (verbose){ printf(" %% finished\n");}
}

void powerupdate(struct power *p,double t,double DT)
{
  int verbose=0;
  int nr=0,nt=0,nv=0,nl=0;
  struct neuron *n=NULL;
  int tab_old=0,tab_new=0,tab=0;
  double *temp=NULL;
  double *ra=NULL;
  double **rara=NULL;
  int lt=0;
  int varname_registry_lfp=0;
  int varname_registry_spike_flag=0;
  if (verbose){ printf(" %% [entering powerupdate] t %f dt %f\n",t,DT);}
  switch (GLOBAL_NEURON_MODEL){
  case 0: case 1: case 2: case 3: case 4: varname_registry_spike_flag=VARNAME_REGISTRY_spike_flag; break;
  case 5: /* mainak style */ varname_registry_spike_flag=VARNAME_REGISTRY_mainak_spike_flag; break;
  case 6: /* wilson style */ varname_registry_spike_flag=VARNAME_REGISTRY_wilson_spike_flag; break;
  case 7: /* snx style */ varname_registry_spike_flag=VARNAME_REGISTRY_snx_spike_flag; break;
  default: break;}
  tab_old = p->tst->tab; strobeupdate_old(p->tst,t,t); tab_new = p->tst->tab;
  for (nt=0;nt<p->indexing_ntype_length;nt++){
    if (verbose>1){ printf(" %% nt %d (%s)\n",nt,GLOBAL_TYPENAMES[p->indexing_ntype_checkout[nt]]);}
    for (nr=0;nr<p->Nra->lengthra[p->indexing_ntype_checkout[nt]];nr++){
      if (verbose>1){ printf(" %% nr %d\n",nr);}
      n = nget(p->Nra,p->indexing_ntype_checkout[nt],nr); 
      for (nv=0;nv<p->indexing_nvar_length;nv++){
	if (verbose>1){ printf(" %% nv %d (%s)\n",nv,GLOBAL_VARNAMES[p->indexing_nvar_checkout[nv]]);}
	if (p->strarara[nt][nv]!=NULL){
	  if (verbose>1){ printf(" %% updating p->strarara[nt %d %s][nv %d %s][nr %d] with %f...\n",nt,GLOBAL_TYPENAMES[p->indexing_ntype_checkout[nt]],nv,GLOBAL_VARNAMES[p->indexing_nvar_checkout[nv]],nr,*(n->vpra[p->indexing_nvar_checkout[nv]]));}
	  if (p->indexing_nvar_checkout[nv]==varname_registry_spike_flag){
	    strobeupdate_old(p->strarara[nt][nv][nr],t,*(n->vpra[p->indexing_nvar_checkout[nv]]));}
	  else /* if (p->indexing_nvar_checkout[nv]!=varname_registry_spike_flag) */{
	    strobeupdate(p->strarara[nt][nv][nr],t,DT,*(n->vpra[p->indexing_nvar_checkout[nv]]));}
	  /* fix later */
	  /* 	  strobeupdate(p->strarara[nt][nv][nr],t,DT,*(n->vpra[p->indexing_nvar_checkout[nv]])); */
	  if (verbose>1){ printf(" %% finished\n");}}}}}
  if (tab_old==p->length-1 && tab_new==0){ 
    if (verbose){ printf(" %% time to update\n");}
    ra = (double *) tcalloc(p->length,sizeof(double));
    lt=0;for (nt=0;nt<p->indexing_ntype_length;nt++){ lt+=p->Nra->lengthra[p->indexing_ntype_checkout[nt]];}
    rara = (double **) tcalloc(lt,sizeof(double *));
    tab=0;
    for (nt=0;nt<p->indexing_ntype_length;nt++){
      for (nr=0;nr<p->Nra->lengthra[p->indexing_ntype_checkout[nt]];nr++){
	for (nv=0;nv<p->indexing_nvar_length;nv++){
	  if (p->strarara[nt][nv]!=NULL && p->vpowrarara[nt][nv]!=NULL){
	    temp = ra2power(NULL,p->strarara[nt][nv][nr]->data,p->length,NULL,0,0); 
	    raplusequals(p->vpowrarara[nt][nv],p->length,temp); 
	    tfree(temp); temp=NULL;}}
	rara[tab+nr] = (double *) tcalloc(p->length,sizeof(double));
	if (p->strarara[nt][varname_registry_lfp]!=NULL){
	  for (nl=0;nl<p->length;nl++){ 
	    rara[tab+nr][nl] = p->strarara[nt][varname_registry_lfp][nr]->data[nl]; 
	    ra[nl] += p->strarara[nt][varname_registry_lfp][nr]->data[nl];}}
	p->update_number += 1;}
      tab += p->Nra->lengthra[p->indexing_ntype_checkout[nt]];}
    if (p->correlation_bother){ 
      rara2corr(NULL,NULL,rara,lt,p->length,p->autocorrelation,p->crosscorrelation);
      if (verbose){ raprintf(rara[0],"double",1,p->length,"rara: "); raprintf(p->autocorrelation,"double",1,p->length,"autocorrelation: "); raprintf(p->crosscorrelation,"double",1,p->length,"crosscorrelation: "); }}
    temp = ra2power(NULL,ra,p->length,NULL,0,0); raplusequals(p->firstvarpower,p->length,temp); tfree(temp); temp=NULL;
    tfree(ra);ra=NULL;
    for (nr=0;nr<lt;nr++){ if (rara[nr]!=NULL){ tfree(rara[nr]);rara[nr]=NULL;}}
    tfree(rara);rara=NULL;}
  if (verbose){ printf(" %% [finishing powerupdate] \n");}
}

void powerdump_trial_average(struct power *p,char *fgvn,int ntrials)
{
  /* presumes one very long multi-trial run of length ntrials*(p->length/ntrials), 
     only useful if storage is not an issue (i.e., small system). 
     produces a set of rasters (as well as plots) where trial and time have been averaged */
  char filename[1024],filename_base[1024],command[2048];
  char gs2[512];
  FILE *fp=NULL;
  int nt=0,nv=0,nr=0,nl=0,nl2=0,nl3=0,nrecords=0,sttab=0,stlength=0,tab=0;
  double *ra=NULL,*ra2=NULL,max=0,min=0,mean=0,stdev=0,max2=0,min2=0;
  int trial_length = p->length/ntrials;
  int lcolor=0,lcolormin=0,lcolormax=minimum(30,ntrials),stdev_use=0;
  double tolerance=0.000001,xpos=0,ypos=0,maxdia = 10000;
  int remove_flag=1,jpeg_flag=(ON_MY_COMPUTER ? 1 : 0);
  int time_binsize=16,time_nbins = trial_length/time_binsize;
  int nglom_nbins = GLOBAL_NCLUSTERS,nglom_binsize = 0;
  if (fgvn==NULL){ sprintf(gs2,"%srecord",GLOBAL_STRING_2);}
  else{ sprintf(gs2,"%srecord_%s",GLOBAL_STRING_2,fgvn);}
  for (nt=0;nt<p->indexing_ntype_length;nt++){ 
    nrecords = p->Nra->lengthra[p->indexing_ntype_checkout[nt]];
    nglom_binsize = nrecords/nglom_nbins;
    for (nv=0;nv<p->indexing_nvar_length;nv++){
      if (p->strarara[nt][nv]!=NULL){
	sttab = p->strarara[nt][nv][0]->tab; stlength = p->strarara[nt][nv][0]->length;
	ra = (double *) tcalloc(time_nbins*nrecords,sizeof(double));
	for (nr=0;nr<nrecords;nr++){ for (nl=0;nl<time_nbins;nl++){
	  ra[nl + nr*time_nbins] = 0;
	  for (nl2=0;nl2<time_binsize;nl2++){ for (nl3=0;nl3<ntrials;nl3++){
	    tab = nl2 + nl*time_binsize + nl3*trial_length;
	    tab = periodize(sttab-trial_length*ntrials + tab,0,stlength);
	    ra[nl + nr*time_nbins] += p->strarara[nt][nv][nr]->data[tab];}}
	  ra[nl + nr*time_nbins] /= (double)(time_binsize*ntrials);}}
	stats("double",ra,time_nbins*nrecords,&max2,&min2,&mean,&stdev);
	if (max2-min2<=tolerance){ /* do nothing */ }
	else /* if (max2>min2+tolerance) */{ 
	  max = p->maxra[p->indexing_nvar_checkout[nv]]; min = p->minra[p->indexing_nvar_checkout[nv]];
	  if (max<=min){ if (stdev_use){ max=mean+STD_VIEW*stdev; min=mean-STD_VIEW*stdev;} else if (!stdev_use){ max=max2;min=min2;}}
	  sprintf(filename,"./power_draw_trial_average_%s_%s_%s.pnm",gs2,GLOBAL_TYPENAMES[p->indexing_ntype_checkout[nt]],GLOBAL_VARNAMES[p->indexing_nvar_checkout[nv]]);
	  WritePNMfile_color(ra,time_nbins,nrecords,max,min,filename,7);
	  ra2 = (double *) tcalloc(time_nbins*nglom_nbins,sizeof(double));
	  for (nl=0;nl<time_nbins;nl++){ for (nl2=0;nl2<nglom_nbins;nl2++){
	    ra2[nl + nl2*time_nbins]=0;
	    for (nl3=0;nl3<nglom_binsize;nl3++){
	      ra2[nl + nl2*time_nbins] += ra[nl + (nl3 + nl2*nglom_binsize)*time_nbins];}
	    ra2[nl + nl2*time_nbins] /= (double)nglom_binsize;}}
	  sprintf(filename_base,"./power_plot_trial_average_%s_%s_%s",gs2,GLOBAL_TYPENAMES[p->indexing_ntype_checkout[nt]],GLOBAL_VARNAMES[p->indexing_nvar_checkout[nv]]);
	  sprintf(filename,"%s.fig",filename_base);
	  if ((fp=fopen(filename,"w"))==NULL){ printf(" %% Warning, cannot create %s in ra2jpg\n",filename); fp=stdout;}
	  fprintf(fp,"%s",FIG_PREAMBLE);
	  for (nr=0;nr<nglom_nbins;nr++){
	    periodify("int",&nr,&lcolormin,&lcolormax,&lcolor);
	    fprintf(fp,"2 1 0 %d %d 7 %d 0 -1 0.000 0 0 -1 0 0 %d\n\t",/*width*/1,/*color*/lcolor,/*depth*/(nr%999)+1,/*npoints*/time_nbins);
	    for (nl2=0;nl2<time_nbins;nl2++){
	      xpos = ((double)nl2+0.5)/(double)time_nbins;
	      ypos = (ra2[nl2 + nr*time_nbins] - min)/(max-min)/(double)nglom_nbins + (double)nr/(double)nglom_nbins;
	      fprintf(fp," %d %d",(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos));}
	    fprintf(fp,"\n");}
	  tfree(ra2);ra2=NULL;
	  if (fp!=stdout){ fclose(fp);}
	  if (jpeg_flag){ sprintf(command,"fig2dev -Ljpeg -q %d %s.fig %s.jpg;",/*quality*/5,filename_base,filename_base); system(command);}
	  if (remove_flag){ sprintf(command,"rm %s.fig;",filename_base); system(command);}}
	tfree(ra);ra=NULL;}}}
}

void powerdump_trial_interlace(struct power *p,char *fgvn,int ntrials)
{
  /* presumes one very long multi-trial run of length ntrials*(p->length/ntrials), 
     only useful if storage is not an issue (i.e., small system). 
     produces a set of rasters (as well as plots) where multiple trials of a single neuron are adjacent. */
  char filename[1024],filename_base[1024],command[2048];
  char gs2[512];
  FILE *fp=NULL;
  int nt=0,nv=0,nr=0,nl=0,nl2=0,nrecords=0,sttab=0,stlength=0;
  double *ra=NULL,max=0,min=0,mean=0,stdev=0,max2=0,min2=0;
  int trial_length = p->length/ntrials;
  int lcolor=0,lcolormin=0,lcolormax=minimum(30,ntrials),stdev_use=0;
  double tolerance=0.000001,xpos=0,ypos=0,maxdia = 10000;
  int remove_flag=1,jpeg_flag=(ON_MY_COMPUTER ? 1 : 0);
  if (fgvn==NULL){ sprintf(gs2,"%srecord",GLOBAL_STRING_2);}
  else{ sprintf(gs2,"%srecord_%s",GLOBAL_STRING_2,fgvn);}
  for (nt=0;nt<p->indexing_ntype_length;nt++){ 
    nrecords = p->Nra->lengthra[p->indexing_ntype_checkout[nt]];
    for (nv=0;nv<p->indexing_nvar_length;nv++){
    if (p->strarara[nt][nv]!=NULL){
      sttab = p->strarara[nt][nv][0]->tab; stlength = p->strarara[nt][nv][0]->length;
      ra = (double *) tcalloc(trial_length*ntrials*nrecords,sizeof(double));
      for (nr=0;nr<nrecords;nr++){ for (nl=0;nl<trial_length*ntrials;nl++){
	nl2 = periodize(sttab-trial_length*ntrials + nl,0,stlength);
	ra[nl + nr*trial_length*ntrials] = p->strarara[nt][nv][nr]->data[nl2];}}
      stats("double",ra,trial_length*ntrials*nrecords,&max2,&min2,&mean,&stdev);
      if (max2-min2<=tolerance){ /* do nothing */ }
      else /* if (max2>min2+tolerance) */{ 
	max = p->maxra[p->indexing_nvar_checkout[nv]]; min = p->minra[p->indexing_nvar_checkout[nv]];
	if (max<=min){ if (stdev_use){ max=mean+STD_VIEW*stdev; min=mean-STD_VIEW*stdev;} else if (!stdev_use){ max=max2;min=min2;}}
	sprintf(filename,"./power_draw_interlace_%s_%s_%s.pnm",gs2,GLOBAL_TYPENAMES[p->indexing_ntype_checkout[nt]],GLOBAL_VARNAMES[p->indexing_nvar_checkout[nv]]);
	WritePNMfile_color(ra,trial_length,ntrials*nrecords,max,min,filename,7);
	sprintf(filename_base,"./power_plot_interlace_%s_%s_%s",gs2,GLOBAL_TYPENAMES[p->indexing_ntype_checkout[nt]],GLOBAL_VARNAMES[p->indexing_nvar_checkout[nv]]);
	sprintf(filename,"%s.fig",filename_base);
	if ((fp=fopen(filename,"w"))==NULL){ printf(" %% Warning, cannot create %s in ra2jpg\n",filename); fp=stdout;}
	fprintf(fp,"%s",FIG_PREAMBLE);
	for (nr=0;nr<nrecords;nr++){
	  for (nl=0;nl<ntrials;nl++){
	    periodify("int",&nl,&lcolormin,&lcolormax,&lcolor);
	    fprintf(fp,"2 1 0 %d %d 7 %d 0 -1 0.000 0 0 -1 0 0 %d\n\t",/*width*/1,/*color*/lcolor,/*depth*/(nr%999)+1,/*npoints*/trial_length);
	    for (nl2=0;nl2<trial_length;nl2++){
	      xpos = ((double)nl2+0.5)/(double)trial_length;
	      ypos = (ra[nl2 + nl*trial_length + nr*trial_length*ntrials] - min)/(max-min)/(double)nrecords + (double)nr/(double)nrecords;
	      fprintf(fp," %d %d",(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos));}
	    fprintf(fp,"\n");}}
	if (fp!=stdout){ fclose(fp);}
	if (jpeg_flag){ sprintf(command,"fig2dev -Ljpeg -q %d %s.fig %s.jpg;",/*quality*/5,filename_base,filename_base); system(command);}
	if (remove_flag){ sprintf(command,"rm %s.fig;",filename_base); system(command);}}
      tfree(ra);ra=NULL;}}}
}

void powerdump(struct power *p,char *fgvn,int dump_type)
{
  char filename[1024];
  char gs2[512];
  FILE *fp=NULL;
  int nt=0,nv=0,nr=0;
  if (fgvn==NULL){ sprintf(gs2,"%srecord",GLOBAL_STRING_2);}
  else{ sprintf(gs2,"%srecord_%s",GLOBAL_STRING_2,fgvn);}
  for (nt=0;nt<p->indexing_ntype_length;nt++){
    for (nv=0;nv<p->indexing_nvar_length;nv++){
      if (p->strarara[nt][nv]!=NULL){
	sprintf(filename,"powerstra_%s_%s_%s",gs2,GLOBAL_TYPENAMES[p->indexing_ntype_checkout[nt]],GLOBAL_VARNAMES[p->indexing_nvar_checkout[nv]]); 
	stradump(p->strarara[nt][nv],p->Nra->lengthra[p->indexing_ntype_checkout[nt]],dump_type,filename);}}}
  if (p->update_number>0){
    if (dump_type==0){
      sprintf(filename,"./power_%s.m",gs2);
      if ((fp=fopen(filename,"w"))!=NULL){
	for (nt=0;nt<p->indexing_ntype_length;nt++){ for (nv=0;nv<p->indexing_nvar_length;nv++){ for (nr=0;nr<p->length;nr++){ fprintf(fp,"power_%s_vpowrarara_%s_%s(%d) = %0.16f;\n",gs2,GLOBAL_TYPENAMES[p->indexing_ntype_checkout[nt]],GLOBAL_VARNAMES[p->indexing_nvar_checkout[nv]],nr+1,p->vpowrarara[nt][nv][nr]/(double)p->update_number);}}}
	for (nr=0;nr<p->length;nr++){ fprintf(fp,"power_%s_firstvarpower(%d) = %0.16f;\n",gs2,nr+1,p->firstvarpower[nr]);}
	if (p->correlation_bother){
	  for (nr=0;nr<p->length;nr++){ fprintf(fp,"power_%s_ac(%d) = %0.16f;\n",gs2,nr+1,p->autocorrelation[nr]);}
	  for (nr=0;nr<p->length;nr++){ fprintf(fp,"power_%s_xc(%d) = %0.16f;\n",gs2,nr+1,p->crosscorrelation[nr]);}}
	fprintf(fp,"figure;clf;hold on;\n");
	fprintf(fp,"subplot(1,1,1);semilogy(power_%s_firstvarpower);title('firstvar');\n",gs2);
	fclose(fp);fp=NULL;}}
    else if (dump_type==1){
      for (nt=0;nt<p->indexing_ntype_length;nt++){
	for (nv=0;nv<p->indexing_nvar_length;nv++){
	  if (p->vpowrarara[nt][nv]!=NULL){
	    sprintf(filename,"./powervpow_%s_%s_%s",gs2,GLOBAL_TYPENAMES[p->indexing_ntype_checkout[nt]],GLOBAL_VARNAMES[p->indexing_nvar_checkout[nv]]);
	    ra2jpg(p->vpowrarara[nt][nv],"double",1,p->length,0,filename,3);}}}
      sprintf(filename,"./power_%s_firstvar",gs2);
      ra2jpg(p->firstvarpower,"double",1,p->length,0,filename,3);
      if (p->correlation_bother){
	sprintf(filename,"./power_%s_autocorrelation",gs2);
	ra2jpg(p->autocorrelation,"double",1,p->length,0,filename,0);
	sprintf(filename,"./power_%s_crosscorrelation",gs2);
	ra2jpg(p->crosscorrelation,"double",1,p->length,0,filename,0);}}}
}

/* Here are the region functions anew */

void regionramake(struct ptree *p,double event_within,int event_threshold,int region_type)
{
  /* intended for orientation domains */
  int verbose=0;
  struct neuronarray *Nra=GLOBAL_Nra;
  int nr=0,ni=0,nt=0;
  if (verbose){ printf(" %% [entering regionramake] event_within %f event_threshold %d region type %d nregions %d\n",event_within,event_threshold,region_type,p->nregions);}
  if (verbose){ printf(" %% making regions\n");}
  for (nr=0;nr<p->nregions;nr++){ 
    p->regionra[nr] = (struct region *) tmalloc(sizeof(struct region)); 
    p->regionra[nr]->region_type = region_type;
    p->regionra[nr]->label = nr; /* critical to have label nr match location in regionra */
    p->regionra[nr]->last_event = GLOBAL_TI-rand01;
    p->regionra[nr]->event_within = event_within;
    p->regionra[nr]->event_threshold = event_threshold;
    p->regionra[nr]->neuronllist = llistmake();
    p->regionra[nr]->pn = NULL;}
  switch (region_type){
  case 0: if (verbose){ printf(" region_type %d, NULL\n",region_type);} break;
  case 1: case 2: /* individual neurons by index */
    if (verbose){ printf(" region_type %d, individual neurons by index\n",region_type);}
    for (nr=0;nr<p->nregions;nr++){
      nt=0; ni=nr; while (nt<Nra->ntypes && ni>=Nra->lengthra[nt]){ ni -= Nra->lengthra[nt]; nt += 1;}
      if (verbose){ printf(" region %d, nt %d, ni %d\n",nr,nt,ni);}
      if (nt<Nra->ntypes){ if (ni<Nra->lengthra[nt]){ litemadd(p->regionra[nr]->neuronllist,nget(Nra,nt,ni));}}}
    break;
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
  switch (r->region_type){
  case 2:
    if (verbose){ printf(" %% region_type implies single neuron's input, although multiple neurons will be handled correctly\n");}
    threshold=minimum(L->length,r->event_threshold);
    within=r->event_within;
    if (verbose){ printf(" %% event_threshold %d event_within %0.2f\n",threshold,within);}
    if (L->length>0){ llistsort(L->first,L->last,L->length,&spikeinput_time_compare);}
    exit_flag=0; depth=0;
    if (L->last!=NULL){
      l=L->last;
      n=(struct neuron *)L->last->item; spikelast = n->spikeinput_time;
      if (n->spikeinput_time>r->last_event && n->spikeinput_time>=t && n->spikeinput_time<=t+DT){
	if (verbose){ printf(" %% at final neuron (%d,%d) with spikeinput_time %f\n",n->type,n->index,n->spikeinput_time);}
	l=l->parent; depth=1;
	while (l!=NULL && !exit_flag){
	  n=(struct neuron *)l->item;
	  if (verbose){ printf(" %% at parent n(%d,%d) with spikeinput_time %f -- depth=%d\n",n->type,n->index,n->spikeinput_time,depth);}
	  if (n->spikeinput_time>r->last_event && (spikelast-n->spikeinput_time)<within){ l=l->parent; depth+=1;}
	  else{ exit_flag=1;}}}
      if (depth>=threshold){
	if (verbose){ printf(" %% depth=%d, event happened\n",depth);}
	r->last_event = spikelast;
	return 1;}
      else{ 
	if (verbose){ printf(" %% depth=%d, event didn't happened\n",depth);}
	return 0;}}
    return 0;
    break;
  case 1: 
    if (verbose){ printf(" %% region_type implies single neuron, although multiple neurons will be handled correctly\n");}
    threshold=minimum(L->length,r->event_threshold);
    within=r->event_within;
    if (verbose){ printf(" %% event_threshold %d event_within %0.2f\n",threshold,within);}
    if (L->length>0){ llistsort(L->first,L->last,L->length,&spikelast_compare);}
    exit_flag=0; depth=0;
    if (L->last!=NULL){
      l=L->last;
      n=(struct neuron *)L->last->item; spikelast = n->spikelast;
      if (n->spikelast>r->last_event && n->spikelast>=t && n->spikelast<=t+DT){
	if (verbose){ printf(" %% at final neuron (%d,%d) with spikelast %f\n",n->type,n->index,n->spikelast);}
	l=l->parent; depth=1;
	while (l!=NULL && !exit_flag){
	  n=(struct neuron *)l->item;
	  if (verbose){ printf(" %% at parent neuron (%d,%d) with spikelast %f -- depth=%d\n",n->type,n->index,n->spikelast,depth);}
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
    break;
  default: return 0; break;}
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

/* Here are the clusterdata functions */

struct clusterdatara * clusterdataramake(struct clusterarray *gli,int power_bother,int power_length,int indexing_nvar_length,int *indexing_nvar_checkout,int indexing_nvar_refile_length,int *indexing_nvar_refile,int trajectory_window_length,int power_window_length,int power_window_update_every,double *maxra,double *minra,int cycle_bother,int ptree_bother,int nlegs,int legtime)
{
  /* assumes a premade cluster *c */  
  int verbose=0;
  int nc=0,nr=0,nv=0;
  struct litem *l=NULL;
  struct clusterdatara *cd=NULL;
  struct ptree *p=NULL;
  if (verbose){ 
    printf(" %% [entering clusterdataramake]\n");
    raprintf(indexing_nvar_checkout,"int",1,indexing_nvar_length,"nvar_checkout");
    raprintf(indexing_nvar_refile,"int",1,indexing_nvar_refile_length,"nvar_refile");
    raprintf(maxra,"double",1,indexing_nvar_refile_length,"maxra");
    raprintf(minra,"double",1,indexing_nvar_refile_length,"minra");}
  cd = (struct clusterdatara *) tcalloc(1,sizeof(struct clusterdatara));
  cd->gli = gli;
  cd->power_bother = power_bother;
  cd->vpow_bother = 0;
  cd->lfp_gammalo = 24;
  cd->lfp_gammahi = 48;
  if (cd->power_bother){
    cd->power_update_every=1;
    cd->trajectory_window = trajectory_window_length/cd->power_update_every;
    cd->power_length=power_length;
    cd->indexing_nvar_length = indexing_nvar_length;
    cd->indexing_nvar_checkout = (int *) tcalloc(cd->indexing_nvar_length,sizeof(int));
    for (nv=0;nv<cd->indexing_nvar_length;nv++){ cd->indexing_nvar_checkout[nv] = indexing_nvar_checkout[nv];}
    cd->indexing_nvar_refile_length = indexing_nvar_refile_length;
    cd->indexing_nvar_refile = (int *) tcalloc(cd->indexing_nvar_refile_length,sizeof(int));
    for (nv=0;nv<cd->indexing_nvar_refile_length;nv++){ cd->indexing_nvar_refile[nv] = indexing_nvar_refile[nv];}
    cd->tst = strobemake(cd->power_length,cd->power_update_every,0,0,0,0);
    cd->strarara = (struct strobe ****) tcalloc(cd->gli->nclusters,sizeof(struct strobe ***));
    if (cd->vpow_bother){ cd->vpowrarara = (double ***) tcalloc(cd->gli->nclusters,sizeof(double **));}
    for (nc=0;nc<cd->gli->nclusters;nc++){
      cd->strarara[nc] = (struct strobe ***) tcalloc(cd->indexing_nvar_length,sizeof(struct strobe **));
      if (cd->vpow_bother){ cd->vpowrarara[nc] = (double **) tcalloc(cd->indexing_nvar_length,sizeof(double *));}
      for (nv=0;nv<cd->indexing_nvar_length;nv++){
	cd->strarara[nc][nv] = (struct strobe **) tcalloc(cd->gli->cra[nc]->LN->length,sizeof(struct strobe *));
	for (nr=0;nr<cd->gli->cra[nc]->LN->length;nr++){
	  cd->strarara[nc][nv][nr] = strobemake(cd->power_length,cd->power_update_every,cycle_bother,0,0,0);}
	if (cd->vpow_bother){ cd->vpowrarara[nc][nv] = (double *) tcalloc(cd->power_length,sizeof(double));}}}
    cd->stra_lfp = (struct strobe **) tcalloc(cd->gli->nclusters,sizeof(struct strobe *));
    for (nc=0;nc<cd->gli->nclusters;nc++){ cd->stra_lfp[nc] = strobemake(cd->power_length,cd->power_update_every,cycle_bother,1,power_window_length,power_window_update_every);}
    cd->maxra= (double *) tcalloc(cd->indexing_nvar_refile_length,sizeof(double));
    cd->minra= (double *) tcalloc(cd->indexing_nvar_refile_length,sizeof(double));
    for (nv=0;nv<indexing_nvar_refile_length;nv++){ cd->maxra[nv] = maxra[nv]; cd->minra[nv] = minra[nv];}
    cd->power_update_number=0;}
  cd->ptree_bother = ptree_bother;
  if (cd->ptree_bother){
    cd->pra = (struct ptree **) tcalloc(cd->gli->nclusters,sizeof(struct ptree *));
    for (nc=0;nc<cd->gli->nclusters;nc++){
      p = (struct ptree *) tcalloc(1,sizeof(struct ptree));
      p->nregions = cd->gli->cra[nc]->LN->length;
      p->regionra = (struct region **) tcalloc(p->nregions,sizeof(struct region *));
      l=cd->gli->cra[nc]->LN->first; nr=0;
      while (l!=NULL){
	p->regionra[nr] = (struct region *) tmalloc(sizeof(struct region)); 
	p->regionra[nr]->label = nr; /* critical to have label nr match location in regionra */
	p->regionra[nr]->last_event = GLOBAL_TI-rand01;
	p->regionra[nr]->event_within = 1;
	p->regionra[nr]->event_threshold = 1;
	p->regionra[nr]->neuronllist = llistmake(); litemadd(p->regionra[nr]->neuronllist,(struct neuron *)l->item);
	p->regionra[nr]->pn = NULL;
	l=l->child;
	nr+=1;}
      p->nlegs = nlegs;
      p->legtime = legtime;
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
      cd->pra[nc] = p;}}
  return cd;
}

void clusterdatatfree(struct clusterdatara *cd)
{
  int verbose=0;
  int nc=0,nr=0,nv=0;
  if (verbose){ printf(" %% [entering clusterdatatfree]\n");}
  if (cd->power_bother){
    if (verbose){ printf(" %% %% tfreeing \n");}
    if (verbose){ printf(" %% %% tfreeing cd->tst\n");}
    strobetfree(cd->tst);cd->tst=NULL;
    for (nc=0;nc<cd->gli->nclusters;nc++){
      if (verbose){ printf(" %% %% cluster %d\n",nc);}
      if (cd->stra_lfp[nc]!=NULL){ 
	if (verbose){ printf(" %% %% tfreeing cd->str_lfp[%d]\n",nc);} 
	strobetfree(cd->stra_lfp[nc]);cd->stra_lfp[nc]=NULL;}
      for (nv=0;nv<cd->indexing_nvar_length;nv++){
	if (verbose){ printf(" %% %% %% variable %d \n",nv);}
	if (cd->strarara[nc][nv]!=NULL){
	  for (nr=0;nr<cd->gli->cra[nc]->LN->length;nr++){
	    if (verbose){ printf(" %% %% %% tfreeing cd->strarara[%d][%d][%d]\n",nc,nv,nr);}
	    strobetfree(cd->strarara[nc][nv][nr]); cd->strarara[nc][nv][nr]=NULL;}
	  if (verbose){ printf(" %% %% %% tfreeing cd->strarara[%d][%d]\n",nc,nv);}
	  tfree(cd->strarara[nc][nv]);cd->strarara[nc][nv]=NULL;}
	if (cd->vpow_bother && cd->vpowrarara[nc][nv]!=NULL){
	  if (verbose){ printf(" %% %% %% tfreeing cd->vpowrarara[%d][%d]\n",nc,nv);}
	  tfree(cd->vpowrarara[nc][nv]);cd->vpowrarara[nc][nv]=NULL;}}
      if (verbose){ printf(" %% %% tfreeing cd->strarara[%d] \n",nc);}
      tfree(cd->strarara[nc]);cd->strarara[nc]=NULL;
      if (cd->vpow_bother){ 
	if (verbose){ printf(" %% %% tfreeing cd->vpowrarara[%d] \n",nc);}
	tfree(cd->vpowrarara[nc]);cd->vpowrarara[nc]=NULL;}}
    if (verbose){ printf(" %% tfreeing cd->strarara\n");}
    tfree(cd->strarara);cd->strarara=NULL;
    if (cd->vpow_bother){
      if (verbose){ printf(" %% tfreeing cd->vpowrarara\n");}
      tfree(cd->vpowrarara);cd->vpowrarara=NULL;}
    if (verbose){ printf(" %% tfreeing cd->stra_lfp\n");}
    tfree(cd->stra_lfp);cd->stra_lfp=NULL;
    if (verbose){ printf(" %% tfreeing checkout and refile\n");}
    tfree(cd->indexing_nvar_checkout); cd->indexing_nvar_checkout=NULL;
    tfree(cd->indexing_nvar_refile); cd->indexing_nvar_refile=NULL;
    if (verbose){ printf(" %% tfreeing maxra and minra\n");}
    tfree(cd->maxra);cd->maxra=NULL;
    tfree(cd->minra);cd->minra=NULL;}
  if (cd->ptree_bother){ for (nc=0;nc<cd->gli->nclusters;nc++){ 
    if (verbose){ printf(" %% tfreeing cd->pra[%d]\n",nc);}
    ptreetfree(cd->pra[nc]); cd->pra[nc]=NULL;} 
  if (verbose){ printf(" %% tfreeing cd->pra\n");}
  tfree(cd->pra); cd->pra=NULL;}
  if (verbose){ printf(" %% tfreeing cd\n");}
  tfree(cd); cd=NULL;
  if (verbose){ printf(" %% [finishing clusterdatatfree]\n");}
}

void clusterdataraupdate(struct clusterdatara *cd,double t,double DT)
{
  int verbose=0;
  int nc=0,nr=0,nv=0;
  struct neuron *n=NULL;
  int tab_old=0,tab_new=0;
  double *temp=NULL,*lfptemp=NULL;
  struct litem *l=NULL;
  int varname_registry_lfp=0;
  if (verbose){ printf(" %% [entering clusterdataraupdate] t %f dt %f\n",t,DT);}
  if (cd->power_bother){
    tab_old = cd->tst->tab; strobeupdate_old(cd->tst,t,t); tab_new = cd->tst->tab;
    for (nc=0;nc<cd->gli->nclusters;nc++){
      if (cd->stra_lfp[nc]!=NULL){
	lfptemp = (double *) tcalloc(cd->gli->cra[nc]->LN->length+1,sizeof(double));
	l=cd->gli->cra[nc]->LN->first;nr=0;
	while (l!=NULL){ n=(struct neuron *)l->item;lfptemp[nr]=*(n->vpra[cd->indexing_nvar_checkout[varname_registry_lfp]]);l=l->child;nr+=1;}
	stats("double",lfptemp,cd->gli->cra[nc]->LN->length,NULL,NULL,&(lfptemp[cd->gli->cra[nc]->LN->length]),NULL);
	strobeupdate(cd->stra_lfp[nc],t,DT,lfptemp[cd->gli->cra[nc]->LN->length]);
	tfree(lfptemp);lfptemp=NULL;}
      for (nv=0;nv<cd->indexing_nvar_length;nv++){
	if (verbose>1){ printf(" %% nv %d (%s)\n",nv,GLOBAL_VARNAMES[cd->indexing_nvar_checkout[nv]]);}
	if (cd->strarara[nc][nv]!=NULL){
	  l=cd->gli->cra[nc]->LN->first; nr=0;
	  while (l!=NULL){ 
	    n=(struct neuron *)l->item; 
	    if (verbose>1){ printf(" %% updating cd->strarara[nc=%d][nv=%d %s][nr=%d] with %f...\n",nc,nv,GLOBAL_VARNAMES[cd->indexing_nvar_checkout[nv]],nr,*(n->vpra[cd->indexing_nvar_checkout[nv]]));}
	    strobeupdate(cd->strarara[nc][nv][nr],t,DT,*(n->vpra[cd->indexing_nvar_checkout[nv]]));
	    if (verbose>1){ printf(" %% finished\n");}
	    l=l->child;nr+=1;}}}}
    if (cd->vpow_bother){
      if (tab_old==cd->power_length-1 && tab_new==0){ 
	if (verbose){ printf(" %% time to update\n");}
	for (nc=0;nc<cd->gli->nclusters;nc++){
	  for (nv=0;nv<cd->indexing_nvar_length;nv++){
	    if (cd->strarara[nc][nv]!=NULL && cd->vpowrarara[nc][nv]!=NULL){
	      for (nr=0;nr<cd->gli->cra[nc]->LN->length;nr++){
		temp = ra2power(NULL,cd->strarara[nc][nv][nr]->data,cd->power_length,NULL,0,0);
		raplusequals(cd->vpowrarara[nc][nv],cd->power_length,temp);
		tfree(temp);temp=NULL;}}}}
	cd->power_update_number += 1;}}}
  if (cd->ptree_bother){ for (nc=0;nc<cd->gli->nclusters;nc++){ ptreeupdate(cd->pra[nc],t,DT,1);}}
}

void clusterdatarareset(struct clusterdatara *cd)
{
  int nc=0,nr=0,nv=0;
  if (cd->power_bother){
    strobereset(cd->tst);
    for (nc=0;nc<cd->gli->nclusters;nc++){
      if (cd->stra_lfp[nc]!=NULL){ strobereset(cd->stra_lfp[nc]);}
      for (nv=0;nv<cd->indexing_nvar_length;nv++){
	if (cd->strarara[nc][nv]!=NULL){
	  for (nr=0;nr<cd->gli->cra[nc]->LN->length;nr++){
	    strobereset(cd->strarara[nc][nv][nr]);}}
	if (cd->vpow_bother && cd->vpowrarara[nc][nv]!=NULL){
	  rareset(cd->vpowrarara[nc][nv],"double",cd->power_length,NULL);}}}
    cd->power_update_number=0;}
  if (cd->ptree_bother){ for (nc=0;nc<cd->gli->nclusters;nc++){ 
    pnodeclear_starter(NULL,cd->pra[nc]->pretree); 
    pnodeclear_starter(NULL,cd->pra[nc]->postree); 
    cd->pra[nc]->total_time=0;}}
}

void clusterdataradump(struct clusterdatara *cd,char *fgvn,int dump_what)
{
  /* assumes *fgvn does not have a "./" at the beginning */
  int verbose=0;
  char filename[1024];
  char gs2[512];
/*   struct strobe **stra=NULL; */
  int nc=0,nv=0,nr=0,length=0,length2=0,length3=0,nt=0,nr2=0;
  double *temp=NULL;
/*   double *temp2=NULL; */
  int wrap_at=1024;
  int tab=0;
  if (fgvn==NULL){ sprintf(gs2,"%srecord",GLOBAL_STRING_2);}
  else{ sprintf(gs2,"%srecord_%s",GLOBAL_STRING_2,fgvn);}
  switch (dump_what){
  case 0: 
    if (verbose){ printf("%% clusterdararadump everything \n");}
    if (cd->power_bother){
      for (nv=0;nv<cd->indexing_nvar_length;nv++){
	length=0; for (nc=0;nc<cd->gli->nclusters;nc++){ if (cd->strarara[nc][nv]!=NULL){ length+=cd->gli->cra[nc]->LN->length;}}
	length2=cd->power_length;length3=length2/wrap_at+(length2%wrap_at>0);if (length3>1){ length2 = wrap_at*length3;}
	temp = (double *) tcalloc(length*length2,sizeof(double));
	/*       stra = (struct strobe **) tcalloc(length,sizeof(struct strobe *)); */
	nr=0;
	for (nc=0;nc<cd->gli->nclusters;nc++){
	  if (cd->strarara[nc][nv]!=NULL){ 
	    for (nr2=0;nr2<cd->gli->cra[nc]->LN->length;nr2++){
	      for (nt=0;nt<cd->power_length;nt++){
		tab = periodize(cd->strarara[nc][nv][nr2]->tab+nt,0,cd->power_length);
		temp[nt/wrap_at+nr*length3 + (nt%wrap_at)*length*length3] = cd->strarara[nc][nv][nr2]->data[tab];}
	      /* 	    stra[nr] = cd->strarara[nc][nv][nr2]; */
	      nr+=1;}}}
	sprintf(filename,"./cluster_all_stra_%s_%s.pnm",GLOBAL_VARNAMES[cd->indexing_nvar_checkout[nv]],gs2);
	WritePNMfile_color(temp,length*length3,minimum(wrap_at,length2),0,0,filename,7);
	tfree(temp);temp=NULL;
	/*       temp = stra2ra(stra,length,cd->trajectory_window,0,NULL); */
	/*       temp2 = ra2pca(temp,cd->power_length/cd->trajectory_window,length,NULL); */
	/*       sprintf(filename,"./cluster_all_pca_%s_%s",GLOBAL_VARNAMES[cd->indexing_nvar_checkout[nv]],gs2); */
	/*       pca2jpg(temp2,cd->power_length/cd->trajectory_window,filename); */
	/*       tfree(stra);stra=NULL; */
	/*       tfree(temp);temp=NULL; */
	/*       tfree(temp2);temp2=NULL; */
      }
      for (nc=0;nc<cd->gli->nclusters;nc++){
	if (cd->stra_lfp[nc]!=NULL){
	  sprintf(filename,"cluster%dof%d_stra_lfp_%s",nc,cd->gli->nclusters,gs2);
	  stradump(&(cd->stra_lfp[nc]),1,1,filename);}
	for (nv=0;nv<cd->indexing_nvar_length;nv++){
	  if (cd->strarara[nc][nv]!=NULL){
	    sprintf(filename,"cluster%dof%d_stra_%s_%s",nc,cd->gli->nclusters,GLOBAL_VARNAMES[cd->indexing_nvar_checkout[nv]],gs2); 
	    stradump(cd->strarara[nc][nv],cd->gli->cra[nc]->LN->length,1,filename);}
	  if (cd->vpow_bother && cd->vpowrarara[nc][nv]!=NULL && cd->power_update_number>0){
	    sprintf(filename,"./cluster%dof%d_vpow_%s_%s",nc,cd->gli->nclusters,GLOBAL_VARNAMES[cd->indexing_nvar_checkout[nv]],gs2);
	    ra2jpg(cd->vpowrarara[nc][nv],"double",1,cd->power_length,0,filename,3);}}}}
    if (cd->ptree_bother){
      for (nc=0;nc<cd->gli->nclusters;nc++){
	sprintf(filename,"cluster%dof%d_%s",nc,cd->gli->nclusters,fgvn==NULL? "_":fgvn); 
	ptreedump_starter(cd->pra[nc],filename,2,0,0,0,+1,-1,0);}}
    break;
  case 1:
    if (verbose){ printf("%% clusterdararadump lfp only \n");}
    if (cd->power_bother){
      for (nc=0;nc<cd->gli->nclusters;nc++){
	if (cd->stra_lfp[nc]!=NULL){
	  sprintf(filename,"cluster%dof%d_stra_lfp_%s",nc,cd->gli->nclusters,gs2);
	  stradump(&(cd->stra_lfp[nc]),1,1,filename);}}}
    break;
  case 2: 
    if (verbose){ printf("%% clusterdararadump ptree only \n");}
    if (cd->ptree_bother){
      for (nc=0;nc<cd->gli->nclusters;nc++){
	sprintf(filename,"cluster%dof%d_%s",nc,cd->gli->nclusters,fgvn==NULL? "_":fgvn); 
	ptreedump_starter(cd->pra[nc],filename,2,0,0,0,+1,-1,0);}}
    break;
  default: break;}
}

/* Here are caicor functions */

struct caicor * caicormake(struct neuronarray *Nra,int indexing_nvar_length,int *indexing_nvar_checkout,int indexing_nvar_refile_length,int *indexing_nvar_refile,double *maxra,double *minra,int nbins){
  struct caicor *c=NULL;
  int nt=0,ni=0,nv=0;
  c = (struct caicor *) tcalloc(1,sizeof(struct caicor));
  c->LN = llistmake(); for (nt=0;nt<Nra->ntypes;nt++){ for (ni=0;ni<Nra->lengthra[nt];ni++){ litemadd(c->LN,nget(Nra,nt,ni));}}
  c->indexing_nvar_length = indexing_nvar_length;
  c->indexing_nvar_checkout = (int *) tcalloc(c->indexing_nvar_length,sizeof(int));
  for (nv=0;nv<c->indexing_nvar_length;nv++){ c->indexing_nvar_checkout[nv] = indexing_nvar_checkout[nv];}
  c->indexing_nvar_refile_length = indexing_nvar_refile_length;
  c->indexing_nvar_refile = (int *) tcalloc(c->indexing_nvar_refile_length,sizeof(int));
  for (nv=0;nv<c->indexing_nvar_refile_length;nv++){ c->indexing_nvar_refile[nv] = indexing_nvar_refile[nv];}
  c->maxra= (double *) tcalloc(c->indexing_nvar_refile_length,sizeof(double));
  c->minra= (double *) tcalloc(c->indexing_nvar_refile_length,sizeof(double));
  for (nv=0;nv<indexing_nvar_refile_length;nv++){ c->maxra[nv] = maxra[nv]; c->minra[nv] = minra[nv];}
  c->hra = (struct hist **) tcalloc(c->indexing_nvar_length,sizeof(struct hist *));
  c->hstra = (struct hist **) tcalloc(c->indexing_nvar_length,sizeof(struct hist *));
  c->nbins = nbins;
  for (nv=0;nv<c->indexing_nvar_length;nv++){ 
    c->hra[nv]=histmake(c->nbins,c->maxra[c->indexing_nvar_checkout[nv]],c->minra[c->indexing_nvar_checkout[nv]]);
    c->hstra[nv]=histmake(c->nbins,c->maxra[c->indexing_nvar_checkout[nv]],c->minra[c->indexing_nvar_checkout[nv]]);}
  c->mrara = (double *) tcalloc(c->nbins*2*c->indexing_nvar_length,sizeof(double));
  c->mrarara = (double *) tcalloc(c->nbins*c->nbins*2*c->indexing_nvar_length,sizeof(double));
  return c;
}

void caicortfree(struct caicor *c)
{
  int nv=0;  
  for (nv=0;nv<c->indexing_nvar_length;nv++){ 
    histtfree(c->hra[nv]); c->hra[nv]=NULL;
    histtfree(c->hstra[nv]); c->hstra[nv]=NULL;} 
  tfree(c->hra); c->hra=NULL;
  tfree(c->hstra); c->hstra=NULL;
  tfree(c->indexing_nvar_checkout);c->indexing_nvar_checkout=NULL;
  tfree(c->indexing_nvar_refile);c->indexing_nvar_refile=NULL;
  tfree(c->maxra);c->maxra=NULL;
  tfree(c->minra);c->minra=NULL;
  tfree(c->mrara);c->mrara=NULL;
  tfree(c->mrarara);c->mrarara=NULL;
}

void caicorupdate(struct caicor *c,double t,double DT,int nspikes)
{
  /* nspikes over interval [t,t+DT] */
  int verbose=0;
  struct litem *l=NULL,*l2=NULL;
  struct neuron *n=NULL,*n2=NULL;
  int nv=0,spike_flag=0,tab=0,tab2;
  double temp=0,temp2=0;
  if (verbose){ printf(" %% [entering caicorupdate] t %f dt %f nspikes %d\n",t,DT,nspikes);}
  l=c->LN->first;
  while (l!=NULL){ 
    n=(struct neuron *)l->item;
    spike_flag = (n->spikelast==n->spiketime && n->spiketime>=t && n->spiketime <= t+DT);
    for (nv=0;nv<c->indexing_nvar_length;nv++){ 
      temp = *(n->vpra[c->indexing_nvar_checkout[nv]]);
      tab = maximum(0,minimum(c->nbins-1,(int)floor(c->nbins*(temp-c->minra[c->indexing_nvar_checkout[nv]])/(c->maxra[c->indexing_nvar_checkout[nv]]-c->minra[c->indexing_nvar_checkout[nv]]))));
      tab2 = tab + 0*c->nbins + nv*2*c->nbins;
      if (verbose>1){ printf(" %% using %f to add to c->mrara[%d + ?*%d + %d*2*%d]\n",temp,tab,c->nbins,nv,c->nbins);}
      c->mrara[tab2] += nspikes;
      tab2 = tab + 1*c->nbins + nv*2*c->nbins;
      c->mrara[tab2] += DT;
      /* fix later, O(pow(LN->length,2)) operations */
      l2=l;
      while (l2!=NULL){
	n2=(struct neuron *)l2->item;
	temp2 = *(n2->vpra[c->indexing_nvar_checkout[nv]]);
	tab = maximum(0,minimum(c->nbins-1,(int)floor(c->nbins*(temp-c->minra[c->indexing_nvar_checkout[nv]])/(c->maxra[c->indexing_nvar_checkout[nv]]-c->minra[c->indexing_nvar_checkout[nv]])))) + c->nbins*maximum(0,minimum(c->nbins-1,(int)floor(c->nbins*(temp2-c->minra[c->indexing_nvar_checkout[nv]])/(c->maxra[c->indexing_nvar_checkout[nv]]-c->minra[c->indexing_nvar_checkout[nv]]))));
	tab2 = tab + 0*c->nbins*c->nbins + nv*2*c->nbins*c->nbins;
	c->mrarara[tab2] += (l==l2?1:2)*nspikes;
	tab2 = tab + 1*c->nbins*c->nbins + nv*2*c->nbins*c->nbins;
	c->mrarara[tab2] += DT;
	l2=l2->child;}
      if (verbose>1){ printf(" %% adding %f to c->h[%s]\n",temp,GLOBAL_VARNAMES[c->indexing_nvar_checkout[nv]]);}
      if (nspikes>0){ histadd(c->hstra[nv],*(n->vpra[c->indexing_nvar_checkout[nv]]),nspikes-spike_flag);}
      histadd(c->hra[nv],*(n->vpra[c->indexing_nvar_checkout[nv]]),1);}
    l=l->child;}  
  if (verbose){ raprintf(c->hra[0]->data,"double",1,c->nbins,"h");}
  if (verbose){ raprintf(c->mrara,"double",1,c->nbins,"m");}
}

void caicordump(struct caicor *c,char *fgvn)
{
  char filename[1024];
  char gs2[512];
  int nv=0,nb=0,nb2=0;
  double *temp=NULL,*temp2=NULL;
  if (fgvn==NULL){ sprintf(gs2,"%srecord",GLOBAL_STRING_2);}
  else{ sprintf(gs2,"%srecord_%s",GLOBAL_STRING_2,fgvn);}
  for (nv=0;nv<c->indexing_nvar_length;nv++){
    sprintf(filename,"./caicor_hstra_%s_%s",GLOBAL_VARNAMES[c->indexing_nvar_checkout[nv]],gs2);
    ra2jpg(c->hstra[nv]->data,"double",1,c->hstra[nv]->nbins,0,filename,0);
    sprintf(filename,"./caicor_hra_%s_%s",GLOBAL_VARNAMES[c->indexing_nvar_checkout[nv]],gs2);
    ra2jpg(c->hra[nv]->data,"double",1,c->hra[nv]->nbins,0,filename,0);
    sprintf(filename,"./caicor_m0_%s_%s",GLOBAL_VARNAMES[c->indexing_nvar_checkout[nv]],gs2);
    ra2jpg(&(c->mrara[0+0*c->nbins+nv*2*c->nbins]),"double",1,c->nbins,0,filename,0);
    sprintf(filename,"./caicor_m1_%s_%s",GLOBAL_VARNAMES[c->indexing_nvar_checkout[nv]],gs2);
    ra2jpg(&(c->mrara[0+1*c->nbins+nv*2*c->nbins]),"double",1,c->nbins,0,filename,0);
    temp = (double *) tcalloc(c->nbins,sizeof(double));
    for (nb=0;nb<c->nbins;nb++){ temp[nb] = c->mrara[nb+0*c->nbins+nv*2*c->nbins]/maximum(1,c->mrara[nb+1*c->nbins+nv*2*c->nbins]);}
    sprintf(filename,"./caicor_C1_%s_%s",GLOBAL_VARNAMES[c->indexing_nvar_checkout[nv]],gs2);
    ra2jpg(temp,"double",1,c->nbins,0,filename,0);
    temp2=(double *)tcalloc(c->nbins*c->nbins,sizeof(double));
    for (nb=0;nb<c->nbins;nb++){ for (nb2=0;nb2<c->nbins;nb2++){ temp2[nb+nb2*c->nbins]=temp[nb]*temp[nb2];}}
    tfree(temp);temp=NULL;
    sprintf(filename,"./caicor_C1C1_%s_%s.pnm",GLOBAL_VARNAMES[c->indexing_nvar_checkout[nv]],gs2);
    WritePNMfile_color(temp2,c->nbins,c->nbins,0,0,filename,7);
    temp = (double *) tcalloc(c->nbins*c->nbins,sizeof(double));
    for (nb=0;nb<c->nbins*c->nbins;nb++){ temp[nb] = c->mrarara[nb+0*c->nbins*c->nbins+nv*2*c->nbins*c->nbins]/maximum(1,c->mrarara[nb+1*c->nbins*c->nbins+nv*2*c->nbins*c->nbins]);}
    sprintf(filename,"./caicor_C2_%s_%s.pnm",GLOBAL_VARNAMES[c->indexing_nvar_checkout[nv]],gs2);
    WritePNMfile_color(temp,c->nbins,c->nbins,0,0,filename,7);
    tfree(temp);temp=NULL;
    tfree(temp2);temp2=NULL;}
}

/* Here are lyapunov functions */

struct isi * isimake(struct neuronarray *Nra,int nbins,double maxisi,double minisi)
{
  struct isi *i=NULL;
  i = (struct isi *) tcalloc(1,sizeof(struct isi));
  i->Nra = Nra;
  i->L = llistmake();
  i->L2 = llistmake();
  i->length_max = 0;
  i->nbins = nbins; if (i->nbins>0){ i->length_max=16;}
  i->maxisi = maxisi;
  i->minisi = minisi;
  i->ra = (double *) tcalloc(nbins*nbins,sizeof(double));
  return i;
}

void isitfree(struct isi *i){ llisttfree2(i->L); i->L=NULL; llisttfree2(i->L2); i->L2=NULL; tfree(i->ra); i->ra=NULL; tfree(i); i=NULL;}

void isiupdate(struct isi *i,double t)
{
  int verbose=0;
  double *temp=NULL;
  int nb1=0,nb2=0;
  if (i->L->length>0 && (*(double *)i->L->last->item)==t){ /* do nothing */}
  else /* if new element */{ 
    if (verbose){ printf(" %% updating isi with %f\n",t);} 
    temp=(double *)tcalloc(1,sizeof(double));*temp=t;litemadd(i->L,temp); 
    if (i->length_max>0 && i->L->length>i->length_max){ llistkillfirst2(i->L);}
    if (i->L->length>1){
      temp=(double *)tcalloc(1,sizeof(double));
      *temp = *(double *)i->L->last->item - *(double *)i->L->last->parent->item;
      litemadd(i->L2,temp);
      if (i->length_max>0 && i->L2->length>i->length_max){ llistkillfirst2(i->L2);}
      if (i->L2->length>1){
	nb1 = crop((int)floor(i->nbins*(*(double *)i->L2->last->item - i->minisi)/(i->maxisi-i->minisi)),0,i->nbins-1);
	nb2 = crop((int)floor(i->nbins*(*(double *)i->L2->last->parent->item - i->minisi)/(i->maxisi-i->minisi)),0,i->nbins-1);
	i->ra[nb1 + nb2*i->nbins] += 1;}}}
}

void isidump(struct isi *i,char *fgvn)
{
  char filename[2048];
  char gs2[1024];
  FILE *fp=NULL;
  int nr1=0,nr2=0;
  if (fgvn==NULL){ sprintf(gs2,"%srecord",GLOBAL_STRING_2);}
  else{ sprintf(gs2,"%srecord_%s",GLOBAL_STRING_2,fgvn);}
  sprintf(filename,"isi_%s.m",gs2);
  if ((fp=fopen(filename,"w"))!=NULL){
    fprintf(fp," %% isi_%s dumped at time %d;\n",gs2,(int)GLOBAL_time);
    for (nr1=0;nr1<i->nbins;nr1++){ for (nr2=0;nr2<i->nbins;nr2++){
      fprintf(fp," isi_%s(%d,%d)=%f;\n",gs2,nr1+1,nr2+1,i->ra[nr1+nr2*i->nbins]);}}
    fclose(fp);fp=NULL;}
}

/* Here are the hhlib functions */

struct hhlib_hist *hhlib_histmake(int indexing_nvar_length,int triggering_index,int *libuseflagra,int *logflagra,double *maxra,double *minra,int nbins,int ncumulants)
{
  /* note that the trigger variable is not binned */
  int verbose=0;
  int nv=0;
  struct hhlib_hist *hlh=NULL;
  hlh = (struct hhlib_hist *) tcalloc(1,sizeof(struct hhlib_hist));
  hlh->nvars=indexing_nvar_length;
  hlh->indexing_trigger = triggering_index;
  hlh->libuseflagra = (int *) tcalloc(hlh->nvars,sizeof(int)); for (nv=0;nv<hlh->nvars;nv++){ hlh->libuseflagra[nv]=libuseflagra[nv];}
  hlh->logflagra = (int *) tcalloc(hlh->nvars,sizeof(int)); for (nv=0;nv<hlh->nvars;nv++){ hlh->logflagra[nv]=logflagra[nv];}
  hlh->maxra = (double *) tcalloc(hlh->nvars,sizeof(double)); 
  for (nv=0;nv<hlh->nvars;nv++){ hlh->maxra[nv]= hlh->logflagra[nv] ? log(maxra[nv]) : maxra[nv];}
  hlh->minra = (double *) tcalloc(hlh->nvars,sizeof(double)); 
  for (nv=0;nv<hlh->nvars;nv++){ hlh->minra[nv]= hlh->logflagra[nv] ? log(minra[nv]) : minra[nv];}
  if (verbose){
    raprintf(hlh->logflagra,"int",1,hlh->nvars,"logflagra: ");
    raprintf(hlh->libuseflagra,"int",1,hlh->nvars,"libuseflagra: ");
    raprintf(hlh->maxra,"double",1,hlh->nvars,"maxra: ");
    raprintf(hlh->minra,"double",1,hlh->nvars,"minra: ");}
  hlh->nbins=nbins;
  hlh->ncumulants=ncumulants;
  hlh->data = (double *) tcalloc(pow(hlh->nbins,hlh->nvars-1)*hlh->nvars*hlh->ncumulants,sizeof(double));
  hlh->data_project = (double *) tcalloc(hlh->nbins*(hlh->nvars-1)*hlh->nvars*hlh->ncumulants,sizeof(double));
  return hlh;
}

void hhlib_histupdate(struct hhlib_hist *hlh,double *valra_start,double *valra_end)
{
  /* valra_start contains all variables, even the trigger variable */
  int verbose=0;
  int nv=0,nv_temp=0,*nbra=NULL,tab=0,nc=0;
  double v=0;
  if (verbose){ printf(" %% [entering hhlib_histupdate]\n");}
  if (verbose){ raprintf(valra_start,"double",1,hlh->nvars," %% %% valra_start: ");}
  if (verbose){ raprintf(valra_end,"double",1,hlh->nvars," %% %% valra_end: ");}
  nbra = (int *) tcalloc(hlh->nvars-1,sizeof(int));
  nv_temp=0;
  for (nv=0;nv<hlh->nvars;nv++){
    if (nv!=hlh->indexing_trigger){
      v = hlh->logflagra[nv] ? log(valra_start[nv]) : valra_start[nv];
      nbra[nv_temp]=maximum(0,minimum(hlh->nbins-1,(int)floor(hlh->nbins*(v-hlh->minra[nv])/(hlh->maxra[nv]-hlh->minra[nv]))));
      nv_temp+=1;}}
  if (verbose){ raprintf(nbra,"int",1,hlh->nvars-1," %% %% nbra: ");}
  tab=0; for (nv=0;nv<hlh->nvars-1;nv++){ tab += nbra[nv]*pow(hlh->nbins,nv);}
  if (verbose){ printf(" %% %% tab %d out of %d\n",tab,(int)pow(hlh->nbins,hlh->nvars-1));}
  if (hlh->data!=NULL){
    for (nv=0;nv<hlh->nvars;nv++){
      for (nc=0;nc<hlh->ncumulants;nc++){
	if (verbose){ printf(" %% %% adding %f to var %d cumulant %d...",pow(valra_end[nv],nc),nv,nc);}
	hlh->data[tab+(nv+nc*hlh->nvars)*(int)pow(hlh->nbins,hlh->nvars-1)] += pow(valra_end[nv],nc);
	if (verbose){ printf(" to obtain %f\n",hlh->data[tab+(nv+nc*hlh->nvars)*(int)pow(hlh->nbins,hlh->nvars-1)]);}}}}
  if (verbose>1){ 
    for (nc=0;nc<hlh->ncumulants;nc++){
      printf(" %% hlh->data cumulant %d\n",nc);
      raprintf(&(hlh->data[nc*(int)pow(hlh->nbins,hlh->nvars-1)*hlh->nvars]),"double",(int)pow(hlh->nbins,hlh->nvars-1),hlh->nvars,"hlh->data: ");}}
  tfree(nbra);nbra=NULL;
  if (verbose){ printf(" %% [exiting hhlib_histupdate]\n");}
}

void hhlib_project_printf(struct hhlib *hl,int nt,int nc)
{
  hhlib_histproject(hl->hlhra[nt]);
  hhlib_histproject_printf(hl->hlhra[nt],nc);
}

void hhlib_histproject(struct hhlib_hist *hlh)
{
  int verbose=0;
  int nv=0,nc=0,nv2=0,tab=0,*nbra=NULL;
  double val=0,val2=0;
  nbra = (int *) tcalloc(hlh->nvars-1,sizeof(int));
  rareset(hlh->data_project,"double",hlh->nbins*(hlh->nvars-1)*hlh->nvars*hlh->ncumulants,NULL);
  for (tab=0;tab<(int)pow(hlh->nbins,hlh->nvars-1);tab++){
    for (nv=0;nv<hlh->nvars-1;nv++){ nbra[nv] = (tab/(int)pow(hlh->nbins,nv))%hlh->nbins;}
    if (verbose){ printf("%% tab %d\n",tab); raprintf(nbra,"int",1,hlh->nvars-1," %% nbra:");}
    for (nv=0;nv<hlh->nvars-1;nv++){ for (nc=0;nc<hlh->ncumulants;nc++){ for (nv2=0;nv2<hlh->nvars;nv2++){
      if (verbose){ 
	val = hlh->data[tab+(nv2+nc*hlh->nvars)*(int)pow(hlh->nbins,hlh->nvars-1)];
	val2 = hlh->data_project[nbra[nv] + nv*hlh->nbins + (nv2+nc*hlh->nvars)*(hlh->nvars-1)*hlh->nbins];
	if (val!=0){ printf(" %% nv %d, nc %d, nv2 %d, adding %f to %f\n",nv,nc,nv2,val,val2);}}
      hlh->data_project[nbra[nv] + nv*hlh->nbins + (nv2+nc*hlh->nvars)*(hlh->nvars-1)*hlh->nbins] += hlh->data[tab+(nv2+nc*hlh->nvars)*(int)pow(hlh->nbins,hlh->nvars-1)];}}}}
  tfree(nbra);nbra=NULL;
}

void hhlib_histproject_printf(struct hhlib_hist *hlh,int nc)
{
  int nv=0,nv2=0;
  for (nv=0;nv<hlh->nvars-1;nv++){ for (nv2=0;nv2<(nc>0?hlh->nvars:1);nv2++){
    printf(" %% projecting variable %d, cumulant %d, output variable %d\n",nv,nc,nv2);
    raprintf(&(hlh->data_project[0 + nv*hlh->nbins + (nv2+nc*hlh->nvars)*(hlh->nvars-1)*hlh->nbins]),"double",1,hlh->nbins," %% hlh->data_project:");}}
}

double hhlib_histN(struct hhlib_hist *hlh,double *valra_start)
{
  int nv=0,nv_temp=0,*nbra=NULL,tab=0;
  double v=0;
  nbra = (int *) tcalloc(hlh->nvars-1,sizeof(int));
  nv_temp=0;
  for (nv=0;nv<hlh->nvars;nv++){
    if (nv!=hlh->indexing_trigger){
      v = hlh->logflagra[nv] ? log(valra_start[nv]) : valra_start[nv];
      nbra[nv_temp]=maximum(0,minimum(hlh->nbins-1,(int)floor(hlh->nbins*(v-hlh->minra[nv])/(hlh->maxra[nv]-hlh->minra[nv]))));
      nv_temp+=1;}}
  tab=0; for (nv=0;nv<hlh->nvars-1;nv++){ tab += nbra[nv]*pow(hlh->nbins,nv);}
  v = hlh->data[tab+(0+0*hlh->nvars)*(int)pow(hlh->nbins,hlh->nvars-1)];
  tfree(nbra);nbra=NULL;
  return v;
}

void hhlib_histmean(struct hhlib_hist *hlh,double *valra_start,double *valra_end)
{
  int verbose=0;
  int nv=0,nv_temp=0,*nbra=NULL,tab=0;
  double v=0;
  if (verbose){ printf(" %% [entering hhlib_histmean]\n"); raprintf(valra_start,"double",1,hlh->nvars," %% valra_start: ");}
  nbra = (int *) tcalloc(hlh->nvars-1,sizeof(int));
  nv_temp=0;
  for (nv=0;nv<hlh->nvars;nv++){
    if (nv!=hlh->indexing_trigger){
      v = hlh->logflagra[nv] ? log(valra_start[nv]) : valra_start[nv];
      nbra[nv_temp]=maximum(0,minimum(hlh->nbins-1,(int)floor(hlh->nbins*(v-hlh->minra[nv])/(hlh->maxra[nv]-hlh->minra[nv]))));
      nv_temp+=1;}}
  tab=0; for (nv=0;nv<hlh->nvars-1;nv++){ tab += nbra[nv]*pow(hlh->nbins,nv);}
  if (verbose){ printf(" %% tab %d\n",tab); raprintf(nbra,"int",1,hlh->nvars-1," %% nbra: ");}
  if (hlh->data!=NULL){
    for (nv=0;nv<hlh->nvars;nv++){
      if (hlh->libuseflagra[nv]){
	if (verbose){ printf(" %% using variable %d, sum %f, number %f\n",nv,hlh->data[tab+(nv+1*hlh->nvars)*(int)pow(hlh->nbins,hlh->nvars-1)],hlh->data[tab+(nv+0*hlh->nvars)*(int)pow(hlh->nbins,hlh->nvars-1)]);}
	valra_end[nv] = hlh->data[tab+(nv+1*hlh->nvars)*(int)pow(hlh->nbins,hlh->nvars-1)]/maximum(1,hlh->data[tab+(nv+0*hlh->nvars)*(int)pow(hlh->nbins,hlh->nvars-1)]);}}}
  if (verbose){ raprintf(valra_end,"double",1,hlh->nvars," %% valra_end: "); printf(" %% [finishing hhlib_histmean]\n"); }
  tfree(nbra);nbra=NULL;
}

void hhlib_histtfree(struct hhlib_hist *hlh)
{ 
  tfree(hlh->libuseflagra); hlh->libuseflagra=NULL;
  tfree(hlh->logflagra); hlh->logflagra=NULL;
  tfree(hlh->maxra); hlh->maxra=NULL;
  tfree(hlh->minra); hlh->minra=NULL;
  if (hlh->data!=NULL){ tfree(hlh->data); hlh->data=NULL;} 
  if (hlh->data_project!=NULL){ tfree(hlh->data_project); hlh->data_project=NULL;}
  tfree(hlh); hlh=NULL;
}

void hhlib_histdump(struct hhlib_hist *hlh,FILE *fp)
{
  fwrite(&(hlh->nvars),sizeof(int),1,fp);
  fwrite(&(hlh->indexing_trigger),sizeof(int),1,fp);
  fwrite(hlh->libuseflagra,sizeof(int),hlh->nvars,fp);
  fwrite(hlh->logflagra,sizeof(int),hlh->nvars,fp);
  fwrite(hlh->maxra,sizeof(double),hlh->nvars,fp);
  fwrite(hlh->minra,sizeof(double),hlh->nvars,fp);
  fwrite(&(hlh->nbins),sizeof(int),1,fp);
  fwrite(&(hlh->ncumulants),sizeof(int),1,fp);
  fwrite(hlh->data,sizeof(double),pow(hlh->nbins,hlh->nvars-1)*hlh->nvars*hlh->ncumulants,fp);
  fwrite(hlh->data_project,sizeof(double),hlh->nbins*(hlh->nvars-1)*hlh->nvars*hlh->ncumulants,fp);
}

struct hhlib_hist *hhlib_histread(FILE *fp)
{
  int verbose=0;
  int no_error=1,tab=0;
  struct hhlib_hist *hlh=NULL;
  if (verbose){ printf(" %% [entering hhlib_histread] \n");}
  hlh = (struct hhlib_hist *) tcalloc(1,sizeof(struct hhlib_hist));
  if (no_error){ no_error = (fread(&(hlh->nvars),sizeof(int),1,fp)==(size_t)1);}
  if (no_error){ no_error = (fread(&(hlh->indexing_trigger),sizeof(int),1,fp)==(size_t)1);}
  hlh->libuseflagra = (int *) tcalloc(hlh->nvars,sizeof(int));
  if (no_error){ no_error = (fread(&(hlh->libuseflagra),sizeof(int),hlh->nvars,fp)==(size_t)(hlh->nvars));}
  hlh->logflagra = (int *) tcalloc(hlh->nvars,sizeof(int));
  if (no_error){ no_error = (fread(&(hlh->logflagra),sizeof(int),hlh->nvars,fp)==(size_t)(hlh->nvars));}
  hlh->maxra = (double *) tcalloc(hlh->nvars,sizeof(double)); 
  if (no_error){ no_error = (fread(&(hlh->maxra),sizeof(double),hlh->nvars,fp)==(size_t)(hlh->nvars));}
  hlh->minra = (double *) tcalloc(hlh->nvars,sizeof(double)); 
  if (no_error){ no_error = (fread(&(hlh->minra),sizeof(double),hlh->nvars,fp)==(size_t)(hlh->nvars));}
  if (no_error){ no_error = (fread(&(hlh->nbins),sizeof(int),1,fp)==(size_t)1);}
  if (no_error){ no_error = (fread(&(hlh->ncumulants),sizeof(int),1,fp)==(size_t)1);}
  hlh->data = (double *) tcalloc(pow(hlh->nbins,hlh->nvars-1)*hlh->nvars*hlh->ncumulants,sizeof(double));
  if (no_error){ 
    tab=pow(hlh->nbins,hlh->nvars-1)*hlh->nvars*hlh->ncumulants; no_error = (fread(hlh->data,sizeof(double),tab,fp)==(size_t)tab);}
  hlh->data_project = (double *) tcalloc(hlh->nbins*(hlh->nvars-1)*hlh->nvars*hlh->ncumulants,sizeof(double));
  if (no_error){ 
    tab=hlh->nbins*(hlh->nvars-1)*hlh->nvars*hlh->ncumulants; no_error = (fread(hlh->data_project,sizeof(double),tab,fp)==(size_t)tab);}
  if (!no_error){ printf(" %% Warning! improperly read hhlib_hist in hhlib_histread\n");}
  return hlh;
}

void hhlibdump(struct hhlib *hl,char *filename)
{
  FILE *fp;
  int nt=0;
  if ((fp=fopen(filename,"w"))!=NULL){ printf(" %% warning, could not open %s in hhlibdump\n",filename); fp=stdout;}
  fwrite(&(hl->tau_skip),sizeof(double),1,fp);
  fwrite(&(hl->trigger_value),sizeof(double),1,fp);
  fwrite(&(hl->indexing_trigger),sizeof(int),1,fp);
  fwrite(&(hl->indexing_nvar_length),sizeof(int),1,fp);
  fwrite(hl->indexing_nvar_checkout,sizeof(int),hl->indexing_nvar_length,fp);
  fwrite(&(hl->indexing_nvar_refile_length),sizeof(int),1,fp);
  fwrite(hl->indexing_nvar_refile,sizeof(int),hl->indexing_nvar_refile_length,fp);
  fwrite(hl->logflagra,sizeof(int),hl->indexing_nvar_refile_length,fp);
  fwrite(hl->libuseflagra,sizeof(int),hl->indexing_nvar_refile_length,fp);
  fwrite(hl->maxra,sizeof(double),hl->indexing_nvar_refile_length,fp);
  fwrite(hl->minra,sizeof(double),hl->indexing_nvar_refile_length,fp);
  fwrite(&(hl->nbins),sizeof(int),1,fp);
  for (nt=0;nt<hl->Nra->ntypes;nt++){ hhlib_histdump(hl->hlhra[nt],fp);}
  fwrite(&(hl->use_after),sizeof(int),1,fp);
  if (fp!=stdout){ fclose(fp);fp=NULL;}
}

struct hhlib * hhlibread(struct neuronarray *Nra,char *filename)
{
  int verbose=0;
  struct hhlib *hl=NULL;
  FILE *fp=NULL;
  int nt=0,ni=0,tab,no_error=1;
  if (verbose){ printf(" %% [entering hhlibread] \n");}
  if ((fp=fopen(filename,"r"))!=NULL){ printf(" %% warning, could not open %s in hhlibread\n",filename); fp=stdin;}
  hl = (struct hhlib *) tcalloc(1,sizeof(struct hhlib));
  hl->Nra = Nra;
  if (no_error){ tab=1; no_error = (fread(&(hl->tau_skip),sizeof(double),tab,fp)==(size_t)tab);}
  if (no_error){ tab=1; no_error = (fread(&(hl->trigger_value),sizeof(double),tab,fp)==(size_t)tab);}
  if (no_error){ tab=1; no_error = (fread(&(hl->indexing_trigger),sizeof(int),tab,fp)==(size_t)tab);}
  if (no_error){ tab=1; no_error = (fread(&(hl->indexing_nvar_length),sizeof(int),tab,fp)==(size_t)tab);}
  hl->indexing_nvar_checkout = (int *) tcalloc(hl->indexing_nvar_length,sizeof(int));
  if (no_error){ tab=hl->indexing_nvar_length; no_error = (fread(hl->indexing_nvar_checkout,sizeof(int),tab,fp)==(size_t)tab);}
  if (no_error){ tab=1; no_error = (fread(&(hl->indexing_nvar_refile_length),sizeof(int),tab,fp)==(size_t)tab);}
  if (no_error){ tab=hl->indexing_nvar_refile_length; no_error = (fread(hl->indexing_nvar_refile,sizeof(int),tab,fp)==(size_t)tab);}
  hl->logflagra = (int *) tcalloc(hl->indexing_nvar_refile_length,sizeof(int));
  hl->libuseflagra = (int *) tcalloc(hl->indexing_nvar_refile_length,sizeof(int));
  hl->maxra= (double *) tcalloc(hl->indexing_nvar_refile_length,sizeof(double));
  hl->minra= (double *) tcalloc(hl->indexing_nvar_refile_length,sizeof(double));
  if (no_error){ tab=hl->indexing_nvar_refile_length; no_error = (fread(hl->logflagra,sizeof(int),tab,fp)==(size_t)tab);}
  if (no_error){ tab=hl->indexing_nvar_refile_length; no_error = (fread(hl->libuseflagra,sizeof(int),tab,fp)==(size_t)tab);}
  if (no_error){ tab=hl->indexing_nvar_refile_length; no_error = (fread(hl->maxra,sizeof(double),tab,fp)==(size_t)tab);}
  if (no_error){ tab=hl->indexing_nvar_refile_length; no_error = (fread(hl->minra,sizeof(double),tab,fp)==(size_t)tab);}
  if (no_error){ tab=1; no_error = (fread(&(hl->nbins),sizeof(int),tab,fp)==(size_t)tab);}
  hl->trigger_rarara = (double ***) tcalloc(hl->Nra->ntypes,sizeof(double **));
  for (nt=0;nt<hl->Nra->ntypes;nt++){
    hl->trigger_rarara[nt] = (double **) tcalloc(hl->Nra->lengthra[nt],sizeof(double *));
    for (ni=0;ni<hl->Nra->lengthra[nt];ni++){
      hl->trigger_rarara[nt][ni] = (double *) tcalloc(hl->indexing_nvar_length+2,sizeof(double));}}
  hl->hlhra = (struct hhlib_hist **) tcalloc(hl->Nra->ntypes,sizeof(struct hhlib_hist *));
  for (nt=0;nt<hl->Nra->ntypes;nt++){ hl->hlhra[nt] = hhlib_histread(fp);}
  if (no_error){ tab=1; no_error = (fread(&(hl->use_after),sizeof(int),tab,fp)==(size_t)tab);}
  if (!no_error){ printf(" %% warning! error reading hhlib in hhlibread\n");}
  if (fp!=stdin){ fclose(fp);fp=NULL;}
  return hl;
}

struct hhlib * hhlibmake(struct neuronarray *Nra,double tau_skip,int indexing_trigger,int indexing_nvar_length,int *indexing_nvar_checkout,int indexing_nvar_refile_length,int *indexing_nvar_refile,int *logflagra,int *libuseflagra,double *maxra,double *minra,int nbins,int use_after)
{
  /* hl->trigger_rarara[nt][ni][2+?] stores variables of interest (including trigger variable), then stores:
     hl->trigger_rarara[nt][ni][0] = refracted_flag (0 normal, 1 waiting to read, 2 waiting to write)
     hl->trigger_rarara[nt][ni][1] = time_last
     hl->hlhra[nt] stores cumulants of variables of interest:
  */
  struct hhlib *hl=NULL;
  int nt=0,ni=0,nv=0;
  int *logflagra_temp=NULL,*libuseflagra_temp=NULL;
  double *maxra_temp=NULL,*minra_temp=NULL;
  hl = (struct hhlib *) tcalloc(1,sizeof(struct hhlib));
  hl->Nra = Nra;
  hl->tau_skip = tau_skip;
  hl->trigger_value = VOLTAGE_THRESHOLD_S;
  hl->indexing_trigger = indexing_trigger;
  hl->indexing_nvar_length = indexing_nvar_length;
  hl->indexing_nvar_checkout = (int *) tcalloc(hl->indexing_nvar_length,sizeof(int));
  for (nv=0;nv<hl->indexing_nvar_length;nv++){ hl->indexing_nvar_checkout[nv] = indexing_nvar_checkout[nv];}
  hl->indexing_nvar_refile_length = indexing_nvar_refile_length;
  hl->indexing_nvar_refile = (int *) tcalloc(hl->indexing_nvar_refile_length,sizeof(int));
  for (nv=0;nv<hl->indexing_nvar_refile_length;nv++){ hl->indexing_nvar_refile[nv] = indexing_nvar_refile[nv];}
  hl->logflagra = (int *) tcalloc(hl->indexing_nvar_refile_length,sizeof(int));
  hl->libuseflagra = (int *) tcalloc(hl->indexing_nvar_refile_length,sizeof(int));
  hl->maxra= (double *) tcalloc(hl->indexing_nvar_refile_length,sizeof(double));
  hl->minra= (double *) tcalloc(hl->indexing_nvar_refile_length,sizeof(double));
  for (nv=0;nv<indexing_nvar_refile_length;nv++){ hl->logflagra[nv] = logflagra[nv]; hl->libuseflagra[nv] = libuseflagra[nv]; hl->maxra[nv] = maxra[nv]; hl->minra[nv] = minra[nv];}
  hl->nbins = nbins;
  hl->trigger_rarara = (double ***) tcalloc(hl->Nra->ntypes,sizeof(double **));
  for (nt=0;nt<hl->Nra->ntypes;nt++){
    hl->trigger_rarara[nt] = (double **) tcalloc(hl->Nra->lengthra[nt],sizeof(double *));
    for (ni=0;ni<hl->Nra->lengthra[nt];ni++){
      hl->trigger_rarara[nt][ni] = (double *) tcalloc(hl->indexing_nvar_length+2,sizeof(double));}}
  hl->hlhra = (struct hhlib_hist **) tcalloc(hl->Nra->ntypes,sizeof(struct hhlib_hist *));
  maxra_temp = (double *) tcalloc(hl->indexing_nvar_length,sizeof(double));
  minra_temp = (double *) tcalloc(hl->indexing_nvar_length,sizeof(double));
  logflagra_temp = (int *) tcalloc(hl->indexing_nvar_length,sizeof(int));
  libuseflagra_temp = (int *) tcalloc(hl->indexing_nvar_length,sizeof(int));
  for (nv=0;nv<hl->indexing_nvar_length;nv++){ 
    maxra_temp[nv] = hl->maxra[hl->indexing_nvar_checkout[nv]];
    minra_temp[nv] = hl->minra[hl->indexing_nvar_checkout[nv]];
    logflagra_temp[nv] = hl->logflagra[hl->indexing_nvar_checkout[nv]];
    libuseflagra_temp[nv] = hl->libuseflagra[hl->indexing_nvar_checkout[nv]];}
  for (nt=0;nt<hl->Nra->ntypes;nt++){
    hl->hlhra[nt] = hhlib_histmake(hl->indexing_nvar_length,hl->indexing_nvar_refile[hl->indexing_trigger],libuseflagra_temp,logflagra_temp,maxra_temp,minra_temp,hl->nbins,3);}
  hl->use_after = use_after;
  tfree(libuseflagra_temp);libuseflagra_temp=NULL;
  tfree(logflagra_temp);logflagra_temp=NULL;
  tfree(maxra_temp);maxra_temp=NULL;
  tfree(minra_temp);minra_temp=NULL;
  return hl;
}

void hhlibupdate(struct hhlib *hl,double t,double DT)
{
  /* assumes that n->vra[] holds the previous values of *(n->vpra[]) */
  int verbose=0;
  int nt=0,ni=0,nv=0,tab=0;
  struct neuron *n=NULL;
  double v_trigger_new=0,v_trigger_old=0,t_temp=0,dt_temp=0;
  double *vra=NULL,*vra_temp=NULL;
  if (verbose>1){ printf(" %% [entering hhlibupdate] t %f DT %f\n",t,DT);}
  vra = (double *) tcalloc(hl->indexing_nvar_length,sizeof(double));
  for (nt=0;nt<hl->Nra->ntypes;nt++){ for (ni=0;ni<hl->Nra->lengthra[nt];ni++){ 
    n=nget(hl->Nra,nt,ni);
    if (hl->trigger_rarara[nt][ni][0]==0){
      v_trigger_old = n->vra[hl->indexing_trigger];
      v_trigger_new = *(n->vpra[hl->indexing_trigger]);
      if (verbose>2){ printf(" %% at neuron (%d,%d) with r_flag %d and %s=[%f,%f]\n",nt,ni,(int)hl->trigger_rarara[nt][ni][0],GLOBAL_VARNAMES[hl->indexing_trigger],v_trigger_old,v_trigger_new);}
      if (v_trigger_new>=hl->trigger_value && v_trigger_old<hl->trigger_value){
	dt_temp = linerootfinder(v_trigger_old,v_trigger_new,hl->trigger_value,DT);
	if (verbose){ printf(" %% found neuron (%d,%d) with %s=[%f,%f], crossed at time %f+%f... ",nt,ni,GLOBAL_VARNAMES[hl->indexing_trigger],v_trigger_old,v_trigger_new,t,dt_temp);}
	for (nv=0;nv<hl->indexing_nvar_length;nv++){
	  tab = hl->indexing_nvar_checkout[nv];
	  hl->trigger_rarara[nt][ni][2+nv] = n->vra[tab]+(dt_temp/DT)*(*(n->vpra[tab]) - n->vra[tab]);}
	if (verbose){ raprintf(&(hl->trigger_rarara[nt][ni][2]),"double",1,hl->indexing_nvar_length," %% %% var:");}
	if (verbose){ printf(" %% already recorded %d for this bin\n",(int)hhlib_histN(hl->hlhra[nt],&(hl->trigger_rarara[nt][ni][2])));}
	if (verbose){ 
	  printf(" %% current mean for this bin is:\n");
	  vra_temp = (double *) tcalloc(hl->indexing_nvar_length,sizeof(double));
	  hhlib_histmean(hl->hlhra[nt],&(hl->trigger_rarara[nt][ni][2]),vra_temp);
	  raprintf(vra_temp,"double",1,hl->indexing_nvar_length," %% vra_temp: ");
	  tfree(vra_temp);vra_temp=NULL;}
	hl->trigger_rarara[nt][ni][0] = (hhlib_histN(hl->hlhra[nt],&(hl->trigger_rarara[nt][ni][2]))>hl->use_after?2:1);
	hl->trigger_rarara[nt][ni][1] = t+dt_temp;
	if (verbose){ printf(" set r_flag %d, time %0.2f\n",(int)hl->trigger_rarara[nt][ni][0],hl->trigger_rarara[nt][ni][1]);}}
      else if (v_trigger_new<hl->trigger_value && v_trigger_old<hl->trigger_value){
	hl->trigger_rarara[nt][ni][2+hl->indexing_nvar_refile[hl->indexing_trigger]]=v_trigger_new;}}
    else if (hl->trigger_rarara[nt][ni][0]==1){
      if (verbose>1){ printf(" %% at neuron (%d,%d) with r_flag %d\n",nt,ni,(int)hl->trigger_rarara[nt][ni][0]);}
      t_temp = hl->trigger_rarara[nt][ni][1]+hl->tau_skip;
      if (t<t_temp && (t+DT)>=t_temp){
	dt_temp = t_temp-t;
	if (verbose){ printf(" %% finished recording neuron (%d,%d) at time %f in [%f,%f]\n",nt,ni,t_temp,t,t+DT);}
	for (nv=0;nv<hl->indexing_nvar_length;nv++){
	  tab = hl->indexing_nvar_checkout[nv];
	  vra[nv] = n->vra[tab]+(dt_temp/DT)*(*(n->vpra[tab]) - n->vra[tab]);}
	if (verbose){ raprintf(vra,"double",1,hl->indexing_nvar_length," %% %% var:");}
	hhlib_histupdate(hl->hlhra[nt],&(hl->trigger_rarara[nt][ni][2]),vra);
	if (verbose){ 
	  printf(" %% just checking to ensure hhlib_histupdate worked correctly:\n");
	  vra_temp = (double *) tcalloc(hl->indexing_nvar_length,sizeof(double));
	  hhlib_histmean(hl->hlhra[nt],&(hl->trigger_rarara[nt][ni][2]),vra_temp);
	  raprintf(vra_temp,"double",1,hl->indexing_nvar_length," %% vra_temp: ");
	  tfree(vra_temp);vra_temp=NULL;}
	hl->trigger_rarara[nt][ni][0]=0;
	hl->trigger_rarara[nt][ni][2+hl->indexing_nvar_refile[hl->indexing_trigger]]=*(n->vpra[hl->indexing_trigger]);
	if (verbose){ printf(" %% resetting trigger value to %f\n",hl->trigger_rarara[nt][ni][2+hl->indexing_nvar_refile[hl->indexing_trigger]]);}}}
    else if (hl->trigger_rarara[nt][ni][0]==2){
      if (verbose>1){ printf(" %% at neuron (%d,%d) with r_flag %d\n",nt,ni,(int)hl->trigger_rarara[nt][ni][0]);}
      t_temp = hl->trigger_rarara[nt][ni][1]+hl->tau_skip;
      if (t+DT>=t_temp){ 
	if (verbose){ printf(" %% used library\n");}
	hl->trigger_rarara[n->type][n->index][0]=0;
	hl->trigger_rarara[nt][ni][2+hl->indexing_nvar_refile[hl->indexing_trigger]] = n->vra[hl->indexing_trigger];}}}}
  tfree(vra);vra=NULL;
  if (verbose>1){ printf(" %% [exiting hhlibupdate]\n");}
}

void hhlibprintf(struct hhlib *hl)
{
  int nt=0,nc=0;
  struct hhlib_hist *hlh=NULL;
  for (nt=0;nt<hl->Nra->ntypes;nt++){
    hlh = hl->hlhra[nt];
    printf(" %% type %d\n",nt);
    raprintf(hlh->libuseflagra,"int",1,hlh->nvars,"libuseflagra: ");
    raprintf(hlh->logflagra,"int",1,hlh->nvars,"logflagra: ");
    raprintf(hlh->maxra,"double",1,hlh->nvars,"maxra: ");
    raprintf(hlh->minra,"double",1,hlh->nvars,"minra: ");
    for (nc=0;nc<hlh->ncumulants;nc++){
      printf(" %% cumulant %d\n",nc);
      raprintf(&(hlh->data[(int)pow(hlh->nbins,hlh->nvars)*nc]),"double",hlh->nbins,pow(hlh->nbins,hlh->nvars-1),"data: ");}}
}

/* Here are the suite functions */

void multidra2jpg(double *dra,int ndim,int *dimra,char *filename_base,int *dim_to_plot,double max,double min)
{
  /* assumes that filename_base contains "./" or equivalent */
  /* (*int) dim_to_plot is an array of size ndim, with code
     0 running variable (plotted along horizontal axis, must be unique)
     1 overlapping variable
     2 nonoverlapping variable (not yet implemented)
  */
  int verbose=0;
  int indexmax=0,nd=0,index=0,*P=NULL,*dimra_P=NULL,*indexra=NULL,*indexra_P=NULL,index_P=0,index_R=0,indexmax_R=0;
  int running_variable=0;
  char filename[4096],command[5120];
  FILE *fp=NULL;  
  int remove_flag=0,jpeg_flag=(ON_MY_COMPUTER ? 1 : 0),eps_flag=(ON_MY_COMPUTER ? 1 : 0);
  double maxdia = 10000;
  double xpos=0,ypos=0;
  double rcolor=0,gcolor=0,bcolor=0;
  int colorcode=0,fontsize=0;
  if (verbose){ printf(" %% [entering multidra2jpg]\n");}
  running_variable=-1; indexmax=1; 
  for (nd=0;nd<ndim;nd++){ indexmax *= dimra[nd]; if (dim_to_plot[nd]==0){ running_variable=dim_to_plot[nd];}}
  if (max<=min){ stats("double",dra,indexmax,&max,&min,NULL,NULL); if (max<=min){ max=min+1;}}
  if (running_variable==-1){ printf(" %% Warning! no running_variable in multidra2jpg\n");}
  if (verbose){ printf(" %% indexmax %d, running_variable %d\n",indexmax,running_variable);}
  P = (int *) tcalloc(ndim,sizeof(int));
  for (nd=0;nd<ndim;nd++){ 
    if (nd==running_variable){ P[0]=nd;} 
    else if (nd<running_variable){ P[1+nd]=nd;}
    else if (nd>running_variable){ P[nd]=nd;}}
  if (verbose){ raprintf(P,"int",1,ndim,"P: ");}
  dimra_P = (int *) tcalloc(ndim,sizeof(int)); indexmax_R=1;
  for (nd=0;nd<ndim;nd++){ dimra_P[nd] = dimra[P[nd]]; if (nd>0){ indexmax_R *= dimra_P[nd];}}
  if (verbose){ 
    printf(" %% indexmax_R=%d\n",indexmax_R);
    raprintf(dimra,"int",1,ndim,"dimra: ");
    raprintf(dimra_P,"int",1,ndim,"dimra_P: ");}
  sprintf(filename,"%s.fig",filename_base);
  if ((fp=fopen(filename,"w"))==NULL){ printf(" %% warning! cannot open %s in multidra2jpg\n",filename); fp=stdout;}
  fprintf(fp,"%s",FIG_PREAMBLE); fprintf(fp,"%s",FIG_PREAMBLE_COLOR_7);
  indexra = (int *) tcalloc(ndim,sizeof(int));
  indexra_P = (int *) tcalloc(ndim,sizeof(int));
  index_R=0;colorcode=0;
  for (index_P=0;index_P<indexmax;index_P++){
    indextract(index_P,ndim,dimra_P,indexra_P);
    for (nd=0;nd<ndim;nd++){ indexra[P[nd]] = indexra_P[nd];}
    index = indexpack(ndim,dimra,indexra);
    if (verbose>1){ 
      printf(" %% index_P %d, \n",index_P); raprintf(indexra_P,"int",1,ndim,"indexra_P: ");
      printf(" %% sent to:\n");
      printf(" %% index %d, \n",index); raprintf(indexra,"int",1,ndim,"indexra: ");}
    if (indexra_P[0]==0){ /* line_start */
      if (ndim>1){ 
	index_R = indexpack(ndim-1,&(dimra_P[1]),&(indexra_P[1]));
	colorscale(0,index_R,indexmax_R-1,0,&rcolor,&gcolor,&bcolor);
	colorcode = crop((int)floor(512*rcolor),0,511);}
      fprintf(fp,"2 1 0 %d %d 7 %d 0 -1 0.000 0 0 -1 0 0 %d\n\t",/*width*/1,/*color*/32+colorcode,/*depth*/(index_R%999)+1,/*npoints*/dimra_P[0]);}
    xpos = (double)(indexra_P[0])/(double)maximum(1,dimra_P[0]-1);
    ypos = (dra[index]-min)/(max-min);
    fprintf(fp," %d %d",(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos));
    if (indexra_P[0]==dimra_P[0]-1){ /* line end */ fprintf(fp,"\n");} }
  fprintf(fp,"2 1 0 %d %d 7 %d 0 -1 0.000 0 0 -1 0 0 %d\n\t",/*width*/1,/*color*/0,/*depth*/999,/*npoints*/5);
  xpos = 0; ypos = 0; fprintf(fp," %d %d",(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos));
  xpos = 1; ypos = 0; fprintf(fp," %d %d",(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos));
  xpos = 1; ypos = 1; fprintf(fp," %d %d",(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos));
  xpos = 0; ypos = 1; fprintf(fp," %d %d",(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos));
  xpos = 0; ypos = 0; fprintf(fp," %d %d",(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos));
  fprintf(fp,"\n");
  xpos = -0.1; ypos = 0.0; fontsize=24;
  fprintf(fp,"4 0 0 %d 0 %d %d 1.5708 4 %d %d %d %d %0.16f\\001\n",/*depth*/1,/*font*/0,/*point*/fontsize,/*textheight*/(int)(fontsize*13.5),/*textwidth*/(int)(fontsize*9*18),(int)(maxdia*xpos),(int)maxdia-(int)(maxdia*ypos),min);
  xpos = -0.1; ypos = 1.0; fontsize=24;
  fprintf(fp,"4 0 0 %d 0 %d %d 1.5708 4 %d %d %d %d %0.16f\\001\n",/*depth*/1,/*font*/0,/*point*/fontsize,/*textheight*/(int)(fontsize*13.5),/*textwidth*/(int)(fontsize*9*18),(int)(maxdia*xpos),(int)maxdia-(int)(maxdia*ypos),max);
  xpos = 0.0; ypos = -0.1; fontsize=24;
  fprintf(fp,"4 0 0 %d 0 %d %d 0.0000 4 %d %d %d %d %d\\001\n",/*depth*/1,/*font*/0,/*point*/fontsize,/*textheight*/(int)(fontsize*13.5),/*textwidth*/(int)(fontsize*9*18),(int)(maxdia*xpos),(int)maxdia-(int)(maxdia*ypos),0);
  xpos = 1.0; ypos = -0.1; fontsize=24;
  fprintf(fp,"4 0 0 %d 0 %d %d 0.0000 4 %d %d %d %d %d\\001\n",/*depth*/1,/*font*/0,/*point*/fontsize,/*textheight*/(int)(fontsize*13.5),/*textwidth*/(int)(fontsize*9*18),(int)(maxdia*xpos),(int)maxdia-(int)(maxdia*ypos),dimra_P[0]-1);
  xpos = 0.0; ypos = 1.1; fontsize=8;
  fprintf(fp,"4 0 0 %d 0 %d %d 0.0000 4 %d %d %d %d %s\\001\n",/*depth*/1,/*font*/0,/*point*/fontsize,/*textheight*/(int)(fontsize*13.5),/*textwidth*/(int)(fontsize*9*strlen(filename_base)),(int)(maxdia*xpos),(int)maxdia-(int)(maxdia*ypos),filename_base);
  if (fp!=stdout){ fclose(fp);fp=NULL;}
  if (jpeg_flag){ sprintf(command,"fig2dev -Ljpeg -q %d %s.fig %s.jpg;",/*quality*/5,filename_base,filename_base); system(command);}
  if (eps_flag){ sprintf(command,"fig2dev -Leps %s.fig %s.eps;",filename_base,filename_base); system(command);}
  if (remove_flag){ sprintf(command,"rm %s.fig;",filename_base); system(command);}
  tfree(P);P=NULL;
  tfree(dimra_P);dimra_P=NULL;
  tfree(indexra);indexra=NULL;
  tfree(indexra_P);indexra_P=NULL;
  if (verbose){ printf(" %% [finishing multidra2jpg]\n");}
}

void checktofind_howmany_ra_general_average(int dim,char **finra,struct llist **nameLra,int *dimra,int *found_type,int *found_rows,int *found_cols,int *dim_to_average_maximum,int *dim_to_average_minimum,double **output,int *output_type,int *output_rows,int *output_cols,int *output_nfound)
{
  /* averages along dimension nd from dim_to_average_minimum[nd] ot dim_to_average_maximum[nd] */
  int verbose=0;
  int nd=0,nd2=0,nr=0;
  struct litem *l=NULL;
  char filename[4096],tempchar[1024];
  int *dimra_local=NULL,*indexra=NULL;
  int index=0,indexmax=0,index_full=0;
  int local_type=0,local_rows=0,local_cols=0;
  int *ira=NULL;
  double *dra=NULL;
  for (nd=0;nd<dim;nd++){ if (dim_to_average_maximum[nd]<dim_to_average_minimum[nd] || dim_to_average_maximum[nd]<0 || dim_to_average_maximum[nd]>=dimra[nd] || dim_to_average_minimum[nd]<0 || dim_to_average_minimum[nd]>=dimra[nd]){ dim_to_average_maximum[nd] = dimra[nd]-1; dim_to_average_minimum[nd]=0;}}
  if (verbose){ raprintf(dim_to_average_minimum,"int",1,dim," %% dim_to_average_minimum: ");}
  if (verbose){ raprintf(dim_to_average_maximum,"int",1,dim," %% dim_to_average_maximum: ");}
  dimra_local = (int *) tcalloc(dim,sizeof(int)); indexmax=1;
  for (nd=0;nd<dim;nd++){ dimra_local[nd]=dim_to_average_maximum[nd]-dim_to_average_minimum[nd]+1; indexmax*=dimra_local[nd];}
  if (verbose){ raprintf(dimra_local,"int",1,dim," %% dimra_local: ");}
  indexra = (int *) tcalloc(dim,sizeof(int)); 
  index=0; indextract(index,dim,dimra_local,indexra); for (nd=0;nd<dim;nd++){ indexra[nd] += dim_to_average_minimum[nd];}
  index_full = indexpack(dim,dimra,indexra);
  *output_type = found_type[index_full]; *output_rows = found_rows[index_full]; *output_cols = found_cols[index_full];
  for (index=0;index<indexmax;index++){
    indextract(index,dim,dimra_local,indexra); for (nd=0;nd<dim;nd++){ indexra[nd] += dim_to_average_minimum[nd];}
    index_full = indexpack(dim,dimra,indexra);
    if (found_type[index_full]!=*output_type){ printf(" %% warning! type mismatch %d!=%din checktofind_howmany_ra_general_average\n",found_type[index_full],*output_type);}
    if (found_rows[index_full]!=*output_rows){ printf(" %% warning! rows mismatch %d!=%din checktofind_howmany_ra_general_average\n",found_rows[index_full],*output_rows);}
    if (found_cols[index_full]!=*output_cols){ printf(" %% warning! cols mismatch %d!=%din checktofind_howmany_ra_general_average\n",found_cols[index_full],*output_cols);}}
  if (verbose){ printf(" %% found output_type %d, output_rows %d, output_cols %d\n",*output_type,*output_rows,*output_cols);}
  *output = (double *) tcalloc(*output_rows**output_cols,sizeof(double));
  *output_nfound = 0;
  for (index=0;index<indexmax;index++){
    indextract(index,dim,dimra_local,indexra); for (nd=0;nd<dim;nd++){ indexra[nd] += dim_to_average_minimum[nd];}
    if (verbose){ printf(" %% index %d, ",index); raprintf(indexra,"int",1,dim,"indexra: ");}
    sprintf(filename,"%s",finra[0]);
    for (nd2=1;nd2<=dim;nd2++){ 
      if (nameLra[nd2-1]->length>=1+indexra[nd2-1]){ 
	l=nameLra[nd2-1]->first; nr=0; while (l!=NULL && nr<indexra[nd2-1]){ l=l->child; nr+=1;}
	sprintf(tempchar,"%s",(char *)l->item);}
      else /* if (nameLra[nd2-1]->length<1+indexra[nd2-1]) */{ sprintf(tempchar,"%d",indexra[nd2-1]);}
      sprintf(filename,"%s%s%s",filename,tempchar,finra[nd2]);}
    if (verbose){ printf(" %% checking for %s... ",filename);}
    if (checktofind(filename)){ 
      if (verbose){ printf("found\n");}
      *output_nfound += 1;
      switch (*output_type){
      case 0: /* int */
	ira = raread(filename,&local_type,&local_rows,&local_cols); for (nr=0;nr<*output_rows**output_cols;nr++){ (*output)[nr] += ira[nr];} tfree(ira);ira=NULL;
	break;
      case 1: /* double */
	dra = raread(filename,&local_type,&local_rows,&local_cols); for (nr=0;nr<*output_rows**output_cols;nr++){ (*output)[nr] += dra[nr];} tfree(dra);dra=NULL;
	break;
      default: printf(" %% Warning, improper type %d in checktofind_howmany_ra_general_average\n",*output_type);}}
    else /*if (!checktofind(filename))*/{ 
      if (verbose){ printf("warning! %s not found\n",filename);}}}
  if (verbose){ printf(" %% %d total records found\n",*output_nfound);}
  ratimesequals(*output,*output_rows**output_cols,1.0/(double)maximum(1,*output_nfound));
  tfree(dimra_local);dimra_local=NULL;
  tfree(indexra);indexra=NULL;
}

int checktofind_howmany_general(char *fpre,struct llist *nameL,int *index,char *fpos)
{
  /* checks for fpre%sfpos, with characters drawn from nameL, if the llist is too short, integers are used instead */
  /* limited with inefficient n^2 litem extraction */
  char filename[4096],tempchar[1024];
  int continue_flag=0;
  int nfound=0;
  struct litem *l=NULL;
  int nr=0;
  continue_flag=1;nfound=0;
  do{
    if (nameL->length>=(1+*index)){ 
      nr=0; l=nameL->first; while (l!=NULL && nr<*index){ l=l->child;nr+=1;} sprintf(tempchar,"%s",(char *)l->item);}
    else /* if (nameL->length<(1+*index)) */{ sprintf(tempchar,"%d",*index);}
    sprintf(filename,"%s%s%s",fpre,tempchar,fpos); 
    continue_flag=checktofind(filename);
    if (continue_flag){ *index += 1; nfound += 1;}}
  while (continue_flag);
  return nfound;
}

void checktofind_howmany_ra_general(int dim,char **finra,struct llist **nameLra,int **dimra,int **found_type,int **found_rows,int **found_cols)
{
  /* if finra = "pre","apple","orange","pos"
     then we search for pre%sapple%sorange%spos with bounds boxed in by the maximum achieved along coordinate axes.
     the jth set of characters used in the string are drawn from llist * nameLra[j].
     if the llist is too short, integers are used instead.
     Note that **finra should be of size dim+1. 
     We assume neither dimra nor found_type is initialized 
     Note that this is limited with inefficient n^2 litem extraction */
  int verbose=0;
  char filename1[4096],filename2[4096],tempchar[1024];
  int nd=0,nd2=0,indexmax=0,index=0,*indexra=NULL;
  struct litem *l=NULL;
  int nr=0;
  FILE *fp=NULL;
  if (verbose){ printf("%% [entering checktofind_howmany_ra_general] dim %d\n",dim);}
  *dimra = (int *) tcalloc(dim,sizeof(int));
  for (nd=1;nd<=dim;nd++){
    sprintf(filename1,"%s",finra[0]);
    for (nd2=1;nd2<nd;nd2++){ 
      if (nameLra[nd2-1]->length>=1){ sprintf(tempchar,"%s",(char *)(nameLra[nd2-1]->first->item));}
      else /* if (nameLra[nd2-1]->length==0) */{ sprintf(tempchar,"%s","0");}
      sprintf(filename1,"%s%s%s",filename1,tempchar,finra[nd2]);}
    sprintf(filename2,"%s",finra[nd]);
    for (nd2=nd+1;nd2<=dim;nd2++){ 
      if (nameLra[nd2-1]->length>=1){ sprintf(tempchar,"%s",(char *)(nameLra[nd2-1]->first->item));}
      else /* if (nameLra[nd2-1]->length==0) */{ sprintf(tempchar,"%s","0");}
      sprintf(filename2,"%s%s%s",filename2,tempchar,finra[nd2]);}
    if (verbose){ printf(" %% checking %s?%s...",filename1,filename2);}
    checktofind_howmany_general(filename1,nameLra[nd-1],&((*dimra)[nd-1]),filename2);
    if (verbose){ printf(" found %d entries\n",(*dimra)[nd-1]);}}
  if (verbose){ raprintf(*dimra,"int",1,dim,"dimra: ");}
  indexmax=1; for (nd=0;nd<dim;nd++){ indexmax *= (*dimra)[nd];}
  *found_type = (int *) tcalloc(indexmax,sizeof(int));
  *found_rows = (int *) tcalloc(indexmax,sizeof(int));
  *found_cols = (int *) tcalloc(indexmax,sizeof(int));
  indexra = (int *) tcalloc(dim,sizeof(int));
  for (index=0;index<indexmax;index++){
    indextract(index,dim,*dimra,indexra);
    if (verbose){ printf(" %% index %d, ",index); raprintf(indexra,"int",1,dim,"indexra: ");}
    sprintf(filename1,"%s",finra[0]);
    for (nd2=1;nd2<=dim;nd2++){ 
      if (nameLra[nd2-1]->length>=1+indexra[nd2-1]){ 
	l=nameLra[nd2-1]->first; nr=0; while (l!=NULL && nr<indexra[nd2-1]){ l=l->child; nr+=1;}
	sprintf(tempchar,"%s",(char *)l->item);}
      else /* if (nameLra[nd2-1]->length<1+indexra[nd2-1]) */{ sprintf(tempchar,"%d",indexra[nd2-1]);}
      sprintf(filename1,"%s%s%s",filename1,tempchar,finra[nd2]);}
    if (verbose){ printf(" %% checking for %s... ",filename1);}
    if (checktofind(filename1)){ 
      if (verbose){ printf("found\n");}
      if ((fp=fopen(filename1,"r"))==NULL){ printf(" %% warning! %s not found\n",filename1);}
      fread(&((*found_type)[index]),sizeof(int),1,fp);
      fread(&((*found_rows)[index]),sizeof(int),1,fp);
      fread(&((*found_cols)[index]),sizeof(int),1,fp);
      fclose(fp);fp=NULL;}
    else /*if (!checktofind(filename1))*/{ 
      if (verbose){ printf("not found\n");} 
      (*found_type)[index]=-1;
      (*found_rows)[index]=-1;
      (*found_cols)[index]=-1;}}
  tfree(indexra);indexra=NULL;
  if (verbose){ raprintf(*found_type,"int",1,indexmax,"found_type: ");}
  if (verbose){ raprintf(*found_rows,"int",1,indexmax,"found_rows: ");}
  if (verbose){ raprintf(*found_cols,"int",1,indexmax,"found_cols: ");}
  if (verbose){ printf(" %% [exiting checktofind_howmany_ra_general]\n");}
}

void checktofind_howmany_ra(int dim,char **finra,int **dimra,int **found_type)
{
  /* if finra = "pre","apple","orange","pos"
     then we search for pre%dapple%dorange%dpos with bounds boxed in by the maximum achieved along coordinate axes.
     Note that **finra should be of size dim+1. 
     We assume neither dimra nor found_type is initialized */
  int verbose=0;
  char filename1[4096],filename2[4096];
  int nd=0,nd2=0,indexmax=0,index=0,*indexra=NULL;
  if (verbose){ printf("%% [entering checktofind_howmany_ra] dim %d\n",dim);}
  *dimra = (int *) tcalloc(dim,sizeof(int));
  for (nd=1;nd<=dim;nd++){
    sprintf(filename1,"%s",finra[0]);
    for (nd2=1;nd2<nd;nd2++){ sprintf(filename1,"%s0%s",filename1,finra[nd2]);}
    sprintf(filename2,"%s",finra[nd]);
    for (nd2=nd+1;nd2<=dim;nd2++){ sprintf(filename2,"%s0%s",filename2,finra[nd2]);}
    if (verbose){ printf(" %% checking %s?%s...",filename1,filename2);}
    checktofind_howmany(filename1,&((*dimra)[nd-1]),filename2);
    if (verbose){ printf(" found %d entries\n",(*dimra)[nd-1]);}}
  if (verbose){ raprintf(*dimra,"int",1,dim,"dimra: ");}
  indexmax=1; for (nd=0;nd<dim;nd++){ indexmax *= (*dimra)[nd];}
  *found_type = (int *) tcalloc(indexmax,sizeof(int));
  indexra = (int *) tcalloc(dim,sizeof(int));
  for (index=0;index<indexmax;index++){
    indextract(index,dim,*dimra,indexra);
    sprintf(filename1,"%s",finra[0]);
    for (nd2=1;nd2<=dim;nd2++){ sprintf(filename1,"%s%d%s",filename1,indexra[nd2-1],finra[nd2]);}
    if (verbose){ printf(" %% checking for %s... ",filename1);}
    if (checktofind(filename1)){ if (verbose){ printf("found\n");} (*found_type)[index]=1;}
    else /*if (!checktofind(filename1))*/{ if (verbose){ printf("not found\n");} (*found_type)[index]=0;}}
  tfree(indexra);indexra=NULL;
  if (verbose){ raprintf(*found_type,"int",1,indexmax,"found_type: ");}
  if (verbose){ printf(" %% [exiting checktofind_howmany_ra]\n");}
}

void suite_8_power_process_helper2(char *gs2,int pulse_start,int pulse_length)
{
  /* simply collects the frequency response curves for the kazama-wilson experiment */
  int verbose=0;
  char dirname[1024],**finra=NULL,*tempchar=NULL,filename_base[1536];
  struct llist **nameLra=NULL;
  struct litem *l=NULL;
  int dim=0,nd=0,*dimra=NULL,*found_type=NULL,*found_rows=NULL,*found_cols=NULL,*dimra_min=NULL,*dimra_max=NULL;
  int total_types=0,total_vars=0,total_trials=0,total_concentrations=0,total_odors=0,total_bicucullines=0,total_records=0;
  int tab=0,tab2=0,tab3=0;
  int nconcentration=0,nodor=0,nbicuculline=0;
  int nt=0;
  double *ra=NULL,*ra_epsc=NULL,*ra_nAch=NULL,*ra_gabaA=NULL,*ra_gabaB=NULL,*ra_full=NULL,*ra_charge_128=NULL,*ra_charge_512=NULL;
  int output_type=0,output_rows=0,output_cols=0,output_nfound=0;
  int *ra_full_dim=NULL,*dim_to_plot=NULL;
  int total_bins=0;
  int total_glomeruli=GLOBAL_NCLUSTERS,ngli=0;
  double voltage_leak=0,conductance_leak=0;
  if (verbose){ printf(" %% [entering suite_8_power_process_helper2] gs2 %s\n",gs2);}
  sprintf(dirname,"./dir_temp_suite_8_%s",gs2);
  if (verbose){ printf(" %% checking for directory %s\n",dirname);}
  if (!checktofind(dirname)){ printf(" %% warning! cannot find %s in suite_8_power_process_helper\n",dirname);}
  else /* directory exists */{
    if (verbose){ printf(" %% found %s\n",dirname);}
    dim=6;
    finra = (char **) tcalloc(dim+1,sizeof(char *));
    for (nd=0;nd<=dim;nd++){ finra[nd]=(char *) tcalloc(512,sizeof(char));}
    sprintf(finra[0],"%s/power_deposit_type_",dirname);
    sprintf(finra[1],"_var_");
    sprintf(finra[2],"_instance_");
    sprintf(finra[3],"_concentration_");
    sprintf(finra[4],"_odor_");
    sprintf(finra[5],"_bicuculline_");
    sprintf(finra[6],"_%s",gs2);
    nameLra = (struct llist **) tcalloc(dim,sizeof(struct llist *)); for (nd=0;nd<dim;nd++){ nameLra[nd]=llistmake();}
    tempchar = (char *) tcalloc(64,sizeof(char));
    sprintf(tempchar,"%s",GLOBAL_TYPENAMES[TYPENAME_REGISTRY_wilson_PN]);
    litemadd(nameLra[0],tempchar);
    tempchar = (char *) tcalloc(64,sizeof(char));
    sprintf(tempchar,"%s",GLOBAL_TYPENAMES[TYPENAME_REGISTRY_wilson_LNe]);
    litemadd(nameLra[0],tempchar);
    tempchar = (char *) tcalloc(64,sizeof(char));
    sprintf(tempchar,"%s",GLOBAL_TYPENAMES[TYPENAME_REGISTRY_wilson_LNi]);
    litemadd(nameLra[0],tempchar);
    tempchar = (char *) tcalloc(64,sizeof(char));
    sprintf(tempchar,"%s",GLOBAL_VARNAMES[VARNAME_REGISTRY_wilson_nAch]);
    litemadd(nameLra[1],tempchar);
    tempchar = (char *) tcalloc(64,sizeof(char));
    sprintf(tempchar,"%s",GLOBAL_VARNAMES[VARNAME_REGISTRY_wilson_gabaA]);
    litemadd(nameLra[1],tempchar);
    tempchar = (char *) tcalloc(64,sizeof(char));
    sprintf(tempchar,"%s",GLOBAL_VARNAMES[VARNAME_REGISTRY_wilson_gabaB]);
    litemadd(nameLra[1],tempchar);
    checktofind_howmany_ra_general(dim,finra,nameLra,&dimra,&found_type,&found_rows,&found_cols);
    total_types = dimra[0];
    total_vars = dimra[1];
    total_trials = dimra[2];
    total_concentrations = dimra[3];
    total_odors = dimra[4];
    total_bicucullines = dimra[5];
    dimra_max = (int *) tcalloc(dim,sizeof(int));
    dimra_min = (int *) tcalloc(dim,sizeof(int));
    for (nbicuculline=0;nbicuculline<total_bicucullines;nbicuculline++){ for (nodor=0;nodor<total_odors;nodor++){ for (nconcentration=0;nconcentration<total_concentrations;nconcentration++){
      nt=0 /* TYPENAME_REGISTRY_wilson_PN */;
      switch (nt){ default: voltage_leak = -64.387; conductance_leak = 0.3; break;}
      dimra_max[0] = nt; dimra_min[0] = nt;
      dimra_max[1] = 0 /* VARNAME_REGISTRY_wilson_nAch */; dimra_min[1] = 0 /* VARNAME_REGISTRY_wilson_nAch */;
      dimra_max[2] = total_trials-1; dimra_min[2] = 0;
      dimra_max[3] = nconcentration; dimra_min[3] = nconcentration;
      dimra_max[4] = nodor; dimra_min[4] = nodor;
      dimra_max[5] = nbicuculline; dimra_min[5] = nbicuculline;
      checktofind_howmany_ra_general_average(dim,finra,nameLra,dimra,found_type,found_rows,found_cols,dimra_max,dimra_min,&ra_nAch,&output_type,&output_rows,&output_cols,&output_nfound);
      dimra_max[1] = 1 /* VARNAME_REGISTRY_wilson_gabaA */; dimra_min[1] = 1 /* VARNAME_REGISTRY_wilson_gabaA */;
      checktofind_howmany_ra_general_average(dim,finra,nameLra,dimra,found_type,found_rows,found_cols,dimra_max,dimra_min,&ra_gabaA,&output_type,&output_rows,&output_cols,&output_nfound);
      dimra_max[1] = 2 /* VARNAME_REGISTRY_wilson_gabaB */; dimra_min[1] = 2 /* VARNAME_REGISTRY_wilson_gabaB */;
      checktofind_howmany_ra_general_average(dim,finra,nameLra,dimra,found_type,found_rows,found_cols,dimra_max,dimra_min,&ra_gabaB,&output_type,&output_rows,&output_cols,&output_nfound);
      ra_epsc = (double *) tcalloc(output_rows*output_cols,sizeof(double));
      raaddequals(ra_epsc,output_rows*output_cols,-conductance_leak*(voltage_leak-voltage_leak));
      ratimesequals(ra_nAch,output_rows*output_cols,-(voltage_leak-VOLTAGE_[VARNAME_REGISTRY_wilson_nAch]));
      ratimesequals(ra_gabaA,output_rows*output_cols,-(voltage_leak-VOLTAGE_[VARNAME_REGISTRY_wilson_gabaA]));
      ratimesequals(ra_gabaB,output_rows*output_cols,-(voltage_leak-VOLTAGE_[VARNAME_REGISTRY_wilson_gabaB]));
      raplusequals(ra_epsc,output_rows*output_cols,ra_nAch);
      raplusequals(ra_epsc,output_rows*output_cols,ra_gabaA);
      raplusequals(ra_epsc,output_rows*output_cols,ra_gabaB);
      total_bins = output_rows;
      total_records = output_cols;
      tab2 = 0; tab3 = GLOBAL_NCLUSTERS; periodify("int",&nodor,&tab2,&tab3,&tab);
      ngli = tab; /* we assume that suite_8 emphasizes only the glomeruli corresponding to odor index */ 
      if (ra_full==NULL){ ra_full = (double *) tcalloc(pulse_length*(total_records/total_glomeruli)*total_concentrations*total_odors*total_bicucullines,sizeof(double));}
      ra = raplugout(ra_epsc,total_bins,total_records,pulse_start,ngli*(total_records/total_glomeruli),pulse_length,total_records/total_glomeruli);
      raplugin(&(ra_full[0 + nconcentration*pulse_length*(total_records/total_glomeruli) + nodor*total_concentrations*pulse_length*(total_records/total_glomeruli) + nbicuculline*total_concentrations*total_odors*pulse_length*(total_records/total_glomeruli)]),pulse_length,total_records/total_glomeruli,ra,pulse_length,total_records/total_glomeruli,0,0);
      tfree(ra);ra=NULL;
      if (ra_charge_128==NULL){ ra_charge_128 = (double *) tcalloc(total_concentrations*total_odors*total_bicucullines,sizeof(double));}
      ra = raplugout(ra_epsc,total_bins,total_records,pulse_start,ngli*(total_records/total_glomeruli),minimum(128,pulse_length),total_records/total_glomeruli);
      stats("double",ra,minimum(128,pulse_length)*total_records/total_glomeruli,NULL,NULL,&(ra_charge_128[nconcentration + nodor*total_concentrations + nbicuculline*total_concentrations*total_odors]),NULL);
      tfree(ra);ra=NULL;
      if (ra_charge_512==NULL){ ra_charge_512 = (double *) tcalloc(total_concentrations*total_odors*total_bicucullines,sizeof(double));}
      ra = raplugout(ra_epsc,total_bins,total_records,pulse_start,ngli*(total_records/total_glomeruli),minimum(512,pulse_length),total_records/total_glomeruli);
      stats("double",ra,minimum(512,pulse_length)*total_records/total_glomeruli,NULL,NULL,&(ra_charge_512[nconcentration + nodor*total_concentrations + nbicuculline*total_concentrations*total_odors]),NULL);
      tfree(ra);ra=NULL;
      tfree(ra_nAch);ra_nAch=NULL;
      tfree(ra_gabaA);ra_gabaA=NULL;
      tfree(ra_gabaB);ra_gabaB=NULL;
      tfree(ra_epsc);ra_epsc=NULL;}}}
    ra_full_dim = (int *) tcalloc(5,sizeof(int));
    ra_full_dim[0] = pulse_length;
    ra_full_dim[1] = total_records/total_glomeruli;
    ra_full_dim[2] = total_concentrations;
    ra_full_dim[3] = total_odors;
    ra_full_dim[4] = total_bicucullines;
    dim_to_plot = (int *) tcalloc(5,sizeof(int));
    dim_to_plot[0] = 0;
    dim_to_plot[1] = 1;
    dim_to_plot[2] = 1;
    dim_to_plot[3] = 1;
    dim_to_plot[4] = 1;
    sprintf(filename_base,"%s/suite_8_power_process_helper2_%s",dirname,gs2);
    if (ra_full!=NULL){ multidra2jpg(ra_full,5,ra_full_dim,filename_base,dim_to_plot,0,0);}
    tfree(dim_to_plot);dim_to_plot=NULL;
    tfree(ra_full_dim);ra_full_dim=NULL;
    ra_full_dim = (int *) tcalloc(3,sizeof(int));
    ra_full_dim[0] = total_concentrations;
    ra_full_dim[1] = total_odors;
    ra_full_dim[2] = total_bicucullines;
    dim_to_plot = (int *) tcalloc(3,sizeof(int));
    dim_to_plot[0] = 0;
    dim_to_plot[1] = 1;
    dim_to_plot[2] = 1;
    sprintf(filename_base,"%s/suite_8_power_process_helper2_charge_128_%s",dirname,gs2);
    if (ra_charge_128!=NULL){ multidra2jpg(ra_charge_128,3,ra_full_dim,filename_base,dim_to_plot,0,0);}
    tfree(dim_to_plot);dim_to_plot=NULL;
    tfree(ra_full_dim);ra_full_dim=NULL;
    ra_full_dim = (int *) tcalloc(3,sizeof(int));
    ra_full_dim[0] = total_concentrations;
    ra_full_dim[1] = total_odors;
    ra_full_dim[2] = total_bicucullines;
    dim_to_plot = (int *) tcalloc(3,sizeof(int));
    dim_to_plot[0] = 0;
    dim_to_plot[1] = 1;
    dim_to_plot[2] = 1;
    sprintf(filename_base,"%s/suite_8_power_process_helper2_charge_512_%s",dirname,gs2);
    if (ra_charge_512!=NULL){ multidra2jpg(ra_charge_512,3,ra_full_dim,filename_base,dim_to_plot,0,0);}
    tfree(dim_to_plot);dim_to_plot=NULL;
    tfree(ra_full_dim);ra_full_dim=NULL;
    if (ra_full!=NULL){ tfree(ra_full);ra_full=NULL;}
    if (ra_charge_128!=NULL){ tfree(ra_charge_128);ra_charge_128=NULL;}
    if (ra_charge_512!=NULL){ tfree(ra_charge_512);ra_charge_512=NULL;}
    tfree(dimra_max);dimra_max=NULL;
    tfree(dimra_min);dimra_min=NULL;
    tfree(dimra);dimra=NULL;
    tfree(found_type);found_type=NULL;
    tfree(found_rows);found_rows=NULL;
    tfree(found_cols);found_cols=NULL;
    for (nd=0;nd<dim;nd++){ l=nameLra[nd]->first; while (l!=NULL){ tfree(l->item);l->item=NULL;l=l->child;} llisttfree(nameLra[nd]); nameLra[nd]=NULL;} tfree(nameLra);nameLra=NULL;
    for (nd=0;nd<=dim;nd++){ tfree(finra[nd]);finra[nd]=NULL;} tfree(finra);finra=NULL;
  }
}

void suite_8_power_process_helper(char *gs2)
{
  int verbose=1;
  char dirname[1024],fpre[512],fpos[512],filename[2048],filename_base[1536],command[3072];
  int total_trials=0,total_concentrations=0,total_odors=0,total_bicucullines=0,total_length=0,total_records=0,temp_flag=0;
  int ntrial=0,nconcentration=0,nodor=0,nbicuculline=0,nlength=0,nrecord=0;
  int nt=0,nv=0,trials_found=0;
  double *ra=NULL,*ra_conductance=NULL,*ra_nAch=NULL,*ra_gabaA=NULL,*ra_gabaB=NULL;
  FILE *fp=NULL;
  double maxdia=10000,xpos=0,ypos=0;
  int remove_flag=SUITE_8_CLEANUP;
  int jpeg_flag=(ON_MY_COMPUTER ? 1 : 0),eps_flag=(ON_MY_COMPUTER ? 1 : 0)*(!remove_flag),pnm_flag=(!remove_flag);
  int conductance_plot_flag=1;
  int bin_size=0,bin_overlap=0,bin_skip=0,total_bins=0,nbin=0;
  int total_glomeruli=GLOBAL_NCLUSTERS,ngli=0;
  double temp_max=0,temp_min=0,temp_mean=0,temp_stdev=0;
  int nv_tab=0,nt_tab=0;
  double voltage_leak=0,conductance_leak=0,temp_epsc=0;
  if (verbose){ printf(" %% [entering suite_8_power_process_helper] gs2 %s\n",gs2);}
  sprintf(dirname,"./dir_temp_suite_8_%s",gs2);
  if (verbose){ printf(" %% checking for directory %s\n",dirname);}
  if (!checktofind(dirname)){ printf(" %% warning! cannot find %s in suite_8_power_process_helper\n",dirname);}
  else /* directory exists */{
    if (verbose){ printf(" %% found %s\n",dirname);}
    for (nt=0;nt<GLOBAL_NTYPES;nt++){ for (nv=0;nv<GLOBAL_NVARS;nv++){
      if (verbose){ printf(" %% checking for %s, %s\n",GLOBAL_TYPENAMES[nt],GLOBAL_VARNAMES[nv]);}
      sprintf(filename,"%s/power_deposit_type_%s_var_%s_instance_%d_concentration_%d_odor_%d_bicuculline_%d_%s",dirname,GLOBAL_TYPENAMES[nt],GLOBAL_VARNAMES[nv],0,0,0,0,gs2);
      if (checktofind(filename)){
	if (verbose>1){ printf(" %% Here we assume that the total_* variables can be obtained by looking along coordinate axes\n");}
	sprintf(fpre,"%s/power_deposit_type_%s_var_%s_instance_",dirname,GLOBAL_TYPENAMES[nt],GLOBAL_VARNAMES[nv]);
	sprintf(fpos,"_concentration_%d_odor_%d_bicuculline_%d_%s",0,0,0,gs2);
	total_trials=0; checktofind_howmany(fpre,&total_trials,fpos);
	sprintf(fpre,"%s/power_deposit_type_%s_var_%s_instance_%d_concentration_",dirname,GLOBAL_TYPENAMES[nt],GLOBAL_VARNAMES[nv],0);
	sprintf(fpos,"_odor_%d_bicuculline_%d_%s",0,0,gs2);
	total_concentrations=0; checktofind_howmany(fpre,&total_concentrations,fpos);
	sprintf(fpre,"%s/power_deposit_type_%s_var_%s_instance_%d_concentration_%d_odor_",dirname,GLOBAL_TYPENAMES[nt],GLOBAL_VARNAMES[nv],0,0);
	sprintf(fpos,"_bicuculline_%d_%s",0,gs2);
	total_odors=0; checktofind_howmany(fpre,&total_odors,fpos);
	sprintf(fpre,"%s/power_deposit_type_%s_var_%s_instance_%d_concentration_%d_odor_%d_bicuculline_",dirname,GLOBAL_TYPENAMES[nt],GLOBAL_VARNAMES[nv],0,0,0);
	sprintf(fpos,"_%s",gs2);
	total_bicucullines=0; checktofind_howmany(fpre,&total_bicucullines,fpos);
	ra = (double *) raread(filename,&temp_flag,&total_length,&total_records);tfree(ra);ra=NULL;
	if (verbose){ printf(" %% found %s,%s, total_trials %d total_concentrations %d total_odors %d total_bicucullines %d, total_length %d, total_records %d\n",GLOBAL_TYPENAMES[nt],GLOBAL_VARNAMES[nv],total_trials,total_concentrations,total_odors,total_bicucullines,total_length,total_records);}
	if (conductance_plot_flag && ((nt==TYPENAME_REGISTRY_wilson_ORN && (nv==VARNAME_REGISTRY_wilson_gabaA || nv==VARNAME_REGISTRY_wilson_gabaB || nv==VARNAME_REGISTRY_wilson_s_ORN || nv==VARNAME_REGISTRY_wilson_vesicle_depletion)) || (nt==TYPENAME_REGISTRY_wilson_PN && (nv==VARNAME_REGISTRY_wilson_nAch || nv==VARNAME_REGISTRY_wilson_gabaA || nv==VARNAME_REGISTRY_wilson_gabaB)) || (nt==TYPENAME_REGISTRY_wilson_LNe && (nv==VARNAME_REGISTRY_wilson_nAch || nv==VARNAME_REGISTRY_wilson_gabaA || nv==VARNAME_REGISTRY_wilson_gabaB)) || (nt==TYPENAME_REGISTRY_wilson_LNi && (nv==VARNAME_REGISTRY_wilson_nAch || nv==VARNAME_REGISTRY_wilson_gabaA || nv==VARNAME_REGISTRY_wilson_gabaB)))){
	  if (verbose){ printf(" %% plotting type %s, var %s\n",GLOBAL_TYPENAMES[nt],GLOBAL_VARNAMES[nv]);}
	  for (nbicuculline=0;nbicuculline<total_bicucullines;nbicuculline++){ for (nodor=0;nodor<total_odors;nodor++){ for (nconcentration=0;nconcentration<total_concentrations;nconcentration++){
	    if (verbose){ printf(" %% currently trying nbicuculline %d nodor %d nconcentration %d... ",nbicuculline,nodor,nconcentration);}
	    bin_size = maximum(1,GLOBAL_SPACE_SMOOTHER); bin_overlap = minimum(0,bin_size-1); bin_skip = bin_size-bin_overlap;
	    total_bins = total_length/bin_skip+1;
	    ra_conductance = (double *) tcalloc(total_bins*total_records,sizeof(double));
	    trials_found=0;
	    for (ntrial=0;ntrial<total_trials;ntrial++){
	      sprintf(filename,"%s/power_deposit_type_%s_var_%s_instance_%d_concentration_%d_odor_%d_bicuculline_%d_%s",dirname,GLOBAL_TYPENAMES[nt],GLOBAL_VARNAMES[nv],ntrial,nconcentration,nodor,nbicuculline,gs2);
	      if (checktofind(filename)){
		ra = (double *) raread(filename,&temp_flag,&nlength,&nrecord);
		assert(nlength==total_length && nrecord==total_records);
		for (nrecord=0;nrecord<total_records;nrecord++){
		  for (nbin=0;nbin<total_bins;nbin++){ ra_conductance[nbin+nrecord*total_bins]=0;}
		  for (nlength=0;nlength<total_length;nlength++){ 
		    for (nbin=nlength/bin_skip;nbin>=maximum(0,minimum(total_bins-1,(int)ceil((nlength+1-bin_size)/(double)bin_skip)));nbin--){ ra_conductance[nbin+nrecord*total_bins] += ra[nlength+nrecord*total_length];}}}
		tfree(ra);ra=NULL;
		trials_found += 1;}}
	    if (verbose){ printf(" %d total trials found\n",trials_found);}
	    for (nrecord=0;nrecord<total_records;nrecord++){ for (nbin=0;nbin<total_bins;nbin++){
	      ra_conductance[nbin+nrecord*total_bins] /= trials_found*(double)bin_size;}}
	    if (trials_found>0){
	      temp_max = GLOBAL_POWER_maxra_[nv]; temp_min = GLOBAL_POWER_minra_[nv];
	      stats("double",ra_conductance,total_bins*total_records,&temp_max,&temp_min,&temp_mean,&temp_stdev);
	      if (temp_max<=temp_min){ 
		stats("double",ra_conductance,total_bins*total_records,NULL,NULL,&temp_mean,&temp_stdev);
		temp_max = temp_mean + STD_VIEW*temp_stdev;
		temp_min = temp_mean - STD_VIEW*temp_stdev;}
	      sprintf(filename_base,"%s/power_deposit_type_%s_var_%s_concentration_%d_odor_%d_bicuculline_%d_%s",dirname,GLOBAL_TYPENAMES[nt],GLOBAL_VARNAMES[nv],nconcentration,nodor,nbicuculline,gs2);
	      if (verbose){ printf(" %% writing ra %s...",filename_base);} 
	      radump(ra_conductance,"double",total_bins,total_records,filename_base);
	      sprintf(filename,"%s.pnm",filename_base);
	      if (pnm_flag){ WritePNMfile_color(&(ra_conductance[0+0*total_bins]),total_bins,total_records,temp_max,temp_min,filename,7);}
	      sprintf(filename,"%s.fig",filename_base);
	      if ((fp=fopen(filename,"w"))==NULL){ printf(" %% warning! couldn't open %s in suite_8_power_process_helper\n",filename); fp=stdout;}
	      fprintf(fp,"%s",FIG_PREAMBLE);
	      for (ngli=0;ngli<total_glomeruli;ngli++){
		fprintf(fp,"2 1 0 1 %d 7 %d 0 -1 0.000 0 0 -1 0 0 %d\n\t",/*color*/11,/*depth*/999,/*npoints*/5);
		xpos = 0; ypos = (double)(ngli+0)/(double)total_glomeruli;
		fprintf(fp," %d %d",(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos));
		xpos = 1; ypos = (double)(ngli+0)/(double)total_glomeruli;
		fprintf(fp," %d %d",(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos));
		xpos = 1; ypos = (double)(ngli+1)/(double)total_glomeruli;
		fprintf(fp," %d %d",(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos));
		xpos = 0; ypos = (double)(ngli+1)/(double)total_glomeruli;
		fprintf(fp," %d %d",(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos));
		xpos = 0; ypos = (double)(ngli+0)/(double)total_glomeruli;
		fprintf(fp," %d %d",(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos));
		fprintf(fp,"\n");}
	      for (nrecord=0;nrecord<total_records;nrecord++){
		ngli = nrecord/(total_records/total_glomeruli);
		fprintf(fp,"2 1 0 1 %d 7 %d 0 -1 0.000 0 0 -1 0 0 %d\n\t",/*color*/(nrecord)%32,/*depth*/((nrecord)%999)+1,/*npoints*/total_bins);
		for (nbin=0;nbin<total_bins;nbin++){
		  xpos = ((double)nbin + 0.5)/(double)total_bins;
		  ypos = (ra_conductance[nbin+nrecord*total_bins]-temp_min)/(temp_max-temp_min);
		  ypos = maximum(0,minimum(1,ypos));
		  ypos *= 1.0/(double)total_glomeruli; ypos += (double)ngli/(double)total_glomeruli;
		  fprintf(fp," %d %d",(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos));}
		fprintf(fp,"\n");}
	      fprintf(fp,"4 0 0 %d 0 %d %d 0.0000 4 %d %d %d %d %0.16f\\001\n",/*depth*/1,/*font*/12,/*point*/24,/*textheight*/(int)(24*13.5),/*textwidth*/(int)(24*9*18),/*xpos*/(int)(0.1*maxdia),/*ypos*/(int)maxdia-(int)(1.0*maxdia),temp_max);
	      fprintf(fp,"4 0 0 %d 0 %d %d 0.0000 4 %d %d %d %d %0.16f\\001\n",/*depth*/1,/*font*/12,/*point*/24,/*textheight*/(int)(24*13.5),/*textwidth*/(int)(24*9*18),/*xpos*/(int)(0.1*maxdia),/*ypos*/(int)maxdia-(int)(0.0*maxdia),temp_min);
	      if (fp!=stdout){ fclose(fp);fp=NULL;
	      if (jpeg_flag){ sprintf(command,"fig2dev -Ljpeg -q %d %s.fig %s.jpg;",/*quality*/5,filename_base,filename_base); system(command);}
	      if (eps_flag){ sprintf(command,"fig2dev -Leps %s.fig %s.eps;",filename_base,filename_base); system(command);}
	      if (remove_flag){ sprintf(command,"rm %s.fig;",filename_base); system(command);}}}
	    tfree(ra_conductance);ra_conductance=NULL;}}}}
      }
    }}
    nt = TYPENAME_REGISTRY_wilson_ORN;
    for (nv_tab=0;nv_tab<2;nv_tab++){
      switch (nv_tab){
      case 0: nv = VARNAME_REGISTRY_wilson_s_ORN; break;
      case 1: nv = VARNAME_REGISTRY_wilson_vesicle_depletion; break;
      default: printf(" warning! nv_tab out of bounds in suite_8_power_process_helper\n"); break;}
      sprintf(filename,"%s/power_deposit_type_%s_var_%s_concentration_%d_odor_%d_bicuculline_%d_%s",dirname,GLOBAL_TYPENAMES[nt],GLOBAL_VARNAMES[nv],0,0,0,gs2);
      if (checktofind(filename)){
	if (verbose>1){ printf(" %% Here we assume that the total_* variables can be obtained by looking along coordinate axes\n");}
	sprintf(fpre,"%s/power_deposit_type_%s_var_%s_concentration_",dirname,GLOBAL_TYPENAMES[nt],GLOBAL_VARNAMES[nv]);
	sprintf(fpos,"_odor_%d_bicuculline_%d_%s",0,0,gs2);
	total_concentrations=0; checktofind_howmany(fpre,&total_concentrations,fpos);
	sprintf(fpre,"%s/power_deposit_type_%s_var_%s_concentration_%d_odor_",dirname,GLOBAL_TYPENAMES[nt],GLOBAL_VARNAMES[nv],0);
	sprintf(fpos,"_bicuculline_%d_%s",0,gs2);
	total_odors=0; checktofind_howmany(fpre,&total_odors,fpos);
	sprintf(fpre,"%s/power_deposit_type_%s_var_%s_concentration_%d_odor_%d_bicuculline_",dirname,GLOBAL_TYPENAMES[nt],GLOBAL_VARNAMES[nv],0,0);
	sprintf(fpos,"_%s",gs2);
	total_bicucullines=0; checktofind_howmany(fpre,&total_bicucullines,fpos);
	ra = (double *) raread(filename,&temp_flag,&total_bins,&total_records);tfree(ra);ra=NULL;
	if (verbose){ printf(" %% found %s,%s, total_concentrations %d total_odors %d total_bicucullines %d, total_bins %d, total_records %d\n",GLOBAL_TYPENAMES[nt],GLOBAL_VARNAMES[nv],total_concentrations,total_odors,total_bicucullines,total_bins,total_records);}
	for (nbicuculline=0;nbicuculline<total_bicucullines;nbicuculline++){ for (nodor=0;nodor<total_odors;nodor++){ for (nconcentration=0;nconcentration<total_concentrations;nconcentration++){
	  sprintf(filename,"%s/power_deposit_type_%s_var_%s_concentration_%d_odor_%d_bicuculline_%d_%s",dirname,GLOBAL_TYPENAMES[nt],GLOBAL_VARNAMES[nv],nconcentration,nodor,nbicuculline,gs2);
	  if (checktofind(filename)){
	    ra_conductance = (double *) tcalloc(total_bins*total_glomeruli,sizeof(double));
	    ra = (double *) raread(filename,&temp_flag,&nbin,&nrecord);
	    assert(nbin==total_bins && nrecord==total_records);
	    for (nbin=0;nbin<total_bins;nbin++){ 
	      for (nrecord=0;nrecord<total_records;nrecord++){
		ngli = nrecord/(total_records/total_glomeruli);
		ra_conductance[nbin+ngli*total_bins]+=ra[nbin+nrecord*total_bins];}
	      for (ngli=0;ngli<total_glomeruli;ngli++){
		ra_conductance[nbin+ngli*total_bins] /= total_records/total_glomeruli;}}
	    tfree(ra);ra=NULL;
	    temp_max = GLOBAL_POWER_maxra_[nv]; temp_min = GLOBAL_POWER_minra_[nv];
	    stats("double",ra_conductance,total_bins*total_glomeruli,&temp_max,&temp_min,&temp_mean,&temp_stdev);
	    if (temp_max<=temp_min){ 
	      stats("double",ra_conductance,total_bins*total_glomeruli,NULL,NULL,&temp_mean,&temp_stdev);
	      temp_max = temp_mean + STD_VIEW*temp_stdev;
	      temp_min = temp_mean - STD_VIEW*temp_stdev;}
	    sprintf(filename_base,"%s/power_gli_average_type_%s_var_%s_concentration_%d_odor_%d_bicuculline_%d_%s",dirname,GLOBAL_TYPENAMES[nt],GLOBAL_VARNAMES[nv],nconcentration,nodor,nbicuculline,gs2);
	    sprintf(filename,"%s.pnm",filename_base);
	    if (pnm_flag){ WritePNMfile_color(&(ra_conductance[0+0*total_bins]),total_bins,total_glomeruli,temp_max,temp_min,filename,7);}
	    sprintf(filename,"%s.fig",filename_base);
	    if ((fp=fopen(filename,"w"))==NULL){ printf(" %% warning! couldn't open %s in suite_8_power_process_helper\n",filename); fp=stdout;}
	    fprintf(fp,"%s",FIG_PREAMBLE);
	    for (ngli=0;ngli<total_glomeruli;ngli++){
	      fprintf(fp,"2 1 0 1 %d 7 %d 0 -1 0.000 0 0 -1 0 0 %d\n\t",/*color*/11,/*depth*/999,/*npoints*/5);
	      xpos = 0; ypos = (double)(ngli+0)/(double)total_glomeruli;
	      fprintf(fp," %d %d",(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos));
	      xpos = 1; ypos = (double)(ngli+0)/(double)total_glomeruli;
	      fprintf(fp," %d %d",(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos));
	      xpos = 1; ypos = (double)(ngli+1)/(double)total_glomeruli;
	      fprintf(fp," %d %d",(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos));
	      xpos = 0; ypos = (double)(ngli+1)/(double)total_glomeruli;
	      fprintf(fp," %d %d",(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos));
	      xpos = 0; ypos = (double)(ngli+0)/(double)total_glomeruli;
	      fprintf(fp," %d %d",(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos));
	      fprintf(fp,"\n");}
	    for (ngli=0;ngli<total_glomeruli;ngli++){
	      fprintf(fp,"2 1 0 1 %d 7 %d 0 -1 0.000 0 0 -1 0 0 %d\n\t",/*color*/(ngli)%32,/*depth*/((ngli)%999)+1,/*npoints*/total_bins);
	      for (nbin=0;nbin<total_bins;nbin++){
		xpos = ((double)nbin + 0.5)/(double)total_bins;
		ypos = (ra_conductance[nbin+ngli*total_bins]-temp_min)/(temp_max-temp_min);
		ypos = maximum(0,minimum(1,ypos));
		ypos *= 1.0/(double)total_glomeruli; ypos += (double)ngli/(double)total_glomeruli;
		fprintf(fp," %d %d",(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos));}
	      fprintf(fp,"\n");}
	    fprintf(fp,"4 0 0 %d 0 %d %d 0.0000 4 %d %d %d %d %0.16f\\001\n",/*depth*/1,/*font*/12,/*point*/24,/*textheight*/(int)(24*13.5),/*textwidth*/(int)(24*9*18),/*xpos*/(int)(0.1*maxdia),/*ypos*/(int)maxdia-(int)(1.0*maxdia),temp_max);
	    fprintf(fp,"4 0 0 %d 0 %d %d 0.0000 4 %d %d %d %d %0.16f\\001\n",/*depth*/1,/*font*/12,/*point*/24,/*textheight*/(int)(24*13.5),/*textwidth*/(int)(24*9*18),/*xpos*/(int)(0.1*maxdia),/*ypos*/(int)maxdia-(int)(0.0*maxdia),temp_min);
	    if (fp!=stdout){ fclose(fp);fp=NULL;}
	    if (jpeg_flag){ sprintf(command,"fig2dev -Ljpeg -q %d %s.fig %s.jpg;",/*quality*/5,filename_base,filename_base); system(command);}
	    if (eps_flag){ sprintf(command,"fig2dev -Leps %s.fig %s.eps;",filename_base,filename_base); system(command);}
	    if (remove_flag){ sprintf(command,"rm %s.fig;",filename_base); system(command);}
	    tfree(ra_conductance);ra_conductance=NULL;}}}}}}
    for (nt_tab=0;nt_tab<3;nt_tab++){
      switch (nt_tab){
      case 0: nt = TYPENAME_REGISTRY_wilson_PN; break;
      case 1: nt = TYPENAME_REGISTRY_wilson_LNe; break;
      case 2: nt = TYPENAME_REGISTRY_wilson_LNi; break;
      default: printf(" warning! nt_tab out of bounds in suite_8_power_process_helper\n"); break;}
      switch (nt){ default: voltage_leak = -64.387; conductance_leak = 0.3; break;}
      sprintf(filename,"%s/power_deposit_type_%s_var_%s_concentration_%d_odor_%d_bicuculline_%d_%s",dirname,GLOBAL_TYPENAMES[nt],GLOBAL_VARNAMES[VARNAME_REGISTRY_wilson_nAch],0,0,0,gs2); if (checktofind(filename)){ sprintf(filename,"%s/power_deposit_type_%s_var_%s_concentration_%d_odor_%d_bicuculline_%d_%s",dirname,GLOBAL_TYPENAMES[nt],GLOBAL_VARNAMES[VARNAME_REGISTRY_wilson_gabaA],0,0,0,gs2); if (checktofind(filename)){ sprintf(filename,"%s/power_deposit_type_%s_var_%s_concentration_%d_odor_%d_bicuculline_%d_%s",dirname,GLOBAL_TYPENAMES[nt],GLOBAL_VARNAMES[VARNAME_REGISTRY_wilson_gabaB],0,0,0,gs2); if (checktofind(filename)){
	if (verbose>1){ printf(" %% Here we assume that var_%s is representative of all total_* variables, and that the total_* variables can be obtained by looking along coordinate axes\n",GLOBAL_VARNAMES[VARNAME_REGISTRY_wilson_nAch]);}
	sprintf(fpre,"%s/power_deposit_type_%s_var_%s_concentration_",dirname,GLOBAL_TYPENAMES[nt],GLOBAL_VARNAMES[VARNAME_REGISTRY_wilson_nAch]); sprintf(fpos,"_odor_%d_bicuculline_%d_%s",0,0,gs2); total_concentrations=0; checktofind_howmany(fpre,&total_concentrations,fpos);
	sprintf(fpre,"%s/power_deposit_type_%s_var_%s_concentration_%d_odor_",dirname,GLOBAL_TYPENAMES[nt],GLOBAL_VARNAMES[VARNAME_REGISTRY_wilson_nAch],0); sprintf(fpos,"_bicuculline_%d_%s",0,gs2); total_odors=0; checktofind_howmany(fpre,&total_odors,fpos);
	sprintf(fpre,"%s/power_deposit_type_%s_var_%s_concentration_%d_odor_%d_bicuculline_",dirname,GLOBAL_TYPENAMES[nt],GLOBAL_VARNAMES[VARNAME_REGISTRY_wilson_nAch],0,0); sprintf(fpos,"_%s",gs2); total_bicucullines=0; checktofind_howmany(fpre,&total_bicucullines,fpos);
	ra = (double *) raread(filename,&temp_flag,&total_bins,&total_records);tfree(ra);ra=NULL;
	if (verbose){ printf(" %% found %s,%s, total_concentrations %d total_odors %d total_bicucullines %d, total_bins %d, total_records %d\n",GLOBAL_TYPENAMES[nt],GLOBAL_VARNAMES[VARNAME_REGISTRY_wilson_nAch],total_concentrations,total_odors,total_bicucullines,total_bins,total_records);}
	for (nbicuculline=0;nbicuculline<total_bicucullines;nbicuculline++){ for (nodor=0;nodor<total_odors;nodor++){ for (nconcentration=0;nconcentration<total_concentrations;nconcentration++){
	  sprintf(filename,"%s/power_deposit_type_%s_var_%s_concentration_%d_odor_%d_bicuculline_%d_%s",dirname,GLOBAL_TYPENAMES[nt],GLOBAL_VARNAMES[VARNAME_REGISTRY_wilson_nAch],nconcentration,nodor,nbicuculline,gs2);
	  ra_nAch = (double *) raread(filename,&temp_flag,&nbin,&nrecord); assert(nbin==total_bins && nrecord==total_records);
	  sprintf(filename,"%s/power_deposit_type_%s_var_%s_concentration_%d_odor_%d_bicuculline_%d_%s",dirname,GLOBAL_TYPENAMES[nt],GLOBAL_VARNAMES[VARNAME_REGISTRY_wilson_gabaA],nconcentration,nodor,nbicuculline,gs2);
	  ra_gabaA = (double *) raread(filename,&temp_flag,&nbin,&nrecord); assert(nbin==total_bins && nrecord==total_records);
	  sprintf(filename,"%s/power_deposit_type_%s_var_%s_concentration_%d_odor_%d_bicuculline_%d_%s",dirname,GLOBAL_TYPENAMES[nt],GLOBAL_VARNAMES[VARNAME_REGISTRY_wilson_gabaB],nconcentration,nodor,nbicuculline,gs2);
	  ra_gabaB = (double *) raread(filename,&temp_flag,&nbin,&nrecord); assert(nbin==total_bins && nrecord==total_records);
	  ra_conductance = (double *) tcalloc(total_bins*total_glomeruli,sizeof(double));
	  for (nbin=0;nbin<total_bins;nbin++){ 
	    for (nrecord=0;nrecord<total_records;nrecord++){
	      ngli = nrecord/(total_records/total_glomeruli);
	      temp_epsc = 0 - conductance_leak*(voltage_leak - voltage_leak) - ra_nAch[nbin+nrecord*total_bins]*(voltage_leak-VOLTAGE_[VARNAME_REGISTRY_wilson_nAch]) - ra_gabaA[nbin+nrecord*total_bins]*(voltage_leak-VOLTAGE_[VARNAME_REGISTRY_wilson_gabaA]) - ra_gabaB[nbin+nrecord*total_bins]*(voltage_leak-VOLTAGE_[VARNAME_REGISTRY_wilson_gabaB]);
	      ra_conductance[nbin+ngli*total_bins]+=temp_epsc;}
	    for (ngli=0;ngli<total_glomeruli;ngli++){
	      ra_conductance[nbin+ngli*total_bins] /= total_records/total_glomeruli;}}
	  tfree(ra_nAch);ra_nAch=NULL;
	  tfree(ra_gabaA);ra_gabaA=NULL;
	  tfree(ra_gabaB);ra_gabaB=NULL;
	  stats("double",ra_conductance,total_bins*total_glomeruli,NULL,NULL,&temp_mean,&temp_stdev);
	  temp_max = temp_mean + STD_VIEW*temp_stdev; temp_min = temp_mean - STD_VIEW*temp_stdev;
	  stats("double",ra_conductance,total_bins*total_glomeruli,&temp_max,&temp_min,&temp_mean,&temp_stdev);
	  sprintf(filename_base,"%s/power_gli_average_type_%s_var_epsc_concentration_%d_odor_%d_bicuculline_%d_%s",dirname,GLOBAL_TYPENAMES[nt],nconcentration,nodor,nbicuculline,gs2);
	  sprintf(filename,"%s.pnm",filename_base);
	  if (pnm_flag){ WritePNMfile_color(&(ra_conductance[0+0*total_bins]),total_bins,total_glomeruli,temp_max,temp_min,filename,7);}
	  sprintf(filename,"%s.fig",filename_base);
	  if ((fp=fopen(filename,"w"))==NULL){ printf(" %% warning! couldn't open %s in suite_8_power_process_helper\n",filename); fp=stdout;}
	  fprintf(fp,"%s",FIG_PREAMBLE);
	  for (ngli=0;ngli<total_glomeruli;ngli++){
	    fprintf(fp,"2 1 0 1 %d 7 %d 0 -1 0.000 0 0 -1 0 0 %d\n\t",/*color*/11,/*depth*/999,/*npoints*/5);
	    xpos = 0; ypos = (double)(ngli+0)/(double)total_glomeruli;
	    fprintf(fp," %d %d",(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos));
	    xpos = 1; ypos = (double)(ngli+0)/(double)total_glomeruli;
	    fprintf(fp," %d %d",(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos));
	    xpos = 1; ypos = (double)(ngli+1)/(double)total_glomeruli;
	    fprintf(fp," %d %d",(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos));
	    xpos = 0; ypos = (double)(ngli+1)/(double)total_glomeruli;
	    fprintf(fp," %d %d",(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos));
	    xpos = 0; ypos = (double)(ngli+0)/(double)total_glomeruli;
	    fprintf(fp," %d %d",(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos));
	    fprintf(fp,"\n");}
	  for (ngli=0;ngli<total_glomeruli;ngli++){
	    fprintf(fp,"2 1 0 1 %d 7 %d 0 -1 0.000 0 0 -1 0 0 %d\n\t",/*color*/(ngli)%32,/*depth*/((ngli)%999)+1,/*npoints*/total_bins);
	    for (nbin=0;nbin<total_bins;nbin++){
	      xpos = ((double)nbin + 0.5)/(double)total_bins;
	      ypos = (ra_conductance[nbin+ngli*total_bins]-temp_min)/(temp_max-temp_min);
	      ypos = maximum(0,minimum(1,ypos));
	      ypos *= 1.0/(double)total_glomeruli; ypos += (double)ngli/(double)total_glomeruli;
	      fprintf(fp," %d %d",(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos));}
	    fprintf(fp,"\n");}
	  fprintf(fp,"4 0 0 %d 0 %d %d 0.0000 4 %d %d %d %d %0.16f\\001\n",/*depth*/1,/*font*/12,/*point*/24,/*textheight*/(int)(24*13.5),/*textwidth*/(int)(24*9*18),/*xpos*/(int)(0.1*maxdia),/*ypos*/(int)maxdia-(int)(1.0*maxdia),temp_max);
	  fprintf(fp,"4 0 0 %d 0 %d %d 0.0000 4 %d %d %d %d %0.16f\\001\n",/*depth*/1,/*font*/12,/*point*/24,/*textheight*/(int)(24*13.5),/*textwidth*/(int)(24*9*18),/*xpos*/(int)(0.1*maxdia),/*ypos*/(int)maxdia-(int)(0.0*maxdia),temp_min);
	  if (fp!=stdout){ fclose(fp);fp=NULL;}
	  if (jpeg_flag){ sprintf(command,"fig2dev -Ljpeg -q %d %s.fig %s.jpg;",/*quality*/5,filename_base,filename_base); system(command);}
	  if (eps_flag){ sprintf(command,"fig2dev -Leps %s.fig %s.eps;",filename_base,filename_base); system(command);}
	  if (remove_flag){ sprintf(command,"rm %s.fig;",filename_base); system(command);}
	  tfree(ra_conductance);ra_conductance=NULL;}}}}}}}
    if (remove_flag){ 
      if (verbose){ printf(" %% cleaning up trial averages\n");}
      for (nt=0;nt<GLOBAL_NTYPES;nt++){ for (nv=0;nv<GLOBAL_NVARS;nv++){
	if (verbose){ printf(" %% checking for %s, %s\n",GLOBAL_TYPENAMES[nt],GLOBAL_VARNAMES[nv]);}
	sprintf(filename,"%s/power_deposit_type_%s_var_%s_instance_%d_concentration_%d_odor_%d_bicuculline_%d_%s",dirname,GLOBAL_TYPENAMES[nt],GLOBAL_VARNAMES[nv],0,0,0,0,gs2);
	if (checktofind(filename)){
	  if (verbose>1){ printf(" %% Here we assume that the total_* variables can be obtained by looking along coordinate axes\n");}
	  sprintf(fpre,"%s/power_deposit_type_%s_var_%s_instance_",dirname,GLOBAL_TYPENAMES[nt],GLOBAL_VARNAMES[nv]);
	  sprintf(fpos,"_concentration_%d_odor_%d_bicuculline_%d_%s",0,0,0,gs2);
	  total_trials=0; checktofind_howmany(fpre,&total_trials,fpos);
	  sprintf(fpre,"%s/power_deposit_type_%s_var_%s_instance_%d_concentration_",dirname,GLOBAL_TYPENAMES[nt],GLOBAL_VARNAMES[nv],0);
	  sprintf(fpos,"_odor_%d_bicuculline_%d_%s",0,0,gs2);
	  total_concentrations=0; checktofind_howmany(fpre,&total_concentrations,fpos);
	  sprintf(fpre,"%s/power_deposit_type_%s_var_%s_instance_%d_concentration_%d_odor_",dirname,GLOBAL_TYPENAMES[nt],GLOBAL_VARNAMES[nv],0,0);
	  sprintf(fpos,"_bicuculline_%d_%s",0,gs2);
	  total_odors=0; checktofind_howmany(fpre,&total_odors,fpos);
	  sprintf(fpre,"%s/power_deposit_type_%s_var_%s_instance_%d_concentration_%d_odor_%d_bicuculline_",dirname,GLOBAL_TYPENAMES[nt],GLOBAL_VARNAMES[nv],0,0,0);
	  sprintf(fpos,"_%s",gs2);
	  total_bicucullines=0; checktofind_howmany(fpre,&total_bicucullines,fpos);
	  ra = (double *) raread(filename,&temp_flag,&total_length,&total_records);tfree(ra);ra=NULL;
	  for (nbicuculline=0;nbicuculline<total_bicucullines;nbicuculline++){ for (nodor=0;nodor<total_odors;nodor++){ for (nconcentration=0;nconcentration<total_concentrations;nconcentration++){
	    sprintf(filename_base,"%s/power_deposit_type_%s_var_%s_concentration_%d_odor_%d_bicuculline_%d_%s",dirname,GLOBAL_TYPENAMES[nt],GLOBAL_VARNAMES[nv],nconcentration,nodor,nbicuculline,gs2);
	    sprintf(command,"rm %s;",filename_base); system(command);}}}}}}}
  }
}

void suite_8_power_deposit_helper(int instance,int concentration,int odor,int bicuculline,int time)
{
  int verbose=1;
  char dirname[1024],filename[1536],gs2[128],command[2048];
  struct power *p=GLOBAL_POWER;
  int nt=0,nrecords=0,nv=0,sttab=0,stlength=0,nr=0,nl=0,nl2=0;
  double *ra=NULL;
  if (verbose){ printf(" %% [entering suite_8_power_deposit_helper] trial %d concentration %d odor %d bicuculline %d time %d\n",instance,concentration,odor,bicuculline,time);}
  sprintf(gs2,"%srecord",GLOBAL_STRING_2);
  sprintf(dirname,"./dir_temp_suite_8_%s",gs2);
  if (verbose){ printf(" %% checking for directory %s\n",dirname);}
  if (!checktofind(dirname)){ if (verbose){ printf(" %% making %s\n",dirname);} sprintf(command,"mkdir %s;",dirname); system(command);}
  for (nt=0;nt<p->indexing_ntype_length;nt++){ 
    nrecords = p->Nra->lengthra[p->indexing_ntype_checkout[nt]];
    for (nv=0;nv<p->indexing_nvar_length;nv++){
      if (p->strarara[nt][nv]!=NULL){
	if (verbose){ printf(" %% depositing p->strarara[%s][%s]\n",GLOBAL_TYPENAMES[p->indexing_ntype_checkout[nt]],GLOBAL_VARNAMES[p->indexing_nvar_checkout[nv]]);}
	sttab = p->strarara[nt][nv][0]->tab; stlength = p->strarara[nt][nv][0]->length;
	ra = (double *) tcalloc(p->length*nrecords,sizeof(double));
	for (nr=0;nr<nrecords;nr++){ for (nl=0;nl<p->length;nl++){
	  nl2 = periodize(sttab-p->length + nl,0,stlength);
	  ra[nl + nr*p->length] = p->strarara[nt][nv][nr]->data[nl2];}}
	sprintf(filename,"%s/power_deposit_type_%s_var_%s_instance_%d_concentration_%d_odor_%d_bicuculline_%d_%s",dirname,GLOBAL_TYPENAMES[p->indexing_ntype_checkout[nt]],GLOBAL_VARNAMES[p->indexing_nvar_checkout[nv]],instance,concentration,odor,bicuculline,gs2);
	if (verbose){ printf(" %% writing ra %s...",filename);}
	radump(ra,"double",p->length,nrecords,filename);
	if (verbose){ printf(" freeing ra\n");}
	tfree(ra);ra=NULL;}}}
}

int spike_flag_to_spike_time(double t1,double t2,double dt,double tau,double increment,double *spiketime)
{
  /* backs out spiketime from spike_flag, intended for wilson_spike_flag, which assumes that spike_flag decays with timescale TAU_[VARNAME_REGISTRY_wilson_spike_flag], and is incremented by 1.0 for every spike */
  int verbose=0;
  int spike_found=0;
  double tolerance=0.000001;
  if (t2<=t1*exp(-dt/tau)+tolerance){ 
    if (verbose){ printf(" %% t2 %f > t1*exp(-dt/tau) %f, no spike\n",t2,t1*exp(-dt/tau));}
    spike_found=0;}
  else /* if (t2>t1*exp(-dt/tau)+tolerance) */{ 
    *spiketime = dt + tau*log((t2-t1*exp(-dt/tau))/increment);
    if (verbose){ printf(" %% t2-t1*exp(-dt/tau) %f, spike found at spiketime %f\n",t2-t1*exp(-dt/tau),*spiketime);}
    if (*spiketime<0 || *spiketime>dt){ printf(" %% Warning! t2 %f, t1 %f, dt %f, tau %f, increment %f, tau*log((t2-t1*exp(-dt/tau))/increment) %f, 0<spiketime %f<dt %f in spike_flag_to_spike_time\n",t2,t1,dt,tau,increment,tau*log((t2-t1*exp(-dt/tau))/increment),*spiketime,dt);}
    spike_found=1;}
  return spike_found;
}

void suite_7_power_process_helper_2(int total_clusters,int total_odors,double *ra1,int total_1_records,double *ra2,int total_2_records,char *fgvn)
{
  /* assumes *fgvn begins with "./" */
  /* the *ra? should be of size total_odors*total_?_records*2 (mean,stdev) */
  /* sets up a scatter plot of mean +/- stdev for ra1 (x-axis) vs ra2 (y-axis) */
  /* attempting to ensure that the points are visible */
  int nr0,nr1=0,nr2=0,nr1_a=0,nr2_a=0,no=0,nl=0;
  char filename_base[1024],filename[2048],command[4096];
  FILE *fp=NULL;
  double maxdia=10000,xpos=0,ypos=0,rcolor=0,gcolor=0,bcolor=0;
  int colorcode=0;
  int remove_flag=SUITE_7_CLEANUP,jpeg_flag=(ON_MY_COMPUTER ? 1 : 0),eps_flag=(ON_MY_COMPUTER ? 1 : 0);
  double ra1max=0,ra2max=0;
  stats("double",ra1,total_odors*total_1_records,&ra1max,NULL,NULL,NULL);
  stats("double",ra2,total_odors*total_2_records,&ra2max,NULL,NULL,NULL);
  sprintf(filename_base,"%s",fgvn);
  sprintf(filename,"%s.fig",filename_base);
  if ((fp=fopen(filename,"w"))==NULL){ printf(" %% warning! not open in %s suite_7_power_process_helper_1\n",filename); fp=stdout;}
  fprintf(fp,"%s",FIG_PREAMBLE); fprintf(fp,"%s",FIG_PREAMBLE_COLOR_7);
  fprintf(fp,"2 1 0 1 %d 7 %d 0 -1 0.000 0 0 -1 0 0 %d\n\t",/*color*/0,/*depth*/999,/*npoints*/2);
  xpos=0;ypos=0; fprintf(fp," %d %d",(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos));
  xpos=0;ypos=1; fprintf(fp," %d %d",(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos));
  fprintf(fp,"\n");
  fprintf(fp,"2 1 0 1 %d 7 %d 0 -1 0.000 0 0 -1 0 0 %d\n\t",/*color*/0,/*depth*/999,/*npoints*/2);
  xpos=0;ypos=0; fprintf(fp," %d %d",(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos));
  xpos=1;ypos=0; fprintf(fp," %d %d",(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos));
  fprintf(fp,"\n");
  for (nr0=0;nr0<total_clusters;nr0++){
    colorscale(0,(double)(nr0+0.5)/(double)total_clusters,1,0,&rcolor,&gcolor,&bcolor);
    colorcode = crop((int)floor(512*rcolor),0,511);
    colorcode = 0;
    for (nr1=0;nr1<total_1_records/total_clusters;nr1++){ for (nr2=0;nr2<total_2_records/total_clusters;nr2++){
      nr1_a = nr1 + nr0*(total_1_records/total_clusters); nr2_a = nr2 + nr0*(total_2_records/total_clusters);
      for (no=0;no<total_odors;no++){
	/* fprintf(fp,"2 1 0 1 %d 7 %d 0 -1 0.000 0 0 -1 0 0 %d\n\t",/\*color*\/colorcode+32,/\*depth*\/998,/\*npoints*\/2);	 */
	/* xpos=ra1[no+nr1_a*total_odors+0*total_odors*total_1_records]-ra1[no+nr1_a*total_odors+1*total_odors*total_1_records];xpos/=ra1max; */
	/* ypos=ra2[no+nr2_a*total_odors+0*total_odors*total_2_records];ypos/=ra2max; */
	/* fprintf(fp," %d %d",(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos)); */
	/* xpos=ra1[no+nr1_a*total_odors+0*total_odors*total_1_records]+ra1[no+nr1_a*total_odors+1*total_odors*total_1_records];xpos/=ra1max; */
	/* ypos=ra2[no+nr2_a*total_odors+0*total_odors*total_2_records];ypos/=ra2max; */
	/* fprintf(fp," %d %d",(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos)); */
	/* fprintf(fp,"\n"); */
	/* fprintf(fp,"2 1 0 1 %d 7 %d 0 -1 0.000 0 0 -1 0 0 %d\n\t",/\*color*\/colorcode+32,/\*depth*\/998,/\*npoints*\/2);	 */
	/* xpos=ra1[no+nr1_a*total_odors+0*total_odors*total_1_records];xpos/=ra1max; */
	/* ypos=ra2[no+nr2_a*total_odors+0*total_odors*total_2_records]-ra2[no+nr2_a*total_odors+1*total_odors*total_2_records];ypos/=ra2max; */
	/* fprintf(fp," %d %d",(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos)); */
	/* xpos=ra1[no+nr1_a*total_odors+0*total_odors*total_1_records];xpos/=ra1max; */
	/* ypos=ra2[no+nr2_a*total_odors+0*total_odors*total_2_records]+ra2[no+nr2_a*total_odors+1*total_odors*total_2_records];ypos/=ra2max; */
	/* fprintf(fp," %d %d",(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos)); */
	/* fprintf(fp,"\n"); */
	for (nl=0;nl<1;nl++){
	  xpos=ra1[no+nr1_a*total_odors+0*total_odors*total_1_records] + 0.5*randn()*ra1[no+nr1_a*total_odors+1*total_odors*total_1_records];xpos/=ra1max;
	  ypos=ra2[no+nr2_a*total_odors+0*total_odors*total_2_records] + 0.5*randn()*ra2[no+nr2_a*total_odors+1*total_odors*total_2_records];ypos/=ra2max;
	  fprintf(fp,"1 3 0 %d %d %d %d 0 20 0.000 1 0.0000 %d %d %d %d %d %d %d %d\n",/*width*/1,/*outercolor*/0,/*innercolor*/colorcode+32,/*depth*/maximum(1,997-nr0),(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos),/*radius*/100,/*radius*/100,(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos),(int)floor(maxdia*xpos)+/*radius*/100,(int)maxdia-(int)floor(maxdia*ypos));
	/* for (nl=0;nl<24;nl++){ } */}
      }}}}
  fprintf(fp,"4 0 0 %d 0 %d %d 0.0000 4 %d %d %d %d %0.16f\\001\n",/*depth*/1,/*font*/12,/*point*/24,/*textheight*/(int)(24*13.5),/*textwidth*/(int)(24*9*18),/*xpos*/(int)(1.0*maxdia),/*ypos*/(int)maxdia-(int)(0.1*maxdia),ra1max);
  fprintf(fp,"4 0 0 %d 0 %d %d 0.0000 4 %d %d %d %d %0.16f\\001\n",/*depth*/1,/*font*/12,/*point*/24,/*textheight*/(int)(24*13.5),/*textwidth*/(int)(24*9*18),/*xpos*/(int)(0.1*maxdia),/*ypos*/(int)maxdia-(int)(1.0*maxdia),ra2max);
  if (fp!=stdout){ fclose(fp);fp=NULL;}
  if (jpeg_flag){ sprintf(command,"fig2dev -Ljpeg -q %d %s.fig %s.jpg;",/*quality*/5,filename_base,filename_base); system(command);}
  if (eps_flag){ sprintf(command,"fig2dev -Leps %s.fig %s.eps;",filename_base,filename_base); system(command);}
  if (remove_flag){ sprintf(command,"rm %s.fig;",filename_base); system(command);}
}

void suite_7_power_process_helper_1(int total_clusters,int total_odors,double *ra1,int total_1_records,double *ra2,int total_2_records,char *fgvn)
{
  /* assumes *fgvn begins with "./" */
  /* the *ra? should be of size total_odors*total_?_records*2 (mean,stdev) */
  /* sets up a scatter plot of mean +/- stdev for ra1 (x-axis) vs ra2 (y-axis) */
  int nr0,nr1=0,nr2=0,nr1_a=0,nr2_a=0,no=0;
  char filename_base[1024],filename[2048],command[4096];
  FILE *fp=NULL;
  double maxdia=10000,xpos=0,ypos=0,rcolor=0,gcolor=0,bcolor=0;
  int colorcode=0;
  int remove_flag=SUITE_7_CLEANUP,jpeg_flag=(ON_MY_COMPUTER ? 1 : 0),eps_flag=(ON_MY_COMPUTER ? 1 : 0);
  double ra1max=0,ra2max=0;
  stats("double",ra1,total_odors*total_1_records,&ra1max,NULL,NULL,NULL);
  stats("double",ra2,total_odors*total_2_records,&ra2max,NULL,NULL,NULL);
  sprintf(filename_base,"%s",fgvn);
  sprintf(filename,"%s.fig",filename_base);
  if ((fp=fopen(filename,"w"))==NULL){ printf(" %% warning! not open in %s suite_7_power_process_helper_1\n",filename); fp=stdout;}
  fprintf(fp,"%s",FIG_PREAMBLE); fprintf(fp,"%s",FIG_PREAMBLE_COLOR_7);
  fprintf(fp,"2 1 0 1 %d 7 %d 0 -1 0.000 0 0 -1 0 0 %d\n\t",/*color*/0,/*depth*/999,/*npoints*/2);
  xpos=0;ypos=0; fprintf(fp," %d %d",(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos));
  xpos=0;ypos=1; fprintf(fp," %d %d",(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos));
  fprintf(fp,"\n");
  fprintf(fp,"2 1 0 1 %d 7 %d 0 -1 0.000 0 0 -1 0 0 %d\n\t",/*color*/0,/*depth*/999,/*npoints*/2);
  xpos=0;ypos=0; fprintf(fp," %d %d",(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos));
  xpos=1;ypos=0; fprintf(fp," %d %d",(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos));
  fprintf(fp,"\n");
  for (nr0=0;nr0<total_clusters;nr0++){
    colorscale(0,(double)(nr0+0.5)/(double)total_clusters,1,0,&rcolor,&gcolor,&bcolor);
    colorcode = crop((int)floor(512*rcolor),0,511);
    for (nr1=0;nr1<total_1_records/total_clusters;nr1++){ for (nr2=0;nr2<total_2_records/total_clusters;nr2++){
      nr1_a = nr1 + nr0*(total_1_records/total_clusters); nr2_a = nr2 + nr0*(total_2_records/total_clusters);
      for (no=0;no<total_odors;no++){
	fprintf(fp,"2 1 0 1 %d 7 %d 0 -1 0.000 0 0 -1 0 0 %d\n\t",/*color*/colorcode+32,/*depth*/998,/*npoints*/2);	
	xpos=ra1[no+nr1_a*total_odors+0*total_odors*total_1_records]-ra1[no+nr1_a*total_odors+1*total_odors*total_1_records];xpos/=ra1max;
	ypos=ra2[no+nr2_a*total_odors+0*total_odors*total_2_records];ypos/=ra2max;
	fprintf(fp," %d %d",(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos));
	xpos=ra1[no+nr1_a*total_odors+0*total_odors*total_1_records]+ra1[no+nr1_a*total_odors+1*total_odors*total_1_records];xpos/=ra1max;
	ypos=ra2[no+nr2_a*total_odors+0*total_odors*total_2_records];ypos/=ra2max;
	fprintf(fp," %d %d",(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos));
	fprintf(fp,"\n");
	fprintf(fp,"2 1 0 1 %d 7 %d 0 -1 0.000 0 0 -1 0 0 %d\n\t",/*color*/colorcode+32,/*depth*/998,/*npoints*/2);	
	xpos=ra1[no+nr1_a*total_odors+0*total_odors*total_1_records];xpos/=ra1max;
	ypos=ra2[no+nr2_a*total_odors+0*total_odors*total_2_records]-ra2[no+nr2_a*total_odors+1*total_odors*total_2_records];ypos/=ra2max;
	fprintf(fp," %d %d",(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos));
	xpos=ra1[no+nr1_a*total_odors+0*total_odors*total_1_records];xpos/=ra1max;
	ypos=ra2[no+nr2_a*total_odors+0*total_odors*total_2_records]+ra2[no+nr2_a*total_odors+1*total_odors*total_2_records];ypos/=ra2max;
	fprintf(fp," %d %d",(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos));
	fprintf(fp,"\n");
	xpos=ra1[no+nr1_a*total_odors+0*total_odors*total_1_records];xpos/=ra1max;
	ypos=ra2[no+nr2_a*total_odors+0*total_odors*total_2_records];ypos/=ra2max;
	fprintf(fp,"1 3 0 %d %d %d %d 0 20 0.000 1 0.0000 %d %d %d %d %d %d %d %d\n",/*width*/2,/*outercolor*/0,/*innercolor*/colorcode+32,/*depth*/maximum(1,997-nr0),(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos),/*radius*/100,/*radius*/100,(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos),(int)floor(maxdia*xpos)+/*radius*/100,(int)maxdia-(int)floor(maxdia*ypos));}}}}
  fprintf(fp,"4 0 0 %d 0 %d %d 0.0000 4 %d %d %d %d %0.16f\\001\n",/*depth*/1,/*font*/12,/*point*/24,/*textheight*/(int)(24*13.5),/*textwidth*/(int)(24*9*18),/*xpos*/(int)(1.0*maxdia),/*ypos*/(int)maxdia-(int)(0.1*maxdia),ra1max);
  fprintf(fp,"4 0 0 %d 0 %d %d 0.0000 4 %d %d %d %d %0.16f\\001\n",/*depth*/1,/*font*/12,/*point*/24,/*textheight*/(int)(24*13.5),/*textwidth*/(int)(24*9*18),/*xpos*/(int)(0.1*maxdia),/*ypos*/(int)maxdia-(int)(1.0*maxdia),ra2max);
  if (fp!=stdout){ fclose(fp);fp=NULL;}
  if (jpeg_flag){ sprintf(command,"fig2dev -Ljpeg -q %d %s.fig %s.jpg;",/*quality*/5,filename_base,filename_base); system(command);}
  if (eps_flag){ sprintf(command,"fig2dev -Leps %s.fig %s.eps;",filename_base,filename_base); system(command);}
  if (remove_flag){ sprintf(command,"rm %s.fig;",filename_base); system(command);}
}

void spearman_rank_correlation_histogram(int total_clusters,int total_odors,double *ra1,int total_1_records,double *ra2,int total_2_records,char *fgvn)
{
  /* assumes *fgvn begins with "./" */
  /* the *ra? should be of size total_odors*total_?_records*2 (mean,stdev) */
  /* assumes that total_?_records are grouped into total_clusters clumps of adjacent records */
  int verbose=0;
  int total_trials = 2048*total_clusters;
  int nt=0,nr0,nr1=0,nr2=0,no=0;
  int nbin=0,total_bins=64;
  struct hist *h=NULL;
  double *ra_temp=NULL,temp=0,hmax=0;
  struct llist *L=NULL;
  struct litem *l=NULL;
  char filename_base[1024],filename[2048],command[4096];
  FILE *fp=NULL;
  double maxdia=10000,xpos=0,ypos=0;
  int remove_flag=SUITE_7_CLEANUP,jpeg_flag=(ON_MY_COMPUTER ? 1 : 0),eps_flag=(ON_MY_COMPUTER ? 1 : 0);
  if (verbose){ printf(" %% [entering spearman_rank_correlation_histogram]\n");}
  h = histmake(total_bins,+1,-1);
  for (nt=0;nt<total_trials;nt++){
    if (verbose>1){ printf(" %% trial %d/%d\n",nt,total_trials);}
    nr0 = total_clusters*rand01; nr0 = periodize(nr0,0,total_clusters);
    nr1=total_1_records/total_clusters*rand01; nr1=periodize(nr1,0,total_1_records/total_clusters); nr1+=nr0*(total_1_records/total_clusters);
    nr2=total_2_records/total_clusters*rand01; nr2=periodize(nr2,0,total_2_records/total_clusters); nr2+=nr0*(total_2_records/total_clusters);
    L=llistmake();
    for (no=0;no<total_odors;no++){
      ra_temp = (double *) tcalloc(4,sizeof(double));
      ra_temp[0]=randn()*ra1[no+nr1*total_odors+1*total_odors*total_1_records]+ra1[no+nr1*total_odors+0*total_odors*total_1_records];
      ra_temp[1]=randn()*ra2[no+nr2*total_odors+1*total_odors*total_2_records]+ra2[no+nr2*total_odors+0*total_odors*total_2_records];
      litemadd(L,ra_temp);}
    if (verbose>1){ printf(" %% L:\n");l=L->first;while(l!=NULL){ raprintf((double *)l->item,"double",1,4,"%%");l=l->child;}}
    llistsort(L->first,L->last,L->length,&double_compare);
    l=L->first;no=0;while(l!=NULL){ 
      ra_temp=(double *)l->item;
      temp=ra_temp[0]; ra_temp[0]=ra_temp[1]; ra_temp[1]=temp;
      ra_temp[2]=no;
      no++;l=l->child;}
    llistsort(L->first,L->last,L->length,&double_compare);
    l=L->first;no=0;while(l!=NULL){ 
      ra_temp=(double *)l->item;
      ra_temp[3]=no;
      no++;l=l->child;}
    if (verbose>1){ printf(" %% L:\n");l=L->first;while(l!=NULL){ raprintf((double *)l->item,"double",1,4,"%%");l=l->child;}}
    temp=0;l=L->first;while(l!=NULL){ ra_temp=(double *)l->item; temp += pow(ra_temp[2]-ra_temp[3],2); l=l->child;}
    llisttfree2(L);
    histadd(h,1-(6.0*temp)/(double)(total_odors*(total_odors*total_odors-1)),1);}
  if (verbose){ histprintf(h,"correlation: ");}
  stats("double",h->data,h->nbins,&hmax,NULL,NULL,NULL);
  sprintf(filename_base,"%s_rank_correlation",fgvn);
  sprintf(filename,"%s.fig",filename_base);
  if ((fp=fopen(filename,"w"))==NULL){ printf(" %% warning! not open in %s spearman_rank_correlation_histogram\n",filename); fp=stdout;}
  fprintf(fp,"%s",FIG_PREAMBLE);
  for (nbin=0;nbin<h->nbins;nbin++){
    fprintf(fp,"2 1 0 1 %d 7 %d 0 -1 0.000 0 0 -1 0 0 %d\n\t",/*color*/0,/*depth*/999,/*npoints*/5);
    xpos = (double)(nbin+0)/(double)h->nbins; ypos = 0;
    fprintf(fp," %d %d",(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos));
    xpos = (double)(nbin+0)/(double)h->nbins; ypos = h->data[nbin]/hmax;
    fprintf(fp," %d %d",(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos));
    xpos = (double)(nbin+1)/(double)h->nbins; ypos = h->data[nbin]/hmax;
    fprintf(fp," %d %d",(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos));
    xpos = (double)(nbin+1)/(double)h->nbins; ypos = 0;
    fprintf(fp," %d %d",(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos));
    xpos = (double)(nbin+0)/(double)h->nbins; ypos = 0;
    fprintf(fp," %d %d",(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos));
    fprintf(fp,"\n");}
  fprintf(fp,"4 0 0 %d 0 %d %d 0.0000 4 %d %d %d %d %d\\001\n",/*depth*/1,/*font*/12,/*point*/24,/*textheight*/(int)(24*13.5),/*textwidth*/(int)(24*9*9),/*xpos*/(int)(0.1*maxdia),/*ypos*/(int)maxdia-(int)(1.0*maxdia),total_trials);
  if (fp!=stdout){ fclose(fp);fp=NULL;}
  if (verbose){ printf(" %% making rank plot...");}
  if (jpeg_flag){ sprintf(command,"fig2dev -Ljpeg -q %d %s.fig %s.jpg;",/*quality*/5,filename_base,filename_base); system(command);}
  if (eps_flag){ sprintf(command,"fig2dev -Leps %s.fig %s.eps;",filename_base,filename_base); system(command);}
  if (remove_flag){ sprintf(command,"rm %s.fig;",filename_base); system(command);}
  histtfree(h);h=NULL;
  if (verbose){ printf(" %% [finishing spearman_rank_correlation_histogram]\n");}
}

void suite_7_power_process_helper(char *gs2)
{
  int verbose=GLOBAL_verbose;
  char dirname[1024],fpre[512],fpos[512],filename[2048],filename_base[1536],command[3072];
  int total_trials=0,total_concentrations=0,total_odors=0,total_bicucullines=0,total_length=0,total_records=0,temp_flag=0;
  int ntrial=0,nconcentration=0,nodor=0,nbicuculline=0,nlength=0,nrecord=0;
  int nt=0,nv=0,trials_found=0;
  double *ra=NULL,*ra_reliability=NULL,*ra_orn_rank=NULL,*ra_pn_rank=NULL,**pra=NULL;
  FILE *fp=NULL;
  double maxdia=10000,xpos=0,ypos=0;
  int remove_flag=SUITE_7_CLEANUP;
  int jpeg_flag=(ON_MY_COMPUTER ? 1 : 0),eps_flag=(ON_MY_COMPUTER ? 1 : 0),pnm_flag=(!remove_flag);
  //int raster_plot_flag=0,reliability_plot_flag=1,rank_plot_flag=1,discriminability_plot_flag=1;
  int raster_plot_flag=0,reliability_plot_flag=0,rank_plot_flag=1,discriminability_plot_flag=0;
  int bin_size=0,bin_overlap=0,bin_skip=0,total_bins=0,nbin=0;
  double bins_per_second=0;
  int total_glomeruli=GLOBAL_NCLUSTERS,ngli=0;
  double temp_max=0,temp_min=0;
  char temp_char[32];
  int temp_tab=0,temp_bin=0,total_orn_records=0,total_pn_records=0;
  struct llist **Lra=NULL;
  double *temp=NULL,temp_dt=0;
  if (verbose){ printf(" %% [entering suite_7_power_process_helper] gs2 %s\n",gs2);}
  sprintf(dirname,"./dir_temp_suite_7_%s",gs2);
  if (verbose){ printf(" %% checking for directory %s\n",dirname);}
  if (!checktofind(dirname)){ printf(" %% warning! cannot find %s in suite_7_power_process_helper\n",dirname);}
  else /* directory exists */{
    if (verbose){ printf(" %% found %s\n",dirname);}
    for (nt=0;nt<GLOBAL_NTYPES;nt++){ for (nv=0;nv<GLOBAL_NVARS;nv++){
      if (verbose){ printf(" %% checking for %s, %s\n",GLOBAL_TYPENAMES[nt],GLOBAL_VARNAMES[nv]);}
      sprintf(filename,"%s/power_deposit_type_%s_var_%s_instance_%d_concentration_%d_odor_%d_bicuculline_%d_%s",dirname,GLOBAL_TYPENAMES[nt],GLOBAL_VARNAMES[nv],0,0,0,0,gs2);
      if (checktofind(filename)){
	sprintf(fpre,"%s/power_deposit_type_%s_var_%s_instance_",dirname,GLOBAL_TYPENAMES[nt],GLOBAL_VARNAMES[nv]);
	sprintf(fpos,"_concentration_%d_odor_%d_bicuculline_%d_%s",0,0,0,gs2);
	total_trials=0; checktofind_howmany(fpre,&total_trials,fpos);
	sprintf(fpre,"%s/power_deposit_type_%s_var_%s_instance_%d_concentration_",dirname,GLOBAL_TYPENAMES[nt],GLOBAL_VARNAMES[nv],0);
	sprintf(fpos,"_odor_%d_bicuculline_%d_%s",0,0,gs2);
	total_concentrations=0; checktofind_howmany(fpre,&total_concentrations,fpos);
	sprintf(fpre,"%s/power_deposit_type_%s_var_%s_instance_%d_concentration_%d_odor_",dirname,GLOBAL_TYPENAMES[nt],GLOBAL_VARNAMES[nv],0,0);
	sprintf(fpos,"_bicuculline_%d_%s",0,gs2);
	total_odors=0; checktofind_howmany(fpre,&total_odors,fpos);
	sprintf(fpre,"%s/power_deposit_type_%s_var_%s_instance_%d_concentration_%d_odor_%d_bicuculline_",dirname,GLOBAL_TYPENAMES[nt],GLOBAL_VARNAMES[nv],0,0,0);
	sprintf(fpos,"_%s",gs2);
	total_bicucullines=0; checktofind_howmany(fpre,&total_bicucullines,fpos);
	ra = (double *) raread(filename,&temp_flag,&total_length,&total_records);tfree(ra);ra=NULL;
	if (verbose){ printf(" %% found %s,%s, total_trials %d total_concentrations %d total_odors %d total_bicucullines %d, total_length %d, total_records %d\n",GLOBAL_TYPENAMES[nt],GLOBAL_VARNAMES[nv],total_trials,total_concentrations,total_odors,total_bicucullines,total_length,total_records);}
	if (raster_plot_flag && nv==VARNAME_REGISTRY_wilson_spike_flag){
	  if (verbose){ printf(" %% plotting rasters for type %s, var %s\n",GLOBAL_TYPENAMES[nt],GLOBAL_VARNAMES[nv]);}
	  for (nbicuculline=0;nbicuculline<total_bicucullines;nbicuculline++){ for (nodor=0;nodor<total_odors;nodor++){ for (nconcentration=0;nconcentration<total_concentrations;nconcentration++){
	    if (verbose){ printf(" %% currently trying nbicuculline %d nodor %d nconcentration %d... ",nbicuculline,nodor,nconcentration);}
	    sprintf(filename_base,"%s/power_deposit_type_%s_var_%s_raster_concentration_%d_odor_%d_bicuculline_%d_%s",dirname,GLOBAL_TYPENAMES[nt],GLOBAL_VARNAMES[nv],nconcentration,nodor,nbicuculline,gs2);
	    sprintf(filename,"%s.fig",filename_base);
	    if ((fp=fopen(filename,"w"))==NULL){ printf(" %% warning! couldn't open %s in suite_7_power_process_helper\n",filename); fp=stdout;}
	    fprintf(fp,"%s",FIG_PREAMBLE);
	    for (nrecord=0;nrecord<total_records;nrecord++){
	      fprintf(fp,"2 1 0 1 %d 7 %d 0 -1 0.000 0 0 -1 0 0 %d\n\t",/*color*/11,/*depth*/999,/*npoints*/5);
	      xpos = 0; ypos = (double)(nrecord+0)/(double)(total_records/total_glomeruli);
	      fprintf(fp," %d %d",(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos));
	      xpos = 1; ypos = (double)(nrecord+0)/(double)(total_records/total_glomeruli);
	      fprintf(fp," %d %d",(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos));
	      xpos = 1; ypos = (double)(nrecord+1)/(double)(total_records/total_glomeruli);
	      fprintf(fp," %d %d",(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos));
	      xpos = 0; ypos = (double)(nrecord+1)/(double)(total_records/total_glomeruli);
	      fprintf(fp," %d %d",(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos));
	      xpos = 0; ypos = (double)(nrecord+0)/(double)(total_records/total_glomeruli);
	      fprintf(fp," %d %d",(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos));
	      fprintf(fp,"\n");}
	    trials_found=0;
	    for (ntrial=0;ntrial<total_trials;ntrial++){
	      sprintf(filename,"%s/power_deposit_type_%s_var_%s_instance_%d_concentration_%d_odor_%d_bicuculline_%d_%s",dirname,GLOBAL_TYPENAMES[nt],GLOBAL_VARNAMES[nv],ntrial,nconcentration,nodor,nbicuculline,gs2);
	      if (checktofind(filename)){
		ra = (double *) raread(filename,&temp_flag,&nlength,&nrecord);
		assert(nlength==total_length && nrecord==total_records);
		for (nlength=0;nlength<total_length;nlength++){ for (nrecord=0;nrecord<total_records;nrecord++){
		  if (ra[nlength+nrecord*total_length]>exp(-1.5)){
		    fprintf(fp,"2 1 0 1 %d 7 %d 0 -1 0.000 0 0 -1 0 0 %d\n\t",/*color*/0,/*depth*/(nrecord%999)+1,/*npoints*/2);
		    xpos = ((double)nlength + 0.5)/(double)total_length;
		    ypos = ((double)(trials_found + nrecord*total_trials)+0.0)/(double)(total_trials*(total_records/total_glomeruli));
		    fprintf(fp," %d %d",(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos));
		    ypos = ((double)(trials_found + nrecord*total_trials)+1.0)/(double)(total_trials*(total_records/total_glomeruli));
		    fprintf(fp," %d %d",(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos));
		    fprintf(fp,"\n");}}}
		tfree(ra);ra=NULL;
		trials_found += 1;}}
	    if (verbose){ printf(" %d total trials found\n",trials_found);}
	    if (fp!=stdout){ fclose(fp); fp=NULL;}
	    if (trials_found>0){
	      if (verbose){ printf(" %% making raster plot...");}
	      if (jpeg_flag){ sprintf(command,"fig2dev -Ljpeg -q %d %s.fig %s.jpg;",/*quality*/5,filename_base,filename_base); system(command);}
	      if (eps_flag){ sprintf(command,"fig2dev -Leps %s.fig %s.eps;",filename_base,filename_base); system(command);}
	      if (remove_flag){ sprintf(command,"rm %s.fig;",filename_base); system(command);}
	      if (verbose){ printf(" done plotting\n");}}}}}}
	if (reliability_plot_flag && nv==VARNAME_REGISTRY_wilson_spike_flag){
	  if (verbose){ printf(" %% plotting reliability for type %s, var %s\n",GLOBAL_TYPENAMES[nt],GLOBAL_VARNAMES[nv]);}
	  for (nbicuculline=0;nbicuculline<total_bicucullines;nbicuculline++){ for (nodor=0;nodor<total_odors;nodor++){ for (nconcentration=0;nconcentration<total_concentrations;nconcentration++){
	    if (verbose){ printf(" %% currently trying nbicuculline %d nodor %d nconcentration %d... ",nbicuculline,nodor,nconcentration);}
	    bin_size = 50; bin_overlap = minimum(25,bin_size-1); bin_skip = bin_size-bin_overlap;
	    total_bins = total_length/bin_skip+1;
	    ra_reliability = (double *) tcalloc(total_bins*total_records*3,sizeof(double));
	    trials_found=0;
	    for (ntrial=0;ntrial<total_trials;ntrial++){
	      sprintf(filename,"%s/power_deposit_type_%s_var_%s_instance_%d_concentration_%d_odor_%d_bicuculline_%d_%s",dirname,GLOBAL_TYPENAMES[nt],GLOBAL_VARNAMES[nv],ntrial,nconcentration,nodor,nbicuculline,gs2);
	      if (checktofind(filename)){
		ra = (double *) raread(filename,&temp_flag,&nlength,&nrecord);
		assert(nlength==total_length && nrecord==total_records);
		for (nrecord=0;nrecord<total_records;nrecord++){
		  for (nbin=0;nbin<total_bins;nbin++){ ra_reliability[nbin+nrecord*total_bins+0*total_bins*total_records]=0;}
		  for (nlength=0;nlength<total_length;nlength++){ 
		    for (nbin=nlength/bin_skip;nbin>=maximum(0,(nlength-bin_size)/bin_skip);nbin--){
		      ra_reliability[nbin+nrecord*total_bins+0*total_bins*total_records] += ra[nlength+nrecord*total_length];}}
		  for (nbin=0;nbin<total_bins;nbin++){ 
		    ra_reliability[nbin+nrecord*total_bins+1*total_bins*total_records] += ra_reliability[nbin+nrecord*total_bins+0*total_bins*total_records];
		    ra_reliability[nbin+nrecord*total_bins+2*total_bins*total_records] += pow(ra_reliability[nbin+nrecord*total_bins+0*total_bins*total_records],2);}}
		tfree(ra);ra=NULL;
		trials_found += 1;}}
	    if (verbose){ printf(" %d total trials found\n",trials_found);}
	    for (nrecord=0;nrecord<total_records;nrecord++){ for (nbin=0;nbin<total_bins;nbin++){
	      ra_reliability[nbin+nrecord*total_bins+1*total_bins*total_records] /= trials_found;
	      ra_reliability[nbin+nrecord*total_bins+2*total_bins*total_records] = sqrt(ra_reliability[nbin+nrecord*total_bins+2*total_bins*total_records]/(double)trials_found - pow(ra_reliability[nbin+nrecord*total_bins+1*total_bins*total_records],2));
	      ra_reliability[nbin+nrecord*total_bins+0*total_bins*total_records] = ra_reliability[nbin+nrecord*total_bins+2*total_bins*total_records]/ra_reliability[nbin+nrecord*total_bins+1*total_bins*total_records];
	      if (!finite(ra_reliability[nbin+nrecord*total_bins+0*total_bins*total_records])){ ra_reliability[nbin+nrecord*total_bins+0*total_bins*total_records] = 0;}}}
	    if (trials_found>0){
	      for (temp_tab=0;temp_tab<3;temp_tab++){
		switch (temp_tab){
		case 0: sprintf(temp_char,"cv");temp_max=4;temp_min=0; break;
		case 1: sprintf(temp_char,"mean");temp_max=(double)bin_size/(double)8/(double)2;temp_min=0; break;
		case 2: sprintf(temp_char,"stdev");temp_max=1.5;temp_min=0; break;
		default: break;}
		sprintf(filename_base,"%s/power_deposit_type_%s_var_%s_firingrate_reliability_%s_concentration_%d_odor_%d_bicuculline_%d_%s",dirname,GLOBAL_TYPENAMES[nt],GLOBAL_VARNAMES[nv],temp_char,nconcentration,nodor,nbicuculline,gs2);
		sprintf(filename,"%s.pnm",filename_base);
		if (pnm_flag){ WritePNMfile_color(&(ra_reliability[0+0*total_bins+temp_tab*total_bins*total_records]),total_bins,total_records,temp_max,temp_min,filename,7);}
		sprintf(filename,"%s.fig",filename_base);
		if ((fp=fopen(filename,"w"))==NULL){ printf(" %% warning! couldn't open %s in suite_7_power_process_helper\n",filename); fp=stdout;}
		fprintf(fp,"%s",FIG_PREAMBLE);
		for (ngli=0;ngli<total_glomeruli;ngli++){
		  fprintf(fp,"2 1 0 1 %d 7 %d 0 -1 0.000 0 0 -1 0 0 %d\n\t",/*color*/11,/*depth*/999,/*npoints*/5);
		  xpos = 0; ypos = (double)(ngli+0)/(double)total_glomeruli;
		  fprintf(fp," %d %d",(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos));
		  xpos = 1; ypos = (double)(ngli+0)/(double)total_glomeruli;
		  fprintf(fp," %d %d",(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos));
		  xpos = 1; ypos = (double)(ngli+1)/(double)total_glomeruli;
		  fprintf(fp," %d %d",(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos));
		  xpos = 0; ypos = (double)(ngli+1)/(double)total_glomeruli;
		  fprintf(fp," %d %d",(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos));
		  xpos = 0; ypos = (double)(ngli+0)/(double)total_glomeruli;
		  fprintf(fp," %d %d",(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos));
		  fprintf(fp,"\n");}
		for (nrecord=0;nrecord<total_records;nrecord++){
		  ngli = nrecord/(total_records/total_glomeruli);
		  fprintf(fp,"2 1 0 1 %d 7 %d 0 -1 0.000 0 0 -1 0 0 %d\n\t",/*color*/(nrecord)%32,/*depth*/((nrecord)%999)+1,/*npoints*/total_bins);
		  for (nbin=0;nbin<total_bins;nbin++){
		    xpos = ((double)nbin + 0.5)/(double)total_bins;
		    ypos = (ra_reliability[nbin+nrecord*total_bins+temp_tab*total_bins*total_records]-temp_min)/(temp_max-temp_min);
		    ypos = maximum(0,minimum(1,ypos));
		    ypos *= 1.0/(double)total_glomeruli; ypos += (double)ngli/(double)total_glomeruli;
		    fprintf(fp," %d %d",(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos));}
		  fprintf(fp,"\n");}
		fprintf(fp,"4 0 0 %d 0 %d %d 0.0000 4 %d %d %d %d %0.16f\\001\n",/*depth*/1,/*font*/12,/*point*/24,/*textheight*/(int)(24*13.5),/*textwidth*/(int)(24*9*18),/*xpos*/(int)(0.1*maxdia),/*ypos*/(int)maxdia-(int)(1.0*maxdia),temp_max);
		fprintf(fp,"4 0 0 %d 0 %d %d 0.0000 4 %d %d %d %d %0.16f\\001\n",/*depth*/1,/*font*/12,/*point*/24,/*textheight*/(int)(24*13.5),/*textwidth*/(int)(24*9*18),/*xpos*/(int)(0.1*maxdia),/*ypos*/(int)maxdia-(int)(0.0*maxdia),temp_min);
		if (fp!=stdout){ fclose(fp);fp=NULL;
		if (jpeg_flag){ sprintf(command,"fig2dev -Ljpeg -q %d %s.fig %s.jpg;",/*quality*/5,filename_base,filename_base); system(command);}
		if (eps_flag){ sprintf(command,"fig2dev -Leps %s.fig %s.eps;",filename_base,filename_base); system(command);}
		if (remove_flag){ sprintf(command,"rm %s.fig;",filename_base); system(command);}}}}
	    tfree(ra_reliability);ra_reliability=NULL;}}}}
      }
    }}
    if (rank_plot_flag){
      for (nbicuculline=0;nbicuculline<total_bicucullines;nbicuculline++){ for (nconcentration=/*0*/1;nconcentration<total_concentrations;nconcentration++){
	if (verbose){ printf(" %% plotting rank correlation for nconcentration %d nbicuculline %d\n",nconcentration,nbicuculline);}
	for (temp_tab=0;temp_tab<=1;temp_tab++){
	  nv = VARNAME_REGISTRY_wilson_spike_flag;
	  switch(temp_tab){ 
	  case 0: nt=TYPENAME_REGISTRY_wilson_ORN; pra = &ra_orn_rank; 
	    sprintf(filename,"%s/power_deposit_type_%s_var_%s_instance_%d_concentration_%d_odor_%d_bicuculline_%d_%s",dirname,GLOBAL_TYPENAMES[nt],GLOBAL_VARNAMES[nv],0,0,0,0,gs2); if (checktofind(filename)){ ra = (double *) raread(filename,&temp_flag,&total_length,&total_records);tfree(ra);ra=NULL; total_orn_records=total_records;}
	    break;
	  case 1: nt=TYPENAME_REGISTRY_wilson_PN; pra = &ra_pn_rank; 
	    sprintf(filename,"%s/power_deposit_type_%s_var_%s_instance_%d_concentration_%d_odor_%d_bicuculline_%d_%s",dirname,GLOBAL_TYPENAMES[nt],GLOBAL_VARNAMES[nv],0,0,0,0,gs2); if (checktofind(filename)){ ra = (double *) raread(filename,&temp_flag,&total_length,&total_records);tfree(ra);ra=NULL; total_pn_records=total_records;}
	    break;
	  default: break;}
	  *pra = (double *) tcalloc(total_odors*total_records*2/*mean,stdev*/,sizeof(double));
	  for (nodor=0;nodor<total_odors;nodor++){
	    if (verbose){ printf(" %% working on odor %d type %s\n",nodor,GLOBAL_TYPENAMES[nt]);}
	    bin_size = 50; bin_overlap = minimum(25,bin_size-1); bin_skip = bin_size-bin_overlap; bins_per_second = 1024/(double)bin_size;
	    total_bins = total_length/bin_skip+1;
	    ra_reliability = (double *) tcalloc(total_bins*total_records*3,sizeof(double));
	    trials_found=0;
	    for (ntrial=0;ntrial<total_trials;ntrial++){
	      sprintf(filename,"%s/power_deposit_type_%s_var_%s_instance_%d_concentration_%d_odor_%d_bicuculline_%d_%s",dirname,GLOBAL_TYPENAMES[nt],GLOBAL_VARNAMES[nv],ntrial,nconcentration,nodor,nbicuculline,gs2);
	      if (checktofind(filename)){
		ra = (double *) raread(filename,&temp_flag,&nlength,&nrecord);
		assert(nlength==total_length && nrecord==total_records);
		for (nrecord=0;nrecord<total_records;nrecord++){
		  for (nbin=0;nbin<total_bins;nbin++){ ra_reliability[nbin+nrecord*total_bins+0*total_bins*total_records]=0;}
		  for (nlength=0;nlength<total_length;nlength++){ 
		    for (nbin=nlength/bin_skip;nbin>=maximum(0,(nlength-bin_size)/bin_skip);nbin--){
		      ra_reliability[nbin+nrecord*total_bins+0*total_bins*total_records] += ra[nlength+nrecord*total_length];}}
		  for (nbin=0;nbin<total_bins;nbin++){ 
		    ra_reliability[nbin+nrecord*total_bins+1*total_bins*total_records] += ra_reliability[nbin+nrecord*total_bins+0*total_bins*total_records];
		    ra_reliability[nbin+nrecord*total_bins+2*total_bins*total_records] += pow(ra_reliability[nbin+nrecord*total_bins+0*total_bins*total_records],2);}}
		tfree(ra);ra=NULL;
		trials_found += 1;}}
	    if (verbose){ printf(" %d total trials found\n",trials_found);}
	    for (nrecord=0;nrecord<total_records;nrecord++){ for (nbin=0;nbin<total_bins;nbin++){
	      ra_reliability[nbin+nrecord*total_bins+1*total_bins*total_records] /= trials_found;
	      ra_reliability[nbin+nrecord*total_bins+2*total_bins*total_records] = sqrt(ra_reliability[nbin+nrecord*total_bins+2*total_bins*total_records]/(double)trials_found - pow(ra_reliability[nbin+nrecord*total_bins+1*total_bins*total_records],2));
	      ra_reliability[nbin+nrecord*total_bins+0*total_bins*total_records] = ra_reliability[nbin+nrecord*total_bins+2*total_bins*total_records]/ra_reliability[nbin+nrecord*total_bins+1*total_bins*total_records];
	      if (!finite(ra_reliability[nbin+nrecord*total_bins+0*total_bins*total_records])){ ra_reliability[nbin+nrecord*total_bins+0*total_bins*total_records] = 0;}}}
	    for (nrecord=0;nrecord<total_records;nrecord++){
	      temp_bin=0; for (nbin=0;nbin<total_bins;nbin++){ if (ra_reliability[nbin+nrecord*total_bins+1*total_bins*total_records]>ra_reliability[temp_bin+nrecord*total_bins+1*total_bins*total_records] || (ra_reliability[nbin+nrecord*total_bins+1*total_bins*total_records]==ra_reliability[temp_bin+nrecord*total_bins+1*total_bins*total_records] && ra_reliability[nbin+nrecord*total_bins+2*total_bins*total_records]<ra_reliability[temp_bin+nrecord*total_bins+2*total_bins*total_records])){ temp_bin = nbin;}}
	      if (verbose>1){ printf(" %% odor %d record %d, max %f found at index %d of %d\n",nodor,nrecord,ra_reliability[temp_bin+nrecord*total_bins+1*total_bins*total_records],temp_bin,total_bins);}
	      (*pra)[nodor+nrecord*total_odors + 0*total_odors*total_records] = bins_per_second*ra_reliability[maximum(0,temp_bin-1)+nrecord*total_bins+1*total_bins*total_records]+ra_reliability[minimum(total_bins-1,temp_bin+1)+nrecord*total_bins+1*total_bins*total_records];
	      (*pra)[nodor+nrecord*total_odors + 1*total_odors*total_records] = sqrt(bins_per_second)*ra_reliability[maximum(0,temp_bin-1)+nrecord*total_bins+2*total_bins*total_records]+ra_reliability[minimum(total_bins-1,temp_bin+1)+nrecord*total_bins+2*total_bins*total_records];}
	    tfree(ra_reliability);ra_reliability=NULL;}}
	/* sprintf(filename,"%s/power_deposit_rank_orn_orn_concentration_%d_bicuculline_%d_%s",dirname,nconcentration,nbicuculline,gs2); spearman_rank_correlation_histogram(total_glomeruli,total_odors,ra_orn_rank,total_orn_records,ra_orn_rank,total_orn_records,filename);  */
	/* sprintf(filename,"%s/power_deposit_rank_orn_pn_concentration_%d_bicuculline_%d_%s",dirname,nconcentration,nbicuculline,gs2); spearman_rank_correlation_histogram(total_glomeruli,total_odors,ra_orn_rank,total_orn_records,ra_pn_rank,total_pn_records,filename);  */
	/* sprintf(filename,"%s/power_deposit_rank_pn_pn_concentration_%d_bicuculline_%d_%s",dirname,nconcentration,nbicuculline,gs2); spearman_rank_correlation_histogram(total_glomeruli,total_odors,ra_pn_rank,total_pn_records,ra_pn_rank,total_pn_records,filename); */
	/* sprintf(filename,"%s/power_deposit_orn_pn_firingrate_concentration_%d_bicuculline_%d_%s",dirname,nconcentration,nbicuculline,gs2); suite_7_power_process_helper_1(total_glomeruli,total_odors,ra_orn_rank,total_orn_records,ra_pn_rank,total_pn_records,filename); */
	sprintf(filename,"%s/power_deposit_orn_pn_firingrate_ver2_concentration_%d_bicuculline_%d_%s",dirname,nconcentration,nbicuculline,gs2); suite_7_power_process_helper_2(total_glomeruli,total_odors,ra_orn_rank,total_orn_records,ra_pn_rank,total_pn_records,filename);
	tfree(ra_orn_rank);ra_orn_rank=NULL;tfree(ra_pn_rank);ra_pn_rank=NULL;pra=NULL;}}}
 
    if (discriminability_plot_flag){
      nt=TYPENAME_REGISTRY_wilson_PN;
      nv = VARNAME_REGISTRY_wilson_spike_flag;
      sprintf(filename,"%s/power_deposit_type_%s_var_%s_instance_%d_concentration_%d_odor_%d_bicuculline_%d_%s",dirname,GLOBAL_TYPENAMES[nt],GLOBAL_VARNAMES[nv],0,0,0,0,gs2); if (checktofind(filename)){ ra = (double *) raread(filename,&temp_flag,&total_length,&total_records);tfree(ra);ra=NULL; total_pn_records=total_records;}
      for (nbicuculline=0;nbicuculline<total_bicucullines;nbicuculline++){ for (nconcentration=0;nconcentration<total_concentrations;nconcentration++){ for (nodor=0;nodor<total_odors;nodor++){ 
	Lra = (struct llist **) tcalloc(total_pn_records,sizeof(struct llist *));
	for (nrecord=0;nrecord<total_pn_records;nrecord++){ Lra[nrecord]=llistmake();}
	for (ntrial=0;ntrial<total_trials;ntrial++){
	  if (verbose){ printf(" %% collecting firing events for discriminability correlation for ntrial %d nodor %d nconcentration %d nbicuculline %d\n",ntrial,nodor,nconcentration,nbicuculline);}
	  sprintf(filename,"%s/power_deposit_type_%s_var_%s_instance_%d_concentration_%d_odor_%d_bicuculline_%d_%s",dirname,GLOBAL_TYPENAMES[nt],GLOBAL_VARNAMES[nv],ntrial,nconcentration,nodor,nbicuculline,gs2);
	  if (checktofind(filename)){
	    ra = (double *) raread(filename,&temp_flag,&nlength,&nrecord);assert(nlength==total_length && nrecord==total_records);
	    for (nrecord=0;nrecord<total_pn_records;nrecord++){ for (nlength=1;nlength<total_length;nlength++){
	      temp_dt=0; 
	      if (spike_flag_to_spike_time(ra[nlength-1 + nrecord*total_length],ra[nlength + nrecord*total_length],1.0,TAU_[VARNAME_REGISTRY_wilson_spike_flag],1.0,&temp_dt)){
		temp = (double *) tcalloc(1,sizeof(double)); 
		*temp = total_length*ntrial + nlength-1 + temp_dt;
		litemadd(Lra[nrecord],temp);}}}
	    tfree(ra);ra=NULL;}}
	sprintf(filename,"%s/Lra_t_%s_c_%d_o_%d_b_%d_%s",dirname,GLOBAL_TYPENAMES[nt],nconcentration,nodor,nbicuculline,gs2);
	if (verbose){ printf(" %% dumping llistra as %s, which is an abbreviation of %s/llistra_deposit_type_%s_var_%s_concentration_%d_odor_%d_bicuculline_%d_%s \n",filename,dirname,GLOBAL_TYPENAMES[nt],GLOBAL_VARNAMES[nv],nconcentration,nodor,nbicuculline,gs2);}
	llistra2dump(Lra,total_records,filename);
	for (nrecord=0;nrecord<total_pn_records;nrecord++){ llisttfree2(Lra[nrecord]);Lra[nrecord]=NULL;}
	tfree(Lra);Lra=NULL;}}
      sprintf(filename,"%s/Lra_t_%s_b_%d_%s.in",dirname,GLOBAL_TYPENAMES[nt],nbicuculline,gs2);
      if ((fp=fopen(filename,"w"))==NULL){ printf(" %% warning! couldn't open %s in suite_7_power_process_helper\n",filename); fp=stdout;}
      fprintf(fp,"GLOBAL_verbose= 0;\n");
      fprintf(fp,"GLOBAL_STRING= Lra_t_%s_b_%d_%s;\n",GLOBAL_TYPENAMES[nt],nbicuculline,gs2);
      fprintf(fp,"GLOBAL_MODE= 4;\n");
      fprintf(fp,"GLOBAL_DATA_FORMAT= 1;\n");
      fprintf(fp,"GLOBAL_NFILES_TO_COMPARE= %d;\n",total_concentrations*total_odors);
      fprintf(fp,"GLOBAL_FILENAMES_TO_COMPARE= ");
      for (nconcentration=0;nconcentration<total_concentrations;nconcentration++){ for (nodor=0;nodor<total_odors;nodor++){ 
	fprintf(fp,"Lra_t_%s_c_%d_o_%d_b_%d_%s",GLOBAL_TYPENAMES[nt],nconcentration,nodor,nbicuculline,gs2);
	fprintf(fp,"%s",(nconcentration==total_concentrations-1 && nodor==total_odors-1) ? ";" : ",");}}
      fprintf(fp,"\n");
      fprintf(fp,"GLOBAL_TF= %d;\n",total_trials*total_length);
      fprintf(fp,"GLOBAL_PTREE_NREGIONS= %d;\n",total_pn_records);
      fprintf(fp,"GLOBAL_PTREE_NLEGS= %d;\n",2);
      fprintf(fp,"GLOBAL_PTREE_LEGTIME= %d;\n",32);
      fprintf(fp,"GLOBAL_PTREE_OBSTIME= %d;\n",total_length);
      fprintf(fp,"END= 0;\n");
      if (fp!=stdout){ fclose(fp);fp=NULL;}
      sprintf(command,"cd %s; nice ./ptree_test < Lra_t_%s_b_%d_%s.in; rm -rf ptree_obsdisthist*.fig; cd ..;",dirname,GLOBAL_TYPENAMES[nt],nbicuculline,gs2); system(command);
      } /* for (nbicuculline=0;nbicuculline<total_bicucullines;nbicuculline++) */ } /* if (discriminability_plot_flag) */

  }
}

void suite_7_power_deposit_helper(int instance,int concentration,int odor,int bicuculline,int time)
{
  int verbose=1;
  char dirname[1024],filename[1536],gs2[128],command[2048];
  struct power *p=GLOBAL_POWER;
  int nt=0,nrecords=0,nv=0,sttab=0,stlength=0,nr=0,nl=0,nl2=0;
  double *ra=NULL;
  if (verbose){ printf(" %% [entering suite_7_power_deposit_helper] trial %d concentration %d odor %d bicuculline %d time %d\n",instance,concentration,odor,bicuculline,time);}
  sprintf(gs2,"%srecord",GLOBAL_STRING_2);
  sprintf(dirname,"./dir_temp_suite_7_%s",gs2);
  if (verbose){ printf(" %% checking for directory %s\n",dirname);}
  if (!checktofind(dirname)){ if (verbose){ printf(" %% making %s\n",dirname);} sprintf(command,"mkdir %s;",dirname); system(command);}
  for (nt=0;nt<p->indexing_ntype_length;nt++){ 
    nrecords = p->Nra->lengthra[p->indexing_ntype_checkout[nt]];
    for (nv=0;nv<p->indexing_nvar_length;nv++){
      if (p->strarara[nt][nv]!=NULL){
	if (verbose){ printf(" %% depositing p->strarara[%s][%s]\n",GLOBAL_TYPENAMES[p->indexing_ntype_checkout[nt]],GLOBAL_VARNAMES[p->indexing_nvar_checkout[nv]]);}
	sttab = p->strarara[nt][nv][0]->tab; stlength = p->strarara[nt][nv][0]->length;
	ra = (double *) tcalloc(p->length*nrecords,sizeof(double));
	for (nr=0;nr<nrecords;nr++){ for (nl=0;nl<p->length;nl++){
	  nl2 = periodize(sttab-p->length + nl,0,stlength);
	  ra[nl + nr*p->length] = p->strarara[nt][nv][nr]->data[nl2];}}
	sprintf(filename,"%s/power_deposit_type_%s_var_%s_instance_%d_concentration_%d_odor_%d_bicuculline_%d_%s",dirname,GLOBAL_TYPENAMES[p->indexing_ntype_checkout[nt]],GLOBAL_VARNAMES[p->indexing_nvar_checkout[nv]],instance,concentration,odor,bicuculline,gs2);
	if (verbose){ printf(" %% writing ra %s...",filename);}
	radump(ra,"double",p->length,nrecords,filename);
	if (verbose){ printf(" freeing ra\n");}
	tfree(ra);ra=NULL;}}}
}

void system_monitor(struct neuronarray *Nra,double t,double DT)
{
  /* govern input and data-structures */
  int verbose=1;
  double time_dumpevery=0;
  char **fnamebase=NULL;
  char *gs2=GLOBAL_STRING_2;
  int nr1=0,nr2=0,nr3=0,nr4=0,nr=0,nc=0,tab=0,tab2=0,tab3=0,nv=0;
  double *timera=NULL;
  double pulse_waitfor=0,pulse_howlong=0,pulse_temp=0,pulse_time=0;
  double pulse_phase_max=0,pulse_phase_min=0,pulse_phase=0;
  int nt=0,ntimes=0,ninputs=0,nowtime=0,nowodor=0,nowconcentration=0,nowinstance=0,nowbicuculline=0,nowinput=0,nowreweight=0;
  int oldtime=0,oldodor=0,oldconcentration=0,oldinstance=0,oldbicuculline=0,oldinput=0,oldreweight=0;
  char tempchar[1024],command[2048];
  double **temprara=NULL,*tempra=NULL,*tempra2=NULL;
  struct clusterdatara *cd=NULL;
  struct ptree *p=NULL;
  int npoints=0;
  struct odor *o1=NULL;
  if (verbose>1){ printf(" %% [enter system_monitor] with t %0.2f and DT %0.2f\n",t,DT);}
  switch (SUITE_BOTHER){
  case 1:
    ntimes = SUITE_NINSTANCES*SUITE_NCONCENTRATIONS*SUITE_NODORS*SUITE_NBICUCULLINES + 1;
    ninputs = SUITE_NCONCENTRATIONS*SUITE_NODORS*SUITE_NBICUCULLINES;
    timera = (double *) tcalloc(ntimes,sizeof(double));
    for (nt=0;nt<ntimes;nt++){ timera[nt] = nt*1024*SUITE_NSECONDS;}
    time_dumpevery = 1024*SUITE_NSECONDS;
    pulse_waitfor = 2048;
    pulse_howlong = 512;
    if (t==0){ /* clear everything */
      if (verbose){ printf(" %% \n");}
      if (verbose){ printf(" %% initializing suite one\n");}
      if (verbose>2){ printf(" times are: ");raprintf(timera,"double",1,ntimes,"timera ");}
      GLOBAL_TF=0;
      SUITE_ONSET_TIME=-1;
      if (POWER_BOTHER!=0){ powertfree(GLOBAL_POWER); GLOBAL_POWER=NULL; POWER_BOTHER=0;}
      if (PTREE_BOTHER!=0){ ptreetfree(GLOBAL_PTREE); GLOBAL_PTREE=NULL; GLOBAL_PTREE_BITBYBIT=0; GLOBAL_PTREE_BITBYBIT_CURRENT_NUMBER=0; PTREE_BOTHER=0;}
      if (CLUSTERDATA_BOTHER!=0){ clusterdatatfree(GLOBAL_CDRA); GLOBAL_CDRA=NULL; CLUSTERDATA_BOTHER=0;}}
    else /* if (t>0) */{
      nowtime = 0; for (nt=0;nt<ntimes-1;nt++){ if (t>=timera[nt] && t<timera[nt+1]){ nowtime=nt;}} 
      if (t>=timera[ntimes-1]){ nowtime=ntimes-1;}
      nowbicuculline = nowtime/(SUITE_NINSTANCES*SUITE_NCONCENTRATIONS*SUITE_NODORS);
      nowodor = (nowtime-nowbicuculline*SUITE_NINSTANCES*SUITE_NCONCENTRATIONS*SUITE_NODORS)/(SUITE_NINSTANCES*SUITE_NCONCENTRATIONS);
      nowconcentration = (nowtime-nowbicuculline*SUITE_NINSTANCES*SUITE_NCONCENTRATIONS*SUITE_NODORS-nowodor*SUITE_NINSTANCES*SUITE_NCONCENTRATIONS)/(SUITE_NINSTANCES);
      nowinstance = (nowtime-nowbicuculline*SUITE_NINSTANCES*SUITE_NCONCENTRATIONS*SUITE_NODORS-nowodor*SUITE_NINSTANCES*SUITE_NCONCENTRATIONS-nowconcentration*SUITE_NINSTANCES)/(1);
      nowinput = nowconcentration + nowodor*SUITE_NCONCENTRATIONS + nowbicuculline*SUITE_NCONCENTRATIONS*SUITE_NODORS;
      oldtime = maximum(0,nowtime-1);
      oldbicuculline = oldtime/(SUITE_NINSTANCES*SUITE_NCONCENTRATIONS*SUITE_NODORS);
      oldodor = (oldtime-oldbicuculline*SUITE_NINSTANCES*SUITE_NCONCENTRATIONS*SUITE_NODORS)/(SUITE_NINSTANCES*SUITE_NCONCENTRATIONS);
      oldconcentration = (oldtime-oldbicuculline*SUITE_NINSTANCES*SUITE_NCONCENTRATIONS*SUITE_NODORS-oldodor*SUITE_NINSTANCES*SUITE_NCONCENTRATIONS)/(SUITE_NINSTANCES);
      oldinstance = (oldtime-oldbicuculline*SUITE_NINSTANCES*SUITE_NCONCENTRATIONS*SUITE_NODORS-oldodor*SUITE_NINSTANCES*SUITE_NCONCENTRATIONS-oldconcentration*SUITE_NINSTANCES)/(1);
      oldinput = oldconcentration + oldodor*SUITE_NCONCENTRATIONS + oldbicuculline*SUITE_NCONCENTRATIONS*SUITE_NODORS;
      if (nowtime<ntimes){
	if (fabs(SUITE_ONSET_TIME-timera[nowtime])>0.1){
	  if (verbose){ printf(" %% time %f, entering %d_th phase of suite one, bicuculline %d odor %d concentration %d instance %d\n",GLOBAL_time,nowtime,nowbicuculline,nowodor,nowconcentration,nowinstance);}
	  SUITE_ONSET_TIME = timera[nowtime]; OUTPUT_DUMP_EVERY = 0; /* don't bother dumping data */
	  if (verbose){ printf(" %% applying bicuculline...");}
	  GLOBAL_CS_SRA_SCALE_[GLOBAL_INDEXING_REFILE_sra[VARNAME_REGISTRY_s_G]]=pow(0.01,nowbicuculline/maximum(1,SUITE_NBICUCULLINES-1));
	  if (verbose){ printf(" s_G scale %f \n",GLOBAL_CS_SRA_SCALE_[GLOBAL_INDEXING_REFILE_sra[VARNAME_REGISTRY_s_G]]);}
	  odorplusequals(GLOBAL_ODORra,nowodor,SUITE_ODORra_BACON);
	  if (verbose){ printf("%% %% GLOBAL_ODORra: ");odorfprintf(stdout,GLOBAL_ODORra);}
	  INPUT_CONTRAST = pow(SUITE_INPUT_CONTRAST,nowconcentration);
	  if (verbose){ printf(" %% %% setting INPUT_CONTRAST=%0.2f\n",INPUT_CONTRAST);}
	  if (SUITE_PTREE_BOTHER){
	    if (GLOBAL_PTREE==NULL){
	      if (verbose){ printf(" %% %% making ptree\n");}
	      PTREE_BOTHER=1; GLOBAL_PTREE_BITBYBIT=1; GLOBAL_PTREE_BITBYBIT_CURRENT_NUMBER=0; GLOBAL_PTREE_REGION_TYPE=1;
	      GLOBAL_PTREE = ptreemake(GLOBAL_PTREE_NREGIONS,GLOBAL_PTREE_EVENT_WITHIN,GLOBAL_PTREE_EVENT_THRESHOLD,GLOBAL_PTREE_REGION_TYPE,GLOBAL_PTREE_NLEGS,GLOBAL_PTREE_LEGTIME);}
	    else /* if (GLOBAL_PTREE!=NULL) */{
	      if (verbose){ printf(" %% %% updating odh->hra[%d]... ",oldinput);}
	      pnode_obs2dist_starter(oldinput,ninputs,0,NULL,GLOBAL_PTREE->postree,-1,0,1);
	      if (SUITE_BITBYBIT_RECORD){
		GLOBAL_PTREE_BITBYBIT_CURRENT_NUMBER=oldinstance;
		if (verbose){ printf("dumping ptree... ");}
		ptreerate(GLOBAL_PTREE); 
		sprintf(tempchar,"input%d",oldinput); 
		ptreedump_starter(GLOBAL_PTREE,tempchar,2,0,0,0,+1,-1,0);}
	      if (verbose){ printf("clearing ptree... ");}
	      pnodeclear_starter(NULL,GLOBAL_PTREE->pretree); pnodeclear_starter(NULL,GLOBAL_PTREE->postree); GLOBAL_PTREE->total_time=0;
	      if (verbose){ printf("\n");}}}
	  if (SUITE_CLUSTERDATA_BOTHER){
	    if (GLOBAL_CDRA==NULL){
	      if (verbose){ printf(" %% %% making clusterdatara\n");}
	      CLUSTERDATA_BOTHER=1; GLOBAL_PTREE_BITBYBIT=1; GLOBAL_PTREE_BITBYBIT_CURRENT_NUMBER=0; GLOBAL_PTREE_REGION_TYPE=1; GLOBAL_CDRA = clusterdataramake(GLOBAL_Nra->gli,1,/* GLOBAL_POWER_LENGTH */1024*SUITE_NSECONDS,GLOBAL_POWER_INDEXING_NVAR_LENGTH,GLOBAL_POWER_INDEXING_NVAR_CHECKOUT,GLOBAL_Nra->nvars,GLOBAL_POWER_INDEXING_NVAR_REFILE,GLOBAL_POWER_TRAJECTORY_LENGTH,GLOBAL_POWER_WINDOW_LENGTH,GLOBAL_POWER_WINDOW_UPDATE_EVERY,GLOBAL_POWER_maxra_,GLOBAL_POWER_minra_,/* GLOBAL_POWER_CYCLE_BOTHER */1,/* ptree_bother */0,GLOBAL_PTREE_NLEGS,GLOBAL_PTREE_LEGTIME);}
	    else /* if (GLOBAL_CDRA!=NULL) */{ 
	      cd = GLOBAL_CDRA;
	      if (cd->ptree_bother){
		if (verbose){ printf(" %% %% updating odh->hra[%d]... ",oldinput);}
		for (nr1=0;nr1<cd->gli->nclusters;nr1++){
		  pnode_obs2dist_starter(oldinput,ninputs,0,NULL,cd->pra[nr1]->postree,-1,0,1);
		  if (verbose){ printf("clearing ptree... ");}
		  pnodeclear_starter(NULL,cd->pra[nr1]->pretree); 
		  pnodeclear_starter(NULL,cd->pra[nr1]->postree); 
		  cd->pra[nr1]->total_time=0;}
		if (verbose){ printf("\n");}}
	      if (cd->ptree_bother){
		if (verbose){ printf(" %% %% dumping clusterdatara->pra\n");}
		sprintf(tempchar,"input%d",oldinput); 
		clusterdataradump(cd,tempchar,2);}
	      if (nowinstance==0 && nowtime>0){
		if (verbose){ printf(" %% %% loading up cd->vtemp*\n");}
		if (cd->vtemp1==NULL){
		  if (verbose){ printf(" %% %% %% creating cd->vtemp1[size %d]\n",cd->gli->nclusters);}
		  temprara = (double **) tcalloc(cd->gli->nclusters,sizeof(double *));
		  for (nc=0;nc<cd->gli->nclusters;nc++){
		    if (verbose){ printf(" %% %% %% creating cd->vtemp1[size %d][size %d]\n",nc,cd->gli->cra[nc]->LN->length*cd->power_length/cd->trajectory_window*ninputs);}
		    temprara[nc] = (double *) tcalloc(cd->gli->cra[nc]->LN->length*cd->power_length/cd->trajectory_window*ninputs,sizeof(double));}
		  cd->vtemp1=temprara;}
		if (cd->vtemp1!=NULL){
		  temprara = (double **) cd->vtemp1;
		  for (nc=0;nc<cd->gli->nclusters;nc++){ 
		    tab = cd->gli->cra[nc]->LN->length*cd->power_length/cd->trajectory_window;
		    tab2 = 0+oldconcentration*tab+oldodor*SUITE_NCONCENTRATIONS*tab+oldbicuculline*SUITE_NCONCENTRATIONS*SUITE_NODORS*tab;
		    if (verbose){ printf(" %% %% %% manipulating cd->vtemp1[%d][%d]\n",nc,tab2);}
		    stra2ra(cd->strarara[nc][cd->indexing_nvar_refile[VARNAME_REGISTRY_spike_flag]],cd->gli->cra[nc]->LN->length,cd->trajectory_window,1,&(temprara[nc][tab2]));}}
		if (cd->vtemp2==NULL){
		  if (verbose){ printf(" %% %% %% creating cd->vtemp2[size %d]\n",cd->gli->nclusters);}
		  temprara = (double **) tcalloc(cd->gli->nclusters,sizeof(double *));
		  for (nc=0;nc<cd->gli->nclusters;nc++){
		    if (verbose){ printf(" %% %% %% creating cd->vtemp2[size %d][size %d]\n",cd->gli->nclusters,cd->stra_lfp[nc]->lpower_length*ninputs);}
		    temprara[nc] = (double *) tcalloc(cd->stra_lfp[nc]->lpower_length*ninputs,sizeof(double));}
		  cd->vtemp2=temprara;}
		if (cd->vtemp2!=NULL){
		  temprara = (double **) cd->vtemp2;
		  for (nc=0;nc<cd->gli->nclusters;nc++){
		    tab = cd->stra_lfp[nc]->lpower_length;
		    tab2 = 0+oldconcentration*tab+oldodor*SUITE_NCONCENTRATIONS*tab+oldbicuculline*SUITE_NCONCENTRATIONS*SUITE_NODORS*tab;
		    if (verbose){ printf(" %% %% %% manipulating cd->vtemp2[%d][%d]\n",nc,tab2);}
		    raplugin(&(temprara[nc][tab2]),1,tab,cd->stra_lfp[nc]->lpowerdata,1,tab,0,0);}}
		if (verbose){ printf(" %% %% resetting clusterdatara\n");}
		clusterdatarareset(GLOBAL_CDRA);}}}}
	INPUT_PULSE = ((GLOBAL_time-timera[nowtime])>=pulse_waitfor && (GLOBAL_time-timera[nowtime])<=(pulse_waitfor+pulse_howlong));
	if (verbose>2){ printf(" %% %% stepping input by adjusting input_pulse=%f\n",INPUT_PULSE);}}
      if (nowtime>=ntimes-1){
	if (verbose){ printf(" %% time %f, entering %d_th phase of suite one, finishing\n",GLOBAL_time,nowtime);}
	cd = GLOBAL_CDRA;
	if (cd->vtemp1!=NULL){
	  temprara = (double **) cd->vtemp1;
	  for (nc=0;nc<cd->gli->nclusters;nc++){ 
	    tab = cd->gli->cra[nc]->LN->length*cd->power_length/cd->trajectory_window;
	    tab3 = cd->power_length/cd->trajectory_window;
	    tempra = (double *) tcalloc(tab*ninputs,sizeof(double));
	    for (nowbicuculline=0;nowbicuculline<SUITE_NBICUCULLINES;nowbicuculline++){ for (nowodor=0;nowodor<SUITE_NODORS;nowodor++){ for (nowconcentration=0;nowconcentration<SUITE_NCONCENTRATIONS;nowconcentration++){
	      sprintf(tempchar,"./strara_concentration%d_odor%d_bcc%d_cluster%d_%srecord.pnm",nowconcentration,nowodor,nowbicuculline,nc,gs2);
	      tab2 = 0+nowconcentration*tab+nowodor*SUITE_NCONCENTRATIONS*tab+nowbicuculline*SUITE_NCONCENTRATIONS*SUITE_NODORS*tab;
	      WritePNMfile_color(&(temprara[nc][tab2]),cd->power_length/cd->trajectory_window,cd->gli->cra[nc]->LN->length,0,0,tempchar,7);
	      for (nr=0;nr<cd->gli->cra[nc]->LN->length;nr++){ raplugin(&(tempra[0+tab2/tab*tab3+nr*tab3*ninputs]),1,tab3,&(temprara[nc][0+nr*tab3+tab2]),1,tab3,0,0);}}}}
	    tempra2 = ra2pca(tempra,tab3*ninputs,cd->gli->cra[nc]->LN->length,NULL);
	    for (nowbicuculline=0;nowbicuculline<SUITE_NBICUCULLINES;nowbicuculline++){ for (nowodor=0;nowodor<SUITE_NODORS;nowodor++){ for (nowconcentration=0;nowconcentration<SUITE_NCONCENTRATIONS;nowconcentration++){
	      sprintf(tempchar,"./strara_pca_concentration%d_odor%d_bcc%d_cluster%d_%srecord",nowconcentration,nowodor,nowbicuculline,nc,gs2);
	      tab2 = 0+nowconcentration*tab3+nowodor*SUITE_NCONCENTRATIONS*tab3+nowbicuculline*SUITE_NCONCENTRATIONS*SUITE_NODORS*tab3;
	      pca2jpg(&(tempra2[tab2]),tab3,tempchar);}}}
	    tfree(tempra2);tempra2=NULL;
	    tfree(tempra);tempra=NULL;
	    tfree(temprara[nc]);temprara[nc]=NULL;}
	  tfree(temprara);temprara=NULL;cd->vtemp1=NULL;}
	if (cd->vtemp2!=NULL){
	  temprara = (double **) cd->vtemp2;
	  for (nc=0;nc<cd->gli->nclusters;nc++){ for (nowbicuculline=0;nowbicuculline<SUITE_NBICUCULLINES;nowbicuculline++){ for (nowodor=0;nowodor<SUITE_NODORS;nowodor++){ for (nowconcentration=0;nowconcentration<SUITE_NCONCENTRATIONS;nowconcentration++){
	    sprintf(tempchar,"./stra_lfp_concentration%d_odor%d_bcc%d_cluster%d_%srecord.pnm",nowconcentration,nowodor,nowbicuculline,nc,gs2);
	    tab = cd->stra_lfp[nc]->lpower_length;
	    tab2 = 0+nowconcentration*tab+nowodor*SUITE_NCONCENTRATIONS*tab+nowbicuculline*SUITE_NCONCENTRATIONS*SUITE_NODORS*tab;
	    tab3 = tab/cd->stra_lfp[nc]->lpower_window_length;
	    WritePNMfile_color(&(temprara[nc][tab2]),cd->stra_lfp[nc]->lpower_window_length,tab3,0,0,tempchar,7);
	    tempra = (double *) tcalloc(tab3,sizeof(double));
	    nr1=maximum(1,floor(cd->stra_lfp[nc]->lpower_window_length*cd->stra_lfp[nc]->update_timestep/cd->lfp_gammalo));
	    nr2=maximum(nr1+1,ceil(cd->stra_lfp[nc]->lpower_window_length*cd->stra_lfp[nc]->update_timestep/cd->lfp_gammahi));
	    if (nr2<cd->stra_lfp[nc]->lpower_window_length){ for (nr=0;nr<tab3;nr++){
	      stats("double",&(temprara[nc][nr1+nr*cd->stra_lfp[nc]->lpower_window_length+tab2]),nr2-nr1,NULL,NULL,&(tempra[nr]),NULL);}}
	    sprintf(tempchar,"./stra_lfpgamma_concentration%d_odor%d_bcc%d_cluster%d_%srecord",nowconcentration,nowodor,nowbicuculline,nc,gs2);
	    ra2jpg(tempra,"double",1,tab3,0,tempchar,0);
	    tfree(tempra);tempra=NULL;}}}
	  tfree(temprara[nc]);temprara[nc]=NULL;}
	  tfree(temprara);temprara=NULL;cd->vtemp2=NULL;}
/* 	if (cd->ptree_bother){ */
/* 	  if (verbose){ printf(" %% %% discriminability for cd->pra\n");} */
/* 	  for (nr1=0;nr1<cd->gli->nclusters;nr1++){ */
/* 	    if (verbose){ printf(" %% %% %% cluster %d of %d\n",nr1,cd->gli->nclusters);} */
/* 	    p = cd->pra[nr1]; */
/* 	    if (verbose){ printf(" %% %% %% calling ptree_obs2dist_starter\n");} */
/* 	    pnode_obs2dist_starter(0,0,SUITE_NINSTANCES,NULL,p->postree,-1,0,3); */
/* 	    if (verbose){ printf(" %% %% %% calling ptree_dumptemp_starter\n");} */
/* 	    ptree_dumptemp_starter(p,NULL,2,&obsdisthist2file,1); */
/* 	    if (verbose){ printf(" %% %% %% calling ptree_obsdisthist_tofig_starter\n");} */
/* 	    pnode_obsdisthist_tofig_starter(p,NULL,ninputs); */
/* 	    if (SUITE_BITBYBIT_RECORD){ */
/* 	      if (verbose){ printf(" %% %% %% making fnamebase\n");} */
/* 	      fnamebase = (char **) tcalloc(ninputs,sizeof(char *)); */
/* 	      for (nr2=0;nr2<ninputs;nr2++){  */
/* 		fnamebase[nr2] = (char *) tcalloc(256,sizeof(char));  */
/* 		sprintf(fnamebase[nr2],"ptree_cluster%dof%d_odor%d_%srecord",nr1,cd->gli->nclusters,nr2,gs2);} */
/* 	      if (verbose){ printf(" %% %% %% calling ptree_classify_starter\n");} */
/* 	      ptree_classify_starter(ninputs,GLOBAL_PTREE_NLEGS+1,2,fnamebase,NULL,SUITE_BITBYBIT_REMOVE); */
/* 	      if (verbose){ printf(" %% freeing fnamebase\n");} */
/* 	      for (nr2=0;nr2<ninputs;nr2++){ tfree(fnamebase[nr2]);fnamebase[nr2]=NULL;} */
	/* 	      tfree(fnamebase);fnamebase=NULL;}}} */
	if (SUITE_PTREE_BOTHER){
	  p = GLOBAL_PTREE;
	  if (verbose){ printf(" %% %% %% calling ptree_obs2dist_starter\n");}
	  pnode_obs2dist_starter(0,0,SUITE_NINSTANCES,NULL,p->postree,-1,0,3);
	  if (verbose){ printf(" %% %% %% calling ptree_dumptemp_starter\n");}
	  ptree_dumptemp_starter(p,NULL,2,&obsdisthist2file,1);
	  if (verbose){ printf(" %% %% %% calling ptree_obsdisthist_tofig_starter\n");}
	  pnode_obsdisthist_tofig_starter(p,NULL,ninputs);
	  if (SUITE_BITBYBIT_RECORD){
	    if (verbose){ printf(" %% %% %% making fnamebase\n");}
	    fnamebase = (char **) tcalloc(ninputs,sizeof(char *));
	    for (nr2=0;nr2<ninputs;nr2++){ 
	      fnamebase[nr2] = (char *) tcalloc(256,sizeof(char)); 
	      sprintf(fnamebase[nr2],"ptree_input%d_%srecord",nr2,gs2);}
	    if (verbose){ printf(" %% %% %% calling ptree_classify_starter\n");}
	    ptree_classify_starter(ninputs,GLOBAL_PTREE_NLEGS+1,2,fnamebase,NULL,SUITE_BITBYBIT_REMOVE);
	    if (verbose){ printf(" %% freeing fnamebase\n");}
	    for (nr2=0;nr2<ninputs;nr2++){ tfree(fnamebase[nr2]);fnamebase[nr2]=NULL;}
	    tfree(fnamebase);fnamebase=NULL;}}
	SUITE_ONSET_TIME = timera[nowtime]; OUTPUT_DUMP_EVERY = 0; /* don't bother dumping data */
	GLOBAL_TF = SUITE_ONSET_TIME+1;}}
    tfree(timera);timera=NULL;
    break;
  case 2: 
    ntimes = 2;
    timera = (double *) tcalloc(ntimes,sizeof(double));
    for (nt=0;nt<ntimes;nt++){ timera[nt] = nt*1024*SUITE_NSECONDS;}
    OUTPUT_DUMP_EVERY = SUITE_DUMPEVERY;
    pulse_waitfor = 128;
    pulse_howlong = 256;
    if (t==0){ /* clear everything */
      if (verbose){ printf(" %% \n");}
      if (verbose){ printf(" %% initializing suite two\n");}
      granule_findint(GLOBAL_ODOR_BASE,GLOBAL_ODOR_BASE,&nr1,&nr2);
      granule_int2input_or_input2int(1,GLOBAL_ODORra,GLOBAL_ODOR_BASE,&nr1,GLOBAL_ODOR_BASE,&nr2);
      GLOBAL_TF=0;
      SUITE_ONSET_TIME=-1;
      if (POWER_BOTHER!=0){ powertfree(GLOBAL_POWER); GLOBAL_POWER=NULL; POWER_BOTHER=0;}
      if (PTREE_BOTHER!=0){ ptreetfree(GLOBAL_PTREE); GLOBAL_PTREE=NULL; GLOBAL_PTREE_BITBYBIT=0; GLOBAL_PTREE_BITBYBIT_CURRENT_NUMBER=0; PTREE_BOTHER=0; if (GLOBAL_INPUT_PTREE!=NULL){ ptreetfree(GLOBAL_INPUT_PTREE); GLOBAL_INPUT_PTREE=NULL;}}
      if (CLUSTERDATA_BOTHER!=0){ clusterdatatfree(GLOBAL_CDRA); GLOBAL_CDRA=NULL; CLUSTERDATA_BOTHER=0;}}
    else /* if (t>0) */{
      nowtime = 0; for (nt=0;nt<ntimes-1;nt++){ if (t>=timera[nt] && t<timera[nt+1]){ nowtime=nt;}} 
      if (t>=timera[ntimes-1]){ nowtime=ntimes-1;}
      if (nowtime<ntimes-1){
	if (SUITE_ONSET_TIME != timera[nowtime]){
	  if (verbose){ printf(" %% time %f, entering %d_th phase of suite two... ",GLOBAL_time,nowtime);}
	  if (SUITE_PTREE_BOTHER){
	    if (GLOBAL_PTREE==NULL){
	      PTREE_BOTHER=1; GLOBAL_PTREE_REGION_TYPE=1;
	      GLOBAL_PTREE = ptreemake(GLOBAL_PTREE_NREGIONS,GLOBAL_PTREE_EVENT_WITHIN,GLOBAL_PTREE_EVENT_THRESHOLD,GLOBAL_PTREE_REGION_TYPE,GLOBAL_PTREE_NLEGS,GLOBAL_PTREE_LEGTIME);}
	    if (GLOBAL_INPUT_PTREE==NULL){
	      GLOBAL_INPUT_PTREE = ptreemake(GLOBAL_PTREE_NREGIONS,GLOBAL_PTREE_EVENT_WITHIN,GLOBAL_PTREE_EVENT_THRESHOLD,2,0,1);}}
	  SUITE_ONSET_TIME = timera[nowtime];}
	if ((OUTPUT_DUMP_EVERY>0 && (int)floor(GLOBAL_time/OUTPUT_DUMP_EVERY)<(int)floor((GLOBAL_time+GLOBAL_DT)/OUTPUT_DUMP_EVERY)) || (GLOBAL_TF > 0 && GLOBAL_time > GLOBAL_TF)){
	  if (PTREE_BOTHER){
	    if (verbose){ printf(" %% %% updating ptree at time %f...",GLOBAL_time);}
	    pnode_obs2dist_starter(0,1,0,NULL,GLOBAL_PTREE->postree,-1,0,1);
	    if (verbose){ printf("clearing ptree... ");}
	    pnodeclear_starter(NULL,GLOBAL_PTREE->pretree); pnodeclear_starter(NULL,GLOBAL_PTREE->postree); GLOBAL_PTREE->total_time=0;
	    if (verbose){ printf("\n");}
	    if (verbose){ printf(" %% %% updating input_ptree at time %f...",GLOBAL_time);}
	    pnode_obs2dist_starter(0,1,0,NULL,GLOBAL_INPUT_PTREE->postree,-1,0,1);
	    if (verbose){ printf("clearing input_ptree... ");}
	    pnodeclear_starter(NULL,GLOBAL_INPUT_PTREE->pretree); pnodeclear_starter(NULL,GLOBAL_INPUT_PTREE->postree); 
	    GLOBAL_INPUT_PTREE->total_time=0; if (verbose){ printf("\n");}}}
	pulse_temp = (GLOBAL_time - OUTPUT_DUMP_EVERY*(int)floor(GLOBAL_time/OUTPUT_DUMP_EVERY));
	INPUT_PULSE = (pulse_temp>pulse_waitfor && pulse_temp<(pulse_waitfor+pulse_howlong));
	if (verbose>2){ printf(" %% %% stepping input by adjusting input_pulse=%f\n",INPUT_PULSE);}}
      else /* if (nowtime==ntimes-1) */{
	if (SUITE_ONSET_TIME != timera[nowtime]){
	  if (verbose){ printf(" %% time %f, entering %d_th phase of suite two -- finishing\n",GLOBAL_time,nowtime);}
	  SUITE_ONSET_TIME=timera[nowtime];
	  GLOBAL_TF = timera[nowtime]+1;
	  if (PTREE_BOTHER){ 
	    if (verbose){ printf(" %% %% dumping ptree at time %f...",GLOBAL_time);}
	    granule_int2input_or_input2int(0,GLOBAL_ODORra,GLOBAL_ODOR_BASE,&nr1,GLOBAL_ODOR_BASE,&nr2);
	    sprintf(tempchar,"i%dj%d",nr1,nr2);
	    ptree_dumptemp_starter(GLOBAL_PTREE,tempchar,2,&obsdisthist2file,1);
	    sprintf(tempchar,"input_i%dj%d",nr1,nr2);
	    ptree_dumptemp_starter(GLOBAL_INPUT_PTREE,tempchar,2,&obsdisthist2file,1);} 
	  if (verbose){ printf(" %% %% calling granule_plot\n");}
	  granule_plot(GLOBAL_ODOR_BASE,GLOBAL_ODOR_BASE);
	  if (verbose){ printf(" %% %% lazy, exiting\n");}
	  exit(EXIT_SUCCESS);}}}
    tfree(timera);timera=NULL;
    break;
  case 3: 
    ntimes = 2;
    timera = (double *) tcalloc(ntimes,sizeof(double));
    for (nt=0;nt<ntimes;nt++){ timera[nt] = nt*1024*SUITE_NSECONDS;}
    OUTPUT_DUMP_EVERY = SUITE_DUMPEVERY;
    pulse_waitfor = 128;
    pulse_howlong = 256;
    npoints=16*16;
    if (t==0){ /* clear everything */
      if (verbose){ printf(" %% \n");}
      if (verbose){ printf(" %% initializing suite three...");}
      GLOBAL_Nra->reweight = ptreemake(GLOBAL_PTREE_NREGIONS,1,1,0,GLOBAL_PTREE_NLEGS,GLOBAL_PTREE_LEGTIME);
      synaptic_sphere_findint2input_or_input2int(1,GLOBAL_Nra->reweight,npoints,&nr1);
      if (verbose>2){ printf(" input %d, reweight:\n",nr1);pnodeprintf(NULL,GLOBAL_Nra->reweight->postree,-1,0);}
      GLOBAL_TF=0;
      SUITE_ONSET_TIME=-1;
      if (POWER_BOTHER!=0){ powertfree(GLOBAL_POWER); GLOBAL_POWER=NULL; POWER_BOTHER=0;}
      if (PTREE_BOTHER!=0){ ptreetfree(GLOBAL_PTREE); GLOBAL_PTREE=NULL; GLOBAL_PTREE_BITBYBIT=0; GLOBAL_PTREE_BITBYBIT_CURRENT_NUMBER=0; PTREE_BOTHER=0; if (GLOBAL_INPUT_PTREE!=NULL){ ptreetfree(GLOBAL_INPUT_PTREE); GLOBAL_INPUT_PTREE=NULL;}}
      if (CLUSTERDATA_BOTHER!=0){ clusterdatatfree(GLOBAL_CDRA); GLOBAL_CDRA=NULL; CLUSTERDATA_BOTHER=0;}}
    else /* if (t>0) */{
      nowtime = 0; for (nt=0;nt<ntimes-1;nt++){ if (t>=timera[nt] && t<timera[nt+1]){ nowtime=nt;}} 
      if (t>=timera[ntimes-1]){ nowtime=ntimes-1;}
      if (nowtime<ntimes-1){
	if (SUITE_ONSET_TIME != timera[nowtime]){
	  if (verbose){ printf(" %% time %f, entering %d_th phase of suite three... ",GLOBAL_time,nowtime);}
	  if (SUITE_PTREE_BOTHER){
	    if (GLOBAL_PTREE==NULL){
	      PTREE_BOTHER=1; GLOBAL_PTREE_REGION_TYPE=1;
	      GLOBAL_PTREE = ptreemake(GLOBAL_PTREE_NREGIONS,GLOBAL_PTREE_EVENT_WITHIN,GLOBAL_PTREE_EVENT_THRESHOLD,GLOBAL_PTREE_REGION_TYPE,GLOBAL_PTREE_NLEGS,GLOBAL_PTREE_LEGTIME);}
	    if (GLOBAL_INPUT_PTREE==NULL){
	      GLOBAL_INPUT_PTREE = ptreemake(GLOBAL_PTREE_NREGIONS,GLOBAL_PTREE_EVENT_WITHIN,GLOBAL_PTREE_EVENT_THRESHOLD,2,0,1);}}
	  SUITE_ONSET_TIME = timera[nowtime];}
	if ((OUTPUT_DUMP_EVERY>0 && (int)floor(GLOBAL_time/OUTPUT_DUMP_EVERY)<(int)floor((GLOBAL_time+GLOBAL_DT)/OUTPUT_DUMP_EVERY)) || (GLOBAL_TF > 0 && GLOBAL_time > GLOBAL_TF)){
	  if (PTREE_BOTHER){
	    if (verbose){ printf(" %% %% updating ptree at time %f...",GLOBAL_time);}
	    pnode_obs2dist_starter(0,1,0,NULL,GLOBAL_PTREE->postree,-1,0,1);
	    if (verbose){ printf("clearing ptree... ");}
	    pnodeclear_starter(NULL,GLOBAL_PTREE->pretree); pnodeclear_starter(NULL,GLOBAL_PTREE->postree); GLOBAL_PTREE->total_time=0;
	    if (verbose){ printf("\n");}
	    if (verbose){ printf(" %% %% updating input_ptree at time %f...",GLOBAL_time);}
	    pnode_obs2dist_starter(0,1,0,NULL,GLOBAL_INPUT_PTREE->postree,-1,0,1);
	    if (verbose){ printf("clearing input_ptree... ");}
	    pnodeclear_starter(NULL,GLOBAL_INPUT_PTREE->pretree); pnodeclear_starter(NULL,GLOBAL_INPUT_PTREE->postree); 
	    GLOBAL_INPUT_PTREE->total_time=0; if (verbose){ printf("\n");}}}
	pulse_temp = (GLOBAL_time - OUTPUT_DUMP_EVERY*(int)floor(GLOBAL_time/OUTPUT_DUMP_EVERY));
	INPUT_PULSE = (pulse_temp>pulse_waitfor && pulse_temp<(pulse_waitfor+pulse_howlong));
	if (verbose>2){ printf(" %% %% stepping input by adjusting input_pulse=%f\n",INPUT_PULSE);}}
      else /* if (nowtime==ntimes-1) */{
	if (SUITE_ONSET_TIME != timera[nowtime]){
	  if (verbose){ printf(" %% time %f, entering %d_th phase of suite three -- finishing\n",GLOBAL_time,nowtime);}
	  SUITE_ONSET_TIME=timera[nowtime];
	  GLOBAL_TF = timera[nowtime]+1;
	  if (PTREE_BOTHER){ 
	    if (verbose){ printf(" %% %% dumping ptree at time %f...",GLOBAL_time);}
	    synaptic_sphere_findint2input_or_input2int(0,GLOBAL_Nra->reweight,npoints,&nr1);
	    sprintf(tempchar,"s%d",nr1);
	    ptree_dumptemp_starter(GLOBAL_PTREE,tempchar,2,&obsdisthist2file,1);
	    sprintf(tempchar,"input_s%d",nr1);
	    ptree_dumptemp_starter(GLOBAL_INPUT_PTREE,tempchar,2,&obsdisthist2file,1);} 
	  if (verbose){ printf(" %% %% calling synaptic_sphere_prediction_plot\n");}
	  synaptic_sphere_prediction_plot(npoints);
	  if (verbose){ printf(" %% %% lazy, exiting\n");}
	  exit(EXIT_SUCCESS);}}}
    tfree(timera);timera=NULL;
    break;    
  case 4: /* run a background run, estimate bi-poo direction, then test granulation and compression */
    ntimes = SUITE_NODORS*2 + 1;
    timera = (double *) tcalloc(ntimes,sizeof(double));
    for (nt=0;nt<ntimes;nt++){ timera[nt] = nt*1024*SUITE_NSECONDS;}
    INPUT_PULSE = 1;
    if (t==0){ /* clear everything */
      if (verbose){ printf(" %% \n");}
      if (verbose){ printf(" %% initializing suite four\n");}
      if (verbose){ printf(" times are: ");raprintf(timera,"double",1,ntimes,"timera ");}
      GLOBAL_TF=0;
      SUITE_ONSET_TIME=-1;
      if (POWER_BOTHER!=0){ powertfree(GLOBAL_POWER); GLOBAL_POWER=NULL; POWER_BOTHER=0;}
      if (PTREE_BOTHER!=0){ ptreetfree(GLOBAL_PTREE); GLOBAL_PTREE=NULL; GLOBAL_PTREE_BITBYBIT=0; GLOBAL_PTREE_BITBYBIT_CURRENT_NUMBER=0; PTREE_BOTHER=0;}
      if (CLUSTERDATA_BOTHER!=0){ clusterdatatfree(GLOBAL_CDRA); GLOBAL_CDRA=NULL; CLUSTERDATA_BOTHER=0;}}
    else /* if (t>0) */{
      nowtime = 0; for (nt=0;nt<ntimes-1;nt++){ if (t>=timera[nt] && t<timera[nt+1]){ nowtime=nt;}} 
      if (t>=timera[ntimes-1]){ nowtime=ntimes-1;}
      nowodor = (nowtime>1? (nowtime-2)/2 : -1); oldodor = ((nowtime-1)>1? (nowtime-1-2)/2 : -1);
      nowreweight = nowtime%2; oldreweight = (nowtime-1)%2;
      if (nowtime<ntimes-1){
	if (fabs(SUITE_ONSET_TIME-timera[nowtime])>0.1){
	  if (verbose){ printf(" %% time %f, entering %d_th phase of suite four\n",GLOBAL_time,nowtime);}
	  SUITE_ONSET_TIME = timera[nowtime]; OUTPUT_DUMP_EVERY = 0; /* don't bother dumping data */
	  if (nowodor>=0){ 
	    if (nowreweight==0){
	      odortfree(GLOBAL_ODORra); GLOBAL_ODORra=NULL;
	      GLOBAL_ODORra = odormake(0,NULL,0,NULL,GLOBAL_ODORra_START);
	      o1 = SUITE_ODORra_BACON; 
	      for (nr1=0;nr1<o1->ntypes;nr1++){ for (nr2=0;nr2<o1->lengthra[nr1];nr2++){ o1->rara[nr1][nr2] = 2*rand01-1;}}
	      odorplusequals(GLOBAL_ODORra,1,o1);}
	    else /* if (nowreweight==1) */{ /* do nothing */ }
	    sprintf(tempchar,"odor_%srecord_odor%d_reweight%d",GLOBAL_STRING_2,nowodor,nowreweight);
	    odor_full_fwrite(tempchar,GLOBAL_ODORra);}
	  else if (nowodor<0){ 
	    sprintf(tempchar,"odor_%srecord_basic_reweight%d",GLOBAL_STRING_2,nowreweight);
	    odor_full_fwrite(tempchar,GLOBAL_ODORra);}
	  if (verbose){ printf("%% %% GLOBAL_ODORra: ");odorfprintf_full(stdout,GLOBAL_ODORra);}
	  if (nowreweight){ GLOBAL_REWEIGHT_STRENGTH=SUITE_4_REWEIGHT_STRENGTH;} else if (!nowreweight){ GLOBAL_REWEIGHT_STRENGTH=0;}
	  if (verbose){ printf(" %% %% GLOBAL_REWEIGHT_STRENGTH %f\n",GLOBAL_REWEIGHT_STRENGTH);}
	  if (verbose){ printf(" %% %% setting INPUT_CONTRAST=%0.2f\n",INPUT_CONTRAST);}
	  if (SUITE_PTREE_BOTHER){
	    if (GLOBAL_PTREE==NULL){
	      if (verbose){ printf(" %% %% making ptree\n");}
	      PTREE_BOTHER=1; GLOBAL_PTREE_BITBYBIT=0; GLOBAL_PTREE_BITBYBIT_CURRENT_NUMBER=0; GLOBAL_PTREE_REGION_TYPE=1;
	      GLOBAL_PTREE = ptreemake(GLOBAL_PTREE_NREGIONS,GLOBAL_PTREE_EVENT_WITHIN,GLOBAL_PTREE_EVENT_THRESHOLD,GLOBAL_PTREE_REGION_TYPE,/* GLOBAL_PTREE_NLEGS */2,GLOBAL_PTREE_LEGTIME);}
	    else /* if (GLOBAL_PTREE!=NULL) */{
	      if (verbose){ printf("dumping ptree... ");}
	      ptreerate(GLOBAL_PTREE); 
	      if (oldodor>=0){ sprintf(tempchar,"./ptree_%srecord_odor%d_reweight%d",GLOBAL_STRING_2,oldodor,oldreweight); }
	      else if (oldodor<0){ sprintf(tempchar,"./ptree_%srecord_basic_reweight%d",GLOBAL_STRING_2,oldreweight); }
	      ptreedump_starter(GLOBAL_PTREE,tempchar,2,1,0,0,+1,-1,0);
	      if (verbose){ printf("clearing ptree... ");}
	      pnodeclear_starter(NULL,GLOBAL_PTREE->pretree); pnodeclear_starter(NULL,GLOBAL_PTREE->postree); GLOBAL_PTREE->total_time=0;
	      if (verbose){ printf("\n");}}
	    if (oldodor<0 && !oldreweight){
	      if (verbose){ printf(" %% obtaining reweight from basic ptree\n");}
	      sprintf(tempchar,"./ptree_%srecord_basic_reweight0",GLOBAL_STRING_2);
	      GLOBAL_Nra->reweight = ptree_obtain_bipoo_reweight(tempchar);
	      sprintf(tempchar,"./ptree_%srecord_reweight",GLOBAL_STRING_2);
	      ptreedump_starter(GLOBAL_Nra->reweight,tempchar,2,1,0,0,0,0,1);}}}}
      if (nowtime>=ntimes-1){
	if (fabs(SUITE_ONSET_TIME-timera[nowtime])>0.1){
	  if (verbose){ printf(" %% time %f, entering %d_th phase of suite four, finishing\n",GLOBAL_time,nowtime);}
	  SUITE_ONSET_TIME = timera[nowtime]; OUTPUT_DUMP_EVERY = 0; /* don't bother dumping data */
	  if (SUITE_PTREE_BOTHER){
	    if (GLOBAL_PTREE!=NULL) {
	      if (verbose){ printf("dumping ptree... ");}
	      ptreerate(GLOBAL_PTREE); 
	      if (oldodor>=0){ sprintf(tempchar,"./ptree_%srecord_odor%d_reweight%d",GLOBAL_STRING_2,oldodor,oldreweight); }
	      else if (oldodor<0){ sprintf(tempchar,"./ptree_%srecord_basic_reweight%d",GLOBAL_STRING_2,oldreweight); }
	      ptreedump_starter(GLOBAL_PTREE,tempchar,2,1,0,0,+1,-1,0);
	      if (verbose){ printf("freeing ptree... ");}
	      ptreetfree(GLOBAL_PTREE); GLOBAL_PTREE=NULL; PTREE_BOTHER=0;}}
	  GLOBAL_TF = SUITE_ONSET_TIME+1;
	  ptree_test_mcpit();}}}
    tfree(timera);timera=NULL;
    break;
  case 5: /* run a simple odor */
    ntimes = SUITE_NODORS + 1;
    timera = (double *) tcalloc(ntimes,sizeof(double));
    for (nt=0;nt<ntimes;nt++){ timera[nt] = nt*1024*SUITE_NSECONDS;}
    pulse_waitfor = 1024;
    pulse_howlong = 2048;
    INPUT_PULSE = 1;
    if (t==0){ /* clear everything */
      if (verbose){ printf(" %% \n");}
      if (verbose){ printf(" %% initializing suite five\n");}
      if (verbose){ printf(" times are: ");raprintf(timera,"double",1,ntimes,"timera ");}
      GLOBAL_TF=0;
      SUITE_ONSET_TIME=-1;
      if (POWER_BOTHER!=0){ powertfree(GLOBAL_POWER); GLOBAL_POWER=NULL; POWER_BOTHER=0;}
      if (PTREE_BOTHER!=0){ ptreetfree(GLOBAL_PTREE); GLOBAL_PTREE=NULL; GLOBAL_PTREE_BITBYBIT=0; GLOBAL_PTREE_BITBYBIT_CURRENT_NUMBER=0; PTREE_BOTHER=0;}
      if (CLUSTERDATA_BOTHER!=0){ clusterdatatfree(GLOBAL_CDRA); GLOBAL_CDRA=NULL; CLUSTERDATA_BOTHER=0;}}
    else /* if (t>0) */{
      nowtime = 0; for (nt=0;nt<ntimes-1;nt++){ if (t>=timera[nt] && t<timera[nt+1]){ nowtime=nt;}} 
      if (t>=timera[ntimes-1]){ nowtime=ntimes-1;}
      nowodor = nowtime;
      if (nowtime<ntimes-1){
	if (fabs(SUITE_ONSET_TIME-timera[nowtime])>0.1){
	  if (verbose){ printf(" %% time %f, entering %d_th phase of suite five\n",GLOBAL_time,nowtime);}
	  SUITE_ONSET_TIME = timera[nowtime]; OUTPUT_DUMP_EVERY = 0; /* don't bother dumping data */
	  if (nowodor>=0){ 
	    odortfree(GLOBAL_ODORra); GLOBAL_ODORra=NULL;
	    GLOBAL_ODORra = odormake(0,NULL,0,NULL,GLOBAL_ODORra_START);
	    o1 = SUITE_ODORra_BACON; 
	    for (nr1=0;nr1<o1->ntypes;nr1++){ for (nr2=0;nr2<o1->lengthra[nr1];nr2++){ o1->rara[nr1][nr2] = 2*rand01-1;}}
	    odorplusequals(GLOBAL_ODORra,1,o1);}
	  if (verbose){ printf("%% %% GLOBAL_ODORra: ");odorfprintf_full(stdout,GLOBAL_ODORra);}}
	pulse_temp = GLOBAL_time - timera[nowtime];
	INPUT_PULSE = (pulse_temp>pulse_waitfor && pulse_temp<(pulse_waitfor+pulse_howlong));}
      if (nowtime>=ntimes-1){
	if (fabs(SUITE_ONSET_TIME-timera[nowtime])>0.1){
	  if (verbose){ printf(" %% time %f, entering %d_th phase of suite five, finishing\n",GLOBAL_time,nowtime);}
	  SUITE_ONSET_TIME = timera[nowtime]; OUTPUT_DUMP_EVERY = 0; /* don't bother dumping data */
	  GLOBAL_TF = SUITE_ONSET_TIME+1;}}}
    tfree(timera);timera=NULL;
    break;
  case 6: /* run a simple odor for ntrials */
    ntimes = SUITE_NINSTANCES + 1;
    SUITE_NSECONDS = maximum(2,SUITE_NSECONDS);
    SUITE_POWER_BOTHER = 1;
    timera = (double *) tcalloc(ntimes,sizeof(double));
    for (nt=0;nt<ntimes;nt++){ timera[nt] = nt*1024*SUITE_NSECONDS;}
    pulse_waitfor = 1024;
    pulse_howlong = 512;
    INPUT_PULSE = 1;
    if (t==0){ /* clear everything */
      if (verbose){ printf(" %% \n");}
      if (verbose){ printf(" %% initializing suite six\n");}
      if (verbose){ printf(" times are: ");raprintf(timera,"double",1,ntimes,"timera ");}
      GLOBAL_TF=0;
      SUITE_ONSET_TIME=-1;
      if (POWER_BOTHER!=0){ powertfree(GLOBAL_POWER); GLOBAL_POWER=NULL; POWER_BOTHER=0;}
      if (PTREE_BOTHER!=0){ ptreetfree(GLOBAL_PTREE); GLOBAL_PTREE=NULL; GLOBAL_PTREE_BITBYBIT=0; GLOBAL_PTREE_BITBYBIT_CURRENT_NUMBER=0; PTREE_BOTHER=0;}
      if (CLUSTERDATA_BOTHER!=0){ clusterdatatfree(GLOBAL_CDRA); GLOBAL_CDRA=NULL; CLUSTERDATA_BOTHER=0;}
      POWER_BOTHER=1;
      GLOBAL_POWER = powermake(GLOBAL_Nra,SUITE_NINSTANCES*SUITE_NSECONDS*1024,GLOBAL_POWER_INDEXING_NTYPE_LENGTH,GLOBAL_POWER_INDEXING_NTYPE_CHECKOUT,GLOBAL_POWER_INDEXING_NTYPE_REFILE,GLOBAL_POWER_INDEXING_NVAR_LENGTH,GLOBAL_POWER_INDEXING_NVAR_CHECKOUT,GLOBAL_POWER_INDEXING_NVAR_REFILE,GLOBAL_POWER_maxra_,GLOBAL_POWER_minra_,0/* GLOBAL_POWER_CORRELATION_BOTHER */);}
    else /* if (t>0) */{
      nowtime = 0; for (nt=0;nt<ntimes-1;nt++){ if (t>=timera[nt] && t<timera[nt+1]){ nowtime=nt;}}
      if (t>=timera[ntimes-1]){ nowtime=ntimes-1;}
      if (nowtime<ntimes-1){
	if (fabs(SUITE_ONSET_TIME-timera[nowtime])>0.1){
	  if (verbose){ printf(" %% time %f, entering %d_th phase of suite six\n",GLOBAL_time,nowtime);}
	  SUITE_ONSET_TIME = timera[nowtime]; OUTPUT_DUMP_EVERY = 0;}
	pulse_temp = GLOBAL_time - timera[nowtime];
	INPUT_PULSE = (pulse_temp>pulse_waitfor && pulse_temp<(pulse_waitfor+pulse_howlong));}
      if (nowtime>=ntimes-1){
	if (fabs(SUITE_ONSET_TIME-timera[nowtime])>0.1){
	  if (verbose){ printf(" %% time %f, entering %d_th phase of suite six, finishing\n",GLOBAL_time,nowtime);}
	  SUITE_ONSET_TIME = timera[nowtime]; OUTPUT_DUMP_EVERY = 0; /* don't bother dumping data */
	  powerdump_trial_average(GLOBAL_POWER,NULL,SUITE_NINSTANCES);
	  powerdump_trial_interlace(GLOBAL_POWER,NULL,SUITE_NINSTANCES);
	  exit(EXIT_SUCCESS);
	  GLOBAL_TF = SUITE_ONSET_TIME+1;}}}
    tfree(timera);timera=NULL;
    break;
  case 7: 
    /* run multiple odors for multiple trials over various concentrations. 
       for use with wilson_type neurons (i.e., GLOBAL_NEURON_MODEL==6)
       ninstances trials, 
       nconcentration concentrations, 
       up to minimum(suite_nodors,global_nclusters*(global_nclusters-1)) combinatorially different odors, indexed by nowodor.
       minimum(5,nbicucullines) alternate architectures.
       nbicuculline==0: standard architecture
       nbicuculline==1: gabaA blocked
       nbicuculline==2: gabaB blocked 
       nbicuculline==3: gabaA and gabaB blocked
       nbicuculline==4: nAch blocked */
    SUITE_NSECONDS = maximum(2,SUITE_NSECONDS);
    SUITE_NODORS = minimum(GLOBAL_NCLUSTERS*maximum(1,GLOBAL_NCLUSTERS-1),SUITE_NODORS);
    SUITE_NBICUCULLINES = minimum(5,SUITE_NBICUCULLINES);
    SUITE_POWER_BOTHER = 1;
    ntimes = SUITE_NINSTANCES*SUITE_NCONCENTRATIONS*SUITE_NODORS*SUITE_NBICUCULLINES + 1;
    ninputs = SUITE_NCONCENTRATIONS*SUITE_NODORS*SUITE_NBICUCULLINES;
    timera = (double *) tcalloc(ntimes,sizeof(double));
    for (nt=0;nt<ntimes;nt++){ timera[nt] = nt*1024*SUITE_NSECONDS;}
    pulse_waitfor = 1024;
    pulse_howlong = 512;
    INPUT_PULSE = 1;
    if (t==0){ /* clear everything */
      if (verbose){ printf(" %% \n");}
      if (verbose){ printf(" %% initializing suite seven\n");}
      if (verbose){ printf(" times are: ");raprintf(timera,"double",1,ntimes,"timera ");}
      GLOBAL_TF=0;
      SUITE_ONSET_TIME=-1;
      if (POWER_BOTHER!=0){ powertfree(GLOBAL_POWER); GLOBAL_POWER=NULL; POWER_BOTHER=0;}
      if (PTREE_BOTHER!=0){ ptreetfree(GLOBAL_PTREE); GLOBAL_PTREE=NULL; GLOBAL_PTREE_BITBYBIT=0; GLOBAL_PTREE_BITBYBIT_CURRENT_NUMBER=0; PTREE_BOTHER=0;}
      if (CLUSTERDATA_BOTHER!=0){ clusterdatatfree(GLOBAL_CDRA); GLOBAL_CDRA=NULL; CLUSTERDATA_BOTHER=0;}
      POWER_BOTHER=1;
      GLOBAL_POWER = powermake(GLOBAL_Nra,SUITE_NSECONDS*1024,GLOBAL_POWER_INDEXING_NTYPE_LENGTH,GLOBAL_POWER_INDEXING_NTYPE_CHECKOUT,GLOBAL_POWER_INDEXING_NTYPE_REFILE,GLOBAL_POWER_INDEXING_NVAR_LENGTH,GLOBAL_POWER_INDEXING_NVAR_CHECKOUT,GLOBAL_POWER_INDEXING_NVAR_REFILE,GLOBAL_POWER_maxra_,GLOBAL_POWER_minra_,0/* GLOBAL_POWER_CORRELATION_BOTHER */);}
    else /* if (t>0) */{
      nowtime = 0; for (nt=0;nt<ntimes-1;nt++){ if (t>=timera[nt] && t<timera[nt+1]){ nowtime=nt;}} 
      if (t>=timera[ntimes-1]){ nowtime=ntimes-1;}
      nowbicuculline = nowtime/(SUITE_NINSTANCES*SUITE_NCONCENTRATIONS*SUITE_NODORS);
      nowodor = (nowtime-nowbicuculline*SUITE_NINSTANCES*SUITE_NCONCENTRATIONS*SUITE_NODORS)/(SUITE_NINSTANCES*SUITE_NCONCENTRATIONS);
      nowconcentration = (nowtime-nowbicuculline*SUITE_NINSTANCES*SUITE_NCONCENTRATIONS*SUITE_NODORS-nowodor*SUITE_NINSTANCES*SUITE_NCONCENTRATIONS)/(SUITE_NINSTANCES);
      nowinstance = (nowtime-nowbicuculline*SUITE_NINSTANCES*SUITE_NCONCENTRATIONS*SUITE_NODORS-nowodor*SUITE_NINSTANCES*SUITE_NCONCENTRATIONS-nowconcentration*SUITE_NINSTANCES)/(1);
      nowinput = nowconcentration + nowodor*SUITE_NCONCENTRATIONS + nowbicuculline*SUITE_NCONCENTRATIONS*SUITE_NODORS;
      oldtime = maximum(0,nowtime-1);
      oldbicuculline = oldtime/(SUITE_NINSTANCES*SUITE_NCONCENTRATIONS*SUITE_NODORS);
      oldodor = (oldtime-oldbicuculline*SUITE_NINSTANCES*SUITE_NCONCENTRATIONS*SUITE_NODORS)/(SUITE_NINSTANCES*SUITE_NCONCENTRATIONS);
      oldconcentration = (oldtime-oldbicuculline*SUITE_NINSTANCES*SUITE_NCONCENTRATIONS*SUITE_NODORS-oldodor*SUITE_NINSTANCES*SUITE_NCONCENTRATIONS)/(SUITE_NINSTANCES);
      oldinstance = (oldtime-oldbicuculline*SUITE_NINSTANCES*SUITE_NCONCENTRATIONS*SUITE_NODORS-oldodor*SUITE_NINSTANCES*SUITE_NCONCENTRATIONS-oldconcentration*SUITE_NINSTANCES)/(1);
      oldinput = oldconcentration + oldodor*SUITE_NCONCENTRATIONS + oldbicuculline*SUITE_NCONCENTRATIONS*SUITE_NODORS;
      if (nowtime<ntimes-1){
	if (fabs(SUITE_ONSET_TIME-timera[nowtime])>0.1){
	  if (verbose){ printf(" %% time %f, entering %d_th phase of suite seven, bicuculline %d odor %d concentration %d instance %d\n",GLOBAL_time,nowtime,nowbicuculline,nowodor,nowconcentration,nowinstance);}
	  SUITE_ONSET_TIME = timera[nowtime]; OUTPUT_DUMP_EVERY = 0; /* don't bother dumping data */
	  if (nowtime>0){ 
	    if (verbose){ printf(" %% depositing old power\n");}
	    suite_7_power_deposit_helper(oldinstance,oldconcentration,oldodor,oldbicuculline,oldtime);
	    if (verbose){ printf(" %% remaking new power\n");}
	    powertfree(GLOBAL_POWER); GLOBAL_POWER=NULL;
	    GLOBAL_POWER = powermake(GLOBAL_Nra,SUITE_NSECONDS*1024,GLOBAL_POWER_INDEXING_NTYPE_LENGTH,GLOBAL_POWER_INDEXING_NTYPE_CHECKOUT,GLOBAL_POWER_INDEXING_NTYPE_REFILE,GLOBAL_POWER_INDEXING_NVAR_LENGTH,GLOBAL_POWER_INDEXING_NVAR_CHECKOUT,GLOBAL_POWER_INDEXING_NVAR_REFILE,GLOBAL_POWER_maxra_,GLOBAL_POWER_minra_,0/* GLOBAL_POWER_CORRELATION_BOTHER */);}
	  odortfree(GLOBAL_ODORra); GLOBAL_ODORra=NULL;
	  GLOBAL_ODORra = odormake(0,NULL,0,NULL,GLOBAL_ODORra_START);
	  o1 = GLOBAL_ODORra;
	  for (nr1=0;nr1<o1->ntypes;nr1++){ for (nr2=0;nr2<o1->lengthra[nr1];nr2++){ o1->rara[nr1][nr2] = 0;}}
	  nr1 = TYPENAME_REGISTRY_wilson_ORN;
	  tab2 = 0; tab3 = GLOBAL_NCLUSTERS*(GLOBAL_NCLUSTERS-1); periodify("int",&nowodor,&tab2,&tab3,&tab);
	  tab2 = tab/maximum(1,GLOBAL_NCLUSTERS-1); tab3 = tab%maximum(1,GLOBAL_NCLUSTERS-1); if (tab3>=tab2){ tab3 += 1;}
	  tab2 = periodize(tab2,0,GLOBAL_NCLUSTERS); tab3 = periodize(tab3,0,GLOBAL_NCLUSTERS);
	  if (verbose>1){ printf(" %% clusters %d,%d chosen by periodize(nowodor,0,GLOBAL_NCLUSTERS*(GLOBAL_NCLUSTERS-1))\n",tab2,tab3);}
	  for (nr2=0;nr2<o1->lengthra[nr1];nr2++){ 
	    nr3 = nr2%(o1->lengthra[nr1]/GLOBAL_NCLUSTERS);
	    nr4 = nr2/(o1->lengthra[nr1]/GLOBAL_NCLUSTERS);
	    if (nr4==tab2){ o1->rara[nr1][nr2] = o1->base*exp(-0*1.75/(double)GLOBAL_NCLUSTERS);}
	    else if (nr4==tab3){ o1->rara[nr1][nr2] = o1->base*exp(-1*1.75/(double)GLOBAL_NCLUSTERS);}
	    else /* if (nr4!=tab2 && nr4!=tab3) */{ o1->rara[nr1][nr2] = o1->base*exp(-(2+nr4-(nr4>tab2)-(nr4>tab3))*1.75/(double)GLOBAL_NCLUSTERS);}
	    o1->rara[nr1][nr2] *= (double)(nowconcentration+1)/(double)SUITE_NCONCENTRATIONS;}
	  if (verbose){ printf("%% %% GLOBAL_ODORra: ");odorfprintf_full(stdout,GLOBAL_ODORra);}
	  switch (nowbicuculline){
	  case 0: /* normal */ 
	    GLOBAL_CS_SRA_SCALE_[GLOBAL_INDEXING_REFILE_sra[VARNAME_REGISTRY_wilson_nAch]] = 1;
	    GLOBAL_CS_SRA_SCALE_[GLOBAL_INDEXING_REFILE_sra[VARNAME_REGISTRY_wilson_gabaA]] = 1;
	    GLOBAL_CS_SRA_SCALE_[GLOBAL_INDEXING_REFILE_sra[VARNAME_REGISTRY_wilson_gabaB]] = 1;
	    break;
	  case 1: /* picrotoxin (no gabaA) */ 
	    GLOBAL_CS_SRA_SCALE_[GLOBAL_INDEXING_REFILE_sra[VARNAME_REGISTRY_wilson_nAch]] = 1;
	    GLOBAL_CS_SRA_SCALE_[GLOBAL_INDEXING_REFILE_sra[VARNAME_REGISTRY_wilson_gabaA]] = 0;
	    GLOBAL_CS_SRA_SCALE_[GLOBAL_INDEXING_REFILE_sra[VARNAME_REGISTRY_wilson_gabaB]] = 1;
	    break;
	  case 2: /* CGP (no gabaB) */ 
	    GLOBAL_CS_SRA_SCALE_[GLOBAL_INDEXING_REFILE_sra[VARNAME_REGISTRY_wilson_nAch]] = 1;
	    GLOBAL_CS_SRA_SCALE_[GLOBAL_INDEXING_REFILE_sra[VARNAME_REGISTRY_wilson_gabaA]] = 1;
	    GLOBAL_CS_SRA_SCALE_[GLOBAL_INDEXING_REFILE_sra[VARNAME_REGISTRY_wilson_gabaB]] = 0;
	    break;
	  case 3: /* CGP+picrotoxin (no gaba) */ 
	    GLOBAL_CS_SRA_SCALE_[GLOBAL_INDEXING_REFILE_sra[VARNAME_REGISTRY_wilson_nAch]] = 1;
	    GLOBAL_CS_SRA_SCALE_[GLOBAL_INDEXING_REFILE_sra[VARNAME_REGISTRY_wilson_gabaA]] = 0;
	    GLOBAL_CS_SRA_SCALE_[GLOBAL_INDEXING_REFILE_sra[VARNAME_REGISTRY_wilson_gabaB]] = 0;
	    break;
	  case 4: /* no nAch */ 
	    GLOBAL_CS_SRA_SCALE_[GLOBAL_INDEXING_REFILE_sra[VARNAME_REGISTRY_wilson_nAch]] = 0;
	    GLOBAL_CS_SRA_SCALE_[GLOBAL_INDEXING_REFILE_sra[VARNAME_REGISTRY_wilson_gabaA]] = 1;
	    GLOBAL_CS_SRA_SCALE_[GLOBAL_INDEXING_REFILE_sra[VARNAME_REGISTRY_wilson_gabaB]] = 1;
	    break;
	  default: break;}
          if (verbose){ raprintf(GLOBAL_CS_SRA_SCALE_,"double",1,GLOBAL_INDEXING_sra_LENGTH,"GLOBAL_CS_SRA_SCALE_ ");}}
	pulse_temp = GLOBAL_time - timera[nowtime];
	INPUT_PULSE = (pulse_temp>pulse_waitfor && pulse_temp<(pulse_waitfor+pulse_howlong));}
      if (nowtime>=ntimes-1){
	if (fabs(SUITE_ONSET_TIME-timera[nowtime])>0.1){
	  if (verbose){ printf(" %% time %f, entering %d_th phase of suite seven, finishing\n",GLOBAL_time,nowtime);}
	  SUITE_ONSET_TIME = timera[nowtime]; OUTPUT_DUMP_EVERY = 0; /* don't bother dumping data */
	  if (nowtime>0){ 
	    if (verbose){ printf(" %% depositing old power\n");}
	    suite_7_power_deposit_helper(oldinstance,oldconcentration,oldodor,oldbicuculline,oldtime);
	    if (verbose){ printf(" %% remaking new power\n");}
	    powertfree(GLOBAL_POWER); GLOBAL_POWER=NULL;
	    GLOBAL_POWER = powermake(GLOBAL_Nra,SUITE_NSECONDS*1024,GLOBAL_POWER_INDEXING_NTYPE_LENGTH,GLOBAL_POWER_INDEXING_NTYPE_CHECKOUT,GLOBAL_POWER_INDEXING_NTYPE_REFILE,GLOBAL_POWER_INDEXING_NVAR_LENGTH,GLOBAL_POWER_INDEXING_NVAR_CHECKOUT,GLOBAL_POWER_INDEXING_NVAR_REFILE,GLOBAL_POWER_maxra_,GLOBAL_POWER_minra_,0/* GLOBAL_POWER_CORRELATION_BOTHER */);}
	  sprintf(tempchar,"%srecord",gs2);
	  suite_7_power_process_helper(tempchar);
	  if (SUITE_7_CLEANUP){ for (nt=0;nt<GLOBAL_POWER->indexing_ntype_length;nt++){ for (nv=0;nv<GLOBAL_POWER->indexing_nvar_length;nv++){ for (nowinstance=0;nowinstance<SUITE_NINSTANCES;nowinstance++){ for (nowconcentration=0;nowconcentration<SUITE_NCONCENTRATIONS;nowconcentration++){ for (nowodor=0;nowodor<SUITE_NODORS;nowodor++){ for (nowbicuculline=0;nowbicuculline<SUITE_NBICUCULLINES;nowbicuculline++){ sprintf(command,"""rm"" ./dir_temp_suite_7_%srecord/power_deposit_type_%s_var_%s_instance_%d_concentration_%d_odor_%d_bicuculline_%d_%srecord;",gs2,GLOBAL_TYPENAMES[GLOBAL_POWER->indexing_ntype_checkout[nt]],GLOBAL_VARNAMES[GLOBAL_POWER->indexing_nvar_checkout[nv]],nowinstance,nowconcentration,nowodor,nowbicuculline,gs2); if (verbose){ printf("%s\n",command);} system(command);}}}}}}}
	  exit(EXIT_SUCCESS);
	  GLOBAL_TF = SUITE_ONSET_TIME+1;}}}
    tfree(timera);timera=NULL;
    break;
  case 8: 
    /* for use with wilson_type neurons (i.e., GLOBAL_NEURON_MODEL==6)
       runs multiple odors for multiple trials over various concentrations. 
       the ORNs are stimulated directly, and there is no true combinatorial odor.
       different odors stimulate (a single) different cluster.
       each concentration corresponds to a different frequency of antennal nerve stimulation.
       during the "stimulus off" periods the frequency is 8Hz (period 128 ms).
       during the "stimulus on" periods the frequency ranges from 8Hz to 256 Hz. (8,16,32,64,128,256), a total of 6 concentrations.
       The oscillation phase is lined up so that a odor peak occurs exactly at pulse_waitfor, regardless of the stimulus frequency.
       ninstances trials, 
       maximum(6,nconcentration) concentrations, 
       up to minimum(suite_nodors,global_nclusters) combinatorially different odors, indexed by nowodor.
       minimum(5,nbicucullines) alternate architectures.
       nbicuculline==0: standard architecture
       nbicuculline==1: gabaA blocked
       nbicuculline==2: gabaB blocked 
       nbicuculline==3: gabaA and gabaB blocked
       nbicuculline==4: nAch blocked */
    ORN_BACKRATE = 0;
    CS_ORN_wilson_stim_[TYPENAME_REGISTRY_wilson_ORN]= SUITE_8_wilson_axon_stim;
    SUITE_NSECONDS = maximum(2,SUITE_NSECONDS);
    SUITE_NODORS = minimum(GLOBAL_NCLUSTERS*maximum(1,GLOBAL_NCLUSTERS-1),SUITE_NODORS);
    SUITE_NCONCENTRATIONS = minimum(6,SUITE_NCONCENTRATIONS);
    SUITE_NBICUCULLINES = minimum(5,SUITE_NBICUCULLINES);
    SUITE_POWER_BOTHER = 1;
    ntimes = SUITE_NINSTANCES*SUITE_NCONCENTRATIONS*SUITE_NODORS*SUITE_NBICUCULLINES + 1;
    ninputs = SUITE_NCONCENTRATIONS*SUITE_NODORS*SUITE_NBICUCULLINES;
    timera = (double *) tcalloc(ntimes,sizeof(double));
    for (nt=0;nt<ntimes;nt++){ timera[nt] = nt*1024*SUITE_NSECONDS;}
    pulse_waitfor = 2048;
    pulse_howlong = 1024;
    INPUT_PULSE = 1;
    if (t==0){ /* clear everything */
      if (verbose){ printf(" %% \n");}
      if (verbose){ printf(" %% initializing suite eight\n");}
      if (verbose){ printf(" times are: ");raprintf(timera,"double",1,ntimes,"timera ");}
      GLOBAL_TF=0;
      SUITE_ONSET_TIME=-1;
      if (POWER_BOTHER!=0){ powertfree(GLOBAL_POWER); GLOBAL_POWER=NULL; POWER_BOTHER=0;}
      if (PTREE_BOTHER!=0){ ptreetfree(GLOBAL_PTREE); GLOBAL_PTREE=NULL; GLOBAL_PTREE_BITBYBIT=0; GLOBAL_PTREE_BITBYBIT_CURRENT_NUMBER=0; PTREE_BOTHER=0;}
      if (CLUSTERDATA_BOTHER!=0){ clusterdatatfree(GLOBAL_CDRA); GLOBAL_CDRA=NULL; CLUSTERDATA_BOTHER=0;}
      POWER_BOTHER=1;
      GLOBAL_POWER = powermake(GLOBAL_Nra,SUITE_NSECONDS*1024,GLOBAL_POWER_INDEXING_NTYPE_LENGTH,GLOBAL_POWER_INDEXING_NTYPE_CHECKOUT,GLOBAL_POWER_INDEXING_NTYPE_REFILE,GLOBAL_POWER_INDEXING_NVAR_LENGTH,GLOBAL_POWER_INDEXING_NVAR_CHECKOUT,GLOBAL_POWER_INDEXING_NVAR_REFILE,GLOBAL_POWER_maxra_,GLOBAL_POWER_minra_,0/* GLOBAL_POWER_CORRELATION_BOTHER */);}
    else /* if (t>0) */{
      nowtime = 0; for (nt=0;nt<ntimes-1;nt++){ if (t>=timera[nt] && t<timera[nt+1]){ nowtime=nt;}} 
      if (t>=timera[ntimes-1]){ nowtime=ntimes-1;}
      nowbicuculline = nowtime/(SUITE_NINSTANCES*SUITE_NCONCENTRATIONS*SUITE_NODORS);
      nowodor = (nowtime-nowbicuculline*SUITE_NINSTANCES*SUITE_NCONCENTRATIONS*SUITE_NODORS)/(SUITE_NINSTANCES*SUITE_NCONCENTRATIONS);
      nowconcentration = (nowtime-nowbicuculline*SUITE_NINSTANCES*SUITE_NCONCENTRATIONS*SUITE_NODORS-nowodor*SUITE_NINSTANCES*SUITE_NCONCENTRATIONS)/(SUITE_NINSTANCES);
      nowinstance = (nowtime-nowbicuculline*SUITE_NINSTANCES*SUITE_NCONCENTRATIONS*SUITE_NODORS-nowodor*SUITE_NINSTANCES*SUITE_NCONCENTRATIONS-nowconcentration*SUITE_NINSTANCES)/(1);
      nowinput = nowconcentration + nowodor*SUITE_NCONCENTRATIONS + nowbicuculline*SUITE_NCONCENTRATIONS*SUITE_NODORS;
      oldtime = maximum(0,nowtime-1);
      oldbicuculline = oldtime/(SUITE_NINSTANCES*SUITE_NCONCENTRATIONS*SUITE_NODORS);
      oldodor = (oldtime-oldbicuculline*SUITE_NINSTANCES*SUITE_NCONCENTRATIONS*SUITE_NODORS)/(SUITE_NINSTANCES*SUITE_NCONCENTRATIONS);
      oldconcentration = (oldtime-oldbicuculline*SUITE_NINSTANCES*SUITE_NCONCENTRATIONS*SUITE_NODORS-oldodor*SUITE_NINSTANCES*SUITE_NCONCENTRATIONS)/(SUITE_NINSTANCES);
      oldinstance = (oldtime-oldbicuculline*SUITE_NINSTANCES*SUITE_NCONCENTRATIONS*SUITE_NODORS-oldodor*SUITE_NINSTANCES*SUITE_NCONCENTRATIONS-oldconcentration*SUITE_NINSTANCES)/(1);
      oldinput = oldconcentration + oldodor*SUITE_NCONCENTRATIONS + oldbicuculline*SUITE_NCONCENTRATIONS*SUITE_NODORS;
      if (nowtime<ntimes-1){
	if (fabs(SUITE_ONSET_TIME-timera[nowtime])>0.1){
	  if (verbose){ printf(" %% time %f, entering %d_th phase of suite eight, bicuculline %d odor %d concentration %d instance %d\n",GLOBAL_time,nowtime,nowbicuculline,nowodor,nowconcentration,nowinstance);}
	  SUITE_ONSET_TIME = timera[nowtime]; OUTPUT_DUMP_EVERY = 0; /* don't bother dumping data */
	  if (nowtime>0){ 
	    if (verbose){ printf(" %% depositing old power\n");}
	    suite_8_power_deposit_helper(oldinstance,oldconcentration,oldodor,oldbicuculline,oldtime);
	    if (verbose){ printf(" %% remaking new power\n");}
	    powertfree(GLOBAL_POWER); GLOBAL_POWER=NULL;
	    GLOBAL_POWER = powermake(GLOBAL_Nra,SUITE_NSECONDS*1024,GLOBAL_POWER_INDEXING_NTYPE_LENGTH,GLOBAL_POWER_INDEXING_NTYPE_CHECKOUT,GLOBAL_POWER_INDEXING_NTYPE_REFILE,GLOBAL_POWER_INDEXING_NVAR_LENGTH,GLOBAL_POWER_INDEXING_NVAR_CHECKOUT,GLOBAL_POWER_INDEXING_NVAR_REFILE,GLOBAL_POWER_maxra_,GLOBAL_POWER_minra_,0/* GLOBAL_POWER_CORRELATION_BOTHER */);}
	  switch (nowbicuculline){
	  case 0: /* normal */ 
	    GLOBAL_CS_SRA_SCALE_[GLOBAL_INDEXING_REFILE_sra[VARNAME_REGISTRY_wilson_nAch]] = 1;
	    GLOBAL_CS_SRA_SCALE_[GLOBAL_INDEXING_REFILE_sra[VARNAME_REGISTRY_wilson_gabaA]] = 1;
	    GLOBAL_CS_SRA_SCALE_[GLOBAL_INDEXING_REFILE_sra[VARNAME_REGISTRY_wilson_gabaB]] = 1;
	    break;
	  case 1: /* picrotoxin (no gabaA) */ 
	    GLOBAL_CS_SRA_SCALE_[GLOBAL_INDEXING_REFILE_sra[VARNAME_REGISTRY_wilson_nAch]] = 1;
	    GLOBAL_CS_SRA_SCALE_[GLOBAL_INDEXING_REFILE_sra[VARNAME_REGISTRY_wilson_gabaA]] = 0;
	    GLOBAL_CS_SRA_SCALE_[GLOBAL_INDEXING_REFILE_sra[VARNAME_REGISTRY_wilson_gabaB]] = 1;
	    break;
	  case 2: /* CGP (no gabaB) */ 
	    GLOBAL_CS_SRA_SCALE_[GLOBAL_INDEXING_REFILE_sra[VARNAME_REGISTRY_wilson_nAch]] = 1;
	    GLOBAL_CS_SRA_SCALE_[GLOBAL_INDEXING_REFILE_sra[VARNAME_REGISTRY_wilson_gabaA]] = 1;
	    GLOBAL_CS_SRA_SCALE_[GLOBAL_INDEXING_REFILE_sra[VARNAME_REGISTRY_wilson_gabaB]] = 0;
	    break;
	  case 3: /* CGP+picrotoxin (no gaba) */ 
	    GLOBAL_CS_SRA_SCALE_[GLOBAL_INDEXING_REFILE_sra[VARNAME_REGISTRY_wilson_nAch]] = 1;
	    GLOBAL_CS_SRA_SCALE_[GLOBAL_INDEXING_REFILE_sra[VARNAME_REGISTRY_wilson_gabaA]] = 0;
	    GLOBAL_CS_SRA_SCALE_[GLOBAL_INDEXING_REFILE_sra[VARNAME_REGISTRY_wilson_gabaB]] = 0;
	    break;
	  case 4: /* no nAch */ 
	    GLOBAL_CS_SRA_SCALE_[GLOBAL_INDEXING_REFILE_sra[VARNAME_REGISTRY_wilson_nAch]] = 0;
	    GLOBAL_CS_SRA_SCALE_[GLOBAL_INDEXING_REFILE_sra[VARNAME_REGISTRY_wilson_gabaA]] = 1;
	    GLOBAL_CS_SRA_SCALE_[GLOBAL_INDEXING_REFILE_sra[VARNAME_REGISTRY_wilson_gabaB]] = 1;
	    break;
	  default: break;}
          if (verbose){ raprintf(GLOBAL_CS_SRA_SCALE_,"double",1,GLOBAL_INDEXING_sra_LENGTH,"GLOBAL_CS_SRA_SCALE_ ");}}
	pulse_temp = GLOBAL_time - timera[nowtime];
	o1 = GLOBAL_ODORra; nr1 = TYPENAME_REGISTRY_wilson_ORN;
	tab2 = 0; tab3 = GLOBAL_NCLUSTERS; periodify("int",&nowodor,&tab2,&tab3,&tab);
	if (verbose>1){ printf(" %% cluster %d chosen by periodize(nowodor,0,GLOBAL_NCLUSTERS)\n",tab);}
	for (nr2=0;nr2<o1->lengthra[nr1];nr2++){ 
	  nr3 = nr2%(o1->lengthra[nr1]/GLOBAL_NCLUSTERS); nr4 = nr2/(o1->lengthra[nr1]/GLOBAL_NCLUSTERS);
	  if (nr4==tab){ o1->rara[nr1][nr2] = o1->base*1;}
	  else /* if (nr4!=tab) */{ o1->rara[nr1][nr2] = o1->base*0;}
	  if (pulse_temp<=pulse_waitfor){ 
	    pulse_time = pulse_waitfor - pulse_temp;
	    pulse_phase_max = 0.5*128.0; pulse_phase_min = -0.5*128.0;
	    periodify("double",&pulse_time,&pulse_phase_min,&pulse_phase_max,&pulse_phase);
	    o1->rara[nr1][nr2] = o1->base*(nr4==tab)*(abs(pulse_phase)<2/*stimulus width in ms*/);}
	  else if (pulse_temp>pulse_waitfor && pulse_temp<=pulse_waitfor+pulse_howlong){
	    pulse_time = pulse_temp - pulse_waitfor;
	    pulse_phase_max = 0.5*128.0/pow(2,nowconcentration); pulse_phase_min = -0.5*128.0/pow(2,nowconcentration);
	    periodify("double",&pulse_time,&pulse_phase_min,&pulse_phase_max,&pulse_phase);
	    o1->rara[nr1][nr2] = o1->base*(nr4==tab)*(abs(pulse_phase)<2/*stimulus width in ms*/);}
	  else if (pulse_temp>pulse_waitfor+pulse_howlong){
	    o1->rara[nr1][nr2] = 0;}}
	if (verbose>1){ printf(" %% cluster %d stimulated with period %f, odor %f\n",tab,128.0/pow(2,nowconcentration),o1->rara[nr1][0]);}}
      if (nowtime>=ntimes-1){
	if (fabs(SUITE_ONSET_TIME-timera[nowtime])>0.1){
	  if (verbose){ printf(" %% time %f, entering %d_th phase of suite eight, finishing\n",GLOBAL_time,nowtime);}
	  SUITE_ONSET_TIME = timera[nowtime]; OUTPUT_DUMP_EVERY = 0; /* don't bother dumping data */
	  if (nowtime>0){ 
	    if (verbose){ printf(" %% depositing old power\n");}
	    suite_8_power_deposit_helper(oldinstance,oldconcentration,oldodor,oldbicuculline,oldtime);
	    if (verbose){ printf(" %% remaking new power\n");}
	    powertfree(GLOBAL_POWER); GLOBAL_POWER=NULL;
	    GLOBAL_POWER = powermake(GLOBAL_Nra,SUITE_NSECONDS*1024,GLOBAL_POWER_INDEXING_NTYPE_LENGTH,GLOBAL_POWER_INDEXING_NTYPE_CHECKOUT,GLOBAL_POWER_INDEXING_NTYPE_REFILE,GLOBAL_POWER_INDEXING_NVAR_LENGTH,GLOBAL_POWER_INDEXING_NVAR_CHECKOUT,GLOBAL_POWER_INDEXING_NVAR_REFILE,GLOBAL_POWER_maxra_,GLOBAL_POWER_minra_,0/* GLOBAL_POWER_CORRELATION_BOTHER */);}
	  sprintf(tempchar,"%srecord",gs2);
	  suite_8_power_process_helper2(tempchar,pulse_waitfor,minimum(pulse_howlong,1024));
	  suite_8_power_process_helper(tempchar);
	  if (SUITE_8_CLEANUP){ for (nt=0;nt<GLOBAL_POWER->indexing_ntype_length;nt++){ for (nv=0;nv<GLOBAL_POWER->indexing_nvar_length;nv++){ for (nowinstance=0;nowinstance<SUITE_NINSTANCES;nowinstance++){ for (nowconcentration=0;nowconcentration<SUITE_NCONCENTRATIONS;nowconcentration++){ for (nowodor=0;nowodor<SUITE_NODORS;nowodor++){ for (nowbicuculline=0;nowbicuculline<SUITE_NBICUCULLINES;nowbicuculline++){ sprintf(command,"""rm"" ./dir_temp_suite_8_%srecord/power_deposit_type_%s_var_%s_instance_%d_concentration_%d_odor_%d_bicuculline_%d_%srecord;",gs2,GLOBAL_TYPENAMES[GLOBAL_POWER->indexing_ntype_checkout[nt]],GLOBAL_VARNAMES[GLOBAL_POWER->indexing_nvar_checkout[nv]],nowinstance,nowconcentration,nowodor,nowbicuculline,gs2); if (verbose){ printf("%s\n",command);} system(command);}}}}}}}
	  exit(EXIT_SUCCESS);
	  GLOBAL_TF = SUITE_ONSET_TIME+1;}}}
    tfree(timera);timera=NULL;
    break;
  default: break;}
}

/* here are lattice3d functions */

struct spack3d * spack3dmake(struct lattice3d *l3d,int index1,int index2,int index3,double xpos,double ypos,double zpos,double value)
{
  struct spack3d *s=NULL;
  s = (struct spack3d *) tcalloc(1,sizeof(struct spack3d));
  s->l3d = l3d;
  s->i = index1; s->j = index2; s->k = index3; s->x = xpos; s->y = ypos; s->z = zpos;
  s->value = value;
  s->L = llistmake();
  s->adjacency_ra_s=NULL;
  s->adjacency_ra_v=NULL;
  s->temp = NULL;
  return s;
}

void spack3dprintf(struct spack3d *s,void (*void_printf)(void *,void *),void *void_parameters)
{
  int nr=0;
  struct spack3d *s2=NULL;
  printf("s(%d,%d,%d) at (%0.1f,%0.1f,%0.1f), value %f\n",s->i,s->j,s->k,s->x,s->y,s->z,s->value);
  if (s->temp!=NULL && void_printf!=NULL){ void_printf(s->temp,void_parameters);}
  if (s->adjacency_ra_s!=NULL && s->adjacency_ra_v!=NULL){
    printf("\t neighbors:\n"); 
    for (nr=0;nr<s->L->length;nr++){ 
      s2 = s->adjacency_ra_s[nr];
      printf("\t\ts(%d,%d,%d), value %f\n",s2->i,s2->j,s2->k,s->adjacency_ra_v[nr]);}}
}

void spack3dtfree(struct spack3d *s,void (*free_function)(void *,void *),void *void_parameters)
{
  if (free_function!=NULL && s->temp!=NULL){ (*free_function)(s->temp,void_parameters); s->temp=NULL;}
  llisttfree(s->L); s->L=NULL;
  if (s->adjacency_ra_s!=NULL){ tfree(s->adjacency_ra_s); s->adjacency_ra_s=NULL;}
  if (s->adjacency_ra_v!=NULL){ tfree(s->adjacency_ra_v); s->adjacency_ra_v=NULL;}
  tfree(s); s=NULL;
}

struct spack3d * lattice3dget(struct lattice3d *l3d,int nr1,int nr2,int nr3)
{
  int verbose=0;
  int tab=0;
  struct spack3d *s=NULL;
  if (verbose){
    if (nr1<0 || nr1 >= l3d->input1max){ printf(" %% warning! nr1 %d out of bounds in lattice3d\n",nr1);}
    if (nr2<0 || nr2 >= l3d->input2max_scale){ printf(" %% warning! nr2 %d out of bounds in lattice3d\n",nr2);}
    if (nr3<0 || nr3 >= l3d->input3max_scale){ printf(" %% warning! nr3 %d out of bounds in lattice3d\n",nr3);}}
  if (nr1>=0 && nr1<l3d->input1max && nr2>=0 && nr2<l3d->input2max_scale && nr3>=0 && nr3<l3d->input3max_scale){
    tab = nr1+nr2*l3d->input1max+nr3*l3d->input1max*l3d->input2max_scale; s=l3d->lattice[tab];}
  return s;
}

void spack3dlink(struct spack3d *s1,struct spack3d *s2)
{
  double d = 0;
  double tolerance=0.000001;
  if (s1!=NULL && s2!=NULL){
    d=sqrt(pow(s1->x-s2->x,2)+pow(s1->y-s2->y,2)+pow(s1->z-s2->z,2));
    if (fabs(d-1)>tolerance){ printf(" %% warning! s(%f,%f,%f) and s(%f,%f,%f) are %f away in spack3dlink\n",s1->x,s1->y,s1->z,s2->x,s2->y,s2->z,d);}
    if (!isin(s1->L,s2)){ litemadd(s1->L,s2);}
    if (!isin(s2->L,s1)){ litemadd(s2->L,s1);}
    if (s1->L->length>12){ printf(" %% warning! s(%f,%f,%f)->L->length=%d\n",s1->x,s1->y,s1->z,s1->L->length);}
    if (s2->L->length>12){ printf(" %% warning! s(%f,%f,%f)->L->length=%d\n",s2->x,s2->y,s2->z,s2->L->length);}}
}

double vra2vra_diff_sym(void *vra1,void *vra2,void *void_parameters)
{
  /* parameter list:
     (void *) (int *) index
     (void *) (int *) length
     (void *) (double *) prefactor_matrix
     (void *) (int *) prefactor_rows
  */
  void ** vrara=NULL;
  int *length=NULL;
  void ** vra=NULL;
  double *dra1=NULL,*dra2=NULL,*prefactor_matrix=NULL,*dra3=NULL,*dra4=NULL,*dra5=NULL;
  int *index=NULL,*prefactor_rows=NULL;
  double output=0;
  if (void_parameters!=NULL){
    vrara = (void **) void_parameters;
    index = vrara[0]; length = vrara[1]; prefactor_matrix = vrara[2]; prefactor_rows = vrara[3];
    vra = (void **)vra1; dra1 = (double *)vra[*index]; vra = (void **)vra2; dra2 = (double *)vra[*index];
    if (prefactor_matrix==NULL){ 
      dra5=ra2ra_minus(dra1,dra2,*length);    
      output = ra_norm(dra5,*length);}
    else if (prefactor_matrix!=NULL){ 
      dra3 = ra2ra_matrix_multiply(prefactor_matrix,*prefactor_rows,*length,0,dra1,*length,1,0);
      dra4 = ra2ra_matrix_multiply(prefactor_matrix,*prefactor_rows,*length,0,dra2,*length,1,0);
      dra5 = ra2ra_minus(dra3,dra4,*length); tfree(dra3);dra3=NULL; tfree(dra4);dra4=NULL;
      output = ra_norm(dra5,*prefactor_rows);}
    tfree(dra5);dra5=NULL;}
  return output;
}

double vra2vra_rq_sym(void *vra1,void *vra2,void *void_parameters)
{
  /* parameter list:
     (void *) (int *) index
     (void *) (int *) length
     (void *) (double *) prefactor_matrix
     (void *) (int *) prefactor_rows
  */
  void ** vrara=NULL;
  int *length=NULL;
  void ** vra=NULL;
  double *dra1=NULL,*dra2=NULL,*prefactor_matrix=NULL;
  int *index=NULL,*prefactor_rows=NULL;
  int nr=0,nc=0;
  double output=1,sum1=0,sum2=0,sum1a=0,sum2a=0;
  if (void_parameters!=NULL){
    vrara = (void **) void_parameters;
    index = vrara[0]; length = vrara[1]; prefactor_matrix = vrara[2]; prefactor_rows = vrara[3];
    vra = (void **)vra1; dra1 = (double *)vra[*index]; vra = (void **)vra2; dra2 = (double *)vra[*index];
    output=0;sum1a=0;sum2a=0;
    for (nr=0;nr<(prefactor_rows==NULL ? *length : *prefactor_rows);nr++){ 
      sum1=0;
      if (prefactor_matrix==NULL){ sum1 = dra1[nr];}
      else{ for (nc=0;nc<*length;nc++){ sum1+=prefactor_matrix[nr+nc*(*prefactor_rows)]*dra1[nc];}}
      sum2=0;
      if (prefactor_matrix==NULL){ sum2 = dra2[nr];}
      else{ for (nc=0;nc<*length;nc++){ sum2+=prefactor_matrix[nr+nc*(*prefactor_rows)]*dra2[nc];}}
      output += sum1*sum2;
      sum1a+=pow(sum1,2); sum2a+=pow(sum2,2);}
    output = output/sqrt(sum1a*sum2a);}
  return 1-output;
}

void lattice3d_edge_eval(struct lattice3d *l3d,double (*void_eval)(void *,void *,void *),void * void_parameters)
{
  int nr1=0,nr2=0,nr3=0,tab=0;
  struct spack3d *s=NULL;
  for (nr1=0;nr1<l3d->input1max;nr1++){ for (nr2=0;nr2<l3d->input2max_scale;nr2++){ for (nr3=0;nr3<l3d->input3max_scale;nr3++){
    s = lattice3dget(l3d,nr1,nr2,nr3);
    if (s->adjacency_ra_s!=NULL && s->adjacency_ra_v!=NULL){
      for (tab=0;tab<s->L->length;tab++){
	s->adjacency_ra_v[tab] = void_eval(s->temp,s->adjacency_ra_s[tab]->temp,void_parameters);}}}}}
}

void latticelink(struct lattice3d *l3d)
{
  int nr1=0,nr2=0,nr3=0,tab=0;
  struct spack3d *s=NULL;
  struct litem *l=NULL;
  if (l3d->tet_vs_cube==1){
    for (nr1=0;nr1<l3d->input1max;nr1++){ for (nr2=0;nr2<l3d->input2max_scale;nr2++){ for (nr3=0;nr3<l3d->input3max_scale;nr3++){
      s = lattice3dget(l3d,nr1,nr2,nr3);
      if ((s->k)%2==0){      
	if ((s->j)%2==0){ spack3dlink(s,lattice3dget(l3d,nr1+1,nr2,nr3)); spack3dlink(s,lattice3dget(l3d,nr1-1,nr2,nr3)); spack3dlink(s,lattice3dget(l3d,nr1,nr2+1,nr3)); spack3dlink(s,lattice3dget(l3d,nr1+1,nr2+1,nr3)); spack3dlink(s,lattice3dget(l3d,nr1,nr2-1,nr3)); spack3dlink(s,lattice3dget(l3d,nr1+1,nr2-1,nr3)); spack3dlink(s,lattice3dget(l3d,nr1,nr2,nr3+1)); spack3dlink(s,lattice3dget(l3d,nr1,nr2-1,nr3+1)); spack3dlink(s,lattice3dget(l3d,nr1+1,nr2-1,nr3+1)); spack3dlink(s,lattice3dget(l3d,nr1,nr2,nr3-1)); spack3dlink(s,lattice3dget(l3d,nr1,nr2-1,nr3-1)); spack3dlink(s,lattice3dget(l3d,nr1+1,nr2-1,nr3-1));}
	else if ((s->j)%2==1){ spack3dlink(s,lattice3dget(l3d,nr1+1,nr2,nr3)); spack3dlink(s,lattice3dget(l3d,nr1-1,nr2,nr3)); spack3dlink(s,lattice3dget(l3d,nr1,nr2+1,nr3)); spack3dlink(s,lattice3dget(l3d,nr1-1,nr2+1,nr3)); spack3dlink(s,lattice3dget(l3d,nr1,nr2-1,nr3)); spack3dlink(s,lattice3dget(l3d,nr1-1,nr2-1,nr3)); spack3dlink(s,lattice3dget(l3d,nr1,nr2,nr3+1)); spack3dlink(s,lattice3dget(l3d,nr1,nr2-1,nr3+1)); spack3dlink(s,lattice3dget(l3d,nr1-1,nr2-1,nr3+1)); spack3dlink(s,lattice3dget(l3d,nr1,nr2,nr3-1)); spack3dlink(s,lattice3dget(l3d,nr1,nr2-1,nr3-1)); spack3dlink(s,lattice3dget(l3d,nr1-1,nr2-1,nr3-1));}}
      else if ((s->k)%2==1){
	if ((s->j)%2==0){ spack3dlink(s,lattice3dget(l3d,nr1+1,nr2,nr3)); spack3dlink(s,lattice3dget(l3d,nr1-1,nr2,nr3)); spack3dlink(s,lattice3dget(l3d,nr1,nr2+1,nr3)); spack3dlink(s,lattice3dget(l3d,nr1+1,nr2+1,nr3)); spack3dlink(s,lattice3dget(l3d,nr1,nr2-1,nr3)); spack3dlink(s,lattice3dget(l3d,nr1+1,nr2-1,nr3)); spack3dlink(s,lattice3dget(l3d,nr1,nr2,nr3+1)); spack3dlink(s,lattice3dget(l3d,nr1,nr2+1,nr3+1)); spack3dlink(s,lattice3dget(l3d,nr1+1,nr2+1,nr3+1)); spack3dlink(s,lattice3dget(l3d,nr1,nr2,nr3-1)); spack3dlink(s,lattice3dget(l3d,nr1,nr2+1,nr3-1)); spack3dlink(s,lattice3dget(l3d,nr1+1,nr2+1,nr3-1));}
	else if ((s->j)%2==1){ spack3dlink(s,lattice3dget(l3d,nr1+1,nr2,nr3)); spack3dlink(s,lattice3dget(l3d,nr1-1,nr2,nr3)); spack3dlink(s,lattice3dget(l3d,nr1,nr2+1,nr3)); spack3dlink(s,lattice3dget(l3d,nr1-1,nr2+1,nr3)); spack3dlink(s,lattice3dget(l3d,nr1,nr2-1,nr3)); spack3dlink(s,lattice3dget(l3d,nr1-1,nr2-1,nr3)); spack3dlink(s,lattice3dget(l3d,nr1,nr2,nr3+1)); spack3dlink(s,lattice3dget(l3d,nr1,nr2+1,nr3+1)); spack3dlink(s,lattice3dget(l3d,nr1-1,nr2+1,nr3+1)); spack3dlink(s,lattice3dget(l3d,nr1,nr2,nr3-1)); spack3dlink(s,lattice3dget(l3d,nr1,nr2+1,nr3-1)); spack3dlink(s,lattice3dget(l3d,nr1-1,nr2+1,nr3-1));}}}}}}
  else /* if (l3d->tet_vs_cube==0) */{
    for (nr1=0;nr1<l3d->input1max;nr1++){ for (nr2=0;nr2<l3d->input2max_scale;nr2++){ for (nr3=0;nr3<l3d->input3max_scale;nr3++){
      s = lattice3dget(l3d,nr1,nr2,nr3);
      spack3dlink(s,lattice3dget(l3d,nr1+1,nr2,nr3)); spack3dlink(s,lattice3dget(l3d,nr1-1,nr2,nr3)); spack3dlink(s,lattice3dget(l3d,nr1,nr2+1,nr3)); spack3dlink(s,lattice3dget(l3d,nr1,nr2-1,nr3)); spack3dlink(s,lattice3dget(l3d,nr1,nr2,nr3+1)); spack3dlink(s,lattice3dget(l3d,nr1,nr2,nr3-1));}}}}
  for (nr1=0;nr1<l3d->input1max;nr1++){ for (nr2=0;nr2<l3d->input2max_scale;nr2++){ for (nr3=0;nr3<l3d->input3max_scale;nr3++){
    s = lattice3dget(l3d,nr1,nr2,nr3);
    if (s->L->length>0){
      s->adjacency_ra_s = (struct spack3d **) tcalloc(s->L->length,sizeof(struct spack3d *));
      s->adjacency_ra_v = (double *) tcalloc(s->L->length,sizeof(double));
      l=s->L->first; tab=0; while (l!=NULL){ s->adjacency_ra_s[tab] = (struct spack3d *) l->item; l=l->child; tab+=1;}}}}}
}

struct lattice3d * lattice3dmake(int tet_vs_cube,int input1max,int input2max,int input3max)
{
  struct lattice3d *l3d=NULL;
  int nr1=0,nr2=0,nr3=0,tab=0;
  double xpos=0,ypos=0,zpos=0;
  l3d = (struct lattice3d *) tcalloc(1,sizeof(struct lattice3d));
  l3d->tet_vs_cube = tet_vs_cube;
  l3d->input1max = input1max;
  l3d->input2max = input2max;
  l3d->input3max = input3max;
  l3d->input2max_scale = l3d->tet_vs_cube ? (int)(l3d->input2max/(sqrt(3)/2.0)) : l3d->input2max;
  l3d->input3max_scale = l3d->tet_vs_cube ? (int)(l3d->input3max/sqrt(2.0/3.0)) : l3d->input3max;
  l3d->total_length = l3d->input1max*l3d->input2max_scale*l3d->input3max_scale;
  l3d->lattice = (struct spack3d **) tcalloc(l3d->total_length,sizeof(struct spack3d *));
  for (nr1=0;nr1<l3d->input1max;nr1++){ for (nr2=0;nr2<l3d->input2max_scale;nr2++){ for (nr3=0;nr3<l3d->input3max_scale;nr3++){
    lattice3d_int2input_or_input2int(1,l3d->tet_vs_cube,l3d->input1max,l3d->input2max,l3d->input3max,&nr1,&nr2,&nr3,&xpos,&ypos,&zpos);
    tab = nr1+nr2*l3d->input1max+nr3*l3d->input1max*l3d->input2max_scale;
    l3d->lattice[tab] = spack3dmake(l3d,nr1,nr2,nr3,xpos,ypos,zpos,0);}}}
  latticelink(l3d);
  l3d->nedges=0;
  l3d->edge_v_max=0;
  l3d->edge_v_min=0;
  l3d->edge_v_mean=0;
  l3d->edge_v_stdev=0;
  return l3d;
}

void lattice3d_edge_v_stats(struct lattice3d *l3d)
{
  int nr1=0,nr2=0,nr3=0,tab=0;
  struct spack3d *s=NULL;
  int max_set=0,min_set=0;
  double val=0;
  l3d->nedges=0;
  l3d->edge_v_max=0;
  l3d->edge_v_min=0;
  l3d->edge_v_mean=0;
  l3d->edge_v_stdev=0;
  for (nr1=0;nr1<l3d->input1max;nr1++){ for (nr2=0;nr2<l3d->input2max_scale;nr2++){ for (nr3=0;nr3<l3d->input3max_scale;nr3++){
    s = lattice3dget(l3d,nr1,nr2,nr3);
    if (s->adjacency_ra_s!=NULL && s->adjacency_ra_v!=NULL){
      l3d->nedges += s->L->length;
      for (tab=0;tab<s->L->length;tab++){
	val = s->adjacency_ra_v[tab];
	if (max_set){ l3d->edge_v_max = maximum(val,l3d->edge_v_max);} else if (!max_set){ l3d->edge_v_max = val; max_set=1;}
	if (min_set){ l3d->edge_v_min = minimum(val,l3d->edge_v_min);} else if (!min_set){ l3d->edge_v_min = val; min_set=1;}
	l3d->edge_v_mean += val;
	l3d->edge_v_stdev += val*val;}}}}}
  l3d->edge_v_mean /= l3d->nedges;
  l3d->edge_v_stdev = sqrt(l3d->edge_v_stdev/l3d->nedges - pow(l3d->edge_v_mean,2));
}

void spack3d_temp_vector_dump(void *v1,FILE *fp,double maxdia,double xside,double yside,double xoffset,double yoffset,void *void_parameters)
{
  /* assumes that dra[*vector_length] holds a relevant magnitude 
     parameter list:
     (void *) (int *) index
     (void *) (int *) vector_length
     (void *) (int *) arrow_flag
     (void *) (double *) max
     (void *) (double *) min
   */
  int *index=NULL,*vector_length=NULL,*arrow_flag=NULL;
  void **vra=NULL,**vra_p=NULL;
  double *dra=NULL,*max=NULL,*min=NULL;
  double xord=0,yord=0,zord=0,xord2=0,yord2=0;
  double rcolor=0,gcolor=0,bcolor=0;
  int depth=1;
  int colorcode=0;
  if (v1!=NULL && void_parameters!=NULL){
    vra_p=(void **) void_parameters; index = (int *) vra_p[0]; vector_length = (int *) vra_p[1]; arrow_flag = (int *) vra_p[2];
    max = (double *) vra_p[3]; min = (double *) vra_p[4];
    vra = (void **) v1; dra = (double *) vra[*index];
    xord = (*vector_length>=0 ? dra[0]: 0); yord = (*vector_length>=1 ? dra[1]: 0); zord = (*vector_length>=2 ? dra[2]: 0);
    xord2 = -yord/4.0; yord2 = xord/4.0;
    colorscale(0,dra[*vector_length],(*max>*min?*max:1),(*max>*min?*min:0),&rcolor,&gcolor,&bcolor); 
    colorcode=crop((int)floor(512*rcolor),0,511);
    if (!*arrow_flag){
      fprintf(fp,"2 1 0 0 0 %d %d 0 %d 0.000 0 0 -1 0 0 %d\n",/*colorcode*/colorcode+32,/*depth*/crop(depth,0,999),/*fill*/20,/*npoints*/4+1); fprintf(fp,"\t"); 
      fprintf(fp,"%d %d ",(int)floor(maxdia*(xoffset+(xord+xord2)*xside)),(int)floor(maxdia*(yoffset+(yord+yord2)*yside))); 
      fprintf(fp,"%d %d ",(int)floor(maxdia*(xoffset+(xord-xord2)*xside)),(int)floor(maxdia*(yoffset+(yord-yord2)*yside))); 
      fprintf(fp,"%d %d ",(int)floor(maxdia*(xoffset-(xord+xord2)*xside)),(int)floor(maxdia*(yoffset-(yord+yord2)*yside))); 
      fprintf(fp,"%d %d ",(int)floor(maxdia*(xoffset-(xord-xord2)*xside)),(int)floor(maxdia*(yoffset-(yord-yord2)*yside))); 
      fprintf(fp,"%d %d ",(int)floor(maxdia*(xoffset+(xord+xord2)*xside)),(int)floor(maxdia*(yoffset+(yord+yord2)*yside))); 
      fprintf(fp,"\n");}
    else if (*arrow_flag){
      fprintf(fp,"2 1 0 0 0 %d %d 0 %d 0.000 0 0 -1 0 0 %d\n",/*colorcode*/colorcode+32,/*depth*/crop(depth,0,999),/*fill*/20,/*npoints*/3+1); fprintf(fp,"\t"); 
      fprintf(fp,"%d %d ",(int)floor(maxdia*(xoffset+xord*xside)),(int)floor(maxdia*(yoffset+yord*yside))); 
      fprintf(fp,"%d %d ",(int)floor(maxdia*(xoffset-(xord+xord2)*xside)),(int)floor(maxdia*(yoffset-(yord+yord2)*yside))); 
      fprintf(fp,"%d %d ",(int)floor(maxdia*(xoffset-(xord-xord2)*xside)),(int)floor(maxdia*(yoffset-(yord-yord2)*yside))); 
      fprintf(fp,"%d %d ",(int)floor(maxdia*(xoffset+xord*xside)),(int)floor(maxdia*(yoffset+yord*yside))); 
      fprintf(fp,"\n");}}
}

void lattice3d_temp_dump(struct lattice3d *l3d,char *fgvn,void (*void_dump)(void *,FILE *,double,double,double,double,double,void *),void *void_parameters)
{
  int remove_flag=1,jpg_flag=(ON_MY_COMPUTER ? 1 : 0),eps_flag=0 /* (ON_MY_COMPUTER ? 1 : 0) */;
  char filename[1024],gs2[1024],filename2[1024],command[2048];
  FILE *fp=NULL;
  int nr1=0,nr2=0,nr3=0;
  struct spack3d *s=NULL;
  int k2x=0,k2y=0;
  double maxdia = 10000,side=1.0/sqrt(2.0);
  if (fgvn==NULL){ sprintf(gs2,"%srecord",GLOBAL_STRING_2);} else{ sprintf(gs2,"%srecord_%s",GLOBAL_STRING_2,fgvn);}
  sprintf(filename,"lattice3d_%s",gs2); sprintf(filename2,"./%s.fig",filename);
  if ((fp=fopen(filename2,"w"))==NULL){ printf(" %% warning! couldn't open %s in lattice3d_temp_dump\n",filename2); fp=stdout;}
  fprintf(fp,"%s",FIG_PREAMBLE); fprintf(fp,"%s",FIG_PREAMBLE_COLOR_7);
  for (nr3=0;nr3<l3d->input3max_scale;nr3++){ 
    k2x = nr3/(int)sqrt(l3d->input3max_scale); k2y = nr3%(int)sqrt(l3d->input3max_scale);
    for (nr1=0;nr1<l3d->input1max;nr1++){ for (nr2=0;nr2<l3d->input2max_scale;nr2++){ 
      s = lattice3dget(l3d,nr1,nr2,nr3);
      void_dump(s->temp,fp,maxdia,side/(double)l3d->input1max,side/(double)l3d->input2max_scale,s->x/(double)l3d->input1max+k2x,s->y/(double)l3d->input2max_scale+k2y,void_parameters);}}
    fprintf(fp,"4 0 0 %d 0 %d %d 0.0000 4 %d %d %d %d k%d\\001\n",/*depth*/999,/*font*/5,/*point*/12,/*textheight*/(int)(12*13.5),/*textwidth*/(int)(12*9*4),/*xpos*/(int)floor(maxdia*k2x),/*ypos*/(int)floor(maxdia*(k2y-0.1)),nr3);}
  if (fp!=stdout){ fclose(fp);fp=NULL;}
  if (jpg_flag){ sprintf(command,"fig2dev -Ljpeg -q %d ./%s.fig ./%s.jpg;",/*quality*/5,filename,filename); system(command);}
  if (eps_flag){ sprintf(command,"fig2dev -Leps ./%s.fig ./%s.eps;",filename,filename); system(command);}
  if (remove_flag){ sprintf(command,"rm ./%s.fig;",filename); system(command);}
}

void lattice3d_edge_v_dump(struct lattice3d *l3d,char *fgvn)
{
  int remove_flag=1,jpg_flag=(ON_MY_COMPUTER ? 1 : 0),eps_flag=0 /* (ON_MY_COMPUTER ? 1 : 0) */;
  char filename[1024],gs2[1024],filename2[1024],command[2048];
  FILE *fp=NULL;
  int nr1=0,nr2=0,nr3=0,tab=0;
  struct spack3d *s=NULL,*s2=NULL;
  double max=0,min=0;
  double rcolor=0,gcolor=0,bcolor=0;
  double xpos=0,ypos=0,xpos1=0,ypos1=0,xpos2=0,ypos2=0,d=0;
  int k2x=0,k2y=0;
  double maxdia = 10000;
  int colorcode=0;
  if (fgvn==NULL){ sprintf(gs2,"%srecord",GLOBAL_STRING_2);} else{ sprintf(gs2,"%srecord_%s",GLOBAL_STRING_2,fgvn);}
  sprintf(filename,"lattice3d_edge_v_%s",gs2); sprintf(filename2,"./%s.fig",filename);
  if ((fp=fopen(filename2,"w"))==NULL){ printf(" %% warning! couldn't open %s in lattice3d_edge_v_dump\n",filename2); fp=stdout;}
  fprintf(fp,"%s",FIG_PREAMBLE); fprintf(fp,"%s",FIG_PREAMBLE_COLOR_7);
  max = l3d->edge_v_mean+STD_VIEW*l3d->edge_v_stdev; min = l3d->edge_v_mean-STD_VIEW*l3d->edge_v_stdev;
  for (nr3=0;nr3<l3d->input3max_scale;nr3++){ 
    k2x = nr3/(int)sqrt(l3d->input3max_scale); k2y = nr3%(int)sqrt(l3d->input3max_scale);
    for (nr1=0;nr1<l3d->input1max;nr1++){ for (nr2=0;nr2<l3d->input2max_scale;nr2++){ 
      s = lattice3dget(l3d,nr1,nr2,nr3);
      if (s->adjacency_ra_s!=NULL && s->adjacency_ra_v!=NULL){
	for (tab=0;tab<s->L->length;tab++){
	  s2 = s->adjacency_ra_s[tab];
	  if (s2->k==s->k){
	    colorscale(0,s->adjacency_ra_v[tab],max,min,&rcolor,&gcolor,&bcolor); colorcode=crop((int)floor(512*rcolor),0,511);	    
	    xpos=(s->x+s2->x)/2.0; ypos=(s->y+s2->y)/2.0; 
	    xpos1 = s->y-s2->y; ypos1 = s2->x-s->x; d = sqrt(pow(xpos1,2)+pow(ypos1,2)); xpos2 = xpos1/d; ypos2 = ypos1/d;
	    if (l3d->tet_vs_cube==1){
	      xpos1 = xpos+xpos2/2/sqrt(3); ypos1 = ypos+ypos2/2/sqrt(3);
	      xpos2 = xpos-xpos2/2/sqrt(3); ypos2 = ypos-ypos2/2/sqrt(3);}
	    else /* if (l3d->tet_vs_cube==0) */{
	      xpos1 = xpos+xpos2/2; ypos1 = ypos+ypos2/2;
	      xpos2 = xpos-xpos2/2; ypos2 = ypos-ypos2/2;}
	    fprintf(fp,"2 1 0 0 0 %d %d 0 %d 0.000 0 0 -1 0 0 %d\n",/*colorcode*/colorcode+32,/*depth*/crop(nr3,0,999),/*fill*/20,/*npoints*/3+1); fprintf(fp,"\t"); fprintf(fp,"%d %d ",(int)floor(maxdia*(s->x/(double)l3d->input1max + k2x)),(int)floor(maxdia*(s->y/(double)l3d->input2max_scale + k2y))); fprintf(fp,"%d %d ",(int)floor(maxdia*(xpos1/(double)l3d->input1max + k2x)),(int)floor(maxdia*(ypos1/(double)l3d->input2max_scale + k2y))); fprintf(fp,"%d %d ",(int)floor(maxdia*(xpos2/(double)l3d->input1max + k2x)),(int)floor(maxdia*(ypos2/(double)l3d->input2max_scale + k2y))); fprintf(fp,"%d %d ",(int)floor(maxdia*(s->x/(double)l3d->input1max + k2x)),(int)floor(maxdia*(s->y/(double)l3d->input2max_scale + k2y))); fprintf(fp,"\n"); }}}}}
    fprintf(fp,"4 0 0 %d 0 %d %d 0.0000 4 %d %d %d %d k%d\\001\n",/*depth*/999,/*font*/5,/*point*/12,/*textheight*/(int)(12*13.5),/*textwidth*/(int)(12*9*4),/*xpos*/(int)floor(maxdia*k2x),/*ypos*/(int)floor(maxdia*(k2y-0.1)),nr3);}
  fprintf(fp,"4 0 0 %d 0 %d %d 0.0000 4 %d %d %d %d min%f_max%f\\001\n",/*depth*/999,/*font*/5,/*point*/12,/*textheight*/(int)(12*13.5),/*textwidth*/(int)(12*9*32),/*xpos*/(int)floor(maxdia*0),/*ypos*/(int)floor(maxdia*(0-0.2)),min,max);
  if (fp!=stdout){ fclose(fp);fp=NULL;}
  if (jpg_flag){ sprintf(command,"fig2dev -Ljpeg -q %d ./%s.fig ./%s.jpg;",/*quality*/5,filename,filename); system(command);}
  if (eps_flag){ sprintf(command,"fig2dev -Leps ./%s.fig ./%s.eps;",filename,filename); system(command);}
  if (remove_flag){ sprintf(command,"rm ./%s.fig;",filename); system(command);}
}

void lattice3dtfree(struct lattice3d *l3d,void (*free_function)(void *,void *),void *void_parameters)
{
  int nr1=0,nr2=0,nr3=0;
  for (nr1=0;nr1<l3d->input1max;nr1++){ for (nr2=0;nr2<l3d->input2max_scale;nr2++){ for (nr3=0;nr3<l3d->input3max_scale;nr3++){
    spack3dtfree(lattice3dget(l3d,nr1,nr2,nr3),free_function,void_parameters);}}}
  tfree(l3d->lattice);l3d->lattice=NULL;
  tfree(l3d);l3d=NULL;
}

struct ptree * binary_projection_s2p_make(double *sra,int length,int depth)
{
  /* actually constructs the appropriate ptree from the statetree array.
     in general we expect {s0->s1->s2...->sD} to be stored in statetree[s_{D} + s_{D-1}*slength + ... + s_{0}*pow(slength,D)]
     in general we expect {0->1->2->...->D} to be stored in ptree[D + (D-1)*length + ... + 0*pow(length,D)]
  */
  int verbose=0;
  int pdepth=0,nr=0,nd=0;
  int slength=(int)pow(2,length);
  double *s2p=NULL,*pra=NULL;
  struct ptree *p=NULL;
  struct llitem *l0=NULL,*l1=NULL;
  struct pnode *pn=NULL;
  int multiplier=1024;
  if (verbose){ printf(" %% [entering binary_projection_s2p_make] length %d depth %d\n",length,depth);}
  p = ptreemake(length,1,1,0,depth,1);
  for (pdepth=1;pdepth<=depth;pdepth++){
    s2p = binary_projection_s2p(0,length,depth,pdepth);
    pra = ra2ra_matrix_multiply(s2p,(int)pow(length,pdepth),(int)pow(slength,depth),0,sra,(int)pow(slength,depth),1,0);
    if (verbose){ printf(" %% ptree at level %d:\n",pdepth); raprintf(pra,"double",(int)pow(length,pdepth),1,"pra: ");}
    for (nr=0;nr<(int)pow(length,pdepth);nr++){
      if (pra[nr]!=0 || pra[nr]==0 /* waste space */){
	l0 = p->postree; pn=NULL;
	for (nd=0;nd<pdepth;nd++){
	  if ((l1=llitemaddorfind(0,l0,p->regionra[nr/(int)pow(length,nd)%length],&region2pnode_compare_label))==NULL){
	    l1=llitemaddorfind(1,l0,pnodemake(pn,p->regionra[nr/(int)pow(length,nd)%length],0,0),&pnode2pnode_compare_label);}
	  pn = (struct pnode *)l1->item;
	  l0 = pn->childllitem;}
	pn->weight = pra[nr]*multiplier; pn->relevance = pra[nr];
	l0 = p->pretree; pn=NULL;
	for (nd=pdepth-1;nd>=0;nd--){
	  if ((l1=llitemaddorfind(0,l0,p->regionra[nr/(int)pow(length,nd)%length],&region2pnode_compare_label))==NULL){
	    l1=llitemaddorfind(1,l0,pnodemake(pn,p->regionra[nr/(int)pow(length,nd)%length],0,0),&pnode2pnode_compare_label);}
	  pn = (struct pnode *)l1->item;
	  l0 = pn->childllitem;}
	pn->weight = pra[nr]*multiplier; pn->relevance = pra[nr];}}
    tfree(s2p);s2p=NULL;tfree(pra);pra=NULL;}
  if (verbose){ printf(" %% postree:\n");pnodeprintf(NULL,p->postree,-1,0);printf(" %% pretree:\n");pnodeprintf(NULL,p->pretree,-1,0);}
  return p;
}

double * binary_projection_s2s(int length,int depth,int sdepth)
{
  /* assumes that slength=(int)pow(2,length);
     assumes that a 2-statetree[ns2+ns*slength] = sra[ns2+ns*slength]*eig[ns + Imaxindex[0]*slength]/sum;
     thus statetree ordering is [s0->s0 s0->s1 s0->s2 ... s0->sN s1->s0 s1->s1 ... s1->sN s2->s0 ... sN->sN]
     in general we expect {s0->s1->s2...->sD} to be stored in statetree[s_{D} + s_{D-1}*slength + ... + s_{0}*pow(slength,D)]
     s2s is pow(slength,sdepth) by pow(slength,depth), */
  double *s2s=NULL;
  int slength = (int)pow(2,length);
  int nr1=0,nr2=0;
  s2s = (double *) tcalloc((int)pow(slength,sdepth)*(int)pow(slength,depth),sizeof(double));
  for (nr1=0;nr1<(int)pow(slength,depth);nr1++){
    nr2 = nr1/(int)pow(slength,depth-sdepth);
    s2s[nr2+nr1*(int)pow(slength,sdepth)] = 1;}
  return s2s;
}

double * binary_projection_s2p(int nonevent_vs_event,int length,int depth,int pdepth)
{
  /* assumes that slength=(int)pow(2,length);
     assumes that a 2-statetree[ns2+ns*slength] = sra[ns2+ns*slength]*eig[ns + Imaxindex[0]*slength]/sum;
     thus statetree ordering is [s0->s0 s0->s1 s0->s2 ... s0->sN s1->s0 s1->s1 ... s1->sN s2->s0 ... sN->sN]
     in general we expect {s0->s1->s2...->sD} to be stored in statetree[s_{D} + s_{D-1}*slength + ... + s_{0}*pow(slength,D)]
     s2p is pow(length,pdepth) by pow(slength,depth),
     in general we expect {0->1->2->...->D} to be stored in ptree[D + (D-1)*length + ... + 0*pow(length,D)] if nonevent_vs_event=0;
     and we expect {0->1->2->...->D-1->~D} to be stored in ptree[D + (D-1)*length + ... + 0*pow(length,D)] if nonevent_vs_event=1;
  */
  int verbose=0;
  double *s2p=NULL;
  int slength = (int)pow(2,length);
  int nr1=0,nr2=0,nr3=0,nr4=0,*ira=NULL;
  struct llist **Lra=NULL;
  int *itemp=NULL;
  struct litem **lra=NULL;
  int continue_flag=0,tab=0,shift_flag=0;
  if (verbose){ printf(" %% [entering binary_project_s2p] (%s),length %d,depth %d,pdepth %d\n",nonevent_vs_event?"nonevent":"event",length,depth,pdepth);}
  s2p = (double *) tcalloc((int)pow(length,pdepth)*(int)pow(slength,depth),sizeof(double));
  ira = (int *) tcalloc(pdepth,sizeof(int));
  Lra = (struct llist **) tcalloc(pdepth,sizeof(struct llist *));
  lra = (struct litem **) tcalloc(pdepth,sizeof(struct litem *));
  for (nr1=0;nr1<(int)pow(slength,depth);nr1++){
    nr2 = nr1/(int)pow(slength,depth-pdepth);
    if (verbose){ printf(" %% state-chain %d collapses to %d, which is composed of:\n",nr1,nr2);}
    for (nr3=0;nr3<pdepth;nr3++){
      ira[nr3] = (nr2/(int)pow(slength,nr3))%slength;
      if (verbose){ printf(" %% %% %d element is state %d\n",nr3,ira[nr3]);}
      if (Lra[nr3]!=NULL){ llisttfree(Lra[nr3]); Lra[nr3]=NULL;} Lra[nr3] = llistmake();
      for (nr4=0;nr4<length;nr4++){
	if ((ira[nr3]/(int)pow(2,nr4))%2 == (nr3==0 && nonevent_vs_event==1 ? 0 : 1)){ 
	  if (verbose){ printf(" %% %% %% neuron %d fired, adding to Lra[%d]\n",nr4,nr3);}
	  itemp = (int *) tcalloc(1,sizeof(int)); *itemp = nr4; litemadd(Lra[nr3],itemp);}}}
    if (verbose){ printf(" %% now stepping through firing neurons of state-chain\n");}
    for (nr3=0;nr3<pdepth;nr3++){ lra[nr3] = Lra[nr3]->first;}
    continue_flag=1; for(nr3=0;nr3<pdepth;nr3++){ continue_flag*=(lra[nr3]!=NULL);}
    while (continue_flag){
      if (verbose){ printf(" %% neuron indices: "); for (nr3=0;nr3<pdepth;nr3++){ printf(" %d",*(int *)lra[nr3]->item);} printf("... ");}
      tab = 0;
      for (nr3=0;nr3<pdepth;nr3++){ tab += *(int *)(lra[nr3]->item)*pow(length,nr3);}
      if (verbose){ printf(" correspond to tab %d\n",tab);}
      s2p[tab + nr1*(int)pow(length,pdepth)] = 1;
      nr3=0; shift_flag=1;
      do{ 
	lra[nr3]=lra[nr3]->child; 
	if (lra[nr3]==NULL){ lra[nr3]=Lra[nr3]->first; nr3+=1; shift_flag=1;}
	else /* if (lra[nr3]!=NULL) */{ shift_flag=0;}}
      while (shift_flag && nr3<pdepth);
      if (nr3>=pdepth){ continue_flag=0;}}}
  tfree(ira);ira=NULL;
  for (nr1=0;nr1<pdepth;nr1++){ llisttfree3(Lra[nr1]);} tfree(Lra); Lra=NULL;
  return s2p;
}

struct lattice3d * threetree_test(int nbins,int tet_vs_cube,double *synapse)
{
  /* assumes that synapse[s+n*nbins] is the connection strength from s to n */
  int verbose=0;
  int input1max=nbins,input2max=nbins,input3max=nbins;
  struct lattice3d *l3d=lattice3dmake(tet_vs_cube,input1max,input2max,input3max);
  struct spack3d *s=NULL;
  int length=3,slength=(int)pow(2,length);
  double *input=NULL,*statetree=NULL,*Gs1=NULL,*Gs2=NULL,*Gp1=NULL,*Gp2=NULL,*Cs1=NULL,*Cs2=NULL,*Cp1=NULL,*Cp2=NULL;
  double input_min = -1, input_max = 1;
  int nr1=0,nr2=0,nr3=0,tab=0,index=0;
  void ** vrara=NULL;
  int depth=0,sdepth=0,pdepth=0,rows=0,cols=0;
  double *s2s=NULL,*s2p=NULL;
  void **vra=NULL;
  int ilength_scale=0;
  double *normra=NULL,norm_max=0,norm_min=0;
  if (verbose){ printf(" %% [entering threetree_test]\n");}
  input=(double *)tcalloc(length,sizeof(double));
  if (verbose){ raprintf(synapse,"double",length,length,"%% synapse: ");}
  ilength_scale = l3d->input1max*l3d->input2max_scale*l3d->input3max_scale;
  normra = (double *) tcalloc(8*ilength_scale,sizeof(double));
  for (nr1=0;nr1<l3d->input1max;nr1++){ for (nr2=0;nr2<l3d->input2max_scale;nr2++){ for (nr3=0;nr3<l3d->input3max_scale;nr3++){
    s=lattice3dget(l3d,nr1,nr2,nr3);
    if (verbose>2){ spack3dprintf(s,NULL,NULL);}
    statetree = (double *)tcalloc(slength*slength,sizeof(double)); 
    Gs1 = (double *)tcalloc(length+1,sizeof(double));
    Gp1 = (double *)tcalloc(length+1,sizeof(double));
    Gs2 = (double *)tcalloc(length+1,sizeof(double));
    Gp2 = (double *)tcalloc(length+1,sizeof(double));
    Cs1 = (double *)tcalloc(length+1,sizeof(double));
    Cp1 = (double *)tcalloc(length+1,sizeof(double));
    Cs2 = (double *)tcalloc(length+1,sizeof(double));
    Cp2 = (double *)tcalloc(length+1,sizeof(double));
    input[0]=input_min+(input_max-input_min)*s->x/(double)l3d->input1max;
    input[1]=input_min+(input_max-input_min)*s->y/(double)l3d->input2max_scale;
    input[2]=input_min+(input_max-input_min)*s->z/(double)l3d->input3max_scale;
    ideal_2_statetree(length,input,synapse,statetree,Gs1,Gp1,Gs2,Gp2,Cs1,Cp1,Cs2,Cp2);
    tab = nr1 + nr2*l3d->input1max + nr3*l3d->input1max*l3d->input2max_scale;
    normra[tab+0*ilength_scale] = Gs1[length];
    normra[tab+1*ilength_scale] = Gp1[length];
    normra[tab+2*ilength_scale] = Gs2[length];
    normra[tab+3*ilength_scale] = Gp2[length];
    normra[tab+4*ilength_scale] = Cs1[length];
    normra[tab+5*ilength_scale] = Cp1[length];
    normra[tab+6*ilength_scale] = Cs2[length];
    normra[tab+7*ilength_scale] = Cp2[length];
    if (verbose>2){ raprintf(input,"double",1,length,"%% input: ");}
    if (verbose>2){ raprintf(statetree,"double",slength,slength,"%% statetree: ");}
    s->temp=(void *)((void **) tcalloc(9,sizeof(void *))); vra=s->temp;
    vra[0] = statetree;
    vra[1] = Gs1; vra[2] = Gp1; vra[3] = Gs2; vra[4] = Gp2;
    vra[5] = Cs1; vra[6] = Cp1; vra[7] = Cs2; vra[8] = Cp2;}}}
  verbose=1;
  if (verbose){ printf(" %% dumping Gs1\n");}
  stats("double",&(normra[0+0*ilength_scale]),ilength_scale,&norm_max,&norm_min,NULL,NULL);
  vra=(void **)tcalloc(5,sizeof(void *)); index=1; tab=0; vra[0]=&index; vra[1] = &length; vra[2]=&tab; vra[3]=&norm_max; vra[4]=&norm_min;
  lattice3d_temp_dump(l3d,"Gs1",&spack3d_temp_vector_dump,vra);
  tfree(vra);vra=NULL;
  if (verbose){ printf(" %% dumping Gp1\n");}
  stats("double",&(normra[0+1*ilength_scale]),ilength_scale,&norm_max,&norm_min,NULL,NULL);
  vra=(void **)tcalloc(5,sizeof(void *)); index=2; tab=0; vra[0]=&index; vra[1] = &length; vra[2]=&tab; vra[3]=&norm_max; vra[4]=&norm_min;
  lattice3d_temp_dump(l3d,"Gp1",&spack3d_temp_vector_dump,vra);
  tfree(vra);vra=NULL;
  if (verbose){ printf(" %% dumping Gs2\n");}
  stats("double",&(normra[0+2*ilength_scale]),ilength_scale,&norm_max,&norm_min,NULL,NULL);
  vra=(void **)tcalloc(5,sizeof(void *)); index=3; tab=0; vra[0]=&index; vra[1] = &length; vra[2]=&tab; vra[3]=&norm_max; vra[4]=&norm_min;
  lattice3d_temp_dump(l3d,"Gs2",&spack3d_temp_vector_dump,vra);
  tfree(vra);vra=NULL;
  if (verbose){ printf(" %% dumping Gp2\n");}
  stats("double",&(normra[0+3*ilength_scale]),ilength_scale,&norm_max,&norm_min,NULL,NULL);
  vra=(void **)tcalloc(5,sizeof(void *)); index=4; tab=0; vra[0]=&index; vra[1] = &length; vra[2]=&tab; vra[3]=&norm_max; vra[4]=&norm_min;
  lattice3d_temp_dump(l3d,"Gp2",&spack3d_temp_vector_dump,vra);
  tfree(vra);vra=NULL;
  if (verbose){ printf(" %% dumping Cs1\n");}
  stats("double",&(normra[0+4*ilength_scale]),ilength_scale,&norm_max,&norm_min,NULL,NULL);
  vra=(void **)tcalloc(5,sizeof(void *)); index=5; tab=1; vra[0]=&index; vra[1] = &length; vra[2]=&tab; vra[3]=&norm_max; vra[4]=&norm_min;
  lattice3d_temp_dump(l3d,"Cs1",&spack3d_temp_vector_dump,vra);
  tfree(vra);vra=NULL;
  if (verbose){ printf(" %% dumping Cp1\n");}
  stats("double",&(normra[0+5*ilength_scale]),ilength_scale,&norm_max,&norm_min,NULL,NULL);
  vra=(void **)tcalloc(5,sizeof(void *)); index=6; tab=1; vra[0]=&index; vra[1] = &length; vra[2]=&tab; vra[3]=&norm_max; vra[4]=&norm_min;
  lattice3d_temp_dump(l3d,"Cp1",&spack3d_temp_vector_dump,vra);
  tfree(vra);vra=NULL;
  if (verbose){ printf(" %% dumping Cs2\n");}
  stats("double",&(normra[0+6*ilength_scale]),ilength_scale,&norm_max,&norm_min,NULL,NULL);
  vra=(void **)tcalloc(5,sizeof(void *)); index=7; tab=1; vra[0]=&index; vra[1] = &length; vra[2]=&tab; vra[3]=&norm_max; vra[4]=&norm_min;
  lattice3d_temp_dump(l3d,"Cs2",&spack3d_temp_vector_dump,vra);
  tfree(vra);vra=NULL;
  if (verbose){ printf(" %% dumping Cp2\n");}
  stats("double",&(normra[0+7*ilength_scale]),ilength_scale,&norm_max,&norm_min,NULL,NULL);
  vra=(void **)tcalloc(5,sizeof(void *)); index=8; tab=1; vra[0]=&index; vra[1] = &length; vra[2]=&tab; vra[3]=&norm_max; vra[4]=&norm_min;
  lattice3d_temp_dump(l3d,"Cp2",&spack3d_temp_vector_dump,vra);
  tfree(vra);vra=NULL;
  if (verbose){ printf(" %% evaluating edges as 2-statetree \n");}
  depth=2; rows = pow(slength,depth); cols = pow(slength,depth);
  vrara = (void **) tcalloc(4,sizeof(void *));
  index=0; vrara[0] = &index; vrara[1] = &cols; vrara[2] = NULL; vrara[3] = NULL;
  lattice3d_edge_eval(l3d,&vra2vra_diff_sym,vrara); lattice3d_edge_v_stats(l3d);
  tfree(vrara);vrara=NULL;
  lattice3d_edge_v_dump(l3d,"s2");
  if (verbose){ printf(" %% evaluating edges as 1-statetree \n");}
  depth=2; sdepth=1; rows = pow(slength,sdepth); cols = pow(slength,depth);
  s2s = binary_projection_s2s(length,depth,sdepth); if (verbose){ raprintf(s2s,"double",rows,cols,"s2s: ");}
  vrara = (void **) tcalloc(4,sizeof(void *));
  index=0; vrara[0] = &index; vrara[1] = &cols; vrara[2] = s2s; vrara[3] = &rows;
  lattice3d_edge_eval(l3d,&vra2vra_diff_sym,vrara); lattice3d_edge_v_stats(l3d);
  tfree(vrara);vrara=NULL;tfree(s2s);s2s=NULL;
  lattice3d_edge_v_dump(l3d,"s1");
  if (verbose){ printf(" %% evaluating edges as 1-ptree \n");}
  depth=2; pdepth=1; rows = pow(length,pdepth); cols = pow(slength,depth);
  s2p = binary_projection_s2p(0,length,depth,pdepth); if (verbose){ raprintf(s2p,"double",rows,cols,"s2p: ");}
  vrara = (void **) tcalloc(4,sizeof(void *));
  index=0; vrara[0] = &index; vrara[1] = &cols; vrara[2] = s2p; vrara[3] = &rows;
  lattice3d_edge_eval(l3d,&vra2vra_diff_sym,vrara); lattice3d_edge_v_stats(l3d);
  tfree(vrara);vrara=NULL;tfree(s2p);s2p=NULL;
  lattice3d_edge_v_dump(l3d,"p1");
  if (verbose){ printf(" %% evaluating edges as 2-ptree \n");}
  depth=2; pdepth=2; rows = pow(length,pdepth); cols = pow(slength,depth);
  s2p = binary_projection_s2p(0,length,depth,pdepth); if (verbose){ raprintf(s2p,"double",rows,cols,"s2p: ");}
  vrara = (void **) tcalloc(4,sizeof(void *));
  index=0; vrara[0] = &index; vrara[1] = &cols; vrara[2] = s2p; vrara[3] = &rows;
  lattice3d_edge_eval(l3d,&vra2vra_diff_sym,vrara); lattice3d_edge_v_stats(l3d);
  tfree(vrara);vrara=NULL;tfree(s2p);s2p=NULL;
  lattice3d_edge_v_dump(l3d,"p2");
  tfree(input);input=NULL;
  tfree(normra);normra=NULL;
  return l3d;
}

void lattice3d_int2input_or_input2int(int direction_flag,int tet_vs_cube,int input1max,int input2max,int input3max,int *input1,int *input2,int *input3,double *output1,double *output2,double *output3)
{
  /* output given on lattice with nearest neighbors 1 unit away */
  int verbose=0;
  int input2max_scale=0,input3max_scale=0,npoints=0;
  double xpos=0,ypos=0,zpos=0;
  int xord=0,yord=0,zord=0;
  int fig_verbose=0;
  char filename1a[1024],filename1b[1024];
  FILE *fp=NULL;
  int nr1=0,nr2=0,nr3=0,tab=0,colorcode=0,ns=0;
  double rcolor=0,gcolor=0,bcolor=0;
  double maxdia = 20000;
  if (verbose){ printf(" %% [entering lattice3d_int2input_or_input2int] direction_flag %d tet_vs_cube %d\n",direction_flag,tet_vs_cube);}
  if (direction_flag==1){
    if (tet_vs_cube==1){
      xpos = *input1 - 0.5*(*input2%2);
      ypos = sqrt(3)/2.0*(*input2) + (*input3%2)*1.0/sqrt(3);
      zpos = sqrt(2.0/3.0)*(*input3);
      if (verbose){ printf(" input %d,%d,%d output %f,%f,%f\n",*input1,*input2,*input3,xpos,ypos,zpos);}
      if (output1!=NULL){ *output1 = xpos;} if (output2!=NULL){ *output2 = ypos;} if (output3!=NULL){ *output3 = zpos;}}
    else /* if (tet_vs_cube==0) */{
      xpos = *input1;
      ypos = *input2;
      zpos = *input3;
      if (verbose){ printf(" input %d,%d,%d output %f,%f,%f\n",*input1,*input2,*input3,xpos,ypos,zpos);}
      if (output1!=NULL){ *output1 = xpos;} if (output2!=NULL){ *output2 = ypos;} if (output3!=NULL){ *output3 = zpos;}}}
  else if (direction_flag==0){
    if (tet_vs_cube==1){
      zord = adi_round((*output3)/sqrt(2.0/3.0));
      yord = adi_round(((*output2)-(zord%2)*1.0/sqrt(3))/(sqrt(3)/2.0));
      xord = adi_round((*output1) + 0.5*(yord%2));
      if (verbose){ printf(" output %f,%f,%f input %d,%d,%d\n",*output1,*output2,*output3,xord,yord,zord);}
      if (input1!=NULL){ *input1=xord;} if (input2!=NULL){ *input2=yord;} if (input3!=NULL){ *input3=zord;}}
    else /* if (tet_vs_cube==0) */{
      zord = adi_round(*output3);
      yord = adi_round(*output2);
      xord = adi_round(*output1);
      if (verbose){ printf(" output %f,%f,%f input %d,%d,%d\n",*output1,*output2,*output3,xord,yord,zord);}
      if (input1!=NULL){ *input1=xord;} if (input2!=NULL){ *input2=yord;} if (input3!=NULL){ *input3=zord;}}}
  else /* if (direction_flag!=1 && direction+flag!=0) */{
    if (verbose){ printf(" %% testing lattice3d_int2input_or_input2int\n");}
    sprintf(filename1a,"./lattice3d_int2input_or_input2int_test"); sprintf(filename1b,"%s.fig",filename1a);
    if ((fp=fopen(filename1b,"w"))==NULL){ printf(" %% Can't create %s in lattice3d_int2input_or_input2int\n",filename1b); fp=stdout;}
    fprintf(fp,"%s",FIG_PREAMBLE); fprintf(fp,"%s",FIG_PREAMBLE_COLOR_7);
    if (tet_vs_cube==1){ input2max_scale = (int)(input2max/(sqrt(3)/2.0)); input3max_scale = (int)(input3max/sqrt(2.0/3.0)); npoints=6;}
    else /* if (tet_vs_cube==0) */{ input2max_scale = input2max; input3max_scale = input3max; npoints=4;}
    for (nr1=0;nr1<input1max;nr1++){ for (nr2=0;nr2<input2max_scale;nr2++){ for (nr3=0;nr3<input3max_scale;nr3++){
      tab = nr1+nr2*input1max+nr3*input1max*input2max_scale;
      lattice3d_int2input_or_input2int(1,tet_vs_cube,input1max,input2max,input3max,&nr1,&nr2,&nr3,&xpos,&ypos,&zpos);
      lattice3d_int2input_or_input2int(0,tet_vs_cube,input1max,input2max,input3max,&xord,&yord,&zord,&xpos,&ypos,&zpos);
      if (nr1!=xord){ printf(" %% warning! 1st coordinate %d improperly inverted to %d\n",nr1,xord);}
      if (nr2!=yord){ printf(" %% warning! 2nd coordinate %d improperly inverted to %d\n",nr2,yord);}
      if (nr3!=zord){ printf(" %% warning! 3rd coordinate %d improperly inverted to %d\n",nr3,zord);}
      colorscale(0,tab,input1max*input2max_scale*input3max_scale,0,&rcolor,&gcolor,&bcolor);colorcode=crop((int)floor(512*rcolor),0,511);
      fprintf(fp,"2 1 0 0 0 %d %d 0 %d 0.000 0 0 -1 0 0 %d\n",/*colorcode*/colorcode+32,/*depth*/crop((int)nr3+1,0,999),/*fill*/20,/*npoints*/npoints+1); fprintf(fp,"\t"); for (ns=0;ns<npoints+1;ns++){ fprintf(fp,"%d %d ",(int)floor(maxdia*(ypos+cos(2*PI*(ns+(npoints-2)/4.0)/(double)npoints)/sqrt(npoints/2))/(double)input2max),(int)floor(maxdia*(xpos+sin(2*PI*(ns+(npoints-2)/4.0)/(double)npoints)/sqrt(npoints/2))/(double)input1max));} fprintf(fp,"\n");if (fig_verbose){ fprintf(fp,"4 0 0 %d 0 %d %d 0.0000 4 %d %d %d %d i%dj%dk%d\\001\n",/*depth*/crop((int)nr3+1,0,999),/*font*/5,/*point*/5,/*textheight*/(int)(5*13.5),/*textwidth*/(int)(5*9*7),/*xpos*/(int)floor(maxdia*ypos/(double)input2max),/*ypos*/(int)floor(maxdia*xpos/(double)input1max),nr1,nr2,nr3);}}}}
    if (fp!=stdout){ fclose(fp);}}
}

double * bp_projection(int length,int depth)
{
  /* assumes depth==2,
     generates standard BP projection from 2-ptree, where {n0->n1} is stored in ptree[n1 + n0*length], 
     and same ordering is used for synapses */
  double *ra=NULL;
  int nr1=0,nr2=0,nr3=0,nr4=0;
  ra = (double *) tcalloc((int)pow(length,2)*(int)pow(length,2),sizeof(double));
  for (nr1=0;nr1<length;nr1++){ for (nr2=0;nr2<length;nr2++){ for (nr3=0;nr3<length;nr3++){ for (nr4=0;nr4<length;nr4++){
    if (nr1!=nr2){
      if (nr4==nr2 && nr3==nr1){ ra[nr2+nr1*length + (nr4+nr3*length)*(int)pow(length,2)] = +1;}
      else if (nr4==nr1 && nr3==nr2){ ra[nr2+nr1*length + (nr4+nr3*length)*(int)pow(length,2)] = -1;}}}}}}
  return ra;
}

double * ds_projection(int rows)
{
  /* generates a unitary matrix such that Q(:,1)=[1 1 ... 1]/sqrt(rows) */
  int verbose=0;
  double *ra=NULL,*ra2=NULL;
  int nr=0,cols=rows,nc=0,nc2=0;
  double d=0;
  ra = (double *) tcalloc(rows*rows,sizeof(double));
  for (nr=0;nr<rows;nr++){ ra[nr]=1.0; ra[nr+nr*rows]=1.0;}
  for (nc=0;nc<cols;nc++){ 
    d=ra2ra_dot(&(ra[0+nc*rows]),&(ra[0+nc*rows]),rows); if (d>0){ for (nr=0;nr<rows;nr++){ ra[nr+nc*rows]/=sqrt(d);}}
    for (nc2=nc+1;nc2<cols;nc2++){
      d=ra2ra_dot(&(ra[0+nc*rows]),&(ra[0+nc2*rows]),rows); for (nr=0;nr<rows;nr++){ ra[nr+nc2*rows]-=d*ra[nr+nc*rows];}}}
  if (verbose){ 
    raprintf(ra,"double",rows,rows,"ds_projection: ");
    ra2 = ra2ra_matrix_multiply(ra,rows,cols,1,ra,rows,cols,0);
    raprintf(ra2,"double",rows,rows,"dsp'*dsp: ");
    tfree(ra2);ra2=NULL;}
  return ra;
}

void ideal_2_statetree_old(int length,double *input,double *synapse,double *statetree_2,double *granulation_direction_input,double *bp_direction_input)
{
  /* assumes that *input is length long, *synapse is length*length long, *statetree_2 is (int)pow((int)pow(2,length),2) long 
     in order for the 1-statetree to be an eigenvector of statetree, we assume that statetree ordering is:
     [ s0->s0 s1->s0 ... sN->s0 ]
     [ s0->s1 s1->s1 ... sN->s1 ]
     [  ...               ...   ]
     [ s0->sN s1->sN ... sN->sN ]
     stored in row dominant form: 
     statetree ordering is [s0->s0 s0->s1 s0->s2 ... s0->sN s1->s0 s1->s1 ... s1->sN s2->s0 ... sN->sN]
     in general we expect {s0->s1->s2...->sD} to be stored in statetree[s_{D} + s_{D-1}*slength + ... + s_{0}*pow(slength,D)]
     ane we expect {n0->n1->...->nD} to be stored in ptree[n_{D} + n_{D-1}*length + ... + n_{0}*pow(length,D)]
     however here we assume that depth D=2 
     we also assume that the relevant probability of state s producing event A is simply 0.5+0.5*erf(\nu_{A}+\sum_{B\in s} syn_{B,A})
     hence the partial derivative with respect to an argument is simply 1.0/sqrt(PI)*exp(-pow(\nu_{A}+\sum_{B\in s} syn_{B,A},2)) */
  int verbose=0;
  int nr=0,ns=0,nr2=0,ns2=0;
  int slength = (int)pow(2,length);
  double *pra=NULL,*sra=NULL,*sraprime_input=NULL,*sraprime_synapse=NULL,*statetree_1=NULL;
  double *eiwr=NULL,*eiwi=NULL,*eig=NULL;
  int *Imaxindex=NULL,Inmax=0,*Iminindex=NULL,Inmin=0;
  double tolerance=0.000001,sum=0;
  char text[32];
  struct matrix *ms1=mmake(),*ms2=mmake(),*mdsp=mmake(),*mdspt=mmake(),*mqtlq=mmake();
  struct matrix *mqtdldi=mmake(),*mqtdlds=mmake(),*mtmp=mmake(),*mtmp2=mmake(),*mrhsi=mmake(),*mrhss=mmake();
  struct matrix *msv=mmake(),*msvr=mmake(),*msvl=mmake();
  if (verbose){ printf(" %% [entering ideal_2_statetree] length %d\n",length);}
  if (verbose){ raprintf(input,"double",1,length," %% input: ");}
  if (verbose){ raprintf(synapse,"double",length,length," %% synapse_n->n: ");}
  pra = (double *) tcalloc(slength*length,sizeof(double));
  for (nr=0;nr<length;nr++){ 
    if (verbose>1){ printf(" %% neuron %d\n",nr);}
    for (ns=0;ns<slength;ns++){
      if (verbose>1){ printf(" %% %% state %d,",ns);}
      pra[ns+nr*slength] = input[nr];
      if (verbose>1){ printf(" starting with input %0.2f,",input[nr]);}
      for (nr2=0;nr2<length;nr2++){ 
	if ((ns/(int)pow(2,nr2))%2==1){ 
	  if (verbose>1){ printf(" neuron %d fired, adding %0.2f",nr2,synapse[nr2+nr*length]);}
	  pra[ns+nr*slength] += synapse[nr2+nr*length];}}
      if (verbose>1){ printf(" ending with input %0.3f\n",pra[ns+nr*slength]);}}}
  sra = (double *) tcalloc(slength*slength,sizeof(double));
  sraprime_input = (double *) tcalloc(slength*slength*length,sizeof(double));
  sraprime_synapse = (double *) tcalloc(slength*slength*(int)pow(length,2),sizeof(double));
  statetree_1 = (double *) tcalloc(slength,sizeof(double));
  for (ns=0;ns<slength;ns++){ for (ns2=0;ns2<slength;ns2++){
    if (verbose>1){ printf(" %% considering link from state %d to state %d\n",ns,ns2);}
    sra[ns2+ns*slength]=1;
    for (nr=0;nr<length;nr++){
      if ((ns2/(int)pow(2,nr))%2==1){
	if (verbose>1){ printf(" %% %% neuron %d fired in state %d, factor %f\n",nr,ns2,0.5+0.5*erf(pra[ns+nr*slength]));}
	sra[ns2+ns*slength] *= 0.5+0.5*erf(pra[ns+nr*slength]);}
      else /* if ((ns2/(int)pow(2,nr))%2!=1) */{ 
	if (verbose>1){ printf(" %% %% neuron %d did not fire in state %d, factor (1-%f)\n",nr,ns2,0.5-0.5*(pra[ns+nr*slength]));}
	sra[ns2+ns*slength] *= 0.5-0.5*erf(pra[ns+nr*slength]);}}
    if (verbose>1){ printf(" %% computed sra[%d+%d*slength]=%f\n",ns2,ns,sra[ns2+ns*slength]);}
    for (nr=0;nr<length;nr++){
      if ((ns2/(int)pow(2,nr))%2==1){ sraprime_input[ns2+ns*slength+nr*slength*slength]=sra[ns2+ns*slength]*(+1.0)*exp(-pow(pra[ns+nr*slength],2))/sqrt(PI)/(0.5+0.5*erf(pra[ns+nr*slength]));}
      else /* if ((ns2/(int)pow(2,nr))%2!=1) */{ sraprime_input[ns2+ns*slength+nr*slength*slength]=sra[ns2+ns*slength]*(-1.0)*exp(-pow(pra[ns+nr*slength],2))/sqrt(PI)/(0.5-0.5*erf(pra[ns+nr*slength]));}
      if (verbose>1){ printf(" %% computed sraprime_input[%d+%d*slength+%d*slength*slength]=%f\n",ns2,ns,nr,sraprime_input[ns2+ns*slength+nr*slength*slength]);}}
    for (nr=0;nr<length;nr++){ for (nr2=0;nr2<length;nr2++){
      if (verbose>1){ printf("ns %d ns2 %d nr%d nr2%d m1 %d m2 %d\n",ns,ns2,nr,nr2,(ns/(int)pow(2,nr))%2==1,(ns2/(int)pow(2,nr2))%2==1);}
      if (((ns/(int)pow(2,nr))%2==1) && ((ns2/(int)pow(2,nr2))%2==1)){ sraprime_synapse[ns2+ns*slength+(nr2+nr*length)*slength*slength]=sra[ns2+ns*slength]*(+1.0)*exp(-pow(pra[ns+nr2*slength],2))/sqrt(PI)/(0.5+0.5*erf(pra[ns+nr2*slength]));}
      else if (((ns/(int)pow(2,nr))%2==1) && ((ns2/(int)pow(2,nr2))%2!=1)){ sraprime_synapse[ns2+ns*slength+(nr2+nr*length)*slength*slength]=sra[ns2+ns*slength]*(-1.0)*exp(-pow(pra[ns+nr2*slength],2))/sqrt(PI)/(0.5-0.5*erf(pra[ns+nr2*slength]));}
      else /* if (((ns/(int)pow(2,nr))%2!=1)) */{ sraprime_synapse[ns2+ns*slength+(nr2+nr*length)*slength*slength]=0;}
      if (verbose>1){ printf(" %% computed sraprime_synapse[%d+%d*slength+(%d+%d*length)*slength*slength]=%f\n",ns2,ns,nr2,nr,sraprime_synapse[ns2+ns*slength+(nr2+nr*length)*slength*slength]);}}}
    if (verbose>1){ printf(" %% finished linking state %d to state %d\n",ns,ns2);}}}
  if (verbose){ 
    raprintf(pra,"double",slength,length,"pra[s->n]: ");
    raprintf(sra,"double",slength,slength,"sra[s->s]: ");
    for (nr=0;nr<length;nr++){
      sprintf(text,"sraprime_input_n%d",nr);
      raprintf(&(sraprime_input[0+nr*slength*slength]),"double",slength,slength,text);}
    for (nr=0;nr<length;nr++){ for (nr2=0;nr2<length;nr2++){
      sprintf(text,"sraprime_synapse_n%d->n%d",nr,nr2);
      raprintf(&(sraprime_synapse[0+(nr2+nr*length)*slength*slength]),"double",slength,slength,text);}}}
  eiwr = (double *) tcalloc(slength,sizeof(double));
  eiwi = (double *) tcalloc(slength,sizeof(double));
  eig = (double *) tcalloc(slength*slength,sizeof(double));
  ra2eig(sra,slength,eiwr,eiwi,eig);
  maxmindex("double",eiwr,slength,&Imaxindex,&Inmax,&Iminindex,&Inmin,tolerance);
  if (verbose){ printf(" %% %d largest eigenvalues, first of which is %f\n",Inmax,eiwr[Imaxindex[0]]);}
  stats("double",&(eig[0+Imaxindex[0]*slength]),slength,NULL,NULL,&sum,NULL); sum*=slength;
  for (ns=0;ns<slength;ns++){ statetree_1[ns] = eig[ns+Imaxindex[0]*slength]/sum;}
  if (verbose){ raprintf(statetree_1,"double",1,slength,"statetree_1: ");}
  if (statetree_2!=NULL){
    for (ns=0;ns<slength;ns++){ for (ns2=0;ns2<slength;ns2++){
      statetree_2[ns2+ns*slength] = sra[ns2+ns*slength]*statetree_1[ns];}}
    if (verbose){ raprintf(statetree_2,"double",slength,slength," %% statetree_2: ");}}

  //  verbose=1;
  mtfree(ms1);minit0(ms1,slength,1);memcpy(ms1->mtrx,statetree_1,slength*sizeof(double));
  mtfree(mdsp);mdsp->rows=slength;mdsp->cols=mdsp->rows;mdsp->mtrx=ds_projection(mdsp->rows); mplugout(mdsp,0,slength-1,1,slength-1,mdsp);
  if (verbose){ raprintf(mdsp->mtrx,"double",mdsp->rows,mdsp->cols,"mdsp: ");}
  mtrans(mdsp,mdspt);
  minit1(mtmp,slength,slength);minit0(mtmp2,slength,slength);memcpy(mtmp2->mtrx,sra,slength*slength*sizeof(double));
  msubtm(mtmp,mtmp2,mqtlq);
  if (verbose){ raprintf(mqtlq->mtrx,"double",mqtlq->rows,mqtlq->cols,"mqtlq: ");}
  mtimesm(mqtlq,mdsp,mtmp);mtimesm(mdspt,mtmp,mqtlq);minv(mqtlq,mqtlq);
  if (verbose){ raprintf(mqtlq->mtrx,"double",mqtlq->rows,mqtlq->cols,"mqtlq: ");}
  minit0(mqtdldi,slength,length);
  for (nr=0;nr<length;nr++){ 
    minit0(mtmp,slength,slength);
    memcpy(mtmp->mtrx,&(sraprime_input[0+0*slength+nr*slength*slength]),slength*slength*sizeof(double));
    mtimesm(mtmp,ms1,mtmp);
    mplugin(mtmp,0,slength-1,0,0,mqtdldi,0,nr);}
  mtimesm(mdspt,mqtdldi,mqtdldi);
  if (verbose){ raprintf(mqtdldi->mtrx,"double",mqtdldi->rows,mqtdldi->cols,"mqtdldi: ");}
  mtimesm(mqtlq,mqtdldi,mrhsi);
  if (verbose){ raprintf(mrhsi->mtrx,"double",mrhsi->rows,mrhsi->cols,"mrhsi: ");}
  minit0(msv,minimum(mrhsi->rows,mrhsi->cols),1);
  minit0(msvl,mrhsi->rows,mrhsi->rows);
  minit0(msvr,mrhsi->cols,mrhsi->cols);
  ra2svd(mrhsi->mtrx,mrhsi->rows,mrhsi->cols,msv->mtrx,msvl->mtrx,msvr->mtrx);
  for (nr=0;nr<length;nr++){ granulation_direction_input[nr] = gentry(msvr,0,nr);} granulation_direction_input[length]=gentry(msv,0,0);
  if (verbose){ raprintf(granulation_direction_input,"double",1,length,"granulation_direction_input: "); printf("expansion %f\n",granulation_direction_input[length]);}
  
  minit0(ms2,(int)pow(slength,2),1); memcpy(ms2->mtrx,statetree_2,slength*slength*sizeof(double));
  if (verbose){ raprintf(ms2->mtrx,"double",ms2->rows,ms2->cols,"ms2: ");}
  minit0(mqtdlds,slength,(int)pow(length,2));
  for (nr=0;nr<length;nr++){ for (nr2=0;nr2<length;nr2++){
    minit0(mtmp,slength,slength);
    memcpy(mtmp->mtrx,&(sraprime_synapse[0+0*slength+(nr2+nr*length)*slength*slength]),slength*slength*sizeof(double));
    mtimesm(mtmp,ms1,mtmp);
    mplugin(mtmp,0,slength-1,0,0,mqtdlds,0,nr2+nr*length);}}
  mtimesm(mdspt,mqtdlds,mqtdlds);
  if (verbose){ raprintf(mqtdlds->mtrx,"double",mqtdlds->rows,mqtdlds->cols,"mqtdlds: ");}
  mtimesm(mqtlq,mqtdlds,mrhss);
  if (verbose){ raprintf(mrhss->mtrx,"double",mrhss->rows,mrhss->cols,"mrhss: ");}
  mtfree(mtmp2);mtmp2->rows=(int)pow(length,2);mtmp2->cols=(int)pow(slength,2);mtmp2->mtrx=binary_projection_s2p(0,length,2,2);
  if (verbose){ raprintf(mtmp2->mtrx,"double",mtmp2->rows,mtmp2->cols,"s2p: ");}
  mtfree(mtmp);mtmp->rows=(int)pow(length,2);mtmp->cols=(int)pow(length,2);mtmp->mtrx=bp_projection(length,2);
  if (verbose){ raprintf(mtmp->mtrx,"double",mtmp->rows,mtmp->cols,"p2syn: ");}
  mtimesm(mtmp,mtmp2,mtmp2);mtimesm(mtmp2,ms2,mtmp);
  if (verbose){ raprintf(mtmp->mtrx,"double",mtmp->rows,mtmp->cols,"p2syn*s2p*s2: ");}
  mtimesm(mrhss,mtmp,mtmp2);  
  if (verbose){ raprintf(mtmp2->mtrx,"double",mtmp2->rows,mtmp2->cols,"mrhss*p2syn*s2p*s2: ");}
  minit0(mtmp,slength-1,1);
  ra2lss(mrhsi->mtrx,slength-1,length,mtmp2->mtrx,1,mtmp->mtrx);
  for (nr=0;nr<length;nr++){ bp_direction_input[nr] = gentry(mtmp,nr,0);} bp_direction_input[length] = 1;
  ratimesequals(bp_direction_input,length,1.0/ra_norm(bp_direction_input,length));
  if (verbose){ raprintf(bp_direction_input,"double",1,length,"bp_direction_input: ");}
  //exit(0);
  
  mtfree(ms1);tfree(ms1); mtfree(ms2);tfree(ms2); mtfree(mdsp);tfree(mdsp); mtfree(mdspt);tfree(mdspt); mtfree(mqtlq);tfree(mqtlq); 
  mtfree(mqtdldi);tfree(mqtdldi); mtfree(mqtdlds);tfree(mqtdlds);
  mtfree(mtmp);tfree(mtmp); mtfree(mtmp2);tfree(mtmp2);
  mtfree(mrhsi);tfree(mrhsi); mtfree(mrhss);tfree(mrhss);
  mtfree(msv);tfree(msv); mtfree(msvl);tfree(msvl); mtfree(msvr);tfree(msvr);
  tfree(Imaxindex);Imaxindex=NULL; tfree(Iminindex);Iminindex=NULL;
  tfree(pra);pra=NULL; tfree(sra);sra=NULL; tfree(sraprime_input);sraprime_input=NULL; tfree(sraprime_synapse);sraprime_synapse=NULL;
  tfree(eiwr);eiwr=NULL; tfree(eiwi);eiwi=NULL; tfree(eig);eig=NULL;
}

void ideal_2_statetree(int length,double *input,double *synapse,double *statetree_2,double *Gs1,double *Gp1,double *Gs2,double *Gp2,double *Cs1,double *Cp1,double *Cs2,double *Cp2)
{
  /* assumes that *input is length long, *synapse is length*length long, *statetree_2 is (int)pow((int)pow(2,length),2) long 
     in order for the 1-statetree to be an eigenvector of statetree, we assume that statetree ordering is:
     [ s0->s0 s1->s0 ... sN->s0 ]
     [ s0->s1 s1->s1 ... sN->s1 ]
     [  ...               ...   ]
     [ s0->sN s1->sN ... sN->sN ]
     stored in row dominant form: 
     statetree ordering is [s0->s0 s0->s1 s0->s2 ... s0->sN s1->s0 s1->s1 ... s1->sN s2->s0 ... sN->sN]
     in general we expect {s0->s1->s2...->sD} to be stored in statetree[s_{D} + s_{D-1}*slength + ... + s_{0}*pow(slength,D)]
     ane we expect {n0->n1->...->nD} to be stored in ptree[n_{D} + n_{D-1}*length + ... + n_{0}*pow(length,D)]
     however here we assume that depth D=2 
     we also assume that the relevant probability of state s producing event A is simply 0.5+0.5*erf(\nu_{A}+\sum_{B\in s} syn_{B,A})
     hence the partial derivative with respect to an argument is simply 1.0/sqrt(PI)*exp(-pow(\nu_{A}+\sum_{B\in s} syn_{B,A},2)) */
  int verbose=0;
  int nr=0,ns=0,nr2=0,ns2=0,nr3=0,ns3=0;
  int slength = (int)pow(2,length);
  double *pra=NULL,*sra=NULL,*sraprime_input=NULL,*sraprime_synapse=NULL,*sraprimeprime_input_synapse=NULL,*statetree_1=NULL;
  double *eiwr=NULL,*eiwi=NULL,*eig=NULL;
  int *Imaxindex=NULL,Inmax=0,*Iminindex=NULL,Inmin=0;
  double tolerance=0.000001,sum=0;
  char text[32];
  struct matrix *mtmp1=mmake(),*mtmp2=mmake(),*mtmp3=mmake();
  struct matrix *ms1=mmake(),*ms2=mmake(),*ms2p2=mmake(),*ms1p1=mmake();
  struct matrix *mL=mmake(),*mLh=mmake();
  struct matrix *mQ=mmake(),*mQt=mmake();
  struct matrix *mM=mmake();
  struct matrix *mdiLs1=mmake(),*mdis1=mmake(),*mdiLhs1=mmake(),*mdis2=mmake(),*mdiLhdss1=mmake();
  struct matrix *mbp=mmake(),*mbps2=mmake();
  struct matrix *mdsLs1=mmake(),*mdss1=mmake(),*mdsLhs1=mmake(),*mdss2=mmake();
  struct matrix *mdiMdsLs1bps2=mmake(),*mMdidsLs1bps2=mmake(),*mMdsLdis1bps2=mmake(),*mMdsLs1bpdis2=mmake();
  struct matrix *mC=mmake(),*mC2=mmake();  
  struct matrix *msv=mmake(),*msvl=mmake(),*msvr=mmake();
  struct matrix *meiwr=mmake(),*meiwi=mmake(),*meig=mmake();
/*   double threshold=0.000001; */
  if (verbose){ printf(" %% [entering ideal_2_statetree] length %d\n",length);}
  if (verbose){ raprintf(input,"double",1,length," %% input: ");}
  if (verbose){ raprintf(synapse,"double",length,length," %% synapse_n->n: ");}
  pra = (double *) tcalloc(slength*length,sizeof(double));
  for (nr=0;nr<length;nr++){ 
    if (verbose>1){ printf(" %% neuron %d\n",nr);}
    for (ns=0;ns<slength;ns++){
      if (verbose>1){ printf(" %% %% state %d,",ns);}
      pra[ns+nr*slength] = input[nr];
      if (verbose>1){ printf(" starting with input %0.2f,",input[nr]);}
      for (nr2=0;nr2<length;nr2++){ 
	if ((ns/(int)pow(2,nr2))%2==1){ 
	  if (verbose>1){ printf(" neuron %d fired, adding %0.2f",nr2,synapse[nr2+nr*length]);}
	  pra[ns+nr*slength] += synapse[nr2+nr*length];}}
      if (verbose>1){ printf(" ending with input %0.3f\n",pra[ns+nr*slength]);}}}
  sra = (double *) tcalloc(slength*slength,sizeof(double));
  sraprime_input = (double *) tcalloc(slength*slength*length,sizeof(double));
  sraprime_synapse = (double *) tcalloc(slength*slength*(int)pow(length,2),sizeof(double));
  sraprimeprime_input_synapse = (double *) tcalloc(slength*slength*(int)pow(length,2)*length,sizeof(double));
  statetree_1 = (double *) tcalloc(slength,sizeof(double));
  for (ns=0;ns<slength;ns++){ for (ns2=0;ns2<slength;ns2++){
    if (verbose>1){ printf(" %% considering link from state %d to state %d\n",ns,ns2);}
    sra[ns2+ns*slength]=1;
    for (nr=0;nr<length;nr++){
      if ((ns2/(int)pow(2,nr))%2==1){
	if (verbose>1){ printf(" %% %% neuron %d fired in state %d, factor %f\n",nr,ns2,0.5+0.5*erf(pra[ns+nr*slength]));}
	sra[ns2+ns*slength] *= 0.5+0.5*erf(pra[ns+nr*slength]);}
      else /* if ((ns2/(int)pow(2,nr))%2!=1) */{ 
	if (verbose>1){ printf(" %% %% neuron %d did not fire in state %d, factor (1-%f)\n",nr,ns2,0.5-0.5*(pra[ns+nr*slength]));}
	sra[ns2+ns*slength] *= 0.5-0.5*erf(pra[ns+nr*slength]);}}
    if (verbose>1){ printf(" %% computed sra[%d+%d*slength]=%f\n",ns2,ns,sra[ns2+ns*slength]);}
    for (nr=0;nr<length;nr++){
      if ((ns2/(int)pow(2,nr))%2==1){ sraprime_input[ns2+ns*slength+nr*slength*slength]=sra[ns2+ns*slength]*(+1.0)*exp(-pow(pra[ns+nr*slength],2))/sqrt(PI)/(0.5+0.5*erf(pra[ns+nr*slength]));}
      else /* if ((ns2/(int)pow(2,nr))%2!=1) */{ sraprime_input[ns2+ns*slength+nr*slength*slength]=sra[ns2+ns*slength]*(-1.0)*exp(-pow(pra[ns+nr*slength],2))/sqrt(PI)/(0.5-0.5*erf(pra[ns+nr*slength]));}
      if (verbose>1){ printf(" %% computed sraprime_input[%d+%d*slength+%d*slength*slength]=%f\n",ns2,ns,nr,sraprime_input[ns2+ns*slength+nr*slength*slength]);}}
    for (nr=0;nr<length;nr++){ for (nr2=0;nr2<length;nr2++){
      if (verbose>1){ printf("ns %d ns2 %d nr%d nr2%d m1 %d m2 %d\n",ns,ns2,nr,nr2,(ns/(int)pow(2,nr))%2==1,(ns2/(int)pow(2,nr2))%2==1);}
      if (((ns/(int)pow(2,nr))%2==1) && ((ns2/(int)pow(2,nr2))%2==1)){ sraprime_synapse[ns2+ns*slength+(nr2+nr*length)*slength*slength]=sra[ns2+ns*slength]*(+1.0)*exp(-pow(pra[ns+nr2*slength],2))/sqrt(PI)/(0.5+0.5*erf(pra[ns+nr2*slength]));}
      else if (((ns/(int)pow(2,nr))%2==1) && ((ns2/(int)pow(2,nr2))%2!=1)){ sraprime_synapse[ns2+ns*slength+(nr2+nr*length)*slength*slength]=sra[ns2+ns*slength]*(-1.0)*exp(-pow(pra[ns+nr2*slength],2))/sqrt(PI)/(0.5-0.5*erf(pra[ns+nr2*slength]));}
      else /* if (((ns/(int)pow(2,nr))%2!=1)) */{ sraprime_synapse[ns2+ns*slength+(nr2+nr*length)*slength*slength]=0;}
      if (verbose>1){ printf(" %% computed sraprime_synapse[%d+%d*slength+(%d+%d*length)*slength*slength]=%f\n",ns2,ns,nr2,nr,sraprime_synapse[ns2+ns*slength+(nr2+nr*length)*slength*slength]);}}}
    for (nr=0;nr<length;nr++){ for (nr2=0;nr2<length;nr2++){ for (nr3=0;nr3<length;nr3++){
      if (((ns/(int)pow(2,nr))%2==1) && ((ns2/(int)pow(2,nr2))%2==1) && nr2!=nr3 && ((ns2/(int)pow(2,nr3))%2==1)){ sraprimeprime_input_synapse[ns2+ns*slength+(nr2+nr*length + nr3*length*length)*slength*slength]=sra[ns2+ns*slength]*(+1.0)*exp(-pow(pra[ns+nr2*slength],2))/sqrt(PI)/(0.5+0.5*erf(pra[ns+nr2*slength]))*(+1.0)*exp(-pow(pra[ns+nr3*slength],2))/sqrt(PI)/(0.5+0.5*erf(pra[ns+nr3*slength]));}
      else if (((ns/(int)pow(2,nr))%2==1) && ((ns2/(int)pow(2,nr2))%2==1) && nr2!=nr3 && ((ns2/(int)pow(2,nr3))%2!=1)){ sraprimeprime_input_synapse[ns2+ns*slength+(nr2+nr*length + nr3*length*length)*slength*slength]=sra[ns2+ns*slength]*(+1.0)*exp(-pow(pra[ns+nr2*slength],2))/sqrt(PI)/(0.5+0.5*erf(pra[ns+nr2*slength]))*(-1.0)*exp(-pow(pra[ns+nr3*slength],2))/sqrt(PI)/(0.5-0.5*erf(pra[ns+nr3*slength]));}
      else if (((ns/(int)pow(2,nr))%2==1) && ((ns2/(int)pow(2,nr2))%2==1) && nr2==nr3){ sraprimeprime_input_synapse[ns2+ns*slength+(nr2+nr*length + nr3*length*length)*slength*slength]=sra[ns2+ns*slength]*(+1.0)*exp(-pow(pra[ns+nr2*slength],2))/sqrt(PI)/(0.5+0.5*erf(pra[ns+nr2*slength]))*(-2.0)*pra[ns+nr2*slength];}
      else if (((ns/(int)pow(2,nr))%2==1) && ((ns2/(int)pow(2,nr2))%2!=1) && nr2!=nr3 && ((ns2/(int)pow(2,nr3))%2==1)){ sraprimeprime_input_synapse[ns2+ns*slength+(nr2+nr*length + nr3*length*length)*slength*slength]=sra[ns2+ns*slength]*(-1.0)*exp(-pow(pra[ns+nr2*slength],2))/sqrt(PI)/(0.5-0.5*erf(pra[ns+nr2*slength]))*(+1.0)*exp(-pow(pra[ns+nr3*slength],2))/sqrt(PI)/(0.5+0.5*erf(pra[ns+nr3*slength]));}
      else if (((ns/(int)pow(2,nr))%2==1) && ((ns2/(int)pow(2,nr2))%2!=1) && nr2!=nr3 && ((ns2/(int)pow(2,nr3))%2!=1)){ sraprimeprime_input_synapse[ns2+ns*slength+(nr2+nr*length + nr3*length*length)*slength*slength]=sra[ns2+ns*slength]*(-1.0)*exp(-pow(pra[ns+nr2*slength],2))/sqrt(PI)/(0.5-0.5*erf(pra[ns+nr2*slength]))*(-1.0)*exp(-pow(pra[ns+nr3*slength],2))/sqrt(PI)/(0.5-0.5*erf(pra[ns+nr3*slength]));}
      else if (((ns/(int)pow(2,nr))%2==1) && ((ns2/(int)pow(2,nr2))%2!=1) && nr2==nr3){ sraprimeprime_input_synapse[ns2+ns*slength+(nr2+nr*length + nr3*length*length)*slength*slength]=sra[ns2+ns*slength]*(-1.0)*exp(-pow(pra[ns+nr2*slength],2))/sqrt(PI)/(0.5-0.5*erf(pra[ns+nr2*slength]))*(-2.0)*pra[ns+nr2*slength];}
      else if (((ns/(int)pow(2,nr))%2!=1)){ sraprimeprime_input_synapse[ns2+ns*slength+(nr2+nr*length + nr3*length*length)*slength*slength]=0;}
      if (verbose>1){ printf(" %% computed sraprimeprime_input_synapse[%d+%d*slength+(%d+%d*length + %d*length*length)*slength*slength]=%f\n",ns2,ns,nr2,nr,nr3,sraprimeprime_input_synapse[ns2+ns*slength+(nr2+nr*length + nr3*length*length)*slength*slength]);}}}}
    if (verbose>1){ printf(" %% finished linking state %d to state %d\n",ns,ns2);}}}
  if (verbose){ 
    raprintf(pra,"double",slength,length,"pra[s->n]: ");
    raprintf(sra,"double",slength,slength,"sra[s->s]: ");
    for (nr=0;nr<length;nr++){
      sprintf(text,"sraprime_input_n%d",nr);
      raprintf(&(sraprime_input[0+nr*slength*slength]),"double",slength,slength,text);}
    for (nr=0;nr<length;nr++){ for (nr2=0;nr2<length;nr2++){
      sprintf(text,"sraprime_synapse_n%d->n%d",nr,nr2);
      raprintf(&(sraprime_synapse[0+(nr2+nr*length)*slength*slength]),"double",slength,slength,text);}}}
  eiwr = (double *) tcalloc(slength,sizeof(double));
  eiwi = (double *) tcalloc(slength,sizeof(double));
  eig = (double *) tcalloc(slength*slength,sizeof(double));
  ra2eig(sra,slength,eiwr,eiwi,eig);
  maxmindex("double",eiwr,slength,&Imaxindex,&Inmax,&Iminindex,&Inmin,tolerance);
  if (verbose){ printf(" %% %d largest eigenvalues, first of which is %f\n",Inmax,eiwr[Imaxindex[0]]);}
  stats("double",&(eig[0+Imaxindex[0]*slength]),slength,NULL,NULL,&sum,NULL); sum*=slength;
  for (ns=0;ns<slength;ns++){ statetree_1[ns] = eig[ns+Imaxindex[0]*slength]/sum;}
  if (verbose){ raprintf(statetree_1,"double",1,slength,"statetree_1: ");}
  tfree(Imaxindex);Imaxindex=NULL; tfree(Iminindex);Iminindex=NULL;
  if (statetree_2!=NULL){
    for (ns=0;ns<slength;ns++){ for (ns2=0;ns2<slength;ns2++){
      statetree_2[ns2+ns*slength] = sra[ns2+ns*slength]*statetree_1[ns];}}
    if (verbose){ raprintf(statetree_2,"double",slength,slength," %% statetree_2: ");}}
  minit0(ms1,slength,1);memcpy(ms1->mtrx,statetree_1,slength*sizeof(double));
  mtfree(ms2p2);ms2p2->rows=length*length;ms2p2->cols=slength*slength;ms2p2->mtrx=binary_projection_s2p(0,length,2,2);
  mtfree(ms1p1);ms1p1->rows=length;ms1p1->cols=slength;ms1p1->mtrx=binary_projection_s2p(0,length,1,1);
  minit0(mL,slength,slength);memcpy(mL->mtrx,sra,slength*slength*sizeof(double));
  minit0(mLh,slength*slength,slength);
  for (ns=0;ns<slength;ns++){ for (ns2=0;ns2<slength;ns2++){ ns3=ns; sentry(mLh,ns2+ns*slength,ns3,gentry(mL,ns2,ns));}}
  mtimesm(mLh,ms1,ms2);
  if (verbose){ 
    printf("ms1:\n"); mprintf(ms1);
    raprintf(ms2p2->mtrx,"double",length*length,slength*slength,"s2p2: ");
    raprintf(ms1p1->mtrx,"double",length,slength,"s1p1: ");
    printf("mL:\n"); mprintf(mL); printf("mLh:\n"); mprintf(mLh); printf("ms2:\n"); mprintf(ms2);}
  mtfree(mQ);mQ->rows=slength;mQ->cols=mQ->rows;mQ->mtrx=ds_projection(mQ->rows); mplugout(mQ,0,slength-1,1,slength-1,mQ); mtrans(mQ,mQt);
  minit1(mM,slength,slength); msubtm(mM,mL,mM); mtimesm(mQt,mM,mM); mtimesm(mM,mQ,mM); minv(mM,mM); mtimesm(mQ,mM,mM); mtimesm(mM,mQt,mM);
  if (verbose){ printf("mQ: \n"); mprintf(mQ); printf("mQt: \n"); mprintf(mQt); printf("mM: \n"); mprintf(mM);}
  minit0(mdiLs1,slength,length);
  for (nr=0;nr<length;nr++){
    minit0(mtmp1,slength,slength);memcpy(mtmp1->mtrx,&(sraprime_input[0+nr*slength*slength]),slength*slength*sizeof(double));
    mtimesm(mtmp1,ms1,mtmp2); mplugin(mtmp2,0,slength-1,0,0,mdiLs1,0,nr);}
  mtimesm(mM,mdiLs1,mdis1);
  minit0(mdiLhs1,slength*slength,length);
  for (nr=0;nr<length;nr++){
    minit0(mtmp1,slength,slength);memcpy(mtmp1->mtrx,&(sraprime_input[0+nr*slength*slength]),slength*slength*sizeof(double));
    minit0(mtmp2,slength*slength,slength);
    for (ns=0;ns<slength;ns++){ mplugin(mtmp1,0,slength-1,ns,ns,mtmp2,ns*slength,ns);}
    mtimesm(mtmp2,ms1,mtmp3); mplugin(mtmp3,0,slength*slength-1,0,0,mdiLhs1,0,nr);}
  mtimesm(mLh,mdis1,mdis2); mplusm(mdis2,mdiLhs1,mdis2);
  if (verbose){
    printf("mdiLs1: \n"); mprintf(mdiLs1); printf("mdiLhs1: \n"); mprintf(mdiLhs1);
    printf("mdis1: \n"); mprintf(mdis1); printf("mdis2: \n"); mprintf(mdis2);}
  mtfree(mbp);mbp->rows=(int)pow(length,2);mbp->cols=(int)pow(length,2);mbp->mtrx=bp_projection(length,2);mtimesm(mbp,ms2p2,mbp);
  mtimesm(mbp,ms2,mbps2); 
  if (verbose){     
    raprintf(mbp->mtrx,"double",length*length,slength*slength,"bp: ");
    printf("mbps2: \n"); mprintf(mbps2);}
  minit0(mdsLs1,slength,length*length);
  for (nr=0;nr<length;nr++){ for (nr2=0;nr2<length;nr2++){
    minit0(mtmp1,slength,slength);memcpy(mtmp1->mtrx,&(sraprime_synapse[0+(nr2+nr*length)*slength*slength]),slength*slength*sizeof(double));
    mtimesm(mtmp1,ms1,mtmp2); mplugin(mtmp2,0,slength-1,0,0,mdsLs1,0,nr2+nr*length);}}
  mtimesm(mdsLs1,mbps2,mdss1);mtimesm(mM,mdss1,mdss1);
  minit0(mdsLhs1,slength*slength,length*length);
  for (nr=0;nr<length;nr++){ for (nr2=0;nr2<length;nr2++){
    minit0(mtmp1,slength,slength);memcpy(mtmp1->mtrx,&(sraprime_synapse[0+(nr2+nr*length)*slength*slength]),slength*slength*sizeof(double));
    minit0(mtmp2,slength*slength,slength);
    for (ns=0;ns<slength;ns++){ mplugin(mtmp1,0,slength-1,ns,ns,mtmp2,ns*slength,ns);}
    mtimesm(mtmp2,ms1,mtmp3); mplugin(mtmp3,0,slength*slength-1,0,0,mdsLhs1,0,nr2+nr*length);}}
  mtimesm(mdsLhs1,mbps2,mtmp1); mtimesm(mLh,mdss1,mtmp2); mplusm(mtmp1,mtmp2,mdss2);
  if (verbose){ 
    printf("mdsLs1: \n"); mprintf(mdsLs1); printf("mdss1: \n"); mprintf(mdss1);
    printf("mdsLhs1: \n"); mprintf(mdsLhs1); printf("mdss2: \n"); mprintf(mdss2);}
  minit0(mdiMdsLs1bps2,slength,length);
  for (nr=0;nr<length;nr++){
    minit0(mtmp1,slength,slength);memcpy(mtmp1->mtrx,&(sraprime_input[0+nr*slength*slength]),slength*slength*sizeof(double));
    mtimesm(mQt,mtmp1,mtmp1); mtimesm(mtmp1,mQ,mtmp1); mtimesm(mQ,mtmp1,mtmp1); mtimesm(mtmp1,mQt,mtmp1); 
    mtimesm(mM,mtmp1,mtmp1); mtimesm(mtmp1,mM,mtmp1);
    mtimesm(mtmp1,mdsLs1,mtmp1); mtimesm(mtmp1,mbps2,mtmp1);
    mplugin(mtmp1,0,slength-1,0,0,mdiMdsLs1bps2,0,nr);}
  minit0(mMdidsLs1bps2,slength,length);
  for (nr=0;nr<length;nr++){
    minit0(mtmp2,slength,1);
    for (nr2=0;nr2<length;nr2++){ for (nr3=0;nr3<length;nr3++){
      minit0(mtmp1,slength,slength);memcpy(mtmp1->mtrx,&(sraprimeprime_input_synapse[0+(nr3+nr2*length + nr*length*length)*slength*slength]),slength*slength*sizeof(double));
      mtimesm(mtmp1,ms1,mtmp1); mtimesd(mtmp1,gentry(mbps2,nr3+nr2*length,0),mtmp1); mplusm(mtmp2,mtmp1,mtmp2);}}
    mplugin(mtmp2,0,slength-1,0,0,mMdidsLs1bps2,0,nr);}
  mtimesm(mM,mMdidsLs1bps2,mMdidsLs1bps2);
  minit0(mMdsLdis1bps2,slength,length);
  for (nr=0;nr<length;nr++){ for (nr2=0;nr2<length;nr2++){
    minit0(mtmp1,slength,slength);memcpy(mtmp1->mtrx,&(sraprime_synapse[0+(nr2+nr*length)*slength*slength]),slength*slength*sizeof(double));
    mtimesm(mtmp1,mM,mtmp2); mtimesm(mtmp2,mdiLs1,mtmp1); mtimesd(mtmp1,gentry(mbps2,nr2+nr*length,0),mtmp1); 
    mplusm(mMdsLdis1bps2,mtmp1,mMdsLdis1bps2);}}
  mtimesm(mM,mMdsLdis1bps2,mMdsLdis1bps2);
  mtimesm(mdsLs1,mbp,mMdsLs1bpdis2); mtimesm(mMdsLs1bpdis2,mdis2,mMdsLs1bpdis2); mtimesm(mM,mMdsLs1bpdis2,mMdsLs1bpdis2);
  mplusm(mdiMdsLs1bps2,mMdidsLs1bps2,mtmp1); mplusm(mMdsLdis1bps2,mMdsLs1bpdis2,mtmp2); mplusm(mtmp1,mtmp2,mC);
  minit0(mdiLhdss1,slength*slength,length);
  for (nr=0;nr<length;nr++){
    minit0(mtmp1,slength,slength);memcpy(mtmp1->mtrx,&(sraprime_input[0+nr*slength*slength]),slength*slength*sizeof(double));
    minit0(mtmp2,slength*slength,slength);
    for (ns=0;ns<slength;ns++){ mplugin(mtmp1,0,slength-1,ns,ns,mtmp2,ns*slength,ns);}
    mtimesm(mtmp2,mdss1,mtmp3); mplugin(mtmp3,0,slength*slength-1,0,0,mdiLhdss1,0,nr);}
  mtimesm(mLh,mC,mC2); mplusm(mC2,mdiLhdss1,mC2);
  if (verbose){ 
    printf("mdiMdsLs1bps2: \n");mprintf(mdiMdsLs1bps2); printf("mMdidsLs1bps2: \n");mprintf(mMdidsLs1bps2);
    printf("mMdsLdis1bps2: \n");mprintf(mMdsLdis1bps2); printf("mMdsLs1bpdis2: \n");mprintf(mMdsLs1bpdis2);
    printf("mC: \n");mprintf(mC);
    printf("mC2: \n");mprintf(mC2);}
  mcopy(mdis1,mtmp1);
  minit0(msv,minimum(mtmp1->rows,mtmp1->cols),1); minit0(msvl,mtmp1->rows,mtmp1->rows); minit0(msvr,mtmp1->cols,mtmp1->cols);
  ra2svd(mtmp1->mtrx,mtmp1->rows,mtmp1->cols,msv->mtrx,msvl->mtrx,msvr->mtrx);
  for (nr=0;nr<length;nr++){ Gs1[nr] = gentry(msvr,0,nr);} Gs1[length]=gentry(msv,0,0);
  if (verbose){ raprintf(Gs1,"double",1,length+1,"Gs1: ");}
  mtimesm(ms1p1,mdis1,mtmp1);
  minit0(msv,minimum(mtmp1->rows,mtmp1->cols),1); minit0(msvl,mtmp1->rows,mtmp1->rows); minit0(msvr,mtmp1->cols,mtmp1->cols);
  ra2svd(mtmp1->mtrx,mtmp1->rows,mtmp1->cols,msv->mtrx,msvl->mtrx,msvr->mtrx);
  for (nr=0;nr<length;nr++){ Gp1[nr] = gentry(msvr,0,nr);} Gp1[length]=gentry(msv,0,0);
  if (verbose){ raprintf(Gp1,"double",1,length+1,"Gp1: ");}
  mcopy(mdis2,mtmp1);
  minit0(msv,minimum(mtmp1->rows,mtmp1->cols),1); minit0(msvl,mtmp1->rows,mtmp1->rows); minit0(msvr,mtmp1->cols,mtmp1->cols);
  ra2svd(mtmp1->mtrx,mtmp1->rows,mtmp1->cols,msv->mtrx,msvl->mtrx,msvr->mtrx);
  for (nr=0;nr<length;nr++){ Gs2[nr] = gentry(msvr,0,nr);} Gs2[length]=gentry(msv,0,0);
  if (verbose){ raprintf(Gs2,"double",1,length+1,"Gs2: ");}
  mtimesm(ms2p2,mdis2,mtmp1);
  minit0(msv,minimum(mtmp1->rows,mtmp1->cols),1); minit0(msvl,mtmp1->rows,mtmp1->rows); minit0(msvr,mtmp1->cols,mtmp1->cols);
  ra2svd(mtmp1->mtrx,mtmp1->rows,mtmp1->cols,msv->mtrx,msvl->mtrx,msvr->mtrx);
  for (nr=0;nr<length;nr++){ Gp2[nr] = gentry(msvr,0,nr);} Gp2[length]=gentry(msv,0,0);
  if (verbose){ raprintf(Gp2,"double",1,length+1,"Gp2: ");}
  mtrans(mdis1,mtmp1); mtimesm(mtmp1,mC,mtmp1); 
  minit0(meiwr,mtmp1->rows,1); minit0(meiwi,mtmp1->rows,1); minit0(meig,mtmp1->rows,mtmp1->rows);
  ra2eig(mtmp1->mtrx,mtmp1->rows,meiwr->mtrx,meiwi->mtrx,meig->mtrx);
  maxmindex("double",meiwr->mtrx,mtmp1->rows,&Imaxindex,&Inmax,&Iminindex,&Inmin,tolerance);
  if (verbose){ printf(" %% %d largest eigenvalues, first of which is %f\n",Inmax,meiwr->mtrx[Imaxindex[0]]);}
  for (nr=0;nr<length;nr++){ Cs1[nr] = gentry(meig,nr,Imaxindex[0]);} Cs1[length] = gentry(meiwr,Imaxindex[0],0);
  if (verbose){ raprintf(Cs1,"double",1,length+1,"Cs1: ");}
  tfree(Imaxindex);Imaxindex=NULL; tfree(Iminindex);Iminindex=NULL;
  mtimesm(ms1p1,mdis1,mtmp1); mtrans(mtmp1,mtmp1); mtimesm(mtmp1,ms1p1,mtmp1); mtimesm(mtmp1,mC,mtmp1); 
  minit0(meiwr,mtmp1->rows,1); minit0(meiwi,mtmp1->rows,1); minit0(meig,mtmp1->rows,mtmp1->rows);
  ra2eig(mtmp1->mtrx,mtmp1->rows,meiwr->mtrx,meiwi->mtrx,meig->mtrx);
  maxmindex("double",meiwr->mtrx,mtmp1->rows,&Imaxindex,&Inmax,&Iminindex,&Inmin,tolerance);
  if (verbose){ printf(" %% %d largest eigenvalues, first of which is %f\n",Inmax,meiwr->mtrx[Imaxindex[0]]);}
  for (nr=0;nr<length;nr++){ Cp1[nr] = gentry(meig,nr,Imaxindex[0]);} Cp1[length] = gentry(meiwr,Imaxindex[0],0);
  if (verbose){ raprintf(Cp1,"double",1,length+1,"Cp1: ");}
  tfree(Imaxindex);Imaxindex=NULL; tfree(Iminindex);Iminindex=NULL;
    mtrans(mdis2,mtmp1); mtimesm(mtmp1,mC2,mtmp1); 
  minit0(meiwr,mtmp1->rows,1); minit0(meiwi,mtmp1->rows,1); minit0(meig,mtmp1->rows,mtmp1->rows);
  ra2eig(mtmp1->mtrx,mtmp1->rows,meiwr->mtrx,meiwi->mtrx,meig->mtrx);
  maxmindex("double",meiwr->mtrx,mtmp1->rows,&Imaxindex,&Inmax,&Iminindex,&Inmin,tolerance);
  if (verbose){ printf(" %% %d largest eigenvalues, first of which is %f\n",Inmax,meiwr->mtrx[Imaxindex[0]]);}
  for (nr=0;nr<length;nr++){ Cs2[nr] = gentry(meig,nr,Imaxindex[0]);} Cs2[length] = gentry(meiwr,Imaxindex[0],0);
  if (verbose){ raprintf(Cs2,"double",1,length+1,"Cs2: ");}
  tfree(Imaxindex);Imaxindex=NULL; tfree(Iminindex);Iminindex=NULL;
  mtimesm(ms2p2,mdis2,mtmp1); mtrans(mtmp1,mtmp1); mtimesm(mtmp1,ms2p2,mtmp1); mtimesm(mtmp1,mC2,mtmp1); 
  minit0(meiwr,mtmp1->rows,1); minit0(meiwi,mtmp1->rows,1); minit0(meig,mtmp1->rows,mtmp1->rows);
  ra2eig(mtmp1->mtrx,mtmp1->rows,meiwr->mtrx,meiwi->mtrx,meig->mtrx);
  maxmindex("double",meiwr->mtrx,mtmp1->rows,&Imaxindex,&Inmax,&Iminindex,&Inmin,tolerance);
  if (verbose){ printf(" %% %d largest eigenvalues, first of which is %f\n",Inmax,meiwr->mtrx[Imaxindex[0]]);}
  for (nr=0;nr<length;nr++){ Cp2[nr] = gentry(meig,nr,Imaxindex[0]);} Cp2[length] = gentry(meiwr,Imaxindex[0],0);
  if (verbose){ raprintf(Cp2,"double",1,length+1,"Cp2: ");}
  tfree(Imaxindex);Imaxindex=NULL; tfree(Iminindex);Iminindex=NULL;

/*   mtrans(mdss1,mtmp1);mtimesm(mtmp1,mC,mtmp1);  */
/*   Cs1[length]=mfrobnorm(mtmp1);  */
/*   if (Cs1[length]>threshold){ for (nr=0;nr<length;nr++){ Cs1[nr] = gentry(mtmp1,0,nr)/Cs1[length];}} */
/*   else if (Cs1[length]<=threshold){ for (nr=0;nr<length;nr++){ Cs1[nr] = 0;}} */
/*   if (verbose){ raprintf(Cs1,"double",1,length+1,"Cs1: ");} */
/*   mtimesm(ms1p1,mdss1,mtmp2);mtrans(mtmp2,mtmp1);mtimesm(ms1p1,mC,mtmp2);mtimesm(mtmp1,mtmp2,mtmp1); */
/*   Cp1[length]=mfrobnorm(mtmp1);  */
/*   if (Cp1[length]>threshold){ for (nr=0;nr<length;nr++){ Cp1[nr] = gentry(mtmp1,0,nr)/Cp1[length];}} */
/*   else if (Cp1[length]<=threshold){ for (nr=0;nr<length;nr++){ Cp1[nr] = 0;}} */
/*   if (verbose){ raprintf(Cp1,"double",1,length+1,"Cp1: ");} */
/*   mtrans(mdss2,mtmp1);mtimesm(mtmp1,mC2,mtmp1); */
/*   Cs2[length]=mfrobnorm(mtmp1);  */
/*   if (Cs2[length]>threshold){ for (nr=0;nr<length;nr++){ Cs2[nr] = gentry(mtmp1,0,nr)/Cs2[length];}} */
/*   else if (Cs2[length]<=threshold){ for (nr=0;nr<length;nr++){ Cs2[nr] = 0;}} */
/*   if (verbose){ raprintf(Cs2,"double",1,length+1,"Cs2: ");} */
/*   mtimesm(ms2p2,mdss2,mtmp2);mtrans(mtmp2,mtmp1);mtimesm(ms2p2,mC2,mtmp2);mtimesm(mtmp1,mtmp2,mtmp1); */
/*   Cp2[length]=mfrobnorm(mtmp1);  */
/*   if (Cp2[length]>threshold){ for (nr=0;nr<length;nr++){ Cp2[nr] = gentry(mtmp1,0,nr)/Cp2[length];}} */
/*   else if (Cp2[length]<=threshold){ for (nr=0;nr<length;nr++){ Cp2[nr] = 0;}} */
/*   if (verbose){ raprintf(Cp2,"double",1,length+1,"Cp2: ");} */

  mtfree(mtmp1);mtfree(mtmp2);mtfree(mtmp3);
  mtfree(ms1);mtfree(ms2);mtfree(ms2p2);mtfree(ms1p1);
  mtfree(mL);mtfree(mLh);
  mtfree(mQ);mtfree(mQt);
  mtfree(mM);
  mtfree(mdiLs1);mtfree(mdis1);mtfree(mdiLhs1);mtfree(mdis2);mtfree(mdiLhdss1);
  mtfree(mbp);mtfree(mbps2);
  mtfree(mdsLs1);mtfree(mdss1);mtfree(mdsLhs1);mtfree(mdss2);
  mtfree(mdiMdsLs1bps2);mtfree(mMdidsLs1bps2);mtfree(mMdsLdis1bps2);mtfree(mMdsLs1bpdis2);
  mtfree(mC);mtfree(mC2);
  mtfree(msv);mtfree(msvl);mtfree(msvr);
  mtfree(meiwr);mtfree(meiwi);mtfree(meig);
  tfree(mtmp1);tfree(mtmp2);tfree(mtmp3);
  tfree(ms1);tfree(ms2);tfree(ms2p2);tfree(ms1p1);
  tfree(mL);tfree(mLh);
  tfree(mQ);tfree(mQt);
  tfree(mM);
  tfree(mdiLs1);tfree(mdis1);tfree(mdiLhs1);tfree(mdis2);tfree(mdiLhdss1);
  tfree(mbp);tfree(mbps2);
  tfree(mdsLs1);tfree(mdss1);tfree(mdsLhs1);tfree(mdss2);
  tfree(mdiMdsLs1bps2);tfree(mMdidsLs1bps2);tfree(mMdsLdis1bps2);tfree(mMdsLs1bpdis2);
  tfree(mC);tfree(mC2);
  tfree(msv);tfree(msvl);tfree(msvr);
  tfree(meiwr);tfree(meiwi);tfree(meig);
  tfree(pra);pra=NULL; tfree(sra);sra=NULL; tfree(sraprime_input);sraprime_input=NULL; tfree(sraprime_synapse);sraprime_synapse=NULL;
  tfree(sraprimeprime_input_synapse);sraprimeprime_input_synapse=NULL;
  tfree(eiwr);eiwr=NULL; tfree(eiwi);eiwi=NULL; tfree(eig);eig=NULL;
}

void ra2lss(double *ra,int rows,int cols,double *rhs,int rhs_cols,double *lss)
{
  /* assumes ra,rhs are row dominant, and that rows>cols. Need to fix later for the case rows<cols */
  int verbose=0;
  integer m=rows,n=cols,nrhs=rhs_cols,lda=m,ldb=maximum(rows,cols),rank=0,lwork=-1,*iwork=NULL,info=0;
  doublereal rcond=0.0000001, *A=NULL,*B=NULL,*S=NULL,*work=NULL;
  int nr=0,nc=0;
  A = (doublereal *) tcalloc(rows*cols,sizeof(doublereal));
  B = (doublereal *) tcalloc((int)ldb*rhs_cols,sizeof(doublereal));
  S = (doublereal *) tcalloc(minimum(rows,cols),sizeof(doublereal));
  work = (doublereal *) tcalloc(1,sizeof(doublereal));
  iwork = (integer *) tcalloc(1,sizeof(integer));
  for (nr=0;nr<rows;nr++){ for (nc=0;nc<cols;nc++){ A[nr+nc*rows] = ra[nr+nc*rows];}}
  for (nr=0;nr<rows;nr++){ for (nc=0;nc<rhs_cols;nc++){ B[nr+nc*rows] = rhs[nr+nc*rows];}}
  if (verbose){
    raprintf(ra,"double",rows,cols,"ra ");
    printf(" A:\n"); for (nr=0;nr<rows;nr++){ for (nc=0;nc<cols;nc++){ printf("%f ",A[nr+nc*rows]);} printf("\n");}
    printf(" B:\n"); for (nr=0;nr<rows;nr++){ for (nc=0;nc<rhs_cols;nc++){ printf("%f ",B[nr+nc*rows]);} printf("\n");}}
  dgelsd_(&m,&n,&nrhs,A,&lda,B,&ldb,S,&rcond,&rank,work,&lwork,iwork,&info); lwork=(integer)work[0];
  tfree(work);work=NULL; tfree(iwork);iwork=NULL; 
  work = (doublereal *) tcalloc(lwork,sizeof(doublereal));
  iwork = (integer *) tcalloc(lwork,sizeof(integer));
  dgelsd_(&m,&n,&nrhs,A,&lda,B,&ldb,S,&rcond,&rank,work,&lwork,iwork,&info);
  if (verbose){
    printf(" info %d\n",(int)info);
    printf(" S:\n"); for (nr=0;nr<minimum(rows,cols);nr++){ printf("%f ",S[nr]);} printf("\n");
    printf(" B:\n"); for (nr=0;nr<(int)ldb;nr++){ for (nc=0;nc<nrhs;nc++){ printf("%f ",B[nr+nc*(int)ldb]);} printf("\n");}}
  if (lss!=NULL){ for (nr=0;nr<(int)ldb;nr++){ for (nc=0;nc<nrhs;nc++){ lss[nr+nc*(int)ldb]=B[nr+nc*(int)ldb];}}}
  tfree(A);A=NULL;
  tfree(B);B=NULL;
  tfree(S);S=NULL;
  tfree(work);work=NULL;
  tfree(iwork);iwork=NULL;
}

void ra2svd(double *ra,int rows,int cols,double *sv,double *svl,double *svr)
{
  /* assumes ra is row dominant. For some reason A must be row dominant.
     ra = svl*diag(sv)*svr */
  int verbose=0;
  char jobz='A';
  integer m=rows,n=cols,lda=m,ldu=m,ldv=n,lwork=2*12*(int)pow(maximum(rows,cols),2),*iwork=NULL,info=0;
  doublereal *A=NULL,*S=NULL,*U=NULL,*VT=NULL,*work=NULL;
  int nr=0,nc=0;
  double *ra2=NULL,*ra3=NULL,*ra4=NULL;
  A = (doublereal *) tcalloc(rows*cols,sizeof(doublereal));
  S = (doublereal *) tcalloc(minimum(rows,cols),sizeof(doublereal));
  U = (doublereal *) tcalloc(rows*rows,sizeof(doublereal));
  VT = (doublereal *) tcalloc(cols*cols,sizeof(doublereal));
  work = (doublereal *) tcalloc(lwork,sizeof(doublereal));
  iwork = (integer *) tcalloc(8*maximum(rows,cols),sizeof(integer));
  for (nr=0;nr<rows;nr++){ for (nc=0;nc<cols;nc++){ A[nr+nc*rows] = ra[nr+nc*rows];}}
  if (verbose){
    raprintf(ra,"double",rows,cols,"ra ");
    printf(" A:\n"); for (nr=0;nr<rows;nr++){ for (nc=0;nc<cols;nc++){ printf("%f ",A[nr+nc*rows]);} printf("\n");}}
  dgesdd_(&jobz,&m,&n,A,&lda,S,U,&ldu,VT,&ldv,work,&lwork,iwork,&info);
  if (verbose){
    printf(" info %d\n",(int)info);
    printf(" S:\n"); for (nr=0;nr<minimum(rows,cols);nr++){ printf("%f ",S[nr]);} printf("\n");
    printf(" U:\n"); for (nr=0;nr<rows;nr++){ for (nc=0;nc<rows;nc++){ printf("%f ",U[nc+nr*rows]);} printf("\n");}
    printf(" VT:\n"); for (nr=0;nr<cols;nr++){ for (nc=0;nc<cols;nc++){ printf("%f ",VT[nc+nr*cols]);} printf("\n");}}
  if (sv!=NULL){ for (nr=0;nr<minimum(rows,cols);nr++){ sv[nr]=S[nr];}}
  if (svl!=NULL){ for (nr=0;nr<rows;nr++){ for (nc=0;nc<rows;nc++){ svl[nr+nc*rows]=U[nr+nc*rows];}}}
  if (svr!=NULL){ for (nr=0;nr<cols;nr++){ for (nc=0;nc<cols;nc++){ svr[nr+nc*cols]=VT[nr+nc*cols];}}}
  if (verbose){
    ra2 = (double *) tcalloc(rows*cols,sizeof(double)); for (nr=0;nr<minimum(rows,cols);nr++){ ra2[nr+nr*rows]=sv[nr];}
    ra3 = ra2ra_matrix_multiply(ra2,rows,cols,0,svr,cols,cols,0);
    ra4 = ra2ra_matrix_multiply(svl,rows,rows,0,ra3,rows,cols,0);
    raprintf(ra4,"double",rows,cols,"A?: ");
    printf("optimal work: %f\n",work[0]);
    tfree(ra2);ra2=NULL; tfree(ra3);ra3=NULL; tfree(ra4);ra4=NULL;}
  tfree(A);A=NULL;
  tfree(S);S=NULL;
  tfree(U);U=NULL;
  tfree(VT);VT=NULL;
  tfree(work);work=NULL;
  tfree(iwork);iwork=NULL;
}

void ra2eig(double *ra,int length,double *eiw_real,double *eiw_imag,double *eig_vect)
{
  /* assumes ra is row dominant 
     for some reason VL ends up storing right eigenvectors
     if eigenvalues come in complex conjugate pairs, the positive imaginary part comes first (i.e., entry j and j+1),
     and the corresponding eigenvectors are eig_vect(:,j)+i*eig_vect(:,j+1) and eig_vect(:,j)-i*eig_vect(:,j+1)
  */
  int verbose=0;
  char jobvl='V',jobvr='N';
  integer N=length,lda=length,ldvl=length,ldvr=length,lwork=maximum(4,length)*length,info=0;
  doublereal *A=NULL,*WR=NULL,*WI=NULL,*VL=NULL,*VR=NULL,*work=NULL;
  int nr=0,nc=0;
  double *ra2=NULL,*ra3=NULL;
  A = (doublereal *) tcalloc(length*length,sizeof(doublereal));
  WR = (doublereal *) tcalloc(length,sizeof(doublereal));
  WI = (doublereal *) tcalloc(length,sizeof(doublereal));
  VL = (doublereal *) tcalloc(length*length,sizeof(doublereal));
  VR = (doublereal *) tcalloc(length*length,sizeof(doublereal));
  work = (doublereal *) tcalloc(maximum(4,length)*length,sizeof(doublereal));
  for (nr=0;nr<length;nr++){ for (nc=0;nc<length;nc++){ A[nc+nr*length] = ra[nr+nc*length];}}
  if (verbose){
    raprintf(ra,"double",length,length,"ra ");
    printf(" A:\n"); for (nr=0;nr<length;nr++){ for (nc=0;nc<length;nc++){ printf("%f ",A[nc+nr*length]);} printf("\n");}}
  dgeev_(&jobvl,&jobvr,&N,A,&lda,WR,WI,VL,&ldvl,VR,&ldvr,work,&lwork,&info);
  if (verbose){
    printf(" WR:\n"); for (nr=0;nr<length;nr++){ printf("%f ",WR[nr]);} printf("\n");
    printf(" WI:\n"); for (nr=0;nr<length;nr++){ printf("%f ",WI[nr]);} printf("\n");
    printf(" VL:\n"); for (nr=0;nr<length;nr++){ for (nc=0;nc<length;nc++){ printf("%f ",VL[nr+nc*length]);} printf("\n");}
    printf(" VR:\n"); for (nr=0;nr<length;nr++){ for (nc=0;nc<length;nc++){ printf("%f ",VR[nr+nc*length]);} printf("\n");}}
  if (verbose){
    ra2 = (double *) tcalloc(length*length,sizeof(double));
    for (nr=0;nr<length;nr++){ for (nc=0;nc<length;nc++){ ra2[nr+nc*length]=VL[nr+nc*length];}}
    ra3 = ra2ra_matrix_multiply(ra,length,length,0,ra2,length,length,0);
    raprintf(ra3,"double",length,length,"A*VR: ");
    tfree(ra2);ra2=NULL; tfree(ra3);ra3=NULL;}
  if (eiw_real!=NULL){ for (nr=0;nr<length;nr++){ eiw_real[nr] = WR[nr];}}
  if (eiw_imag!=NULL){ for (nr=0;nr<length;nr++){ eiw_imag[nr] = WI[nr];}}
  if (eig_vect!=NULL){ for (nr=0;nr<length;nr++){ for (nc=0;nc<length;nc++){ eig_vect[nr+nc*length] = VL[nr+nc*length];}}}
  tfree(A);A=NULL; tfree(WR);WR=NULL; tfree(WI);WI=NULL; tfree(VL);VL=NULL; tfree(VR);VR=NULL; tfree(work);work=NULL;
}

/* Here are the snxdata functions */

void snxdatadump(struct snxdata *s)
{
  /* this function writes a matlab script which carries out the subnetwork expansion to derive the firing rates, 
     assuming GLOBAL_SNX_VERSION==1 */
  int verbose=1;
  int version_flag = GLOBAL_SNX_VERSION;
  int total_Vs_states = 0, total_states = 0;
  char filename[1024];
  FILE *fp=NULL;
  double *L=NULL,*DL=NULL;
  int nr=0,nr1=0,nr2=0,ns=0,ns1=0,ns2=0;
  switch (version_flag){
  case 1: if (verbose){ printf(" %% snxdatadump designed for version %d\n",version_flag);}
    total_Vs_states = GLOBAL_SNX_NSTATES_[VARNAME_REGISTRY_snx_Vs]; 
    total_states = total_Vs_states;
    sprintf(filename,"snxdata_chain_0_%srecord.m",GLOBAL_STRING_2);
    if ((fp=fopen(filename,"w"))==NULL){ printf(" %% Can't create %s in snxdatadump\n",filename); fp=stdout;}
    fprintf(fp,";\n");
    fprintf(fp," %% matlab script for snxdata_chain_0_%srecord.m;\n",GLOBAL_STRING_2);
    fprintf(fp," L = zeros(%d,%d,%d);\n",total_states,total_states,s->length);
    fprintf(fp," DL = zeros(%d,%d,%d);\n",total_states,total_states,s->length);
    L=(double *) tcalloc(total_states*total_states,sizeof(double)); DL=(double *) tcalloc(total_states*total_states,sizeof(double));
    for (nr=0;nr<s->length;nr++){ 
      for (ns=0;ns<total_Vs_states;ns++){
	snx_statematrix_1(ns,total_Vs_states,s->input[nr],&(L[0+ns*total_states]),&(DL[0+ns*total_states]));}
      for (ns1=0;ns1<total_Vs_states;ns1++){ for (ns2=0;ns2<total_Vs_states;ns2++){
	fprintf(fp," L(%d,%d,%d)=%0.16f;\n",ns1+1,ns2+1,nr+1,L[ns1+ns2*total_states]);
	fprintf(fp," DL(%d,%d,%d)=%0.16f;\n",ns1+1,ns2+1,nr+1,DL[ns1+ns2*total_states]);}}}
    tfree(L);L=NULL; tfree(DL);DL=NULL;
    fprintf(fp," number_of_neurons = %d;\n",s->length);
    fprintf(fp," connectivity = zeros(%d,%d); %% connectivity(1,2) is connection from 1 --> 2; \n",s->length,s->length);
    for (nr1=0;nr1<s->length;nr1++){ for (nr2=0;nr2<s->length;nr2++){ 
      fprintf(fp," connectivity(%d,%d)=%0.16f;\n",nr1+1,nr2+1,s->connectivity[nr1 + nr2*s->length]);}}
    fprintf(fp," input = zeros(%d,1);\n",s->length);
    for (nr=0;nr<s->length;nr++){ 
      fprintf(fp," input(%d)=%0.16f;\n",nr+1,s->input[nr]);}
    fprintf(fp," chain_0 = zeros(%d,1);\n",s->length);
    for (nr=0;nr<s->length;nr++){ 
      fprintf(fp," chain_0(%d)=%0.16f;\n",nr+1,s->chain_0[nr]/s->total_time);}
    fprintf(fp,";\n");
    fprintf(fp,";\n");
    fprintf(fp,";\n");
    fprintf(fp,";\n");
    fprintf(fp,";\n");
    fprintf(fp,";\n");
    fprintf(fp,";\n");
    fprintf(fp,";\n");
    fprintf(fp,";\n");
    fprintf(fp,";\n");
    fprintf(fp,";\n");
    if (fp!=stdout){ fclose(fp); fp=NULL;}
    break;
  case 2: if (verbose){ printf(" %% snxdatadump designed for version %d\n",version_flag);}
    sprintf(filename,"snxdata_abc_bca_%srecord.txt",GLOBAL_STRING_2);
    if ((fp=fopen(filename,"w"))==NULL){ printf(" %% Can't create %s in snxdatadump\n",filename); fp=stdout;}
    fprintf(fp,"%0.16f,%0.16f\n",s->abc/s->total_time,s->bca/s->total_time); 
    if (fp!=stdout){ fclose(fp); fp=NULL;}
    break;
  default: printf(" %% warning! snxdatadump not designed for version %d\n",version_flag); break;}
}

struct snxdata *snxdatamake(struct neuronarray *Nra,int lookback)
{
  int verbose=0;
  int connectivity_type=2;
  struct snxdata *s=NULL;
  int nr1=0,nr2=0;
  s = (struct snxdata *) tcalloc(1,sizeof(struct snxdata));
  s->Nra = Nra;
  s->length = s->Nra->lengthra[0];
  s->lookback = maximum(1,lookback);
  s->stra = (struct strobe **) tcalloc(s->length,sizeof(struct strobe *));
  for (nr1=0;nr1<s->length;nr1++){ s->stra[nr1] = strobemake(s->lookback+3,1,0,0,0,0);}
  s->chain_0 = (double *) tcalloc(s->length,sizeof(double));
  s->chain_1_ = (double *) tcalloc(s->length*s->lookback,sizeof(double));
  s->total_time = 0;
  s->input = (double *) tcalloc(s->length,sizeof(double));
  s->connectivity = (double *) tcalloc(s->length*s->length,sizeof(double));
  s->abc = 0; s->bca = 0;
  switch(GLOBAL_SNX_VERSION){
  case 1:
    switch (connectivity_type){
    case 0: /* random */
      for (nr1=0;nr1<s->length;nr1++){ s->input[nr1] = (1 + 0.1*(2*rand01-1))*CS_ORN_[TYPENAME_REGISTRY_snx_PN];}
      for (nr1=0;nr1<s->length;nr1++){ for (nr2=0;nr2<s->length;nr2++){ if (nr1!=nr2){ s->connectivity[nr1+nr2*s->length] = CS__[TYPENAME_REGISTRY_snx_PN+TYPENAME_REGISTRY_snx_PN*GLOBAL_NTYPES+0/* GLOBAL_INDEXING_REFILE_sra[VARNAME_REGISTRY_snx_nAch] */*GLOBAL_NTYPES*GLOBAL_NTYPES]*(2*rand01-1);}}}
      break;
    case 1: /* ring */
      for (nr1=0;nr1<s->length;nr1++){ s->input[nr1] = CS_ORN_[TYPENAME_REGISTRY_snx_PN];}
      for (nr1=0;nr1<s->length;nr1++){ for (nr2=0;nr2<s->length;nr2++){ if (nr1==periodize(nr2-1,0,s->length)){ s->connectivity[nr1+nr2*s->length] = CS__[TYPENAME_REGISTRY_snx_PN+TYPENAME_REGISTRY_snx_PN*GLOBAL_NTYPES+0/* GLOBAL_INDEXING_REFILE_sra[VARNAME_REGISTRY_snx_nAch] */*GLOBAL_NTYPES*GLOBAL_NTYPES];}}}
      break;
    case 2: /* ring with multiple positive weak connections leading to nr==1 */
      for (nr1=0;nr1<s->length;nr1++){ s->input[nr1] = CS_ORN_[TYPENAME_REGISTRY_snx_PN];}
      for (nr1=0;nr1<s->length;nr1++){ for (nr2=0;nr2<s->length;nr2++){ 
	if (nr1==periodize(nr2-1,0,s->length)){ s->connectivity[nr1+nr2*s->length] = CS__[TYPENAME_REGISTRY_snx_PN+TYPENAME_REGISTRY_snx_PN*GLOBAL_NTYPES+0/* GLOBAL_INDEXING_REFILE_sra[VARNAME_REGISTRY_snx_nAch] */*GLOBAL_NTYPES*GLOBAL_NTYPES];}
	else /* if (nr1!=periodize(nr2-1,0,s->length)) */{ if (nr2==1){ s->connectivity[nr1+nr2*s->length] = pow(CS__[TYPENAME_REGISTRY_snx_PN+TYPENAME_REGISTRY_snx_PN*GLOBAL_NTYPES+0/* GLOBAL_INDEXING_REFILE_sra[VARNAME_REGISTRY_snx_nAch] */*GLOBAL_NTYPES*GLOBAL_NTYPES],2);}}}}
      break;
    default: break;}
    break;
  case 2: 
    switch (connectivity_type){
    case 2: /* two way ring */
      for (nr1=0;nr1<s->length;nr1++){ s->input[nr1] = CS_ORN_[TYPENAME_REGISTRY_snx_PN];}
      for (nr1=0;nr1<s->length;nr1++){ for (nr2=0;nr2<s->length;nr2++){ if (nr1==periodize(nr2-1,0,s->length) || nr1==periodize(nr2+1,0,s->length)){ s->connectivity[nr1+nr2*s->length] = CS__[TYPENAME_REGISTRY_snx_PN+TYPENAME_REGISTRY_snx_PN*GLOBAL_NTYPES+0/* GLOBAL_INDEXING_REFILE_sra[VARNAME_REGISTRY_snx_nAch] */*GLOBAL_NTYPES*GLOBAL_NTYPES];}}}
      break;
    default: break;}
  default: break;}
  if (verbose){ raprintf(s->connectivity,"double",s->length,s->length,"s->connectivity: ");}
  return s;
}

void snxdatatfree(struct snxdata *s)
{
  int nr1=0;
  for (nr1=0;nr1<s->length;nr1++){ strobetfree(s->stra[nr1]); s->stra[nr1]=NULL;}
  tfree(s->stra); s->stra=NULL;
  tfree(s->chain_0); s->chain_0=NULL;
  tfree(s->chain_1_); s->chain_1_=NULL;
  tfree(s->input); s->input=NULL;
  tfree(s->connectivity); s->connectivity=NULL;
  tfree(s);s=NULL;
}

void snxdataupdate(struct snxdata *s,double t,double DT)
{
  int verbose=0;
  int nr=0,nl=0,nr1=0,nr2=0,nr3=0,tab=0,tab1=0,tab2=0,tab3=0;
  double d1=0,d2=0,d3=0;
  struct neuron *n=NULL;
  if (verbose){ printf(" %% [entering snxdataupdate] t %f dt %f\n",t,DT);}
  if (s->stra!=NULL){
  for (nr=0;nr<s->length;nr++){
    n=nget(s->Nra,0,nr);
    strobeupdate(s->stra[nr],t,DT,*(n->vpra[VARNAME_REGISTRY_snx_spike_flag]));
    if (verbose>2){ printf("nr %d tab %d",nr,s->stra[nr]->tab); raprintf(s->stra[nr]->data,"double",1,s->stra[nr]->length,"data: ");}}}
  for (nr=0;nr<s->length;nr++){
    tab = periodize(s->stra[nr]->tab-1,0,s->stra[nr]->length);
    if (s->stra[nr]->data[tab] >= 0.5){ 
      if (verbose){ printf(" %% n(%d) fired",nr);}
      s->chain_0[nr] += 1;
      nr2 = periodize(nr-1,0,s->length);
      for (nl=0;nl<s->lookback;nl++){
	tab2 = periodize(tab-1-nl,0,s->stra[nr]->length);
	if (s->stra[nr2]->data[tab2]>=0.5){
	  if (verbose){ printf(" n(%d) fired %d ago,",nr2,nl+1);}
	  s->chain_1_[nr + nl*s->length] += 1;}}
      if (verbose){ printf("\n");}}}
  for (nr=0;nr<s->length;nr++){
    nr1 = periodize(nr-0,0,s->length); tab1 = periodize(s->stra[nr1]->tab-1,0,s->stra[nr1]->length); d1 = s->stra[nr1]->data[tab1];
    if (d1 >= 0.5){      
      nr2 = periodize(nr-1,0,s->length); tab2 = periodize(s->stra[nr2]->tab-2,0,s->stra[nr2]->length); d2 = s->stra[nr2]->data[tab2];
      if (d2 >= 0.5){
	nr3 = periodize(nr-2,0,s->length); tab3 = periodize(s->stra[nr3]->tab-3,0,s->stra[nr3]->length); d3 = s->stra[nr3]->data[tab3];
	if (d3 >= 0.5){ s->abc += 1; if (verbose){ printf(" %% n(%d),n(%d),n(%d) fired in sequence",nr3,nr2,nr1);}}}
      nr2 = periodize(nr+1,0,s->length); tab2 = periodize(s->stra[nr2]->tab-2,0,s->stra[nr2]->length); d2 = s->stra[nr2]->data[tab2];
      if (d2 >= 0.5){
	nr3 = periodize(nr+2,0,s->length); tab3 = periodize(s->stra[nr3]->tab-3,0,s->stra[nr3]->length); d3 = s->stra[nr3]->data[tab3];
	if (d3 >= 0.5){ s->abc += 1; if (verbose){ printf(" %% n(%d),n(%d),n(%d) fired in sequence",nr3,nr2,nr1);}}}      
      nr2 = periodize(nr-2,0,s->length); tab2 = periodize(s->stra[nr2]->tab-2,0,s->stra[nr2]->length); d2 = s->stra[nr2]->data[tab2];
      if (d2 >= 0.5){
	nr3 = periodize(nr-1,0,s->length); tab3 = periodize(s->stra[nr3]->tab-3,0,s->stra[nr3]->length); d3 = s->stra[nr3]->data[tab3];
	if (d3 >= 0.5){ s->bca += 1; if (verbose){ printf(" %% n(%d),n(%d),n(%d) fired out of sequence",nr3,nr2,nr1);}}}
      nr2 = periodize(nr+2,0,s->length); tab2 = periodize(s->stra[nr2]->tab-2,0,s->stra[nr2]->length); d2 = s->stra[nr2]->data[tab2];
      if (d2 >= 0.5){
	nr3 = periodize(nr+1,0,s->length); tab3 = periodize(s->stra[nr3]->tab-3,0,s->stra[nr3]->length); d3 = s->stra[nr3]->data[tab3];
	if (d3 >= 0.5){ s->bca += 1; if (verbose){ printf(" %% n(%d),n(%d),n(%d) fired out of sequence",nr3,nr2,nr1);}}}}}
  s->total_time += DT;
  if (verbose){ printf(" %% [finished with snxdataupdate]\n");}
}

/* Here are the rho functions */

struct rho * rhomake(struct neuronarray *Nra,int length,int indexing_nvar_length,int *indexing_nvar_checkout,int *indexing_nvar_refile,double *maxra,double *minra,int *nbinra)
{
  int verbose=0;
  int nr=0,nv=0;
  struct rho *r=NULL;
  if (verbose){ 
    printf(" %% [entering rhomake]\n");
    raprintf(indexing_nvar_checkout,"int",1,indexing_nvar_length,"nvar_checkout");
    raprintf(indexing_nvar_refile,"int",1,Nra->nvars,"nvar_refile");
    raprintf(maxra,"double",1,Nra->nvars,"maxra");
    raprintf(minra,"double",1,Nra->nvars,"minra");}
  r = (struct rho *) tmalloc(sizeof(struct rho));
  r->Nra=Nra;
  r->length=length;
  r->neuronllistra = (struct llist **) tcalloc(r->length,sizeof(struct llist *));
  for (nr=0;nr<r->length;nr++){ r->neuronllistra[nr] = llistmake(); rho_neuronllistmake(r,nr);}  
  r->indexing_nvar_length = indexing_nvar_length;
  r->indexing_nvar_checkout = (int *) tcalloc(r->indexing_nvar_length,sizeof(int));
  for (nv=0;nv<r->indexing_nvar_length;nv++){ r->indexing_nvar_checkout[nv] = indexing_nvar_checkout[nv];}
  r->indexing_nvar_refile = (int *) tcalloc(r->Nra->nvars,sizeof(int));
  for (nv=0;nv<r->Nra->nvars;nv++){ r->indexing_nvar_refile[nv] = indexing_nvar_refile[nv];}
  r->maxra= (double *) tcalloc(r->Nra->nvars,sizeof(double));
  r->minra= (double *) tcalloc(r->Nra->nvars,sizeof(double));
  for (nv=0;nv<Nra->nvars;nv++){ r->maxra[nv] = maxra[nv]; r->minra[nv] = minra[nv];}
  r->nbinra = (int *) tcalloc(r->Nra->nvars,sizeof(int));
  for (nv=0;nv<Nra->nvars;nv++){ r->nbinra[nv] = nbinra[nv];}
  r->nbins = 1; for (nv=0;nv<r->indexing_nvar_length;nv++){ r->nbins *= r->nbinra[r->indexing_nvar_checkout[nv]];}
  r->rhora = (double **) tcalloc(r->length,sizeof(double *));
  for (nr=0;nr<r->length;nr++){ r->rhora[nr] = (double *) tcalloc(r->nbins,sizeof(double));}
  r->total_time = 0;
  return r;
}

void rhotfree(struct rho *r)
{
  int nr=0;
  for (nr=0;nr<r->length;nr++){ llisttfree(r->neuronllistra[nr]);r->neuronllistra[nr]=NULL;} tfree(r->neuronllistra);r->neuronllistra=NULL;
  tfree(r->indexing_nvar_checkout);r->indexing_nvar_checkout=NULL;
  tfree(r->indexing_nvar_refile);r->indexing_nvar_refile=NULL;
  tfree(r->maxra);r->maxra=NULL;
  tfree(r->minra);r->minra=NULL;
  tfree(r->nbinra);r->nbinra=NULL;
  for (nr=0;nr<r->length;nr++){ tfree(r->rhora[nr]);r->rhora[nr]=NULL;} tfree(r->rhora);r->rhora=NULL;
  tfree(r);r=NULL;
}

int rhoupdate_getbins(struct rho *r,struct neuron *n,int *ira)
{
  /* we presume that *ira is appropriately initialized with r->indexing_nvar_length elements 
     returns the appropriate nested multi-index */
  int verbose=0;
  int nv=0,nb=0,tab=0;
  if (verbose){ printf(" %% [entering rhoupdate_getbins]\n");}
  for (nv=0;nv<r->indexing_nvar_length;nv++){
    nb = crop((int)floor(r->nbinra[r->indexing_nvar_checkout[nv]]*(*(n->vpra[r->indexing_nvar_checkout[nv]])-r->minra[r->indexing_nvar_checkout[nv]])/(r->maxra[r->indexing_nvar_checkout[nv]]-r->minra[r->indexing_nvar_checkout[nv]])),0,r->nbinra[r->indexing_nvar_checkout[nv]]);
    ira[nv] = nb;}
  tab=0;
  for (nv=r->indexing_nvar_length-1;nv>=0;nv--){
    tab += ira[nv];
    if (nv>0){ tab *= r->nbinra[r->indexing_nvar_checkout[nv-1]];}}
  if (verbose){ raprintf(ira,"int",1,r->indexing_nvar_length," %% ira: "); printf(" %% multi-index %d\n",tab);}
  if (verbose){ printf(" %% [finishing rhoupdate_getbins]\n");}
  return tab;
}

void rhoupdate(struct rho *r,double t,double DT)
{
  int verbose=0;
  int nr=0;
  struct litem *l=NULL;
  struct neuron *n=NULL;
  int *ira=NULL;
  if (verbose){ printf(" %% [entering rhoupdate] t %f DT %f\n",t,DT);}
  ira = (int *) tcalloc(r->indexing_nvar_length,sizeof(int));
  for (nr=0;nr<r->length;nr++){ 
    if (verbose){ printf(" %% updating r->rhora[%d]\n",nr);}
    l = r->neuronllistra[nr]->first;
    while (l!=NULL){
      n = (struct neuron *) l->item;
      if (verbose){ printf(" %% currently considering neuron (%d,%d)\n",n->type,n->index);}
      r->rhora[nr][rhoupdate_getbins(r,n,ira)] += DT;
      l=l->child;}}
  r->total_time += DT;
  tfree(ira);ira=NULL;
  if (verbose){ printf(" %% [finishing rhoupdate]\n");}
}

void rho_neuronllistmake(struct rho *r,int nl)
{
  /* very simple function which takes the first nr/r->length neurons and puts them into r->neuronllistra[nr] */
  int verbose=0;
  int nr=0,nt=0;
  if (verbose){ printf(" %% [entering rho_neuronllistmake] with index %d\n",nl);}
  nt=0;
  for (nr=nl*(r->Nra->lengthra[nt]/r->length);nr<(nl+1)*(r->Nra->lengthra[nt]/r->length);nr++){
    if (nr>=0 && nr<r->Nra->lengthra[nt]){ 
      if (verbose){ printf(" %% adding neuron (%d,%d)\n",nt,nr);}
      litemadd(r->neuronllistra[nl],nget(r->Nra,nt,nr));}}
  if (verbose){ printf(" %% [finishing rho_neuronllistmake]\n");}
}

void rhodump(struct rho *r,char *fgvn,int dump_type)
{
  /* fix later */
}
