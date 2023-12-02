#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define maximum(A,B) ((A) > (B) ? (A) : (B))
#define minimum(A,B) ((A) < (B) ? (A) : (B))
#define periodize(A,B,C) ((A) < (B) ? (A) + (C) - (B) : ((A) >= (C) ? (A) - (C) + (B) : (A)))
#define crop(A,B,C) ((A) < (B) ? (B) : ((A) > (C) ? (C) : (A)))
#define rand01 ((double)rand()/(double)RAND_MAX)

static char PBS_PREAMBLE_1[128]="#!/bin/bash\n#PBS -N ";
static char PBS_PREAMBLE_2[512]="\n#PBS -l nodes=1:ppn=2\n#PBS -l walltime=99:00:00\n#PBS -o /home/sdc4/adi_runs/adi_output/\n#PBS -e /home/sdc4/adi_runs/adi_error/\nARCH=$(uname -m)\n#path to your compiler - may need to add gnu\nIBMCMP=/opt/ibmcmp\nPATH=$MPICH/bin:$PATH\nPATH=$IBMCMP/vacpp/7.0/bin:$PATH\nPATH=$IBMCMP/xlf/9.1/bin:$PATH\nulimit -s unlimited\ncd /home/sdc4/adi_runs\n";

void indextract(int dindex,int nvar,int *dlength,int *indexra)
{
  /* extract individual indices from a bound multi-index */
  /* for example, (127,4,[2 5 3 7]) yields [a + b*2 + c*2*5 + d*2*5*3] yields [1 + 3*2 + 0*2*5 + 4*2*5*3 ] yields [1 3 0 4] */
  /* both int *dlength and int *indexra should be initialized */
  int curvar=0;
  if (nvar > 0){
    curvar = 0;
    do{
      indexra[curvar] = dindex%dlength[curvar];
      dindex = (dindex-indexra[curvar])/dlength[curvar];
      curvar += 1;}
    while (curvar<nvar);}
}

void schain(char *output,int nvar,int *indexra)
{
  /* given int *indexra [1 2 3 4] we create char *output "_1_2_3_4" */
  /* we assume output has the space to handle itself */
  char chartemp[256];
  int vindex=0;
  sprintf(output,"");
  for  (vindex=0;vindex<nvar;vindex++){ sprintf(chartemp,"%s_%d",output,indexra[vindex]); sprintf(output,"%s",chartemp);}
}

void ping(){printf("ping\n");}
void pong(){printf("pong\n");}

int main(int argc, char **argv)
{
  int verbose=2;
  char inputname[128],outputname[256],vname[128],helpfile[1024],programname[64],tuningcurvevsrate[64];
  char equals[1024],comma_vs_semicolon[1024];
  int still_have_options=argc-1;
  FILE *ip=NULL;
  FILE **fp=NULL;
  int charsize=64,nvar=0;
  char *chartemp=NULL,**keyvarra=NULL,**keyvartypera=NULL;
  int *keyvarindexra=NULL;
  int varindex=0;
  int match_index_and_name=0,matching_vindex=0;
  int rtemp=0,*rangera=NULL,*scaletypera=NULL,*multorsumra=NULL;
  double ftemp=0,*factorra=NULL;
  double **rangefactorsra=NULL;
  int dim=0,vindex=0,rindex=0,tindex=0,r2=0,*indexra=NULL,dindex=0,found_flag=0;
  FILE *hp=NULL;
  char gstring[64];
  int grecord=0,maxrecord=1;
  int seqorsim=0,perprocessor=0;
  int runfilenumber=0;
  sprintf(helpfile,"-Finputname [-Pperprocessor] [-Mmaxrecord] [-Nprogramname] {repeat the following} -Z -Vkeyvarname, -Iindex -Tkeyvartype, -Rrange -Xfactor\n R- balanced range, R+ increasing range\n X+ factor X- sum\n");
  sprintf(programname,"**PROGRAM_NAME**");
  sprintf(tuningcurvevsrate,"**TUNINGCURVEVSRATE**");
  if (verbose>1){ printf("starting inputcreate...\n");}
  if (still_have_options==0){ printf(helpfile); exit(EXIT_SUCCESS);}
  while (still_have_options){
    switch(getopt(argc,argv,"F:f:P:p:M:m:N:n:V:v:I:i:T:t:R:r:X:x:Z")){
    case 'F': case 'f': if (verbose>1){ printf("option Ff\n");}
      sprintf(inputname,optarg); if (verbose>1){ printf(" inputname %s\n",inputname);}
      break;
    case 'P': case 'p': if (verbose>1){ printf("option Pp\n");} 
      seqorsim = atoi(optarg); /* ; or & */
      perprocessor = abs(seqorsim);
      break;
    case 'M': case 'm': if (verbose>1){ printf("option Mm\n");} 
      grecord = atoi(optarg); maxrecord = maximum(1,grecord);
      break;
    case 'N': case 'n': if (verbose>1){ printf("option Nn\n");} 
      sprintf(programname,optarg); if (verbose>1){ printf(" programname %s\n",programname);}
      break;
    case 'Z': if (verbose>1){ printf("option Z\n");}
      nvar += 1; 
      keyvarra = (char **) realloc(keyvarra,nvar*sizeof(char *));
      keyvarindexra = (int *) realloc(keyvarindexra,nvar*sizeof(int));
      keyvartypera = (char **) realloc(keyvartypera,nvar*sizeof(char *));
      scaletypera = (int *) realloc(scaletypera,nvar*sizeof(int));
      multorsumra = (int *) realloc(multorsumra,nvar*sizeof(int));
      rangera = (int *) realloc(rangera,nvar*sizeof(int));
      factorra = (double *) realloc(factorra,nvar*sizeof(double));
      break;
    case 'V': case 'v': if (verbose>1){ printf("option Vv\n");}
      keyvarra[nvar-1] = (char *) malloc(charsize*sizeof(char));
      sprintf(keyvarra[nvar-1],optarg); 
      break;
    case 'I': case 'i': if (verbose>1){ printf("option Ii\n");}
      keyvarindexra[nvar-1] = atoi(optarg);
      break;
    case 'T': case 't': if (verbose>1){ printf("option Tt\n");}
      keyvartypera[nvar-1] = (char *) malloc(charsize*sizeof(char));
      sprintf(keyvartypera[nvar-1],optarg);
      break;
    case 'R': case 'r': if (verbose>1){ printf("option Rf\n");} 
      rtemp = atoi(optarg); 
      scaletypera[nvar-1] = (rtemp>0); rangera[nvar-1] = abs(rtemp); 
      break;
    case 'X': case 'x': if (verbose>1){ printf("option Xx\n");} 
      ftemp = atof(optarg);
      multorsumra[nvar-1] = (ftemp>0); factorra[nvar-1] = fabs(ftemp); 
      break;
    default: printf(helpfile); exit(EXIT_SUCCESS); break;}
    still_have_options--;}
  if (verbose>1){ 
    printf("finished reading options:\n");
    for (vindex=0;vindex<nvar;vindex++){
      printf("keyvar %s of type %s index %d read with range %d (%s) and factor %f (%s)\n",keyvarra[vindex],keyvartypera[vindex],keyvarindexra[vindex],rangera[vindex],scaletypera[vindex]?"increasing":"balanced",factorra[vindex],multorsumra[vindex]?"mult":"sum");}}
  dim = 1;
  rangefactorsra = (double **) calloc(nvar,sizeof(double *));
  for (vindex=0;vindex<nvar;vindex++){ 
    dim *= rangera[vindex];
    rangefactorsra[vindex] = (double *) calloc(rangera[vindex],sizeof(double));
    for (rindex=0;rindex<rangera[vindex];rindex++){
      if (multorsumra[vindex]>0){
	ftemp = scaletypera[vindex]>0 ? rindex : rindex-(double)(rangera[vindex]-1.0)/2.0;
	rangefactorsra[vindex][rindex] = pow(factorra[vindex],ftemp);}
      else /* if (multorsumra[vindex]<=0) */{ 
	ftemp = scaletypera[vindex]>0 ? rindex : rindex-(double)(rangera[vindex]-1.0)/2.0;
	rangefactorsra[vindex][rindex] = factorra[vindex]*ftemp;}}}
  if (verbose>1){ 
    printf("rangefactorsra are...\n");
    for (vindex=0;vindex<nvar;vindex++){ 
      printf("\t variable %d = %s of type %s -- ",vindex,keyvarra[vindex],keyvartypera[vindex]);
      for (rindex=0;rindex<rangera[vindex];rindex++){
	printf("%0.3f ",rangefactorsra[vindex][rindex]);}
      printf("\n");}}
  sprintf(outputname,"%s.in",inputname);
  if ((ip=fopen(outputname,"r"))==NULL){ printf("cannot read file %s.in, aborting\n",inputname); exit(EXIT_SUCCESS);}
  indexra = (int *) calloc(nvar,sizeof(int));
  fp = (FILE **) calloc(dim,sizeof(FILE *));
  chartemp = (char *) calloc(charsize,sizeof(char));
  for (dindex=0;dindex<dim;dindex++){
    indextract(dindex,nvar,rangera,indexra); schain(chartemp,nvar,indexra); sprintf(outputname,"%s%s.in",inputname,chartemp);
    if (verbose){ printf("creating output_file %s\n",outputname);}
    if ((fp[dindex]=fopen(outputname,"w"))==NULL){ printf("cannot create file %s, aborting\n",outputname); exit(EXIT_SUCCESS);}}
  do{
    fscanf(ip,"%[^=]",vname);fscanf(ip,"%s",equals);fscanf(ip,"%c",equals);
    if (strcmp(vname,"GLOBAL_STRING")==0){
      fscanf(ip,"%[^,;]",equals);
      if (verbose){ printf("found GLOBAL_STRING variable %s= %s, recording and changing to:\n",vname,equals);}
      if (verbose){ printf("%s= %s%s;\n",vname,equals,"_?");}
      sprintf(gstring,"%s",equals);
      for (dindex=0;dindex<dim;dindex++){
	indextract(dindex,nvar,rangera,indexra); schain(chartemp,nvar,indexra);
	fprintf(fp[dindex],"%s= %s%s;",vname,equals,chartemp);}
      fscanf(ip,"%c",equals);}
    else /* if (strcmp(vname,"GLOBAL_STRING")!=0) */{
      found_flag=-1; for (vindex=0;vindex<nvar;vindex++){ if (strcmp(vname,keyvarra[vindex])==0){ found_flag = vindex;}}
      if (found_flag<0){
	if (verbose){ printf("found regular variable %s= ",vname);}
	for (dindex=0;dindex<dim;dindex++){ fprintf(fp[dindex],"%s= ",vname);}
	do{
	  fscanf(ip,"%[^,;]",equals);
	  if (verbose){ printf("%s",equals);}
	  for (dindex=0;dindex<dim;dindex++){ fprintf(fp[dindex],"%s",equals);}
	  fscanf(ip,"%c",comma_vs_semicolon);
	  if (strcmp(comma_vs_semicolon,",")==0){
	    if (verbose){ printf(",");}
	    for (dindex=0;dindex<dim;dindex++){ fprintf(fp[dindex],",",equals);}}
	  else if (strcmp(comma_vs_semicolon,";")==0){
	    if (verbose){ printf(";\n");}
	    for (dindex=0;dindex<dim;dindex++){ fprintf(fp[dindex],";\n");}}
	  else{ printf(" error! poor divider %d\n",comma_vs_semicolon);}}
	while (strcmp(comma_vs_semicolon,",")==0);}
      else /* if (found_flag>=0) */{
	if (verbose){ printf("found our keyname variable %s at variable index %d\n",keyvarra[found_flag],found_flag);}
	if (verbose){ printf("%s= ",keyvarra[found_flag]);}
	for (dindex=0;dindex<dim;dindex++){ fprintf(fp[dindex],"%s= ",keyvarra[found_flag]);}
	varindex=0;
	do{
	  fscanf(ip,"%[^,;]",equals);
	  if (verbose){ printf("%s",equals);}
	  match_index_and_name=0; matching_vindex=-1;
	  for (vindex=0;vindex<nvar;vindex++){ 
	    if (varindex==keyvarindexra[vindex] && strcmp(vname,keyvarra[vindex])==0){ matching_vindex=vindex;}}
	  match_index_and_name = (matching_vindex!=-1);
	  for (dindex=0;dindex<dim;dindex++){
	    if (match_index_and_name){ 
	      if (verbose){ printf("changing...");}
	      indextract(dindex,nvar,rangera,indexra);
	      if (0){}
	      else if (strcmp(keyvartypera[matching_vindex],"int")==0){
		if (multorsumra[matching_vindex]>0){
		  fprintf(fp[dindex],"%d",(int)(rangefactorsra[matching_vindex][indexra[matching_vindex]]*atoi(equals)));}
		else /* if (multorsumra[matching_vindex]<=0) */{
		  fprintf(fp[dindex],"%d",(int)(rangefactorsra[matching_vindex][indexra[matching_vindex]]+atoi(equals)));}}
	      else if (strcmp(keyvartypera[matching_vindex],"double")==0){
		if (multorsumra[matching_vindex]>0){
		  fprintf(fp[dindex],"%0.16lf",(double)(rangefactorsra[matching_vindex][indexra[matching_vindex]]*atof(equals)));}
		else /* if (multorsumra[matching_vindex]<=0) */{
		  fprintf(fp[dindex],"%0.16lf",(double)(rangefactorsra[matching_vindex][indexra[matching_vindex]]+atof(equals)));}}
	      else{
		printf(" error! wackty type %s\n",keyvartypera[matching_vindex]);
		fprintf(fp[dindex],"%s",equals);}}
	    else /* if (!match_index_and_name) */{
	      fprintf(fp[dindex],"%s",equals);}}
	  fscanf(ip,"%c",comma_vs_semicolon);
	  if (strcmp(comma_vs_semicolon,",")==0){
	    if (verbose){ printf(",");}
	    for (dindex=0;dindex<dim;dindex++){ fprintf(fp[dindex],",",equals);}}
	  else if (strcmp(comma_vs_semicolon,";")==0){
	    if (verbose){ printf(";\n");}
	    for (dindex=0;dindex<dim;dindex++){ fprintf(fp[dindex],";\n");}}
	  else{ printf(" error! poor divider %d\n",comma_vs_semicolon);}
	  varindex += 1;}
	while (strcmp(comma_vs_semicolon,",")==0);}
      fscanf(ip,"%c",equals);}}
  while (strcmp(vname,"END")!=0);
  fclose(ip);
  if (verbose){ printf("now create helper run files\n");}
  sprintf(outputname,"""rm"" helper_run_%s_*.sh;\n",inputname); system(outputname);
  for (dindex=0;dindex<dim;dindex++){
    runfilenumber = (perprocessor>0 ? dindex/perprocessor : 0);
    sprintf(outputname,"helper_run_%s_%d.sh",inputname,runfilenumber);
    if ((hp=fopen(outputname,"a"))==NULL){ printf("cannot create file %s, aborting\n",outputname); exit(EXIT_SUCCESS);}
    indextract(dindex,nvar,rangera,indexra); schain(chartemp,nvar,indexra);
    for (rindex=0;rindex<maxrecord;rindex++){ 
      if (verbose>1){ printf(" writing: nohup nice ./%s < %s%s.in %s\n",programname,inputname,chartemp,seqorsim>0?";":"&");}
      fprintf(hp,"nohup nice ./%s < %s%s.in %s\n",programname,inputname,chartemp,seqorsim>0?";":"&");}
    fclose(hp);}
  if (verbose){ printf("now create helper unpack file\n");}
  sprintf(outputname,"""rm"" helper_unpack_%s.sh;\n",inputname); system(outputname); 
  sprintf(outputname,"helper_unpack_%s.sh",inputname);
  if ((hp=fopen(outputname,"w"))==NULL){ printf("cannot create file %s, aborting\n",outputname); exit(EXIT_SUCCESS);}
  fprintf(hp,"gunzip *.gz;\n");  
  for (dindex=0;dindex<dim;dindex++){
    indextract(dindex,nvar,rangera,indexra); schain(chartemp,nvar,indexra);
    for (rindex=0;rindex<maxrecord;rindex++){ fprintf(hp,"tar xvf %s%s_%d_bundle.tar;\n",gstring,chartemp,rindex+(grecord ? 1 : 0));}}
  fclose(hp);
  for (dindex=0;dindex<dim;dindex++){ fclose(fp[dindex]);} free(fp);
  for (vindex=0;vindex<nvar;vindex++){ free(rangefactorsra[vindex]);} free(rangefactorsra);
  for (vindex=0;vindex<nvar;vindex++){ free(keyvarra[vindex]);} free(keyvarra);
  for (vindex=0;vindex<nvar;vindex++){ free(keyvartypera[vindex]);} free(keyvartypera);
  free(keyvarindexra);
  free(indexra);
  free(scaletypera);
  free(rangera);
  free(factorra);
  return 1;
}

  
    
