#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int main(int argc, char **argv)
{
  int verbose=2;
  char keyvarname[64],keyvartype[16],inputname[64],outputname[128],vname[64],equals[64];
  int range=0,arange=0,rindex=0;
  double factor=2;
  double *rangefactors=NULL;
  int still_have_options=argc-1;
  FILE *ip=NULL;
  FILE **fp=NULL;
  if (verbose>1){ printf("starting inputcreate...\n");}
  if (still_have_options==0){ printf("-Finputname[char*], -Vkeyvarname[char*], -Tkeyvartype{'int','double','float'}, -Rrange[int] -Dfactor[double]\n"); exit(EXIT_SUCCESS);}
  while (still_have_options){
    switch(getopt(argc,argv,"F:f:V:v:T:t:R:r:D:d:")){
    case 'F': case 'f': sprintf(inputname,optarg); break;
    case 'V': case 'v': sprintf(keyvarname,optarg); break;
    case 'T': case 't': sprintf(keyvartype,optarg); break;
    case 'R': case 'r': range = atoi(optarg); break;
    case 'D': case 'd': factor = atof(optarg); break;
    default: printf("-Finputname, -Vkeyvarname, -Tkeyvartype, -Rrange -Dfactor\n"); exit(EXIT_SUCCESS); break;}
    still_have_options--;}
  if (verbose){ printf("inputfile %s, variable %s of type %s, range %d\n",inputname,keyvarname,keyvartype,range);}
  arange = abs(range);
  rangefactors = (double *) calloc(arange,sizeof(double));
  for (rindex=0;rindex<arange;rindex++){ rangefactors[rindex] = pow(factor,(range > 0 ? 1 : -1)*(rindex+1));}
  sprintf(outputname,"%s.in",inputname);
  if ((ip=fopen(outputname,"r"))==NULL){ printf("cannot read file %s.in, aborting\n",inputname); exit(EXIT_SUCCESS);}
  fp = (FILE **) calloc(arange,sizeof(FILE *));
  for (rindex=0;rindex<arange;rindex++){
    sprintf(outputname,"%s_%s_%d.in",inputname,keyvarname,rindex);
    if (verbose){ printf("creating output_file %s\n",outputname);}
    if ((fp[rindex]=fopen(outputname,"w"))==NULL){ printf("cannot create file %s, aborting\n",outputname); exit(EXIT_SUCCESS);}}
  do{
    fscanf(ip,"%[^=]",vname);fscanf(ip,"%s",equals);
    if (strcmp(vname,"GLOBAL_STRING")==0){
      fscanf(ip,"%s",equals);
      if (verbose){ printf("found GLOBAL_STRING variable %s= %s, changing!\n",vname,equals);}
      for (rindex=0;rindex<arange;rindex++){
	sprintf(outputname,"%s_%s_%d",equals,keyvarname,rindex);
	fprintf(fp[rindex],"%s= %s\n",vname,outputname);}
      fscanf(ip,"%c",equals);}
    else if (strcmp(vname,keyvarname)!=0){
      fscanf(ip,"%s",equals);
      if (verbose){ printf("found regular variable %s= %s...\n",vname,equals);}
      for (rindex=0;rindex<arange;rindex++){ fprintf(fp[rindex],"%s= %s\n",vname,equals);}
      fscanf(ip,"%c",equals);}
    else /* if (strcmp(vname,keyvarname)==0) */{
      fscanf(ip,"%s",equals);
      if (verbose){ printf("found our key %s variable %s= %s changing!\n",keyvarname,vname,equals);}
      if (0){}
      else if (strcmp(keyvartype,"int")==0){
	for (rindex=0;rindex<arange;rindex++){ fprintf(fp[rindex],"%s= %d\n",vname,(int)(rangefactors[rindex]*atoi(equals)));}}
      else if (strcmp(keyvartype,"double")==0 || strcmp(keyvartype,"float")==0){
	for (rindex=0;rindex<arange;rindex++){ fprintf(fp[rindex],"%s= %0.16lf\n",vname,(double)(rangefactors[rindex]*atof(equals)));}}
      else{ 
	printf("wacky type %s not yet consistent\n",keyvartype);
	for (rindex=0;rindex<arange;rindex++){ fprintf(fp[rindex],"%s= %s\n",vname,equals);}}
      fscanf(ip,"%c",equals);}}
  while (strcmp(vname,"END")!=0);
  fclose(ip);
  for (rindex=0;rindex<arange;rindex++){ fclose(fp[rindex]);}
  free(fp);
  free(rangefactors);
}

  
    
