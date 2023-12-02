#include <stdio.h>     /* for printf */
#include <stdlib.h>    /* for exit */
#include <math.h>
#include <getopt.h>

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

int GRAYSCALE=0;

void rgb2hsv(double r,double g,double b,double *h,double *s,double *v)
{
  /* shamelessly stolen from http://cs.haifa.ac.il/hagit/courses/ist/Lectures/Demos/ColorApplet2/t_convert.html
     r,g,b values are from 0 to 1
     h = [0,360], s = [0,1], v = [0,1]
     if s == 0, then h = -1 (undefined) */
  double min=0,max=0,delta=0;
  min = minimum(r,minimum(g,b));
  max = maximum(r,maximum(g,b));
  *v = max;
  delta = max - min;
  if (max!=0){ *s = delta/max;}
  else{ *s = 0; *h = -1; return;}
  if (r==max){ /* between yellow and magenta */ *h = (g - b)/delta;}
  else if (g==max){ /* between cyan and yellow */ *h = 2+(b-r)/delta;}
  else /* if (b==max) */{ /* between magenta and cyan */ *h = 4+(r-g)/delta;}
  *h *= 60; /* degrees */
  if (*h<0){ *h+=360;}
}

void hsv2rgb(double h,double s,double v,double *r,double *g,double *b)
{
  int i;
  double f=0,p=0,q=0,t=0;
  if (s==0){ /* greyscale */ *r=v;*g=v;*b=v; return;}
  h /= 60; i=floor(h); /* six sectors */
  f = h-i;
  p = v*(1-s);
  q = v*(1-s*f);
  t = v*(1-s*(1-f));
  switch (i){
  case 0: *r=v;*g=t;*b=p; break;
  case 1: *r=q;*g=v;*b=p; break;
  case 2: *r=p;*g=v;*b=t; break;
  case 3: *r=p;*g=q;*b=v; break;
  case 4: *r=t;*g=p;*b=v; break;
  default: /* case 5 */ *r=v;*g=p;*b=q; break;}
}

void colorscale(int coloringtype,double val,double valmax,double valmin,double *rcolor,double *gcolor,double *bcolor)
{
  double gamma=(coloringtype==5?1.0:0.5);
  int norm_flag=1;
  double norm=1;
  double v=0;
  double h=0;
  switch (coloringtype){
  case 0: /* greyscale coloring */
    if (valmax<=valmin){ valmax=valmin+1;}
    v = (val-valmin)/(valmax-valmin); v = maximum(0,minimum(1,v));
    *rcolor = v; *gcolor = v; *bcolor = v;
    break;
  case 2: /* cyclical coloring */
    if (valmax<=valmin){ valmax=valmin+1;}
    v = 3*((val-valmin)/(valmax-valmin)-0.5);
    if (v < -0.5 || v > 0.5){ *rcolor = sqrt(fabs(v)-0.5);} else{ *rcolor=0;}
    if (v > -0.5){ *gcolor = sqrt(1-fabs(v-0.5));} else{ *gcolor=0;}
    if (v < 0.5){ *bcolor = sqrt(1-fabs(v+0.5));} else{ *bcolor=0;}
    break;
  case 3: /* IC IS EC ES coloring */
    switch ((int)val){
    case 0: /* IC */ *rcolor=0.25; *gcolor=0.75; *bcolor=1.00; break;
    case 1: /* IS */ *rcolor=0.00; *gcolor=0.00; *bcolor=1.00; break;
    case 2: /* EC */ *rcolor=1.00; *gcolor=0.75; *bcolor=0.25; break;
    case 3: /* ES */ *rcolor=1.00; *gcolor=0.00; *bcolor=0.00; break;
    default: *rcolor=1; *gcolor=1; *bcolor=1; break;}
    break;
  case 4: /* coloring with gamma */
    if (valmax==valmin){ valmax++;}
    v = (val-valmin)/(valmax-valmin);
    if (v<0){ v=0;}
    if (v>1){ v=1;}
    v = pow(v,gamma);
    if (!GRAYSCALE){
      if (v<0.5){
	v = asin(2*v);
	*rcolor = 0;
	*gcolor = sin(v);
	*bcolor = cos(v);}
      else if (v>=0.5){
	v = asin(2*(v-0.5));
	*rcolor = sin(v);
	*gcolor = cos(v);
	*bcolor = 0;}}
    else if (GRAYSCALE){
      *rcolor = v;
      *gcolor = v;
      *bcolor = v;}
    break;
  case 5: /* super coloring with blue-cyan-green-yellow-red and gamma value */
    if (valmax==valmin){ valmax++;}
    v = (val-valmin)/(valmax-valmin);
    if (v<0){ v=0;}
    if (v>1){ v=1;}
    v = pow(v,gamma);
    if (!GRAYSCALE){
      if (v>=0.0 && v<0.25){
	*rcolor = 0;
	*gcolor = (v-0.0)/0.25;
	*bcolor = 1;}
      if (v>=0.25 && v<0.5){
	*rcolor = 0;
	*gcolor = 1;
	*bcolor = 1-(v-0.25)/0.25;}
      if (v>=0.5 && v<0.75){
	*rcolor = (v-0.5)/0.25;
	*gcolor = 1;
	*bcolor = 0;}
      if (v>=0.75 && v<=1){
	*rcolor = 1;
	*gcolor = 1-(v-0.75)/0.25;
	*bcolor = 0;}
      if (norm_flag){ norm = sqrt(pow(*rcolor,2)+pow(*gcolor,2)+pow(*bcolor,2)); *rcolor /= norm; *gcolor /= norm; *bcolor /= norm;}}
    else if (GRAYSCALE){
      *rcolor = v;
      *gcolor = v;
      *bcolor = v;}
    break;
  case 6: /* use hsv2rgb with h in [0,360] */
    if (valmax<=valmin){ valmax=valmin+1;}
    h = 360*(val-valmin)/(valmax-valmin);
    hsv2rgb(h,1,1,rcolor,gcolor,bcolor);
    break;
  case 7: /* use hsv2rgb with h in [240,-60] */
    if (valmax<=valmin){ valmax=valmin+1;}
    h = 300*(valmax-val)/(valmax-valmin)-60;
    hsv2rgb(h,1,1,rcolor,gcolor,bcolor);
    break;
  default:
    if (valmax==valmin){ valmax++;}
    v = (val-valmin)/(valmax-valmin);
    if (v<0){ v=0;}
    if (v>1){ v=1;}
    if (!GRAYSCALE){
      if (v<0.5){
	v = asin(2*v);
	*rcolor = 0;
	*gcolor = sin(v);
	*bcolor = cos(v);}
      else if (v>=0.5){
	v = asin(2*(v-0.5));
	*rcolor = sin(v);
	*gcolor = cos(v);
	*bcolor = 0;}}
    else if (GRAYSCALE){
      *rcolor = v;
      *gcolor = v;
      *bcolor = v;}
    break;}
}

int main (int argc, char **argv) 
{
  int nr=0;
  double rcolor=0,gcolor=0,bcolor=0;
  int rcode=0,gcode=0,bcode=0;
  int xord=0,yord=0,side=32;
  char filename[512],filename_txt[1024],filename_fig[1024];
  FILE *fp=NULL;
  int cc=0;
  for (cc=0;cc<=7;cc++){
    sprintf(filename,"preamble_colorscale_%d",cc);sprintf(filename_txt,"./%s.txt",filename);sprintf(filename_fig,"./%s.fig");
    if ((fp=fopen(filename_txt,"w"))==NULL){ printf(" %% error, can't open %s\n",filename_txt); exit(EXIT_SUCCESS);};
    fprintf(fp,"#FIG 3.2\\nLandscape\\nCenter\\nInches\\nLetter\\n100.00\\nSingle\\n-2\\n1200 2\\n");
    for (nr=0;nr<512;nr++){
      colorscale(cc,nr,511,0,&rcolor,&gcolor,&bcolor);
      rcode = crop((int)floor(256*rcolor),0,255);
      gcode = crop((int)floor(256*gcolor),0,255);
      bcode = crop((int)floor(256*bcolor),0,255);
      fprintf(fp,"0 %d #%0.2x%0.2x%0.2x\\n",nr+32,rcode,gcode,bcode);}
    fclose(fp);
    if ((fp=fopen(filename_fig,"w"))==NULL){ printf(" %% error, can't open %s\n",filename_fig); exit(EXIT_SUCCESS);};
    fprintf(fp,"#FIG 3.2\nLandscape\nCenter\nInches\nLetter\n100.00\nSingle\n-2\n1200 2\n");
    for (nr=0;nr<512;nr++){
      colorscale(cc,nr,511,0,&rcolor,&gcolor,&bcolor);
      rcode = crop((int)floor(256*rcolor),0,255);
      gcode = crop((int)floor(256*gcolor),0,255);
      bcode = crop((int)floor(256*bcolor),0,255);
      fprintf(fp,"0 %d #%0.2x%0.2x%0.2x\n",nr+32,rcode,gcode,bcode);}
    for (nr=0;nr<512;nr++){
      fprintf(fp,"2 2 0 0 -1 %d 50 0 %d 0.000 0 0 -1 0 0 5\n",/*colorcode*/nr+32,/*fill*/20);
      xord = side*(nr/16); yord = side*(nr%16);
      fprintf(fp,"\t %d %d %d %d %d %d %d %d %d %d\n",xord,yord,xord+side,yord,xord+side,yord+side,xord,yord+side,xord,yord);}
    fclose(fp);}
}
    
