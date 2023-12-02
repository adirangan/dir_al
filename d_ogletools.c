/* #ifndef NOTGRAPH */

/* int WritePPMFile(const char *, GLubyte *,int,int,int); */
/* int DumpWindow(const char *,int,int); */
/* void ftexto(float,float,float,char *); */
/* GLvoid InitGL(GLsizei,GLsizei); */
/* GLvoid ReSizeGLScene(GLsizei,GLsizei); */
/* GLvoid DrawNra(double,double,double); */
/* GLvoid Drawlgn(double,double,double,double); */
/* GLvoid glbox(double,double,double,double,double,double); */
/* GLvoid Drawcolorbar(int,double,double,double); */
/* GLvoid Drawmenu(double,double,double); */
/* void keyPressed(unsigned char,int,int);  */
/* void specialKeyPressed(int,int,int);  */
/* GLvoid DrawGLScene(GLvoid); */

/* #endif /\* NOTGRAPH *\/ */

#ifndef NOTGRAPH

int WritePPMFile(const char *filename, GLubyte *buffer, int win_width, int win_height, int mode)
{
  int i, j, k, q;
  unsigned char *ibuffer;
  FILE *fp;
  int pixelSize = mode==GL_RGBA?4:3;
  ibuffer = (unsigned char *) malloc(win_width*win_height*RGB3);
  if ( (fp = fopen(filename, "wb")) == NULL ) {
    printf(" %% Warning: cannot open %s in WritePPMFile\n", filename);
    return 0;}
  fprintf(fp, "P6\n# CREATOR: glReadPixel()\n%d %d\n%d\n",win_width, win_height, UCHAR_MAX);
  q = 0;
  for (i = 0; i < win_height; i++)
    for (j = 0; j < win_width; j++)
      for (k = 0; k < RGB3; k++)
	ibuffer[q++] = (unsigned char)*(buffer + (pixelSize*((win_height-1-i)*win_width+j)+k));
  fwrite(ibuffer, sizeof(unsigned char), RGB3*win_width*win_height, fp);
  fclose(fp);
  free(ibuffer);
  printf(" %% wrote file %s, (%d x %d pixels, %d bytes)\n",filename,win_width, win_height, RGB3*win_width*win_height);
  return 1;
}

int DumpWindow(const char *filename, int win_width, int win_height) 
{
  GLubyte *buffer;
  int result;
  buffer = (GLubyte *) malloc(win_width*win_height*RGBA);
  /* read window contents from color buffer with glReadPixels */
  glFinish();
  glReadPixels(0, 0, win_width, win_height,GL_RGBA, GL_UNSIGNED_BYTE, buffer);
  result = WritePPMFile( filename, buffer, win_width, win_height,GL_RGBA );
  free(buffer);
  return result;
}

void ftexto(float x, float y, float z, char *text)
{
  /* thanks to alex */
  char *p;
  glRasterPos3f(x,y,z);
  for (p = text; *p; p++){ glutBitmapCharacter(GLUT_BITMAP_8_BY_13, *p);}
}

GLvoid InitGL(GLsizei Width, GLsizei Height)
{
  /* A general OpenGL initialization function.  Sets all of the initial parameters. */
  glClearColor(0.0f, 0.0f, 0.0f, 0.0f);	// This Will Clear The Background Color To Black
  glClearDepth(1.0);				// Enables Clearing Of The Depth Buffer
  glDepthFunc(GL_LESS);			// The Type Of Depth Test To Do
  glEnable(GL_DEPTH_TEST);			// Enables Depth Testing
  glShadeModel(GL_SMOOTH);			// Enables Smooth Color Shading
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();				// Reset The Projection Matrix
  gluPerspective(45.0f,(GLfloat)Width/(GLfloat)Height,0.1f,100.0f);	// Calculate The Aspect Ratio Of The Window
  glMatrixMode(GL_MODELVIEW);
}

GLvoid ReSizeGLScene(GLsizei Width, GLsizei Height)
{
  /* The function called when our window is resized */
  if (Height==0){ Height=1;} /* don't divide by zero */
  glViewport(0, 0, Width, Height); /* reset viewport & perspective */
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(45.0f,(GLfloat)Width/(GLfloat)Height,0.1f,100.0f);
  glMatrixMode(GL_MODELVIEW);
}

GLvoid loglogra(void *ra,char *type,int length,double max,double min,double xside,double yside,double xoffset,double yoffset,double rcolor,double gcolor,double bcolor)
{
  int nr=0;
  double max2=0,min2=0,mean=0,stdev=0;
  double xord=0,yord=0;
  double *dra=NULL;
  int *ira=NULL;
  if (strcmp(type,"double")==0){ 
    dra = (double *)ra;
    if (max>min){ max2=max;min2=min;} 
    else if (max<min){ stats(type,dra,length,&max2,&min2,&mean,&stdev);} 
    else if (max==min){ stats(type,dra,length,&max2,&min2,&mean,&stdev); min2=0;}
    if (min2==0){ min2=0.01*max2;}
    glColor3f(rcolor,gcolor,bcolor);
    glBegin(GL_LINE_STRIP);
    for (nr=0;nr<length;nr++){
      if (dra[nr]>0){
	yord = (log(dra[nr])-log(min2))/(log(max2)-log(min2));
	xord = log((double)nr+0.5)/log((double)length);
	glVertex3f(xord*xside+xoffset,yord*yside+yoffset,0);}}
    glEnd();}
  else if (strcmp(type,"int")==0){ 
    ira = (int *)ra;
    if (max>min){ max2=max;min2=min;} 
    else if (max<min){ stats(type,ira,length,&max2,&min2,&mean,&stdev);} 
    else if (max==min){ stats(type,ira,length,&max2,&min2,&mean,&stdev); min2=0;}
    if (min2==0){ min2=0.01*max2;}
    glColor3f(rcolor,gcolor,bcolor);
    glBegin(GL_LINE_STRIP);
    for (nr=0;nr<length;nr++){
      if (ira[nr]>0){
	yord = (log(ira[nr])-log(min2))/(log(max2)-log(min2));
	xord = log((double)nr+0.5)/log((double)length);
	glVertex3f(xord*xside+xoffset,yord*yside+yoffset,0);}}
    glEnd();}
}

GLvoid Plotra(void *ra,char *type,int length,double max,double min,double xside,double yside,double xoffset,double yoffset,double rcolor,double gcolor,double bcolor)
{
  int nr=0;
  double max2=0,min2=0,mean=0,stdev=0;
  double xord=0,yord=0;
  double *dra=NULL;
  int *ira=NULL;
  if (strcmp(type,"double")==0){ 
    dra = (double *)ra;
    if (max>min){ max2=max;min2=min;} 
    else if (max<min){ stats(type,dra,length,&max2,&min2,&mean,&stdev);} 
    else if (max==min){ stats(type,dra,length,&max2,&min2,&mean,&stdev); min2=0;}
    glColor3f(rcolor,gcolor,bcolor);
    glBegin(GL_LINE_STRIP);
    for (nr=0;nr<length;nr++){
      yord = (dra[nr]-min2)/(max2-min2);
      xord = ((double)nr+0.5)/(double)length;
      glVertex3f(xord*xside+xoffset,yord*yside+yoffset,0);}
    glEnd();}
  else if (strcmp(type,"int")==0){ 
    ira = (int *)ra;
    if (max>min){ max2=max;min2=min;} 
    else if (max<min){ stats(type,ira,length,&max2,&min2,&mean,&stdev);} 
    else if (max==min){ stats(type,ira,length,&max2,&min2,&mean,&stdev); min2=0;}
    glColor3f(rcolor,gcolor,bcolor);
    glBegin(GL_LINE_STRIP);
    for (nr=0;nr<length;nr++){
      yord = (ira[nr]-min2)/(max2-min2);
      xord = ((double)nr+0.5)/(double)length;
      glVertex3f(xord*xside+xoffset,yord*yside+yoffset,0);}
    glEnd();}
}

GLvoid Drawra(void *ra,char *type,int rows,int cols,int coloringtype,int stduse,int gsuse,char *text,double xside,double yside,double xoffset,double yoffset)
{
  /* draws a standard *ra */
  char text2[64];
  int i=0,j=0;
  double varmax=0,varmin=0,varmean=0,varstd=0,vMAX=0,vMIN=0;
  double *dra=NULL,*dra2=NULL;
  int *ira=NULL,*ira2=NULL;
  double side=maximum(xside,yside),bigside=maximum(rows,cols);
  double xord=0,yord=0;
  double rcolor=0,gcolor=0,bcolor=0; 
  int clear_flag=0;
/*   int i2=0,j2=0,i3=0,j3=0; */
/*   double mask=0,dist2=0,angle2=0,dist3=0,angle3=0; */
/*   int highlight2=0,highlight3=0; */
  if (strcmp(type,"double")==0){ 
    dra = (double *)ra;
    if (gsuse){ dra2 = spacesmear(dra,rows,cols,gsuse*GLOBAL_SPACE_SMOOTHER); clear_flag=1;} else{ dra2 = dra; clear_flag=0;}
    stats("double",dra2,rows*cols,&varmax,&varmin,&varmean,&varstd);
    if (stduse){ vMAX=varmean+STD_VIEW*varstd; vMIN=varmean-STD_VIEW*varstd;} else{ vMAX=varmax;vMIN=varmin;}
    for (i=0;i<rows;i++){ for (j=0;j<cols;j++){
      xord = (j+0.5)/bigside; yord = -(i+0.5)/bigside;
      colorscale(coloringtype,dra2[i+j*rows],vMAX,vMIN,&rcolor,&gcolor,&bcolor);
/*       /\* fix later *\/ */
/*       i2 = periodize(SYSTEM_ROW_SIZE/4+ARBOR_DIA-i,-SYSTEM_ROW_SIZE/2,SYSTEM_ROW_SIZE/2); */
/*       j2 = periodize(SYSTEM_COL_SIZE/4+ARBOR_DIA-j,-SYSTEM_COL_SIZE/2,SYSTEM_COL_SIZE/2); */
/*       i2 = SYSTEM_ROW_SIZE/4+ARBOR_DIA-i; */
/*       j2 = SYSTEM_COL_SIZE/4+ARBOR_DIA-j; */
/*       dist2 = sqrt(pow(i2,2)+pow(j2,2)); */
/*       angle2 = atan2(j2,i2)-PI/6; */
/*       highlight2=0; */
/*       i3 = periodize(SYSTEM_ROW_SIZE/4+ARBOR_DIA-i,-SYSTEM_ROW_SIZE/2,SYSTEM_ROW_SIZE/2); */
/*       j3 = periodize(SYSTEM_COL_SIZE/4-ARBOR_DIA-j,-SYSTEM_COL_SIZE/2,SYSTEM_COL_SIZE/2); */
/*       i3 = SYSTEM_ROW_SIZE/4+ARBOR_DIA-i; */
/*       j3 = SYSTEM_COL_SIZE/4-ARBOR_DIA-j; */
/*       dist3 = sqrt(pow(i3,2)+pow(j3,2)); */
/*       angle3 = atan2(j3,i3)+PI/6; */
/*       highlight3=0; */
/*       if (fabs(periodize(angle2,-PI/3,PI/3))<(5*PI/18) && dist2<0.9*ARBOR_DIA){ highlight2=1;} */
/*       if (fabs(periodize(angle3,-PI/3,PI/3))<(5*PI/18) && dist3<0.9*ARBOR_DIA){ highlight3=1;} */
/*       mask = 1-0.8*exp(-maximum(dist2,dist3)/8/ARBOR_DIA); */
/*       if (highlight2 || highlight3){ mask=1;} */
/*       mask = maximum(0,minimum(1,mask)); */
/*       rcolor *= mask; gcolor *= mask; bcolor *= mask;       */
      glbox(xord*xside+xoffset,yord*yside+yoffset,side/bigside,rcolor,gcolor,bcolor);}}
    if (clear_flag){ tfree(dra2); dra2=NULL;}}
  else if (strcmp(type,"int")==0){ 
    ira = (int *)ra; ira2 = ira; clear_flag=0;
    stats("int",ira2,rows*cols,&varmax,&varmin,&varmean,&varstd);
    if (stduse){ vMAX=varmean+STD_VIEW*varstd; vMIN=varmean-STD_VIEW*varstd;} else{ vMAX=varmax;vMIN=varmin;}
    for (i=0;i<rows;i++){ for (j=0;j<cols;j++){
      xord = j/bigside; yord = -i/bigside;
      colorscale(coloringtype,ira2[i+j*rows],vMAX,vMIN,&rcolor,&gcolor,&bcolor);
      glbox(xord*xside+xoffset,yord*yside+yoffset,side/bigside,rcolor,gcolor,bcolor);}}
    if (clear_flag){ tfree(ira2); ira2=NULL;}}
  glColor3f(1,1,1);
  if (text!=NULL){ 
    if (maximum(fabs(varmin),fabs(varmax))<0.005){ sprintf(text2,"%s-[%0.3e..%0.3e]",text,varmin,varmax);} 
    else{ sprintf(text2,"%s-[%0.2f..%0.2f]",text,varmin,varmax);}} 
  else{ 
    if (maximum(fabs(varmin),fabs(varmax))<0.005){ sprintf(text2,"[%0.3e..%0.3e]",varmin,varmax);} 
    else{ sprintf(text2,"[%0.2f..%0.2f]",varmin,varmax);}} 
  ftexto(xoffset,yoffset+0.1,0,text2);
}

GLvoid DrawLR(double side,double xoffset,double yoffset)
{
  /* draws the lr_kernels and puts them to use */
  char text[32];
  int nr=0,nc=0,tab=0;
  struct neuron *n=NULL;
  struct neuronarray *Nra=GLOBAL_Nra;
  int rows=Nra->rows,cols=Nra->cols;
  double bigside = maximum(rows,cols);
  double xord=0,yord=0,temp=0;
  double rcolor=0,gcolor=0,bcolor=0; 
  fftw_complex **lrx=GLOBAL_LRKERNEL_X;
  int nl=DRAW_FLAG%NLRKERNELS;
  int *gs = GLOBAL_SHUFFLES;
  for (nr=0;nr<rows;nr++){ for (nc=0;nc<cols;nc++){
    n=nget(Nra,nr,nc);
    tab = n->ang + n->pierow*PIE_ROW_DIA*PIE_COL_DIA + n->piecol*PIE_ROW_DIA*PIE_COL_DIA*NPIEROWS;
    xord = nc/bigside; yord = -nr/bigside;
    temp = lrx[nl][n->pierow + n->piecol*NPIEROWS + gs[tab]*NPIEROWS*NPIECOLS][0];
    // temp = lrx[nl][n->pierow + n->piecol*NPIEROWS + n->ang*NPIEROWS*NPIECOLS][0];
    colorscale(7,temp,1.0/LR_PCONNECT,0,&rcolor,&gcolor,&bcolor);
    glbox(xord*side+xoffset,yord*side+yoffset,side/bigside,rcolor,gcolor,bcolor);}}
  glColor3f(1,1,1);
  sprintf(text,"lrk");
  ftexto(xoffset,yoffset,0,text);
}

GLvoid DrawNra(double side,double xoffset,double yoffset)
{
  /*
     DRAW_FLAG==
     -1: warning
     0: pierow
     1: piecol
     2: inpierow
     3: inpiecol
     4: t2s
     5: inputrate
     6: rad
     7: ang
     8: lgnangle
     9: lgnphase
     10: 
     11: sog
     12: spiketime
     13:
     14:
     15:
     16: n->V
     17: n->sA
     18: n->sN
     19: n->sG
     20: n->VS
  */
  char type[32],text[32];
  int varp=0,i=0,j=0,DdiI=0; /* (D)ouble array, (d)ouble, (i)nteger or (I)nteger array */
  double *vra=NULL;
  struct neuron *n=NULL;
  struct neuronarray *Nra=GLOBAL_Nra;
  int rows=Nra->rows,cols=Nra->cols;
  int coloringtype=0, stduse=0, gsuse=0;
  int vra_clear_flag=0;
  n = nget(Nra,0,0);
  switch (periodize(DRAW_FLAG,-1,21)){
  case -1: sprintf(text,"warning"); varp = (int)&(n->warning)-(int)n; DdiI=+1; coloringtype=7; stduse=0; gsuse=0; break;
  case 0: sprintf(text,"pierow"); varp = (int)&(n->pierow)-(int)n; DdiI=-1; coloringtype=7; stduse=0; gsuse=0; break;
  case 1: sprintf(text,"piecol"); varp = (int)&(n->piecol)-(int)n; DdiI=-1; coloringtype=7; stduse=0; gsuse=0; break;
  case 2: sprintf(text,"inpierow"); varp = (int)&(n->inpierow)-(int)n; DdiI=-1; coloringtype=7; stduse=0; gsuse=0; break;
  case 3: sprintf(text,"inpiecol"); varp = (int)&(n->inpiecol)-(int)n; DdiI=-1; coloringtype=7; stduse=0; gsuse=0; break;
  case 4: sprintf(text,"t2s"); varp = (int)&(Nra->t2sra)-(int)Nra; DdiI=+2; coloringtype=3; stduse=0; gsuse=0; break;
  case 5: sprintf(text,"inputrate"); varp = (int)&(n->inputrate)-(int)n; DdiI=+1; coloringtype=7; stduse=1; gsuse=1; break;
  case 6: sprintf(text,"rad"); varp = (int)&(n->rad)-(int)n; DdiI=+1; coloringtype=7; stduse=0; gsuse=0; break;
  case 7: sprintf(text,"ang"); varp = (int)&(n->ang)-(int)n; DdiI=-1; coloringtype=6; stduse=0; gsuse=0; break;
  case 8: sprintf(text,"lgnangle"); varp = (int)&(n->lgnangle)-(int)n; DdiI=+1; coloringtype=6; stduse=0; gsuse=0; break;
  case 9: sprintf(text,"lgnphase"); varp = (int)&(n->lgnphase)-(int)n; DdiI=+1; coloringtype=6; stduse=0; gsuse=0; break;
  case 11: sprintf(text,"sog"); varp = (int)&(n->sog)-(int)n; DdiI=-1; coloringtype=1; stduse=0; gsuse=0; break;
  case 12: sprintf(text,"spiketime"); varp = (int)&(n->spiketime)-(int)n; DdiI=+1; coloringtype=7; stduse=0; gsuse=0; break;
  case 16: sprintf(text,"Vra"); varp = (int)&(Nra->Vra)-(int)Nra; DdiI=+2; coloringtype=7; stduse=1; gsuse=1; break;
  case 17: sprintf(text,"sAra"); varp = (int)&(Nra->sAra)-(int)Nra; DdiI=+2; coloringtype=7; stduse=1; gsuse=1; break;
  case 18: sprintf(text,"sNra"); varp = (int)&(Nra->sNra)-(int)Nra; DdiI=+2; coloringtype=7; stduse=1; gsuse=1; break;
  case 19: sprintf(text,"sGra"); varp = (int)&(Nra->sGra)-(int)Nra; DdiI=+2; coloringtype=7; stduse=1; gsuse=1; break;
  case 20: sprintf(text,"VSra"); varp = (int)&(Nra->VSra)-(int)Nra; DdiI=+2; coloringtype=7; stduse=1; gsuse=1; break;
  default: sprintf(text,"lgnangle"); varp = (int)&(n->lgnangle)-(int)n; DdiI=+1; coloringtype=6; stduse=0; gsuse=0; break;}
  if (abs(DdiI)>1){ 
    if (DdiI==+2){ vra_clear_flag=0;vra = *(double **)((int)Nra+varp); sprintf(type,"double");}}
  else{ 
    vra_clear_flag=1; vra = (double *) tcalloc(rows*cols,sizeof(double));
    if (DdiI==+1){
      for (i=0;i<rows;i++){ for (j=0;j<cols;j++){ 
	n = nget(Nra,i,j);
	vra[i+j*rows] = *(double *)((int)n+varp);}}
      sprintf(type,"double");}
    else if (DdiI==-1){
      for (i=0;i<rows;i++){ for (j=0;j<cols;j++){ 
	n = nget(Nra,i,j);
	vra[i+j*rows] = (double)(*(int *)((int)n+varp));}}
      sprintf(type,"double");}}
  Drawra(vra,type,rows,cols,coloringtype,stduse,gsuse,text,side,side,xoffset,yoffset);
  if (vra_clear_flag){ tfree(vra);}
}

GLvoid Drawlgn(double xside,double yside,double xoffset,double yoffset)
{
  int verbose=1;
  char text[64],s[32];
  double rcolor=0,gcolor=0,bcolor=0;
  struct lgn *p=GLOBAL_LGN;
  int lag=0,tab=1,leap=0,rshift=0,cshift=0;
  int angle_flag=0,phase_flag=0;
  double *ra=NULL,max=0,min=0,mean=0,std=0;
  int r=0,c=0,nr2=0,nc2=0;
  double xord=0,yord=0;
  double temp=0,*temp2=0;
  double side = maximum(xside,yside);
  int bigside = maximum(p->rows,p->cols);
  if (verbose){ temp2 = (double *) tmalloc(sizeof(double)*p->cols*p->rows);}
  if (DRAW_FLAG>0){
    switch ((DRAW_FLAG-1)%5){
    case 1: ra = p->onshape; max = STD_VIEW*(1); min = STD_VIEW*(-1); lag = 0; tab = 1; rshift = p->rows/2; cshift = p->cols/2; 
      sprintf(s,"onshape"); break;
    case 2: ra = p->onin; max = STD_VIEW*(1); min = STD_VIEW*(-1); lag = 0; tab = 1; rshift = 0; cshift = 0;
      sprintf(s,"oninput"); break; 
    case 3: ra = p->onrate; max = STD_VIEW*p->peakrate; min = -STD_VIEW*p->peakrate; lag = 4; tab = 5; rshift = 0; cshift = 0;
      sprintf(s,"onrate"); break;
    case 4: ra = p->offrate; max = STD_VIEW*p->peakrate; min = -STD_VIEW*p->peakrate; lag = 4; tab = 5; rshift = 0; cshift = 0;
      sprintf(s,"offrate"); break;
    default: ra = p->pnm; max = STD_VIEW*(1); min = STD_VIEW*(-1); lag = 0; tab = 1; 
      sprintf(s,"pnm"); break; }}
  else{
    leap = -DRAW_FLAG/(NSLICES*NPHASES);
    //angle_flag = (-DRAW_FLAG-leap*NSLICES*NPHASES)%NSLICES; phase_flag = (-DRAW_FLAG-leap*NSLICES*NPHASES)/NSLICES; 
    phase_flag = (-DRAW_FLAG-leap*NSLICES*NPHASES)%NPHASES; angle_flag = (-DRAW_FLAG-leap*NSLICES*NPHASES)/NPHASES; 
    switch (leap%3){
    case 0: 
      ra = p->angleshape[angle_flag + phase_flag*NSLICES]; 
      max = STD_VIEW*(1); min = STD_VIEW*(-1); lag = 0; tab = 1; rshift = p->rows/2; cshift = p->cols/2; 
      sprintf(s,"angle_%d_phase_%d_shape",angle_flag,phase_flag);
      break;
    case 1:
      ra = p->anglein[angle_flag + phase_flag*NSLICES]; 
      max = STD_VIEW*(1); min = STD_VIEW*(-1); lag = 0; tab = 1; rshift = 0; cshift = 0;
      sprintf(s,"angle_%d_phase_%d_in",angle_flag,phase_flag);
      break;
    case 2:
      ra = p->anglerate[angle_flag + phase_flag*NSLICES]; 
      max = STD_VIEW*p->peakrate; min = -STD_VIEW*p->peakrate; lag = 4; tab = 5; rshift = 0; cshift = 0;
      sprintf(s,"angle_%d_phase_%d_rate",angle_flag,phase_flag);
      break;
    default: break;}}
  for (r=0;r<p->rows;r++){ for (c=0;c<p->cols;c++){
    nr2=periodize(r+rshift,0,p->rows); nc2 = periodize(c+cshift,0,p->cols);
    temp = ra[lag + nr2*tab + nc2*p->rows*tab];
    if (verbose){ temp2[nr2+nc2*p->rows] = temp;}
    colorscale(1,temp,max,min,&rcolor,&gcolor,&bcolor);
    xord = (double)(c+0.5)/(double)bigside; yord = (double)(r+0.5)/(double)bigside;
    glbox(xord*xside+xoffset,yord*yside+yoffset,side/bigside,rcolor,gcolor,bcolor);}}
  glColor3f(1,1,1);
  sprintf(text,"%s, min = %0.1f to %0.1f = max",s,min,max);
  ftexto(xoffset,yoffset-.1,0,text);
  if (verbose){ 
    stats("double",temp2,p->rows*p->cols,&max,&min,&mean,&std);
    sprintf(text,"max=%0.3f, min = %0.3f, mean = %0.3f, std = %0.3f",max,min,mean,std);
    ftexto(xoffset,yoffset-.2,0,text);
    tfree(temp2);}
}

GLvoid Drawrtc(double xside,double yside,double xoffset,double yoffset)
{
  char text[32],text2[64];
  int nl=0,na=0,nt2s=0;
  struct rtc *r=GLOBAL_RTC;
  double max=0,min=0,xshift=0;
  double rcolor=0,gcolor=0,bcolor=0;
  double *ra=NULL;
  double *vra=NULL;
  int varp=0;
  int drawvsplot=DRAW_FLAG2%2;
  ra = (double *) tcalloc(r->nangles*r->length,sizeof(double));
  if (DRAW_FLAG<0){
    for (na=0;na<r->nangles;na++){ for (nl=0;nl<r->length;nl++){
      ra[na+nl*r->nangles] = r->anglera[nl]==na;}}
    sprintf(text,"angle"); 
    Drawra(ra,"double",r->nangles,r->length,7,1,0,text,xside,yside,xoffset,yoffset-1.5*yside);
    switch (periodize(DRAW_FLAG,0,7)){
    case 1: sprintf(text,"V"); varp = (int)&(r->Vra)-(int)r; break;
    case 2: sprintf(text,"sA"); varp = (int)&(r->sAra)-(int)r; break;
    case 3: sprintf(text,"sN"); varp = (int)&(r->sNra)-(int)r; break;
    case 4: sprintf(text,"sG"); varp = (int)&(r->sGra)-(int)r; break;
    case 5: sprintf(text,"VS"); varp = (int)&(r->VSra)-(int)r; break;
    default: sprintf(text,"m"); varp = (int)&(r->mra)-(int)r; break;}
    vra = *(double **)((int)r+varp);
    for (nt2s=0;nt2s<4;nt2s++){ 
      sprintf(text2,"%s%d",text,nt2s);
      Drawra(vra+nt2s*r->nangles*r->length,"double",r->nangles,r->length,7,1,0,text2,xside,yside,xoffset,yoffset+1.5*nt2s*yside);}}
  else{
    for (nt2s=0;nt2s<4;nt2s++){
      colorscale(3,nt2s,1,0,&rcolor,&gcolor,&bcolor);
      xshift=0;
      for (na=0;na<r->nangles-1;na++){ for (nl=0;nl<r->length;nl++){ 
	ra[na+nl*(r->nangles-1)] = log(r->rtcmra[na+nl*r->nangles+nt2s*r->length*r->nangles]/r->rtcmra[r->nangles-1+nl*r->nangles+nt2s*r->length*r->nangles]);}}
      max=1; min=-1;
      if (drawvsplot==0){ for (nl=0;nl<r->length;nl++){ Plotra(&(ra[0+nl*(r->nangles-1)]),"double",r->nangles-1,max,min,xside,yside,xoffset+xshift*xside,yoffset+nl*yside,rcolor,gcolor,bcolor);}}
      else{ sprintf(text,"R(m)%d",nt2s); Drawra(ra,"double",r->nangles-1,r->length,7,0,0,text,xside,yside,xoffset+xshift*1.25*xside,yoffset+nt2s*1.25*yside);}
      sprintf(text,"R(m)"); ftexto(xoffset+xshift*xside,yoffset-.2,0,text);
      xshift=1;
      for (na=0;na<r->nangles;na++){ for (nl=0;nl<r->length;nl++){ 
	ra[na+nl*r->nangles] = r->rtcmra[na+nl*r->nangles+nt2s*r->length*r->nangles]/maximum(1,r->rtcra[nl]);}}
      stats("double",ra,r->nangles*r->length,&max,&min,NULL,NULL); min=0;
      if (drawvsplot==0){ for (nl=0;nl<r->length;nl++){ Plotra(&(ra[0+nl*r->nangles]),"double",r->nangles,max,min,xside,yside,xoffset+xshift*xside,yoffset+nl*yside,rcolor,gcolor,bcolor);}}
      else{ sprintf(text,"m%d",nt2s); Drawra(ra,"double",r->nangles,r->length,7,0,0,text,xside,yside,xoffset+xshift*1.25*xside,yoffset+nt2s*1.25*yside);}
      sprintf(text,"m"); ftexto(xoffset+xshift*xside,yoffset-.2,0,text);
      xshift=2;
      for (na=0;na<r->nangles;na++){ for (nl=0;nl<r->length;nl++){ 
	ra[na+nl*r->nangles] = r->rtcVra[na+nl*r->nangles+nt2s*r->length*r->nangles]/maximum(1,r->rtcra[nl]);}}
      stats("double",ra,r->nangles*r->length,&max,&min,NULL,NULL);
      if (drawvsplot==0){ for (nl=0;nl<r->length;nl++){ Plotra(&(ra[0+nl*r->nangles]),"double",r->nangles,max,min,xside,yside,xoffset+xshift*xside,yoffset+nl*yside,rcolor,gcolor,bcolor);}}
      else{ sprintf(text,"V%d",nt2s); Drawra(ra,"double",r->nangles,r->length,7,0,0,text,xside,yside,xoffset+xshift*1.25*xside,yoffset+nt2s*1.25*yside);}
      sprintf(text,"V"); ftexto(xoffset+xshift*xside,yoffset-.2,0,text);
      xshift=3;
      for (na=0;na<r->nangles;na++){ for (nl=0;nl<r->length;nl++){ 
	ra[na+nl*r->nangles] = r->rtcsAra[na+nl*r->nangles+nt2s*r->length*r->nangles]/maximum(1,r->rtcra[nl]);}}
      stats("double",ra,r->nangles*r->length,&max,&min,NULL,NULL); min=0;
      if (drawvsplot==0){ for (nl=0;nl<r->length;nl++){ Plotra(&(ra[0+nl*r->nangles]),"double",r->nangles,max,min,xside,yside,xoffset+xshift*xside,yoffset+nl*yside,rcolor,gcolor,bcolor);}}
      else{ sprintf(text,"sA%d",nt2s); Drawra(ra,"double",r->nangles,r->length,7,0,0,text,xside,yside,xoffset+xshift*1.25*xside,yoffset+nt2s*1.25*yside);}
      sprintf(text,"sA"); ftexto(xoffset+xshift*xside,yoffset-.2,0,text);
      xshift=4;
      for (na=0;na<r->nangles;na++){ for (nl=0;nl<r->length;nl++){ 
	ra[na+nl*r->nangles] = r->rtcsNra[na+nl*r->nangles+nt2s*r->length*r->nangles]/maximum(1,r->rtcra[nl]);}}
      stats("double",ra,r->nangles*r->length,&max,&min,NULL,NULL); min=0;
      if (drawvsplot==0){ for (nl=0;nl<r->length;nl++){ Plotra(&(ra[0+nl*r->nangles]),"double",r->nangles,max,min,xside,yside,xoffset+xshift*xside,yoffset+nl*yside,rcolor,gcolor,bcolor);}}
      else{ sprintf(text,"sN%d",nt2s); Drawra(ra,"double",r->nangles,r->length,7,0,0,text,xside,yside,xoffset+xshift*1.25*xside,yoffset+nt2s*1.25*yside);}
      sprintf(text,"sN"); ftexto(xoffset+xshift*xside,yoffset-.2,0,text);
      xshift=5;
      for (na=0;na<r->nangles;na++){ for (nl=0;nl<r->length;nl++){ 
	ra[na+nl*r->nangles] = r->rtcsGra[na+nl*r->nangles+nt2s*r->length*r->nangles]/maximum(1,r->rtcra[nl]);}}
      stats("double",ra,r->nangles*r->length,&max,&min,NULL,NULL); min=0;
      if (drawvsplot==0){ for (nl=0;nl<r->length;nl++){ Plotra(&(ra[0+nl*r->nangles]),"double",r->nangles,max,min,xside,yside,xoffset+xshift*xside,yoffset+nl*yside,rcolor,gcolor,bcolor);}}
      else{ sprintf(text,"sG%d",nt2s); Drawra(ra,"double",r->nangles,r->length,7,0,0,text,xside,yside,xoffset+xshift*1.25*xside,yoffset+nt2s*1.25*yside);}
      sprintf(text,"sG"); ftexto(xoffset+xshift*xside,yoffset-.2,0,text);
      xshift=6;
      for (na=0;na<r->nangles;na++){ for (nl=0;nl<r->length;nl++){ 
	ra[na+nl*r->nangles] = r->rtcVSra[na+nl*r->nangles+nt2s*r->length*r->nangles]/maximum(1,r->rtcra[nl]);}}
      stats("double",ra,r->nangles*r->length,&max,&min,NULL,NULL); min=0;
      if (drawvsplot==0){ for (nl=0;nl<r->length;nl++){ Plotra(&(ra[0+nl*r->nangles]),"double",r->nangles,max,min,xside,yside,xoffset+xshift*xside,yoffset+nl*yside,rcolor,gcolor,bcolor);}}
      else{ sprintf(text,"VS%d",nt2s); Drawra(ra,"double",r->nangles,r->length,7,0,0,text,xside,yside,xoffset+xshift*1.25*xside,yoffset+nt2s*1.25*yside);}
      sprintf(text,"VS"); ftexto(xoffset+xshift*xside,yoffset-.2,0,text);}
    glColor3f(1,1,1);
    for (nl=0;nl<r->length;nl++){ sprintf(text,"%d",(int)r->rtcra[nl]); ftexto(xoffset+.2*nl,yoffset-.4,0,text);}}
  tfree(ra);
}

GLvoid Drawstra(struct strobe **stra,int ralength,double side,double xoffset,double yoffset,double max,double min,double rcolor,double gcolor,double bcolor,int drawvsplot)
{
  int sttab=stra[0]->tab,stlength=stra[0]->length,stcycle_bother=stra[0]->cycle_bother,stcyclenum=stra[0]->cyclenum;
  int na=0,nt=0,nt2=0;
  double xord=0,yord=0;
  double max2=max,min2=min,stdscale=0;
  double *meanra=NULL,*stdevra=NULL,mean=0,stdev=0,*ra=NULL;
  int coloringtype = rcolor >= 0 ? 0 : (int)gcolor;
  double rcolor2=0,gcolor2=0,bcolor2=0;
  int automatic_color=0;
  if (drawvsplot==0){
    if (max<=min){ 
      stdscale=17;
      meanra = (double *) tcalloc(ralength,sizeof(double));
      stdevra = (double *) tcalloc(ralength,sizeof(double));
      for (na=0;na<ralength;na++){ stats("double",stra[na]->data,stlength,NULL,NULL,meanra+na,stdevra+na);}
      stats("double",meanra,ralength,NULL,NULL,&mean,NULL);
      stats("double",stdevra,ralength,NULL,NULL,&stdev,NULL);
      tfree(meanra);tfree(stdevra);
      max2 = mean+stdscale*STD_VIEW*stdev; min2 = mean-stdscale*STD_VIEW*stdev;}
    for (na=0;na<ralength;na++){ 
      switch (coloringtype){
      case -1: /* automatic coloring */ automatic_color=1; break;
      case 0: /* specified coloring */ rcolor2=rcolor;gcolor2=gcolor;bcolor2=bcolor; break;
      case 1: colorscale(7,na,ralength-1,0,&rcolor2,&gcolor2,&bcolor2); break;
      case 2: colorscale(6,na,ralength-1,0,&rcolor2,&gcolor2,&bcolor2); break;
      case 3: colorscale(3,na,3,1,&rcolor2,&gcolor2,&bcolor2); break;
      default: rcolor2=rcolor;gcolor2=gcolor;bcolor2=bcolor; break;}
      for (nt=sttab;nt<sttab+stlength;nt++){
	nt2 = periodize(nt,0,stlength);
	xord = (double)(nt-sttab+0.5)/(double)stlength;
	if (automatic_color){ 
	  colorscale(7,stra[na]->data[nt2],max2,min2,&rcolor2,&gcolor2,&bcolor2);
	  yord = (double)(-na+0.5)/(double)ralength;
	  glbox(xord*side+xoffset,yord*side+yoffset,side/(double)minimum(ralength,stlength),rcolor2,gcolor2,bcolor2);}
	else{
	  yord = (stra[na]->data[nt2]-min2)/(max2-min2);
	  glbox(xord*side+xoffset,yord*side+yoffset,side/(double)minimum(100,stlength),rcolor2,gcolor2,bcolor2);}}}
    if (stcycle_bother && stcyclenum>0){
      if (max<=min){
	stdscale=17;
	meanra = (double *) tcalloc(ralength,sizeof(double));
	stdevra = (double *) tcalloc(ralength,sizeof(double));
	for (na=0;na<ralength;na++){ stats("double",stra[na]->cycledata,stlength,NULL,NULL,meanra+na,stdevra+na);}
	stats("double",meanra,ralength,NULL,NULL,&mean,NULL);
	stats("double",stdevra,ralength,NULL,NULL,&stdev,NULL);
	tfree(meanra);tfree(stdevra);
	max2 = mean+stdscale*STD_VIEW*stdev; min2 = mean-stdscale*STD_VIEW*stdev;}
      else{ max2 = max*stcyclenum; min2 = min*stcyclenum;}
      for (na=0;na<ralength;na++){ 
	switch (coloringtype){
	case -1: automatic_color=1; break;
	case 0: rcolor2=rcolor;gcolor2=gcolor;bcolor2=bcolor; break;
	case 1: colorscale(7,na,ralength-1,0,&rcolor2,&gcolor2,&bcolor2); break;
	case 2: colorscale(6,na,ralength-1,0,&rcolor2,&gcolor2,&bcolor2); break;
	case 3: colorscale(3,na,3,1,&rcolor2,&gcolor2,&bcolor2); break;
	default: rcolor2=rcolor;gcolor2=gcolor;bcolor2=bcolor; break;}
	for (nt=0;nt<stlength;nt++){
	  xord = (double)(nt+0.5)/(double)stlength;
	  if (automatic_color){ 
	    colorscale(7,stra[na]->cycledata[nt]/stcyclenum,max2,min2,&rcolor2,&gcolor2,&bcolor2);
	    yord = (double)(-na+0.5)/(double)ralength;
	    glbox(xord*side+xoffset,yord*side+yoffset,side/(double)minimum(ralength,stlength),rcolor2,gcolor2,bcolor2);}
	  else{
	    yord = (stra[na]->cycledata[nt]-min2)/(max2-min2);
	    glbox(xord*side+xoffset,yord*side+yoffset+side,side/(double)minimum(100,stlength),rcolor2,gcolor2,bcolor2);}}}}}
  else if (drawvsplot==1){
    ra = (double *) tcalloc(ralength*stlength,sizeof(double));
    for (na=0;na<ralength;na++){ for (nt=sttab;nt<sttab+stlength;nt++){
      nt2 = periodize(nt,0,stlength);
      ra[na + (nt-sttab)*ralength] = stra[na]->data[nt2];}}
    Drawra(ra,"double",ralength,stlength,7,1,0,NULL,side,side,xoffset,yoffset);
    if (stcycle_bother && stcyclenum>0){
      for (na=0;na<ralength;na++){ for (nt=0;nt<stlength;nt++){
	ra[na + nt*ralength] = stra[na]->cycledata[nt]/stcyclenum;}}
      Drawra(ra,"double",ralength,stlength,7,1,0,NULL,side,side,xoffset,yoffset+side);}
    tfree(ra);}
}

GLvoid Drawstrobetrace(double side,double xoffset,double yoffset)
{
  char text[32],text2[32];
  int local_draw_flag = abs(DRAW_FLAG)%7;
  int na=0,nt=0,nt2=0,nr=0,nc=0,nd=0,d=0,nt2s=0;
  int rows=0,cols=0;
  double rcolor=0,gcolor=0,bcolor=0;
  double xord=0,yord=0,max=0,min=0,*maxra=NULL,*minra=NULL,max_cycle=0,min_cycle=0;
  struct strobetrace *st=GLOBAL_STROBETRACE;
  double *ra=NULL;
  double sumsn=0,sumsnsn=0,sumvs=0,sumvsvs=0,sumsnvs=0,sumt=0;
  int drawvsplot = DRAW_FLAG2%2;
  if (DRAW_FLAG<0){
    rows = st->pcspiesacross*PIE_ROW_DIA;
    cols = st->pcspiestall*PIE_COL_DIA;
    Drawra(&(st->pcsra[0+0*rows+(abs(DRAW_FLAG)%(st->nangles))*rows*cols]),"double",rows,cols,7,1,1,"pcs",side,side,xoffset,yoffset);
    rows = NARBORS_TALL;
    cols = NARBORS_WIDE;
    ra = (double *) tcalloc(rows*cols,sizeof(double));
    for (nr=0;nr<rows;nr++){ for (nc=0;nc<cols;nc++){
      sumsn = st->snvstc[0 + nr*5 + nc*rows*5]; sumsnsn = st->snvstc[1 + nr*5 + nc*rows*5];
      sumvs = st->snvstc[2 + nr*5 + nc*rows*5]; sumvsvs = st->snvstc[3 + nr*5 + nc*rows*5];
      sumsnvs = st->snvstc[4 + nr*5 + nc*rows*5]; sumt = st->total_time;
      ra[nr + nc*rows] = (sumt*sumsnvs - sumsn*sumvs)/sqrt((sumt*sumsnsn-sumsn*sumsn)*(sumt*sumvsvs-sumvs*sumvs));}}
    Drawra(ra,"double",rows,cols,7,1,0,"snvstc",side,side,xoffset-side,yoffset);
    tfree(ra);
    for (nd=st->avalanche->logdmin;nd<st->avalanche->logdmax;nd++){ d=(int)pow(2,nd); for (nt2s=0;nt2s<5;nt2s++){
      colorscale(3,nt2s,4,0,&rcolor,&gcolor,&bcolor);
      Plotra(st->avalanche->mhist[nd+nt2s*st->avalanche->logdmax],"double",(st->avalanche->Nra->rows/d+1)*(st->avalanche->Nra->cols/d+1),0,0,side/2,side/2,xoffset+side+nd*side/2,yoffset-side+nt2s*side/2,rcolor,gcolor,bcolor);}}}
  else /*if (DRAW_FLAG>=0) */{
    switch (local_draw_flag){
    case 0: /* single neuron */
      Drawstra(&(st->Vst),1,side,xoffset,yoffset,VOLTAGE_THRESHOLD,VOLTAGE_IN,1,1,1,0);
      Drawstra(&(st->VSst),1,side,xoffset,yoffset,VOLTAGE_THRESHOLD,VOLTAGE_IN,1,1,0,0);
      Drawstra(&(st->sAst),1,side,xoffset,yoffset,0,0,1,0,0,0);
      Drawstra(&(st->sNst),1,side,xoffset,yoffset,0,0,0,1,0,0);
      Drawstra(&(st->sGst),1,side,xoffset,yoffset,0,0,0,0,1,0);
      break;
    case 1:
      sprintf(text,"m"); Drawstra(st->mstra,4,side,xoffset,yoffset,0,0,-1,3,0,drawvsplot); break;
    case 2: /* Voltage of angular subpopulations */
      sprintf(text,"VS"); Drawstra(st->VSstra,st->nangles,side,xoffset,yoffset,VOLTAGE_THRESHOLD,VOLTAGE_IN,-1,2,0,drawvsplot); break;
    case 3: /* Conductance of angular subpopulations */
      sprintf(text,"sA"); Drawstra(st->sAstra,st->nangles,side,xoffset,yoffset,0,0,-1,2,0,drawvsplot); break;
    case 4: 
      sprintf(text,"sN"); Drawstra(st->sNstra,st->nangles,side,xoffset,yoffset,0,0,-1,2,0,drawvsplot); break;
    case 5: 
      sprintf(text,"sG"); Drawstra(st->sGstra,st->nangles,side,xoffset,yoffset,0,0,-1,2,0,drawvsplot); break;
    case 6: 
      sprintf(text,"pattern");
      /* pcs similarity & firing rates */
      maxra = (double *) tcalloc(st->nangles,sizeof(double)); minra = (double *) tcalloc(st->nangles,sizeof(double));
      for (na=0;na<st->nangles;na++){ stats("double",st->pmstra[na]->data,st->length,maxra+na,minra+na,NULL,NULL);};
      stats("double",maxra,st->nangles,&max,NULL,NULL,NULL); stats("double",minra,st->nangles,NULL,&min,NULL,NULL);
      tfree(maxra);tfree(minra);
      if (st->cycle_bother){
	maxra = (double *) tcalloc(st->nangles,sizeof(double)); minra = (double *) tcalloc(st->nangles,sizeof(double));
	for (na=0;na<st->nangles;na++){ stats("double",st->pmstra[na]->cycledata,st->length,maxra+na,minra+na,NULL,NULL);};
	stats("double",maxra,st->nangles,&max_cycle,NULL,NULL,NULL); stats("double",minra,st->nangles,NULL,&min_cycle,NULL,NULL);
	tfree(maxra);tfree(minra);}
      for (na=0;na<st->nangles;na++){
	for (nt=st->tst->tab;nt<st->tst->tab+st->length;nt++){
	  /* and faith here */
	  nt2 = periodize(nt,0,st->length);
	  xord = (double)(nt-st->tst->tab+0.5)/(double)st->length;
	  yord = (double)(na+0.5)/(double)st->nangles;
	  colorscale(7,st->patternstra_VS[na]->data[nt2],0.5,-0.5,&rcolor,&gcolor,&bcolor);
	  glbox(xord*side+xoffset,yord*side+yoffset-1*side,side/(double)st->nangles,rcolor,gcolor,bcolor);
	  colorscale(7,st->patternstra_sN[na]->data[nt2],0.5,-0.5,&rcolor,&gcolor,&bcolor);
	  glbox(xord*side+xoffset,yord*side+yoffset+0*side,side/(double)st->nangles,rcolor,gcolor,bcolor);
	  colorscale(7,st->pmstra[na]->data[nt2],max,min,&rcolor,&gcolor,&bcolor);
	  glbox(xord*side+xoffset,yord*side+yoffset+1*side,side/(double)st->nangles,rcolor,gcolor,bcolor);
	  if (st->cycle_bother){
	    colorscale(7,st->patternstra_VS[na]->cycledata[nt-st->tst->tab],0.5,-0.5,&rcolor,&gcolor,&bcolor);
	    glbox(xord*side+xoffset,yord*side+yoffset+2*side,side/(double)st->nangles,rcolor,gcolor,bcolor);
	    colorscale(7,st->patternstra_sN[na]->cycledata[nt-st->tst->tab],0.5,-0.5,&rcolor,&gcolor,&bcolor);
	    glbox(xord*side+xoffset,yord*side+yoffset+3*side,side/(double)st->nangles,rcolor,gcolor,bcolor);
	    colorscale(7,st->pmstra[na]->cycledata[nt-st->tst->tab],max_cycle,min_cycle,&rcolor,&gcolor,&bcolor);
	    glbox(xord*side+xoffset,yord*side+yoffset+4*side,side/(double)st->nangles,rcolor,gcolor,bcolor);}}}
      glColor3f(1,1,1);
      sprintf(text2,"max=%0.1e,min=%0.1e",max,min);
      ftexto(xoffset,yoffset+2*side,0,text2);
      sprintf(text,"patternc");
      Drawstra(&(st->patternc),1,side,xoffset+1*side,yoffset+0*side,+1,-1.0,1,1,1,0);
      glBegin(GL_LINES);
      glColor3f(1,1,1);
      xord = 0; yord = 0; glVertex3f(xord*side+xoffset+1*side,yord*side+yoffset+0*side,0);
      xord = 1; glVertex3f(xord*side+xoffset+1*side,yord*side+yoffset+0*side,0);
      xord = 0; yord = 0.5; glVertex3f(xord*side+xoffset+1*side,yord*side+yoffset+0*side,0);
      xord = 1; glVertex3f(xord*side+xoffset+1*side,yord*side+yoffset+0*side,0);
      xord = 0; yord = 1; glVertex3f(xord*side+xoffset+1*side,yord*side+yoffset+0*side,0);
      xord = 1; glVertex3f(xord*side+xoffset+1*side,yord*side+yoffset+0*side,0);
      xord = 0; yord = 0; glVertex3f(xord*side+xoffset+1*side,yord*side+yoffset+0*side,0);
      yord = 1; glVertex3f(xord*side+xoffset+1*side,yord*side+yoffset+0*side,0);
      xord = 1; yord = 0; glVertex3f(xord*side+xoffset+1*side,yord*side+yoffset+0*side,0);
      yord = 1; glVertex3f(xord*side+xoffset+1*side,yord*side+yoffset+0*side,0);
      glEnd();
      break;
    default: break;}
    if (!drawvsplot){
      glBegin(GL_LINES);
      glColor3f(1,1,1);
      xord = 0; yord = (VOLTAGE_REST-VOLTAGE_IN)/(VOLTAGE_THRESHOLD-VOLTAGE_IN); glVertex3f(xord*side+xoffset,yord*side+yoffset,0);
      xord = 1; glVertex3f(xord*side+xoffset,yord*side+yoffset,0);
      glColor3f(1,0,0);
      xord = 0; yord = (VOLTAGE_THRESHOLD-VOLTAGE_IN)/(VOLTAGE_THRESHOLD-VOLTAGE_IN); glVertex3f(xord*side+xoffset,yord*side+yoffset,0);
      xord = 1; glVertex3f(xord*side+xoffset,yord*side+yoffset,0);
      glColor3f(0,0,1);
      xord = 0; yord = (VOLTAGE_RESET-VOLTAGE_IN)/(VOLTAGE_THRESHOLD-VOLTAGE_IN); glVertex3f(xord*side+xoffset,yord*side+yoffset,0);
      xord = 1; glVertex3f(xord*side+xoffset,yord*side+yoffset,0);
      if (st->cycle_bother){
	glColor3f(1,1,1);
	xord = 0; yord = (VOLTAGE_REST-VOLTAGE_IN)/(VOLTAGE_THRESHOLD-VOLTAGE_IN); glVertex3f(xord*side+xoffset,yord*side+yoffset+side,0);
	xord = 1; glVertex3f(xord*side+xoffset,yord*side+yoffset+side,0);
	glColor3f(1,0,0);
	xord = 0; yord = (VOLTAGE_THRESHOLD-VOLTAGE_IN)/(VOLTAGE_THRESHOLD-VOLTAGE_IN); glVertex3f(xord*side+xoffset,yord*side+yoffset+side,0);
	xord = 1; glVertex3f(xord*side+xoffset,yord*side+yoffset+side,0);
	glColor3f(0,0,1);
	xord = 0; yord = (VOLTAGE_RESET-VOLTAGE_IN)/(VOLTAGE_THRESHOLD-VOLTAGE_IN); glVertex3f(xord*side+xoffset,yord*side+yoffset+side,0);
	xord = 1; glVertex3f(xord*side+xoffset,yord*side+yoffset+side,0);}
      glEnd();}
    glColor3f(1,1,1);
    ftexto(xoffset,yoffset+side,0,text);}
}

GLvoid Drawtuningcurve(double xside,double yside,double xoffset,double yoffset)
{
  struct tuningcurve *T = GLOBAL_TUNINGCURVE;
  int nr=0,nt2s=0;
  double rcolor=0,gcolor=0,bcolor=0;
  char text[32];
  for (nr=0;nr<T->nradius;nr++){
    for (nt2s=0;nt2s<4;nt2s++){
      colorscale(3,nt2s,1,0,&rcolor,&gcolor,&bcolor);
      Plotra(&(T->mra[0 + nr*T->nangles + nt2s*T->nangles*T->nradius]),"double",T->nangles,0,0,xside,yside,xoffset+0*xside,yoffset+nr*yside,rcolor,gcolor,bcolor);
      Plotra(&(T->Vra[0 + nr*T->nangles + nt2s*T->nangles*T->nradius]),"double",T->nangles,-1,0,xside,yside,xoffset+1*xside,yoffset+nr*yside,rcolor,gcolor,bcolor);
      Plotra(&(T->sAra[0 + nr*T->nangles + nt2s*T->nangles*T->nradius]),"double",T->nangles,0,0,xside,yside,xoffset+2*xside,yoffset+nr*yside,rcolor,gcolor,bcolor);
      Plotra(&(T->sNra[0 + nr*T->nangles + nt2s*T->nangles*T->nradius]),"double",T->nangles,0,0,xside,yside,xoffset+3*xside,yoffset+nr*yside,rcolor,gcolor,bcolor);
      Plotra(&(T->sGra[0 + nr*T->nangles + nt2s*T->nangles*T->nradius]),"double",T->nangles,0,0,xside,yside,xoffset+4*xside,yoffset+nr*yside,rcolor,gcolor,bcolor);
      Plotra(&(T->VSra[0 + nr*T->nangles + nt2s*T->nangles*T->nradius]),"double",T->nangles,-1,0,xside,yside,xoffset+5*xside,yoffset+nr*yside,rcolor,gcolor,bcolor);}}
  glColor3f(1,1,1);
  sprintf(text,"mra");
  ftexto(xoffset+0*xside,yoffset-0.1*yside,0,text);
  sprintf(text,"Vra");
  ftexto(xoffset+1*xside,yoffset-0.1*yside,0,text);
  sprintf(text,"sAra");
  ftexto(xoffset+2*xside,yoffset-0.1*yside,0,text);
  sprintf(text,"sNra");
  ftexto(xoffset+3*xside,yoffset-0.1*yside,0,text);
  sprintf(text,"sGra");
  ftexto(xoffset+4*xside,yoffset-0.1*yside,0,text);
  sprintf(text,"VSra");
  ftexto(xoffset+5*xside,yoffset-0.1*yside,0,text);
}

GLvoid Drawlmitri(double xside,double yside,double xoffset,double yoffset)
{
  struct lmitri *lt=GLOBAL_LMITRI;
  int nt2s=0,tab=0;
  char text[32];
  for (nt2s=0;nt2s<5;nt2s++){
    tab = 0 + 0*lt->space_length + nt2s*lt->space_length*lt->time_length;
    sprintf(text,"m%d",nt2s);
    Drawra(&(lt->mra[tab]),"double",lt->space_length,lt->time_length,7,1,0,text,xside,yside,xoffset+0*xside*1.25,yoffset+nt2s*yside*1.25);
    sprintf(text,"V%d",nt2s);
    Drawra(&(lt->Vra[tab]),"double",lt->space_length,lt->time_length,7,1,0,text,xside,yside,xoffset+1*xside*1.25,yoffset+nt2s*yside*1.25);
    sprintf(text,"sA%d",nt2s);
    Drawra(&(lt->sAra[tab]),"double",lt->space_length,lt->time_length,7,1,0,text,xside,yside,xoffset+2*xside*1.25,yoffset+nt2s*yside*1.25);
    sprintf(text,"sN%d",nt2s);
    Drawra(&(lt->sNra[tab]),"double",lt->space_length,lt->time_length,7,1,0,text,xside,yside,xoffset+3*xside*1.25,yoffset+nt2s*yside*1.25);
    sprintf(text,"sG%d",nt2s);
    Drawra(&(lt->sGra[tab]),"double",lt->space_length,lt->time_length,7,1,0,text,xside,yside,xoffset+4*xside*1.25,yoffset+nt2s*yside*1.25);
    sprintf(text,"VS%d",nt2s);
    Drawra(&(lt->VSra[tab]),"double",lt->space_length,lt->time_length,7,1,0,text,xside,yside,xoffset+5*xside*1.25,yoffset+nt2s*yside*1.25);}}

GLvoid Draweventra(struct ptree *p,double side,double xoffset,double yoffset)
{
  struct llist *L=NULL;
  struct litem *l=NULL;
  struct region *r=NULL;
  int nt=0,tab=0;
  double xord=0,yord=0;
  for (nt=0;nt<p->length;nt++){
    tab = periodize(periodize(p->tab+1+nt,0,p->length),0,p->length);
    L = llistmake();
    llistgrowllitem(L,p->eventra[tab]);
    l=L->first;
    while (l!=NULL){
      r = (struct region *) l->item;
      xord = +1*(double)nt/(double)p->length*side+xoffset;
      yord = -1*(double)r->label/(double)p->nregions*side+yoffset;
      glbox(xord,yord,side/maximum((double)p->nregions,(double)p->length),1,1,1);
      l=l->child;}
    llisttfree(L);L=NULL;}
}

double norm(double x,double y,double z){ return sqrt(x*x+y*y+z*z);}

void crossproduct(double x1,double y1,double z1,double x2,double y2,double z2,double *x3,double *y3,double *z3)
{
  *x3 = y1*z2-y2*z1;
  *y3 = x2*z1-x1*z2;
  *z3 = x1*y2-x2*y1;
}

void rotate100(double *x,double *y,double *z,double c,double s)
{
  /* overwrites *x,*y,*z with rotated coordinates about x-axis with cos(angle)=c and sin(angle)=s */
  double ynew = c**y-s**z;
  double znew = s**y+c**z;
  *y = ynew;
  *z = znew;
}

void rotate010(double *x,double *y,double *z,double c,double s)
{
  /* overwrites *x,*y,*z with rotated coordinates about y-axis with cos(angle)=c and sin(angle)=s */
  double xnew = c**x-s**z;
  double znew = s**x+c**z;
  *x = xnew;
  *z = znew;
}

GLvoid glgrid(int inorout,int xlines,int ylines,double xside,double yside,double xoffset,double yoffset,double rcolor,double gcolor,double bcolor)
{
  /* draws a grid */
  int nl=0;
  double x1=0,x2=0,y1=0,y2=0;
  glColor3f(rcolor,gcolor,bcolor);
  glBegin(GL_LINES);
  if (inorout){
    x1 = 0*xside+xoffset;
    x2 = 1*xside+xoffset;
    for (nl=0;nl<xlines;nl++){
      y1 = (double)(-nl)/(double)xlines*yside+yoffset;
      glVertex3f(x1,y1,0);
      glVertex3f(x2,y1,0);}
    y1 = -0*yside+yoffset;
    y2 = -1*yside+yoffset;
    for (nl=0;nl<ylines;nl++){
      x1 = (double)(nl)/(double)ylines*xside+xoffset;
      glVertex3f(x1,y1,0);
      glVertex3f(x1,y2,0);}}
  else /* if (!inorout) */{
    x1 = (0-0.5/(double)xlines)*xside+xoffset;
    x2 = (1-0.5/(double)xlines)*xside+xoffset;
    for (nl=0;nl<xlines+1;nl++){
      y1 = (double)(-nl+0.5)/(double)xlines*yside+yoffset;
      glVertex3f(x1,y1,0);
      glVertex3f(x2,y1,0);}
    y1 = -(0-0.5/(double)ylines)*yside+yoffset;
    y2 = -(1-0.5/(double)ylines)*yside+yoffset;
    for (nl=0;nl<ylines+1;nl++){
      x1 = (double)(nl-0.5)/(double)ylines*xside+xoffset;
      glVertex3f(x1,y1,0);
      glVertex3f(x1,y2,0);}}
  glEnd();
}

GLvoid Drawllitem_starter_ring(int posorpre,int depth,struct ptree *p,struct pnode *parent,struct llitem *l0,struct llist *L,double wmax,double side,double xoffset,double yoffset)
{
  int firstentry=0;
  struct pnode *pn=NULL;
  double rcolor=0,gcolor=0,bcolor=0;
  double xord=0,yord=0,rord=side/p->nregions;
  struct litem *l=NULL;
  int nl=0,nr=0,nr2=0,nr3=0,nr4=0;
  double epsilon=0,vx=0,vy=0,vrad=0,bigangle=0,smlangle=0;
  if (parent!=NULL && l0->parent==NULL){ /* first descent into parent->childllitem */ 
    assert(parent->childllitem==l0); firstentry=1;}
  if (firstentry){ litemadd(L,parent);}
  if (firstentry){ 
    if (DRAW_FLAG2>0){
      if (p->nregions<=8){
	epsilon=0.5;
	l=L->last;nl=0;xord=0,yord=0,vx=0;vy=0;vrad=PI/(double)p->nregions;vrad=0.5*sin(vrad)/(1+sin(vrad));
	while (l!=NULL){
	  pn=(struct pnode *)l->item;
	  vx = cos(2*PI*(pn->region->label+0.5)/(double)p->nregions)*(0.5-vrad);
	  vy = sin(2*PI*(pn->region->label+0.5)/(double)p->nregions)*(0.5-vrad);
	  xord += pow(3+epsilon,nl)*vx; yord += pow(3+epsilon,nl)*vy;
	  nl+=1;
	  l=l->parent;}
	xord /= pow(3+epsilon,p->nlegs); yord /= pow(3+epsilon,p->nlegs); rord = vrad*pow(3,-p->nlegs);
	pn=(struct pnode *)L->last->item;
	if (wmax>0){ colorscale(7,log(pn->weight),log(wmax),0,&rcolor,&gcolor,&bcolor);}
	else /* if (wmax<=0) */{ colorscale(7,pn->relevance,+1,-1,&rcolor,&gcolor,&bcolor);}
	gldot(p->nregions,xord*side+xoffset+side/2*posorpre,yord*side+yoffset,rord*side,rcolor,gcolor,bcolor);}
      else if (p->nregions>8){
	epsilon=4;
	l=L->last;nl=0;xord=0,yord=0,vx=0;vy=0;
	bigangle=2*PI/(double)((p->nregions-1)/4+1);smlangle=bigangle/2;vrad=0.5*sin(bigangle)/(1+sin(bigangle));rord=vrad/3;
	while (l!=NULL){ 
	  pn=(struct pnode *)l->item;
	  nr = pn->region->label/4; nr2 = pn->region->label%4; nr3 = (nr2>=2?+1:-1); nr4 = (nr2%2?+1:-1);
	  vx = (4.0/6.0)*((1.0/3.0+nr4*vrad/2)*cos(nr*bigangle)+vrad*cos(nr*bigangle+nr3*smlangle));
	  vy = (4.0/6.0)*((1.0/3.0+nr4*vrad/2)*sin(nr*bigangle)+vrad*sin(nr*bigangle+nr3*smlangle));
	  xord += pow(3+epsilon,nl)*vx; yord += pow(3+epsilon,nl)*vy;
	  nl+=1;l=l->parent;}
	xord /= pow(3+epsilon,p->nlegs); yord /= pow(3+epsilon,p->nlegs); rord /= pow(3+epsilon,p->nlegs);
	pn=(struct pnode *)L->last->item;
	if (wmax>0){ colorscale(7,log(pn->weight),log(wmax),0,&rcolor,&gcolor,&bcolor);}
	else /* if (wmax<=0) */{ colorscale(7,pn->relevance,+1,-1,&rcolor,&gcolor,&bcolor);}
	//gldot(p->nregions,xord*side+xoffset+side/2*posorpre,yord*side+yoffset,rord*side,rcolor,gcolor,bcolor);
	glBegin(GL_QUADS);glColor3f(rcolor,gcolor,bcolor);
	nr = pn->region->label/4; nr2 = pn->region->label%4; nr3 = (nr2>=2?+1:-1); nr4 = (nr2%2?+1:-1);
	glVertex3f((xord+rord*cos(nr*bigangle+nr3*smlangle))*side+xoffset+side/2*posorpre,(yord+rord*sin(nr*bigangle+nr3*smlangle))*side+yoffset,0);
	glVertex3f((xord+rord/2*cos(nr*bigangle+nr3*smlangle+PI/2))*side+xoffset+side/2*posorpre,(yord+rord/2*sin(nr*bigangle+nr3*smlangle+PI/2))*side+yoffset,0);
	glVertex3f((xord+rord*cos(nr*bigangle+nr3*smlangle+PI))*side+xoffset+side/2*posorpre,(yord+rord*sin(nr*bigangle+nr3*smlangle+PI))*side+yoffset,0);
	glVertex3f((xord+rord/2*cos(nr*bigangle+nr3*smlangle-PI/2))*side+xoffset+side/2*posorpre,(yord+rord/2*sin(nr*bigangle+nr3*smlangle-PI/2))*side+yoffset,0);
	glEnd();
}}
    else /* if DRAW_FLAG2<=0 */{
      epsilon=0.1;
      l=L->last;nl=0;xord=0;yord=0;vrad=PI/(double)p->nregions;vrad=0.5*sin(vrad)/(1+sin(vrad));
      while (l!=NULL){
	pn=(struct pnode *)l->item;
	switch(pn->region->label){ case 0:case 1:case 2: vy=1.0/3.0;break;case 5:case 6:case 7:vy=-1.0/3.0;break;default:vy=0;break;}
	switch(pn->region->label){ case 0:case 3:case 5: vx=-1.0/3.0;break;case 2:case 4:case 7:vx=1.0/3.0;break;default:vx=0;break;}
	xord += pow(3+epsilon,nl)*vx; yord += pow(3+epsilon,nl)*vy;
	nl+=1;
	l=l->parent;}
      xord /= pow(3+epsilon,p->nlegs); yord /= pow(3+epsilon,p->nlegs); rord = vrad*pow(3,-p->nlegs);
      pn=(struct pnode *)L->last->item;
      if (wmax>0){ colorscale(7,log(pn->weight),log(wmax),0,&rcolor,&gcolor,&bcolor);}
      else /* if (wmax<=0) */{ colorscale(7,pn->relevance,+1,-1,&rcolor,&gcolor,&bcolor);}
      glbox(xord*side+xoffset+side/2*posorpre,yord*side+yoffset,rord*side,rcolor,gcolor,bcolor);}}
  if (l0->kidl!=NULL){ Drawllitem_starter_ring(posorpre,depth,p,parent,l0->kidl,L,wmax,side,xoffset,yoffset);}
  if (l0->item!=NULL /* && (L->length-depth)<3 */){ pn=(struct pnode *)l0->item; Drawllitem_starter_ring(posorpre,depth,p,pn,pn->childllitem,L,wmax,side,xoffset,yoffset);}
  if (l0->kidr!=NULL){ Drawllitem_starter_ring(posorpre,depth,p,parent,l0->kidr,L,wmax,side,xoffset,yoffset);}
  if (firstentry){ llistkilllast(L);}
}

GLvoid Drawllitem_starter(int posorpre,int depth,struct ptree *p,struct pnode *parent,struct llitem *l0,struct llist *L,double side,double xoffset,double yoffset)
{
  int firstentry=0;
  struct pnode *pn=NULL,*pn1=NULL,*pn2=NULL,*pn3=NULL;
  double rcolor=0,gcolor=0,bcolor=0;
  double xord=0,yord=0,zord=0,rord=side/p->nregions;
  if (parent!=NULL && l0->parent==NULL){ /* first descent into parent->childllitem */ 
    assert(parent->childllitem==l0); firstentry=1;}
  if (firstentry){ litemadd(L,parent);}
  if (firstentry){ 
    switch(L->length-depth){
    case 0: break;
    case 1: 
      pn1 = (struct pnode *) L->last->item;
      xord = (0.5+0.1*posorpre)*side + xoffset;
      yord = -1*(double)pn1->region->label*rord + yoffset;
      colorscale(7,pn1->relevance,+1,-1,&rcolor,&gcolor,&bcolor);
      //glbox(xord,yord,rord,rcolor,gcolor,bcolor);
      glnum(xord,yord,0,pn1->weight,rcolor,gcolor,bcolor);
      //glnum(xord+0.1,yord+0.1,0,pn1->relevance,1,1,1);
      break;
    case 2: 
      pn1 = (struct pnode *) L->last->parent->item;
      pn2 = (struct pnode *) L->last->item;
      xord = +1*(double)pn1->region->label*rord + posorpre*1.0*side + xoffset;
      yord = -1*(double)pn2->region->label*rord + yoffset;
      colorscale(7,pn2->relevance,+1,-1,&rcolor,&gcolor,&bcolor);
      glbox(xord,yord,rord,rcolor,gcolor,bcolor);
      glnum(xord-0.1,yord,0.01,pn2->weight,0,0,0);
      //glnum(xord+0.1,yord+0.1,0,pn2->relevance,1,1,1);
      break;
    case 3: 
      pn1 = (struct pnode *) L->last->parent->parent->item;
      pn2 = (struct pnode *) L->last->parent->item;
      pn3 = (struct pnode *) L->last->item;

/*       xord = +1*(double)(pn1->region->label + 0.5)/(double)p->nregions - 0.5; */
/*       yord = -1*(double)(pn2->region->label + 0.5)/(double)p->nregions + 0.5; */
/*       zord = -1*(double)(pn3->region->label + 0.5)/(double)p->nregions + 0.5; */
/*       rotate100(&xord,&yord,&zord,cos(RXangle),sin(RXangle)); */
/*       rotate010(&xord,&yord,&zord,cos(RYangle),sin(RYangle)); */
/*       colorscale(7,pn3->relevance,+1,-1,&rcolor,&gcolor,&bcolor); */
/*       glBegin(GL_POINTS); */
/*       glColor3f(rcolor,gcolor,bcolor); */
/*       glVertex3f((xord+0.5)*side + posorpre*2.0*side+xoffset,(yord-0.5)*side + yoffset,(zord-0.5)*side + 0); */
/*       glEnd(); */

      xord = +1*(double)pn1->region->label*rord + posorpre*2.0*side + xoffset;
      yord = -1*(double)pn2->region->label*rord + yoffset;
      zord = -1*(double)pn3->region->label*rord + 0;
      colorscale(7,pn3->relevance,+1,-1,&rcolor,&gcolor,&bcolor);
      rord = (double)(1 + pn3->region->label)/(double)p->nregions*side/(double)p->nregions;
      zord = -0.001*(double)pn3->region->label/(double)p->nregions;
      glBegin(GL_QUADS);
      glColor3f(rcolor,gcolor,bcolor);
      glVertex3f(xord+rord/2,yord+rord/2,zord);
      glVertex3f(xord+rord/2,yord-rord/2,zord);
      glVertex3f(xord-rord/2,yord-rord/2,zord);
      glVertex3f(xord-rord/2,yord+rord/2,zord);
      glEnd();

      //glcube(xord,yord,zord,rord/2,rcolor,gcolor,bcolor);
      //glnum(xord,yord,zord,pn3->weight,rcolor,gcolor,bcolor);
      break;
    default: break;}}
  if (l0->kidl!=NULL){ Drawllitem_starter(posorpre,depth,p,parent,l0->kidl,L,side,xoffset,yoffset);}
  if (l0->item!=NULL && (L->length-depth)<3){ pn=(struct pnode *)l0->item; Drawllitem_starter(posorpre,depth,p,pn,pn->childllitem,L,side,xoffset,yoffset);}
  if (l0->kidr!=NULL){ Drawllitem_starter(posorpre,depth,p,parent,l0->kidr,L,side,xoffset,yoffset);}
  if (firstentry){ llistkilllast(L);}
}

GLvoid Drawptree(struct ptree *p,double side,double xoffset,double yoffset)
{
  struct neuronarray *Nra=GLOBAL_Nra;
  struct region *r=NULL,*r2=NULL;
  struct llitem *l0=NULL;
  struct pnode *pn=NULL;
  char text[128];
  int nr=0,nc=0,nl=0;
  double bigside=maximum(Nra->rows,Nra->cols);
  double rcolor=0,gcolor=0,bcolor=0;
  int fr=0,frmin=0,frmax=0;
  int fc=0,fcmin=0,fcmax=0;
  int rsides=7;
  double cA=0,cN=0,cG=0,mA=0,mN=0,mG=0;
  struct neuron *s=NULL,*n=NULL;
  int bign = NPIEROWS*NPIECOLS*PIE_ROW_DIA*PIE_COL_DIA;
  int tab=0;
  double max=0,temp=0;
  struct llist *L=NULL;
  struct litem *l=NULL;
  double *connections=NULL;
  double xord=0,yord=0;
  double wmax=0;
  if (p!=NULL){
    switch(abs(DRAW_FLAG)%5){
    case 4: /* relevance histogram */
      Plotra(p->wh->data,"double",p->wh->nbins,pow(p->nregions,p->nlegs+1)/p->wh->nbins,0,side,side,xoffset,yoffset,0,1,1);
      Plotra(p->rh->data,"double",p->rh->nbins,pow(p->nregions,p->nlegs+1)/p->rh->nbins,0,side,side,xoffset,yoffset,1,1,0);
      glBegin(GL_LINES);
      glColor3f(0.6,0.6,0.6);
      xord = 0.0; yord = 0.0; glVertex3f(xord*side+xoffset+0*side,yord*side+yoffset+0*side,0);
      xord = 1.0; yord = 0.0; glVertex3f(xord*side+xoffset+0*side,yord*side+yoffset+0*side,0);
      xord = 0.5; yord = 0.0; glVertex3f(xord*side+xoffset+0*side,yord*side+yoffset+0*side,0);
      xord = 0.5; yord = 1.0; glVertex3f(xord*side+xoffset+0*side,yord*side+yoffset+0*side,0);
      xord = 0.0; yord = 1.0; glVertex3f(xord*side+xoffset+0*side,yord*side+yoffset+0*side,0);
      xord = 1.0; yord = 1.0; glVertex3f(xord*side+xoffset+0*side,yord*side+yoffset+0*side,0);
      xord = 0.0; yord = 0.0; glVertex3f(xord*side+xoffset+0*side,yord*side+yoffset+0*side,0);
      xord = 0.0; yord = 1.0; glVertex3f(xord*side+xoffset+0*side,yord*side+yoffset+0*side,0);
      xord = 1.0; yord = 0.0; glVertex3f(xord*side+xoffset+0*side,yord*side+yoffset+0*side,0);
      xord = 1.0; yord = 1.0; glVertex3f(xord*side+xoffset+0*side,yord*side+yoffset+0*side,0);
      glEnd();
      break;
    case 3: /* region->neuronllist */
      for (nr=0;nr<p->nregions;nr++){
	colorscale(7,nr,p->nregions-1,0,&rcolor,&gcolor,&bcolor);
	l = p->regionra[nr]->neuronllist->first;
	while (l!=NULL){
	  n=(struct neuron *)l->item;
	  glbox(n->col/bigside*side+xoffset,-n->row/bigside*side+yoffset,side/bigside,rcolor,gcolor,bcolor);
	  l=l->child;}}
      break;
    case 2: /* raster */
      Draweventra(p,side,xoffset,yoffset);
      break;
    case 1: /* ptree */
      glgrid(0,p->nregions,p->nregions,side,side,xoffset+side,yoffset,1,1,1);
      frmin = -1; frmax = p->nregions;
      fcmin = -1; fcmax = p->nregions;
      periodify("int",&FIDDLE_ROW,&frmin,&frmax,&fr);
      periodify("int",&FIDDLE_COL,&fcmin,&fcmax,&fc);
      if (fr==-1){
	L=llistmake(); if (DRAW_FLAG2%2){ pnodestats_starter(NULL,p->postree,-1,-1,0,NULL,&wmax,NULL,NULL,NULL,NULL,NULL,NULL,NULL); Drawllitem_starter_ring(+1,0,p,NULL,p->postree,L,wmax,side,xoffset,yoffset);} else{ Drawllitem_starter(+1,0,p,NULL,p->postree,L,side,xoffset,yoffset);} llisttfree(L); L=NULL;
	L=llistmake(); if (DRAW_FLAG2%2){ pnodestats_starter(NULL,p->pretree,-1,-1,0,NULL,&wmax,NULL,NULL,NULL,NULL,NULL,NULL,NULL); Drawllitem_starter_ring(-1,0,p,NULL,p->pretree,L,wmax,side,xoffset,yoffset);} else{ Drawllitem_starter(-1,0,p,NULL,p->pretree,L,side,xoffset,yoffset);} llisttfree(L); L=NULL;}
      else /* if (fr>=0) */{
	if (fc==-1){
	  r = p->regionra[fr];
	  if ((l0 = llitemaddorfind(0,p->postree,r,&region2pnode_compare_label))!=NULL){
	    pn = (struct pnode *) l0->item;
	    L=llistmake();
	    Drawllitem_starter(+1,1,p,pn,pn->childllitem,L,side,xoffset,yoffset);
	    llisttfree(L);L=NULL;}
	  if ((l0 = llitemaddorfind(0,p->pretree,r,&region2pnode_compare_label))!=NULL){
	    pn = (struct pnode *) l0->item;
	    L=llistmake();
	    Drawllitem_starter(-1,1,p,pn,pn->childllitem,L,side,xoffset,yoffset);
	    llisttfree(L);L=NULL;}}
	else /* if (fc>=0) */{
	  r = p->regionra[fr];
	  r2 = p->regionra[fc];
	  if ((l0 = llitemaddorfind(0,p->postree,r,&region2pnode_compare_label))!=NULL){
	    pn = (struct pnode *) l0->item;
	    L=llistmake(); litemadd(L,pn);
	    if ((l0 = llitemaddorfind(0,pn->childllitem,r2,&region2pnode_compare_label))!=NULL){
	      pn = (struct pnode *) l0->item;
	      Drawllitem_starter(+1,2,p,pn,pn->childllitem,L,side,xoffset,yoffset);}
	    llisttfree(L);L=NULL;}
	  if ((l0 = llitemaddorfind(0,p->pretree,r,&region2pnode_compare_label))!=NULL){
	    pn = (struct pnode *) l0->item;
	    L=llistmake(); litemadd(L,pn);
	    if ((l0 = llitemaddorfind(0,pn->childllitem,r2,&region2pnode_compare_label))!=NULL){
	      pn = (struct pnode *) l0->item;
	      Drawllitem_starter(-1,2,p,pn,pn->childllitem,L,side,xoffset,yoffset);}
	    llisttfree(L);L=NULL;}}}
      break;
    case 0: /* cortical connections */
      switch (DRAW_FLAG2%2){
      case 0: /* network representation */
	frmin = 0; frmax = Nra->rows;
	fcmin = 0; fcmax = Nra->cols;
	periodify("int",&FIDDLE_ROW,&frmin,&frmax,&fr);
	periodify("int",&FIDDLE_COL,&fcmin,&fcmax,&fc);
	for (nr=0;nr<Nra->rows;nr++){ for (nc=0;nc<Nra->cols;nc++){
	  colorscale(3,Nra->t2sra[nr+nc*Nra->rows],1,0,&rcolor,&gcolor,&bcolor);
	  gldot(rsides,nc/bigside*side+xoffset,-nr/bigside*side+yoffset,1*side/bigside/8,rcolor,gcolor,bcolor);}}
	s=nget(Nra,fr,fc);
	glring(rsides,s->col/bigside*side+xoffset,-s->row/bigside*side+yoffset,3*side/bigside/8,4*side/bigside/8,1,1,1);
	for (nr=0;nr<Nra->rows;nr++){ for (nc=0;nc<Nra->cols;nc++){
	  n=nget(Nra,nr,nc);
	  if (slink(s,n,&cA,&cN,&cG,&mA,&mN,&mG)){
	    colorscale(0,maximum(cA,cG),maximum(mA,mG),0,&rcolor,&gcolor,&bcolor);
	    glring(rsides,n->col/bigside*side+xoffset,-n->row/bigside*side+yoffset,1*side/bigside/8,2.5*side/bigside/8,(*(s->t2s)>=2)*rcolor,0,(*(s->t2s)<=1)*bcolor);}}}
	n=s; nl=n->lr_kernel_type;
	switch ((int)*(n->t2s)){
	case 0: max=CS_ICLR; break;
	case 1: max=CS_ISLR; break;
	case 2: max=CS_ECLR; break;
	case 3: max=CS_ESLR; break;
	default: printf("error! wrong (int)*(n->t2s)=%d in Drawptree\n",(int)*(n->t2s)); break;}
	for (nr=0;nr<bign;nr++){ GLOBAL_LRSWAP_X[nl][nr][0]=0; GLOBAL_LRSWAP_X[nl][nr][1]=0;}
	GLOBAL_LRSWAP_X[nl][n->pierow+n->piecol*NPIEROWS+n->ang*NPIEROWS*NPIECOLS][0] += max;
	fftw_execute_dft(GLOBAL_FFTWPLAN_SWAP_FORWARD[nl],GLOBAL_LRSWAP_X[nl],GLOBAL_LRSWAP_K[nl]);
	for (nr=0;nr<bign;nr++){ GLOBAL_LRSWAP_X[nl][nr][0] = GLOBAL_LRKERNEL_K[nl][nr][0]*GLOBAL_LRSWAP_K[nl][nr][0] - GLOBAL_LRKERNEL_K[nl][nr][1]*GLOBAL_LRSWAP_K[nl][nr][1]; GLOBAL_LRSWAP_X[nl][nr][1] = GLOBAL_LRKERNEL_K[nl][nr][0]*GLOBAL_LRSWAP_K[nl][nr][1] + GLOBAL_LRKERNEL_K[nl][nr][1]*GLOBAL_LRSWAP_K[nl][nr][0];}
	fftw_execute_dft(GLOBAL_FFTWPLAN_SWAP_BACKWARD[nl],GLOBAL_LRSWAP_X[nl],GLOBAL_LRSWAP_K[nl]);
	for (nr=0;nr<Nra->rows;nr++){ for (nc=0;nc<Nra->cols;nc++){
	  n = nget(Nra,nr,nc);
	  tab = n->ang + n->pierow*PIE_ROW_DIA*PIE_COL_DIA + n->piecol*PIE_ROW_DIA*PIE_COL_DIA*NPIEROWS;
	  temp = GLOBAL_LRSWAP_K[nl][n->pierow+n->piecol*NPIEROWS+GLOBAL_SHUFFLES[tab]*NPIEROWS*NPIECOLS][0]/(double)bign;
	  switch((int)*(n->t2s)){
	  case 0: /* IC */ temp *= CS_LRIC; break;
	  case 1: /* IS */ temp *= CS_LRIS; break;
	  case 2: /* EC */ temp *= CS_LREC; break;
	  case 3: /* ES */ temp *= CS_LRES; break;
	  default: printf("error! wrong (int)*(n->t2s)=%d in Drawptree\n",(int)*(n->t2s)); break;}
	  if (temp!=0){
	    colorscale(0,temp,max,0,&rcolor,&gcolor,&bcolor);
	    glring(rsides,n->col/bigside*side+xoffset,-n->row/bigside*side+yoffset,2.5*side/bigside/8,4*side/bigside/8,0,gcolor,0);}}}
	glColor3f(1,1,1);
	sprintf(text,"(%d,%d) connections",s->row,s->col);
	ftexto(xoffset,yoffset+0.1,0,text);
	break;
      case 1: /* matrix representation */
	connections = (double *) tcalloc((Nra->rows*Nra->cols)*(Nra->rows*Nra->cols+2),sizeof(double));
	for (nr=0;nr<Nra->rows*Nra->cols;nr++){ for (nc=0;nc<Nra->rows*Nra->cols;nc++){
	  s = Nra->N[nr]; n = Nra->N[nc];
	  if (slink(s,n,&cA,&cN,&cG,&mA,&mN,&mG)){ 
	    if (*(s->t2s)<=1){ connections[nr+nc*Nra->rows*Nra->cols] = -cG;}
	    else if (*(s->t2s)>1){ connections[nr+nc*Nra->rows*Nra->cols] = cA+cN;}}}}
	stats("double",connections,(Nra->rows*Nra->cols)*(Nra->rows*Nra->cols),&mA,&mG,NULL,NULL); mN=maximum(fabs(mA),fabs(mG));
	for (nr=0;nr<Nra->rows*Nra->cols;nr++){
	  connections[nr+(Nra->rows*Nra->cols+1)*(Nra->rows*Nra->cols)] = -mN + 2*mN*(double)nr/(double)(Nra->rows*Nra->cols-1);}
	Drawra(connections,"double",Nra->rows*Nra->cols,Nra->rows*Nra->cols+2,7,0,0,"connections",side,side,xoffset,yoffset);
	tfree(connections);connections=NULL;
	break;
      default: break;}
      break;
    default: break;}}
}

void Drawlyapunov(double side,double xoffset,double yoffset)
{
  struct lyapunov *y=GLOBAL_LYAPUNOV;
  int nbins=32,nb=0;
  double mean=0,stdev=0,max=0,min=0,*ra=NULL;
  struct litem *l=NULL;
  char text[32];
  lliststats(y->L,NULL,NULL,&mean,&stdev);
  max=mean+STD_VIEW*stdev;
  min=mean-STD_VIEW*stdev;
  ra = (double *) tcalloc(nbins,sizeof(double));
  l=y->L->first;
  while (l!=NULL){
    nb=(int)floor(nbins*(*(double *)l->item-min)/(max-min));
    if (nb>=0 && nb<nbins){ ra[nb]+=1;}
    l=l->child;}
  Plotra(ra,"double",nbins,0,0,side,side,xoffset,yoffset,0.8,0.8,0.8);
  tfree(ra);
  glColor3f(1,1,1);
  sprintf(text,"L[%0.2f,%0.2f]",min,max);
  ftexto(xoffset,yoffset+0.1,0,text);
  lliststats(y->pL,NULL,NULL,&mean,&stdev);
  max=mean+STD_VIEW*stdev;
  min=mean-STD_VIEW*stdev;
  ra = (double *) tcalloc(nbins,sizeof(double));
  l=y->pL->first;
  while (l!=NULL){
    nb=(int)floor(nbins*(*(double *)l->item-min)/(max-min));
    if (nb>=0 && nb<nbins){ ra[nb]+=1;}
    l=l->child;}
  Plotra(ra,"double",nbins,0,0,side,side,xoffset,yoffset,1,0.5,1);
  tfree(ra);
  glColor3f(1,1,1);
  sprintf(text,"pL[%0.2f,%0.2f]",min,max);
  ftexto(xoffset,yoffset+0.2,0,text);
}

void Drawpower(double side,double xoffset,double yoffset)
{
  struct power *p=GLOBAL_POWER;
  char text[32];
  int drawvsplot = DRAW_FLAG2%2;
  switch (DRAW_FLAG%8){
  case 0: 
    sprintf(text,"V"); Drawstra(p->Vstra,p->N,side,xoffset,yoffset,VOLTAGE_THRESHOLD,VOLTAGE_IN,-1,2,0,drawvsplot); 
    if (drawvsplot){ loglogra(p->Vpower,"double",p->length,0,0,side,side,xoffset,yoffset+side,1,1,1);}
    else{ Plotra(p->Vpower,"double",p->length,0,0,side,side,xoffset,yoffset+side,1,1,1);} break;
  case 1: 
    sprintf(text,"sA");  Drawstra(p->sAstra,p->N,side,xoffset,yoffset,0,0,-1,2,0,drawvsplot);
    if (drawvsplot){ loglogra(p->sApower,"double",p->length,0,0,side,side,xoffset,yoffset+side,1,1,1);}
    else{ Plotra(p->sApower,"double",p->length,0,0,side,side,xoffset,yoffset+side,1,1,1);} break;
  case 2: 
    sprintf(text,"sN");  Drawstra(p->sNstra,p->N,side,xoffset,yoffset,0,0,-1,2,0,drawvsplot);
    if (drawvsplot){ loglogra(p->sNpower,"double",p->length,0,0,side,side,xoffset,yoffset+side,1,1,1);}
    else{ Plotra(p->sNpower,"double",p->length,0,0,side,side,xoffset,yoffset+side,1,1,1);} break;
  case 3: 
    sprintf(text,"sG");  Drawstra(p->sGstra,p->N,side,xoffset,yoffset,0,0,-1,2,0,drawvsplot);
    if (drawvsplot){ loglogra(p->sGpower,"double",p->length,0,0,side,side,xoffset,yoffset+side,1,1,1);}
    else{ Plotra(p->sGpower,"double",p->length,0,0,side,side,xoffset,yoffset+side,1,1,1);} break;
  case 4: 
    sprintf(text,"VS");  Drawstra(p->VSstra,p->N,side,xoffset,yoffset,VOLTAGE_THRESHOLD,VOLTAGE_IN,-1,2,0,drawvsplot);
    if (drawvsplot){ loglogra(p->VSpower,"double",p->length,0,0,side,side,xoffset,yoffset+side,1,1,1);}
    else{ Plotra(p->VSpower,"double",p->length,0,0,side,side,xoffset,yoffset+side,1,1,1);} break;
  case 5: sprintf(text,"ac"); 
    if (drawvsplot){ loglogra(p->autocorrelation,"double",p->length,0,0,side,side,xoffset,yoffset,1,1,1);}
    else{ Plotra(p->autocorrelation,"double",p->length,0,0,side,side,xoffset,yoffset,1,1,1);} break;
  case 6: sprintf(text,"xc"); 
    if (drawvsplot){ loglogra(p->crosscorrelation,"double",p->length,0,0,side,side,xoffset,yoffset,1,1,1);}
    else{ Plotra(p->crosscorrelation,"double",p->length,0,0,side,side,xoffset,yoffset,1,1,1);} break;
  case 7: sprintf(text,"raster"); Drawstra(p->raster,p->N,side,xoffset,yoffset,1,0,-1,2,0,drawvsplot); break;
  default: break;}
  glColor3f(1,1,1);
  ftexto(xoffset,yoffset+0.1,0,text);
}

void Drawtaof(double side,double xoffset,double yoffset)
{
  struct taof *tf=GLOBAL_TAOF;
  char text[32];
  int drawvsplot = DRAW_FLAG2%2;
  switch (DRAW_FLAG%5){
  case 0: sprintf(text,"input_contrast"); Drawstra(&(tf->input_contrast),1,side,xoffset,yoffset,0,0,-1,2,0,drawvsplot); break;
  case 1: sprintf(text,"input_spaceangle"); Drawstra(&(tf->input_spaceangle),1,side,xoffset,yoffset,0,0,-1,2,0,drawvsplot); break;
  case 2: sprintf(text,"f1r"); Drawstra(tf->f1r,tf->Nra->rows*tf->Nra->cols,side,xoffset,yoffset,0,0,-1,2,0,drawvsplot); break;
  case 3: sprintf(text,"f1i"); Drawstra(tf->f1i,tf->Nra->rows*tf->Nra->cols,side,xoffset,yoffset,0,0,-1,2,0,drawvsplot); break;
  case 4: sprintf(text,"f0"); Drawstra(tf->f0,tf->Nra->rows*tf->Nra->cols,side,xoffset,yoffset,0,0,-1,2,0,drawvsplot); break;
  default: break;}
  glColor3f(1,1,1);
  ftexto(xoffset,yoffset-0.1,0,text);
}

GLvoid glbox(double x,double y,double r,double rcolor,double gcolor,double bcolor)
{
  glBegin(GL_QUADS);
  glColor3f(rcolor,gcolor,bcolor);
  glVertex3f(x+r/2,y+r/2,0);
  glVertex3f(x+r/2,y-r/2,0);
  glVertex3f(x-r/2,y-r/2,0);
  glVertex3f(x-r/2,y+r/2,0);
  glEnd();
}

GLvoid glnum(double x,double y,double z,double weight,double rcolor,double gcolor,double bcolor)
{
  char text[16];
  glColor3f(rcolor,gcolor,bcolor);
  sprintf(text,"%0.0f",weight);
  ftexto(x,y,z,text);
}

GLvoid gldot(int sides,double x,double y,double r,double rcolor,double gcolor,double bcolor)
{
  int ns=0;
  glBegin(GL_POLYGON);
  glColor3f(rcolor,gcolor,bcolor);
  for (ns=0;ns<sides;ns++){
    glVertex3f(x+r*cos(2*PI*((double)ns+0.5)/(double)sides),y+r*sin(2*PI*((double)ns+0.5)/(double)sides),0);}
  glEnd();
}

GLvoid glring(int sides,double x,double y,double rin,double rout,double rcolor,double gcolor,double bcolor)
{
  int ns=0;
  glBegin(GL_QUADS);
  glColor3f(rcolor,gcolor,bcolor);
  for (ns=0;ns<sides;ns++){
    glVertex3f(x+rin*cos(2*PI*((double)ns+0.5)/(double)sides),y+rin*sin(2*PI*((double)ns+0.5)/(double)sides),0);
    glVertex3f(x+rout*cos(2*PI*((double)ns+0.5)/(double)sides),y+rout*sin(2*PI*((double)ns+0.5)/(double)sides),0);
    glVertex3f(x+rout*cos(2*PI*((double)ns-0.5)/(double)sides),y+rout*sin(2*PI*((double)ns-0.5)/(double)sides),0);
    glVertex3f(x+rin*cos(2*PI*((double)ns-0.5)/(double)sides),y+rin*sin(2*PI*((double)ns-0.5)/(double)sides),0);}
  glEnd();
}

GLvoid glarrow(int nheads,double x1,double y1,double x2,double y2,double rmin,double rcolor,double gcolor,double bcolor)
{
  double theta = atan2(y2-y1,x2-x1),theta1=theta+PI-PI/8,theta2=theta+PI+PI/8;
  int nh=0;
  double x3=(x1+x2)/2.0,y3=(y1+y2)/2.0;
  double x4=0,y4=0;
  glBegin(GL_LINES);
  glColor3f(rcolor,gcolor,bcolor);
  glVertex3f(x1,y1,0); glVertex3f(x2,y2,0);
  for (nh=0;nh<nheads;nh++){
    x4 = x3 + (double)nh*rmin/5.0*cos(theta);
    y4 = y3 + (double)nh*rmin/5.0*sin(theta);
    glVertex3f(x4,y4,0);
    glVertex3f(x4+rmin*cos(theta1),y4+rmin*sin(theta1),0);
    glVertex3f(x4,y4,0);
    glVertex3f(x4+rmin*cos(theta2),y4+rmin*sin(theta2),0);}
  glEnd();
}

GLvoid glcube(double x,double y,double z,double r,double rcolor,double gcolor,double bcolor)
{
  glBegin(GL_QUADS);
  glColor3f(rcolor,gcolor,bcolor);
  /* top side */
  glVertex3f(x+r/2,y+r/2,z+r/2);
  glVertex3f(x+r/2,y-r/2,z+r/2);
  glVertex3f(x-r/2,y-r/2,z+r/2);
  glVertex3f(x-r/2,y+r/2,z+r/2);
  /* bottom side */
  glVertex3f(x+r/2,y+r/2,z-r/2);
  glVertex3f(x+r/2,y-r/2,z-r/2);
  glVertex3f(x-r/2,y-r/2,z-r/2);
  glVertex3f(x-r/2,y+r/2,z-r/2);
  /* right side */
  glVertex3f(x+r/2,y+r/2,z+r/2);
  glVertex3f(x+r/2,y+r/2,z-r/2);
  glVertex3f(x+r/2,y-r/2,z-r/2);
  glVertex3f(x+r/2,y-r/2,z+r/2);
  /* left side */
  glVertex3f(x-r/2,y+r/2,z+r/2);
  glVertex3f(x-r/2,y+r/2,z-r/2);
  glVertex3f(x-r/2,y-r/2,z-r/2);
  glVertex3f(x-r/2,y-r/2,z+r/2);
  /* back side */
  glVertex3f(x+r/2,y+r/2,z+r/2);
  glVertex3f(x+r/2,y+r/2,z-r/2);
  glVertex3f(x-r/2,y+r/2,z-r/2);
  glVertex3f(x-r/2,y+r/2,z+r/2);
  /* front side */
  glVertex3f(x+r/2,y-r/2,z+r/2);
  glVertex3f(x+r/2,y-r/2,z-r/2);
  glVertex3f(x-r/2,y-r/2,z-r/2);
  glVertex3f(x-r/2,y-r/2,z+r/2);
  glEnd();
}

GLvoid Drawcolorbar(int coloringtype,double side,double xoffset,double yoffset)
{
  int nv=500,nc=0;
  double xord=0,yord=0;
  double val=0,dval=0.5/(double)nv;
  double rcolor=0,gcolor=0,bcolor=0;
  for (nc=0;nc<nv;nc++){
    val = (double)(nc+0.5)/(double)nv;
    colorscale(coloringtype,val,1,0,&rcolor,&gcolor,&bcolor);
    glBegin(GL_QUADS);
    glColor3f(rcolor,gcolor,bcolor);
    xord = val-dval; yord = 0; glVertex3f(xord*side+xoffset,yord*side+yoffset,0);
    xord = val+dval; yord = 0; glVertex3f(xord*side+xoffset,yord*side+yoffset,0);
    xord = val+dval; yord = 1; glVertex3f(xord*side+xoffset,yord*side+yoffset,0);
    xord = val-dval; yord = 1; glVertex3f(xord*side+xoffset,yord*side+yoffset,0);
    glEnd();}
}

GLvoid Drawmenu(double side,double xoffset,double yoffset)
{
  char text[64];
  int menupos=0;
  double tcolor=0;
  struct neuronarray *Nra=GLOBAL_Nra;
  menupos=-6;
  glColor3f(1,1,1);
  sprintf(text,"%d,%d neurons, %d,%d pies, lr_range=%0.1f",SYSTEM_ROW_SIZE,SYSTEM_COL_SIZE,NPIEROWS,NPIECOLS,LR_DIST_DIA_INPIES*minimum(PIE_ROW_DIA,PIE_COL_DIA));
  ftexto(0*side+xoffset,menupos*side+yoffset,0,text);
  menupos=-5;
  glColor3f(1,1,1);
  sprintf(text,"t=%0.1f,DT=%0.3f",GLOBAL_time,GLOBAL_DT);
  ftexto(0*side+xoffset,menupos*side+yoffset,0,text);
  if (CORTEX_BOTHER){ 
    menupos=-4;
    tcolor = .25+.75*(FIDDLE_PARAMETER==menupos);glColor3f(tcolor,tcolor,tcolor);
    sprintf(text,"reset? ES%0.1f,IS%0.1f,EC%0.1f,IC%0.1f",Nra->mES,Nra->mIS,Nra->mEC,Nra->mIC);
    ftexto(0*side+xoffset,menupos*side+yoffset,0,text);}
  menupos=-3;
  tcolor = .25+.75*(FIDDLE_PARAMETER==menupos);glColor3f(tcolor,tcolor,tcolor);
  sprintf(text,"GLOBAL_ODOR=%d,BASE=%d",GLOBAL_ODOR,GLOBAL_ODOR_BASE);
  ftexto(0*side+xoffset,menupos*side+yoffset,0,text);
  menupos=-2;
  tcolor = .25+.75*(FIDDLE_PARAMETER==menupos);glColor3f(tcolor,tcolor,tcolor);
  sprintf(text,"LGN_STRENGTH=%0.1f=e^%0.1f,OTHERLAYER_STRENGTH=%0.1f=e^%0.1f",LGN_STRENGTH,log(LGN_STRENGTH),OTHERLAYER_STRENGTH,log(OTHERLAYER_STRENGTH));
  ftexto(0*side+xoffset,menupos*side+yoffset,0,text);
  menupos=-1;
  tcolor = .25+.75*(FIDDLE_PARAMETER==menupos);glColor3f(tcolor,tcolor,tcolor);
  sprintf(text,"LGN_BACKRATE=%0.1f=e^%0.1f,OTHERLAYER_BACKRATE=%0.1f=e^%0.1f",LGN_BACKRATE,log(LGN_BACKRATE),OTHERLAYER_BACKRATE,log(OTHERLAYER_BACKRATE));
  ftexto(0*side+xoffset,menupos*side+yoffset,0,text);
  menupos=0;
  tcolor = .25+.75*(FIDDLE_PARAMETER==menupos);glColor3f(tcolor,tcolor,tcolor);
  sprintf(text,"CS_ESLR=%0.1f=e^%0.1f,CS_ECLR=%0.1f=e^%0.1f",CS_ESLR,log(CS_ESLR),CS_ECLR,log(CS_ECLR));
  ftexto(0*side+xoffset,menupos*side+yoffset,0,text);
  menupos=1;
  tcolor = .25+.75*(FIDDLE_PARAMETER==menupos);glColor3f(tcolor,tcolor,tcolor);
  sprintf(text,"CS_ISLR=%0.1f=e^%0.1f,CS_ICLR=%0.1f=e^%0.1f",CS_ISLR,log(CS_ISLR),CS_ICLR,log(CS_ICLR));
  ftexto(0*side+xoffset,menupos*side+yoffset,0,text);
  menupos=2;
  tcolor = .25+.75*(FIDDLE_PARAMETER==menupos);glColor3f(tcolor,tcolor,tcolor);
  sprintf(text,"CS_LRES=%0.1f=e^%0.1f,CS_LREC=%0.1f=e^%0.1f",CS_LRES,log(CS_LRES),CS_LREC,log(CS_LREC));
  ftexto(0*side+xoffset,menupos*side+yoffset,0,text);
  menupos=3;
  tcolor = .25+.75*(FIDDLE_PARAMETER==menupos);glColor3f(tcolor,tcolor,tcolor);
  sprintf(text,"CS_LRIS=%0.1f=e^%0.1f,CS_LRIC=%0.1f=e^%0.1f",CS_LRIS,log(CS_LRIS),CS_LRIC,log(CS_LRIC));
  ftexto(0*side+xoffset,menupos*side+yoffset,0,text);
  tcolor = .25+.75*(FIDDLE_PARAMETER==4);glColor3f(tcolor,tcolor,tcolor);
  sprintf(text,"CS_ESES_A=%0.1f, CS_ECES_A=%0.1f",CS_ESES_A,CS_ECES_A);
  ftexto(0*side+xoffset,4*side+yoffset,0,text);
  tcolor = .25+.75*(FIDDLE_PARAMETER==5);glColor3f(tcolor,tcolor,tcolor);
  sprintf(text,"CS_ESIS_A=%0.1f, CS_ECIS_A=%0.1f",CS_ESIS_A,CS_ECIS_A);
  ftexto(0*side+xoffset,5*side+yoffset,0,text);
  tcolor = .25+.75*(FIDDLE_PARAMETER==6);glColor3f(tcolor,tcolor,tcolor);
  sprintf(text,"CS_ESEC_A=%0.1f, CS_ECEC_A=%0.1f",CS_ESEC_A,CS_ECEC_A);
  ftexto(0*side+xoffset,6*side+yoffset,0,text);
  tcolor = .25+.75*(FIDDLE_PARAMETER==7);glColor3f(tcolor,tcolor,tcolor);
  sprintf(text,"CS_ESIC_A=%0.1f, CS_ECIC_A=%0.1f",CS_ESIC_A,CS_ECIC_A);
  ftexto(0*side+xoffset,7*side+yoffset,0,text);
  tcolor = .25+.75*(FIDDLE_PARAMETER==8);glColor3f(tcolor,tcolor,tcolor);
  sprintf(text,"CS_ISES_G=%0.1f, CS_ICES_G=%0.1f",CS_ISES_G,CS_ICES_G);
  ftexto(0*side+xoffset,8*side+yoffset,0,text);
  tcolor = .25+.75*(FIDDLE_PARAMETER==9);glColor3f(tcolor,tcolor,tcolor);
  sprintf(text,"CS_ISIS_G=%0.1f, CS_ICIS_G=%0.1f",CS_ISIS_G,CS_ICIS_G);
  ftexto(0*side+xoffset,9*side+yoffset,0,text);
  tcolor = .25+.75*(FIDDLE_PARAMETER==10);glColor3f(tcolor,tcolor,tcolor);
  sprintf(text,"CS_ISEC_G=%0.1f, CS_ICEC_G=%0.1f",CS_ISEC_G,CS_ICEC_G);
  ftexto(0*side+xoffset,10*side+yoffset,0,text);
  tcolor = .25+.75*(FIDDLE_PARAMETER==11);glColor3f(tcolor,tcolor,tcolor);
  sprintf(text,"CS_ISIC_G=%0.1f, CS_ICIC_G=%0.1f",CS_ISIC_G,CS_ICIC_G);
  ftexto(0*side+xoffset,11*side+yoffset,0,text);
  tcolor = .25 + .75*(FIDDLE_PARAMETER==12);glColor3f(tcolor,tcolor,tcolor);
  sprintf(text,"AMPA_DIA=%0.1f,P_AMPA=%0.1f",AMPA_DIA,P_AMPA);
  ftexto(0*side+xoffset,12*side+yoffset,0,text);
  tcolor = .25 + .75*(FIDDLE_PARAMETER==13);glColor3f(tcolor,tcolor,tcolor);
  sprintf(text,"NMDA_DIA=%0.1f,P_NMDA=%0.1f",NMDA_DIA,P_NMDA);
  ftexto(0*side+xoffset,13*side+yoffset,0,text);
  tcolor = .25 + .75*(FIDDLE_PARAMETER==14);glColor3f(tcolor,tcolor,tcolor);
  sprintf(text,"GABA_DIA=%0.1f,P_GABA=%0.1f",GABA_DIA,P_GABA);
  ftexto(0*side+xoffset,14*side+yoffset,0,text);
  tcolor = .25 + .75*(FIDDLE_PARAMETER==15);glColor3f(tcolor,tcolor,tcolor);
  sprintf(text,"STD_VIEW=%0.1f,PTREE_VIEW=%0.1f",STD_VIEW,PTREE_VIEW);
  ftexto(0*side+xoffset,15*side+yoffset,0,text);
  menupos=17;
  tcolor = .25+.75*(FIDDLE_PARAMETER==menupos);glColor3f(tcolor,tcolor,tcolor);
  sprintf(text,"diarowmin=%d,diacolmin=%d",GLOBAL_BLOCK_ROW_DIA_MIN,GLOBAL_BLOCK_COL_DIA_MIN);
  ftexto(0*side+xoffset,menupos*side+yoffset,0,text);
  menupos=18;
  tcolor = .25+.75*(FIDDLE_PARAMETER==menupos);glColor3f(tcolor,tcolor,tcolor);
  sprintf(text,"dtmax=%0.2f",GLOBAL_DTmax);
  ftexto(0*side+xoffset,menupos*side+yoffset,0,text);
  menupos=19;
  tcolor = .25+.75*(FIDDLE_PARAMETER==menupos);glColor3f(tcolor,tcolor,tcolor);
  sprintf(text,"LRtoAMPA=%0.5f",LR_TO_AMPA);
  ftexto(0*side+xoffset,menupos*side+yoffset,0,text);
  menupos=20;
  tcolor = .25+.75*(FIDDLE_PARAMETER==menupos);glColor3f(tcolor,tcolor,tcolor);
  sprintf(text,"STEPS_PER_DRAW=%d",STEPS_PER_DRAW);
  ftexto(0*side+xoffset,menupos*side+yoffset,0,text);
  menupos=21;
  tcolor = .25+.75*(FIDDLE_PARAMETER==menupos);glColor3f(tcolor,tcolor,tcolor);
  sprintf(text,"DRAW_FLAG=%d,%d",DRAW_FLAG,DRAW_FLAG2);
  ftexto(0*side+xoffset,menupos*side+yoffset,0,text);
  menupos=22;
  tcolor = .25+.75*(FIDDLE_PARAMETER==menupos);glColor3f(tcolor,tcolor,tcolor);
  sprintf(text,"SMOOTHER=%d",GLOBAL_SPACE_SMOOTHER);
  ftexto(0*side+xoffset,menupos*side+yoffset,0,text);
  menupos=23;
  tcolor = .25+.75*(FIDDLE_PARAMETER==menupos);glColor3f(tcolor,tcolor,tcolor);
  sprintf(text,"FIDDLE_ROW=%d, FIDDLE_COL=%d, %s",FIDDLE_ROW,FIDDLE_COL,INPUT_IGNORE_FLAG ? "input ignored" : "");
  ftexto(0*side+xoffset,menupos*side+yoffset,0,text);
  menupos=24;
  tcolor = .25+.75*(FIDDLE_PARAMETER==menupos);glColor3f(tcolor,tcolor,tcolor);
  sprintf(text,"INPUT_CONTRAST=%0.3f,OTHERLAYER_INPUTRATE=%0.3f",INPUT_CONTRAST,OTHERLAYER_INPUTRATE);
  ftexto(0*side+xoffset,menupos*side+yoffset,0,text);
}

void keyPressed(unsigned char key, int x, int y) 
{
  int nr=0;int nc=0;
  struct neuron *n=NULL;
  /* The function called whenever a normal key is pressed. */
  usleep(100); /* still not sure why this is called */
  switch (key) {    
  case ESCAPE_KEY: /* stop everything */
    glutDestroyWindow(GLUTWINDOWNUMBER); /* shut down our window */
    printf("(x,y,z) = (%f,%f,%f)\n",xdepth,ydepth,zdepth);
    exit(EXIT_SUCCESS); /* exit the program...normal termination. */
    break; /* just in case */
  case SPACE_KEY: /* stop/restart temporarily */
    RYangle = 0; RXangle = 0; RYdangle = 0; RXdangle = 0;
    STEPS_PER_DRAW=8*!STEPS_PER_DRAW;if (!STEPS_PER_DRAW){ printf("PROCESS HALTED\n");}
    break;
  case 'w': FIDDLE_PARAMETER++; if (FIDDLE_PARAMETER > 24){ FIDDLE_PARAMETER = -4;} break;
  case 's': FIDDLE_PARAMETER--; if (FIDDLE_PARAMETER < -4){ FIDDLE_PARAMETER = 24;} break;
  case 'd': 
    switch (FIDDLE_PARAMETER){
    case -5: GRAYSCALE = !GRAYSCALE; break;
    case -4:     
      for (nr=0;nr<SYSTEM_ROW_SIZE;nr++){ for (nc=0;nc<SYSTEM_COL_SIZE;nc++){ 
	n = nget(GLOBAL_Nra,nr,nc);
	*(n->sA)=0;*(n->sN)=0;*(n->sG)=0;}}
      break;
    case -3: GLOBAL_ODOR += GLOBAL_ODOR_BACON;break;
    case -2: LGN_STRENGTH *= 1.1; break;
    case -1: LGN_BACKRATE *= 1.1; break;
    case 0: CS_ESLR *= 1.1; break;
    case 1: CS_ISLR *= 1.1; break;
    case 2: CS_LRES *= 1.1; break;
    case 3: CS_LRIS *= 1.1; break;
    case 4: CS_ESES_A *= 1.1; break;
    case 5: CS_ESIS_A *= 1.1; break;
    case 6: CS_ESEC_A *= 1.1; break;
    case 7: CS_ESIC_A *= 1.1; break;
    case 8: CS_ISES_G *= 1.1; break;
    case 9: CS_ISIS_G *= 1.1; break;
    case 10: CS_ISEC_G *= 1.1; break;
    case 11: CS_ISIC_G *= 1.1; break;
    case 12: AMPA_DIA += 1; if (AMPA_DIA < 0){ AMPA_DIA = 0;} break;
    case 13: NMDA_DIA += 1; if (NMDA_DIA < 0){ NMDA_DIA = 0;} break;
    case 14: GABA_DIA += 1; if (GABA_DIA < 0){ GABA_DIA = 0;} break;
    case 15: STD_VIEW *= 1.1; break;
    case 17: GLOBAL_BLOCK_ROW_DIA_MIN *= 2; break;
    case 18: GLOBAL_DTmax *= 2; break;
    case 19: if (LR_TO_AMPA < 0.5){ LR_TO_AMPA *= 2;} else{ LR_TO_AMPA = 1 - (1-LR_TO_AMPA)/2.0;} break;
    case 20: STEPS_PER_DRAW *= 2; if (STEPS_PER_DRAW==0){ STEPS_PER_DRAW=1;} break;
    case 21: DRAW_FLAG++; break;
    case 22: GLOBAL_SPACE_SMOOTHER++; break;
    case 23: FIDDLE_ROW += 1; break;
    case 24: INPUT_CONTRAST*=1.1; break;
    default: break;}
    break;
  case 'D':
    switch (FIDDLE_PARAMETER){
    case -3: break;
    case -2: OTHERLAYER_STRENGTH *= 1.1; break;
    case -1: OTHERLAYER_BACKRATE *= 1.1; break;
    case 0: CS_ECLR *= 1.1; break;
    case 1: CS_ICLR *= 1.1; break;
    case 2: CS_LREC *= 1.1; break;
    case 3: CS_LRIC *= 1.1; break;
    case 4: CS_ECES_A *= 1.1; break;
    case 5: CS_ECIS_A *= 1.1; break;
    case 6: CS_ECEC_A *= 1.1; break;
    case 7: CS_ECIC_A *= 1.1; break;
    case 8: CS_ICES_G *= 1.1; break;
    case 9: CS_ICIS_G *= 1.1; break;
    case 10: CS_ICEC_G *= 1.1; break;
    case 11: CS_ICIC_G *= 1.1; break;
    case 12: P_AMPA += 0.1; if (P_AMPA>1){ P_AMPA=1;} break;
    case 13: P_NMDA += 0.1; if (P_NMDA>1){ P_NMDA=1;} break;
    case 14: P_GABA += 0.1; if (P_GABA>1){ P_GABA=1;} break;
    case 15: PTREE_VIEW *= 1.1; break;
    case 17: GLOBAL_BLOCK_COL_DIA_MIN *= 2; break;
    case 21: DRAW_FLAG2++; break;
    case 23: FIDDLE_COL += 1; break;
    case 24: OTHERLAYER_INPUTRATE*=1.1; break;
    default: break;}
    break;
  case 'a':
    switch (FIDDLE_PARAMETER){
    case -5: GRAYSCALE = !GRAYSCALE; break;
    case -4: break;
    case -3: GLOBAL_ODOR -= GLOBAL_ODOR_BACON;break;
    case -2: LGN_STRENGTH /= 1.1; break;
    case -1: LGN_BACKRATE /= 1.1; break;
    case 0: CS_ESLR /= 1.1; break;
    case 1: CS_ISLR /= 1.1; break;
    case 2: CS_LRES /= 1.1; break;
    case 3: CS_LRIS /= 1.1; break;
    case 4: CS_ESES_A /= 1.1; break;
    case 5: CS_ESIS_A /= 1.1; break;
    case 6: CS_ESEC_A /= 1.1; break;
    case 7: CS_ESIC_A /= 1.1; break;
    case 8: CS_ISES_G /= 1.1; break;
    case 9: CS_ISIS_G /= 1.1; break;
    case 10: CS_ISEC_G /= 1.1; break;
    case 11: CS_ISIC_G /= 1.1; break;
    case 12: AMPA_DIA -= 1; if (AMPA_DIA < 0){ AMPA_DIA = 0;} break;
    case 13: NMDA_DIA -= 1; if (NMDA_DIA < 0){ NMDA_DIA = 0;} break;
    case 14: GABA_DIA -= 1; if (GABA_DIA < 0){ GABA_DIA = 0;} break;
    case 15: STD_VIEW /= 1.1; break;
    case 17: GLOBAL_BLOCK_ROW_DIA_MIN /= 2; break;
    case 18: GLOBAL_DTmax /= 2; break;
    case 19: if (LR_TO_AMPA <= 0.5){ LR_TO_AMPA /= 2;} else{ LR_TO_AMPA = 1 - (1-LR_TO_AMPA)*2.0;} break;
    case 20: STEPS_PER_DRAW /= 2; if (STEPS_PER_DRAW < 1){ STEPS_PER_DRAW=0;} break;
    case 21: DRAW_FLAG--; break;
    case 22: GLOBAL_SPACE_SMOOTHER--; if (GLOBAL_SPACE_SMOOTHER<0){ GLOBAL_SPACE_SMOOTHER=0;} break;
    case 23: FIDDLE_ROW -= 1; break;
    case 24: INPUT_CONTRAST/=1.1; break;
    default: break;}
    break;
  case 'A':
    switch (FIDDLE_PARAMETER){
    case -3: break;
    case -2: OTHERLAYER_STRENGTH /= 1.1; break;
    case -1: OTHERLAYER_BACKRATE /= 1.1; break;
    case 0: CS_ECLR /= 1.1; break;
    case 1: CS_ICLR /= 1.1; break;
    case 2: CS_LREC /= 1.1; break;
    case 3: CS_LRIC /= 1.1; break;
    case 4: CS_ECES_A /= 1.1; break;
    case 5: CS_ECIS_A /= 1.1; break;
    case 6: CS_ECEC_A /= 1.1; break;
    case 7: CS_ECIC_A /= 1.1; break;
    case 8: CS_ICES_G /= 1.1; break;
    case 9: CS_ICIS_G /= 1.1; break;
    case 10: CS_ICEC_G /= 1.1; break;
    case 11: CS_ICIC_G /= 1.1; break;
    case 12: P_AMPA -= 0.1; if (P_AMPA<0){ P_AMPA=0;} break;
    case 13: P_NMDA -= 0.1; if (P_NMDA<0){ P_NMDA=0;} break;
    case 14: P_GABA -= 0.1; if (P_GABA<0){ P_GABA=0;} break;
    case 15: PTREE_VIEW /= 1.1; break;
    case 17: GLOBAL_BLOCK_COL_DIA_MIN /= 2; break;
    case 21: DRAW_FLAG2--; break;
    case 23: FIDDLE_COL -= 1; break;
    case 24: OTHERLAYER_INPUTRATE/=1.1; break;
    default: break;}
    break;
  default:
    break;}
}

void specialKeyPressed(int key, int x, int y) 
{
  /* The function called whenever a special key is pressed. */
  int mod=0;
  usleep(100);
  switch (key) {    
  case GLUT_KEY_PAGE_UP:
    zdepth /= 1.05;
    break;
  case GLUT_KEY_PAGE_DOWN:
    zdepth *= 1.05;
    break;
  case GLUT_KEY_UP:
    mod = glutGetModifiers();
    if (mod==GLUT_ACTIVE_SHIFT){ RXdangle -= 0.1;}
    else{ ydepth+= -0.05*zdepth;}
    break;
  case GLUT_KEY_DOWN:
    mod = glutGetModifiers();
    if (mod==GLUT_ACTIVE_SHIFT){ RXdangle += 0.1;}
    else{ ydepth-= -0.05*zdepth;}
    break;
  case GLUT_KEY_LEFT:
    mod = glutGetModifiers();
    if (mod==GLUT_ACTIVE_SHIFT){ RYdangle -= 0.1;}
    else{ xdepth-= -0.05*zdepth;}
    break;
  case GLUT_KEY_RIGHT:
    mod = glutGetModifiers();
    if (mod==GLUT_ACTIVE_SHIFT){ RYdangle += 0.1;}
    else{ xdepth+= -0.05*zdepth;}
    break;
  case GLUT_KEY_END:
    mod = glutGetModifiers();
    if (mod==GLUT_ACTIVE_SHIFT){ 
      if (PTREE_BOTHER){ ptreereset(GLOBAL_PTREE); /* ptreerate(GLOBAL_PTREE); ptreedump_starter(GLOBAL_PTREE,NULL,2,0,0,0,+1,-1); */} 
      if (CLOSET_BOTHER){  /* ptreerate(GLOBAL_CLOSET->p); */ ptreedump_starter(GLOBAL_CLOSET->p,NULL,2,0,0,0,+1,-1);}
      if (YGGDRASIL_BOTHER){ 
	if (GLOBAL_YGGDRASIL->p!=NULL){
	  ptreerate(GLOBAL_YGGDRASIL->p); ptreedump_starter(GLOBAL_YGGDRASIL->p,"yggdrasil_p",2,0,0,0,+1,-1);}
	if (GLOBAL_YGGDRASIL->pp!=NULL){ 
	  ptreerate(GLOBAL_YGGDRASIL->pp); ptreedump_starter(GLOBAL_YGGDRASIL->pp,"yggdrasil_pp",2,0,0,0,+1,-1);}}
      if (BONSAI_BOTHER){ bonsaidump(GLOBAL_BONSAI,"bonsai",0);bonsaidump(GLOBAL_BONSAI,"bonsai",1);bonsaidump(GLOBAL_BONSAI,"bonsai",2);}
      dumpoutput("d_output");}
    else{ if (!STEPS_PER_DRAW){ computestep();}}
    break;
  case GLUT_KEY_HOME:
    mod = glutGetModifiers();
    if (mod==GLUT_ACTIVE_SHIFT){INPUT_IGNORE_FLAG=!INPUT_IGNORE_FLAG; if (INPUT_IGNORE_FLAG && LGN_BOTHER){ lgnresetrates(GLOBAL_LGN);}}
    else{ GLOBAL_DRAW_FLAG = !GLOBAL_DRAW_FLAG;}
    break;
  default:
    break;}
}

GLvoid DrawGLScene(GLvoid)
{
  int index=0;
  for (index=0;index<STEPS_PER_DRAW;index++){ if (!RUN_DONE){ computestep();}}
  /* draw stuff */
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // Clear The Screen And The Depth Buffer
  glLoadIdentity(); // Reset The View
  glTranslatef(xdepth,ydepth,zdepth);
/*   glRotatef(RYangle,0,1,0); */
/*   glRotatef(RXangle,1,0,0); */
  if (GLOBAL_DRAW_FLAG){
    if (CORTEX_BOTHER){ /* DrawLR(2,1,6.0); */ DrawNra(2,1,3.25);}
    if (LGN_BOTHER && !CORTEX_BOTHER){ Drawlgn(3,3,1,0);}
    if (RTC_BOTHER){ Drawrtc(0.5,0.5,3,-0.5);}
    if (STROBETRACE_BOTHER){ Drawstrobetrace(1,1.9,-.5);}
    if (TUNINGCURVE_BOTHER){ Drawtuningcurve(0.5,0.5,2,-0.5);}
    if (LMITRI_BOTHER){ Drawlmitri(0.5,0.5,2,-0.5);}
    if (PTREE_BOTHER){ Drawptree(GLOBAL_PTREE,2,1,0);}
    if (CLOSET_BOTHER){ Drawptree(GLOBAL_CLOSET->p,2,1,0);}
    if (YGGDRASIL_BOTHER){ Drawptree(GLOBAL_YGGDRASIL->p,2,1,0); Drawptree(GLOBAL_YGGDRASIL->pp,2,1,-3);}
    if (BONSAI_BOTHER){ Drawptree(GLOBAL_BONSAI->p[0],2,1,0); Drawptree(GLOBAL_BONSAI->p[1],2,1,-3);}
    if (HYDRA_BOTHER){ Drawptree(GLOBAL_HYDRA->pjuston_1,2,1,0); Drawptree(GLOBAL_HYDRA->pstayon_1,2,1,-3); Drawptree(GLOBAL_HYDRA->pjustoff_1,2,1,-6); Drawptree(GLOBAL_HYDRA->pjuston_2,2,1,-9); Drawptree(GLOBAL_HYDRA->pstayon_2,2,1,-12); Drawptree(GLOBAL_HYDRA->pjustoff_2,2,1,-15); Drawptree(GLOBAL_HYDRA->pstayoff,2,1,-18);}
    if (LYAPUNOV_BOTHER){ Drawlyapunov(2,1,-1);}
    if (POWER_BOTHER){ Drawpower(2,3,0);}
    if (TAOF_BOTHER){ Drawtaof(2,3,0);}}
  Drawmenu(.1,-1,1);
  glutSwapBuffers(); /* since double buffered */
  RYangle += RYdangle;
  RXangle += RXdangle;
}

#endif /* NOTGRAPH */

