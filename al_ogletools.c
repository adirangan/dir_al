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

GLvoid Drawra(void *ra,char *type,int rows,int cols,int coloringtype,int stduse,double max,double min,int gsuse,char *text,double xside,double yside,double xoffset,double yoffset)
{
  /* draws a standard *ra */
  char text2[64];
  int i=0,j=0;
  double varmax=0,varmin=0,varmean=0,varstd=0,vMAX=0,vMIN=0;
  double *dra=NULL,*dra2=NULL;
  int *ira=NULL,*ira2=NULL;
  double xord=0,yord=0;
  double rcolor=0,gcolor=0,bcolor=0; 
  int clear_flag=0;
  if (strcmp(type,"double")==0){ 
    dra = (double *)ra;
    if (gsuse){ dra2 = spacesmear(dra,rows,cols,gsuse*GLOBAL_SPACE_SMOOTHER); clear_flag=1;} else{ dra2 = dra; clear_flag=0;}
    stats("double",dra2,rows*cols,&varmax,&varmin,&varmean,&varstd);
    if (stduse){ vMAX=varmean+STD_VIEW*varstd; vMIN=varmean-STD_VIEW*varstd;} else{ vMAX=varmax;vMIN=varmin;}
    if (max>min){ vMAX = max; vMIN = min;}
    for (i=0;i<rows;i++){ for (j=0;j<cols;j++){
      xord = (j+0.5)/(double)cols; yord = -(i+0.5)/(double)rows;
      colorscale(coloringtype,dra2[i+j*rows],vMAX,vMIN,&rcolor,&gcolor,&bcolor);
      glrect(xord*xside+xoffset,yord*yside+yoffset,xside/(double)cols,yside/(double)rows,rcolor,gcolor,bcolor);}}
    if (clear_flag){ tfree(dra2); dra2=NULL;}}
  else if (strcmp(type,"int")==0){ 
    ira = (int *)ra; ira2 = ira; clear_flag=0;
    stats("int",ira2,rows*cols,&varmax,&varmin,&varmean,&varstd);
    if (stduse){ vMAX=varmean+STD_VIEW*varstd; vMIN=varmean-STD_VIEW*varstd;} else{ vMAX=varmax;vMIN=varmin;}
    if (max>min){ vMAX = max; vMIN = min;}
    for (i=0;i<rows;i++){ for (j=0;j<cols;j++){
      xord = j/(double)cols; yord = -i/(double)rows;
      colorscale(coloringtype,ira2[i+j*rows],vMAX,vMIN,&rcolor,&gcolor,&bcolor);
      glrect(xord*xside+xoffset,yord*yside+yoffset,xside/(double)cols,yside/(double)rows,rcolor,gcolor,bcolor);}}
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

GLvoid Drawconnectivity(struct neuronarray *Nra,double xside,double yside,double xoffset,double yoffset)
{
  int nt=0,nr=0,nt2=0,nr2=0,tab=0,tab2=0;
  double *cra=NULL,*maxsra=NULL;
  struct neuron *s=NULL,*n=NULL;
  cra = (double *) tcalloc(Nra->lt*Nra->lt,sizeof(double));
  maxsra = (double *) tcalloc(GLOBAL_INDEXING_sra_LENGTH,sizeof(double));
  tab=0;
  for (nt=0;nt<Nra->ntypes;nt++){ 
    if (Nra->lengthra[nt]>0){
      for (nr=0;nr<Nra->lengthra[nt];nr++){
	s = nget(Nra,nt,nr);
	tab2=0;
	for (nt2=0;nt2<Nra->ntypes;nt2++){ 
	  if (Nra->lengthra[nt2]>0){
	    for (nr2=0;nr2<Nra->lengthra[nt2];nr2++){
	      n = nget(Nra,nt2,nr2);
	      if (slink(s,n,NULL,maxsra)){ cra[tab+s->index + (tab2+n->index)*Nra->lt]=1;}
	      else /* if not connected */{ cra[tab+s->index + (tab2+n->index)*Nra->lt]=0;}}}
	  tab2 += Nra->lengthra[nt2];}}}
    tab += Nra->lengthra[nt];}
  Drawra(&(cra[0]),"double",Nra->lt,Nra->lt,7,1,0,0,0,"connectivity",xside,yside,xoffset,yoffset);
  tfree(maxsra);maxsra=NULL;
  tfree(cra);cra=NULL;
}

GLvoid DrawNra(struct neuronarray *Nra,double xside,double yside,double xoffset,double yoffset)
{
  int nt=0,nv=0;
  char text[128];
  for (nt=0;nt<Nra->ntypes;nt++){
    nv = periodize(DRAW_FLAG,0,Nra->nvars); sprintf(text,"%s_%s",GLOBAL_TYPENAMES[nt],GLOBAL_VARNAMES[nv]);
    Drawra(Nra->vrarara[nt][nv],"double",1,Nra->lengthra[nt],7,0,0,0,0,text,xside,yside,xoffset,yoffset-yside*nt);}
}

GLvoid Drawsnx(double side,double xoffset,double yoffset)
{
  double *ra=NULL,*ra2=NULL;
  int total_Vs_states=0, total_nAch_states=0, total_states=0;
  double input=0;
  int ns=0,ns1=0,ns2=0,nv1=0,ng1=0,nv2=0,ng2=0;
  switch (GLOBAL_SNX_VERSION){
  case 0: case 2:
    total_Vs_states=GLOBAL_SNX_NSTATES_[VARNAME_REGISTRY_snx_Vs];
    total_nAch_states=GLOBAL_SNX_NSTATES_[VARNAME_REGISTRY_snx_nAch]; 
    total_states=total_Vs_states*total_nAch_states;
    input=CS_ORN_[TYPENAME_REGISTRY_snx_PN];
    ra = (double *) tcalloc(total_states*total_states,sizeof(double));
    for (ns1=0;ns1<total_Vs_states;ns1++){ for (ns2=0;ns2<total_nAch_states;ns2++){
      ns = ns1+ns2*total_Vs_states;
      snx_statematrix_0(ns1,total_Vs_states,ns2,total_nAch_states,input,&(ra[0+ns*total_states]));}}
    switch (abs(DRAW_FLAG2%2)){
    case 0: 
      Drawra(ra,"double",total_states,total_states,7,0,2.0/(double)total_states,0.0,0,"snx",side,side,xoffset,yoffset); 
      break;
    case 1: 
      ra2 = (double *) tcalloc(total_states*total_states,sizeof(double));
      for (nv1=0;nv1<total_Vs_states;nv1++){ for (ng1=0;ng1<total_nAch_states;ng1++){
	ns1 = nv1 + ng1*total_Vs_states;
	for (nv2=0;nv2<total_Vs_states;nv2++){ for (ng2=0;ng2<total_nAch_states;ng2++){
	  ns2 = nv2 + ng2*total_Vs_states;
	  ra2[(ng1+nv1*total_nAch_states) + (ng2+nv2*total_nAch_states)*total_states] = ra[ns1+ns2*total_states];}}}}
      Drawra(ra2,"double",total_states,total_states,7,0,2.0/(double)total_states,0.0,0,"snx",side,side,xoffset,yoffset); 
      tfree(ra2);ra2=NULL;
      break;
    default: break;}
    tfree(ra);ra=NULL;
    break;
  case 1: 
    total_Vs_states=GLOBAL_SNX_NSTATES_[VARNAME_REGISTRY_snx_Vs]; total_states=total_Vs_states;
    input=CS_ORN_[TYPENAME_REGISTRY_snx_PN];
    ra = (double *) tcalloc(total_states*total_states,sizeof(double));
    for (ns1=0;ns1<total_Vs_states;ns1++){
      ns = ns1;
      snx_statematrix_1(ns1,total_Vs_states,input,&(ra[0+ns*total_states]),NULL);}
    Drawra(ra,"double",total_states,total_states,7,0,2.0/(double)total_states,0,0,"snx",side,side,xoffset,yoffset);
    break;
  default: break;}
}

GLvoid Drawcra(struct neuronarray *Nra,double side,double xoffset,double yoffset)
{
  char text[32];
  int nt=0,nr=0;
  double rcolor=0,gcolor=0,bcolor=0;
  double xord=0,yord=0,xord2=0,yord2=0;
  double arc = 2*PI/(double)Nra->gli->nclusters;
  double scale = STD_VIEW*0.01;
  double ascale = 1;
  struct neuron *n=NULL,*n2=NULL;
  struct llist *L=NULL;
  struct litem *l=NULL;
  int nv1 =   abs(DRAW_FLAG)%GLOBAL_NVARS, nv2 = 0, log_flag=1;
  double max1=0,min1=0;
  switch (GLOBAL_NEURON_MODEL){
  case 0: case 1: case 2: case 3: case 4: nv2 = VARNAME_REGISTRY_spike_flag; break;
  case 5: nv2 = VARNAME_REGISTRY_mainak_spike_flag; break;
  case 6: nv2 = VARNAME_REGISTRY_wilson_spike_flag; break;
  case 7: nv2 = VARNAME_REGISTRY_snx_spike_flag; break;
  default: break;}
  max1 = GLOBAL_POWER_maxra_[nv1], min1 = GLOBAL_POWER_minra_[nv1];
  for (nt=0;nt<Nra->ntypes;nt++){
    for (nr=0;nr<Nra->lengthra[nt];nr++){
      n=nget(Nra,nt,nr);
      xord=Nra->gli->n2x[nt][nr]; yord=Nra->gli->n2y[nt][nr];
      colorscale(7,*(n->vpra[nv1]),max1,min1,&rcolor,&gcolor,&bcolor);
      gldot(3+n->type,xord*side+xoffset,yord*side+yoffset,side*arc*scale,rcolor,gcolor,bcolor);
      if (DRAW_FLAG2%2){
      if (log_flag==0){ colorscale(7,*(n->vpra[nv2]),GLOBAL_POWER_maxra_[nv2],GLOBAL_POWER_minra_[nv2],&rcolor,&gcolor,&bcolor);}
      else if (log_flag==1){ colorscale(7,log(*(n->vpra[nv2])),0,-32,&rcolor,&gcolor,&bcolor);}
      glring(3+n->type,xord*side+xoffset,yord*side+yoffset,side*arc*scale,2*side*arc*scale,rcolor,gcolor,bcolor);
	L=llistmake();
	llitem2llist(n->sparse_link,L); l=L->first;
	while(l!=NULL){
	  n2=(struct neuron *)l->item;
	  xord2=Nra->gli->n2x[n2->type][n2->index];
	  yord2=Nra->gli->n2y[n2->type][n2->index];
	  glarc(8,xord*side+xoffset,yord*side+yoffset,xord2*side+xoffset,yord2*side+yoffset,-0.01,ascale*side*((n->sparse_out-n2->sparse_in)%Nra->gli->nclusters+0.5)/(double)Nra->gli->nclusters,rcolor,gcolor,bcolor);
	  l=l->child;}
	llisttfree(L); L=NULL;}}}
  glColor3f(1,1,1); sprintf(text,"%s,%s",GLOBAL_VARNAMES[nv1],GLOBAL_VARNAMES[nv2]); ftexto(xoffset,yoffset-0.5*side,0,text);
}

GLvoid Drawstra(struct strobe **stra,int ralength,double xside,double yside,double xoffset,double yoffset,double max,double min,double rcolor,double gcolor,double bcolor,int drawvsplot)
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
      if (!automatic_color){
	glBegin(GL_LINES);
	glColor3f(1,1,1);
	glVertex3f(xoffset,0.5*yside+yoffset,0);
	glVertex3f(xoffset+xside,0.5*yside+yoffset,0);
	glVertex3f(xoffset,0*yside+yoffset,0);
	glVertex3f(xoffset+xside,0*yside+yoffset,0);
	glVertex3f(xoffset,1*yside+yoffset,0);
	glVertex3f(xoffset+xside,1*yside+yoffset,0);
	glVertex3f(xoffset,0*yside+yoffset,0);
	glVertex3f(xoffset,1*yside+yoffset,0);
	glVertex3f(xoffset+xside,0*yside+yoffset,0);
	glVertex3f(xoffset+xside,1*yside+yoffset,0);
	glEnd();}
      glBegin(GL_LINE_STRIP);
      for (nt=sttab;nt<sttab+stlength;nt++){
	nt2 = periodize(nt,0,stlength);
	xord = (double)(nt-sttab+0.5)/(double)stlength;
	if (automatic_color){ 
	  colorscale(7,stra[na]->data[nt2],max2,min2,&rcolor2,&gcolor2,&bcolor2);
	  yord = (double)(-na+0.5)/(double)ralength;}
	else{ yord = (stra[na]->data[nt2]-min2)/(max2-min2);}
	glColor3f(rcolor2,gcolor2,bcolor2);
	glVertex3f(xord*xside+xoffset,yord*yside+yoffset,0);}
      glEnd();}
    if (stra[0]->lpower_bother){
      Drawra(stra[0]->lpowerdata,"double",stra[0]->lpower_window_length,stra[0]->lpower_length/stra[0]->lpower_window_length,7,1,0,0,0,NULL,xside,yside,xoffset,yoffset+yside);}
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
	    glrect(xord*xside+xoffset,yord*yside+yoffset,xside/(double)stlength,yside/(double)ralength,rcolor2,gcolor2,bcolor2);}
	  else{
	    yord = (stra[na]->cycledata[nt]-min2)/(max2-min2);
	    glrect(xord*xside+xoffset,yord*yside+yoffset+yside,xside/(double)stlength,yside/(double)ralength,rcolor2,gcolor2,bcolor2);}}}}}
  else if (drawvsplot==1){
    ra = (double *) tcalloc(ralength*stlength,sizeof(double));
    for (na=0;na<ralength;na++){ for (nt=sttab;nt<sttab+stlength;nt++){
      nt2 = periodize(nt,0,stlength);
      ra[na + (nt-sttab)*ralength] = stra[na]->data[nt2];}}
    Drawra(ra,"double",ralength,stlength,7,1,max,min,0,NULL,xside,yside,xoffset,yoffset);
    if (stcycle_bother && stcyclenum>0){
      for (na=0;na<ralength;na++){ for (nt=0;nt<stlength;nt++){
	ra[na + nt*ralength] = stra[na]->cycledata[nt]/stcyclenum;}}
      Drawra(ra,"double",ralength,stlength,7,1,0,0,0,NULL,xside,yside,xoffset,yoffset+yside);}
    if (stra[0]->lpower_bother){
      Drawra(stra[0]->lpowerdata,"double",stra[0]->lpower_window_length,stra[0]->lpower_length/stra[0]->lpower_window_length,7,1,0,0,0,NULL,xside,yside,xoffset,yoffset+yside);}
    tfree(ra);}
}

GLvoid Drawpower(struct power *p,double xside,double yside,double xoffset,double yoffset)
{
  int nt=0,nv=0;
  int drawvsplot = DRAW_FLAG2%2;
  int plot_flag = DRAW_FLAG%(p->indexing_nvar_length + 3);
  char text[128];
  if (plot_flag<p->indexing_nvar_length){
    for (nt=0;nt<p->indexing_ntype_length;nt++){
      nv = maximum(0,minimum(p->indexing_nvar_length-1,plot_flag)); 
      sprintf(text,"%s_%s tab %d",GLOBAL_TYPENAMES[p->indexing_ntype_checkout[nt]],GLOBAL_VARNAMES[p->indexing_nvar_checkout[nv]],p->strarara[nt][nv][0]->tab);
      glColor3f(1,1,1);ftexto(xoffset,yoffset-yside*nt*1.5+0.2,0,text);
      if (p->strarara[nt][nv]!=NULL){
	Drawstra(p->strarara[nt][nv],p->Nra->lengthra[p->indexing_ntype_checkout[nt]],xside,yside,xoffset,yoffset-yside*nt*1.5,p->maxra[p->indexing_nvar_checkout[nv]],p->minra[p->indexing_nvar_checkout[nv]],-1,1,0,drawvsplot);
	if (drawvsplot){ loglogra(p->vpowrarara[nt][nv],"double",p->length,0,0,xside,yside,xoffset+xside,yoffset-yside*nt*1.5,1,1,1);}
	else{ Plotra(p->vpowrarara[nt][nv],"double",p->length,0,0,xside,yside,xoffset+xside,yoffset-yside*nt*1.5,1,1,1);}}}}
  else{
    switch (plot_flag-p->indexing_nvar_length){
    case 0: sprintf(text,"firstvar");
      if (drawvsplot){ loglogra(p->firstvarpower,"double",p->length,0,0,xside,yside,xoffset,yoffset,1,1,1);}
      else{ Plotra(p->firstvarpower,"double",p->length,0,0,xside,yside,xoffset,yoffset,1,1,1);} break;
    case 1: sprintf(text,"ac"); 
      if (p->autocorrelation!=NULL){
	if (drawvsplot){ loglogra(&(p->autocorrelation[minimum(1,p->length-1)]),"double",p->length-1,0,0,xside,yside,xoffset,yoffset,1,1,1);}
	else{ Plotra(&(p->autocorrelation[minimum(1,p->length-1)]),"double",p->length-1,0,0,xside,yside,xoffset,yoffset,1,1,1);}}
      break;
    case 2: sprintf(text,"xc"); 
      if (p->crosscorrelation!=NULL){
	if (drawvsplot){ loglogra(&(p->crosscorrelation[minimum(1,p->length-1)]),"double",p->length-1,0,0,xside,yside,xoffset,yoffset,1,1,1);}
	else{ Plotra(&(p->crosscorrelation[minimum(1,p->length-1)]),"double",p->length-1,0,0,xside,yside,xoffset,yoffset,1,1,1);}}
      break;
    default: break;}
    ftexto(xoffset,yoffset+0.1,0,text);}
}

GLvoid Drawrho(struct rho *r,double xside,double yside,double xoffset,double yoffset)
{
  int verbose=0;
  char text[64];
  int nr=0,nb=0;
  int nv_draw=abs(DRAW_FLAG)%r->indexing_nvar_length;
  int tab_draw = r->indexing_nvar_checkout[nv_draw];
  int *ira=NULL,*nbinra=NULL;
  double *ra=NULL,*ra_diff=NULL;
  double rcolor=0,gcolor=0,bcolor=0;
  if (verbose){ printf(" %% [entering Drawrho], nv_draw %d, tab_draw %d\n",nv_draw,tab_draw);}
  ira = (int *) tcalloc(r->indexing_nvar_length,sizeof(int));
  nbinra = (int *) tcalloc(r->indexing_nvar_length,sizeof(int));
  ra_diff = (double *) tcalloc(r->nbinra[tab_draw],sizeof(double));
  for (nb=0;nb<r->indexing_nvar_length;nb++){ nbinra[nb] = r->nbinra[r->indexing_nvar_checkout[nb]];}
  for (nr=0;nr<r->length;nr++){
    if (verbose){ printf(" %% \n");}
    if (verbose){ printf(" %% drawing r->rhora[%d]\n",nr); raprintf(r->rhora[nr],"double",1,r->nbins,"rhora: ");}
    ra = (double *) tcalloc(r->nbinra[tab_draw],sizeof(double));
    for (nb=0;nb<r->nbins;nb++){
      if (verbose){ printf(" %% currently testing bin %d... ",nb);}
      indextract(nb,r->indexing_nvar_length,nbinra,ira);
      if (verbose){ printf(" %% nv %d, nvar %d, ",nb,r->indexing_nvar_length); raprintf(nbinra,"int",1,r->indexing_nvar_length,"nbinra: "); raprintf(ira,"int",1,r->indexing_nvar_length," ira: ");}
      ra[ira[nv_draw]] += r->rhora[nr][nb];}
    ratimesequals(ra,r->nbinra[tab_draw],1.0/maximum(1.0,(double)r->neuronllistra[nr]->length)/maximum(1.0,(double)r->total_time));
    colorscale(7,nr,r->length-1,0,&rcolor,&gcolor,&bcolor);
    Plotra(ra,"double",r->nbinra[tab_draw],1.0/(double)r->nbinra[tab_draw],0,xside,yside,xoffset,yoffset,rcolor,gcolor,bcolor);
    if (nr==0){ for (nb=0;nb<r->nbinra[tab_draw];nb++){ ra_diff[nb] += ra[nb];}}
    if (nr==1){ for (nb=0;nb<r->nbinra[tab_draw];nb++){ ra_diff[nb] -= ra[nb];}}
    tfree(ra);ra=NULL;}
  tfree(ira);ira=NULL;
  tfree(nbinra);nbinra=NULL;
  Plotra(ra_diff,"double",r->nbinra[tab_draw],0,0,xside,yside,xoffset,yoffset,1,1,1);
  tfree(ra_diff);ra_diff=NULL;
  sprintf(text,"%s,%f->-%d->-%f",GLOBAL_VARNAMES[tab_draw],r->minra[tab_draw],r->nbinra[tab_draw],r->maxra[tab_draw]);
  glColor3f(1,1,1);
  ftexto(xoffset,yoffset-0.1,0,text);
  glgrid(1,3,3,xside,yside,xoffset,yoffset+yside,0.7,0.5,0.3);
} 

GLvoid Draweventra(struct ptree *p,double side,double xoffset,double yoffset)
{
  struct llist *L=NULL;
  struct litem *l=NULL;
  struct region *r=NULL;
  int nt=0,tab=0,nr=0;
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
  glBegin(GL_LINES);
  glColor3f(1,1,1);
  for (nr=0;nr<p->nregions;nr++){
    xord = 0*side+xoffset; yord = -1*nr/(double)p->nregions*side+yoffset; glVertex3f(xord,yord,0);
    xord = 1*side+xoffset; yord = -1*nr/(double)p->nregions*side+yoffset; glVertex3f(xord,yord,0);}
  glEnd();
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
    for (nl=0;nl<=xlines;nl++){
      y1 = (double)(-nl)/(double)xlines*yside+yoffset;
      glVertex3f(x1,y1,0);
      glVertex3f(x2,y1,0);}
    y1 = -0*yside+yoffset;
    y2 = -1*yside+yoffset;
    for (nl=0;nl<=ylines;nl++){
      x1 = (double)(nl)/(double)ylines*xside+xoffset;
      glVertex3f(x1,y1,0);
      glVertex3f(x1,y2,0);}}
  else /* if (!inorout) */{
    x1 = (0-0.5/(double)xlines)*xside+xoffset;
    x2 = (1-0.5/(double)xlines)*xside+xoffset;
    for (nl=0;nl<=xlines+1;nl++){
      y1 = (double)(-nl+0.5)/(double)xlines*yside+yoffset;
      glVertex3f(x1,y1,0);
      glVertex3f(x2,y1,0);}
    y1 = -(0-0.5/(double)ylines)*yside+yoffset;
    y2 = -(1-0.5/(double)ylines)*yside+yoffset;
    for (nl=0;nl<=ylines+1;nl++){
      x1 = (double)(nl-0.5)/(double)ylines*xside+xoffset;
      glVertex3f(x1,y1,0);
      glVertex3f(x1,y2,0);}}
  glEnd();
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

GLvoid glrect(double x,double y,double rx,double ry,double rcolor,double gcolor,double bcolor)
{
  glBegin(GL_QUADS);
  glColor3f(rcolor,gcolor,bcolor);
  glVertex3f(x+rx/2,y+ry/2,0);
  glVertex3f(x+rx/2,y-ry/2,0);
  glVertex3f(x-rx/2,y-ry/2,0);
  glVertex3f(x-rx/2,y+ry/2,0);
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

GLvoid glarc(int nsides,double x1,double y1,double x2,double y2,double z,double r,double rcolor,double gcolor,double bcolor)
{
  int na=0;
  double d1=0,d2=0,d3=0;
  double x3=0,y3=0,x4=0,y4=0;
  double theta1=0,theta2=0,theta3=0;
  double dim=0;
  d1=sqrt(pow(x1-x2,2)+pow(y1-y2,2));
  x3 = (x1+x2)/2; y3 = (y1+y2)/2;
  x4 = (y1-y2)/2; y4 = (x2-x1)/2; d2 = sqrt(pow(x4,2)+pow(y4,2)); if (d2>0){ x4/=d2; y4/=d2;}
  if (r>=0){ r=maximum(d2,r);} else /* if (r<0) */{ r=minimum(-d2,r);}
  d3 = (r>0?+1:-1)*sqrt(pow(r,2)-pow(d2,2));
  x3 += d3*x4; y3 += d3*y4;
  theta1 = atan2(y1-y3,x1-x3); theta2 = atan2(y2-y3,x2-x3); 
  while (theta2<theta1){ theta2 += 2*PI;} while (theta2-theta1>2*PI){ theta2-=2*PI;}
  /* smaller of two arcs */
  if (theta2-theta1>PI){ 
    x3 -= 2*d3*x4; y3 -= 2*d3*y4;
    theta1 = atan2(y1-y3,x1-x3); theta2 = atan2(y2-y3,x2-x3); 
    while (theta2<theta1){ theta2 += 2*PI;} while (theta2-theta1>2*PI){ theta2-=2*PI;}}
  glBegin(GL_LINE_STRIP);
  for (na=0;na<=nsides;na++){
    dim = 1-(double)na/(double)nsides;
    theta3 = theta1+(double)na/(double)nsides*(theta2-theta1);
    glColor3f(rcolor*dim,gcolor*dim,bcolor*dim);
    glVertex3f(x3+fabs(r)*cos(theta3),y3+fabs(r)*sin(theta3),z);}
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

GLvoid Drawptree(int drawflag,struct ptree *p,double side,double xoffset,double yoffset)
{
  struct neuronarray *Nra=GLOBAL_Nra;
  struct region *r=NULL,*r2=NULL;
  struct llitem *l0=NULL;
  struct pnode *pn=NULL;
  int nr=0;
  double bigside=0;
  double rcolor=0,gcolor=0,bcolor=0;
  int fr=0,frmin=0,frmax=0;
  int fc=0,fcmin=0,fcmax=0;
  struct neuron *n=NULL;
  struct llist *L=NULL;
  struct litem *l=NULL;
  double xord=0,yord=0;
  double wmax=0;
  stats("int",Nra->lengthra,Nra->ntypes,NULL,NULL,&bigside,NULL); bigside*=Nra->ntypes; bigside=maximum(Nra->ntypes,bigside);
  if (p!=NULL){
    switch(abs(drawflag)%5){
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
	  glbox(n->index/bigside*side+xoffset,-n->type/bigside*side+yoffset,side/bigside,rcolor,gcolor,bcolor);
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
      Drawconnectivity(Nra,side,side,xoffset,yoffset);
      break;
    default: break;}}
}

GLvoid Drawcd(struct clusterdatara *cd,double xside,double yside,double xoffset,double yoffset)
{
  int nc=0,nv=0;
  int drawvsplot = DRAW_FLAG2%2;
  int plot_flag = DRAW_FLAG;
  double aord=0,xord=0,yord=0;
  char text[128];
  int varname_registry_lfp=0;
  for (nc=0;nc<cd->gli->nclusters;nc++){
    aord = 2*PI*(nc+0.5)/(double)cd->gli->nclusters;
    xord = xoffset + 2*xside*cos(aord);
    yord = yoffset + 2*yside*sin(aord);
    if (plot_flag>=0 && plot_flag<cd->indexing_nvar_length && cd->power_bother){
      nv = maximum(0,minimum(cd->indexing_nvar_length-1,plot_flag)); 
      sprintf(text,"c%d_%s",nc,GLOBAL_VARNAMES[cd->indexing_nvar_checkout[nv]]); 
      glColor3f(1,1,1);ftexto(xord,yord+0.2,0,text);
      if (cd->strarara[nc][nv]!=NULL){
	Drawstra(cd->strarara[nc][nv],cd->gli->cra[nc]->LN->length,xside,yside,xord,yord,cd->maxra[cd->indexing_nvar_checkout[nv]],cd->minra[cd->indexing_nvar_checkout[nv]],-1,1,0,drawvsplot);
	if (cd->vpow_bother){
	  if (drawvsplot){ loglogra(cd->vpowrarara[nc][nv],"double",cd->power_length,0,0,xside,yside,xord,yord-yside,1,1,1);}
	  else{ Plotra(cd->vpowrarara[nc][nv],"double",cd->power_length,0,0,xside,yside,xord,yord-yside,1,1,1);}}}}
    else if (plot_flag<0){
      sprintf(text,"lfp_%s",GLOBAL_VARNAMES[cd->indexing_nvar_checkout[varname_registry_lfp]]); 
      glColor3f(1,1,1);ftexto(xord,yord+0.2,0,text);
      Drawstra(&(cd->stra_lfp[nc]),1,xside,yside,xord,yord,cd->maxra[cd->indexing_nvar_checkout[varname_registry_lfp]],cd->minra[cd->indexing_nvar_checkout[varname_registry_lfp]],-1,1,0,drawvsplot);}
    else if (plot_flag>=cd->indexing_nvar_length && cd->ptree_bother){
      Drawptree(plot_flag-cd->indexing_nvar_length,cd->pra[nc],minimum(xside,yside),xord,yord);}}
}

GLvoid Drawhist(double *data,int length,double datamax,double datamin,char *text,double xside,double yside,double xoffset,double yoffset)
{
  /* actually can draw any array */
  int nr=0;
  char text2[64];
  double binmax=0,binmin=0;
  double xord1=0,yord1=0,xord2=0,yord2=0;
  stats("double",data,length,&binmax,&binmin,NULL,NULL); if (binmax<=binmin){ binmax=binmin+1;}
  sprintf(text2,"%s[%0.3f,%0.3f],[%0.3f,%0.3f]",text,datamin,datamax,binmin,binmax);
  for (nr=0;nr<length;nr++){
    xord1 = (double)(nr+0)/(double)length; xord2 = (double)(nr+1)/(double)length;
    yord1 = 0; yord2 = (data[nr]-binmin)/(binmax-binmin);
    glBegin(GL_QUADS);
    glVertex3f(xord1*xside+xoffset,yord1*yside+yoffset,0);
    glVertex3f(xord1*xside+xoffset,yord2*yside+yoffset,0);
    glVertex3f(xord2*xside+xoffset,yord2*yside+yoffset,0);
    glVertex3f(xord2*xside+xoffset,yord1*yside+yoffset,0);
    glEnd();}
  ftexto(xoffset,yoffset+0.1,0,text2);
}

GLvoid Drawcaicor(struct caicor *c,double xside,double yside,double xoffset,double yoffset)
{
  int nv=abs(DRAW_FLAG%c->indexing_nvar_length),nb=0,nb2=0;
  char text[32];
  double *temp=NULL,*temp2=NULL;
  sprintf(text,"hst_%s",GLOBAL_VARNAMES[c->indexing_nvar_checkout[nv]]); 
  glColor3f(1,1,1); Drawhist(c->hstra[nv]->data,c->hstra[nv]->nbins,c->hstra[nv]->max,c->hstra[nv]->min,text,xside,yside,xoffset,yoffset);
  sprintf(text,"h_%s",GLOBAL_VARNAMES[c->indexing_nvar_checkout[nv]]); 
  glColor3f(1,1,1); Drawhist(c->hra[nv]->data,c->hra[nv]->nbins,c->hra[nv]->max,c->hra[nv]->min,text,xside,yside,xoffset,yoffset-yside);
  temp = (double *) tcalloc(c->nbins,sizeof(double));
  for (nb=0;nb<c->nbins;nb++){ temp[nb] = c->mrara[nb+0*c->nbins+nv*2*c->nbins]/maximum(1,c->mrara[nb+1*c->nbins+nv*2*c->nbins]);}
  sprintf(text,"m0_%s",GLOBAL_VARNAMES[c->indexing_nvar_checkout[nv]]); 
  glColor3f(1,1,1); Drawhist(&(c->mrara[0+0*c->nbins+nv*2*c->nbins]),c->nbins,c->maxra[nv],c->minra[nv],text,xside,yside,xoffset+xside,yoffset-0*yside);
  sprintf(text,"m1_%s",GLOBAL_VARNAMES[c->indexing_nvar_checkout[nv]]); 
  glColor3f(1,1,1); Drawhist(&(c->mrara[0+1*c->nbins+nv*2*c->nbins]),c->nbins,c->maxra[nv],c->minra[nv],text,xside,yside,xoffset+xside,yoffset-1*yside);
  sprintf(text,"m0/m1_%s",GLOBAL_VARNAMES[c->indexing_nvar_checkout[nv]]); 
  glColor3f(1,1,1); Drawhist(temp,c->nbins,c->maxra[nv],c->minra[nv],text,xside,yside,xoffset+xside,yoffset-2*yside);
  temp2=(double *)tcalloc(c->nbins*c->nbins,sizeof(double));
  for (nb=0;nb<c->nbins;nb++){ for (nb2=0;nb2<c->nbins;nb2++){ temp2[nb+nb2*c->nbins]=temp[nb]*temp[nb2];}}
  tfree(temp);temp=NULL;
  temp = (double *) tcalloc(c->nbins*c->nbins,sizeof(double));
  for (nb=0;nb<c->nbins*c->nbins;nb++){ temp[nb] = c->mrarara[nb+0*c->nbins*c->nbins+nv*2*c->nbins*c->nbins]/maximum(1,c->mrarara[nb+1*c->nbins*c->nbins+nv*2*c->nbins*c->nbins]);}
  sprintf(text,"mra0_%s",GLOBAL_VARNAMES[c->indexing_nvar_checkout[nv]]);
  Drawra(&(c->mrarara[0+0*c->nbins*c->nbins+nv*2*c->nbins*c->nbins]),"double",c->nbins,c->nbins,7,1,0,0,0,text,xside,yside,xoffset+2*xside,yoffset+1*yside);
  sprintf(text,"mra1_%s",GLOBAL_VARNAMES[c->indexing_nvar_checkout[nv]]);
  Drawra(&(c->mrarara[0+1*c->nbins*c->nbins+nv*2*c->nbins*c->nbins]),"double",c->nbins,c->nbins,7,1,0,0,0,text,xside,yside,xoffset+2*xside,yoffset-0*yside);
  sprintf(text,"mra0/mra1_%s",GLOBAL_VARNAMES[c->indexing_nvar_checkout[nv]]);
  Drawra(temp,"double",c->nbins,c->nbins,7,1,0,0,0,text,xside,yside,xoffset+2*xside,yoffset-1*yside);
  sprintf(text,"m0/m1*m0/m1_%s",GLOBAL_VARNAMES[c->indexing_nvar_checkout[nv]]);
  Drawra(temp2,"double",c->nbins,c->nbins,7,1,0,0,0,text,xside,yside,xoffset+3*xside,yoffset-1*yside);
  tfree(temp);temp=NULL;
  tfree(temp2);temp2=NULL;
}

GLvoid Drawisi(struct isi *i,double xside,double yside,double xoffset,double yoffset)
{
  double max=0,min=0,mean=0,stdev=0;
  char text[512];
  struct litem *l=NULL;
  double xord=0,yord=0;
  int l2orra=1;
  if (l2orra==0){
    lliststats(i->L2,&max,&min,&mean,&stdev); max = mean+STD_VIEW*stdev; min=mean-STD_VIEW*stdev;
    l=i->L2->first;
    while (l!=NULL){
      if (l->child!=NULL){ 
	xord = (*(double *)l->item - min)/(max-min);
	yord = (*(double *)l->child->item - min)/(max-min);
	glbox(xord*xside+xoffset,yord*yside+yoffset,xside/100,1,1,1);}
      l=l->child;}
    sprintf(text,"isi (%0.3f,%0.3f)",min,max);
    glColor3f(1,1,1); ftexto(xoffset,yoffset+0.1,0,text);}
  else if (l2orra==1){
    sprintf(text,"isi (%0.3f,%0.3f)",i->minisi,i->maxisi);
    Drawra(i->ra,"double",i->nbins,i->nbins,7,1,0,0,0,text,xside,yside,xoffset,yoffset);}
}

void spack3d_temp_vector_draw(void *v1,int draw_flag,double xside,double yside,double zside,double xoffset,double yoffset,double yzoffset,void *void_parameters)
{
  /* assumes that dra[*vector_length] holds a relevant magnitude 
     parameter list:
     (void *) (int *) index
     (void *) (int *) vector_length
   */
  int *index=NULL,*vector_length=NULL;
  void **vra=NULL,**vra_p=NULL;
  double *dra=NULL;
  double xord=0,yord=0,zord=0;
  double rcolor=0,gcolor=0,bcolor=0,kcolor=0;
  if (v1!=NULL && void_parameters!=NULL){
    vra_p=(void **) void_parameters; index = (int *) vra_p[0]; vector_length = (int *) vra_p[1];
    vra = (void **) v1; dra = (double *) vra[*index];
    switch (draw_flag){
    case 1: 
      if (DRAW_FLAG2==0){
	xord = (*vector_length>=0 ? dra[0]: 0); yord = (*vector_length>=1 ? dra[1]: 0); zord = (*vector_length>=2 ? dra[2]: 0);
	colorscale(7,dra[*vector_length],1*STD_VIEW,0,&rcolor,&gcolor,&bcolor);
	glarrow(1,xoffset-xord*xside,yoffset-yord*yside,xoffset+xord*xside,yoffset+yord*yside,0,rcolor,gcolor,bcolor);}
      if (DRAW_FLAG2==1){
	xord = (*vector_length>=0 ? dra[0]: 0); yord = (*vector_length>=1 ? dra[1]: 0); zord = (*vector_length>=2 ? dra[2]: 0);
	colorscale(2,periodize(atan2(yord,xord),-PI/2,PI/2),PI/2,-PI/2,&rcolor,&gcolor,&bcolor);
	colorscale(0,ra_norm(dra,minimum(*vector_length,2)),ra_norm(dra,minimum(*vector_length,3)),0,&kcolor,&kcolor,&kcolor);
	gldot(6,xoffset,yoffset,xside,rcolor*kcolor,gcolor*kcolor,bcolor*kcolor);}
      break;
    default: break;}}
}

GLvoid Drawlattice3d(struct lattice3d *l3d,int draw_flag,double side,double xoffset,double yoffset,void (*void_draw)(void *,int,double,double,double,double,double,double,void *), void *void_parameters)
{
  char text[64];
  int nr1=0,nr2=0,nr3=0;
  struct spack3d *s=NULL;
  int k2x=0,k2y=0;
  double scaleside=side/l3d->input1max/sqrt(2);
  switch (draw_flag){
  case 0:
    for (nr1=0;nr1<l3d->input1max;nr1++){ for (nr2=0;nr2<l3d->input2max_scale;nr2++){ for (nr3=0;nr3<l3d->input3max_scale;nr3++){ s = lattice3dget(l3d,nr1,nr2,nr3); void_draw(s->temp,draw_flag,scaleside,scaleside,scaleside,s->x/l3d->input1max*side+xoffset,s->y/l3d->input2max_scale*side+yoffset,s->z/l3d->input3max_scale*side,void_parameters);}}}
    sprintf(text,"O"); glColor3f(1,1,1); ftexto(0*side+xoffset,0*side+yoffset,0*side,text);
    sprintf(text,"x"); glColor3f(1,0,0); ftexto(1*side+xoffset,0*side+yoffset,0*side,text);
    sprintf(text,"y"); glColor3f(0,1,0); ftexto(0*side+xoffset,1*side+yoffset,0*side,text);
    sprintf(text,"z"); glColor3f(0,0,1); ftexto(0*side+xoffset,0*side+yoffset,1*side,text);
    break;
  case 1: 
    for (nr3=0;nr3<l3d->input3max_scale;nr3++){ 
      k2x = nr3/(int)sqrt(l3d->input3max_scale); k2y = nr3%(int)sqrt(l3d->input3max_scale);
      for (nr1=0;nr1<l3d->input1max;nr1++){ for (nr2=0;nr2<l3d->input2max_scale;nr2++){ s = lattice3dget(l3d,nr1,nr2,nr3); void_draw(s->temp,draw_flag,scaleside,scaleside,scaleside,s->x/l3d->input1max*side+k2x*side+xoffset,s->y/l3d->input2max_scale*side+k2y*side+yoffset,0,void_parameters);}}
      sprintf(text,"k%d",nr3); glColor3f(1,1,1); ftexto(k2x*side+xoffset,k2y*side+yoffset-0.1,0,text);}
    break;
  default: break;}
}

GLvoid Drawlattice3d_edge(struct lattice3d *l3d,double side,double xoffset,double yoffset)
{
  char text[64];
  int nr1=0,nr2=0,nr3=0,tab=0;
  struct spack3d *s0=NULL,*s=NULL,*s2=NULL;
  double max=0,min=0;
  double rcolor=0,gcolor=0,bcolor=0;
  double xpos=0,ypos=0,zpos=0,xpos1=0,ypos1=0,xpos2=0,ypos2=0,d=0;
  int k2x=0,k2y=0;
  int length=GLOBAL_LATTICE3D_LENGTH,depth=GLOBAL_LATTICE3D_DEPTH,zero=0,slength=(int)pow(2,length),pdepth=0;
  int nx=0,ny=0,nz=0;
  double *s2p=NULL,*pra=NULL;
  void **vra=NULL;
  periodify("int",&FIDDLE_ROW,&zero,&(l3d->input1max),&nx);
  periodify("int",&FIDDLE_COL,&zero,&(l3d->input2max_scale),&ny);
  nz = l3d->input3max_scale/2;
  s0=lattice3dget(l3d,nx,ny,nz); vra=(void **)s0->temp;
  sprintf(text,"s2");Drawra((double *)vra[0],"double",slength,slength,7,1,0,0,0,text,side,side,xoffset+2*side,yoffset+3*side);
  pdepth=depth; s2p = binary_projection_s2p(0,length,depth,pdepth);
  pra = ra2ra_matrix_multiply(s2p,(int)pow(length,pdepth),(int)pow(slength,depth),0,(double *)vra[0],(int)pow(slength,depth),1,0);
  sprintf(text,"p2");Drawra(pra,"double",length,length,7,1,0,0,0,text,side,side,xoffset+3*side,yoffset+3*side);
  tfree(pra);pra=NULL;tfree(s2p);s2p=NULL;
  max = l3d->edge_v_mean+MEAN_VIEW*l3d->edge_v_stdev+STD_VIEW*l3d->edge_v_stdev;
  min = l3d->edge_v_mean+MEAN_VIEW*l3d->edge_v_stdev-STD_VIEW*l3d->edge_v_stdev;
  switch (abs(DRAW_FLAG%2)){
  case 0:
    glcube(s0->x*side/l3d->input1max+xoffset,s0->y*side/l3d->input2max_scale+yoffset,s0->z*side/l3d->input3max_scale,0.5*side/l3d->input1max,1,1,1);
    glBegin(GL_LINES);
    for (nr1=0;nr1<l3d->input1max;nr1++){ for (nr2=0;nr2<l3d->input2max_scale;nr2++){ for (nr3=0;nr3<l3d->input3max_scale;nr3++){
      s = lattice3dget(l3d,nr1,nr2,nr3);
      if (s->adjacency_ra_s!=NULL && s->adjacency_ra_v!=NULL){
	for (tab=0;tab<s->L->length;tab++){
	  colorscale(7,s->adjacency_ra_v[tab],max,min,&rcolor,&gcolor,&bcolor); glColor3f(rcolor,gcolor,bcolor);
	  s2 = s->adjacency_ra_s[tab];
	  xpos=(s->x+s2->x)/2.0; ypos=(s->y+s2->y)/2.0; zpos=(s->z+s2->z)/2.0;
	  glVertex3f(s->x*side/l3d->input1max+xoffset,s->y*side/l3d->input2max_scale+yoffset,s->z*side/l3d->input3max_scale);
	  glVertex3f(xpos*side/l3d->input1max+xoffset,ypos*side/l3d->input2max_scale+yoffset,zpos*side/l3d->input3max_scale);}}}}}
    glEnd();
    sprintf(text,"edge [%f,%f]",min,max); glColor3f(1,1,1); ftexto(xoffset,yoffset-0.2,0,text);
    sprintf(text,"O"); glColor3f(1,1,1); ftexto(0*side+xoffset,0*side+yoffset,0*side,text);
    sprintf(text,"x"); glColor3f(1,0,0); ftexto(1*side+xoffset,0*side+yoffset,0*side,text);
    sprintf(text,"y"); glColor3f(0,1,0); ftexto(0*side+xoffset,1*side+yoffset,0*side,text);
    sprintf(text,"z"); glColor3f(0,0,1); ftexto(0*side+xoffset,0*side+yoffset,1*side,text);
    break;
  case 1: 
    k2x = nz/(int)sqrt(l3d->input3max_scale); k2y = nz%(int)sqrt(l3d->input3max_scale);
    glbox(s0->x*side/l3d->input1max+k2x*side+xoffset,s0->y*side/l3d->input2max_scale+k2y*side+yoffset,0.5*side/l3d->input1max,1,1,1);
    for (nr3=0;nr3<l3d->input3max_scale;nr3++){ 
      k2x = nr3/(int)sqrt(l3d->input3max_scale); k2y = nr3%(int)sqrt(l3d->input3max_scale);
      for (nr1=0;nr1<l3d->input1max;nr1++){ for (nr2=0;nr2<l3d->input2max_scale;nr2++){ 
	s = lattice3dget(l3d,nr1,nr2,nr3);
	if (s->adjacency_ra_s!=NULL && s->adjacency_ra_v!=NULL){
	  for (tab=0;tab<s->L->length;tab++){
	    s2 = s->adjacency_ra_s[tab];
	    if (s2->k==s->k){
	      colorscale(7,s->adjacency_ra_v[tab],max,min,&rcolor,&gcolor,&bcolor); glColor3f(rcolor,gcolor,bcolor);
	      xpos=(s->x+s2->x)/2.0; ypos=(s->y+s2->y)/2.0; 
	      xpos1 = s->y-s2->y; ypos1 = s2->x-s->x; d = sqrt(pow(xpos1,2)+pow(ypos1,2)); xpos2 = xpos1/d; ypos2 = ypos1/d;
	      if (l3d->tet_vs_cube==1){
		xpos1 = xpos+xpos2/2/sqrt(3); ypos1 = ypos+ypos2/2/sqrt(3);
		xpos2 = xpos-xpos2/2/sqrt(3); ypos2 = ypos-ypos2/2/sqrt(3);}
	      else /* if (l3d->tet_vs_cube==0) */{
		xpos1 = xpos+xpos2/2; ypos1 = ypos+ypos2/2;
		xpos2 = xpos-xpos2/2; ypos2 = ypos-ypos2/2;}
	      glBegin(GL_POLYGON);
	      glVertex3f(s->x*side/l3d->input1max+k2x*side+xoffset,s->y*side/l3d->input2max_scale+k2y*side+yoffset,0);
	      glVertex3f(xpos1*side/l3d->input1max+k2x*side+xoffset,ypos1*side/l3d->input2max_scale+k2y*side+yoffset,0);
	      glVertex3f(xpos2*side/l3d->input1max+k2x*side+xoffset,ypos2*side/l3d->input2max_scale+k2y*side+yoffset,0);
	      glEnd();}}}}}
      sprintf(text,"k%d",nr3); glColor3f(1,1,1); ftexto(k2x*side+xoffset,k2y*side+yoffset-0.1,0,text);}
    sprintf(text,"edge [%f,%f]",min,max); glColor3f(1,1,1); ftexto(xoffset,yoffset-0.2,0,text);
    break;
  default: break;}
}

GLvoid Drawsnxdata(struct snxdata *s,double side,double xoffset,double yoffset)
{
  char text[1024];
  int nr=0,nl=0;;
  switch (GLOBAL_SNX_VERSION){
  case 1:
    glColor3f(1,1,1);
    sprintf(text,"chain_0 "); for (nr=0;nr<s->length;nr++){ sprintf(text,"%s %f",text,s->chain_0[nr]/s->total_time);}
    ftexto(xoffset,yoffset-0.2,0,text);
    for (nl=0;nl<s->lookback;nl++){
      sprintf(text,"chain_1_%d ",nl+1); 
      for (nr=0;nr<s->length;nr++){ sprintf(text,"%s %f",text,s->chain_1_[nr+nl*s->length]/s->total_time);}
      ftexto(xoffset,yoffset+0.2*nl,0,text);}
    break;
  case 2: 
    Drawra(s->chain_0,"double",1,s->length,7,1,0,0,0,"chain_0",side,side,xoffset,yoffset);
    Drawra(s->chain_1_,"double",s->length,s->lookback,7,1,0,0,0,"chain_1_",side,side,xoffset+side,yoffset);
    glColor3f(1,1,1);
    sprintf(text,"abc %f,bca %f,abc-bca rate %0.16f ",s->abc,s->bca,(s->abc-s->bca)/s->total_time); 
    ftexto(xoffset,yoffset+0.2,0,text);
    break;
  default: break;}    
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
  char text[256];
  int menupos=0;
  double tcolor=0;
  int menupos_cs__start=0;
  int menupos_pf__start=-8;
  int nt1=0,nt2=0,nv=0,tab1=0,tab2=0;
  switch (GLOBAL_NEURON_MODEL){
  case 0: case 1: case 2: case 3: case 4:
    menupos=-10;
    glColor3f(1,1,1);
    sprintf(text,"%s, t=%0.1f,DT=%0.3f",GLOBAL_STRING_2,GLOBAL_time,GLOBAL_DT);
    ftexto(0*side+xoffset,menupos*side+yoffset,0,text);
    for (nt1=0;nt1<GLOBAL_NTYPES;nt1++){ for (nv=0;nv<GLOBAL_INDEXING_sra_LENGTH;nv++){
      nt2=0; tab1=nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES;
      nt2=1; tab2=nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES;
      menupos=menupos_cs__start+nt1+nv*GLOBAL_NTYPES;
      tcolor = .25 + .75*(FIDDLE_PARAMETER==menupos);glColor3f(tcolor,tcolor,tcolor);
      sprintf(text,"CS__[%s->%s,%s]=%0.6f CS__[%s->%s,%s]=%0.6f",GLOBAL_TYPENAMES[nt1],GLOBAL_TYPENAMES[0],GLOBAL_VARNAMES[GLOBAL_INDEXING_CHECKOUT_sra[nv]],CS__[tab1],GLOBAL_TYPENAMES[nt1],GLOBAL_TYPENAMES[1],GLOBAL_VARNAMES[GLOBAL_INDEXING_CHECKOUT_sra[nv]],CS__[tab2]);
      ftexto(0*side+xoffset,menupos*side+yoffset,0,text);}}
    for (nt1=0;nt1<GLOBAL_NTYPES;nt1++){ for (nv=0;nv<GLOBAL_INDEXING_sra_LENGTH;nv++){
      nt2=0; tab1=nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES;
      nt2=1; tab2=nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES;
      menupos=menupos_pf__start+nt1+nv*GLOBAL_NTYPES;
      tcolor = .25 + .75*(FIDDLE_PARAMETER==menupos);glColor3f(tcolor,tcolor,tcolor);
      sprintf(text,"PF__[%s->%s,%s]=%0.4f PF__[%s->%s,%s]=%0.4f",GLOBAL_TYPENAMES[nt1],GLOBAL_TYPENAMES[0],GLOBAL_VARNAMES[GLOBAL_INDEXING_CHECKOUT_sra[nv]],P_FAIL__[tab1],GLOBAL_TYPENAMES[nt1],GLOBAL_TYPENAMES[1],GLOBAL_VARNAMES[GLOBAL_INDEXING_CHECKOUT_sra[nv]],P_FAIL__[tab2]);
      ftexto(0*side+xoffset,menupos*side+yoffset,0,text);}}  
    menupos=9;
    menupos=10;
    tcolor = .25 + .75*(FIDDLE_PARAMETER==menupos);glColor3f(tcolor,tcolor,tcolor);
    sprintf(text,"CS_ORN_SCALE=%0.1f",GLOBAL_CS_ORN_SCALE);
    for (nt1=0;nt1<GLOBAL_NTYPES;nt1++){ sprintf(text,"%s,CS_ORN_[%d]=%0.6f",text,nt1,CS_ORN_[nt1]);}
    ftexto(0*side+xoffset,menupos*side+yoffset,0,text);
    menupos=11;
    tcolor = .25 + .75*(FIDDLE_PARAMETER==menupos);glColor3f(tcolor,tcolor,tcolor);
    sprintf(text,"CONDUCTANCE_SD=%0.2f,CONDUCTANCE_DS=%0.2f",CONDUCTANCE_SD,CONDUCTANCE_DS);
    ftexto(0*side+xoffset,menupos*side+yoffset,0,text);
    menupos=12;
    tcolor = .25 + .75*(FIDDLE_PARAMETER==menupos);glColor3f(tcolor,tcolor,tcolor);
    sprintf(text,"CURRENT_INJECTION_S=%0.2f,CURRENT_INJECTION_D=%0.2f",CURRENT_INJECTION_S,CURRENT_INJECTION_D);
    ftexto(0*side+xoffset,menupos*side+yoffset,0,text); 
    menupos=13;
    tcolor = .25 + .75*(FIDDLE_PARAMETER==menupos);glColor3f(tcolor,tcolor,tcolor);
    sprintf(text,"s_Ca=%0.4f,s_KCa=%0.4f",CONDUCTANCE_[VARNAME_REGISTRY_s_Ca],CONDUCTANCE_[VARNAME_REGISTRY_s_KCa]);
    ftexto(0*side+xoffset,menupos*side+yoffset,0,text); 
    menupos=14;
    tcolor = .25 + .75*(FIDDLE_PARAMETER==menupos);glColor3f(tcolor,tcolor,tcolor);
    sprintf(text,"TAU_Ca=%0.2f",TAU_[VARNAME_REGISTRY_Ca]);
    ftexto(0*side+xoffset,menupos*side+yoffset,0,text); 
    menupos=15;
    tcolor = .25 + .75*(FIDDLE_PARAMETER==menupos);glColor3f(tcolor,tcolor,tcolor);
    sprintf(text,"STD_VIEW=%0.1f,MEAN_VIEW=%0.1f",STD_VIEW,MEAN_VIEW);
    ftexto(0*side+xoffset,menupos*side+yoffset,0,text);
    menupos=18;
    tcolor = .25+.75*(FIDDLE_PARAMETER==menupos);glColor3f(tcolor,tcolor,tcolor);
    sprintf(text,"dtmax=%0.2f",GLOBAL_DTmax);
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
    sprintf(text,"FIDDLE_ROW=%d, FIDDLE_COL=%d",FIDDLE_ROW,FIDDLE_COL);
    ftexto(0*side+xoffset,menupos*side+yoffset,0,text);
    menupos=24;
    tcolor = .25+.75*(FIDDLE_PARAMETER==menupos);glColor3f(tcolor,tcolor,tcolor);
    sprintf(text,"INPUT_CONTRAST=%0.3f,INPUT_PULSE=%0.1f",INPUT_CONTRAST,INPUT_PULSE);
    ftexto(0*side+xoffset,menupos*side+yoffset,0,text);
    break;
  case 5: /* mainak style */
    menupos=-10;
    glColor3f(1,1,1);
    sprintf(text,"t=%0.1f,DT=%0.3f",GLOBAL_time,GLOBAL_DT);
    ftexto(0*side+xoffset,menupos*side+yoffset,0,text);
    for (nt1=0;nt1<GLOBAL_NTYPES;nt1++){ for (nv=0;nv<GLOBAL_INDEXING_sra_LENGTH;nv++){
      nt2=0; tab1=nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES;
      nt2=1; tab2=nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES;
      menupos=menupos_cs__start+nt1+nv*GLOBAL_NTYPES;
      tcolor = .25 + .75*(FIDDLE_PARAMETER==menupos);glColor3f(tcolor,tcolor,tcolor);
      sprintf(text,"CS__[%s->%s,%s]=%0.6f CS__[%s->%s,%s]=%0.6f",GLOBAL_TYPENAMES[nt1],GLOBAL_TYPENAMES[0],GLOBAL_VARNAMES[GLOBAL_INDEXING_CHECKOUT_sra[nv]],CS__[tab1],GLOBAL_TYPENAMES[nt1],GLOBAL_TYPENAMES[1],GLOBAL_VARNAMES[GLOBAL_INDEXING_CHECKOUT_sra[nv]],CS__[tab2]);
      ftexto(0*side+xoffset,menupos*side+yoffset,0,text);}}
    for (nt1=0;nt1<GLOBAL_NTYPES;nt1++){ for (nv=0;nv<GLOBAL_INDEXING_sra_LENGTH;nv++){
      nt2=0; tab1=nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES;
      nt2=1; tab2=nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES;
      menupos=menupos_pf__start+nt1+nv*GLOBAL_NTYPES;
      tcolor = .25 + .75*(FIDDLE_PARAMETER==menupos);glColor3f(tcolor,tcolor,tcolor);
      sprintf(text,"PF__[%s->%s,%s]=%0.4f PF__[%s->%s,%s]=%0.4f",GLOBAL_TYPENAMES[nt1],GLOBAL_TYPENAMES[0],GLOBAL_VARNAMES[GLOBAL_INDEXING_CHECKOUT_sra[nv]],P_FAIL__[tab1],GLOBAL_TYPENAMES[nt1],GLOBAL_TYPENAMES[1],GLOBAL_VARNAMES[GLOBAL_INDEXING_CHECKOUT_sra[nv]],P_FAIL__[tab2]);
      ftexto(0*side+xoffset,menupos*side+yoffset,0,text);}}  
    menupos=9;
    menupos=10;
    tcolor = .25 + .75*(FIDDLE_PARAMETER==menupos);glColor3f(tcolor,tcolor,tcolor);
    sprintf(text,"CS_ORN_SCALE=%0.1f",GLOBAL_CS_ORN_SCALE);
    for (nt1=0;nt1<GLOBAL_NTYPES;nt1++){ sprintf(text,"%s,CS_ORN_[%d]=%0.6f",text,nt1,CS_ORN_[nt1]);}
    ftexto(0*side+xoffset,menupos*side+yoffset,0,text);
    menupos=11;
    tcolor = .25 + .75*(FIDDLE_PARAMETER==menupos);glColor3f(tcolor,tcolor,tcolor);
    sprintf(text," ");
    for (nt1=0;nt1<GLOBAL_NTYPES;nt1++){ sprintf(text,"%s,CS_ORN_mainak_stim_[%d]=%0.6f",text,nt1,CS_ORN_mainak_stim_[nt1]);}
    ftexto(0*side+xoffset,menupos*side+yoffset,0,text);
    menupos=12;
    tcolor = .25 + .75*(FIDDLE_PARAMETER==menupos);glColor3f(tcolor,tcolor,tcolor);
    sprintf(text,"CURRENT_INJECTION_S=%0.2f,CURRENT_INJECTION_D=%0.2f",CURRENT_INJECTION_S,CURRENT_INJECTION_D);
    ftexto(0*side+xoffset,menupos*side+yoffset,0,text); 
    menupos=13;
    menupos=14;
    menupos=15;
    tcolor = .25 + .75*(FIDDLE_PARAMETER==menupos);glColor3f(tcolor,tcolor,tcolor);
    sprintf(text,"STD_VIEW=%0.1f,MEAN_VIEW=%0.1f",STD_VIEW,MEAN_VIEW);
    ftexto(0*side+xoffset,menupos*side+yoffset,0,text);
    menupos=18;
    tcolor = .25+.75*(FIDDLE_PARAMETER==menupos);glColor3f(tcolor,tcolor,tcolor);
    sprintf(text,"dtmax=%0.2f",GLOBAL_DTmax);
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
    sprintf(text,"FIDDLE_ROW=%d, FIDDLE_COL=%d",FIDDLE_ROW,FIDDLE_COL);
    ftexto(0*side+xoffset,menupos*side+yoffset,0,text);
    menupos=24;
    tcolor = .25+.75*(FIDDLE_PARAMETER==menupos);glColor3f(tcolor,tcolor,tcolor);
    sprintf(text,"INPUT_CONTRAST=%0.3f,INPUT_PULSE=%0.1f",INPUT_CONTRAST,INPUT_PULSE);
    ftexto(0*side+xoffset,menupos*side+yoffset,0,text);
    break;
  case 6: /* wilson style */
    menupos=-10;
    glColor3f(1,1,1);
    sprintf(text,"t=%0.1f,DT=%0.3f",GLOBAL_time,GLOBAL_DT);
    ftexto(0*side+xoffset,menupos*side+yoffset,0,text);
    menupos=0;
    tcolor = .25 + .75*(FIDDLE_PARAMETER==menupos);glColor3f(tcolor,tcolor,tcolor);
    nv=0;nt1=0,nt2=1; tab2 = nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES; sprintf(text,"CS__[%d], %s->%s,%s %f",tab2,GLOBAL_TYPENAMES[nt1],GLOBAL_TYPENAMES[nt2],GLOBAL_VARNAMES[GLOBAL_INDEXING_CHECKOUT_sra[nv]],CS__[tab2]);
    nv=0;nt1=1,nt2=1; tab2 = nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES; sprintf(text,"%s, CS__[%d], %s->%s,%s %f",text,tab2,GLOBAL_TYPENAMES[nt1],GLOBAL_TYPENAMES[nt2],GLOBAL_VARNAMES[GLOBAL_INDEXING_CHECKOUT_sra[nv]],CS__[tab2]);
    ftexto(0*side+xoffset,menupos*side+yoffset,0,text);
    menupos=1;
    tcolor = .25 + .75*(FIDDLE_PARAMETER==menupos);glColor3f(tcolor,tcolor,tcolor);
    nv=0;nt1=0,nt2=2; tab2 = nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES; sprintf(text,"CS__[%d], %s->%s,%s %f",tab2,GLOBAL_TYPENAMES[nt1],GLOBAL_TYPENAMES[nt2],GLOBAL_VARNAMES[GLOBAL_INDEXING_CHECKOUT_sra[nv]],CS__[tab2]);
    nv=0;nt1=1,nt2=2; tab2 = nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES; sprintf(text,"%s, CS__[%d], %s->%s,%s %f",text,tab2,GLOBAL_TYPENAMES[nt1],GLOBAL_TYPENAMES[nt2],GLOBAL_VARNAMES[GLOBAL_INDEXING_CHECKOUT_sra[nv]],CS__[tab2]);
    ftexto(0*side+xoffset,menupos*side+yoffset,0,text);
    menupos=2;
    tcolor = .25 + .75*(FIDDLE_PARAMETER==menupos);glColor3f(tcolor,tcolor,tcolor);
    nv=0;nt1=0,nt2=3; tab2 = nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES; sprintf(text,"CS__[%d], %s->%s,%s %f",tab2,GLOBAL_TYPENAMES[nt1],GLOBAL_TYPENAMES[nt2],GLOBAL_VARNAMES[GLOBAL_INDEXING_CHECKOUT_sra[nv]],CS__[tab2]);
    nv=0;nt1=1,nt2=3; tab2 = nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES; sprintf(text,"%s, CS__[%d], %s->%s,%s %f",text,tab2,GLOBAL_TYPENAMES[nt1],GLOBAL_TYPENAMES[nt2],GLOBAL_VARNAMES[GLOBAL_INDEXING_CHECKOUT_sra[nv]],CS__[tab2]);
    ftexto(0*side+xoffset,menupos*side+yoffset,0,text);
    menupos=3;
    tcolor = .25 + .75*(FIDDLE_PARAMETER==menupos);glColor3f(tcolor,tcolor,tcolor);
    nv=1;nt1=3,nt2=0; tab2 = nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES; sprintf(text,"CS__[%d], %s->%s,%s %f",tab2,GLOBAL_TYPENAMES[nt1],GLOBAL_TYPENAMES[nt2],GLOBAL_VARNAMES[GLOBAL_INDEXING_CHECKOUT_sra[nv]],CS__[tab2]);
    nv=2;nt1=3,nt2=0; tab2 = nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES; sprintf(text,"%s, CS__[%d], %s->%s,%s %f",text,tab2,GLOBAL_TYPENAMES[nt1],GLOBAL_TYPENAMES[nt2],GLOBAL_VARNAMES[GLOBAL_INDEXING_CHECKOUT_sra[nv]],CS__[tab2]);
    ftexto(0*side+xoffset,menupos*side+yoffset,0,text);
    menupos=4;
    tcolor = .25 + .75*(FIDDLE_PARAMETER==menupos);glColor3f(tcolor,tcolor,tcolor);
    nv=1;nt1=3,nt2=1; tab2 = nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES; sprintf(text,"CS__[%d], %s->%s,%s %f",tab2,GLOBAL_TYPENAMES[nt1],GLOBAL_TYPENAMES[nt2],GLOBAL_VARNAMES[GLOBAL_INDEXING_CHECKOUT_sra[nv]],CS__[tab2]);
    nv=2;nt1=3,nt2=1; tab2 = nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES; sprintf(text,"%s, CS__[%d], %s->%s,%s %f",text,tab2,GLOBAL_TYPENAMES[nt1],GLOBAL_TYPENAMES[nt2],GLOBAL_VARNAMES[GLOBAL_INDEXING_CHECKOUT_sra[nv]],CS__[tab2]);
    ftexto(0*side+xoffset,menupos*side+yoffset,0,text);
    menupos=5;
    tcolor = .25 + .75*(FIDDLE_PARAMETER==menupos);glColor3f(tcolor,tcolor,tcolor);
    nv=1;nt1=3,nt2=2; tab2 = nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES; sprintf(text,"CS__[%d], %s->%s,%s %f",tab2,GLOBAL_TYPENAMES[nt1],GLOBAL_TYPENAMES[nt2],GLOBAL_VARNAMES[GLOBAL_INDEXING_CHECKOUT_sra[nv]],CS__[tab2]);
    nv=2;nt1=3,nt2=2; tab2 = nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES; sprintf(text,"%s, CS__[%d], %s->%s,%s %f",text,tab2,GLOBAL_TYPENAMES[nt1],GLOBAL_TYPENAMES[nt2],GLOBAL_VARNAMES[GLOBAL_INDEXING_CHECKOUT_sra[nv]],CS__[tab2]);
    ftexto(0*side+xoffset,menupos*side+yoffset,0,text);
    menupos=6;
    tcolor = .25 + .75*(FIDDLE_PARAMETER==menupos);glColor3f(tcolor,tcolor,tcolor);
    nv=1;nt1=3,nt2=3; tab2 = nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES; sprintf(text,"CS__[%d], %s->%s,%s %f",tab2,GLOBAL_TYPENAMES[nt1],GLOBAL_TYPENAMES[nt2],GLOBAL_VARNAMES[GLOBAL_INDEXING_CHECKOUT_sra[nv]],CS__[tab2]);
    nv=2;nt1=3,nt2=3; tab2 = nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES; sprintf(text,"%s, CS__[%d], %s->%s,%s %f",text,tab2,GLOBAL_TYPENAMES[nt1],GLOBAL_TYPENAMES[nt2],GLOBAL_VARNAMES[GLOBAL_INDEXING_CHECKOUT_sra[nv]],CS__[tab2]);
    ftexto(0*side+xoffset,menupos*side+yoffset,0,text);
    menupos=7;
    tcolor = .25 + .75*(FIDDLE_PARAMETER==menupos);glColor3f(tcolor,tcolor,tcolor);
    nv=0;nt1=2,nt2=1; tab2 = nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES; sprintf(text,"CS__[%d], %s->%s,%s %f",tab2,GLOBAL_TYPENAMES[nt1],GLOBAL_TYPENAMES[nt2],GLOBAL_VARNAMES[GLOBAL_INDEXING_CHECKOUT_sra[nv]],CS__[tab2]);
    ftexto(0*side+xoffset,menupos*side+yoffset,0,text);
    menupos=8;
    tcolor = .25 + .75*(FIDDLE_PARAMETER==menupos);glColor3f(tcolor,tcolor,tcolor);
    nv=0;nt1=2,nt2=2; tab2 = nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES; sprintf(text,"CS__[%d], %s->%s,%s %f",tab2,GLOBAL_TYPENAMES[nt1],GLOBAL_TYPENAMES[nt2],GLOBAL_VARNAMES[GLOBAL_INDEXING_CHECKOUT_sra[nv]],CS__[tab2]);
    ftexto(0*side+xoffset,menupos*side+yoffset,0,text);
    menupos=9;
    tcolor = .25 + .75*(FIDDLE_PARAMETER==menupos);glColor3f(tcolor,tcolor,tcolor);
    nv=0;nt1=2,nt2=3; tab2 = nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES; sprintf(text,"CS__[%d], %s->%s,%s %f",tab2,GLOBAL_TYPENAMES[nt1],GLOBAL_TYPENAMES[nt2],GLOBAL_VARNAMES[GLOBAL_INDEXING_CHECKOUT_sra[nv]],CS__[tab2]);
    ftexto(0*side+xoffset,menupos*side+yoffset,0,text);
    menupos=10;
    menupos=11;
    menupos=12;
    tcolor = .25 + .75*(FIDDLE_PARAMETER==menupos);glColor3f(tcolor,tcolor,tcolor);
    sprintf(text,"CURRENT_INJECTION_S=%0.2f,CURRENT_INJECTION_D=%0.2f",CURRENT_INJECTION_S,CURRENT_INJECTION_D);
    ftexto(0*side+xoffset,menupos*side+yoffset,0,text); 
    menupos=13;
    menupos=14;
    menupos=15;
    tcolor = .25 + .75*(FIDDLE_PARAMETER==menupos);glColor3f(tcolor,tcolor,tcolor);
    sprintf(text,"STD_VIEW=%0.1f,MEAN_VIEW=%0.1f",STD_VIEW,MEAN_VIEW);
    ftexto(0*side+xoffset,menupos*side+yoffset,0,text);
    menupos=18;
    tcolor = .25+.75*(FIDDLE_PARAMETER==menupos);glColor3f(tcolor,tcolor,tcolor);
    sprintf(text,"dtmax=%0.2f",GLOBAL_DTmax);
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
    sprintf(text,"FIDDLE_ROW=%d, FIDDLE_COL=%d",FIDDLE_ROW,FIDDLE_COL);
    ftexto(0*side+xoffset,menupos*side+yoffset,0,text);
    menupos=24;
    tcolor = .25+.75*(FIDDLE_PARAMETER==menupos);glColor3f(tcolor,tcolor,tcolor);
    sprintf(text,"INPUT_CONTRAST=%0.3f,INPUT_PULSE=%0.1f",INPUT_CONTRAST,INPUT_PULSE);
    ftexto(0*side+xoffset,menupos*side+yoffset,0,text);
    break;
  case 7: /* snx style */
    menupos=-10;
    glColor3f(1,1,1);
    sprintf(text,"t=%0.1f,DT=%0.3f",GLOBAL_time,GLOBAL_DT);
    ftexto(0*side+xoffset,menupos*side+yoffset,0,text);
    menupos=0;menupos=1;menupos=2;
    menupos=3;menupos=4;menupos=5;
    menupos=6;menupos=7;menupos=8;
    menupos=9;menupos=10;menupos=11;
    menupos=12;
    menupos=13;
    tcolor = .25 + .75*(FIDDLE_PARAMETER==menupos);glColor3f(tcolor,tcolor,tcolor);
    nv=0;nt1=0,nt2=0; tab2 = nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES; sprintf(text,"CS__[%d], %s->%s,%s %f",tab2,GLOBAL_TYPENAMES[nt1],GLOBAL_TYPENAMES[nt2],GLOBAL_VARNAMES[GLOBAL_INDEXING_CHECKOUT_sra[nv]],CS__[tab2]);
    ftexto(0*side+xoffset,menupos*side+yoffset,0,text);
    menupos=14;
    tcolor = .25 + .75*(FIDDLE_PARAMETER==menupos);glColor3f(tcolor,tcolor,tcolor);
    sprintf(text,"CS_ORN_[%s]=%f",GLOBAL_TYPENAMES[TYPENAME_REGISTRY_snx_PN],CS_ORN_[TYPENAME_REGISTRY_snx_PN]);
    ftexto(0*side+xoffset,menupos*side+yoffset,0,text);
    menupos=15;
    tcolor = .25 + .75*(FIDDLE_PARAMETER==menupos);glColor3f(tcolor,tcolor,tcolor);
    sprintf(text,"STD_VIEW=%0.1f,MEAN_VIEW=%0.1f",STD_VIEW,MEAN_VIEW);
    ftexto(0*side+xoffset,menupos*side+yoffset,0,text);
    menupos=18;
    tcolor = .25+.75*(FIDDLE_PARAMETER==menupos);glColor3f(tcolor,tcolor,tcolor);
    sprintf(text,"dtmax=%0.2f",GLOBAL_DTmax);
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
    sprintf(text,"FIDDLE_ROW=%d, FIDDLE_COL=%d",FIDDLE_ROW,FIDDLE_COL);
    ftexto(0*side+xoffset,menupos*side+yoffset,0,text);
    menupos=24;
    tcolor = .25+.75*(FIDDLE_PARAMETER==menupos);glColor3f(tcolor,tcolor,tcolor);
    sprintf(text,"INPUT_CONTRAST=%0.3f,INPUT_PULSE=%0.1f",INPUT_CONTRAST,INPUT_PULSE);
    ftexto(0*side+xoffset,menupos*side+yoffset,0,text);
    break;
  default: break;}

}

void keyPressed(unsigned char key, int x, int y) 
{
  int menupos_cs__start=0;
  int menupos_pf__start=-8;
  int nt1=0,nt2=0,nv=0,tab1=0,tab2=0;
  /* The function called whenever a normal key is pressed. */
  usleep(100); /* still not sure why this is called */
  switch (GLOBAL_NEURON_MODEL){
  case 0: case 1: case 2: case 3: case 4: 
    switch (key) {    
    case ESCAPE_KEY: /* stop everything */
      glutDestroyWindow(GLUTWINDOWNUMBER); /* shut down our window */
      printf("(x,y,z) = (%f,%f,%f)\n",xdepth,ydepth,zdepth);
      exit(EXIT_SUCCESS); /* exit the program...normal termination. */
      break; /* just in case */
    case SPACE_KEY: /* stop/restart temporarily */
      RYangle = 0; RXangle = 0; RYdangle = 0; RXdangle = 0;
      STEPS_PER_DRAW=1024*!STEPS_PER_DRAW;if (!STEPS_PER_DRAW){ printf("PROCESS HALTED\n");}
      break;
    case 'w': FIDDLE_PARAMETER++; if (FIDDLE_PARAMETER > 25){ FIDDLE_PARAMETER = -10;} break;
    case 's': FIDDLE_PARAMETER--; if (FIDDLE_PARAMETER < -10){ FIDDLE_PARAMETER = 25;} break;
    case 'd': 
      if (FIDDLE_PARAMETER>=menupos_cs__start && FIDDLE_PARAMETER<menupos_cs__start+GLOBAL_NTYPES*GLOBAL_INDEXING_sra_LENGTH){
	nt1=(FIDDLE_PARAMETER-menupos_cs__start)%GLOBAL_NTYPES;
	nt2=0;
	nv=(FIDDLE_PARAMETER-menupos_cs__start)/GLOBAL_NTYPES;
	tab1=nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES;
	CS__[tab1] *=pow(2,1.0/8.0);}
      if (FIDDLE_PARAMETER>=menupos_pf__start && FIDDLE_PARAMETER<menupos_pf__start+GLOBAL_NTYPES*GLOBAL_INDEXING_sra_LENGTH){
	nt1=(FIDDLE_PARAMETER-menupos_pf__start)%GLOBAL_NTYPES;
	nt2=0;
	nv=(FIDDLE_PARAMETER-menupos_pf__start)/GLOBAL_NTYPES;
	tab1=nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES;
	P_FAIL__[tab1] *=pow(2,1.0/8.0);}
      switch (FIDDLE_PARAMETER){
      case -10: GRAYSCALE = !GRAYSCALE; break;
      case 9: break;
      case 10: CS_ORN_[0] *= pow(2,1.0/8.0); break;
      case 11: CONDUCTANCE_SD *= pow(2,1.0/8.0); break;
      case 12: CURRENT_INJECTION_S *= pow(2,1.0/8.0); break;
      case 13: CONDUCTANCE_[VARNAME_REGISTRY_s_Ca] *= pow(2,1.0/8.0); break;
      case 14: TAU_[VARNAME_REGISTRY_Ca] *= pow(2,1.0/8.0); break;
      case 15: STD_VIEW *= pow(2,1.0/8.0); break;
      case 18: GLOBAL_DTmax *= 2; break;
      case 20: STEPS_PER_DRAW *= 2; if (STEPS_PER_DRAW==0){ STEPS_PER_DRAW=1;} break;
      case 21: DRAW_FLAG++; break;
      case 22: GLOBAL_SPACE_SMOOTHER++; break;
      case 23: FIDDLE_ROW += 1; break;
      case 24: INPUT_CONTRAST+=0.05; break;
      case 25: break;
      default: break;}
      break;
    case 'D':
      if (FIDDLE_PARAMETER>=menupos_cs__start && FIDDLE_PARAMETER<menupos_cs__start+GLOBAL_NTYPES*GLOBAL_INDEXING_sra_LENGTH){
	nt1=(FIDDLE_PARAMETER-menupos_cs__start)%GLOBAL_NTYPES;
	nt2=1;
	nv=(FIDDLE_PARAMETER-menupos_cs__start)/GLOBAL_NTYPES;
	tab1=nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES;
	CS__[tab1] *=pow(2,1.0/8.0);}
      if (FIDDLE_PARAMETER>=menupos_pf__start && FIDDLE_PARAMETER<menupos_pf__start+GLOBAL_NTYPES*GLOBAL_INDEXING_sra_LENGTH){
	nt1=(FIDDLE_PARAMETER-menupos_pf__start)%GLOBAL_NTYPES;
	nt2=1;
	nv=(FIDDLE_PARAMETER-menupos_pf__start)/GLOBAL_NTYPES;
	tab1=nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES;
	P_FAIL__[tab1] *=pow(2,1.0/8.0);}
      switch (FIDDLE_PARAMETER){
      case 9: break;
      case 10: CS_ORN_[1] *= pow(2,1.0/8.0); break;
      case 11: CONDUCTANCE_DS *= pow(2,1.0/8.0); break;
      case 12: CURRENT_INJECTION_D *= pow(2,1.0/8.0); break;
      case 13: CONDUCTANCE_[VARNAME_REGISTRY_s_KCa] *= pow(2,1.0/8.0); break;
      case 15: MEAN_VIEW += 0.1; break;
      case 21: DRAW_FLAG2++; break;
      case 23: FIDDLE_COL += 1; break;
      case 24: INPUT_PULSE=1; break;
      default: break;}
      break;
    case 'a':
      if (FIDDLE_PARAMETER>=menupos_cs__start && FIDDLE_PARAMETER<menupos_cs__start+GLOBAL_NTYPES*GLOBAL_INDEXING_sra_LENGTH){
	nt1=(FIDDLE_PARAMETER-menupos_cs__start)%GLOBAL_NTYPES;
	nt2=0;
	nv=(FIDDLE_PARAMETER-menupos_cs__start)/GLOBAL_NTYPES;
	tab1=nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES;
	CS__[tab1] /=pow(2,1.0/8.0);}
      if (FIDDLE_PARAMETER>=menupos_pf__start && FIDDLE_PARAMETER<menupos_pf__start+GLOBAL_NTYPES*GLOBAL_INDEXING_sra_LENGTH){
	nt1=(FIDDLE_PARAMETER-menupos_pf__start)%GLOBAL_NTYPES;
	nt2=0;
	nv=(FIDDLE_PARAMETER-menupos_pf__start)/GLOBAL_NTYPES;
	tab1=nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES;
	P_FAIL__[tab1] /=pow(2,1.0/8.0);}
      switch (FIDDLE_PARAMETER){
      case -10: GRAYSCALE = !GRAYSCALE; break;
      case -4: break;
      case 9: break;
      case 10: CS_ORN_[0] /= pow(2,1.0/8.0); break;
      case 11: CONDUCTANCE_SD /= pow(2,1.0/8.0); break;
      case 12: CURRENT_INJECTION_S /= pow(2,1.0/8.0); break;
      case 13: CONDUCTANCE_[VARNAME_REGISTRY_s_Ca] /= pow(2,1.0/8.0); break;
      case 14: TAU_[VARNAME_REGISTRY_Ca] /= pow(2,1.0/8.0); break;
      case 15: STD_VIEW /= pow(2,1.0/8.0); break;
      case 18: GLOBAL_DTmax /= 2; break;
      case 20: STEPS_PER_DRAW /= 2; if (STEPS_PER_DRAW < 1){ STEPS_PER_DRAW=0;} break;
      case 21: DRAW_FLAG--; break;
      case 22: GLOBAL_SPACE_SMOOTHER--; if (GLOBAL_SPACE_SMOOTHER<0){ GLOBAL_SPACE_SMOOTHER=0;} break;
      case 23: FIDDLE_ROW -= 1; break;
      case 24: INPUT_CONTRAST-=0.05; break;
      case 25: break;
      default: break;}
      break;
    case 'A':
      if (FIDDLE_PARAMETER>=menupos_cs__start && FIDDLE_PARAMETER<menupos_cs__start+GLOBAL_NTYPES*GLOBAL_INDEXING_sra_LENGTH){
	nt1=(FIDDLE_PARAMETER-menupos_cs__start)%GLOBAL_NTYPES;
	nt2=1;
	nv=(FIDDLE_PARAMETER-menupos_cs__start)/GLOBAL_NTYPES;
	tab1=nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES;
	CS__[tab1] /=pow(2,1.0/8.0);}
      if (FIDDLE_PARAMETER>=menupos_pf__start && FIDDLE_PARAMETER<menupos_pf__start+GLOBAL_NTYPES*GLOBAL_INDEXING_sra_LENGTH){
	nt1=(FIDDLE_PARAMETER-menupos_pf__start)%GLOBAL_NTYPES;
	nt2=1;
	nv=(FIDDLE_PARAMETER-menupos_pf__start)/GLOBAL_NTYPES;
	tab1=nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES;
	P_FAIL__[tab1] /=pow(2,1.0/8.0);}
      switch (FIDDLE_PARAMETER){
      case 9: break;
      case 10: CS_ORN_[1] /= pow(2,1.0/8.0); break;
      case 11: CONDUCTANCE_DS /= pow(2,1.0/8.0); break;
      case 12: CURRENT_INJECTION_D /= pow(2,1.0/8.0); break;
      case 13: CONDUCTANCE_[VARNAME_REGISTRY_s_KCa] /= pow(2,1.0/8.0); break;
      case 15: MEAN_VIEW -= 0.1; break;
      case 21: DRAW_FLAG2--; break;
      case 23: FIDDLE_COL -= 1; break;
      case 24: INPUT_PULSE=0; break;
      default: break;}
      break;
    default: break;}
    break;
  case 5: /* mainak style */
    switch (key) {    
    case ESCAPE_KEY: /* stop everything */
      glutDestroyWindow(GLUTWINDOWNUMBER); /* shut down our window */
      printf("(x,y,z) = (%f,%f,%f)\n",xdepth,ydepth,zdepth);
      exit(EXIT_SUCCESS); /* exit the program...normal termination. */
      break; /* just in case */
    case SPACE_KEY: /* stop/restart temporarily */
      RYangle = 0; RXangle = 0; RYdangle = 0; RXdangle = 0;
      STEPS_PER_DRAW=1024*!STEPS_PER_DRAW;if (!STEPS_PER_DRAW){ printf("PROCESS HALTED\n");}
      break;
    case 'w': FIDDLE_PARAMETER++; if (FIDDLE_PARAMETER > 25){ FIDDLE_PARAMETER = -10;} break;
    case 's': FIDDLE_PARAMETER--; if (FIDDLE_PARAMETER < -10){ FIDDLE_PARAMETER = 25;} break;
    case 'd': 
      if (FIDDLE_PARAMETER>=menupos_cs__start && FIDDLE_PARAMETER<menupos_cs__start+GLOBAL_NTYPES*GLOBAL_INDEXING_sra_LENGTH){
	nt1=(FIDDLE_PARAMETER-menupos_cs__start)%GLOBAL_NTYPES;
	nt2=0;
	nv=(FIDDLE_PARAMETER-menupos_cs__start)/GLOBAL_NTYPES;
	tab1=nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES;
	CS__[tab1] *=pow(2,1.0/8.0);}
      if (FIDDLE_PARAMETER>=menupos_pf__start && FIDDLE_PARAMETER<menupos_pf__start+GLOBAL_NTYPES*GLOBAL_INDEXING_sra_LENGTH){
	nt1=(FIDDLE_PARAMETER-menupos_pf__start)%GLOBAL_NTYPES;
	nt2=0;
	nv=(FIDDLE_PARAMETER-menupos_pf__start)/GLOBAL_NTYPES;
	tab1=nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES;
	P_FAIL__[tab1] *=pow(2,1.0/8.0);}
      switch (FIDDLE_PARAMETER){
      case -10: GRAYSCALE = !GRAYSCALE; break;
      case 9: break;
      case 10: CS_ORN_[0] *= pow(2,1.0/8.0); break;
      case 11: CS_ORN_mainak_stim_[0] *= pow(2,1.0/8.0); break;
      case 12: CURRENT_INJECTION_S *= pow(2,1.0/8.0); break;
      case 13: break;
      case 14: break;
      case 15: STD_VIEW *= pow(2,1.0/8.0); break;
      case 18: GLOBAL_DTmax *= 2; break;
      case 20: STEPS_PER_DRAW *= 2; if (STEPS_PER_DRAW==0){ STEPS_PER_DRAW=1;} break;
      case 21: DRAW_FLAG++; break;
      case 22: GLOBAL_SPACE_SMOOTHER++; break;
      case 23: FIDDLE_ROW += 1; break;
      case 24: INPUT_CONTRAST+=0.05; break;
      case 25: break;
      default: break;}
      break;
    case 'D':
      if (FIDDLE_PARAMETER>=menupos_cs__start && FIDDLE_PARAMETER<menupos_cs__start+GLOBAL_NTYPES*GLOBAL_INDEXING_sra_LENGTH){
	nt1=(FIDDLE_PARAMETER-menupos_cs__start)%GLOBAL_NTYPES;
	nt2=1;
	nv=(FIDDLE_PARAMETER-menupos_cs__start)/GLOBAL_NTYPES;
	tab1=nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES;
	CS__[tab1] *=pow(2,1.0/8.0);}
      if (FIDDLE_PARAMETER>=menupos_pf__start && FIDDLE_PARAMETER<menupos_pf__start+GLOBAL_NTYPES*GLOBAL_INDEXING_sra_LENGTH){
	nt1=(FIDDLE_PARAMETER-menupos_pf__start)%GLOBAL_NTYPES;
	nt2=1;
	nv=(FIDDLE_PARAMETER-menupos_pf__start)/GLOBAL_NTYPES;
	tab1=nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES;
	P_FAIL__[tab1] *=pow(2,1.0/8.0);}
      switch (FIDDLE_PARAMETER){
      case 9: break;
      case 10: CS_ORN_[1] *= pow(2,1.0/8.0); break;
      case 11: CS_ORN_mainak_stim_[1] *= pow(2,1.0/8.0); break;
      case 12: CURRENT_INJECTION_D *= pow(2,1.0/8.0); break;
      case 13: break;
      case 14: break;
      case 15: MEAN_VIEW += 0.1; break;
      case 21: DRAW_FLAG2++; break;
      case 23: FIDDLE_COL += 1; break;
      case 24: INPUT_PULSE=1; break;
      default: break;}
      break;
    case 'a':
      if (FIDDLE_PARAMETER>=menupos_cs__start && FIDDLE_PARAMETER<menupos_cs__start+GLOBAL_NTYPES*GLOBAL_INDEXING_sra_LENGTH){
	nt1=(FIDDLE_PARAMETER-menupos_cs__start)%GLOBAL_NTYPES;
	nt2=0;
	nv=(FIDDLE_PARAMETER-menupos_cs__start)/GLOBAL_NTYPES;
	tab1=nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES;
	CS__[tab1] /=pow(2,1.0/8.0);}
      if (FIDDLE_PARAMETER>=menupos_pf__start && FIDDLE_PARAMETER<menupos_pf__start+GLOBAL_NTYPES*GLOBAL_INDEXING_sra_LENGTH){
	nt1=(FIDDLE_PARAMETER-menupos_pf__start)%GLOBAL_NTYPES;
	nt2=0;
	nv=(FIDDLE_PARAMETER-menupos_pf__start)/GLOBAL_NTYPES;
	tab1=nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES;
	P_FAIL__[tab1] /=pow(2,1.0/8.0);}
      switch (FIDDLE_PARAMETER){
      case -10: GRAYSCALE = !GRAYSCALE; break;
      case -4: break;
      case 9: break;
      case 10: CS_ORN_[0] /= pow(2,1.0/8.0); break;
      case 11: CS_ORN_mainak_stim_[0] /= pow(2,1.0/8.0); break;
      case 12: CURRENT_INJECTION_S /= pow(2,1.0/8.0); break;
      case 13: break;
      case 14: break;
      case 15: STD_VIEW /= pow(2,1.0/8.0); break;
      case 18: GLOBAL_DTmax /= 2; break;
      case 20: STEPS_PER_DRAW /= 2; if (STEPS_PER_DRAW < 1){ STEPS_PER_DRAW=0;} break;
      case 21: DRAW_FLAG--; break;
      case 22: GLOBAL_SPACE_SMOOTHER--; if (GLOBAL_SPACE_SMOOTHER<0){ GLOBAL_SPACE_SMOOTHER=0;} break;
      case 23: FIDDLE_ROW -= 1; break;
      case 24: INPUT_CONTRAST-=0.05; break;
      case 25: break;
      default: break;}
      break;
    case 'A':
      if (FIDDLE_PARAMETER>=menupos_cs__start && FIDDLE_PARAMETER<menupos_cs__start+GLOBAL_NTYPES*GLOBAL_INDEXING_sra_LENGTH){
	nt1=(FIDDLE_PARAMETER-menupos_cs__start)%GLOBAL_NTYPES;
	nt2=1;
	nv=(FIDDLE_PARAMETER-menupos_cs__start)/GLOBAL_NTYPES;
	tab1=nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES;
	CS__[tab1] /=pow(2,1.0/8.0);}
      if (FIDDLE_PARAMETER>=menupos_pf__start && FIDDLE_PARAMETER<menupos_pf__start+GLOBAL_NTYPES*GLOBAL_INDEXING_sra_LENGTH){
	nt1=(FIDDLE_PARAMETER-menupos_pf__start)%GLOBAL_NTYPES;
	nt2=1;
	nv=(FIDDLE_PARAMETER-menupos_pf__start)/GLOBAL_NTYPES;
	tab1=nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES;
	P_FAIL__[tab1] /=pow(2,1.0/8.0);}
      switch (FIDDLE_PARAMETER){
      case 9: break;
      case 10: CS_ORN_[1] /= pow(2,1.0/8.0); break;
      case 11: CS_ORN_mainak_stim_[1] /= pow(2,1.0/8.0); break;
      case 12: CURRENT_INJECTION_D /= pow(2,1.0/8.0); break;
      case 13: break;
      case 14: break;
      case 15: MEAN_VIEW -= 0.1; break;
      case 21: DRAW_FLAG2--; break;
      case 23: FIDDLE_COL -= 1; break;
      case 24: INPUT_PULSE=0; break;
      default: break;}
      break;
    default: break;}
    break;
  case 6: /* wilson style */
    switch (key) {    
    case ESCAPE_KEY: /* stop everything */
      glutDestroyWindow(GLUTWINDOWNUMBER); /* shut down our window */
      printf("(x,y,z) = (%f,%f,%f)\n",xdepth,ydepth,zdepth);
      exit(EXIT_SUCCESS); /* exit the program...normal termination. */
      break; /* just in case */
    case SPACE_KEY: /* stop/restart temporarily */
      RYangle = 0; RXangle = 0; RYdangle = 0; RXdangle = 0;
      STEPS_PER_DRAW=1024*!STEPS_PER_DRAW;if (!STEPS_PER_DRAW){ printf("PROCESS HALTED\n");}
      break;
    case 'w': FIDDLE_PARAMETER++; if (FIDDLE_PARAMETER > 25){ FIDDLE_PARAMETER = -10;} break;
    case 's': FIDDLE_PARAMETER--; if (FIDDLE_PARAMETER < -10){ FIDDLE_PARAMETER = 25;} break;
    case 'd': 
      switch (FIDDLE_PARAMETER){
      case -10: GRAYSCALE = !GRAYSCALE; break;
      case 0: nv=0;nt1=0,nt2=1; tab2 = nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES; CS__[tab2] *= pow(2,1.0/8.0); break;
      case 1: nv=0;nt1=0,nt2=2; tab2 = nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES; CS__[tab2] *= pow(2,1.0/8.0); break;
      case 2: nv=0;nt1=0,nt2=3; tab2 = nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES; CS__[tab2] *= pow(2,1.0/8.0); break;
      case 3: nv=1;nt1=3,nt2=0; tab2 = nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES; CS__[tab2] *= pow(2,1.0/8.0); break;
      case 4: nv=1;nt1=3,nt2=1; tab2 = nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES; CS__[tab2] *= pow(2,1.0/8.0); break;
      case 5: nv=1;nt1=3,nt2=2; tab2 = nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES; CS__[tab2] *= pow(2,1.0/8.0); break;
      case 6: nv=1;nt1=3,nt2=3; tab2 = nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES; CS__[tab2] *= pow(2,1.0/8.0); break;
      case 7: nv=0;nt1=2,nt2=1; tab2 = nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES; CS__[tab2] *= pow(2,1.0/8.0); break;
      case 8: nv=0;nt1=2,nt2=2; tab2 = nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES; CS__[tab2] *= pow(2,1.0/8.0); break;
      case 9: nv=0;nt1=2,nt2=3; tab2 = nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES; CS__[tab2] *= pow(2,1.0/8.0); break;
      case 10: break;
      case 11: break;
      case 12: CURRENT_INJECTION_S *= pow(2,1.0/8.0); break;
      case 13: break;
      case 14: break;
      case 15: STD_VIEW *= pow(2,1.0/8.0); break;
      case 18: GLOBAL_DTmax *= 2; break;
      case 20: STEPS_PER_DRAW *= 2; if (STEPS_PER_DRAW==0){ STEPS_PER_DRAW=1;} break;
      case 21: DRAW_FLAG++; break;
      case 22: GLOBAL_SPACE_SMOOTHER++; break;
      case 23: FIDDLE_ROW += 1; break;
      case 24: INPUT_CONTRAST+=0.05; break;
      case 25: break;
      default: break;}
      break;
    case 'D':
      switch (FIDDLE_PARAMETER){
      case 0: nv=0;nt1=1,nt2=1; tab2 = nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES; CS__[tab2] *= pow(2,1.0/8.0); break;
      case 1: nv=0;nt1=1,nt2=2; tab2 = nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES; CS__[tab2] *= pow(2,1.0/8.0); break;
      case 2: nv=0;nt1=1,nt2=3; tab2 = nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES; CS__[tab2] *= pow(2,1.0/8.0); break;
      case 3: nv=2;nt1=3,nt2=0; tab2 = nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES; CS__[tab2] *= pow(2,1.0/8.0); break;
      case 4: nv=2;nt1=3,nt2=1; tab2 = nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES; CS__[tab2] *= pow(2,1.0/8.0); break;
      case 5: nv=2;nt1=3,nt2=2; tab2 = nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES; CS__[tab2] *= pow(2,1.0/8.0); break;
      case 6: nv=2;nt1=3,nt2=3; tab2 = nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES; CS__[tab2] *= pow(2,1.0/8.0); break;
      case 9: break;
      case 10: break;
      case 11: break;
      case 12: CURRENT_INJECTION_D *= pow(2,1.0/8.0); break;
      case 13: break;
      case 14: break;
      case 15: MEAN_VIEW += 0.1; break;
      case 21: DRAW_FLAG2++; break;
      case 23: FIDDLE_COL += 1; break;
      case 24: INPUT_PULSE=1; break;
      default: break;}
      break;
    case 'a':
      switch (FIDDLE_PARAMETER){
      case -10: GRAYSCALE = !GRAYSCALE; break;
      case 0: nv=0;nt1=0,nt2=1; tab2 = nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES; CS__[tab2] /= pow(2,1.0/8.0); break;
      case 1: nv=0;nt1=0,nt2=2; tab2 = nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES; CS__[tab2] /= pow(2,1.0/8.0); break;
      case 2: nv=0;nt1=0,nt2=3; tab2 = nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES; CS__[tab2] /= pow(2,1.0/8.0); break;
      case 3: nv=1;nt1=3,nt2=0; tab2 = nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES; CS__[tab2] /= pow(2,1.0/8.0); break;
      case 4: nv=1;nt1=3,nt2=1; tab2 = nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES; CS__[tab2] /= pow(2,1.0/8.0); break;
      case 5: nv=1;nt1=3,nt2=2; tab2 = nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES; CS__[tab2] /= pow(2,1.0/8.0); break;
      case 6: nv=1;nt1=3,nt2=3; tab2 = nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES; CS__[tab2] /= pow(2,1.0/8.0); break;
      case 7: nv=0;nt1=2,nt2=1; tab2 = nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES; CS__[tab2] /= pow(2,1.0/8.0); break;
      case 8: nv=0;nt1=2,nt2=2; tab2 = nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES; CS__[tab2] /= pow(2,1.0/8.0); break;
      case 9: nv=0;nt1=2,nt2=3; tab2 = nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES; CS__[tab2] /= pow(2,1.0/8.0); break;
      case 10: break;
      case 11: break;
      case 12: CURRENT_INJECTION_S /= pow(2,1.0/8.0); break;
      case 13: break;
      case 14: break;
      case 15: STD_VIEW /= pow(2,1.0/8.0); break;
      case 18: GLOBAL_DTmax /= 2; break;
      case 20: STEPS_PER_DRAW /= 2; if (STEPS_PER_DRAW < 1){ STEPS_PER_DRAW=0;} break;
      case 21: DRAW_FLAG--; break;
      case 22: GLOBAL_SPACE_SMOOTHER--; if (GLOBAL_SPACE_SMOOTHER<0){ GLOBAL_SPACE_SMOOTHER=0;} break;
      case 23: FIDDLE_ROW -= 1; break;
      case 24: INPUT_CONTRAST-=0.05; break;
      case 25: break;
      default: break;}
      break;
    case 'A':
      switch (FIDDLE_PARAMETER){
      case 0: nv=0;nt1=1,nt2=1; tab2 = nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES; CS__[tab2] /= pow(2,1.0/8.0); break;
      case 1: nv=0;nt1=1,nt2=2; tab2 = nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES; CS__[tab2] /= pow(2,1.0/8.0); break;
      case 2: nv=0;nt1=1,nt2=3; tab2 = nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES; CS__[tab2] /= pow(2,1.0/8.0); break;
      case 3: nv=2;nt1=3,nt2=0; tab2 = nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES; CS__[tab2] /= pow(2,1.0/8.0); break;
      case 4: nv=2;nt1=3,nt2=1; tab2 = nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES; CS__[tab2] /= pow(2,1.0/8.0); break;
      case 5: nv=2;nt1=3,nt2=2; tab2 = nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES; CS__[tab2] /= pow(2,1.0/8.0); break;
      case 6: nv=2;nt1=3,nt2=3; tab2 = nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES; CS__[tab2] /= pow(2,1.0/8.0); break;
      case 9: break;
      case 10: break;
      case 11: break;
      case 12: CURRENT_INJECTION_D /= pow(2,1.0/8.0); break;
      case 13: break;
      case 14: break;
      case 15: MEAN_VIEW -= 0.1; break;
      case 21: DRAW_FLAG2--; break;
      case 23: FIDDLE_COL -= 1; break;
      case 24: INPUT_PULSE=0; break;
      default: break;}
      break;
    default: break;}
    break;
  case 7: /* snx style */
    switch (key) {    
    case ESCAPE_KEY: /* stop everything */
      glutDestroyWindow(GLUTWINDOWNUMBER); /* shut down our window */
      printf("(x,y,z) = (%f,%f,%f)\n",xdepth,ydepth,zdepth);
      exit(EXIT_SUCCESS); /* exit the program...normal termination. */
      break; /* just in case */
    case SPACE_KEY: /* stop/restart temporarily */
      RYangle = 0; RXangle = 0; RYdangle = 0; RXdangle = 0;
      STEPS_PER_DRAW=1*!STEPS_PER_DRAW;if (!STEPS_PER_DRAW){ printf("PROCESS HALTED\n");}
      break;
    case 'w': FIDDLE_PARAMETER++; if (FIDDLE_PARAMETER > 25){ FIDDLE_PARAMETER = -10;} break;
    case 's': FIDDLE_PARAMETER--; if (FIDDLE_PARAMETER < -10){ FIDDLE_PARAMETER = 25;} break;
    case 'd': 
      switch (FIDDLE_PARAMETER){
      case -10: GRAYSCALE = !GRAYSCALE; break;
      case 13: nv=0;nt1=0,nt2=0; tab2 = nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES; CS__[tab2] *= pow(2,1.0/8.0); if (CS__[tab2]==0){ CS__[tab2]+=1;} break;
      case 14: CS_ORN_[TYPENAME_REGISTRY_snx_PN] *= pow(2,1.0/8.0); if (CS_ORN_[TYPENAME_REGISTRY_snx_PN]==0){ CS_ORN_[TYPENAME_REGISTRY_snx_PN] += 1;} break;
      case 15: STD_VIEW *= pow(2,1.0/8.0); break;
      case 18: GLOBAL_DTmax *= 2; break;
      case 20: STEPS_PER_DRAW *= 2; if (STEPS_PER_DRAW==0){ STEPS_PER_DRAW=1;} break;
      case 21: DRAW_FLAG++; break;
      case 22: GLOBAL_SPACE_SMOOTHER++; break;
      case 23: FIDDLE_ROW += 1; break;
      case 24: INPUT_CONTRAST+=0.05; break;
      case 25: break;
      default: break;}
      break;
    case 'D':
      switch (FIDDLE_PARAMETER){
      case 15: MEAN_VIEW += 0.1; break;
      case 21: DRAW_FLAG2++; break;
      case 23: FIDDLE_COL += 1; break;
      case 24: INPUT_PULSE=1; break;
      default: break;}
      break;
    case 'a':
      switch (FIDDLE_PARAMETER){
      case -10: GRAYSCALE = !GRAYSCALE; break;
      case 13: nv=0;nt1=0,nt2=0; tab2 = nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES; CS__[tab2] /= pow(2,1.0/8.0); break;
      case 14: CS_ORN_[TYPENAME_REGISTRY_snx_PN] /= pow(2,1.0/8.0); break;
      case 15: STD_VIEW /= pow(2,1.0/8.0); break;
      case 18: GLOBAL_DTmax /= 2; break;
      case 20: STEPS_PER_DRAW /= 2; if (STEPS_PER_DRAW < 1){ STEPS_PER_DRAW=0;} break;
      case 21: DRAW_FLAG--; break;
      case 22: GLOBAL_SPACE_SMOOTHER--; if (GLOBAL_SPACE_SMOOTHER<0){ GLOBAL_SPACE_SMOOTHER=0;} break;
      case 23: FIDDLE_ROW -= 1; break;
      case 24: INPUT_CONTRAST-=0.05; break;
      case 25: break;
      default: break;}
      break;
    case 'A':
      switch (FIDDLE_PARAMETER){
      case 15: MEAN_VIEW -= 0.1; break;
      case 21: DRAW_FLAG2--; break;
      case 23: FIDDLE_COL -= 1; break;
      case 24: INPUT_PULSE=0; break;
      default: break;}
      break;
    default: break;}
    break;
  default: break;}
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
      if (CLUSTERDATA_BOTHER){ clusterdatarareset(GLOBAL_CDRA);}
      if (PTREE_BOTHER){ ptreereset(GLOBAL_PTREE); /* ptreerate(GLOBAL_PTREE); ptreedump_starter(GLOBAL_PTREE,NULL,2,0,0,0,+1,-1,0); */} 
      if (HHLIB_BOTHER){ hhlib_project_printf(GLOBAL_HHLIB,0,0);}
      dumpoutput("al_output");}
    else{ if (!STEPS_PER_DRAW){ computestep();}}
    break;
  case GLUT_KEY_HOME:
    mod = glutGetModifiers();
    if (mod==GLUT_ACTIVE_SHIFT){ INPUT_PULSE = (INPUT_PULSE==1?0:1);}
    else{ GLOBAL_DRAW_FLAG = !GLOBAL_DRAW_FLAG;}
    break;
  default:
    break;}
}

GLvoid DrawGLScene(GLvoid)
{
  int index=0;
  void **vra=NULL;
  for (index=0;index<STEPS_PER_DRAW;index++){ if (!RUN_DONE){ computestep();}}
  /* draw stuff */
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // Clear The Screen And The Depth Buffer
  glLoadIdentity(); // Reset The View
  glTranslatef(xdepth,ydepth,zdepth);
  glRotatef(RYangle,0,1,0);
  glRotatef(RXangle,1,0,0);
  if (GLOBAL_DRAW_FLAG){
    if (PTREE_BOTHER){ Drawptree(DRAW_FLAG,GLOBAL_PTREE,2,1,0);}
    if (GLOBAL_INPUT_PTREE!=NULL){ Drawptree(DRAW_FLAG,GLOBAL_INPUT_PTREE,2,1,-1);}
    if (POWER_BOTHER){ Drawpower(GLOBAL_POWER,4,2,2,2.5);}
    if (RHO_BOTHER){ Drawrho(GLOBAL_RHO,1,1,3,0);}
    if (CLUSTERDATA_BOTHER){ Drawcd(GLOBAL_CDRA,1,0.5,3,0);}
    if (CAICOR_BOTHER){ Drawcaicor(GLOBAL_CAICOR,1,0.5,3,0);}
    if (SNXDATA_BOTHER){ Drawsnxdata(GLOBAL_SNXDATA,1,1,0);}
    if (ISI_BOTHER){ Drawisi(GLOBAL_ISI,1,1,3,0);}
/*     DrawNra(GLOBAL_Nra,1,1,1,3.5); */
    switch(GLOBAL_NEURON_MODEL){
    case 0: case 1: case 2: case 3: case 4: Drawcra(GLOBAL_Nra,1,1,-1); break;
    case 5: /* mainak style */ Drawcra(GLOBAL_Nra,1,1,-1); break;
    case 6: /* wilson style */ /* Drawcra(GLOBAL_Nra,1,1,-1); */ break;
    case 7: /* snx style */ Drawsnx(2,2,-1); break;
    default: break;}
  }
  if (LATTICE3D_BOTHER){ 
    if (DRAW_FLAG<2){ Drawlattice3d_edge(GLOBAL_LATTICE3D,1,0,0);}
    else if (DRAW_FLAG>=2){ 
      vra=(void **)tcalloc(2,sizeof(void *)); index=1; vra[0]=&index; vra[1] = &GLOBAL_LATTICE3D_LENGTH;
      Drawlattice3d(GLOBAL_LATTICE3D,DRAW_FLAG%2,1,0,0,&spack3d_temp_vector_draw,vra);
      tfree(vra);vra=NULL;
      vra=(void **)tcalloc(2,sizeof(void *)); index=2; vra[0]=&index; vra[1] = &GLOBAL_LATTICE3D_LENGTH;
      Drawlattice3d(GLOBAL_LATTICE3D,DRAW_FLAG%2,1,5,0,&spack3d_temp_vector_draw,vra);
      tfree(vra);vra=NULL;}}
  Drawmenu(.1,-1,1);
  glutSwapBuffers(); /* since double buffered */
  RYangle += RYdangle;
  RXangle += RXdangle;
}

#endif /* NOTGRAPH */

