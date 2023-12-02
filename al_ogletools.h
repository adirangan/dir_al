#ifndef NOTGRAPH

int WritePPMFile(const char *,GLubyte *,int,int,int);
int DumpWindow(const char *,int,int); 
void ftexto(float,float,float,char *);
GLvoid InitGL(GLsizei,GLsizei);
GLvoid ReSizeGLScene(GLsizei,GLsizei);
GLvoid loglogra(void *,char *,int,double,double,double,double,double,double,double,double,double);
GLvoid Plotra(void *,char *,int,double,double,double,double,double,double,double,double,double);
GLvoid Drawra(void *,char *,int,int,int,int,double,double,int,char *,double,double,double,double);
GLvoid Drawconnectivity(struct neuronarray *,double,double,double,double);
GLvoid DrawNra(struct neuronarray *,double,double,double,double);
GLvoid Drawpower(struct power *,double,double,double,double);
GLvoid Drawstra(struct strobe **,int,double,double,double,double,double,double,double,double,double,int);
double norm(double,double,double);
void crossproduct(double,double,double,double,double,double,double *,double *,double *);
void rotate100(double *,double *,double *,double,double);
void rotate010(double *,double *,double *,double,double);
GLvoid glgrid(int,int,int,double,double,double,double,double,double,double);
GLvoid glbox(double,double,double,double,double,double);
GLvoid glrect(double,double,double,double,double,double,double);
GLvoid glnum(double,double,double,double,double,double,double);
GLvoid gldot(int,double,double,double,double,double,double);
GLvoid glring(int,double,double,double,double,double,double,double);
GLvoid glarrow(int,double,double,double,double,double,double,double,double);
GLvoid glcube(double,double,double,double,double,double,double);
GLvoid glarc(int,double,double,double,double,double,double,double,double,double);
GLvoid Drawcolorbar(int,double,double,double);
GLvoid Drawmenu(double,double,double);
void keyPressed(unsigned char,int,int); 
void specialKeyPressed(int,int,int); 
GLvoid DrawGLScene(GLvoid);

#endif /* NOTGRAPH */

