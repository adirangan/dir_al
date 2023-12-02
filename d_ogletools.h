#ifndef NOTGRAPH

int WritePPMFile(const char *, GLubyte *,int,int,int);
int DumpWindow(const char *,int,int);
void ftexto(float,float,float,char *);
GLvoid InitGL(GLsizei,GLsizei);
GLvoid ReSizeGLScene(GLsizei,GLsizei);
GLvoid Plotra(void *,char *,int,double,double,double,double,double,double,double,double,double);
GLvoid loglogra(void *,char *,int,double,double,double,double,double,double,double,double,double);
GLvoid Drawra(void *,char *,int,int,int,int,int,char *,double,double,double,double);
GLvoid DrawLR(double,double,double);
GLvoid DrawNra(double,double,double);
GLvoid Drawlgn(double,double,double,double);
GLvoid Drawrtc(double,double,double,double);
GLvoid Drawstrobetrace(double,double,double);
GLvoid Drawstra(struct strobe **,int,double,double,double,double,double,double,double,double,int);
GLvoid Drawtuningcurve(double,double,double,double);
GLvoid Drawlmitri(double,double,double,double);
GLvoid Drawptree(struct ptree *,double,double,double);
GLvoid glbox(double,double,double,double,double,double);
GLvoid glnum(double,double,double,double,double,double,double);
GLvoid gldot(int,double,double,double,double,double,double);
GLvoid glring(int,double,double,double,double,double,double,double);
GLvoid glarrow(int,double,double,double,double,double,double,double,double);
GLvoid glcube(double,double,double,double,double,double,double);
GLvoid Drawcolorbar(int,double,double,double);
GLvoid Drawmenu(double,double,double);
void keyPressed(unsigned char,int,int); 
void specialKeyPressed(int,int,int); 
GLvoid DrawGLScene(GLvoid);

#endif /* NOTGRAPH */

