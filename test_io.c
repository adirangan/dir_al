/* This is a tutorial which will show you how to read and write simple input files.

   You can compile it by the command:
   gcc -o test_io.out test_io.c;
   and you can run it easily by typing the command
   nice ./test_io.out < test_io_input.in;

   The input file is "test_io_input.in".
   Each line should contain a variable name immediately followed by an "=" equals sign, and then a " " space and then the value.
   Note that each line is read in sequence, so the final value for any given variable will be read from the last line it appears in.
   The function "readinput" terminates when the variablename "END_OF_INPUT" is reached.

   The program generates an outputfile named "test_io_output.out" which is a list of global variables ready to be read back in.
   For example, the command
   nice ./test_io.out < test_io_output.out
   should work just fine (after you create "test_io_output.out" of course).

   The whole point of this is that you can change the input file without having to recompile the original "test_io.out" executable.
   If you would like to decipher things a little more easily, just set GLOBAL_verbose to 1 in the input file.
*/

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>

/* Here are the global variables */
int GLOBAL_verbose=0;
char GLOBAL_STRING[64]="default";
int GLOBAL_INT=0;
double GLOBAL_DOUBLE=0;

/* Here are the input reading and output dumping functions */
void readinput();
int updateglobals(char *);
void dumpoutput(char *);


void readinput()
{
  /* This reads the piped input file, or standard input.
     the variable names should not be longer than 128 characters */
  int verbose=0;
  char vname[128],equals[128],space[128],semicolon[128];
  do{
    scanf("%[^=]",vname);scanf("%s",equals);scanf("%c",space);updateglobals(vname);scanf("%c",semicolon);
    if (verbose){ printf("At this point variable name is (%s), equals is (%s), semicolon is (%s)\n",vname,equals,semicolon);} 
    scanf("%c",semicolon);}
  while (strcmp(vname,"END")!=0);
}

int updateglobals(char *vname)
{
  /* given a variable name, scan the appropriate format into the appropriate global variable */
  int if_we_could_understand_vname=1;
  if (GLOBAL_verbose){ printf(" %% %% %% [entering updateglobals] variablename=%s\n",vname);}
  if (strcmp(vname,"GLOBAL_verbose")==0){ 
    scanf("%d",&GLOBAL_verbose); 
    if (GLOBAL_verbose){ printf(" %% %% %% %s read to be %d\n",vname,GLOBAL_verbose);}}
  else if (strcmp(vname,"GLOBAL_STRING")==0){ 
    scanf("%[^,;]",GLOBAL_STRING);
    if (GLOBAL_verbose){ printf(" %% %% %% %s read to be %s\n",vname,GLOBAL_STRING);}}
  else if (strcmp(vname,"GLOBAL_INT")==0){ 
    scanf("%d",&GLOBAL_INT); 
    if (GLOBAL_verbose){ printf(" %% %% %% %s read to be %d\n",vname,GLOBAL_INT);}}
  else if (strcmp(vname,"GLOBAL_DOUBLE")==0){ 
    scanf("%lf",&GLOBAL_DOUBLE); 
    if (GLOBAL_verbose){ printf(" %% %% %% %s read to be %lf\n",vname,GLOBAL_DOUBLE);}}
  else if (strcmp(vname,"END_OF_INPUT")==0){ /* do nothing */ 
    if (GLOBAL_verbose){ printf("end of input reached\n");}}
  else{ 
    if (GLOBAL_verbose){ printf(" %% %% %% vname %s is not understandable\n",vname);} 
    if_we_could_understand_vname=0;}
  return if_we_could_understand_vname;
}

void dumpoutput(char *filename)
{
  /* prints an output file */
  char text[256];
  char vname[128];
  FILE *fp=NULL;
  sprintf(text,"./%s",filename);
  if (filename==NULL || (fp = fopen(text, "w")) == NULL){ printf("dumpoutput to stdout\n"); fp = stdout;}
  sprintf(vname,"GLOBAL_verbose"); fprintf(fp,"%s= %d;\n",vname,GLOBAL_verbose);
  sprintf(vname,"GLOBAL_STRING"); fprintf(fp,"%s= %s;\n",vname,GLOBAL_STRING);
  sprintf(vname,"GLOBAL_INT"); fprintf(fp,"%s= %d;\n",vname,GLOBAL_INT);
  sprintf(vname,"GLOBAL_DOUBLE"); fprintf(fp,"%s= %0.16lf;\n",vname,GLOBAL_DOUBLE);
  sprintf(vname,"END"); fprintf(fp,"%s= 0;\n",vname);
  if (fp!=stdout){ fclose(fp);}
}

int main(int argc, char **argv) 
{
  readinput();
  printf(" %% main: GLOBAL_STRING %s, GLOBAL_DOUBLE %+0.6f, GLOBAL_INT %+0.0d\n",GLOBAL_STRING,GLOBAL_DOUBLE,GLOBAL_INT);
  dumpoutput("test_io.out");
  return 1;
}
