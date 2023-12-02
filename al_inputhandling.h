/* This holds a odor */
struct odor
{
  int ntypes;
  int *lengthra;
  int base;
  int basedigits;
  double **rara;
};

/* Here are the input functions */ 
void readinput();
void updateglobals(char *);
void dumpoutput(char *);
void sparse_link_dump_ascii(struct neuronarray *,char *);
void sparse_link_read_ascii(struct neuronarray *,char *);
