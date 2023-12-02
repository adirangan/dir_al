/* This holds a rtc */
struct rtc
{
  int length;
  int nback;
  int nangles;
  int nphases;
  double total_time;
  double update_every;
  double update_last;
  int tab;
  int *anglera;
  int *phasera;
  /* stored as [ang_index + tab*NANGLES + t2s*length*NANGLES] */
  double *mra;
  double *Vra;
  double *sAra;
  double *sNra;
  double *sGra;
  double *VSra;
  /* stored as [ang_index + tab*NANGLES + t2s*length*NANGLES] */
  double *rtcmra;
  double *rtcVra;
  double *rtcsAra;
  double *rtcsNra;
  double *rtcsGra;
  double *rtcVSra;
  double *rtcra;
};

/* This holds a strobe */
struct strobe
{
  int length;
  int tab;
  double last_time;
  double update_timestep;
  double *data;
  int cycle_bother;
  int cyclenum;
  double *cycledata;
};

/* This holds a strobetrace */
struct strobetrace
{
  struct neuron *n; 
  int nangles;
  int *t2stotal;
  int *natotal;
  double timelength;
  int length;
  double update_timestep;
  int cycle_bother;
  struct strobe *tst;
  struct strobe *Vst;
  struct strobe *sAst;
  struct strobe *sNst;
  struct strobe *sGst;
  struct strobe *VSst;
  struct strobe **Vstra;
  struct strobe **sAstra;
  struct strobe **sNstra;
  struct strobe **sGstra;
  struct strobe **VSstra;
  struct strobe **mstra;
  int pcspiesacross; /* region around n */
  int pcspiestall; /* region around n */
  double *pcsra; /* preferred cortical states */
  struct strobe **patternstra_VS; /* VS-pcs spatial correlation */
  struct strobe **patternstra_sN; /* sN-pcs spatial correlation */
  struct strobe **pmstra; /* orientation-wise raster */
  struct strobe *patternc; /* instantaneous correlation between VS and sN spatial patterns */
  double *snvstc; /* sN-VS temporal correlation */
  struct avalanche *avalanche;
  double total_time;
};

/* this holds a tuningcurve */
struct tuningcurve
{
  int nangles;
  int nradius;
  double *mra;
  double *Vra;
  double *sAra;
  double *sNra;
  double *sGra;
  double *VSra;
};

/* this holds a lmitri */
struct lmitri
{
  int time_start;
  int time_length;
  int space_length;
  int rmax;
  int rmin;
  double *t2stotal;
  double *nupdates;
  double *mra;
  double *Vra;
  double *sAra;
  double *sNra;
  double *sGra;
  double *VSra;
};


/* this holds a region */
struct region
{
  int label;
  double last_event;
  double event_within;
  int event_threshold;
  struct llist *neuronllist;
  struct pnode *pn; /* for yggdrasil->pp regions */
};

/* this holds a pnode */
struct pnode
{
  struct pnode *parent;
  struct region *region;
  double weight;
  double relevance;
  struct llitem *childllitem;
  int broodsize;
  void *temp; /* for histograms and such */
};

/* this holds a ptree */
struct ptree
{
  int nregions;
  struct region **regionra;
  int nlegs;
  int legtime;
  int length;
  struct llitem **eventra;
  int tab;
  double update_every;
  double update_last;
  double total_time;
  int gate_flag;
  struct llitem *pretree;
  struct llitem *postree;
  struct hist *wh;
  struct hist *rh;
};

/* this holds a closet */
struct closet
{
  struct llist *neuronllist;
  int nlegs;
  int legtime;
  double minworth;
  struct ptree *p;
};

/* this holds a yggdrasil */
struct yggdrasil
{
  struct ptree *p;
  int ppnregions;
  int ppnlegs;
  int pplegtime;
  struct ptree *pp;
  int still_observing;
  int weight_minimum;
};

/* this holds a bonsai */
struct bonsai
{
  struct neuronarray *Nra;
  struct ptree **p;
  int nvsf;
  double switch_last;
  double switch_every;
};

/* this holds a hydra */
struct hydra
{
  double switch_every_two;
  double justontime;
  double stayontime;
  double last_switch;
  double dump_every;
  struct ptree *pjuston_1;
  struct ptree *pstayon_1;
  struct ptree *pjustoff_1;
  struct ptree *pjuston_2;
  struct ptree *pstayon_2;
  struct ptree *pjustoff_2;
  struct ptree *pstayoff;
  void (*switchon)(int);
  void (*switchoff)(int);
};

/* this holds a lyapunov */
struct lyapunov
{
  struct neuronarray *Nra;
  double total_time;
  double update_last;
  double update_every;
  int original_vs_perturbed;
  double perturbnorm;
  double jiggle;
  struct llist *L;
  struct llist *pL;
};

/* this holds a avalanche */
struct avalanche
{
  struct neuronarray *Nra;
  double update_last;
  double update_every;
  int logdmin;
  int logdmax;
  double **mraprev;
  double **mranext;
  double **mhist;
};

/* this holds a power */
struct power
{
  double update_last;
  double update_every;
  int N;
  int length;
  struct neuron **n;
  struct strobe *tst;
  struct strobe **raster;
  struct strobe **Vstra;
  struct strobe **sAstra;
  struct strobe **sNstra;
  struct strobe **sGstra;
  struct strobe **VSstra;
  double *Vpower;
  double *sApower;
  double *sNpower;
  double *sGpower;
  double *VSpower;
  double *VSbarpower;
  int update_number;
  double *autocorrelation;
  double *crosscorrelation;
};

/* this holds a taof */
struct taof
{
  struct neuronarray *Nra;
  int length;
  int step_every;
  struct strobe *input_contrast;
  struct strobe *input_spaceangle;
  struct strobe **f1r;
  struct strobe **f1i;
  struct strobe **f0;
};

/* this holds a seidcorr */
struct seidcorr
{
  char *filename_base;
  struct neuronarray *Nra;
  int space_bin_size;
  int time_bin_size;
  int length;
  double time_start;
  double time_end;
  int rows;
  int cols;
  struct strobe **Vstra;
  struct strobe **sAstra;
  struct strobe **sNstra;
  struct strobe **sGstra;
  struct strobe **VSstra;
};

int findrecordnumber();
int lgnpnmdump(struct lgn *,int,int);
int cortexdump(int,struct neuronarray *);
void connectionsdump(struct neuronarray *,int,char *);
int cleanupoutput();
int trialaverage(int,char **);

/* Here are the rtc functions */
struct rtc * rtcmake(int,int,int,int);
void rtctfree(struct rtc *);
void rtcupdate(struct rtc *,double,double);
void rtcdump(struct rtc *,int);

/* Here are the strobe and strobetrace functions */
struct strobe * strobemake(int,double,int);
void strobeupdate(struct strobe *,double,double,double);
void strobeupdate_sum(struct strobe *,double,double,double);
void strobeupdate_bkp(struct strobe *,double,double,double);
void strobeupdate_old(struct strobe *,double,double);
void stradump(struct strobe **,int,int,char *);
void strobetfree(struct strobe *);
void strobetraceupdate(struct strobetrace *,double,double);
struct strobetrace * strobetracemake(struct neuron *,int,double,double,int,double,int,int);
void strobetracedump(struct strobetrace *,int);
void strobetracetfree(struct strobetrace *);

/* Here are the tuningcurve functions */
struct tuningcurve * tuningcurvemake(int,int);
void tuningcurvetfree(struct tuningcurve *);
void tuningcurveupdate(struct tuningcurve *,double,double);
void tuningcurvedump(struct tuningcurve *,int);

/* Here are the lmitri functions */
struct lmitri * lmitrimake(int,int,double,double);
void lmitriupdate(struct lmitri *,double,double);
void lmitridump(struct lmitri *,int);
void lmitritfree(struct lmitri *);

/* Here are the ptree functions anew */
void regionramake(struct ptree *,double,int,int);
void regionratfree(struct ptree *);
int region_has_event(struct region *,double,double);
int yggdrasil_region_has_event(struct region *,struct ptree *,double,double);
struct pnode * pnodemake(struct pnode *,struct region *,double,double);
void pnodetfree(void *);
struct ptree * ptreemake(int,double,int,int,int,int);
void ptreemakepospre(struct ptree *);
void ptreereset(struct ptree *);
void ptreetfree(struct ptree *);
int region2region_compare_last_event(void *,void *);
int region2region_compare_label(void *,void *);
int region2pnode_compare_label(void *,void *);
int pnode2pnode_compare_label(void *,void *);
int pnode2pnode_compare_relevance(void *,void *);
int pnode2pnode_compare_weight(void *,void *);
void ptreeupdate_helper(struct ptree *,double,double,double,int);
void ptreeupdate(struct ptree *,double,double,double);
void pstrengthen_starter(struct ptree *,struct llist *,double);
void pstrengthen_helper(int,struct ptree *,struct pnode *,struct llitem *,struct llist *,int,double);
void pnodeprintf(struct pnode *,struct llitem *,int,int);
void pnodeprune_starter(int,struct ptree *,struct pnode *,struct llitem *,struct llist *);
void pnodeprune_helper(int,struct ptree *,struct pnode *,struct llist *);
double pnode_shear(double,double);
void pnodeZZZ_starter(int,struct ptree *,struct pnode *,struct llitem *,double,int);
void pnodeZZZ_helper(int,struct ptree *,struct pnode *,double);
void pnodebalance_starter(struct pnode *,struct llitem *);
void ptreerate(struct ptree *);
struct llist * ptreextract_starter(struct ptree *,double,int);
void ptreextract_helper(struct pnode *,struct llitem *,struct region **,double);
double pnode_worth(double,double,double);
void ptree2jpg_starter(struct ptree *,char *,double,double,double,double);
void ptree2jpg_helper(int,int,double,double,double,double,struct pnode *,struct llitem *,FILE *,FILE *,struct llist *);
void ptree2star_helper(int,int,int,int,double,double,double,double,struct pnode *,struct llitem *,FILE *,FILE *,struct llist *);
void ptreedump_starter(struct ptree *,char *,int,int,double,double,double,double);
void ptreedump_helper(int,struct pnode *,struct llitem *,FILE *,struct llist *);
struct ptree * ptreadback(char *);
void ptreeplusequals_starter(struct ptree *,struct ptree *);
void ptreeplusequals_helper(struct ptree *,struct pnode *,struct llitem *,struct ptree *,struct pnode *,struct llitem *);
struct ptree * ptreesetequalto_starter(struct ptree *);
void ptreesetequalto_helper(struct ptree *,struct pnode *,struct llitem *,struct ptree *,struct pnode *,struct llitem *);
struct ptree * ptreesubtptree_starter(struct ptree *,struct ptree *,int,int);
void ptreesubtptree_helper(struct ptree *,struct pnode *,struct llitem *,struct ptree *,struct pnode *,struct llitem *,double,int);
void ptreecomplement_helper(struct ptree *,struct pnode *,struct llitem *,struct ptree *,struct pnode *,struct llitem *,double,double);
void pnodefun_starter(struct pnode *,struct llitem *,int,int,int,void (*)(double *,double *));
void pnodefrob_starter(struct pnode *,struct llitem *,int,int,double *,double *);
void ptreex2(struct ptree *,struct pnode *,struct llitem *,struct ptree *,struct pnode *,struct llitem *,int,int,int,double *,double *,double *,double *);
void ptreerelent_breadth(struct ptree *,struct pnode *,struct llitem *,struct ptree *,struct pnode *,struct llitem *,int,int,int,int,double *,double *,double *);
void ptreerelent_depth(struct ptree *,struct pnode *,struct llitem *,struct ptree *,struct pnode *,struct llitem *,int,int,int,double *);
void pnodeclear_starter(struct pnode *,struct llitem *);
void pnodetimesd_starter(struct pnode *,struct llitem *,double,double);
void pnodestats_starter(struct pnode *,struct llitem *,int,int,int,int *,double *,double *,double *,double *,double *,double *,double *,double *);
void pnode2llist_starter(struct pnode *,struct llitem *,int,int,struct llist *);
void pnode2hist_starter(struct pnode *,struct llitem *,int,int,struct hist *,struct hist *);
void pnodehist_starter(struct pnode *,struct llitem *,int,int,struct hist *,struct hist *);
double ptreex2_match(int,int,double *,int *);
double ptreex1_smatch(int,int,double *,double,int *);
double ptreex1_xmatch(int,int,double *,double,int *);
void ptreex1(int,struct ptree **,struct pnode **,struct llitem **,int,struct ptree *,struct pnode *,struct llitem *,struct ptree *,struct pnode *,struct llitem *,int,int);
void ptreex1sbias(int ,struct ptree **,struct pnode **,struct llitem **,int,struct ptree *,struct pnode *,struct llitem *,struct ptree *,struct pnode *,struct llitem *,int,int,int,double *,double *);
void pnodeprintf_label(void *);
void pnodeshizuffle_starter(struct pnode *,struct llitem *,struct ptree *);
void ptreex2_disthist(double *,double *,int);
double logfactorial(double);
double binomial(double,double,double);
double Pmode(double,double,int,double *);
double multithresh_single(double *,int,double);
double multithresh(double *,int,double);
void discdata_threshold(double *,int,int,int ,char *);
void llistmedian(struct llist *,struct llist *,int (*)(void *,void *),int,void **,int *,int *);
void ptreeobsdistend(int,int,struct ptree *,struct pnode *,struct llitem *,int);
void obsdistdisc(int,int,struct llist **,double *,double *);
void ptreeobsdist(int,struct ptree *,struct pnode *,struct llitem *,int,struct ptree *,struct pnode *,struct llitem *,int,int);
void ptree_trialaverage_helper(int,int,int,int,int,int,int,int,char **,char *);
void ptree_trialaverage(int,char **);
void pnode2distribution_starter(struct ptree *,struct pnode *,struct llitem *,struct ptree *,struct pnode *,struct llitem *,int,int);
void pnode2distribution_helper(struct ptree *,struct pnode *,struct llitem *,int,int,struct llist **,int,int,struct hist **);
void pnode2distribution_ender(struct ptree *,struct pnode *,struct llitem *,int,int);
void ptree_trialdistribution_helper(int,int,int,int,char **,char *);
void ptree_trialdistribution(int,char **);


/* Here are the closet functions */
struct closet * closetmake(int,int,int,int,int,double);
void closettfree(struct closet *);
void closetremake(struct closet *,int);

/* Here are the yggdrasil functions */
struct yggdrasil * yggdrasilmake(int,double,int,int,int,int,int,int,int,int);
void yggdrasiltfree(struct yggdrasil *);
int pnode2pnode_labels_in_common(struct pnode *,struct pnode *);
void yggdrasilupdate(struct yggdrasil *,double,double,double);

/* Here are the bonsai functions */
struct bonsai * bonsaimake(struct neuronarray *,int,double,int,int,int,int);
void bonsaitfree(struct bonsai *);
void bonsaiupdate(struct bonsai *,double,double,double);
void bonsaidump(struct bonsai *,char *,int);

/* Here are the hydra functions */
void ptree_inputswitchon_dumb();
void ptree_inputswitchoff_dumb();
struct hydra * hydramake(double,double,double,void (*)(int),void (*)(int));
void hydratfree(struct hydra *);
void hydraupdate(struct hydra *,double,double,double);
void hydradump(struct hydra *);

/* Here are the lyapunov functions */
void nrasave(struct neuronarray *,char *);
void nraload(struct neuronarray *,char *);
double nraperturb(struct neuronarray *,double);
struct lyapunov * lyapunovmake(struct neuronarray *,double,double);
void lyapunovtfree(struct lyapunov *);
void lyapunovupdate(struct lyapunov *,double,double);
void lyapunovdump(struct lyapunov *,int);

/* Here are the avalanche functions */
struct avalanche * avalanchemake(struct neuronarray *,double,int,int);
void avalanchetfree(struct avalanche *);
void avalancheupdate(struct avalanche *,double,double);
void avalanchedump(struct avalanche *,int,char *);

/* Here are the power function */
struct power * powermake(struct neuronarray *,int,int);
void powertfree(struct power *);
void powerupdate(struct power *,double,double);
void powerdump(struct power *,char *,int);

/* Here are the taof functions */
struct taof * taofmake(struct neuronarray *,double,double);
void taoftfree(struct taof *);
void taofupdate(struct taof *,double,double);
double * taofdump_helper(double *,int,int,int,double *,double *);
void taofdump(struct taof *,char *,int);

/* Here are the seidcorr functions */
struct seidcorr * seidcorrmake(struct neuronarray *,int,int,int,double);
void seidcorrtfree(struct seidcorr *);
void seidcorrupdate(struct seidcorr *,double,double);
void seidcorrdump(struct seidcorr *);
void seidcorr_compile_helper(int,int,int,int);
void seidcorr_compile(int,char **);

/* here are the suite functions */
void cortex_data_tfree();
void lgnswitch(struct lgn *,double,double);
