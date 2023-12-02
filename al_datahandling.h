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
  int lpower_bother;
  int lpower_window_length;
  int lpower_update_every;
  int lpower_length;
  double *lpowerdata;
};

/* this holds a power */
struct power
{
  struct neuronarray *Nra;
  int indexing_ntype_length;
  int *indexing_ntype_checkout;
  int *indexing_ntype_refile;
  int indexing_nvar_length;
  int *indexing_nvar_checkout;
  int *indexing_nvar_refile;
  double update_every;
  int length;
  int update_number;
  struct strobe *tst;
  struct strobe ****strarara;
  double ***vpowrarara;
  double *firstvarpower;
  int correlation_bother;
  double *autocorrelation;
  double *crosscorrelation;
  double *maxra;
  double *minra;
};

/* this holds a rho */
struct rho
{
  struct neuronarray *Nra;
  int length;
  struct llist **neuronllistra;
  int indexing_nvar_length;
  int *indexing_nvar_checkout;
  int *indexing_nvar_refile;
  double *maxra;
  double *minra;
  int *nbinra;
  int nbins;
  double **rhora;
  double total_time;
};

/* this holds a region */
struct region
{
  int region_type;
  int label;
  double last_event;
  double event_within;
  int event_threshold;
  struct llist *neuronllist;
  struct pnode *pn; /* for yggdrasil->pp regions */
};

/* this holds a clusterdatara */
struct clusterdatara
{
  struct clusterarray *gli;
  int power_bother;
  int indexing_nvar_length;
  int *indexing_nvar_checkout;
  int indexing_nvar_refile_length;
  int *indexing_nvar_refile;
  int power_update_every;
  int power_length;
  int power_update_number;
  struct strobe *tst;
  struct strobe ****strarara;
  int vpow_bother;
  double ***vpowrarara;
  struct strobe **stra_lfp;
  double *maxra;
  double *minra;
  int ptree_bother;
  struct ptree **pra;
  int trajectory_window;
  int lfp_gammalo;
  int lfp_gammahi;
  void *vtemp1;
  void *vtemp2;
  void *vtemp3;
};

/* this holds a caicor */
struct caicor
{
  struct llist *LN;
  int indexing_nvar_length;
  int *indexing_nvar_checkout;
  int indexing_nvar_refile_length;
  int *indexing_nvar_refile;
  double power_update_every;
  double *maxra;
  double *minra;
  int nbins;
  struct hist **hra;
  struct hist **hstra;
  double *mrara;
  double *mrarara;
};

/* this holds a hhlib_hist */
struct hhlib_hist
{
  int nvars;
  int indexing_trigger;
  int *libuseflagra;
  int *logflagra;
  double *maxra;
  double *minra;
  int nbins;
  int ncumulants;
  double *data;
  double *data_project;
};

/* this holds a hhlib */
struct hhlib
{
  struct neuronarray *Nra;
  double tau_skip;
  double trigger_value;
  int indexing_trigger;
  int indexing_nvar_length;
  int *indexing_nvar_checkout;
  int indexing_nvar_refile_length;
  int *indexing_nvar_refile;
  int *logflagra;
  int *libuseflagra;
  double *maxra;
  double *minra;
  int nbins;
  int use_after;
  double ***trigger_rarara;
  struct hhlib_hist **hlhra;
};

/* this holds a isi */
struct isi
{
  struct neuronarray *Nra;
  struct llist *L;
  struct llist *L2;
  int length_max;
  double maxisi;
  double minisi;
  int nbins;
  double *ra;
};

/* this holds a spack3d */
struct spack3d
{
  struct lattice3d *l3d;
  int i; int j; int k; double x; double y; double z; double value;
  struct llist *L;
  struct spack3d **adjacency_ra_s;
  double *adjacency_ra_v;
  void *temp;
};

/* this holds a lattice3d */
struct lattice3d
{
  int tet_vs_cube;
  int input1max;
  int input2max;
  int input3max;
  int input2max_scale;
  int input3max_scale;
  int total_length;
  struct spack3d **lattice;
  int nedges;
  double edge_v_max;
  double edge_v_min;
  double edge_v_mean;
  double edge_v_stdev;
};

/* this holds a snxdata */
struct snxdata
{
  struct neuronarray *Nra;
  struct strobe **stra;
  int length;
  int lookback;
  double *chain_0;
  double *chain_1_;
  double total_time;
  double *input;
  double *connectivity;
  double abc;
  double bca;
};

/* Here are the odor functions */
struct odor * odormake(int,int *,int,char **,struct odor *);
void odortfree(struct odor *);
void odorplusequals(struct odor *,double,struct odor *);
void odorfprintf(FILE *,struct odor *);
void odorfprintf_full(FILE *,struct odor *);
void odor_full_fwrite(char *,struct odor *);
struct odor * odor_full_fread(char *);
void granule_local_findint(char *,int *,int *);
void granule_findint(int,int,int *,int *);
void granule_int2input(struct odor *,int,int,int,int);
void granule_input2int(struct odor *,int,int *,int,int *);
void granule_int2input_or_input2int(int,struct odor *,int,int *,int,int *);
void granule_local_plot();
void granule_plot(int,int);
double granule_obtain_reweight_vs_bipoo_angle(struct ptree *,double *);
double * granule_obtain_bipoo_matrix(char *,double,double);
struct ptree * ptree_obtain_bipoo_reweight(char *);
void ptree_test_mcpit();
struct ptree * granule_obtain_reweight(char *,char *,double,double);
/* int regionra_neuron2label(struct ptree *,struct neuronarray *,struct neuron *); */
double reweight2adder(struct ptree *,struct neuronarray *,struct neuron *,struct neuron *,double);
double reweight2multiplier(struct ptree *,struct neuronarray *,struct neuron *,struct neuron *,double);
double * nsphere_index_read(int,int,int);
void nsphere_index_generate(int,int);
void synaptic_sphere_findint2input_or_input2int(int,struct ptree *,int,int *);
void synaptic_sphere_local_findint(char *,int *);
void synaptic_sphere_prediction_check(int,int,char *,char *,double,double,double *,double *);
void synaptic_sphere_prediction_plot(int);

/* Here are system functions */
int cleanupoutput();

/* Here are the strobe and strobetrace functions */
struct strobe * strobemake(int,double,int,int,int,int);
void strobeupdate(struct strobe *,double,double,double);
void strobeupdate_bkp(struct strobe *,double,double,double);
void strobeupdate_old(struct strobe *,double,double);
void stradump(struct strobe **,int,int,char *);
void strobetfree(struct strobe *);
void strobereset(struct strobe *);
double * stra2ra(struct strobe **,int,int,int,double *);
double * ra2pca(double *,int,int,double *);
void pca2jpg(double *,int,char *);

/* Here are the power functions */
struct power * powermake(struct neuronarray *,int,int,int *,int *,int,int *,int *,double *,double *,int);
void powertfree(struct power *);
void powerupdate(struct power *,double,double);
void powerdump(struct power *,char *,int);

/* Here are the rho functions */
struct rho * rhomake(struct neuronarray *,int,int,int *,int *,double *,double *,int *);
void rhotfree(struct rho *);
int rhoupdate_getbins(struct rho *,struct neuron *,int *);
void rhoupdate(struct rho *,double,double);
void rho_neuronllistmake(struct rho *,int);
void rhodump(struct rho *,char *,int);

/* Here are the region functions */
void regionramake(struct ptree *,double,int,int);
void regionratfree(struct ptree *);
int region_has_event(struct region *,double,double);
int yggdrasil_region_has_event(struct region *,struct ptree *,double,double);

/* Here are the clusterdata functions */
struct clusterdatara * clusterdataramake(struct clusterarray *,int,int,int,int *,int,int *,int,int,int,double *,double *,int,int,int,int);
void clusterdataratfree(struct clusterdatara *);
void clusterdataraupdate(struct clusterdatara *,double,double);
void clusterdatarareset(struct clusterdatara *);
void clusterdataradump(struct clusterdatara *,char *,int);

/* Here are the caicor functions */
struct caicor * caicormake(struct neuronarray *,int,int *,int,int *,double *,double *,int);
void caicortfree(struct caicor *);
void caicorupdate(struct caicor *,double,double,int);
void caicordump(struct caicor *,char *);

/* Here are the hhlib functions */
struct hhlib_hist *hhlib_histmake(int,int,int *,int *,double *,double *,int,int);
void hhlib_histupdate(struct hhlib_hist *,double *,double *);
void hhlib_project_printf(struct hhlib *,int,int);
void hhlib_histproject(struct hhlib_hist *);
void hhlib_histproject_printf(struct hhlib_hist *,int);
double hhlib_histN(struct hhlib_hist *,double *);
void hhlib_histmean(struct hhlib_hist *,double *,double *);
void hhlib_histtfree(struct hhlib_hist *);
struct hhlib * hhlibmake(struct neuronarray *,double,int,int,int *,int,int *,int *,int *,double *,double *,int,int);
void hhlibupdate(struct hhlib *,double,double);
void hhlibprintf(struct hhlib *);
void hhlib_histdump(struct hhlib_hist *,FILE *);
struct hhlib_hist *hhlib_histread(FILE *);
void hhlibdump(struct hhlib *,char *);
struct hhlib * hhlibread(struct neuronarray *,char *);

/* Here are the isi functions */
struct isi * isimake(struct neuronarray *,int,double,double);
void isitfree(struct isi *);
void isiupdate(struct isi *,double);
void isidump(struct isi *,char *);

/* Here are the suite functions */

void suite_7_power_process_helper(char *);
void suite_7_power_deposit_helper(int,int,int,int,int);
void system_monitor(struct neuronarray *,double,double);

/* Here are the lattice3d functions */
struct spack3d * spack3dmake(struct lattice3d *,int,int,int,double,double,double,double);
void spack3dprintf(struct spack3d *,void (*)(void *,void *),void *);
void spack3dtfree(struct spack3d *,void (*)(void *,void *),void *);
struct spack3d * lattice3dget(struct lattice3d *,int,int,int);
void spack3dlink(struct spack3d *,struct spack3d *);
double vra2vra_diff_sym(void *,void *,void *);
double vra2vra_rq_sym(void *,void *,void *);
void lattice3d_edge_eval(struct lattice3d *,double (*)(void *,void *,void *),void *);
void latticelink(struct lattice3d *);
struct lattice3d * lattice3dmake(int,int,int,int);
void lattice3d_edge_v_stats(struct lattice3d *);
void spack3d_temp_vector_dump(void *,FILE *,double,double,double,double,double,void *);
void lattice3d_temp_dump(struct lattice3d *,char *,void (*)(void *,FILE *,double,double,double,double,double,void *),void *);
void lattice3d_edge_v_dump(struct lattice3d *,char *);
void lattice3dtfree(struct lattice3d *,void (*)(void *,void *),void *);
struct ptree * binary_projection_s2p_make(double *,int,int);
double * binary_projection_s2s(int,int,int);
double * binary_projection_s2p(int,int,int,int);
struct lattice3d * threetree_test(int,int,double *);
void lattice3d_int2input_or_input2int(int,int,int,int,int,int *,int *,int *,double *,double *,double *);
double * bp_projection(int,int);
double * ds_projection(int);
void ideal_2_statetree(int,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *);
void ra2lss(double *,int,int,double *,int,double *);
void ra2svd(double *,int,int,double *,double *,double *);
void ra2eig(double *,int,double *,double *,double *);

/* Here are the snxdata functions */
struct snxdata *snxdatamake(struct neuronarray *,int);
void snxdatatfree(struct snxdata *);
void snxdataupdate(struct snxdata *,double,double);


