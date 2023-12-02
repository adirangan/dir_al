#define TYPENAME_REGISTRY_PN (0)
#define TYPENAME_REGISTRY_LN (1)
#define TYPENAME_REGISTRY_wilson_ORN (0)
#define TYPENAME_REGISTRY_wilson_PN (1)
#define TYPENAME_REGISTRY_wilson_LNe (2)
#define TYPENAME_REGISTRY_wilson_LNi (3)
#define TYPENAME_REGISTRY_snx_PN (0)

/* These definitions are for use with 
   GLOBAL_NEURON_MODEL==0 --> vraevolve_if (conductance based integrate and fire)
   GLOBAL_NEURON_MODEL==1 --> vraevolve_adi (HH)
   GLOBAL_NEURON_MODEL==2 --> vraevolve_yisun (HH)
   GLOBAL_NEURON_MODEL==3 --> vraevolve_morrislecar
   GLOBAL_NEURON_MODEL==4 --> nllistevolve_cif (current based integrate and fire)
*/
#define VARNAME_REGISTRY_s_A (0)
#define VARNAME_REGISTRY_s_N (1)
#define VARNAME_REGISTRY_s_G (2)
#define VARNAME_REGISTRY_s_ORN (3)
#define VARNAME_REGISTRY_s_NMDA (4)
#define VARNAME_REGISTRY_m_Na (5)
#define VARNAME_REGISTRY_h_Na (6)
#define VARNAME_REGISTRY_m_K (7)
#define VARNAME_REGISTRY_Ca (8)
#define VARNAME_REGISTRY_Vs (9)
#define VARNAME_REGISTRY_Vd (10)
#define VARNAME_REGISTRY_spike_flag (11)
#define VARNAME_REGISTRY_inputrate (12)
#define VARNAME_REGISTRY_s_Na (13)
#define VARNAME_REGISTRY_s_K (14)
#define VARNAME_REGISTRY_s_Ca (15)
#define VARNAME_REGISTRY_s_KCa (16)
#define VARNAME_REGISTRY_s_LEAK_S (17)
#define VARNAME_REGISTRY_s_LEAK_D (18)
#define VARNAME_REGISTRY_lyapunov (19)
#define VARNAME_REGISTRY_m_Ca_slow (20)
#define VARNAME_REGISTRY_s_Ca_slow (21)

/* These definitions are for use with 
   GLOBAL_NEURON_MODEL==5 (Mainak's Code)
*/
#define VARNAME_REGISTRY_mainak_nAch_open (0)
#define VARNAME_REGISTRY_mainak_gaba_open (1)
#define VARNAME_REGISTRY_mainak_slow_inh_g (2)
#define VARNAME_REGISTRY_mainak_s_ORN (3)
#define VARNAME_REGISTRY_mainak_m_Na (4)
#define VARNAME_REGISTRY_mainak_h_Na (5)
#define VARNAME_REGISTRY_mainak_m_K (6)
#define VARNAME_REGISTRY_mainak_m_KA (7)
#define VARNAME_REGISTRY_mainak_h_KA (8)
#define VARNAME_REGISTRY_mainak_m_Ca (9)
#define VARNAME_REGISTRY_mainak_h_Ca (10)
#define VARNAME_REGISTRY_mainak_m_CaK (11)
#define VARNAME_REGISTRY_mainak_Ca (12)
#define VARNAME_REGISTRY_mainak_slow_inh_r (13)
#define VARNAME_REGISTRY_mainak_Vs (14)
#define VARNAME_REGISTRY_mainak_spike_flag (15)
#define VARNAME_REGISTRY_mainak_inputrate (16)
#define VARNAME_REGISTRY_mainak_nAch_open_local (17)
#define VARNAME_REGISTRY_mainak_gaba_open_local (18)
#define VARNAME_REGISTRY_mainak_slow_inh_r_local (19)

/* These definitions are for use with 
   GLOBAL_NEURON_MODEL==6 (Wilson Fly Code)
*/
#define VARNAME_REGISTRY_wilson_nAch (0)
#define VARNAME_REGISTRY_wilson_gabaA (1)
#define VARNAME_REGISTRY_wilson_gabaB (2)
#define VARNAME_REGISTRY_wilson_s_ORN (3)
#define VARNAME_REGISTRY_wilson_m_Na (4)
#define VARNAME_REGISTRY_wilson_h_Na (5)
#define VARNAME_REGISTRY_wilson_m_K (6)
#define VARNAME_REGISTRY_wilson_Vs (7)
#define VARNAME_REGISTRY_wilson_gabaBp (8)
#define VARNAME_REGISTRY_wilson_nAch_local (9)
#define VARNAME_REGISTRY_wilson_gabaA_local (10)
#define VARNAME_REGISTRY_wilson_gabaBp_local (11)
#define VARNAME_REGISTRY_wilson_inputrate (12)
#define VARNAME_REGISTRY_wilson_spike_flag (13)
#define VARNAME_REGISTRY_wilson_Vs_clip (14)
#define VARNAME_REGISTRY_wilson_vesicle_depletion (15)

/* These definitions are for use with
   GLOBAL_NEURON_MODEL==6 (Wilson Fly Code)
*/
#define VARNAME_REGISTRY_snx_nAch (0)
#define VARNAME_REGISTRY_snx_gabaA (1)
#define VARNAME_REGISTRY_snx_Vs (2)
#define VARNAME_REGISTRY_snx_spike_flag (3)

/* This holds a orn structure */
struct orn 
{
  double *rate;
  double rise1;
  double rise2;
  double fall1;
  double fall2;
  double normalizer;
};

/* This holds a pulse_synapse */
struct pulse_synapse
{
  struct llist *pos;
  struct llist *pre;
  struct llist **spiketime_llistra;
  struct llist **value_llistra;
};

/* This holds a neuron state */
struct neuron
{
  int ntypes;
  int type;
  int index;
  int nvars;
  int nsval;
  double **vpra;
  double *vra;
  double spiketime_guess;
  int spiketime_guess_flag;
  double inputrate;
  int spikeinput_flag;
  double spikeinput_multiplicity;
  long long int spikeinput_rseed;
  double spikeinput_time;
  double spikelast;
  double spiketime;
  double spikenext;
  int sparse_out;
  int sparse_in;
  struct llitem *sparse_link;
  double **dense_link;
  int sog;
  int microstep;
  struct pulse_synapse *pulse_synapse;
};

/* This holds multiple neurons */
struct neuronarray
{
  int ntypes;
  int nvars;
  int nsval;
  int *lengthra;
  int lt;
  double ***vrarara;
  struct neuron ***N; 
  struct clusterarray *gli;
  struct ptree *reweight;
};

/* This holds a cluster of neurons */
struct cluster
{
  int index;
  struct llitem *N;
  struct llist *LN;
};

/* This holds multiple clusters */
struct clusterarray
{
  struct neuronarray *Nra;
  int nclusters;
  struct cluster **cra;
  struct llist ***n2c;
  double **n2x;
  double **n2y;
};

/* here are cortex-dependent llists functions */
int spikeinput_time_compare(void *,void *);
int spikelast_compare(void *,void *);
int spiketime_compare(void *,void *);

/* Here are the orn functions */
struct orn * ornmake(double,double,double,double);
void orntfree(struct orn *);
void ornevolve(struct orn *,double,double);

/* Here are the neuron/neuronarray functions */
void neuronmake(struct neuronarray *,int,int);
void neurontfree(struct neuron *);
struct neuron * nget(struct neuronarray *,int,int);
void nset(struct neuronarray *,int,int,void *);
struct neuronarray * neuronarraymake(int,int *,int,int,int);
void neuronarraytfree(struct neuronarray *);
void setinputrate(struct neuronarray *,double);
void spikeinput(struct neuronarray *,double,double);
double spikeguess(struct neuron *,double,double);
int spikescan(struct neuronarray *,struct llist *,struct llist *,double,double);
int ilink(struct neuron *,double *);
int slink(struct neuron *,struct neuron *,double *,double *);
void slavecalc(int,double *,double *,double *,double);
void mainak_synapse_make(struct neuronarray *);
void mainak_synapse_tfree(struct neuronarray *);
void mainak_continuous_synapse_update(struct neuron *,int,double,double);
void mainak_continuous_synapse_evolve(struct neuron *,double,double);
void mainak_continuous_evolve_helper_rhs_j_inv_ln(double *,double *);
void mainak_continuous_evolve_helper_rhs_ln(double *,double *,double *);
void mainak_continuous_evolve_helper_rhs_j_inv_pn(double *,double *);
void mainak_continuous_evolve_helper_rhs_pn(double *,double *,double *);
void mainak_continuous_evolve(struct neuronarray *,double,double);
void wilson_synapse_make(struct neuronarray *);
void wilson_input_synapse_update(int,int,struct neuron *,double *,double);
void wilson_synapse_update(struct neuron *,double*,double);
void wilson_synapse_evolve_helper(int,int,int,struct llist **,struct llist **,double,double **,int **,int *);
void wilson_synapse_evolve(struct neuron *,double,double);
void wilson_evolve_helper_rhs_j_inv(double *,double *);
void wilson_evolve_helper_rhs(double *,double **,double *,double *,int);
void wilson_evolve(struct neuronarray *,double,double);
void snx_statematrix_0(int,int,int,int,double,double *);
void snx_statematrix_1(int,int,double,double *,double *);
void snx_evolve(struct neuronarray *,double,double);
void nllistevolve_cif_helper(struct llist *,struct neuron *,double,int *);
void nllistevolve_cif(struct llist *,struct llist *,double,double);
void vraevolve(struct neuron *,double,double);
void gizmointegrateeiforif(struct neuron *,double,double,double *);
void clumpcorrect(struct llist *,struct llist *,double,double);
void spikecorrect(struct llist *,double,double);
void spikeconduct(struct llist *,double,double,int *,int *);
double linerootfinder(double,double,double,double);
int gizmospiked(struct llist *,double,double);
void gizmoconduct(struct neuronarray *,struct llist *,struct llist *,double,double);

/* compte definitions */
/* #define TYPENAME_REGISTRY_PN (0) */
/* #define TYPENAME_REGISTRY_LN (1) */
/* #define VARNAME_REGISTRY_s_A (0) */
/* #define VARNAME_REGISTRY_s_N (1) */
/* #define VARNAME_REGISTRY_s_G (2) */
/* #define VARNAME_REGISTRY_s_ORN (3) */
/* #define VARNAME_REGISTRY_s_NMDA (4) */
/* #define VARNAME_REGISTRY_h_Na (5) */
/* #define VARNAME_REGISTRY_m_K (6) */
/* #define VARNAME_REGISTRY_h_KA (7) */
/* #define VARNAME_REGISTRY_m_KS (8) */
/* #define VARNAME_REGISTRY_Na (9) */
/* #define VARNAME_REGISTRY_Ca (10) */
/* #define VARNAME_REGISTRY_Vs (11) */
/* #define VARNAME_REGISTRY_Vd (12) */
/* #define VARNAME_REGISTRY_spike_flag (13) */
/* #define VARNAME_REGISTRY_inputrate (14) */
/* #define VARNAME_REGISTRY_s_Na (15) */
/* #define VARNAME_REGISTRY_s_K (16) */
/* #define VARNAME_REGISTRY_s_KA (17) */
/* #define VARNAME_REGISTRY_s_KS (18) */
/* #define VARNAME_REGISTRY_s_KNa (19) */
/* #define VARNAME_REGISTRY_s_Ca (20) */
/* #define VARNAME_REGISTRY_s_NaP (21) */
/* #define VARNAME_REGISTRY_s_AR (22) */
/* #define VARNAME_REGISTRY_s_KCa (23) */
/* #define VARNAME_REGISTRY_s_LEAK_S (24) */
/* #define VARNAME_REGISTRY_s_LEAK_D (25) */

/* pinsky definitions */
/* #define TYPENAME_REGISTRY_PN (0) */
/* #define TYPENAME_REGISTRY_LN (1) */
/* #define VARNAME_REGISTRY_s_A (0) */
/* #define VARNAME_REGISTRY_s_N (1) */
/* #define VARNAME_REGISTRY_s_G (2) */
/* #define VARNAME_REGISTRY_s_ORN (3) */
/* #define VARNAME_REGISTRY_s_NMDA (4) */
/* #define VARNAME_REGISTRY_h_Na (5) */
/* #define VARNAME_REGISTRY_m_K (6) */
/* #define VARNAME_REGISTRY_m_KAHP (7) */
/* #define VARNAME_REGISTRY_m_Ca (8) */
/* #define VARNAME_REGISTRY_m_KCa (9) */
/* #define VARNAME_REGISTRY_Ca (10) */
/* #define VARNAME_REGISTRY_Vs (11) */
/* #define VARNAME_REGISTRY_Vd (12) */
/* #define VARNAME_REGISTRY_spike_flag (13) */
/* #define VARNAME_REGISTRY_inputrate (14) */
/* #define VARNAME_REGISTRY_s_Na (15) */
/* #define VARNAME_REGISTRY_s_K (16) */
/* #define VARNAME_REGISTRY_s_KAHP (17) */
/* #define VARNAME_REGISTRY_s_Ca (18) */
/* #define VARNAME_REGISTRY_s_KCa (19) */
/* #define VARNAME_REGISTRY_s_LEAK_S (20) */
/* #define VARNAME_REGISTRY_s_LEAK_D (21) */

