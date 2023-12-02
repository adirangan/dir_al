/* here are cortex-dependent llists functions */

int vra_Vs_compare(void *item1,void *item2)
{
  /* if n1->spikeinput_time > n2->spikeinput_time */
  int comparison=0;
  struct neuron *n1=(struct neuron *) item1;
  struct neuron *n2=(struct neuron *) item2;
  if (n1->vra[VARNAME_REGISTRY_Vs] < n2->vra[VARNAME_REGISTRY_Vs]){ comparison = -1;}
  else if (n1->vra[VARNAME_REGISTRY_Vs] > n2->vra[VARNAME_REGISTRY_Vs]){ comparison = +1;}
  else if (n1->vra[VARNAME_REGISTRY_Vs] == n2->vra[VARNAME_REGISTRY_Vs] && n1 < n2 ){ comparison = -1;}
  else if (n1->vra[VARNAME_REGISTRY_Vs] == n2->vra[VARNAME_REGISTRY_Vs] && n1 > n2 ){ comparison = +1;}
  return comparison;
}

int spikeinput_time_compare(void *item1,void *item2)
{
  /* if n1->spikeinput_time > n2->spikeinput_time */
  int comparison=0;
  struct neuron *n1=(struct neuron *) item1;
  struct neuron *n2=(struct neuron *) item2;
  if (n1->spikeinput_time < n2->spikeinput_time){ comparison = -1;}
  else if (n1->spikeinput_time > n2->spikeinput_time){ comparison = +1;}
  else if (n1->spikeinput_time == n2->spikeinput_time && n1 < n2 ){ comparison = -1;}
  else if (n1->spikeinput_time == n2->spikeinput_time && n1 > n2 ){ comparison = +1;}
  return comparison;
}

int spikelast_compare(void *item1,void *item2)
{
  /* if n1->spikelast > n2->spikelast */
  int comparison=0;
  struct neuron *n1=(struct neuron *) item1;
  struct neuron *n2=(struct neuron *) item2;
  if (n1->spikelast < n2->spikelast){ comparison = -1;}
  else if (n1->spikelast > n2->spikelast){ comparison = +1;}
  else if (n1->spikelast == n2->spikelast && n1 < n2 ){ comparison = -1;}
  else if (n1->spikelast == n2->spikelast && n1 > n2 ){ comparison = +1;}
  return comparison;
}

int spiketime_compare(void *item1,void *item2)
{
  /* if n1->spiketime > n2->spiketime */
  int comparison=0;
  struct neuron *n1=(struct neuron *) item1;
  struct neuron *n2=(struct neuron *) item2;
  if (n1->spiketime < n2->spiketime){ comparison = -1;}
  else if (n1->spiketime > n2->spiketime){ comparison = +1;}
  else if (n1->spiketime == n2->spiketime && n1 < n2 ){ comparison = -1;}
  else if (n1->spiketime == n2->spiketime && n1 > n2 ){ comparison = +1;}
  return comparison;
}

/* Here are the neuron/neuronarray functions */

void neuronmake(struct neuronarray *Nra,int type,int index)
{
  int verbose=0;
  struct neuron *n=NULL;
  int nv=0;
  if (verbose){ printf(" %% [entering neuronmake] with type %d index %d\n",type,index);}
  n = (struct neuron *) tmalloc(sizeof(struct neuron)); 
  n->ntypes=Nra->ntypes;
  n->type=type;
  n->index=index;
  n->nvars=Nra->nvars;
  n->nsval=Nra->nsval;
  n->vpra = (double **) tcalloc(n->nvars,sizeof(double *));
  for (nv=0;nv<n->nvars;nv++){ n->vpra[nv] = &(Nra->vrarara[n->type][nv][n->index]);}
  n->vra = (double *) tcalloc(n->nvars,sizeof(double));
  n->spiketime_guess=0;
  n->spiketime_guess_flag=0;
  n->inputrate=0;
  n->spikeinput_flag=0;
  n->spikeinput_time=0;
  n->spikeinput_multiplicity=0;
  n->spikeinput_rseed= GLOBAL_SPIKEINPUT_RSEED + n->type + n->index*Nra->ntypes + GLOBAL_RECORD_NUMBER*Nra->lengthra[n->type]*Nra->ntypes;
  n->spikelast=0;
  n->spiketime=0;
  n->spikenext=0;
  n->sparse_out=0;
  n->sparse_in=0;
  switch (GLOBAL_NEURON_MODEL){
  case 0: case 1: case 2: case 3: case 4: 
    *(n->vpra[VARNAME_REGISTRY_Vs]) = VOLTAGE_[VARNAME_REGISTRY_Vs] + rand01*(VOLTAGE_THRESHOLD_S - VOLTAGE_[VARNAME_REGISTRY_Vs]);
    *(n->vpra[VARNAME_REGISTRY_Vd]) = *(n->vpra[VARNAME_REGISTRY_Vs]);
    *(n->vpra[VARNAME_REGISTRY_s_G]) = TAU_[VARNAME_REGISTRY_s_G]*rand01;
    *(n->vpra[VARNAME_REGISTRY_Ca]) = 0.01;
    break;
  case 5: /* mainak style */
    *(n->vpra[VARNAME_REGISTRY_mainak_Vs]) = -80;
    break;
  case 6: /* wilson style */
    *(n->vpra[VARNAME_REGISTRY_wilson_Vs]) = -80;
    break;
  case 7: /* snx style */
    *(n->vpra[VARNAME_REGISTRY_snx_nAch]) = 0;
    *(n->vpra[VARNAME_REGISTRY_snx_Vs]) = 0;
    break;
  default: break;}
  n->spikeinput_time=-2*TAU_REF;
  n->spikelast=-2*TAU_REF;
  n->spiketime=-2*TAU_REF;
  n->spikenext=-TAU_REF;
  n->sparse_out=rand();
  n->sparse_in=rand();
  n->sparse_link=NULL;
  n->dense_link=NULL;
  n->sog=0;
  n->microstep=1;
  n->pulse_synapse=NULL;
  nset(Nra,n->type,n->index,n);
}

void neurontfree(struct neuron *n)
{
  int nt=0;
  if (n->sparse_link!=NULL){ llitemtfree(n->sparse_link,NULL); n->sparse_link=NULL;}
  if (n->dense_link!=NULL){ for (nt=0;nt<GLOBAL_NTYPES;nt++){ tfree(n->dense_link[nt]);n->dense_link[nt]=NULL;}}
  tfree(n->dense_link);n->dense_link=NULL;
  if (n->vpra!=NULL){ tfree(n->vpra); n->vpra=NULL;}
  if (n->vra!=NULL){ tfree(n->vra); n->vra=NULL;}
  tfree(n); n = NULL;
}

struct neuron * nget(struct neuronarray *Nra,int type,int index)
{
  if (type<0 || type>=Nra->ntypes){ printf(" %% error! %d out of bounds %d for nget\n",type,Nra->ntypes);}
  if (index<0 || index>=Nra->lengthra[type]){ printf(" %% error! %d out of bounds %d for nget\n",index,Nra->lengthra[type]);}
  return Nra->N[type][index];
}

void nset(struct neuronarray *Nra,int type,int index,void *n)
{
  if (type<0 || type>=Nra->ntypes){ printf(" %% error! %d out of bounds %d for nset\n",type,Nra->ntypes);}
  if (index<0 || index>=Nra->lengthra[type]){ printf(" %% error! %d out of bounds %d for nset\n",index,Nra->lengthra[type]);}
  if (n==NULL){ Nra->N[type][index] = NULL;}
  else if (n!=NULL){ Nra->N[type][index] = (struct neuron *) n;}
}

struct cluster *clustermake(int index,struct neuron *n)
{
  struct cluster *c=NULL;
  c = (struct cluster *) tmalloc(sizeof(struct cluster));
  c->index = index;
  c->N = llitemmake(); c->N->item = n;
  c->LN = llistmake(); 
  return c;
}

void clustertfree(struct cluster *c){ llitemtfree(c->N,NULL);c->N=NULL;llisttfree(c->LN);c->LN=NULL;tfree(c);c=NULL;}

struct clusterarray *clusterarraymake(struct neuronarray *Nra,int nclusters)
{
  int verbose=0;
  int nc=0,nt=0,nr=0,nc2=0;
  double *pra=GLOBAL_CLUSTER_PRA_;
  struct clusterarray *C=NULL;
  struct neuron *n=NULL;
  struct llist *L=NULL;
  struct litem *l=NULL;
  struct cluster *c=NULL;
  double arc=0,arc2=0,arc3=0;
  if (verbose){ printf(" %% [entering clusterarraymake] nclusters %d\n",nclusters);}
  C = (struct clusterarray *) tmalloc(sizeof(struct clusterarray));
  C->Nra = Nra;
  C->nclusters = nclusters;
  C->cra = (struct cluster **) tcalloc(C->nclusters,sizeof(struct cluster *));
  for (nc=0;nc<C->nclusters;nc++){ 
    C->cra[nc] = clustermake(nc,NULL);
    for (nt=0;nt<C->Nra->ntypes;nt++){
      for (nr=nc*(C->Nra->lengthra[nt]/C->nclusters);nr<(nc+1)*(C->Nra->lengthra[nt]/C->nclusters);nr++){
	n = nget(C->Nra,nt,nr);
	if (C->cra[nc]->N->item==NULL){ C->cra[nc]->N->item=n;}
	else /* if (C->cra[nc]->N->item!=NULL) */{
	  if (llitemaddorfind(0,C->cra[nc]->N,n,&void_compare)==NULL){ 
	    if (verbose){ printf(" %% adding n_%d_(%d,%d) to cluster %d\n",(int)n,nt,nr,nc);}
	    llitemaddorfind(1,C->cra[nc]->N,n,&void_compare);}}}}
    llitembalance(C->cra[nc]->N); C->cra[nc]->N=llitemclimb(C->cra[nc]->N);
    llitem2llist(C->cra[nc]->N,C->cra[nc]->LN);
    if (verbose){ printf(" %% cluster %d:\n",C->cra[nc]->index); llitemprintf(C->cra[nc]->N,NULL);}}
  C->n2c = (struct llist ***) tcalloc(C->Nra->ntypes,sizeof(struct llist **));
  for (nt=0;nt<C->Nra->ntypes;nt++){
    C->n2c[nt] = (struct llist **) tcalloc(C->Nra->lengthra[nt],sizeof(struct llist *));
    for (nr=0;nr<C->Nra->lengthra[nt];nr++){
      C->n2c[nt][nr] = llistmake();
      n = nget(C->Nra,nt,nr);
      if (C->nclusters>0){
	nc = (nr*C->nclusters)/C->Nra->lengthra[nt];
	if (verbose){ printf(" %% %% adding cluster_%d_(%d) to C->n2c[%d][%d]\n",(int)C->cra[nc],nc,nt,nr);}
	litemadd(C->n2c[nt][nr],C->cra[nc]);
	for (nc2=0;nc2<C->nclusters;nc2++){
	  if (nc2!=nc && rand01<pra[nt]){ 
	    if (verbose){ printf(" %% %% adding cluster_%d_(%d) to C->n2c[%d][%d]\n",(int)C->cra[nc],nc2,nt,nr);}
	    litemadd(C->n2c[nt][nr],C->cra[nc2]);}}}}}
  arc = 2*PI/(double)C->nclusters;
  C->n2x = (double **) tcalloc(C->Nra->ntypes,sizeof(double *));
  C->n2y = (double **) tcalloc(C->Nra->ntypes,sizeof(double *));
  for (nt=0;nt<C->Nra->ntypes;nt++){
    C->n2x[nt] = (double *) tcalloc(C->Nra->lengthra[nt],sizeof(double));
    C->n2y[nt] = (double *) tcalloc(C->Nra->lengthra[nt],sizeof(double));
    arc2 = 2*PI/(double)(C->Nra->lengthra[nt]/maximum(1,C->nclusters));
    arc3 = 0.75*(1-(double)nt/(double)C->Nra->ntypes);
    for (nr=0;nr<C->Nra->lengthra[nt];nr++){
      n = nget(C->Nra,nt,nr);
      C->n2x[nt][nr]=0; C->n2y[nt][nr]=0;
      L = C->n2c[nt][nr]; l=L->first; 
      if (l!=NULL){ c=(struct cluster *)l->item; C->n2x[nt][nr] += cos(arc*(c->index+0.5)); C->n2y[nt][nr] += sin(arc*(c->index+0.5));}
      C->n2x[nt][nr] += arc3*cos(arc2*((n->index%(C->Nra->lengthra[nt]/maximum(1,C->nclusters))) + 0.5*n->type)); 
      C->n2y[nt][nr] += arc3*sin(arc2*((n->index%(C->Nra->lengthra[nt]/maximum(1,C->nclusters))) + 0.5*n->type));
      C->n2x[nt][nr] = (C->n2x[nt][nr]+1)/2;
      C->n2y[nt][nr] = (C->n2y[nt][nr]+1)/2;}}
  return C;
}

void clusterarraytfree(struct clusterarray *C)
{
  int nc=0,nt=0,nr=0;
  for (nt=0;nr<C->Nra->ntypes;nt++){ 
    for (nr=0;nr<C->Nra->lengthra[nt];nr++){
      llisttfree(C->n2c[nt][nr]); C->n2c[nt][nr]=NULL;}
    tfree(C->n2c[nt]); C->n2c[nt]=NULL;}
  tfree(C->n2c); C->n2c=NULL;
  for (nt=0;nr<C->Nra->ntypes;nt++){ 
    tfree(C->n2x[nt]); C->n2x[nt]=NULL;
    tfree(C->n2y[nt]); C->n2y[nt]=NULL;}
  tfree(C->n2x); C->n2x=NULL;
  tfree(C->n2y); C->n2y=NULL;
  for (nc=0;nc<C->nclusters;nc++){ clustertfree(C->cra[nc]); C->cra[nc]=NULL;}
  tfree(C->cra); C->cra=NULL;
  tfree(C);C=NULL;
}

struct neuronarray * neuronarraymake(int ntypes,int *lengthra,int nvars,int nsval,int nclusters)
{
  int verbose=0;
  int nr=0,nr2=0,nt=0,nt2=0,nv=0,tab=0;
  struct neuronarray *Nra=NULL;
  struct litem *l1=NULL,*l2=NULL;
  struct cluster *c=NULL;
  struct neuron *n=NULL,*n2=NULL;
  int sparsehit=0,sparse=0,link_flag=0;
  if (verbose){ printf(" %% [entering neuronarraymake] ntypes %d nvars %d nsval %d\n",ntypes,nvars,nsval); raprintf(lengthra,"int",1,ntypes,"lengthra: ");}
  Nra = (struct neuronarray *) tmalloc(sizeof(struct neuronarray));
  Nra->ntypes = ntypes;
  Nra->nvars = nvars;
  Nra->nsval = nsval;
  Nra->lengthra = (int *) tcalloc(Nra->ntypes,sizeof(int));
  Nra->vrarara = (double ***) tcalloc(Nra->ntypes,sizeof(double **));
  Nra->N = (struct neuron ***) tcalloc(Nra->ntypes,sizeof(struct neuron **));
  Nra->lt=0;
  for (nt=0;nt<Nra->ntypes;nt++){
    Nra->lengthra[nt] = lengthra[nt]; Nra->lt += Nra->lengthra[nt];
    Nra->vrarara[nt] = (double **) tcalloc(Nra->nvars,sizeof(double *));
    Nra->N[nt] = (struct neuron **) tcalloc(Nra->lengthra[nt],sizeof(struct neuron *));
    for (nv=0;nv<Nra->nvars;nv++){
      Nra->vrarara[nt][nv] = (double *) tcalloc(Nra->lengthra[nt],sizeof(double));}
    for (nr=0;nr<Nra->lengthra[nt];nr++){ 
      neuronmake(Nra,nt,nr);}}
  if (verbose){ raprintf(Nra->lengthra,"int",1,Nra->ntypes,"lengthra ");}
  Nra->gli=clusterarraymake(Nra,nclusters);
  if (GLOBAL_LINK_SPARSE_OR_DENSE==1){
    if (GLOBAL_CLUSTER_MAKE_OR_READ==+1){
      for (nt=0;nt<Nra->ntypes;nt++){ for (nr=0;nr<Nra->lengthra[nt];nr++){
	n = nget(Nra,nt,nr);
	if (verbose){ printf(" %% considering n_%d_(%d,%d)\n",(int)n,nt,nr);}
	n->sparse_link = llitemmake(); n->sparse_link->item=n;
	l1 = Nra->gli->n2c[nt][nr]->first; 
	while (l1!=NULL){
	  c = (struct cluster *)l1->item;
	  if (verbose){ printf(" %% %% linking to cluster %d\n",c->index);}
	  l2 = c->LN->first;
	  while (l2!=NULL){
	    n2 = (struct neuron *)l2->item;
	    if (verbose){ printf(" %% %% %% n_%d_(%d,%d)...",(int)n2,n2->type,n2->index);}
	    if (l1==Nra->gli->n2c[nt][nr]->first){
	      sparsehit = SPARSEHIT__[n->type + n2->type*n2->ntypes]; sparse = SPARSE__[n->type + n2->type*n2->ntypes];}
	    else /* if (l1!=Nra->gli->n2c[nt][nr]->first) */{
	      sparsehit = SPARSEHIT_OGLI__[n->type + n2->type*n2->ntypes]; sparse = SPARSE_OGLI__[n->type + n2->type*n2->ntypes];}
	    if (periodize((n->sparse_out-n2->sparse_in)%sparsehit,0,sparsehit)<sparse){
	      if (verbose){ printf(" linking\n");}
	      if (llitemaddorfind(0,n->sparse_link,n2,&void_compare)==NULL){ llitemaddorfind(1,n->sparse_link,n2,&void_compare);}}
	    else /* if (periodize((n->sparse_out-n2->sparse_in)%sparsehit,0,sparsehit)>=sparse) */{ if (verbose){ printf("\n");}}
	    l2=l2->child;}
	  l1=l1->child;}
	llitembalance(n->sparse_link); n->sparse_link=llitemclimb(n->sparse_link);}}}
    else if (GLOBAL_CLUSTER_MAKE_OR_READ==-1){ 
      sparse_link_read_ascii(Nra,strcmp(GLOBAL_CLUSTER_READ_FILE,"default")==0 ? NULL : GLOBAL_CLUSTER_READ_FILE);}}
  else if (GLOBAL_LINK_SPARSE_OR_DENSE==-1){
    for (nt2=0;nt2<Nra->ntypes;nt2++){ for (nr2=0;nr2<Nra->lengthra[nt2];nr2++){
      n=nget(Nra,nt2,nr2);
      n->dense_link = (double **) tcalloc(Nra->ntypes,sizeof(double *));
      for (nt=0;nt<Nra->ntypes;nt++){ 
	n->dense_link[nt] = (double *) tcalloc(Nra->lengthra[nt]*GLOBAL_INDEXING_sra_LENGTH,sizeof(double));
	for (nr=0;nr<Nra->lengthra[nt];nr++){
	  n2 = nget(Nra,nt,nr);
	  sparsehit = SPARSEHIT__[n->type + n2->type*n2->ntypes]; sparse = SPARSE__[n->type + n2->type*n2->ntypes];
	  if (periodize((n->sparse_out-n2->sparse_in)%sparsehit,0,sparsehit)<sparse){ link_flag=1;} else{ link_flag=0;}
	  for (nv=0;nv<GLOBAL_INDEXING_sra_LENGTH;nv++){
	    if (verbose){ printf(" considering nt %d nr %d nv %d...",nt,nr,nv);}
	    tab = n->type + n2->type*n->ntypes + nv*n->ntypes*n2->ntypes;
	    n->dense_link[nt][nr+nv*Nra->lengthra[nt]] = (0.25+(0.25+rand01)*link_flag)*CS__[tab];
	    if (verbose){ printf(" setting strength %f\n",n->dense_link[nt][nr+nv*Nra->lengthra[nt]]);}
	    if (AUTAPSES_OFF && n2==n){ n->dense_link[nt][nr+nv*Nra->lengthra[nt]] = 0;}}}}}}}
  if (strcmp(GLOBAL_REWEIGHT_READ_FILE,"default")!=0){ Nra->reweight = ptreadback(GLOBAL_REWEIGHT_READ_FILE,0);}
  else /* if (strcmp(GLOBAL_REWEIGHT_READ_FILE,"default")==0) */{ Nra->reweight = NULL;}
  if (GLOBAL_NEURON_MODEL==5){ mainak_synapse_make(Nra);}
  if (GLOBAL_NEURON_MODEL==6){ wilson_synapse_make(Nra);}
  if (verbose){ printf(" %% [leaving neuronarraymake]\n");}
  return Nra;
}

void neuronarraytfree(struct neuronarray *Nra)
{
  int nr=0,nt=0,nv=0;
  clusterarraytfree(Nra->gli);
  for (nt=0;nt<Nra->ntypes;nt++){
    for (nr=0;nr<Nra->lengthra[nt];nr++){ neurontfree(nget(Nra,nt,nr)); nset(Nra,nt,nr,NULL);}
    for (nv=0;nv<Nra->nvars;nv++){ tfree(Nra->vrarara[nt][nv]); Nra->vrarara[nt][nv]=NULL;}
    tfree(Nra->N[nt]);Nra->N[nt]=NULL;
    tfree(Nra->vrarara[nt]);Nra->vrarara[nt]=NULL;}
  tfree(Nra->N);Nra->N=NULL;
  tfree(Nra->vrarara);Nra->vrarara=NULL;
  tfree(Nra->lengthra);Nra->lengthra=NULL;
  tfree(Nra);Nra=NULL;
}

void setinputrate(struct neuronarray *Nra,double t)
{
  /* sets input spiketimes */
  int nr=0,nt=0;
  double rateplus=0;
  struct neuron *n=NULL;
  double local_odor=0;
  for (nt=0;nt<Nra->ntypes;nt++){ for (nr=0;nr<Nra->lengthra[nt];nr++){
    n=nget(Nra,nt,nr);
    switch (GLOBAL_NEURON_MODEL){
    case 0: case 1: case 2: case 3: case 4:
      local_odor = (double)GLOBAL_ODORra->rara[nt][nr]/(double)GLOBAL_ODORra->base;
      rateplus = local_odor;
      rateplus *= INPUT_CONTRAST;
      rateplus *= (ORN_BOTHER ? GLOBAL_ORN->rate[4] : INPUT_PULSE);
      n->inputrate = maximum(0,ORN_BACKRATE + rateplus);
      break;
    case 5: /* mainak */
      local_odor = (double)GLOBAL_ODORra->rara[nt][nr]/(double)GLOBAL_ODORra->base;
      rateplus = local_odor;
      rateplus *= INPUT_CONTRAST;
      rateplus *= (ORN_BOTHER ? GLOBAL_ORN->rate[4] : INPUT_PULSE);
      n->inputrate = maximum(0,ORN_BACKRATE + rateplus);
      break;
    default: break;}}}
}

void spikeinput(struct neuronarray *Nra,double t,double DT)
{
  /* sets input spiketimes */
  int verbose=0;
  int nr=0,nt=0;
  struct neuron *n=NULL;
  double d=0,temp=0,threshold_time=1024;
  for (nt=0;nt<Nra->ntypes;nt++){ for (nr=0;nr<Nra->lengthra[nt];nr++){ 
    n = nget(Nra,nt,nr);
    switch (GLOBAL_NEURON_MODEL){
    case 0: case 1: case 2: case 3: case 4:
      if (n->inputrate>0){ 
	if (n->spikeinput_time < t){
	  if (nr==0 && nt==0 && verbose){ printf("t=%0.2f DT=%0.2f n(%d)->spikeinput_time=%0.2f, inputrate=%0.2f\n",t,DT,nr,n->spikeinput_time,n->inputrate);}
	  n->spikeinput_multiplicity=0; d=0;
	  while (n->spikeinput_time < t){ 
	    /* fix later! */
	    d = RISIGET(&(n->spikeinput_rseed),n->inputrate); //printf("%f,",d);
	    //d = (0.33+0.003*n->index)/(double)n->inputrate; //printf("%f\n",d);
	    if (d<threshold_time){ n->spikeinput_time +=d; n->spikeinput_multiplicity += 1;}
	    else{ n->spikeinput_time += 1;}
	    if (nr==0 && nt==0 && verbose){ printf("\tinpttm=%0.2f, mltplcty=%f, ",n->spikeinput_time,n->spikeinput_multiplicity);}}}
	if (n->spikeinput_time >= t){
	  if (n->spikeinput_time <= t+DT && n->spikeinput_multiplicity>0){ n->spikeinput_flag=1;}
	  else{ n->spikeinput_flag=0;}
	  if (nr==0 && nt==0 && verbose){ printf("spikeinput_flag %d\n",n->spikeinput_flag);}}}
      else{ n->spikeinput_flag=0;}
      break;
    case 5: /* mainak */
      *(n->vpra[VARNAME_REGISTRY_mainak_s_ORN])+=(rand01<ORN_BACKRATE*DT)*GLOBAL_CS_ORN_SCALE*CS_ORN_[n->type];
      d = INPUT_CONTRAST*(double)GLOBAL_ODORra->rara[nt][nr]/(double)GLOBAL_ODORra->base;
      d *= (ORN_BOTHER ? GLOBAL_ORN->rate[4] : INPUT_PULSE)*0.03572/*ORN_INPUT_RATE*/; 
      d = maximum(0,minimum(1,d*DT));
      temp = (maximum(0,randn()*sqrt(d*(1-d)*200/* NUMORN */)) + d*200/* NUMORN */);
      temp *= GLOBAL_CS_ORN_SCALE*CS_ORN_mainak_stim_[n->type];
      *(n->vpra[VARNAME_REGISTRY_mainak_s_ORN]) += temp;
      break;
    case 6: /* wilson */
      if (CS_ORN_[n->type]!=0 || CS_ORN_wilson_stim_[n->type]!=0){
	*(n->vpra[VARNAME_REGISTRY_wilson_s_ORN])+=(R01GET(&(n->spikeinput_rseed))<ORN_BACKRATE*DT)*GLOBAL_CS_ORN_SCALE*CS_ORN_[n->type];
	d = INPUT_CONTRAST*(double)GLOBAL_ODORra->rara[nt][nr]/(double)GLOBAL_ODORra->base;
	d *= (ORN_BOTHER ? GLOBAL_ORN->rate[4] : INPUT_PULSE)*0.03572/*INPUT_RATE_to_ORNs*/; 
	d = maximum(0,minimum(1,d*DT));
	temp = (maximum(0,RNGET(&(n->spikeinput_rseed))*sqrt(d*(1-d)*50/* NUM_to_ORNs */)) + d*50/* NUM_to_ORNs */);
	temp *= GLOBAL_CS_ORN_SCALE*CS_ORN_wilson_stim_[n->type];
	*(n->vpra[VARNAME_REGISTRY_wilson_s_ORN]) += temp;}
      else /* if (CS_ORN_[n->type]==0 && CS_ORN_wilson_stim[n->type]==0) */{ /* do nothing */}
      break;
    default: break;}}}
}

double spikeguess(struct neuron *n,double t,double DT)
{
  /* use clumpcorrect (with test excitation), fix later*/
  struct llist *L1=llistmake(),*L2=llistmake();
  litemadd(L1,n);
  clumpcorrect(L1,L2,t,DT);
  n->spiketime = n->spiketime_guess_flag ? n->spiketime_guess : t+2*DT;
  llisttfree(L1);llisttfree(L2);
  //if (n->spiketime_guess_flag && n->spiketime_guess >=t && n->spiketime_guess <= t+DT){ printf("spike at %f\n",n->spiketime_guess);}
  return n->spiketime;
}

int spikescan(struct neuronarray *Nra,struct llist *LS,struct llist *LG,double t,double DT)
{
  int nr=0,nt=0,pspikes=0;
  double spiketime_guess=0;
  struct neuron *n=NULL;
  for (nt=0;nt<Nra->ntypes;nt++){ for (nr=0;nr<Nra->lengthra[nt];nr++){
    n = nget(Nra,nt,nr);
    spiketime_guess = spikeguess(n,t,DT);
    if (spiketime_guess >= t && spiketime_guess <= t+DT){ n->sog=1; litemadd(LS,n); pspikes++;}
    else if (spiketime_guess < t || spiketime_guess > t+DT){ n->sog=-1; litemadd(LG,n);}}}
  return pspikes;
}

int ilink(struct neuron *n,double *sra){ 
  int nv=0;
  if (n!=NULL){ if (sra!=NULL){
    if (GLOBAL_NEURON_MODEL==0 || GLOBAL_NEURON_MODEL==1 || GLOBAL_NEURON_MODEL==2 || GLOBAL_NEURON_MODEL==3){
      sra[GLOBAL_INDEXING_REFILE_sra[VARNAME_REGISTRY_s_ORN]] = n->spikeinput_multiplicity*GLOBAL_CS_ORN_SCALE*CS_ORN_[n->type]; 
      return 1;}
    else if (GLOBAL_NEURON_MODEL==4){
      sra[GLOBAL_INDEXING_REFILE_sra[VARNAME_REGISTRY_Vs]] = n->spikeinput_multiplicity*GLOBAL_CS_ORN_SCALE*CS_ORN_[n->type];
      return 1;}}}
  else /* if (n==NULL) */{
    for (nv=0;nv<GLOBAL_INDEXING_sra_LENGTH;nv++){ if (sra!=NULL){ sra[nv]=0;}} return 0;}
  return 0;
}

int slink(struct neuron *s,struct neuron *n,double *sra,double *maxsra)
{
  /* returns synaptic links for s and n 
     Note that with the addition of n->dense_link this function is no longer 'clean', in that it references GLOBAL_Nra */
  int connected=0;
  int nv=0,sparsehit=0,sparse=0;
  int tab=s->type+n->type*n->ntypes,tab2=0;
  double multiplier=1,adder=0,sra_temp=0;
  ilink(NULL,sra);ilink(NULL,maxsra);
  if (AUTAPSES_OFF[s->type] && s==n){ /* do nothing */ }
  else{
    if (GLOBAL_LINK_SPARSE_OR_DENSE==1){
      if (s->sparse_link==NULL){
	sparsehit = SPARSEHIT__[tab]; sparse = SPARSE__[tab];
	if (periodize((s->sparse_out-n->sparse_in)%sparsehit,0,sparsehit)<sparse){ connected=1;}}
      else /* if (s->sparse_link!=NULL) */{ 
	if (llitemaddorfind(0,s->sparse_link,(void *)n,&void_compare)!=NULL){ connected=1;}}
      if (connected && (sra!=NULL || maxsra!=NULL)){ 
	for (nv=0;nv<GLOBAL_INDEXING_sra_LENGTH;nv++){
	  tab2 = s->type + n->type*s->ntypes + nv*n->ntypes*s->ntypes;
	  multiplier = GLOBAL_CS_SCALE*(GLOBAL_CS_PRETYPE_SCALE_==NULL ? 1 : GLOBAL_CS_PRETYPE_SCALE_[s->type])*(GLOBAL_CS_POSTYPE_SCALE_==NULL ? 1 : GLOBAL_CS_POSTYPE_SCALE_[n->type])*(GLOBAL_CS_SRA_SCALE_==NULL ? 1 : GLOBAL_CS_SRA_SCALE_[nv]);
	  if (GLOBAL_Nra->reweight!=NULL && GLOBAL_REWEIGHT_STRENGTH!=0){
	    if (GLOBAL_REWEIGHT_MULTIPLY_OR_ADD==1){
	      multiplier *= reweight2multiplier(GLOBAL_Nra->reweight,GLOBAL_Nra,s,n,GLOBAL_REWEIGHT_STRENGTH);
	      adder = 0;}
	    else if (GLOBAL_REWEIGHT_MULTIPLY_OR_ADD==-1){
	      multiplier *= 1;
	      adder = reweight2adder(GLOBAL_Nra->reweight,GLOBAL_Nra,s,n,GLOBAL_REWEIGHT_STRENGTH);}}
	  if (CS__[tab2]==0){ sra_temp=0;}
	  else if (CS__[tab2]!=0){ sra_temp=maximum(0,multiplier*CS__[tab2]+adder);}
	  if (P_FAIL__[tab2]==1){ if (sra!=NULL){ sra[nv] = sra_temp;}} 
	  else{ if (sra!=NULL){ sra[nv] = (rand01<P_FAIL__[tab2])*sra_temp/P_FAIL__[tab2];}}
	  if (maxsra!=NULL){ maxsra[nv] = sra_temp/P_FAIL__[tab2];}}}}
    else if (GLOBAL_LINK_SPARSE_OR_DENSE==-1){
      connected = 1;
      for (nv=0;nv<GLOBAL_INDEXING_sra_LENGTH;nv++){
	tab2 = s->type + n->type*s->ntypes + nv*n->ntypes*s->ntypes;
	multiplier = GLOBAL_CS_SCALE*(GLOBAL_CS_PRETYPE_SCALE_==NULL ? 1 : GLOBAL_CS_PRETYPE_SCALE_[s->type])*(GLOBAL_CS_POSTYPE_SCALE_==NULL ? 1 : GLOBAL_CS_POSTYPE_SCALE_[n->type])*(GLOBAL_CS_SRA_SCALE_==NULL ? 1 : GLOBAL_CS_SRA_SCALE_[nv]);
	if (GLOBAL_Nra->reweight!=NULL && GLOBAL_REWEIGHT_STRENGTH!=0){
	  if (GLOBAL_REWEIGHT_MULTIPLY_OR_ADD==1){
	    multiplier *= reweight2multiplier(GLOBAL_Nra->reweight,GLOBAL_Nra,s,n,GLOBAL_REWEIGHT_STRENGTH);
	    adder = 0;}
	  else if (GLOBAL_REWEIGHT_MULTIPLY_OR_ADD==-1){
	    multiplier *= 1;
	    adder = reweight2adder(GLOBAL_Nra->reweight,GLOBAL_Nra,s,n,GLOBAL_REWEIGHT_STRENGTH);}}
	sra_temp = s->dense_link[n->type][n->index+nv*GLOBAL_Nra->lengthra[n->type]];
	if (sra_temp==0){ /* do nothing */ }
	else if (sra_temp!=0){ sra_temp = maximum(0,multiplier*sra_temp+adder);}
	if (P_FAIL__[tab2]==1){ if (sra!=NULL){ sra[nv] = sra_temp;}} 
	else{ if (sra!=NULL){ sra[nv] = (rand01<P_FAIL__[tab2])*sra_temp/P_FAIL__[tab2];}}
	if (maxsra!=NULL){ maxsra[nv] = sra_temp/P_FAIL__[tab2];}}}}
  return connected;
}

void threshmaker(double x,double a,double b,double x_flat,double *out,double *outprime)
{
  /* min(ax,b)
   *out = b*(atan((2*x*a/b - 1)*x_flat)/PI + 1/2);
   *outprime = (1/(1+pow((2*x*a/b - 1)*x_flat,2)))*2*x_flat*a/PI; */
  double temp1 = (2*x*a/b - 1)*x_flat;
  if (out!=NULL){ *out = b*(atan(temp1)/PI+1/2);}
  if (outprime!=NULL){ *outprime = 2*x_flat*a/PI/(1+pow(temp1,2));}
}

void hump(double stepat,double height,double width,double x,double *out,double *outprime)
{
  double temp1 = 4*(x-stepat)/width;
  if (out!=NULL){ *out=0.5*height*(1+erf(temp1));}
  if (outprime!=NULL){ *outprime=0.5*height*exp(-pow(temp1,2))*4/width*2/sqrt(PI);}
}

void sra_add(struct neuron *n,double *sra)
{
  int verbose=0;
  int nv=0;
  if (verbose){ printf(" %% [entering sra_add] with:"); raprintf(sra,"double",1,GLOBAL_INDEXING_sra_LENGTH," sra ");}
  if (verbose){ for (nv=0;nv<GLOBAL_INDEXING_sra_LENGTH;nv++){ printf("start %f\n",n->vra[GLOBAL_INDEXING_CHECKOUT_sra[nv]]);}}
  for (nv=0;nv<GLOBAL_INDEXING_sra_LENGTH;nv++){ n->vra[GLOBAL_INDEXING_CHECKOUT_sra[nv]] += sra[nv];}
  if (verbose){ for (nv=0;nv<GLOBAL_INDEXING_sra_LENGTH;nv++){ printf("end %f\n",n->vra[GLOBAL_INDEXING_CHECKOUT_sra[nv]]);}}
}

struct orn * ornmake(double rise1,double rise2,double fall1,double fall2)
{
  int verbose=0;
  struct orn *p=NULL;
  double templ=0,tempr=0,tempc=0,tempf=0;
  if (verbose){ printf(" %% [entering ornmake] rise1 %f rise2 %f fall1 %f fall2 %f\n",rise1,rise2,fall1,fall2);}
  p = (struct orn *) tcalloc(1,sizeof(struct orn));
  p->rate = (double *) tcalloc(5,sizeof(double));
  p->rise1 = rise1; p->rise2 = rise2;
  p->fall1 = fall1; p->fall2 = fall2;
  templ = minimum(rise1,minimum(rise2,minimum(fall1,fall2))); tempr = maximum(rise1,maximum(rise2,maximum(fall1,fall2)));
  while ((((exp(-templ/rise2)-exp(-templ/rise1))/(rise2-rise1) - (exp(-templ/fall2)-exp(-templ/fall1))/(fall2-fall1)))<0){ templ/=2;}
  while ((((exp(-tempr/rise2)-exp(-tempr/rise1))/(rise2-rise1) - (exp(-tempr/fall2)-exp(-tempr/fall1))/(fall2-fall1)))>0){ tempr*=2;}
  do{
    tempc = (templ+tempr)/2;
    tempf = (exp(-tempc/rise2)-exp(-tempc/rise1))/(rise2-rise1) - (exp(-tempc/fall2)-exp(-tempc/fall1))/(fall2-fall1);
    if (verbose){ printf(" templ %f, tempc %f, tempr %f, tempf %f\n",templ,tempc,tempr,tempf);}
    if (tempf>=0){ templ=tempc;} else /* if (tempf<0) */{ tempr=tempc;}}
  while (fabs(tempr-templ)>0.00000001);
  if (verbose){ printf(" f0 found at %f\n",tempc);}
  p->normalizer = ((rise2*(1-exp(-tempc/rise2))-rise1*(1-exp(-tempc/rise1)))/(rise2-rise1) - (fall2*(1-exp(-tempc/fall2))-fall1*(1-exp(-tempc/fall1)))/(fall2-fall1));
  if (verbose){ printf(" normalizer %f\n",p->normalizer);}
  return p;
}

void orntfree(struct orn *p){ tfree(p->rate);p->rate=NULL;tfree(p);p=NULL;}

void ornevolve(struct orn *p,double t,double DT)
{
  /* this evolves the orn structure by a timestep DT */
  /* if linear_flag==1 --> rise1=50;rise2=75; t=0:.1:300; x=(exp(-t/rise1)-exp(-t/rise2))/(1/rise2-1/rise1)/rise2/rise1; fall1=125;fall2=150; y=(exp(-t/fall1)-exp(-t/fall2))/(1/fall2-1/fall1)/fall2/fall1; plot(t,x,t,y,t,x-y); */
  /* if linear_flag==0 --> */
  int verbose=0;
  int linear_flag=0;
  double e1=0,e2=0,e3=0,e4=0,tmp=0;
  if (verbose){ printf(" %% [entering ornevolve] with t=%f,DT=%f\n",t,DT);}
  switch (linear_flag){
  case 1:
    e1=exp(-DT/p->rise1);e2=exp(-DT/p->rise2);e3=exp(-DT/p->fall1);e4=exp(-DT/p->fall2);
    tmp = INPUT_PULSE;
    p->rate[1] = e2*p->rate[1] + 1.0/(p->rise1-p->rise2)*(e1-e2)*(p->rate[0] - p->rise1*tmp) + (1-e2)*tmp;
    p->rate[0] = e1*p->rate[0] + p->rise1*(1-e1)*tmp;
    p->rate[3] = e4*p->rate[3] + 1.0/(p->fall1-p->fall2)*(e3-e4)*(p->rate[2] - p->fall1*tmp) + (1-e4)*tmp;
    p->rate[2] = e3*p->rate[2] + p->fall1*(1-e3)*tmp;
    p->rate[4] = (p->rate[1] - p->rate[3])/p->normalizer;
    break;
  case 0: 
    tmp = INPUT_PULSE;
    if (tmp>p->rate[4]){ p->rate[4] = tmp - exp(-DT/p->rise1)*(tmp-p->rate[4]);}
    else if (tmp<p->rate[4]){ p->rate[4] = tmp + exp(-DT/p->fall2)*(p->rate[4]-tmp);}
    break;
  default: break;}
}

void mainak_synapse_make(struct neuronarray *Nra)
{
  int verbose=0;
  int nt=0,nt2=0,nr=0,nr2=0,tab=0,tab_max=3;
  struct neuron *n=NULL,*s=NULL;
  if (verbose){ printf(" %% [entering mainak_synapse_make]\n");}
  for (nt=0;nt<Nra->ntypes;nt++){ for (nr=0;nr<Nra->lengthra[nt];nr++){
    n = nget(Nra,nt,nr);
    if (n->pulse_synapse!=NULL){ printf(" warning! mainak_synapse!=NULL in mainak_synapse_make\n");}
    n->pulse_synapse = (struct pulse_synapse *) tcalloc(1,sizeof(struct pulse_synapse));
    n->pulse_synapse->pre = llistmake();
    n->pulse_synapse->pos = llistmake();
    n->pulse_synapse->spiketime_llistra = (struct llist **) tcalloc(tab_max,sizeof(struct llist *));
    n->pulse_synapse->value_llistra = (struct llist **) tcalloc(tab_max,sizeof(struct llist *));
    for (tab=0;tab<tab_max;tab++){
      n->pulse_synapse->spiketime_llistra[tab] = llistmake();
      n->pulse_synapse->value_llistra[tab] = llistmake();}}}
  for (nt=0;nt<Nra->ntypes;nt++){ for (nr=0;nr<Nra->lengthra[nt];nr++){
    s = nget(Nra,nt,nr);
    for (nt2=0;nt2<Nra->ntypes;nt2++){ for (nr2=0;nr2<Nra->lengthra[nt2];nr2++){
      n = nget(Nra,nt2,nr2);
      if (verbose){ printf(" %% checking connection between (%d,%d)->(%d,%d)\n",s->type,s->index,n->type,n->index);}
      if (slink(s,n,NULL,NULL)){ 
	if (verbose){ printf(" %% connected\n");}
	litemadd(s->pulse_synapse->pos,n);
	litemadd(n->pulse_synapse->pre,s);}
      else{ if (verbose){ printf(" %% not connected\n");}}}}}}
  if (verbose){ printf(" %% [finishing mainak_synapse_make]\n");}
}

void mainak_synapse_tfree(struct neuronarray *Nra)
{
  int nt=0,nr=0,tab=0,tab_max=3;
  struct neuron *n=NULL;
  for (nt=0;nt<Nra->ntypes;nt++){ for (nr=0;nr<Nra->lengthra[nt];nr++){
    n=nget(Nra,nt,nr);
    llisttfree(n->pulse_synapse->pre); n->pulse_synapse->pre=NULL; llisttfree(n->pulse_synapse->pos); n->pulse_synapse->pos=NULL;
    for (tab=0;tab<tab_max;tab++){ 
      llisttfree2(n->pulse_synapse->spiketime_llistra[tab]); llisttfree2(n->pulse_synapse->value_llistra[tab]);}
    tfree(n->pulse_synapse->spiketime_llistra);n->pulse_synapse->spiketime_llistra=NULL;
    tfree(n->pulse_synapse->value_llistra);n->pulse_synapse->value_llistra=NULL;
    tfree(n->pulse_synapse); n->pulse_synapse=NULL;}}
}

void mainak_continuous_synapse_update(struct neuron *s,int tab,double value,double spiketime)
{
  int verbose=0;//(s->type==0 && s->index==1);
  struct litem *l=NULL,*l_spiketime=NULL,*l_value=NULL,*l3=NULL;
  struct neuron *n=NULL;
  double *temp=NULL;
  if (verbose){ printf(" %% [entering mainak_continuous_synapse_update] with n(%d,%d),tab %d value %f spiketime %f\n",s->type,s->index,tab,value,spiketime);}
  l = s->pulse_synapse->pos->first;
  while (l!=NULL){
    n = (struct neuron *) l->item;
    if ((s->type==0 && tab==0) || (s->type==1 && tab==2) || (s->type==1 && n->type==0 && tab==1)){
      if (verbose>2){
	printf(" %% starting neuron (%d,%d):\n",n->type,n->index);
	llistprintf2(n->pulse_synapse->spiketime_llistra[tab]);
	llistprintf2(n->pulse_synapse->value_llistra[tab]);}
      l_spiketime = n->pulse_synapse->spiketime_llistra[tab]->last;
      l_value = n->pulse_synapse->value_llistra[tab]->last;
      while (l_spiketime!=NULL && *(double *)l_spiketime->item > spiketime){ l_spiketime = l_spiketime->parent; l_value = l_value->parent;}
      if (l_spiketime==NULL && l_value==NULL){
	temp = (double *) tcalloc(1,sizeof(double)); *temp = value;
	if (n->pulse_synapse->value_llistra[tab]->length==0){ litemadd(n->pulse_synapse->value_llistra[tab],temp);} 
	else if (n->pulse_synapse->value_llistra[tab]->length>0){ l3 = liteminsertbefore(n->pulse_synapse->value_llistra[tab],n->pulse_synapse->value_llistra[tab]->first); l3->item = temp;}
	temp = (double *) tcalloc(1,sizeof(double)); *temp = spiketime;
	if (n->pulse_synapse->spiketime_llistra[tab]->length==0){ litemadd(n->pulse_synapse->spiketime_llistra[tab],temp);} 
	else if (n->pulse_synapse->spiketime_llistra[tab]->length>0){ l3 = liteminsertbefore(n->pulse_synapse->spiketime_llistra[tab],n->pulse_synapse->spiketime_llistra[tab]->first); l3->item = temp;}}
      else if (l_spiketime!=NULL && l_value!=NULL){
	temp = (double *) tcalloc(1,sizeof(double)); *temp = value;
	l3 = liteminsertafter(n->pulse_synapse->value_llistra[tab],l_value);
	l3->item = temp;	
	temp = (double *) tcalloc(1,sizeof(double)); *temp = spiketime;
	l3 = liteminsertafter(n->pulse_synapse->spiketime_llistra[tab],l_spiketime);
	l3->item = temp;}
      else if ((l_spiketime!=NULL && l_value==NULL) || (l_spiketime==NULL && l_value!=NULL)){ printf(" %% warning! spiketime_llistra[%d] and value_llistra[%d] of different length in mainak_continuous_synapse_update\n",tab,tab);}
      if (verbose>2){
	printf(" %% finishing neuron (%d,%d):\n",n->type,n->index);
	llistprintf2(n->pulse_synapse->spiketime_llistra[tab]);
	llistprintf2(n->pulse_synapse->value_llistra[tab]);}}
    else{ /* do nothing */}
    l=l->child;}
  if (verbose){ printf(" %% [finishing mainak_continuous_synapse_update]\n");}
}

void mainak_continuous_synapse_evolve(struct neuron *n,double t,double DT)
{
  int verbose=0;//(n->type==0 && n->index==0);
  int tab=0;
  double temp_time=0,temp_spiketime=0;
  struct litem *l_spiketime=NULL,*l_value=NULL;
  if (verbose){ printf(" %% [entering mainak_continuous_synapse_evolve] with n(%d,%d),time %f->%f\n",n->type,n->index,t,t+DT);}
  tab=0;
  if (verbose>2){ 
    printf(" %% examining nAch llist\n"); 
    llistprintf2(n->pulse_synapse->spiketime_llistra[tab]);
    llistprintf2(n->pulse_synapse->value_llistra[tab]);
    printf(" %% current nAch_open given by %f\n",*(n->vpra[VARNAME_REGISTRY_mainak_nAch_open]));}
  temp_time = t;
  l_spiketime = n->pulse_synapse->spiketime_llistra[tab]->first; l_value = n->pulse_synapse->value_llistra[tab]->first;
  while (l_spiketime!=NULL && l_value!=NULL && *(double *)l_spiketime->item<t+DT){
    temp_spiketime = *(double *)l_spiketime->item;
    if (verbose>2){ printf(" %% found temp_spiketime %f\n",temp_spiketime);}
    *(n->vpra[VARNAME_REGISTRY_mainak_nAch_open]) *= exp(-(temp_spiketime-temp_time)*0.2/* /TAU_[VARNAME_REGISTRY_mainak_nAch_open_local] */);
    *(n->vpra[VARNAME_REGISTRY_mainak_nAch_open]) += *(double *)l_value->item;
    temp_time = temp_spiketime;
    llistkillfirst2(n->pulse_synapse->spiketime_llistra[tab]); llistkillfirst2(n->pulse_synapse->value_llistra[tab]);
    l_spiketime = n->pulse_synapse->spiketime_llistra[tab]->first; l_value = n->pulse_synapse->value_llistra[tab]->first;}
  *(n->vpra[VARNAME_REGISTRY_mainak_nAch_open]) *= exp(-(t+DT-temp_time)*0.2/* /TAU_[VARNAME_REGISTRY_mainak_nAch_open_local] */);
  if (verbose>2){ 
    printf(" %% now nAch llist\n"); 
    llistprintf2(n->pulse_synapse->spiketime_llistra[tab]);
    llistprintf2(n->pulse_synapse->value_llistra[tab]);
    printf(" %% current nAch_open given by %f\n",*(n->vpra[VARNAME_REGISTRY_mainak_nAch_open]));}
  tab=1;
  temp_time = t;
  l_spiketime = n->pulse_synapse->spiketime_llistra[tab]->first; l_value = n->pulse_synapse->value_llistra[tab]->first;
  while (l_spiketime!=NULL && l_value!=NULL && *(double *)l_spiketime->item<t+DT){
    temp_spiketime = *(double *)l_spiketime->item;
    *(n->vpra[VARNAME_REGISTRY_mainak_slow_inh_g]) = *(n->vpra[VARNAME_REGISTRY_mainak_slow_inh_g])*exp(-(temp_spiketime-temp_time)*0.06) + 0.1**(n->vpra[VARNAME_REGISTRY_mainak_slow_inh_r])/(0.06-0.0025)*(exp(-(temp_spiketime-temp_time)*0.0025) - exp(-(temp_spiketime-temp_time)*0.06));
    *(n->vpra[VARNAME_REGISTRY_mainak_slow_inh_r]) *= exp(-(temp_spiketime-temp_time)*0.0025);
    *(n->vpra[VARNAME_REGISTRY_mainak_slow_inh_r]) += *(double *)l_value->item;
    temp_time = temp_spiketime;
    llistkillfirst2(n->pulse_synapse->spiketime_llistra[tab]); llistkillfirst2(n->pulse_synapse->value_llistra[tab]);
    l_spiketime = n->pulse_synapse->spiketime_llistra[tab]->first; l_value = n->pulse_synapse->value_llistra[tab]->first;}
  *(n->vpra[VARNAME_REGISTRY_mainak_slow_inh_g]) = *(n->vpra[VARNAME_REGISTRY_mainak_slow_inh_g])*exp(-(t+DT-temp_time)*0.06) + 0.1**(n->vpra[VARNAME_REGISTRY_mainak_slow_inh_r])/(0.06-0.0025)*(exp(-(t+DT-temp_time)*0.0025) - exp(-(t+DT-temp_time)*0.06));
  *(n->vpra[VARNAME_REGISTRY_mainak_slow_inh_r]) *= exp(-(t+DT-temp_time)*0.0025);
  tab=2;
  temp_time = t;
  *(n->vpra[VARNAME_REGISTRY_mainak_gaba_open]) = 0;
  l_spiketime = n->pulse_synapse->spiketime_llistra[tab]->first; l_value = n->pulse_synapse->value_llistra[tab]->first;
  while (l_spiketime!=NULL && l_value!=NULL && *(double *)l_spiketime->item<t+DT){
    *(n->vpra[VARNAME_REGISTRY_mainak_gaba_open]) += *(double *)l_value->item;
    llistkillfirst2(n->pulse_synapse->spiketime_llistra[tab]); llistkillfirst2(n->pulse_synapse->value_llistra[tab]);
    l_spiketime = n->pulse_synapse->spiketime_llistra[tab]->first; l_value = n->pulse_synapse->value_llistra[tab]->first;}
  if (verbose){ printf(" %% [finishing mainak_continuous_synapse_evolve]\n");}
}

void mainak_continuous_synapse_evolve_old(struct neuron *n,double t,double DT)
{
  int verbose=0;//(n->type==0 && n->index==0);
  int tab=0,tab_max=3;
  /* double temp=0; */
  struct litem *l=NULL;
  struct neuron *s=NULL;
  if (verbose){ printf(" %% [entering mainak_continuous_synapse_evolve] with n(%d,%d),time %f->%f\n",n->type,n->index,t,t+DT);}
  *(n->vpra[VARNAME_REGISTRY_mainak_nAch_open])=0;
  *(n->vpra[VARNAME_REGISTRY_mainak_gaba_open])=0;
  *(n->vpra[VARNAME_REGISTRY_mainak_slow_inh_r])=0;
  l = n->pulse_synapse->pre->first;
  while (l!=NULL){
    s = (struct neuron *) l->item;
    if (s->type==0){ *(n->vpra[VARNAME_REGISTRY_mainak_nAch_open]) += *(s->vpra[VARNAME_REGISTRY_mainak_nAch_open_local]);}
    if (s->type==1){ *(n->vpra[VARNAME_REGISTRY_mainak_gaba_open]) += *(s->vpra[VARNAME_REGISTRY_mainak_gaba_open_local]);}
    if (n->type==0 && s->type==1){ 
      *(n->vpra[VARNAME_REGISTRY_mainak_slow_inh_r]) += *(s->vpra[VARNAME_REGISTRY_mainak_slow_inh_r_local]);}
    l=l->child;}
  if (n->type==0){
     *(n->vpra[VARNAME_REGISTRY_mainak_slow_inh_g]) = *(n->vpra[VARNAME_REGISTRY_mainak_slow_inh_g])*exp(-DT*0.06) + 0.1**(n->vpra[VARNAME_REGISTRY_mainak_slow_inh_r])/(0.06-0.0025)*(exp(-DT*0.0025) - exp(-DT*0.06));}
    //temp = *(n->vpra[VARNAME_REGISTRY_mainak_slow_inh_g]); *(n->vpra[VARNAME_REGISTRY_mainak_slow_inh_g]) = temp + DT*(0.1**(n->vpra[VARNAME_REGISTRY_mainak_slow_inh_r]) - 0.06*temp);}
  for (tab=0;tab<tab_max;tab++){
    llisttfree2(n->pulse_synapse->spiketime_llistra[tab]); 
    llisttfree2(n->pulse_synapse->value_llistra[tab]);
    n->pulse_synapse->spiketime_llistra[tab] = llistmake();
    n->pulse_synapse->value_llistra[tab] = llistmake();}
  if (verbose){ printf(" %% [finishing mainak_continuous_synapse_evolve]\n");}
}

void mainak_continuous_evolve_helper_rhs_j_inv_ln(double *ra_rhs,double *ra_j)
{
  /*
    assumes that ra_rhs has component ordering
    |  V  |  Ca | mCaK| m_Ca| h_Ca| m_K |
    assumes that ra_j has the structure
    |  V  |  Ca | mCaK| m_Ca| h_Ca| m_K |
    +-----+-----+-----+-----+-----+-----+
    |  0  |     |  1  |  2  |  3  |  4  |
    +-----+-----+-----+-----+-----+-----+
    |  5  |  6  |     |  7  |  8  |     |
    +-----+-----+-----+-----+-----+-----+
    |     |  9  |  10 |     |     |     |
    +-----+-----+-----+-----+-----+-----+
    |  11 |     |     |  12 |     |     |
    +-----+-----+-----+-----+-----+-----+
    |  13 |     |     |     |  14 |     |
    +-----+-----+-----+-----+-----+-----+
    |  15 |     |     |     |     |  16 |
    +-----+-----+-----+-----+-----+-----+
   */
  int verbose=0;
  double temp=0,j_01=0,j_20=0;
  double *ra1=NULL,*ra2=NULL;
  if (verbose){ 
    ra1=(double *) tcalloc(6,sizeof(double)); raplusequals(ra1,6,ra_rhs);
    ra2=(double *) tcalloc(17,sizeof(double)); raplusequals(ra2,17,ra_j);}
  temp = -ra_j[4]/ra_j[16]; ra_j[0] += temp*ra_j[15]; ra_rhs[0] += temp*ra_rhs[5];
  temp = -ra_j[8]/ra_j[14]; ra_j[5] += temp*ra_j[13]; ra_rhs[1] += temp*ra_rhs[4];
  temp = -ra_j[3]/ra_j[14]; ra_j[0] += temp*ra_j[13]; ra_rhs[0] += temp*ra_rhs[4];
  temp = -ra_j[7]/ra_j[12]; ra_j[5] += temp*ra_j[11]; ra_rhs[1] += temp*ra_rhs[3];
  temp = -ra_j[2]/ra_j[12]; ra_j[0] += temp*ra_j[11]; ra_rhs[0] += temp*ra_rhs[3];
  temp = -ra_j[1]/ra_j[10]; ra_rhs[0] += temp*ra_rhs[2]; j_01 = temp*ra_j[9];
  temp = -j_01/ra_j[6]; ra_j[0] += temp*ra_j[5]; ra_rhs[0] += temp*ra_rhs[1];
  temp = -ra_j[9]/ra_j[6]; j_20 = temp*ra_j[5]; ra_rhs[2] += temp*ra_rhs[1];  
  temp = -ra_j[5]/ra_j[0];  ra_rhs[1] += temp*ra_rhs[0];
  temp = -j_20/ra_j[0];  ra_rhs[2] += temp*ra_rhs[0];
  temp = -ra_j[11]/ra_j[0]; ra_rhs[3] += temp*ra_rhs[0];
  temp = -ra_j[13]/ra_j[0]; ra_rhs[4] += temp*ra_rhs[0];
  temp = -ra_j[15]/ra_j[0]; ra_rhs[5] += temp*ra_rhs[0];
  ra_rhs[0] /= ra_j[0];
  ra_rhs[1] /= ra_j[6];
  ra_rhs[2] /= ra_j[10];
  ra_rhs[3] /= ra_j[12];
  ra_rhs[4] /= ra_j[14];
  ra_rhs[5] /= ra_j[16];
  if (verbose){
    printf("ra2*ra_rhs =  %f, %f, %f, %f, %f, %f\n",ra2[0]*ra_rhs[0]+ra2[1]*ra_rhs[2]+ra2[2]*ra_rhs[3]+ra2[3]*ra_rhs[4]+ra2[4]*ra_rhs[5],ra2[5]*ra_rhs[0]+ra2[6]*ra_rhs[1]+ra2[7]*ra_rhs[3]+ra2[8]*ra_rhs[4],ra2[9]*ra_rhs[1]+ra2[10]*ra_rhs[2],ra2[11]*ra_rhs[0]+ra2[12]*ra_rhs[3],ra2[13]*ra_rhs[0]+ra2[14]*ra_rhs[4],ra2[15]*ra_rhs[0]+ra2[16]*ra_rhs[5]);
    raprintf(ra1,"double",1,6,"ra1 = ");
    tfree(ra1);ra1=NULL;tfree(ra2);ra2=NULL;}
}

void mainak_continuous_evolve_helper_rhs_ln(double *ra,double *ra_rhs,double *ra_j)
{
  /* presumes that ra has component ordering
     V,Ca,m_CaK,m_Ca,h_Ca,m_K,nAch_open,gaba_open,s_ORN
     ra_rhs and ra_j are similarly ordered */
  int nt1=0,nt2=0,nv=0;
  double V=0,Ca=0,m_CaK=0,m_K=0,m_Ca=0,h_Ca=0,nAch_open=0,gaba_open=0,s_ORN=0;
  double alpha=0,alpha_V=0,beta=0,beta_V=0;
  double temp0=0,temp0p=0,temp1=0,temp1p=0,temp2=0,temp2p=0,temp3=0,temp3p=0;
  double I_LEAK_S=0;
  double I_LEAK_S_V=0;
  double I_CaK=0;
  double I_CaK_V=0;
  double I_CaK_m_CaK=0;
  double I_K=0;
  double I_K_V=0;
  double I_K_m_K=0;
  double I_Ca=0;
  double I_Ca_V=0;
  double I_Ca_m_Ca=0;
  double I_Ca_h_Ca=0;
  double I_nAch=0;
  double I_nAch_V=0;
  double I_nAch_nAch_open=0;
  double I_gaba=0;
  double I_gaba_V=0;
  double I_gaba_gaba_open=0;
  double I_orn=0;
  double I_orn_V=0;
  double rhs_m_CaK=0;
  double rhs_m_CaK_Ca=0;
  double rhs_m_CaK_m_CaK=0;
  double rhs_m_K=0;
  double rhs_m_K_V=0;
  double rhs_m_K_m_K=0;
  double rhs_m_Ca=0;
  double rhs_m_Ca_V=0;
  double rhs_m_Ca_m_Ca=0;
  double rhs_h_Ca=0;
  double rhs_h_Ca_V=0;
  double rhs_h_Ca_h_Ca=0;
  double rhs_Ca=0;
  double rhs_Ca_Ca=0;
  double rhs_Ca_V=0;
  double rhs_Ca_m_Ca=0;
  double rhs_Ca_h_Ca=0;
  V=ra[0];
  Ca=ra[1];
  m_CaK=ra[2];
  m_Ca=ra[3];
  h_Ca=ra[4];
  m_K=ra[5];
  nAch_open=ra[6];
  gaba_open=ra[7];
  s_ORN=ra[8];
  I_LEAK_S = -/* CONDUCTANCE_[VARNAME_REGISTRY_mainak_LEAK_LN] */0.3*(V - -50.0/* VOLTAGE_[VARNAME_REGISTRY_mainak_LEAK_LN] */);
  I_LEAK_S_V = -/* CONDUCTANCE_[VARNAME_REGISTRY_mainak_LEAK_LN] */0.3;
  I_Ca = -/*CONDUCTANCE_[VARNAME_REGISTRY_mainak_m_Ca]*/5.0*pow(m_Ca,2)*pow(h_Ca,1)*(V-/*VOLTAGE_[VARNAME_REGISTRY_mainak_m_Ca]*/140.0);
  I_Ca_V = -/*CONDUCTANCE_[VARNAME_REGISTRY_mainak_m_Ca]*/5.0*pow(m_Ca,2)*pow(h_Ca,1);
  I_Ca_m_Ca = -/*CONDUCTANCE_[VARNAME_REGISTRY_mainak_m_Ca]*/5.0*2*pow(m_Ca,1)*pow(h_Ca,1)*(V-/*VOLTAGE_[VARNAME_REGISTRY_mainak_m_Ca]*/140.0);
  I_Ca_h_Ca = -/*CONDUCTANCE_[VARNAME_REGISTRY_mainak_m_Ca]*/5.0*pow(m_Ca,2)*(V-/*VOLTAGE_[VARNAME_REGISTRY_mainak_m_Ca]*/140.0);
  I_CaK = -/*CONDUCTANCE_[VARNAME_REGISTRY_mainak_m_CaK]*/0.045*pow(m_CaK,1)*(V- /*VOLTAGE_[VARNAME_REGISTRY_mainak_m_CaK]*/-95.0);
  I_CaK_V = -/*CONDUCTANCE_[VARNAME_REGISTRY_mainak_m_CaK]*/0.045*pow(m_CaK,1);
  I_CaK_m_CaK = -/*CONDUCTANCE_[VARNAME_REGISTRY_mainak_m_CaK]*/0.045*(V- /*VOLTAGE_[VARNAME_REGISTRY_mainak_m_CaK]*/-95.0);
  I_K = -/*CONDUCTANCE_[VARNAME_REGISTRY_mainak_m_K]*/36.0*pow(m_K,4)*(V- /*VOLTAGE_[VARNAME_REGISTRY_mainak_m_K]*/-95.0);
  I_K_V = -/*CONDUCTANCE_[VARNAME_REGISTRY_mainak_m_K]*/36.0*pow(m_K,4);
  I_K_m_K = -/*CONDUCTANCE_[VARNAME_REGISTRY_mainak_m_K]*/36.0*4*pow(m_K,3)*(V- /*VOLTAGE_[VARNAME_REGISTRY_mainak_m_K]*/-95.0);
  rhs_Ca = 0.0002*I_Ca - (Ca - 0.00024)/150;
  rhs_Ca_Ca = -Ca/150;
  rhs_Ca_V = 0.0002*I_Ca_V;
  rhs_Ca_m_Ca = 0.0002*I_Ca_m_Ca;
  rhs_Ca_h_Ca = 0.0002*I_Ca_h_Ca;
  temp0 = Ca/(Ca+2.0); temp0p = 2.0/pow(Ca+2.0,2); temp1 = 100/(Ca+2.0); temp1p = -100/pow(Ca+2.0,2);
  rhs_m_CaK = (temp0 - m_CaK)/temp1;
  rhs_m_CaK_m_CaK = -1.0/temp1;
  rhs_m_CaK_Ca = (temp0p*temp1 - (temp0-m_CaK)*temp1p)/pow(temp1,2);
  temp0 = V+20.0; temp1 = exp(-temp0/6.5); temp1p = -1.0/6.5*temp1; temp2 = 1.0/(1.0+temp1); temp2p = -pow(temp2,2)*temp1p;
  temp3 = 1 + (V+30)*0.014; temp3p = 0.014;
  rhs_m_Ca = (temp2 - m_Ca)/temp3;
  rhs_m_Ca_m_Ca = -1.0/temp3;
  rhs_m_Ca_V = (temp2p*temp3 - (temp2-m_Ca)*temp3p)/pow(temp3,2);
  temp0 = V+25.0; temp1 = exp(temp0/12); temp1p = 1.0/12.0*temp1; temp2 = 1.0/(1.0+temp1); temp2p = -pow(temp2,2)*temp1p;
  temp3 = 0.3*exp((V-40)/13) + 0.002*exp(-(V-60)/29); temp3p = 0.3/30.0*exp((V-40)/13) - 0.002/29.0*exp(-(V-60)/29);
  rhs_h_Ca = (temp2 - h_Ca)/temp3;
  rhs_h_Ca_h_Ca = -1.0/temp3;
  rhs_h_Ca_V = (temp2p*temp3 - (temp2-m_Ca)*temp3p)/pow(temp3,2);
  temp0 = V+35; temp1 = -exp(-0.2*temp0); temp3p = -0.2*temp1; temp3 = 1.0+temp1;
  alpha = 0.02*temp0/temp3; alpha_V = 0.02*(1.0/temp3 - temp0/pow(temp3,2)*temp3p);
  temp0 = V+40; temp1 = 0.5*exp(-temp0/40); temp1p = -1.0/40.0*temp1; 
  beta = temp1; beta_V = temp1p;
  rhs_m_K = (alpha*(1.0 - m_K) - beta*m_K)/4.65;
  rhs_m_K_V = (alpha_V*(1.0 - m_K) - beta_V*m_K)/4.65;
  rhs_m_K_m_K = (-alpha-beta)/4.65;
  nt1 = TYPENAME_REGISTRY_PN;
  nt2 = TYPENAME_REGISTRY_LN;
  nv = GLOBAL_INDEXING_REFILE_sra[VARNAME_REGISTRY_mainak_nAch_open];
  temp0 = CS__[nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES];
  temp0 *= GLOBAL_CS_PRETYPE_SCALE_[nt1]*GLOBAL_CS_POSTYPE_SCALE_[nt2]*GLOBAL_CS_SRA_SCALE_[GLOBAL_INDEXING_REFILE_sra[VARNAME_REGISTRY_mainak_nAch_open]];
  I_nAch = -temp0*nAch_open*(V-VOLTAGE_[VARNAME_REGISTRY_mainak_nAch_open]);
  I_nAch_V = -temp0*nAch_open;
  I_nAch_nAch_open = -temp0*(V-VOLTAGE_[VARNAME_REGISTRY_mainak_nAch_open]);
  nt1 = TYPENAME_REGISTRY_LN;
  nt2 = TYPENAME_REGISTRY_LN;
  nv = GLOBAL_INDEXING_REFILE_sra[VARNAME_REGISTRY_mainak_gaba_open];
  temp0 = CS__[nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES];
  temp0 *= GLOBAL_CS_PRETYPE_SCALE_[nt1]*GLOBAL_CS_POSTYPE_SCALE_[nt2]*GLOBAL_CS_SRA_SCALE_[GLOBAL_INDEXING_REFILE_sra[VARNAME_REGISTRY_mainak_gaba_open]];
  I_gaba = -temp0*gaba_open*(V-VOLTAGE_[VARNAME_REGISTRY_mainak_gaba_open]);
  I_gaba_V = -temp0*gaba_open;
  I_gaba_gaba_open = -temp0*(V-VOLTAGE_[VARNAME_REGISTRY_mainak_gaba_open]);
  I_orn = s_ORN;
  I_orn_V = 0;
  ra_rhs[0] = I_LEAK_S + I_Ca + I_CaK + I_K + I_nAch + I_gaba + I_orn + CURRENT_INJECTION_S;
  ra_rhs[1] = rhs_Ca;
  ra_rhs[2] = rhs_m_CaK;
  ra_rhs[3] = rhs_m_Ca;
  ra_rhs[4] = rhs_h_Ca;
  ra_rhs[5] = rhs_m_K;
  ra_j[0] = I_LEAK_S_V + I_Ca_V + I_CaK_V + I_K_V + I_nAch_V + I_gaba_V + I_orn_V;
  ra_j[1] = I_CaK_m_CaK;
  ra_j[2] = I_Ca_m_Ca;
  ra_j[3] = I_Ca_h_Ca;
  ra_j[4] = I_K_m_K;
  ra_j[5] = rhs_Ca_V;
  ra_j[6] = rhs_Ca_Ca;
  ra_j[7] = rhs_Ca_m_Ca;
  ra_j[8] = rhs_Ca_h_Ca;
  ra_j[9] = rhs_m_CaK_Ca;
  ra_j[10] = rhs_m_CaK_m_CaK;
  ra_j[11] = rhs_m_Ca_V;
  ra_j[12] = rhs_m_Ca_m_Ca;
  ra_j[13] = rhs_h_Ca_V;
  ra_j[14] = rhs_h_Ca_h_Ca;
  ra_j[15] = rhs_m_K_V;
  ra_j[16] = rhs_m_K_m_K;
}

void mainak_continuous_evolve_helper_rhs_j_inv_pn(double *ra_rhs,double *ra_j)
{
  /*
    assumes that ra_rhs has component ordering
    |  V  | m_Na| h_Na| m_K | m_KA| h_KA|
    assumes that ra_j has the structure
    |  V  | m_Na| h_Na| m_K | m_KA| h_KA|
    +-----+-----+-----+-----+-----+-----+
    |  0  |  1  |  2  |  3  |  4  |  5  |
    +-----+-----+-----+-----+-----+-----+
    |  6  |  7  |     |     |     |     |
    +-----+-----+-----+-----+-----+-----+
    |  8  |     |  9  |     |     |     |
    +-----+-----+-----+-----+-----+-----+
    |  10 |     |     |  11 |     |     |
    +-----+-----+-----+-----+-----+-----+
    |  12 |     |     |     |  13 |     |
    +-----+-----+-----+-----+-----+-----+
    |  14 |     |     |     |     |  15 |
    +-----+-----+-----+-----+-----+-----+
   */
  int verbose=0;
  double temp=0;
  double *ra1=NULL,*ra2=NULL;
  if (verbose){ 
    ra1=(double *) tcalloc(6,sizeof(double)); raplusequals(ra1,6,ra_rhs);
    ra2=(double *) tcalloc(16,sizeof(double)); raplusequals(ra2,16,ra_j);}
  temp = -ra_j[5]/ra_j[15]; ra_j[0] += temp*ra_j[14]; ra_rhs[0] += temp*ra_rhs[5];
  temp = -ra_j[4]/ra_j[13]; ra_j[0] += temp*ra_j[12]; ra_rhs[0] += temp*ra_rhs[4];
  temp = -ra_j[3]/ra_j[11]; ra_j[0] += temp*ra_j[10]; ra_rhs[0] += temp*ra_rhs[3];
  temp = -ra_j[2]/ra_j[9];  ra_j[0] += temp*ra_j[8];  ra_rhs[0] += temp*ra_rhs[2];
  temp = -ra_j[1]/ra_j[7];  ra_j[0] += temp*ra_j[6];  ra_rhs[0] += temp*ra_rhs[1];
  temp = -ra_j[6]/ra_j[0];  ra_rhs[1] += temp*ra_rhs[0];
  temp = -ra_j[8]/ra_j[0];  ra_rhs[2] += temp*ra_rhs[0];
  temp = -ra_j[10]/ra_j[0]; ra_rhs[3] += temp*ra_rhs[0];
  temp = -ra_j[12]/ra_j[0]; ra_rhs[4] += temp*ra_rhs[0];
  temp = -ra_j[14]/ra_j[0]; ra_rhs[5] += temp*ra_rhs[0];
  ra_rhs[0] /= ra_j[0];
  ra_rhs[1] /= ra_j[7];
  ra_rhs[2] /= ra_j[9];
  ra_rhs[3] /= ra_j[11];
  ra_rhs[4] /= ra_j[13];
  ra_rhs[5] /= ra_j[15];
  if (verbose){
    printf("ra2*ra_rhs =  %f, %f, %f, %f, %f, %f\n",ra2[0]*ra_rhs[0]+ra2[1]*ra_rhs[1]+ra2[2]*ra_rhs[2]+ra2[3]*ra_rhs[3]+ra2[4]*ra_rhs[4]+ra2[5]*ra_rhs[5],ra2[6]*ra_rhs[0]+ra2[7]*ra_rhs[1],ra2[8]*ra_rhs[0]+ra2[9]*ra_rhs[2],ra2[10]*ra_rhs[0]+ra2[11]*ra_rhs[3],ra2[12]*ra_rhs[0]+ra2[13]*ra_rhs[4],ra2[14]*ra_rhs[0]+ra2[15]*ra_rhs[5]);
    raprintf(ra1,"double",1,6,"ra1 = ");
    tfree(ra1);ra1=NULL;tfree(ra2);ra2=NULL;}
}

void mainak_continuous_evolve_helper_rhs_pn(double *ra,double *ra_rhs,double *ra_j)
{
  /* presumes that ra has component ordering
     V,m_Na,h_Na,m_K,m_KA,h_KA,nAch_open,gaba_open,slow_inh_g,s_ORN
     ra_rhs and ra_j are similarly ordered */
  int nt1=0,nt2=0,nv=0;
  double V=0,m_Na=0,h_Na=0,m_K=0,m_KA=0,h_KA=0,nAch_open=0,gaba_open=0,slow_inh_g=0,s_ORN=0;
  double alpha=0,alpha_V=0,beta=0,beta_V=0;
  double tempx=0,tempw=0,temp0=0,temp0p=0,temp1=0,temp1p=0,temp2=0,temp2p=0,temp3=0,temp3p=0;
  double I_LEAK_S=0;
  double I_LEAK_S_V=0;
  double I_Na=0;
  double I_Na_V=0;
  double I_Na_m_Na=0;
  double I_Na_h_Na=0;
  double I_K=0;
  double I_K_V=0;
  double I_K_m_K=0;
  double I_KA=0;
  double I_KA_V=0;
  double I_KA_m_KA=0;
  double I_KA_h_KA=0;
  double I_nAch=0;
  double I_nAch_V=0;
  double I_nAch_nAch_open=0;
  double I_gaba=0;
  double I_gaba_V=0;
  double I_gaba_gaba_open=0;
  double I_slow_inh=0;
  double I_slow_inh_V=0;
  double I_slow_inh_slow_inh_g=0;
  double I_orn=0;
  double I_orn_V=0;
  double rhs_m_Na=0;
  double rhs_m_Na_V=0;
  double rhs_m_Na_m_Na=0;
  double rhs_h_Na=0;
  double rhs_h_Na_V=0;
  double rhs_h_Na_h_Na=0;
  double rhs_m_K=0;
  double rhs_m_K_V=0;
  double rhs_m_K_m_K=0;
  double rhs_m_KA=0;
  double rhs_m_KA_V=0;
  double rhs_m_KA_m_KA=0;
  double rhs_h_KA=0;
  double rhs_h_KA_V=0;
  double rhs_h_KA_h_KA=0;
  V=ra[0];
  m_Na=ra[1];
  h_Na=ra[2];
  m_K=ra[3];
  m_KA=ra[4];
  h_KA=ra[5];
  nAch_open=ra[6];
  gaba_open=ra[7];
  slow_inh_g=ra[8];
  s_ORN=ra[9];
  I_LEAK_S = -/* CONDUCTANCE_[VARNAME_REGISTRY_mainak_LEAK_PN] */0.3*(V - -64.387/* VOLTAGE_[VARNAME_REGISTRY_mainak_LEAK_PN] */);
  I_LEAK_S_V = -/* CONDUCTANCE_[VARNAME_REGISTRY_mainak_LEAK_PN] */0.3;
  I_Na = -/* CONDUCTANCE_[VARNAME_REGISTRY_mainak_m_Na] */120*pow(m_Na,3)*pow(h_Na,1)*(V-/* VOLTAGE_[VARNAME_REGISTRY_mainak_m_Na] */ 40);
  I_Na_V = -/* CONDUCTANCE_[VARNAME_REGISTRY_mainak_m_Na] */120*pow(m_Na,3)*pow(h_Na,1);
  I_Na_m_Na = -/* CONDUCTANCE_[VARNAME_REGISTRY_mainak_m_Na] */120*3*pow(m_Na,2)*pow(h_Na,1)*(V-/* VOLTAGE_[VARNAME_REGISTRY_mainak_m_Na] */ 40);
  I_Na_h_Na = -/* CONDUCTANCE_[VARNAME_REGISTRY_mainak_m_Na] */120*pow(m_Na,3)*(V-/* VOLTAGE_[VARNAME_REGISTRY_mainak_m_Na] */ 40);
  I_K = -/* CONDUCTANCE_[VARNAME_REGISTRY_mainak_m_K] */3.6*pow(m_K,4)*(V-/* VOLTAGE_[VARNAME_REGISTRY_mainak_m_K] */ -87);
  I_K_V = -/* CONDUCTANCE_[VARNAME_REGISTRY_mainak_m_K] */3.6*pow(m_K,4);
  I_K_m_K = -/* CONDUCTANCE_[VARNAME_REGISTRY_mainak_m_K] */3.6*4*pow(m_K,3)*(V-/* VOLTAGE_[VARNAME_REGISTRY_mainak_m_K] */ -87);
  I_KA = -/* CONDUCTANCE_[VARNAME_REGISTRY_mainak_m_KA] */1.43*pow(m_KA,4)*pow(h_KA,1)*(V-/* VOLTAGE_[VARNAME_REGISTRY_mainak_m_KA] */ -87);
  I_KA_V = -/* CONDUCTANCE_[VARNAME_REGISTRY_mainak_m_KA] */1.43*pow(m_KA,4)*pow(h_KA,1);
  I_KA_m_KA = -/* CONDUCTANCE_[VARNAME_REGISTRY_mainak_m_KA] */1.43*4*pow(m_KA,3)*pow(h_KA,1)*(V-/* VOLTAGE_[VARNAME_REGISTRY_mainak_m_KA] */ -87);
  I_KA_h_KA = -/* CONDUCTANCE_[VARNAME_REGISTRY_mainak_m_KA] */1.43*pow(m_KA,4)*(V-/* VOLTAGE_[VARNAME_REGISTRY_mainak_m_KA] */ -87);
  temp0 = V+35.0; temp1 = -exp(-0.1*temp0); temp3p = -0.1*temp1; temp3 = 1.0+temp1;
  alpha = 0.1*temp0/temp3; alpha_V = 0.1*(1.0/temp3 - temp0/pow(temp3,2)*temp3p);
  temp0 = V+60.0; temp1 = 4.0*exp(-temp0/18.0); temp1p = -1.0/18.0*temp1;
  beta = temp1; beta_V = temp1p;
  rhs_m_Na = alpha*(1.0 - m_Na) - beta*m_Na;
  rhs_m_Na_V = alpha_V*(1.0 - m_Na) - beta_V*m_Na;
  rhs_m_Na_m_Na = -alpha-beta;
  temp0 = V+60.0; temp1 = 0.07*exp(-temp0/20.0); temp1p = -1.0/20.0*temp1;
  alpha = temp1; alpha_V = temp1p;
  temp0 = V+30.0; temp1 = exp(-0.1*temp0); temp3p = -0.1*temp1; temp3 = 1.0+temp1;
  beta = 1.0/temp3; beta_V = -1.0/pow(temp3,2)*temp3p;
  rhs_h_Na = alpha*(1.0 - h_Na) - beta*h_Na;
  rhs_h_Na_V = alpha_V*(1.0 - h_Na) - beta_V*h_Na;
  rhs_h_Na_h_Na = -alpha-beta;
  temp0 = V+50.0; temp1 = -exp(-0.1*temp0); temp3p = -0.1*temp1; temp3 = 1.0+temp1;
  alpha = 0.01*temp0/temp3; alpha_V = 0.01*(1.0/temp3 - temp0/pow(temp3,2)*temp3p);
  temp0 = V+60.0; temp1 = 0.125*exp(-temp0/80.0); temp1p = -1.0/80.0*temp1;
  beta = temp1; beta_V = temp1p;
  rhs_m_K = alpha*(1.0 - m_K) - beta*m_K;
  rhs_m_K_V = alpha_V*(1.0 - m_K) - beta_V*m_K;
  rhs_m_K_m_K = -alpha-beta;
  temp0 = V+60.0; temp1 = exp(-temp0/8.5); temp3p = -1.0/8.5*temp1; temp3 = 1.0+temp1;
  alpha = 1.0/temp3; alpha_V = -1.0/pow(temp3,2)*temp3p;
  temp0 = 0.27*exp(-(V+35.8)/19.7); temp0p = -1/19.7*temp0;
  temp2 = exp(-(V+79.7)/12.7); temp2p = -1/12.7*temp2;
  beta = temp0+temp2+0.1; beta_V = temp0p+temp2p;
  rhs_m_KA = (alpha - m_KA)/beta;
  rhs_m_KA_V = (alpha_V*beta - (alpha-m_KA)*beta_V)/pow(beta,2);
  rhs_m_KA_m_KA = -1.0/beta;
  temp0 = V+78.0; temp1 = exp(temp0/6); temp3p = 1.0/6.0*temp1; temp3 = 1.0+temp1;
  alpha = 1.0/temp3; alpha_V = -1.0/pow(temp3,2)*temp3p;
  tempw = 2; tempx = -61.6325;
  temp0 = 0.27/(exp((V+46)/5)+exp(-(V+238)/37.5));
  temp0p = -0.27/pow(exp((V+46)/5)+exp(-(V+238)/37.5),2)*(1/5.0*exp((V+46)/5) - 1.0/37.5*exp(-(V+238)/37.5));
  temp1 = 1 - 0.5*(1+erf((V-tempx)/tempw));
  temp1p = -0.5*2/sqrt(PI)*exp(-pow((V-tempx)/tempw,2))*1.0/tempw;
  temp2 = 0.5*(erf((V-tempx)/tempw)+1);
  temp2p = 0.5*2/sqrt(PI)*exp(-pow((V-tempx)/tempw,2))*1.0/tempw;
  temp3 = tempw*exp(-(V-tempx)/tempw) + 0.27/(exp((tempx+46)/5)+exp(-(tempx+238)/37.5));
  temp3p = exp(-(V-tempx)/tempw);
  beta = temp0*temp1 + temp2*temp3;
  beta_V = temp0p*temp1 + temp0*temp1p + temp2p*temp3 + temp2*temp3p;
  rhs_h_KA = (alpha - h_KA)/beta;
  rhs_h_KA_V = (alpha_V*beta - (alpha-h_KA)*beta_V)/pow(beta,2);
  rhs_h_KA_h_KA = -1.0/beta;
  nt1 = TYPENAME_REGISTRY_PN;
  nt2 = TYPENAME_REGISTRY_PN;
  nv = GLOBAL_INDEXING_REFILE_sra[VARNAME_REGISTRY_mainak_nAch_open];
  temp0 = CS__[nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES];
  temp0 *= GLOBAL_CS_PRETYPE_SCALE_[nt1]*GLOBAL_CS_POSTYPE_SCALE_[nt2]*GLOBAL_CS_SRA_SCALE_[GLOBAL_INDEXING_REFILE_sra[VARNAME_REGISTRY_mainak_nAch_open]];
  I_nAch = -temp0*nAch_open*(V-VOLTAGE_[VARNAME_REGISTRY_mainak_nAch_open]);
  I_nAch_V = -temp0*nAch_open;
  I_nAch_nAch_open = -temp0*(V-VOLTAGE_[VARNAME_REGISTRY_mainak_nAch_open]);
  nt1 = TYPENAME_REGISTRY_LN;
  nt2 = TYPENAME_REGISTRY_PN;
  nv = GLOBAL_INDEXING_REFILE_sra[VARNAME_REGISTRY_mainak_gaba_open];
  temp0 = CS__[nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES];
  temp0 *= GLOBAL_CS_PRETYPE_SCALE_[nt1]*GLOBAL_CS_POSTYPE_SCALE_[nt2]*GLOBAL_CS_SRA_SCALE_[GLOBAL_INDEXING_REFILE_sra[VARNAME_REGISTRY_mainak_gaba_open]];
  I_gaba = -temp0*gaba_open*(V-VOLTAGE_[VARNAME_REGISTRY_mainak_gaba_open]);
  I_gaba_V = -temp0*gaba_open;
  I_gaba_gaba_open = -temp0*(V-VOLTAGE_[VARNAME_REGISTRY_mainak_gaba_open]);
  nt1 = TYPENAME_REGISTRY_LN;
  nt2 = TYPENAME_REGISTRY_PN;
  nv = GLOBAL_INDEXING_REFILE_sra[VARNAME_REGISTRY_mainak_slow_inh_g];
  temp0 = CS__[nt1+nt2*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES];
  temp0 *= GLOBAL_CS_PRETYPE_SCALE_[nt1]*GLOBAL_CS_POSTYPE_SCALE_[nt2]*GLOBAL_CS_SRA_SCALE_[GLOBAL_INDEXING_REFILE_sra[VARNAME_REGISTRY_mainak_slow_inh_g]];
  temp1 = pow(slow_inh_g,4)/(pow(slow_inh_g,4)+100);
  temp2 = 4*pow(slow_inh_g,3)/(pow(slow_inh_g,4)+100) - pow(slow_inh_g,4)/pow(pow(slow_inh_g,4)+100,2)*4*pow(slow_inh_g,3);
  I_slow_inh = -temp0*temp1*(V-VOLTAGE_[VARNAME_REGISTRY_mainak_slow_inh_g]);
  I_slow_inh_V = -temp0*temp1;
  I_slow_inh_slow_inh_g = -temp0*temp2*(V-VOLTAGE_[VARNAME_REGISTRY_mainak_slow_inh_g]);
  I_orn = s_ORN;
  I_orn_V = 0;
  ra_rhs[0] = I_LEAK_S + I_Na + I_K + I_KA + I_nAch + I_gaba + I_slow_inh + I_orn + CURRENT_INJECTION_S;
  ra_rhs[1] = rhs_m_Na;
  ra_rhs[2] = rhs_h_Na;
  ra_rhs[3] = rhs_m_K;
  ra_rhs[4] = rhs_m_KA;
  ra_rhs[5] = rhs_h_KA;
  ra_j[0] = I_LEAK_S_V + I_Na_V + I_K_V + I_KA_V + I_nAch_V + I_gaba_V + I_slow_inh_V + I_orn_V;
  ra_j[1] = I_Na_m_Na;
  ra_j[2] = I_Na_h_Na;
  ra_j[3] = I_K_m_K;
  ra_j[4] = I_KA_m_KA;
  ra_j[5] = I_KA_h_KA;
  ra_j[6] = rhs_m_Na_V;
  ra_j[7] = rhs_m_Na_m_Na;
  ra_j[8] = rhs_h_Na_V;
  ra_j[9] = rhs_h_Na_h_Na;
  ra_j[10] = rhs_m_K_V;
  ra_j[11] = rhs_m_K_m_K;
  ra_j[12] = rhs_m_KA_V;
  ra_j[13] = rhs_m_KA_m_KA;
  ra_j[14] = rhs_h_KA_V;
  ra_j[15] = rhs_h_KA_h_KA;
}

void mainak_continuous_evolve(struct neuronarray *Nra,double t,double DT)
{
  /* presumes type 0 is PN and type 1 (if it exists) is LN */
  int verbose=GLOBAL_verbose;
  int nt=0,nr=0,ns=0,iteration=0,iteration_max=0;
  double dt=0,norm=0,tolerance=0.000001;
  double temp_time=0,temp_dt=0,temp_spiketime=0,temp0=0;
  double *ra=NULL,*ra_rhs=NULL,*ra_j=NULL,*ra_temp1=NULL,*ra_temp2=NULL,*ra_temp3=NULL;
  struct neuron *n=NULL;
  int rhs_size=11,j_size=17;
  ra = (double *) tcalloc(rhs_size,sizeof(double));
  ra_rhs = (double *) tcalloc(rhs_size,sizeof(double));
  ra_j = (double *) tcalloc(j_size,sizeof(double));
  ra_temp1 = (double *) tcalloc(j_size,sizeof(double));
  ra_temp2 = (double *) tcalloc(j_size,sizeof(double));
  ra_temp3 = (double *) tcalloc(j_size,sizeof(double));
  for (nt=0;nt<Nra->ntypes;nt++){ for (nr=0;nr<Nra->lengthra[nt];nr++){
    n=nget(Nra,nt,nr);
    *(n->vpra[VARNAME_REGISTRY_mainak_inputrate]) = n->inputrate;
    *(n->vpra[VARNAME_REGISTRY_mainak_spike_flag]) *= exp(-DT);
    switch (n->type){
    case 0: 
      if (verbose){ printf(" %% n(%d,%d) presumed a PN, %d microsteps\n",nt,nr,n->microstep);}
      temp_time = t;
      n->microstep = 1;
      dt = DT/n->microstep;
      for (ns=0;ns<n->microstep;ns++){
	*(n->vpra[VARNAME_REGISTRY_mainak_s_ORN]) *= exp(-dt/TAU_[VARNAME_REGISTRY_mainak_s_ORN]);
	mainak_continuous_synapse_evolve(n,temp_time,dt);
	ra[0] = *(n->vpra[VARNAME_REGISTRY_mainak_Vs]);
	ra[1] = *(n->vpra[VARNAME_REGISTRY_mainak_m_Na]);
	ra[2] = *(n->vpra[VARNAME_REGISTRY_mainak_h_Na]);
	ra[3] = *(n->vpra[VARNAME_REGISTRY_mainak_m_K]);
	ra[4] = *(n->vpra[VARNAME_REGISTRY_mainak_m_KA]);
	ra[5] = *(n->vpra[VARNAME_REGISTRY_mainak_h_KA]);
	ra[6] = *(n->vpra[VARNAME_REGISTRY_mainak_nAch_open]);
	ra[7] = *(n->vpra[VARNAME_REGISTRY_mainak_gaba_open]);
	ra[8] = *(n->vpra[VARNAME_REGISTRY_mainak_slow_inh_g]);
	ra[9] = *(n->vpra[VARNAME_REGISTRY_mainak_s_ORN]);
	raplugin(ra_temp2,10,1,ra,10,1,0,0);
	iteration=0; iteration_max=10;
	do{
	  mainak_continuous_evolve_helper_rhs_pn(ra_temp2,ra_rhs,ra_j);
	  ratimesequals(ra_j,16,-dt);
	  ra_j[0]  += 1; ra_j[7]  += 1; ra_j[9]  += 1; ra_j[11] += 1; ra_j[13] += 1; ra_j[15] += 1;
	  ra_temp1[0] = dt*ra_rhs[0]+ra[0]-ra_temp2[0];
	  ra_temp1[1] = dt*ra_rhs[1]+ra[1]-ra_temp2[1];
	  ra_temp1[2] = dt*ra_rhs[2]+ra[2]-ra_temp2[2];
	  ra_temp1[3] = dt*ra_rhs[3]+ra[3]-ra_temp2[3];
	  ra_temp1[4] = dt*ra_rhs[4]+ra[4]-ra_temp2[4];
	  ra_temp1[5] = dt*ra_rhs[5]+ra[5]-ra_temp2[5];
	  mainak_continuous_evolve_helper_rhs_j_inv_pn(ra_temp1,ra_j);
	  norm = ra_norm(ra_temp1,6);
	  ra_temp2[0] += ra_temp1[0];
	  ra_temp2[1] += ra_temp1[1];
	  ra_temp2[2] += ra_temp1[2];
	  ra_temp2[3] += ra_temp1[3];
	  ra_temp2[4] += ra_temp1[4];
	  ra_temp2[5] += ra_temp1[5];
	  if (verbose){ printf(" %% iteration %d norm %f\n",iteration,norm);}
	  iteration += 1;}
	while (norm>tolerance && iteration<iteration_max);
	if (iteration>=iteration_max){ if (verbose){ printf(" %% warning! newton not converging in mainak_continuous_evolve\n");}}
	if (*(n->vpra[VARNAME_REGISTRY_mainak_Vs]) > 0 /* mainak_pn_threshold */&& ra_temp2[0] <= 0 /* mainak_pn_threshold */){
	  *(n->vpra[VARNAME_REGISTRY_mainak_spike_flag]) += 1;
	  temp_dt = linerootfinder(-*(n->vpra[VARNAME_REGISTRY_mainak_Vs]),-ra_temp2[0],-0/* mainak_pn_threshold */,dt); 
	  temp_spiketime = temp_time + temp_dt;
	  *(n->vpra[VARNAME_REGISTRY_mainak_nAch_open_local]) *= exp(-temp_dt*0.2/* /TAU_[VARNAME_REGISTRY_mainak_nAch_open_local] */);
	  temp0 = (*(n->vpra[VARNAME_REGISTRY_mainak_nAch_open_local]) - 1.0)*(exp(-1.5) - 1.0);
	  mainak_continuous_synapse_update(n,0,temp0,temp_spiketime+0.125 /* GLOBAL_mainak_synapse_delay */);
	  *(n->vpra[VARNAME_REGISTRY_mainak_nAch_open_local]) += temp0;
	  *(n->vpra[VARNAME_REGISTRY_mainak_nAch_open_local]) *= exp(-(dt-temp_dt)*0.2/* /TAU_[VARNAME_REGISTRY_mainak_nAch_open_local] */);}
	else /* if (*(n->vpra[VARNAME_REGISTRY_mainak_Vs]) <= 0 || ra_temp2[0] > 0) */{
	  *(n->vpra[VARNAME_REGISTRY_mainak_nAch_open_local]) *= exp(-dt*0.2/* /TAU_[VARNAME_REGISTRY_mainak_nAch_open_local] */);}
	*(n->vpra[VARNAME_REGISTRY_mainak_Vs]) = ra_temp2[0];
	*(n->vpra[VARNAME_REGISTRY_mainak_m_Na]) = ra_temp2[1];
	*(n->vpra[VARNAME_REGISTRY_mainak_h_Na]) = ra_temp2[2];
	*(n->vpra[VARNAME_REGISTRY_mainak_m_K]) = ra_temp2[3];
	*(n->vpra[VARNAME_REGISTRY_mainak_m_KA]) = ra_temp2[4];
	*(n->vpra[VARNAME_REGISTRY_mainak_h_KA]) = ra_temp2[5];
	temp_time += dt;}
      break;
    case 1:
      if (verbose){ printf(" n(%d,%d) presumed a LN\n",nt,nr);}
      temp_time = t;
      n->microstep = 1;
      dt = DT/n->microstep;
      for (ns=0;ns<n->microstep;ns++){
	*(n->vpra[VARNAME_REGISTRY_mainak_s_ORN]) *= exp(-dt/TAU_[VARNAME_REGISTRY_mainak_s_ORN]);
	mainak_continuous_synapse_evolve(n,temp_time,dt);
	ra[0] = *(n->vpra[VARNAME_REGISTRY_mainak_Vs]);
	ra[1] = *(n->vpra[VARNAME_REGISTRY_mainak_Ca]);
	ra[2] = *(n->vpra[VARNAME_REGISTRY_mainak_m_CaK]);
	ra[3] = *(n->vpra[VARNAME_REGISTRY_mainak_m_Ca]);
	ra[4] = *(n->vpra[VARNAME_REGISTRY_mainak_h_Ca]);
	ra[5] = *(n->vpra[VARNAME_REGISTRY_mainak_m_K]);
	ra[6] = *(n->vpra[VARNAME_REGISTRY_mainak_nAch_open]);
	ra[7] = *(n->vpra[VARNAME_REGISTRY_mainak_gaba_open]);
	ra[8] = *(n->vpra[VARNAME_REGISTRY_mainak_s_ORN]);
	raplugin(ra_temp2,9,1,ra,9,1,0,0);
	iteration=0; iteration_max=10;
	do{
	  mainak_continuous_evolve_helper_rhs_ln(ra_temp2,ra_rhs,ra_j);
	  ratimesequals(ra_j,17,-dt);
	  ra_j[0]  += 1; ra_j[6]  += 1; ra_j[10] += 1; ra_j[12] += 1; ra_j[14] += 1; ra_j[16] += 1;
	  ra_temp1[0] = dt*ra_rhs[0]+ra[0]-ra_temp2[0];
	  ra_temp1[1] = dt*ra_rhs[1]+ra[1]-ra_temp2[1];
	  ra_temp1[2] = dt*ra_rhs[2]+ra[2]-ra_temp2[2];
	  ra_temp1[3] = dt*ra_rhs[3]+ra[3]-ra_temp2[3];
	  ra_temp1[4] = dt*ra_rhs[4]+ra[4]-ra_temp2[4];
	  ra_temp1[5] = dt*ra_rhs[5]+ra[5]-ra_temp2[5];
	  mainak_continuous_evolve_helper_rhs_j_inv_ln(ra_temp1,ra_j);
	  norm = ra_norm(ra_temp1,6);
	  ra_temp2[0] += ra_temp1[0];
	  ra_temp2[1] += ra_temp1[1];
	  ra_temp2[2] += ra_temp1[2];
	  ra_temp2[3] += ra_temp1[3];
	  ra_temp2[4] += ra_temp1[4];
	  ra_temp2[5] += ra_temp1[5];
	  if (verbose){ printf(" %% iteration %d norm %f\n",iteration,norm);}
	  iteration += 1;}
	while (norm>tolerance && iteration<iteration_max);
	if (iteration>=iteration_max){ if (verbose){ printf(" %% warning! newton not converging in mainak_continuous_evolve\n");}}
	if (*(n->vpra[VARNAME_REGISTRY_mainak_Vs]) > -20 /* mainak_ln_threshold */ && ra_temp2[0] <= -20 /* mainak_ln_threshold */){
	  temp_dt = linerootfinder(-*(n->vpra[VARNAME_REGISTRY_mainak_Vs]),-ra_temp2[0],- -20 /* mainak_ln_threshold */,dt);
	  temp_spiketime = temp_time + temp_dt;
	  *(n->vpra[VARNAME_REGISTRY_mainak_slow_inh_r_local]) *= exp(-temp_dt*0.0025/* /TAU_[VARNAME_REGISTRY_mainak_slow_inh_r_local] */);
	  temp0 = (*(n->vpra[VARNAME_REGISTRY_mainak_slow_inh_r_local]) - 1.0)*(exp(-0.15) - 1.0);
	  mainak_continuous_synapse_update(n,1,temp0,temp_spiketime+0.125 /* GLOBAL_mainak_synapse_delay */);
	  *(n->vpra[VARNAME_REGISTRY_mainak_spike_flag]) += 1;
	  *(n->vpra[VARNAME_REGISTRY_mainak_slow_inh_r_local]) += temp0;
	  *(n->vpra[VARNAME_REGISTRY_mainak_slow_inh_r_local]) *= exp(-(dt-temp_dt)*0.0025/* /TAU_[VARNAME_REGISTRY_mainak_slow_inh_r_local] */);}
	else /* if (*(n->vpra[VARNAME_REGISTRY_mainak_Vs]) <= 0 || ra_temp2[0] > 0) */{
	  *(n->vpra[VARNAME_REGISTRY_mainak_slow_inh_r_local]) *= exp(-dt*0.0025/* /TAU_[VARNAME_REGISTRY_mainak_slow_inh_r_local] */);}
	temp0 = ra_temp2[0];
	temp0 = 1.0/(1.0+exp(-(temp0+20)/1.5/* TAU_[VARNAME_REGISTRY_mainak_gaba_open_local] */));
	*(n->vpra[VARNAME_REGISTRY_mainak_gaba_open_local]) = (*(n->vpra[VARNAME_REGISTRY_mainak_gaba_open_local]) + dt*10*temp0)/(1 + dt*0.16 + dt*10*temp0);
	temp0 = *(n->vpra[VARNAME_REGISTRY_mainak_gaba_open_local]);
	mainak_continuous_synapse_update(n,2,temp0,t+dt/2+0.125 /* GLOBAL_mainak_synapse_delay */);
	*(n->vpra[VARNAME_REGISTRY_mainak_Vs]) = ra_temp2[0];
	*(n->vpra[VARNAME_REGISTRY_mainak_Ca]) = ra_temp2[1];
	*(n->vpra[VARNAME_REGISTRY_mainak_m_CaK]) = ra_temp2[2];
	*(n->vpra[VARNAME_REGISTRY_mainak_m_Ca]) = ra_temp2[3];
	*(n->vpra[VARNAME_REGISTRY_mainak_h_Ca]) = ra_temp2[4];
	*(n->vpra[VARNAME_REGISTRY_mainak_m_K]) = ra_temp2[5];
	temp_time += dt;}
      break;
    default: break;}}}
  tfree(ra);ra=NULL;
  tfree(ra_rhs);ra_rhs=NULL;
  tfree(ra_j);ra_j=NULL;
  tfree(ra_temp1);ra_temp1=NULL;
  tfree(ra_temp2);ra_temp2=NULL;
  tfree(ra_temp3);ra_temp3=NULL;
}

void wilson_synapse_make(struct neuronarray *Nra)
{
  int verbose=0;
  int nt=0,nt2=0,nr=0,nr2=0,tab=0,tab_max=GLOBAL_NTYPES*GLOBAL_INDEXING_sra_LENGTH;
  struct neuron *n=NULL,*s=NULL;
  if (verbose){ printf(" %% [entering wilson_synapse_make]\n");}
  for (nt=0;nt<Nra->ntypes;nt++){ for (nr=0;nr<Nra->lengthra[nt];nr++){
    n = nget(Nra,nt,nr);
    if (n->pulse_synapse!=NULL){ printf(" warning! wilson_synapse!=NULL in wilson_synapse_make\n");}
    n->pulse_synapse = (struct pulse_synapse *) tcalloc(1,sizeof(struct pulse_synapse));
    n->pulse_synapse->pre = llistmake();
    n->pulse_synapse->pos = llistmake();
    n->pulse_synapse->spiketime_llistra = (struct llist **) tcalloc(tab_max,sizeof(struct llist *));
    n->pulse_synapse->value_llistra = (struct llist **) tcalloc(tab_max,sizeof(struct llist *));
    for (tab=0;tab<tab_max;tab++){
      n->pulse_synapse->spiketime_llistra[tab] = llistmake();
      n->pulse_synapse->value_llistra[tab] = llistmake();}}}
  for (nt=0;nt<Nra->ntypes;nt++){ for (nr=0;nr<Nra->lengthra[nt];nr++){
    s = nget(Nra,nt,nr);
    for (nt2=0;nt2<Nra->ntypes;nt2++){ for (nr2=0;nr2<Nra->lengthra[nt2];nr2++){
      n = nget(Nra,nt2,nr2);
      if (verbose){ printf(" %% checking connection between (%d,%d)->(%d,%d)\n",s->type,s->index,n->type,n->index);}
      if (slink(s,n,NULL,NULL)){ 
	if (verbose){ printf(" %% connected\n");}
	litemadd(s->pulse_synapse->pos,n);
	litemadd(n->pulse_synapse->pre,s);}
      else{ if (verbose){ printf(" %% not connected\n");}}}}}}
  if (verbose){ printf(" %% [finishing wilson_synapse_make]\n");}
}

void wilson_input_synapse_update(int verbose,int s_type,struct neuron *n,double *value,double spiketime)
{
  int nv=0,tab=0;
  struct litem *l_spiketime=NULL,*l_value=NULL,*l3=NULL;
  double *temp=NULL,value_temp=0;
  if (verbose){ printf(" %% [entering wilson_input_synapse_update] with n(%d,%d), s_type %d, spiketime %f\n",n->type,n->index,s_type,spiketime); raprintf(value,"double",1,GLOBAL_INDEXING_sra_LENGTH,"value: ");}
  for (nv=0;nv<GLOBAL_INDEXING_sra_LENGTH;nv++){
    tab = s_type+nv*GLOBAL_NTYPES;
    value_temp = value[nv];
    value_temp *= (CS__[s_type+n->type*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES] != 0);
    value_temp *= (GLOBAL_CS_PRETYPE_SCALE_[s_type]*GLOBAL_CS_POSTYPE_SCALE_[n->type]*GLOBAL_CS_SRA_SCALE_[nv] != 0);
    if (verbose>1){ printf(" %% s(%d,?)-->n(%d,%d), nv %d, tab %d, value_temp %f --> %f\n",s_type,n->type,n->index,nv,tab,value[nv],value_temp);}
    if (value_temp != 0){
      if (verbose>2){
	printf(" %% starting neuron (%d,%d), value_temp %f:\n",n->type,n->index,value_temp);
	llistprintf2(n->pulse_synapse->spiketime_llistra[tab]);
	llistprintf2(n->pulse_synapse->value_llistra[tab]);}
      l_spiketime = n->pulse_synapse->spiketime_llistra[tab]->last;
      l_value = n->pulse_synapse->value_llistra[tab]->last;
      while (l_spiketime!=NULL && *(double *)l_spiketime->item > spiketime){ l_spiketime = l_spiketime->parent; l_value = l_value->parent;}
      if (l_spiketime==NULL && l_value==NULL){
	temp = (double *) tcalloc(1,sizeof(double)); *temp = value_temp;
	if (n->pulse_synapse->value_llistra[tab]->length==0){ litemadd(n->pulse_synapse->value_llistra[tab],temp);} 
	else if (n->pulse_synapse->value_llistra[tab]->length>0){ l3 = liteminsertbefore(n->pulse_synapse->value_llistra[tab],n->pulse_synapse->value_llistra[tab]->first); l3->item = temp;}
	temp = (double *) tcalloc(1,sizeof(double)); *temp = spiketime;
	if (n->pulse_synapse->spiketime_llistra[tab]->length==0){ litemadd(n->pulse_synapse->spiketime_llistra[tab],temp);} 
	else if (n->pulse_synapse->spiketime_llistra[tab]->length>0){ l3 = liteminsertbefore(n->pulse_synapse->spiketime_llistra[tab],n->pulse_synapse->spiketime_llistra[tab]->first); l3->item = temp;}}
      else if (l_spiketime!=NULL && l_value!=NULL){
	temp = (double *) tcalloc(1,sizeof(double)); *temp = value_temp;
	l3 = liteminsertafter(n->pulse_synapse->value_llistra[tab],l_value);
	l3->item = temp;	
	temp = (double *) tcalloc(1,sizeof(double)); *temp = spiketime;
	l3 = liteminsertafter(n->pulse_synapse->spiketime_llistra[tab],l_spiketime);
	l3->item = temp;}
      else if ((l_spiketime!=NULL && l_value==NULL) || (l_spiketime==NULL && l_value!=NULL)){ printf(" %% warning! spiketime_llistra[%d] and value_llistra[%d] of different length in wilson_input_synapse_update\n",tab,tab);}
      if (verbose>2){
	printf(" %% finishing neuron (%d,%d):\n",n->type,n->index);
	llistprintf2(n->pulse_synapse->spiketime_llistra[tab]);
	llistprintf2(n->pulse_synapse->value_llistra[tab]);}}
    else /* if (value_temp==0) */{ /* do nothing */ }}
  if (verbose){ printf(" %% [finishing wilson_input_synapse_update]\n");}
}

void wilson_synapse_update(struct neuron *s,double* value,double spiketime)
{
  int verbose=0;//(s->type==0 && s->index==1);
  struct litem *l=NULL;
  struct neuron *n=NULL;
  if (verbose){ printf(" %% [entering wilson_synapse_update] with n(%d,%d), spiketime %f\n",s->type,s->index,spiketime); raprintf(value,"double",1,GLOBAL_INDEXING_sra_LENGTH,"value: ");}
  l = s->pulse_synapse->pos->first;
  while (l!=NULL){
    n = (struct neuron *) l->item;
    wilson_input_synapse_update(verbose,s->type,n,value,spiketime);
    l=l->child;}
  if (verbose){ printf(" %% [finishing wilson_synapse_update]\n");}
/*   int nv=0,tab=0; */
/*   struct litem *l=NULL,*l_spiketime=NULL,*l_value=NULL,*l3=NULL; */
/*   struct neuron *n=NULL; */
/*   double *temp=NULL,value_temp=0; */
/*   if (verbose){ printf(" %% [entering wilson_synapse_update] with n(%d,%d), spiketime %f\n",s->type,s->index,spiketime); raprintf(value,"double",1,GLOBAL_INDEXING_sra_LENGTH,"value: ");} */
/*   l = s->pulse_synapse->pos->first; */
/*   while (l!=NULL){ */
/*     n = (struct neuron *) l->item; */
/*     for (nv=0;nv<GLOBAL_INDEXING_sra_LENGTH;nv++){ */
/*       tab = s->type+nv*GLOBAL_NTYPES; */
/*       value_temp = value[nv]; */
/*       value_temp *= (CS__[s->type+n->type*GLOBAL_NTYPES+nv*GLOBAL_NTYPES*GLOBAL_NTYPES] != 0); */
/*       value_temp *= (GLOBAL_CS_PRETYPE_SCALE_[s->type]*GLOBAL_CS_POSTYPE_SCALE_[n->type]*GLOBAL_CS_SRA_SCALE_[nv] != 0); */
/*       if (verbose>1){ printf(" %% s(%d,%d)-->n(%d,%d), nv %d, tab %d, value_temp %f --> %f\n",s->type,s->index,n->type,n->index,nv,tab,value[nv],value_temp);} */
/*       if (value_temp != 0){ */
/* 	if (verbose>2){ */
/* 	  printf(" %% starting neuron (%d,%d), value_temp %f:\n",n->type,n->index,value_temp); */
/* 	  llistprintf2(n->pulse_synapse->spiketime_llistra[tab]); */
/* 	  llistprintf2(n->pulse_synapse->value_llistra[tab]);} */
/* 	l_spiketime = n->pulse_synapse->spiketime_llistra[tab]->last; */
/* 	l_value = n->pulse_synapse->value_llistra[tab]->last; */
/* 	while (l_spiketime!=NULL && *(double *)l_spiketime->item > spiketime){ l_spiketime = l_spiketime->parent; l_value = l_value->parent;} */
/* 	if (l_spiketime==NULL && l_value==NULL){ */
/* 	  temp = (double *) tcalloc(1,sizeof(double)); *temp = value_temp; */
/* 	  if (n->pulse_synapse->value_llistra[tab]->length==0){ litemadd(n->pulse_synapse->value_llistra[tab],temp);} */
/* 	  else if (n->pulse_synapse->value_llistra[tab]->length>0){ l3 = liteminsertbefore(n->pulse_synapse->value_llistra[tab],n->pulse_synapse->value_llistra[tab]->first); l3->item = temp;} */
/* 	  temp = (double *) tcalloc(1,sizeof(double)); *temp = spiketime; */
/* 	  if (n->pulse_synapse->spiketime_llistra[tab]->length==0){ litemadd(n->pulse_synapse->spiketime_llistra[tab],temp);} */
/* 	  else if (n->pulse_synapse->spiketime_llistra[tab]->length>0){ l3 = liteminsertbefore(n->pulse_synapse->spiketime_llistra[tab],n->pulse_synapse->spiketime_llistra[tab]->first); l3->item = temp;}} */
/* 	else if (l_spiketime!=NULL && l_value!=NULL){ */
/* 	  temp = (double *) tcalloc(1,sizeof(double)); *temp = value_temp; */
/* 	  l3 = liteminsertafter(n->pulse_synapse->value_llistra[tab],l_value); */
/* 	  l3->item = temp; */
/* 	  temp = (double *) tcalloc(1,sizeof(double)); *temp = spiketime; */
/* 	  l3 = liteminsertafter(n->pulse_synapse->spiketime_llistra[tab],l_spiketime); */
/* 	  l3->item = temp;} */
/* 	else if ((l_spiketime!=NULL && l_value==NULL) || (l_spiketime==NULL && l_value!=NULL)){ printf(" %% warning! spiketime_llistra[%d] and value_llistra[%d] of different length in wilson_synapse_update\n",tab,tab);} */
/* 	if (verbose>2){ */
/* 	  printf(" %% finishing neuron (%d,%d):\n",n->type,n->index); */
/* 	  llistprintf2(n->pulse_synapse->spiketime_llistra[tab]); */
/* 	  llistprintf2(n->pulse_synapse->value_llistra[tab]);}} */
/*       else /\* if (value_temp==0) *\/{ /\* do nothing *\/ }} */
/*     l=l->child;} */
/*   if (verbose){ printf(" %% [finishing wilson_synapse_update]\n");} */
}

void wilson_synapse_evolve_helper(int verbose,int length,int skip,struct llist **s,struct llist **v,double t_cutoff,double **dra,int **ira,int *length_output)
{
  /* we assume (*dra) and (*ira) are NULL.
     given input s {1.5, 2.5, 3.0, 4.0} {2.0, 3.0, 3.5} (spiketimes) and v {2, 3, 5, 6} {7, 11, 13} (values), 
     and a cutoff time of 4.0, we return 
     (*dra):
     [ 1.5  2.0  2.5  3.0  3.0 3.5  ]
     [ 2    7    3    5    11  13   ]
     (*ira):
     [ 0    1    0    0    1   1    ]
     with length_output = 6.
  */
  int nl=0,continue_flag=0,no_more_spikes=0,next_spike_index=0,total_number_of_spikes=0;
  double next_spike_time=0;
  struct litem *l=NULL,**lsra=NULL,**lvra=NULL;
  if (verbose){
    printf(" %% [entering wilson_synapse_evolve_helper]\n");
    for (nl=0;nl<length;nl++){ printf(" %% index %d\n",nl*skip); llistprintf2(s[nl*skip]); llistprintf2(v[nl*skip]);}}
  *length_output=0;
  for (nl=0;nl<length;nl++){
    l = s[nl*skip]->first; while (l!=NULL && *(double *)l->item < t_cutoff){ *length_output += 1; l=l->child;}}
  if (*length_output>0){
    if (verbose){ printf(" %% creating lsra, lvra of length %d\n",length);}
    lsra = (struct litem **) tcalloc(length,sizeof(struct litem *));
    lvra = (struct litem **) tcalloc(length,sizeof(struct litem *));
    if (verbose){ printf(" %% creating *dra, *ira of length %d\n",*length_output);}
    (*dra) = (double *) tcalloc(2**length_output,sizeof(double));
    (*ira) = (int *) tcalloc(*length_output,sizeof(int));
    continue_flag=1; total_number_of_spikes = 0;
    for (nl=0;nl<length;nl++){ lsra[nl] = s[nl*skip]->first; lvra[nl] = v[nl*skip]->first;}
    if (verbose){ for (nl=0;nl<length;nl++){ printf(" %% lsra[%d] is %f, lvra[%d] is %f\n",nl,lsra[nl]==NULL ? -1:*(double *)lsra[nl]->item,nl,lvra[nl]==NULL ? -1:*(double *)lvra[nl]->item);}}
    while (continue_flag){
      if (verbose){ printf(" %% currently total_number_of_spikes = %d\n",total_number_of_spikes);}
      no_more_spikes = 1; next_spike_index = -1; next_spike_time = t_cutoff;
      for (nl=0;nl<length;nl++){
	if (lsra[nl]!=NULL && *(double *)lsra[nl]->item < next_spike_time){
	  no_more_spikes = 0; next_spike_index = nl; next_spike_time = *(double *)lsra[nl]->item;}}
      if (verbose){ printf(" %% total_number_of_spikes %d, no_more_spikes %d, next_spike_index %d, next_spike_time %f\n",total_number_of_spikes,no_more_spikes,next_spike_index,next_spike_time);}
      if (!no_more_spikes){
	if (verbose){ printf(" %% found at least one spike, updating *dra, *ira\n");}
	(*dra)[total_number_of_spikes + 0**length_output] = *(double *)lsra[next_spike_index]->item;
	(*dra)[total_number_of_spikes + 1**length_output] = *(double *)lvra[next_spike_index]->item;
	(*ira)[total_number_of_spikes] = next_spike_index;
	total_number_of_spikes += 1;
	if (verbose){ printf(" %% killing s[%d]->first, v[%d]->first\n",next_spike_index*skip,next_spike_index*skip);}
	llistkillfirst2(s[next_spike_index*skip]),llistkillfirst2(v[next_spike_index*skip]);
	lsra[next_spike_index] = s[next_spike_index*skip]->first; lvra[next_spike_index] = v[next_spike_index*skip]->first;
	if (verbose){ printf(" %% resetting lsra[%d] to %f, lvra[%d] to %f\n",next_spike_index,lsra[next_spike_index]==NULL ? -1:*(double *)lsra[next_spike_index]->item,next_spike_index,lvra[next_spike_index]==NULL ? -1:*(double *)lvra[next_spike_index]->item);}}
      if (verbose){ printf(" %% now total_number_of_spikes %d\n",total_number_of_spikes); raprintf((*dra),"double",*length_output,2," %% (*dra): ");raprintf((*ira),"int",*length_output,1," %% (*ira): ");}
      continue_flag = !no_more_spikes;}
    tfree(lsra);lsra=NULL; tfree(lvra);lvra=NULL;}
  if (verbose){
    printf(" %% [finishing wilson_synapse_evolve_helper] total_number_of_spikes = %d\n",total_number_of_spikes);}
}

void wilson_synapse_evolve(struct neuron *n,double t,double DT)
{
  int verbose=0;//(n->type==0 && n->index==0);
  int nt=0,nv=0,tab=0,tab2=0,ns=0;
  double temp_time=0,temp_spiketime=0,temp_scale=0;
  double *dra=NULL;
  int *ira=NULL,nspikes=0;
  if (verbose){ printf(" %% [entering wilson_synapse_evolve] with n(%d,%d),time %f->%f\n",n->type,n->index,t,t+DT);}
  nv = GLOBAL_INDEXING_REFILE_sra[VARNAME_REGISTRY_wilson_nAch];
  wilson_synapse_evolve_helper(verbose,GLOBAL_NTYPES,1,&(n->pulse_synapse->spiketime_llistra[0+nv*GLOBAL_NTYPES]),&(n->pulse_synapse->value_llistra[0+nv*GLOBAL_NTYPES]),t+DT,&dra,&ira,&nspikes);
  if (verbose>2){
    for (nt=0;nt<GLOBAL_NTYPES;nt++){
      printf(" %% examining nAch llist for input type %d\n",nt); 
      tab = nt+nv*GLOBAL_NTYPES; 
      llistprintf2(n->pulse_synapse->spiketime_llistra[tab]);
      llistprintf2(n->pulse_synapse->value_llistra[tab]);}
    printf(" %% current nAch given by %f\n",*(n->vpra[VARNAME_REGISTRY_wilson_nAch]));
    if (nspikes>0){
      raprintf(dra,"double",nspikes,2," %% dra: ");
      raprintf(ira,"int",nspikes,1," %% ira: ");}}
  temp_time=t;
  for (ns=0;ns<nspikes;ns++){
    temp_spiketime = dra[ns+0*nspikes];
    tab2 = ira[ns] + n->type*GLOBAL_NTYPES + nv*GLOBAL_NTYPES*GLOBAL_NTYPES;
    temp_scale = CS__[tab2]*GLOBAL_CS_PRETYPE_SCALE_[ira[ns]]*GLOBAL_CS_POSTYPE_SCALE_[n->type]*GLOBAL_CS_SRA_SCALE_[nv];
    *(n->vpra[VARNAME_REGISTRY_wilson_nAch]) *= exp(-(temp_spiketime-temp_time)/TAU_[VARNAME_REGISTRY_wilson_nAch]);
    *(n->vpra[VARNAME_REGISTRY_wilson_nAch]) += temp_scale*dra[ns+1*nspikes];
    temp_time = temp_spiketime;}
  *(n->vpra[VARNAME_REGISTRY_wilson_nAch]) *= exp(-(t+DT-temp_time)/TAU_[VARNAME_REGISTRY_wilson_nAch]);
  if (dra!=NULL){ tfree(dra);dra=NULL;} if (ira!=NULL){ tfree(ira);ira=NULL;} nspikes=0;
  if (verbose>2){
    for (nt=0;nt<GLOBAL_NTYPES;nt++){
      printf(" %% now examining nAch llist for input type %d\n",nt); 
      tab = nt+nv*GLOBAL_NTYPES; 
      llistprintf2(n->pulse_synapse->spiketime_llistra[tab]);
      llistprintf2(n->pulse_synapse->value_llistra[tab]);}
    printf(" %% current nAch given by %f\n",*(n->vpra[VARNAME_REGISTRY_wilson_nAch]));}
  nv = GLOBAL_INDEXING_REFILE_sra[VARNAME_REGISTRY_wilson_gabaA];
  wilson_synapse_evolve_helper(verbose,GLOBAL_NTYPES,1,&(n->pulse_synapse->spiketime_llistra[0+nv*GLOBAL_NTYPES]),&(n->pulse_synapse->value_llistra[0+nv*GLOBAL_NTYPES]),t+DT,&dra,&ira,&nspikes);
  if (verbose>2){
    for (nt=0;nt<GLOBAL_NTYPES;nt++){
      printf(" %% examining gabaA llist for input type %d\n",nt); 
      tab = nt+nv*GLOBAL_NTYPES; 
      llistprintf2(n->pulse_synapse->spiketime_llistra[tab]);
      llistprintf2(n->pulse_synapse->value_llistra[tab]);}
    printf(" %% current gabaA given by %f\n",*(n->vpra[VARNAME_REGISTRY_wilson_gabaA]));
    if (nspikes>0){
      raprintf(dra,"double",nspikes,2," %% dra: ");
      raprintf(ira,"int",nspikes,1," %% ira: ");}}
  temp_time=t;
  for (ns=0;ns<nspikes;ns++){
    temp_spiketime = dra[ns+0*nspikes];
    tab2 = ira[ns] + n->type*GLOBAL_NTYPES + nv*GLOBAL_NTYPES*GLOBAL_NTYPES;
    temp_scale = CS__[tab2]*GLOBAL_CS_PRETYPE_SCALE_[ira[ns]]*GLOBAL_CS_POSTYPE_SCALE_[n->type]*GLOBAL_CS_SRA_SCALE_[nv];
    *(n->vpra[VARNAME_REGISTRY_wilson_gabaA]) *= exp(-(temp_spiketime-temp_time)/TAU_[VARNAME_REGISTRY_wilson_gabaA]);
    *(n->vpra[VARNAME_REGISTRY_wilson_gabaA]) += temp_scale*dra[ns+1*nspikes];
    temp_time = temp_spiketime;}
  *(n->vpra[VARNAME_REGISTRY_wilson_gabaA]) *= exp(-(t+DT-temp_time)/TAU_[VARNAME_REGISTRY_wilson_gabaA]);
  if (dra!=NULL){ tfree(dra);dra=NULL;} if (ira!=NULL){ tfree(ira);ira=NULL;} nspikes=0;
  if (verbose>2){
    for (nt=0;nt<GLOBAL_NTYPES;nt++){
      printf(" %% now examining gabaA llist for input type %d\n",nt); 
      tab = nt+nv*GLOBAL_NTYPES; 
      llistprintf2(n->pulse_synapse->spiketime_llistra[tab]);
      llistprintf2(n->pulse_synapse->value_llistra[tab]);}
    printf(" %% current gabaA given by %f\n",*(n->vpra[VARNAME_REGISTRY_wilson_gabaA]));}
  nv = GLOBAL_INDEXING_REFILE_sra[VARNAME_REGISTRY_wilson_gabaB];
  wilson_synapse_evolve_helper(verbose,GLOBAL_NTYPES,1,&(n->pulse_synapse->spiketime_llistra[0+nv*GLOBAL_NTYPES]),&(n->pulse_synapse->value_llistra[0+nv*GLOBAL_NTYPES]),t+DT,&dra,&ira,&nspikes);
  if (verbose>2){
    for (nt=0;nt<GLOBAL_NTYPES;nt++){
      printf(" %% examining gabaBp llist for input type %d\n",nt); 
      tab = nt+nv*GLOBAL_NTYPES; 
      llistprintf2(n->pulse_synapse->spiketime_llistra[tab]);
      llistprintf2(n->pulse_synapse->value_llistra[tab]);}
    printf(" %% current gabaBp given by %f\n",*(n->vpra[VARNAME_REGISTRY_wilson_gabaBp]));
    if (nspikes>0){
      raprintf(dra,"double",nspikes,2," %% dra: ");
      raprintf(ira,"int",nspikes,1," %% ira: ");}}
  temp_time=t;
  for (ns=0;ns<nspikes;ns++){
    temp_spiketime = dra[ns+0*nspikes];
    tab2 = ira[ns] + n->type*GLOBAL_NTYPES + nv*GLOBAL_NTYPES*GLOBAL_NTYPES;
    temp_scale = CS__[tab2]*GLOBAL_CS_PRETYPE_SCALE_[ira[ns]]*GLOBAL_CS_POSTYPE_SCALE_[n->type]*GLOBAL_CS_SRA_SCALE_[nv];
    *(n->vpra[VARNAME_REGISTRY_wilson_gabaB]) = *(n->vpra[VARNAME_REGISTRY_wilson_gabaB])*exp(-(temp_spiketime-temp_time)/TAU_[VARNAME_REGISTRY_wilson_gabaB]) + *(n->vpra[VARNAME_REGISTRY_wilson_gabaBp])/(1.0/TAU_[VARNAME_REGISTRY_wilson_gabaB] - 1.0/TAU_[VARNAME_REGISTRY_wilson_gabaBp])*(exp(-(temp_spiketime-temp_time)/TAU_[VARNAME_REGISTRY_wilson_gabaBp]) - exp(-(temp_spiketime-temp_time)/TAU_[VARNAME_REGISTRY_wilson_gabaB]));
    *(n->vpra[VARNAME_REGISTRY_wilson_gabaBp]) *= exp(-(temp_spiketime-temp_time)/TAU_[VARNAME_REGISTRY_wilson_gabaBp]);
    *(n->vpra[VARNAME_REGISTRY_wilson_gabaBp]) += temp_scale*dra[ns+1*nspikes];
    temp_time = temp_spiketime;}
  *(n->vpra[VARNAME_REGISTRY_wilson_gabaB]) = *(n->vpra[VARNAME_REGISTRY_wilson_gabaB])*exp(-(t+DT-temp_time)/TAU_[VARNAME_REGISTRY_wilson_gabaB]) + *(n->vpra[VARNAME_REGISTRY_wilson_gabaBp])/(1.0/TAU_[VARNAME_REGISTRY_wilson_gabaB]-1.0/TAU_[VARNAME_REGISTRY_wilson_gabaBp])*(exp(-(t+DT-temp_time)/TAU_[VARNAME_REGISTRY_wilson_gabaBp]) - exp(-(t+DT-temp_time)/TAU_[VARNAME_REGISTRY_wilson_gabaB]));
  *(n->vpra[VARNAME_REGISTRY_wilson_gabaBp]) *= exp(-(t+DT-temp_time)/TAU_[VARNAME_REGISTRY_wilson_gabaBp]);
  if (dra!=NULL){ tfree(dra);dra=NULL;} if (ira!=NULL){ tfree(ira);ira=NULL;} nspikes=0;
  if (verbose>2){
    for (nt=0;nt<GLOBAL_NTYPES;nt++){
      printf(" %% now examining gabaBp llist for input type %d\n",nt); 
      tab = nt+nv*GLOBAL_NTYPES; 
      llistprintf2(n->pulse_synapse->spiketime_llistra[tab]);
      llistprintf2(n->pulse_synapse->value_llistra[tab]);}
    printf(" %% current gabaBp given by %f\n",*(n->vpra[VARNAME_REGISTRY_wilson_gabaBp]));}
  if (verbose){ printf(" %% [finishing wilson_synapse_evolve]\n");}
}

void wilson_evolve_helper_rhs_j_inv(double *ra_rhs,double *ra_j)
{
  /*
    assumes that ra_rhs has component ordering
    |  V  | m_Na| h_Na| m_K |
    assumes that ra_j has the structure
    |  V  | m_Na| h_Na| m_K |
    +-----+-----+-----+-----+
    |  0  |  1  |  2  |  3  |
    +-----+-----+-----+-----+
    |  4  |  5  |     |     |
    +-----+-----+-----+-----+
    |  6  |     |  7  |     |
    +-----+-----+-----+-----+
    |  8  |     |     |  9  |
    +-----+-----+-----+-----+
   */
  int verbose=0;
  double temp=0;
  double *ra1=NULL,*ra2=NULL;
  if (verbose){ 
    ra1=(double *) tcalloc(4,sizeof(double)); raplusequals(ra1,4,ra_rhs);
    ra2=(double *) tcalloc(10,sizeof(double)); raplusequals(ra2,10,ra_j);}
  temp = -ra_j[3]/ra_j[9]; ra_j[0] += temp*ra_j[8]; ra_rhs[0] += temp*ra_rhs[3];
  temp = -ra_j[2]/ra_j[7]; ra_j[0] += temp*ra_j[6]; ra_rhs[0] += temp*ra_rhs[2];
  temp = -ra_j[1]/ra_j[5]; ra_j[0] += temp*ra_j[4]; ra_rhs[0] += temp*ra_rhs[1];  
  temp = -ra_j[4]/ra_j[0];  ra_rhs[1] += temp*ra_rhs[0];
  temp = -ra_j[6]/ra_j[0];  ra_rhs[2] += temp*ra_rhs[0];
  temp = -ra_j[8]/ra_j[0]; ra_rhs[3] += temp*ra_rhs[0];
  ra_rhs[0] /= ra_j[0];
  ra_rhs[1] /= ra_j[5];
  ra_rhs[2] /= ra_j[7];
  ra_rhs[3] /= ra_j[9];
  if (verbose){
    printf("ra2*ra_rhs =  %f, %f, %f, %f \n",ra2[0]*ra_rhs[0]+ra2[1]*ra_rhs[1]+ra2[2]*ra_rhs[2]+ra2[3]*ra_rhs[3],ra2[4]*ra_rhs[0]+ra2[5]*ra_rhs[1],ra2[6]*ra_rhs[0]+ra2[7]*ra_rhs[2],ra2[8]*ra_rhs[0]+ra2[9]*ra_rhs[3]);
    raprintf(ra1,"double",1,4,"ra1 = ");
    tfree(ra1);ra1=NULL;tfree(ra2);ra2=NULL;}
}

void wilson_evolve_helper_rhs(double *ra,double **vpra_temp,double *ra_rhs,double *ra_j,int type)
{
  /* presumes that ra has component ordering
     V,m_Na,h_Na,m_K
     ra_rhs and ra_j are similarly ordered 
     type==TYPENAME_REGISTRY_wilson_ORN are special cased not to receive inhibition */
  double V=0,m_Na=0,h_Na=0,m_K=0,nAch=0,gabaA=0,gabaB=0,s_ORN=0;
  double conductance_leak=0,voltage_leak=0;
  double alpha=0,alpha_V=0,beta=0,beta_V=0;
  double temp0=0,temp1=0,temp1p=0,temp3=0,temp3p=0;
  double I_LEAK_S=0;
  double I_LEAK_S_V=0;
  double I_Na=0;
  double I_Na_V=0;
  double I_Na_m_Na=0;
  double I_Na_h_Na=0;
  double I_K=0;
  double I_K_V=0;
  double I_K_m_K=0;
  double I_nAch=0;
  double I_nAch_V=0;
  double I_gabaA=0;
  double I_gabaA_V=0;
  double I_gabaB=0;
  double I_gabaB_V=0;
  double I_orn=0;
  double I_orn_V=0;
  double rhs_m_Na=0;
  double rhs_m_Na_V=0;
  double rhs_m_Na_m_Na=0;
  double rhs_h_Na=0;
  double rhs_h_Na_V=0;
  double rhs_h_Na_h_Na=0;
  double rhs_m_K=0;
  double rhs_m_K_V=0;
  double rhs_m_K_m_K=0;
  V=ra[0];
  m_Na=ra[1];
  h_Na=ra[2];
  m_K=ra[3];
  switch (type){ default: voltage_leak = -64.387; conductance_leak = 0.3; break;}
  I_LEAK_S = -conductance_leak*(V - voltage_leak);
  I_LEAK_S_V = -conductance_leak;
  I_Na = -/* CONDUCTANCE_[VARNAME_REGISTRY_wilson_m_Na] */120*pow(m_Na,3)*pow(h_Na,1)*(V-/* VOLTAGE_[VARNAME_REGISTRY_wilson_m_Na] */ 40);
  I_Na_V = -/* CONDUCTANCE_[VARNAME_REGISTRY_wilson_m_Na] */120*pow(m_Na,3)*pow(h_Na,1);
  I_Na_m_Na = -/* CONDUCTANCE_[VARNAME_REGISTRY_wilson_m_Na] */120*3*pow(m_Na,2)*pow(h_Na,1)*(V-/* VOLTAGE_[VARNAME_REGISTRY_wilson_m_Na] */ 40);
  I_Na_h_Na = -/* CONDUCTANCE_[VARNAME_REGISTRY_wilson_m_Na] */120*pow(m_Na,3)*(V-/* VOLTAGE_[VARNAME_REGISTRY_wilson_m_Na] */ 40);
  I_K = -/* CONDUCTANCE_[VARNAME_REGISTRY_wilson_m_K] */3.6*pow(m_K,4)*(V-/* VOLTAGE_[VARNAME_REGISTRY_wilson_m_K] */ -87);
  I_K_V = -/* CONDUCTANCE_[VARNAME_REGISTRY_wilson_m_K] */3.6*pow(m_K,4);
  I_K_m_K = -/* CONDUCTANCE_[VARNAME_REGISTRY_wilson_m_K] */3.6*4*pow(m_K,3)*(V-/* VOLTAGE_[VARNAME_REGISTRY_wilson_m_K] */ -87);
  temp0 = V+35.0; temp1 = -exp(-0.1*temp0); temp3p = -0.1*temp1; temp3 = 1.0+temp1;
  alpha = 0.1*temp0/temp3; alpha_V = 0.1*(1.0/temp3 - temp0/pow(temp3,2)*temp3p);
  temp0 = V+60.0; temp1 = 4.0*exp(-temp0/18.0); temp1p = -1.0/18.0*temp1;
  beta = temp1; beta_V = temp1p;
  rhs_m_Na = alpha*(1.0 - m_Na) - beta*m_Na;
  rhs_m_Na_V = alpha_V*(1.0 - m_Na) - beta_V*m_Na;
  rhs_m_Na_m_Na = -alpha-beta;
  temp0 = V+60.0; temp1 = 0.07*exp(-temp0/20.0); temp1p = -1.0/20.0*temp1;
  alpha = temp1; alpha_V = temp1p;
  temp0 = V+30.0; temp1 = exp(-0.1*temp0); temp3p = -0.1*temp1; temp3 = 1.0+temp1;
  beta = 1.0/temp3; beta_V = -1.0/pow(temp3,2)*temp3p;
  rhs_h_Na = alpha*(1.0 - h_Na) - beta*h_Na;
  rhs_h_Na_V = alpha_V*(1.0 - h_Na) - beta_V*h_Na;
  rhs_h_Na_h_Na = -alpha-beta;
  temp0 = V+50.0; temp1 = -exp(-0.1*temp0); temp3p = -0.1*temp1; temp3 = 1.0+temp1;
  alpha = 0.01*temp0/temp3; alpha_V = 0.01*(1.0/temp3 - temp0/pow(temp3,2)*temp3p);
  temp0 = V+60.0; temp1 = 0.125*exp(-temp0/80.0); temp1p = -1.0/80.0*temp1;
  beta = temp1; beta_V = temp1p;
  rhs_m_K = alpha*(1.0 - m_K) - beta*m_K;
  rhs_m_K_V = alpha_V*(1.0 - m_K) - beta_V*m_K;
  rhs_m_K_m_K = -alpha-beta;
  nAch = *(vpra_temp[VARNAME_REGISTRY_wilson_nAch]);
  I_nAch = -nAch*(V-VOLTAGE_[VARNAME_REGISTRY_wilson_nAch]);
  I_nAch_V = nAch;
  temp0 = 0;
  switch (type){ case TYPENAME_REGISTRY_wilson_ORN: gabaA = 0; break; default: gabaA = *(vpra_temp[VARNAME_REGISTRY_wilson_gabaA]); break;}
  I_gabaA = -gabaA*(V-VOLTAGE_[VARNAME_REGISTRY_wilson_gabaA]);
  I_gabaA_V = gabaA;
  temp0 = 0;
  switch (type){ case TYPENAME_REGISTRY_wilson_ORN: gabaB = 0; break; default: gabaB = *(vpra_temp[VARNAME_REGISTRY_wilson_gabaB]); break;}
  I_gabaB = -gabaB*(V-VOLTAGE_[VARNAME_REGISTRY_wilson_gabaB]);
  I_gabaB_V = gabaB;
  temp0 = 0;
  s_ORN = *(vpra_temp[VARNAME_REGISTRY_wilson_s_ORN]);
  I_orn = s_ORN;
  I_orn_V = 0;
  ra_rhs[0] = I_LEAK_S + I_Na + I_K + I_nAch + I_gabaA + I_gabaB + I_orn + CURRENT_INJECTION_S;
  ra_rhs[1] = rhs_m_Na;
  ra_rhs[2] = rhs_h_Na;
  ra_rhs[3] = rhs_m_K;
  ra_j[0] = I_LEAK_S_V + I_Na_V + I_K_V + I_nAch_V + I_gabaA_V + I_gabaB_V + I_orn_V;
  ra_j[1] = I_Na_m_Na;
  ra_j[2] = I_Na_h_Na;
  ra_j[3] = I_K_m_K;
  ra_j[4] = rhs_m_Na_V;
  ra_j[5] = rhs_m_Na_m_Na;
  ra_j[6] = rhs_h_Na_V;
  ra_j[7] = rhs_h_Na_h_Na;
  ra_j[8] = rhs_m_K_V;
  ra_j[9] = rhs_m_K_m_K;
}

void wilson_evolve(struct neuronarray *Nra,double t,double DT)
{
  /* type==TYPENAME_REGISTRY_wilson_ORN are special cased to experience inhibition on their axons, 
     currently a hack job, must fix later */
  int verbose=GLOBAL_verbose;
  int nt=0,nr=0,nv=0,ns=0,iteration=0,iteration_max=0;
  double dt=0,norm=0,tolerance=0.000001;
  double temp_time=0,temp_dt=0,temp_spiketime=0,*temp=NULL;
  double *ra=NULL,*ra_rhs=NULL,*ra_j=NULL,*ra_temp1=NULL,*ra_temp2=NULL,*ra_temp3=NULL;
  struct neuron *n=NULL;
  int rhs_size=4,j_size=10;
  double temp_gabaA_1=0,temp_gabaB_1=0,temp_gabaA_2=0,temp_gabaB_2=0,open_channel_scale=1;
  if (verbose){ printf(" %% [entering wilson_evolve] t %f\n",t);}
  ra = (double *) tcalloc(rhs_size,sizeof(double));
  ra_rhs = (double *) tcalloc(rhs_size,sizeof(double));
  ra_j = (double *) tcalloc(j_size,sizeof(double));
  ra_temp1 = (double *) tcalloc(j_size,sizeof(double));
  ra_temp2 = (double *) tcalloc(j_size,sizeof(double));
  ra_temp3 = (double *) tcalloc(j_size,sizeof(double));
  for (nt=0;nt<Nra->ntypes;nt++){ for (nr=0;nr<Nra->lengthra[nt];nr++){
    n=nget(Nra,nt,nr);
    *(n->vpra[VARNAME_REGISTRY_wilson_inputrate]) = n->inputrate;
    if (verbose>1){ printf(" %% n(%d,%d), %d microsteps\n",nt,nr,n->microstep);}
    temp_time = t;
    n->microstep = 1;
    dt = DT/n->microstep;
    for (ns=0;ns<n->microstep;ns++){
      *(n->vpra[VARNAME_REGISTRY_wilson_s_ORN]) *= exp(-dt/TAU_[VARNAME_REGISTRY_wilson_s_ORN]);
      if (n->type==TYPENAME_REGISTRY_wilson_ORN){ temp_gabaA_1 = *(n->vpra[VARNAME_REGISTRY_wilson_gabaA]); temp_gabaB_1 = *(n->vpra[VARNAME_REGISTRY_wilson_gabaB]);}
      wilson_synapse_evolve(n,temp_time,dt);
      if (n->type==TYPENAME_REGISTRY_wilson_ORN){ temp_gabaA_2 = *(n->vpra[VARNAME_REGISTRY_wilson_gabaA]); temp_gabaB_2 = *(n->vpra[VARNAME_REGISTRY_wilson_gabaB]);}
      ra[0] = *(n->vpra[VARNAME_REGISTRY_wilson_Vs]);
      ra[1] = *(n->vpra[VARNAME_REGISTRY_wilson_m_Na]);
      ra[2] = *(n->vpra[VARNAME_REGISTRY_wilson_h_Na]);
      ra[3] = *(n->vpra[VARNAME_REGISTRY_wilson_m_K]);
      raplugin(ra_temp2,4,1,ra,4,1,0,0);
      iteration=0; iteration_max=10;
      do{
	if (verbose>3){ printf(" %% iteration %d\n",iteration);}
	wilson_evolve_helper_rhs(ra_temp2,n->vpra,ra_rhs,ra_j,n->type);
	ratimesequals(ra_j,10,-dt);
	ra_j[0]  += 1; ra_j[5] += 1; ra_j[7]  += 1; ra_j[9] += 1;
	ra_temp1[0] = dt*ra_rhs[0]+ra[0]-ra_temp2[0];
	ra_temp1[1] = dt*ra_rhs[1]+ra[1]-ra_temp2[1];
	ra_temp1[2] = dt*ra_rhs[2]+ra[2]-ra_temp2[2];
	ra_temp1[3] = dt*ra_rhs[3]+ra[3]-ra_temp2[3];
	wilson_evolve_helper_rhs_j_inv(ra_temp1,ra_j);
	norm = ra_norm(ra_temp1,4);
	ra_temp2[0] += ra_temp1[0];
	ra_temp2[1] += ra_temp1[1];
	ra_temp2[2] += ra_temp1[2];
	ra_temp2[3] += ra_temp1[3];
	if (verbose>3){ printf(" %% iteration %d norm %f\n",iteration,norm);}
	iteration += 1;}
      while (norm>tolerance && iteration<iteration_max);
      if (iteration>=iteration_max){ if (verbose){ printf(" %% warning! newton not converging in wilson_evolve\n");}}
      if (*(n->vpra[VARNAME_REGISTRY_wilson_Vs]) > 0 /* wilson_threshold */&& ra_temp2[0] <= 0 /* wilson_threshold */){
	temp_dt = linerootfinder(-*(n->vpra[VARNAME_REGISTRY_wilson_Vs]),-ra_temp2[0],-0/* wilson_threshold */,dt); 
	temp_spiketime = temp_time + temp_dt;
	*(n->vpra[VARNAME_REGISTRY_wilson_nAch_local]) *= exp(-temp_dt/TAU_[VARNAME_REGISTRY_wilson_nAch_local]);
	*(n->vpra[VARNAME_REGISTRY_wilson_gabaA_local]) *= exp(-temp_dt/TAU_[VARNAME_REGISTRY_wilson_gabaA_local]);
	*(n->vpra[VARNAME_REGISTRY_wilson_gabaBp_local]) *= exp(-temp_dt/TAU_[VARNAME_REGISTRY_wilson_gabaBp_local]);
	*(n->vpra[VARNAME_REGISTRY_wilson_vesicle_depletion]) *= exp(-temp_dt/TAU_[VARNAME_REGISTRY_wilson_vesicle_depletion]);
	*(n->vpra[VARNAME_REGISTRY_wilson_spike_flag]) *= exp(-temp_dt/TAU_[VARNAME_REGISTRY_wilson_spike_flag]);
	temp = (double *) tcalloc(GLOBAL_INDEXING_sra_LENGTH,sizeof(double));
	open_channel_scale = 1;
	if (n->type==TYPENAME_REGISTRY_wilson_ORN){ 
	  open_channel_scale *= (1-*(n->vpra[VARNAME_REGISTRY_wilson_vesicle_depletion]))*exp(-((dt-temp_dt)/dt*(temp_gabaA_1+temp_gabaB_1) + temp_dt/dt*(temp_gabaA_2+temp_gabaB_2)));
	  *(n->vpra[VARNAME_REGISTRY_wilson_vesicle_depletion]) += (*(n->vpra[VARNAME_REGISTRY_wilson_vesicle_depletion]) - 1.0)*(exp(-abs(GLOBAL_VESICLE_DEPLETION)) - 1.0);}
	nv = GLOBAL_INDEXING_REFILE_sra[VARNAME_REGISTRY_wilson_nAch];
	temp[nv] = (*(n->vpra[VARNAME_REGISTRY_wilson_nAch_local]) - 1.0)*(exp(-1.5*open_channel_scale) - 1.0);
	nv = GLOBAL_INDEXING_REFILE_sra[VARNAME_REGISTRY_wilson_gabaA];
	temp[nv] = (*(n->vpra[VARNAME_REGISTRY_wilson_gabaA_local]) - 1.0)*(exp(-1.5*open_channel_scale) - 1.0);
	nv = GLOBAL_INDEXING_REFILE_sra[VARNAME_REGISTRY_wilson_gabaB];
	temp[nv] = (*(n->vpra[VARNAME_REGISTRY_wilson_gabaBp_local]) - 1.0)*(exp(-1.5*open_channel_scale) - 1.0);
	nv = GLOBAL_INDEXING_REFILE_sra[VARNAME_REGISTRY_wilson_s_ORN];
	temp[nv] = 0;
	wilson_synapse_update(n,temp,temp_spiketime+0.125 /* GLOBAL_wilson_synapse_delay */);
	*(n->vpra[VARNAME_REGISTRY_wilson_spike_flag]) += 1;
	nv = GLOBAL_INDEXING_REFILE_sra[VARNAME_REGISTRY_wilson_nAch];
	*(n->vpra[VARNAME_REGISTRY_wilson_nAch_local]) += temp[nv];
	nv = GLOBAL_INDEXING_REFILE_sra[VARNAME_REGISTRY_wilson_gabaA];
	*(n->vpra[VARNAME_REGISTRY_wilson_gabaA_local]) += temp[nv];
	nv = GLOBAL_INDEXING_REFILE_sra[VARNAME_REGISTRY_wilson_gabaB];
	*(n->vpra[VARNAME_REGISTRY_wilson_gabaBp_local]) += temp[nv];
	*(n->vpra[VARNAME_REGISTRY_wilson_nAch_local]) *= exp(-(dt-temp_dt)/TAU_[VARNAME_REGISTRY_wilson_nAch_local]);
	*(n->vpra[VARNAME_REGISTRY_wilson_gabaA_local]) *= exp(-(dt-temp_dt)/TAU_[VARNAME_REGISTRY_wilson_gabaA_local]);
	*(n->vpra[VARNAME_REGISTRY_wilson_gabaBp_local]) *= exp(-(dt-temp_dt)/TAU_[VARNAME_REGISTRY_wilson_gabaBp_local]);
	*(n->vpra[VARNAME_REGISTRY_wilson_vesicle_depletion]) *= exp(-(dt-temp_dt)/TAU_[VARNAME_REGISTRY_wilson_vesicle_depletion]);
	*(n->vpra[VARNAME_REGISTRY_wilson_spike_flag]) *= exp(-(dt-temp_dt)/TAU_[VARNAME_REGISTRY_wilson_spike_flag]);
	tfree(temp);temp=NULL;}
      else /* if (*(n->vpra[VARNAME_REGISTRY_wilson_Vs]) <= 0 || ra_temp2[0] > 0) */{
	*(n->vpra[VARNAME_REGISTRY_wilson_nAch_local]) *= exp(-dt/TAU_[VARNAME_REGISTRY_wilson_nAch_local]);
	*(n->vpra[VARNAME_REGISTRY_wilson_gabaA_local]) *= exp(-dt/TAU_[VARNAME_REGISTRY_wilson_gabaA_local]);
	*(n->vpra[VARNAME_REGISTRY_wilson_gabaBp_local]) *= exp(-dt/TAU_[VARNAME_REGISTRY_wilson_gabaBp_local]);
	*(n->vpra[VARNAME_REGISTRY_wilson_vesicle_depletion]) *= exp(-dt/TAU_[VARNAME_REGISTRY_wilson_vesicle_depletion]);
	*(n->vpra[VARNAME_REGISTRY_wilson_spike_flag]) *= exp(-dt/TAU_[VARNAME_REGISTRY_wilson_spike_flag]);}
      *(n->vpra[VARNAME_REGISTRY_wilson_Vs]) = ra_temp2[0];
      *(n->vpra[VARNAME_REGISTRY_wilson_m_Na]) = ra_temp2[1];
      *(n->vpra[VARNAME_REGISTRY_wilson_h_Na]) = ra_temp2[2];
      *(n->vpra[VARNAME_REGISTRY_wilson_m_K]) = ra_temp2[3];
      *(n->vpra[VARNAME_REGISTRY_wilson_Vs_clip]) = minimum(GLOBAL_POWER_maxra_[VARNAME_REGISTRY_wilson_Vs_clip],maximum(GLOBAL_POWER_minra_[VARNAME_REGISTRY_wilson_Vs_clip],*(n->vpra[VARNAME_REGISTRY_wilson_Vs])));
      temp_time += dt;}}}
  tfree(ra);ra=NULL;
  tfree(ra_rhs);ra_rhs=NULL;
  tfree(ra_j);ra_j=NULL;
  tfree(ra_temp1);ra_temp1=NULL;
  tfree(ra_temp2);ra_temp2=NULL;
  tfree(ra_temp3);ra_temp3=NULL;
}

void snx_statematrix_1(int Vs_state,int total_Vs_states,double input,double *ra,double *ra_prime)
{
  /* attempting simple current-based i&f 
     ra holds the output probabilities P(current_state --> state_i) 
     it is assumed that ra is already initialized with size total_Vs_states */
  int verbose=0;
  double S=0.5;
  double voltage_spread_1 = 0.025, voltage_spread_2 = 0.025, tau_v = 0.1, v=0, v_temp=0, x0=0;
  double *ra_y1=NULL,*ra_y2=NULL,sum_y1=0,sum_y2=0;
  int index=0;
  double input_rectified=0;
  input_rectified = maximum(0,minimum(1,input));
  if (Vs_state>=total_Vs_states-1 || Vs_state<0){
    ra[0]=1; for (index=1;index<total_Vs_states;index++){ ra[index]=0;}
    if (ra_prime!=NULL){ for (index=1;index<total_Vs_states;index++){ ra_prime[index]=0;}}}
  else /* if (Vs_state>=0 && Vs_state<total_Vs_states-1) */{ 
    v = (double)Vs_state/(double)(total_Vs_states-1);
    x0 = v*exp(-tau_v);
    ra_y1=(double *) tcalloc(total_Vs_states,sizeof(double));sum_y1=0; ra_y2=(double *) tcalloc(total_Vs_states,sizeof(double));sum_y2=0;
    for (index=0;index<total_Vs_states;index++){
      v_temp = (double)index/(double)(total_Vs_states-1);
      ra_y1[index] = exp(-pow(v_temp-x0,2)/2/voltage_spread_1); sum_y1 += ra_y1[index];
      ra_y2[index] = exp(-pow(v_temp-(x0+S),2)/2/voltage_spread_2); sum_y2 += ra_y2[index];}
    if (verbose){ raprintf(ra_y1,"double",1,total_Vs_states,"ra_y1, raw: "); raprintf(ra_y2,"double",1,total_Vs_states,"ra_y2, raw: ");}
    ratimesequals(ra_y1,total_Vs_states,1.0/sum_y1); ratimesequals(ra_y2,total_Vs_states,1.0/sum_y2);
    if (verbose){ raprintf(ra_y1,"double",1,total_Vs_states,"ra_y1: ");raprintf(ra_y2,"double",1,total_Vs_states,"ra_y2: ");}
    for (index=0;index<total_Vs_states;index++){
      ra[index] = (1-input_rectified)*ra_y1[index] + input_rectified*ra_y2[index];}
    if (verbose){ printf(" input_rectified %f\n",input_rectified); raprintf(ra,"double",1,total_Vs_states,"ra: ");}
    if (ra_prime!=NULL){ 
      for (index=0;index<total_Vs_states;index++){ ra_prime[index] = (input_rectified==input ? ra_y2[index]-ra_y1[index] : 0);}}
    tfree(ra_y1);ra_y1=NULL;tfree(ra_y2);ra_y2=NULL;}
}

void snx_statematrix_0(int Vs_state,int total_Vs_states,int nAch_state,int total_nAch_states,double input,double *ra)
{
  /* first attempting simple (sloppy) conductance based i&f 
     ra holds the output probabilities P(current_state --> state_i) 
     it is assumed that ra is already initialized with size total_Vs_states*total_nAch_states */
  int verbose=0;
  int ns=0,ns1=0,ns2=0;
  double DT=0,V=0,gE=0,VS=0,gT=0,Vbar=0,gEbar=0,ramean=0;
  double V_sigma=0,gE_sigma=0;
  double V_min=0,V_max=1,gE_min=0,gE_max=0.125;
  if (verbose){ printf(" %% [entering snx_statematrix_0] Vs_state %d/%d, nAch_state %d/%d, input %f\n",Vs_state,total_Vs_states,nAch_state,total_nAch_states,input);}
  DT = GLOBAL_DT;
  V_sigma = pow(2,-3)*V_max/(double)total_Vs_states; gE_sigma = pow(2,-3)*gE_max/(double)total_nAch_states;
  V = V_min + V_max*(double)Vs_state/(double)maximum(1,total_Vs_states-1);
  gE = gE_min + gE_max*(double)nAch_state/(double)maximum(1,total_nAch_states-1);
  gEbar = gE*exp(-DT/TAU_[VARNAME_REGISTRY_snx_nAch]);
  gEbar += input*((1-exp(-DT/TAU_[VARNAME_REGISTRY_snx_nAch]))*TAU_[VARNAME_REGISTRY_snx_nAch]/DT);
  gT = gEbar + 1.0/(double)TAU_[VARNAME_REGISTRY_snx_Vs];
  VS = gEbar*VOLTAGE_[VARNAME_REGISTRY_snx_nAch] + 1.0/(double)TAU_[VARNAME_REGISTRY_snx_Vs]*VOLTAGE_[VARNAME_REGISTRY_snx_Vs]; VS/=gT;
  Vbar = (V-VS)*exp(-DT*gT)+VS;
  Vbar = maximum(V_min,minimum(V_max,Vbar));
  gEbar = maximum(gE_min,minimum(gE_max,gEbar));
  if (verbose){ printf(" %% V %f Vbar %f VS %f gE %f gEbar %f gT %f \n",V,Vbar,VS,gE,gEbar,gT);}
  for (ns1=0;ns1<total_Vs_states;ns1++){ for (ns2=0;ns2<total_nAch_states;ns2++){
    ns = ns1+ns2*total_Vs_states;
    V = V_min + V_max*(double)ns1/(double)maximum(1,total_Vs_states-1);
    gE = gE_min + gE_max*(double)ns2/(double)maximum(1,total_nAch_states-1);
    ra[ns] = (Vs_state==maximum(1,total_Vs_states-1) ? ns1==0 : exp(-fabs(V-Vbar)/(2*V_sigma)));
    ra[ns] *= exp(-fabs(gE-gEbar)/(2*gE_sigma));}}
  stats("double",ra,total_Vs_states*total_nAch_states,NULL,NULL,&ramean,NULL); 
  if (verbose){ printf("V %d, gE %d\n",Vs_state,nAch_state); raprintf(ra,"double",total_Vs_states,total_nAch_states,"ra");}
  ratimesequals(ra,total_Vs_states*total_nAch_states,1.0/ramean/(double)(total_Vs_states*total_nAch_states));
  if (verbose){ printf(" %% [finishing snx_statematrix_0]\n");}
}

void snx_evolve(struct neuronarray *Nra,double t,double DT)
{
  /* first attempting simple all-to-all coupling */
  int verbose=0;
  int nr=0,nr1=0,nt=0,ns=0,ns1=0,ns2=0;
  struct neuron *n=NULL;
  int total_Vs_states=0,total_nAch_states=0,total_states=0;
  double coupling_strength=0, coupling_strength_2=0, total_coupling_strength=0,input=0,temp=0;
  double *ra=NULL,*ra_coupling=NULL;
  int version_flag = GLOBAL_SNX_VERSION;
  if (verbose){ printf(" %% [entering snx_evolve] t %f\n",t);}
  switch (version_flag){
  case 0:
    if (verbose){ printf(" %% version %d conductance based integrate-and-fire with all-to-all coupling\n",version_flag);}
    total_Vs_states=GLOBAL_SNX_NSTATES_[VARNAME_REGISTRY_snx_Vs]; 
    total_nAch_states=GLOBAL_SNX_NSTATES_[VARNAME_REGISTRY_snx_nAch]; 
    total_states=total_Vs_states*total_nAch_states;
    coupling_strength = CS__[TYPENAME_REGISTRY_snx_PN+TYPENAME_REGISTRY_snx_PN*GLOBAL_NTYPES+GLOBAL_INDEXING_REFILE_sra[VARNAME_REGISTRY_snx_nAch]*GLOBAL_NTYPES*GLOBAL_NTYPES];
    input = CS_ORN_[TYPENAME_REGISTRY_snx_PN];
    for (nt=0;nt<Nra->ntypes;nt++){ for (nr=0;nr<Nra->lengthra[nt];nr++){
      n=nget(Nra,nt,nr); 
      *(n->vpra[VARNAME_REGISTRY_snx_spike_flag]) = ((int)*(n->vpra[VARNAME_REGISTRY_snx_Vs])==maximum(1,total_Vs_states-1));
      if ((int)*(n->vpra[VARNAME_REGISTRY_snx_spike_flag])){
	n->spiketime_guess = t+/* rand01 */0.5*DT; n->spiketime = n->spiketime_guess; n->spiketime_guess_flag = 0; 
	n->spikelast = n->spiketime; n->spikenext = n->spikelast;}
      total_coupling_strength += coupling_strength*((int)*(n->vpra[VARNAME_REGISTRY_snx_spike_flag]));}}
    ra = (double *) tcalloc(total_states,sizeof(double));
    for (nt=0;nt<Nra->ntypes;nt++){ for (nr=0;nr<Nra->lengthra[nt];nr++){
      n=nget(Nra,nt,nr); 
      ns1 = (int)*(n->vpra[VARNAME_REGISTRY_snx_Vs]); ns1 = periodize(ns1,0,total_Vs_states);
      ns2 = (int)*(n->vpra[VARNAME_REGISTRY_snx_nAch]); ns2 = periodize(ns2,0,total_nAch_states);
      snx_statematrix_0(ns1,total_Vs_states,ns2,total_nAch_states,input+total_coupling_strength,ra);
      if (verbose){ printf(" %% neuron %d,%d ",nt,nr); raprintf(ra,"double",1,total_states,"ra: ");}
      temp=R01GET(&(n->spikeinput_rseed));
      if (verbose){ printf(" %% temp %f",temp);}
      ns=0;
      while (temp>ra[ns] && ns<(total_states-1)){ temp -= ra[ns]; ns += 1;}
      if (verbose){ printf(" ns %d",ns);}
      ns1 = ns%total_Vs_states; ns2 = ns/total_Vs_states;
      if (verbose){ printf(" ns1 %d ns2 %d\n",ns1,ns2);}
      *(n->vpra[VARNAME_REGISTRY_snx_Vs]) = ns1;
      *(n->vpra[VARNAME_REGISTRY_snx_nAch]) = ns2;}}
    tfree(ra);ra=NULL;
    break;
  case 1:
    if (verbose){ printf(" %% version %d current based integrate-and-fire with skeleton ring architecture\n",version_flag);}
    total_Vs_states=GLOBAL_SNX_NSTATES_[VARNAME_REGISTRY_snx_Vs]; 
    total_states=total_Vs_states;
    coupling_strength = CS__[TYPENAME_REGISTRY_snx_PN+TYPENAME_REGISTRY_snx_PN*GLOBAL_NTYPES+GLOBAL_INDEXING_REFILE_sra[VARNAME_REGISTRY_snx_nAch]*GLOBAL_NTYPES*GLOBAL_NTYPES];
    coupling_strength_2 = pow(coupling_strength,2);
    input = CS_ORN_[TYPENAME_REGISTRY_snx_PN];
    nt=0;
    if (verbose){ printf(" %% if SNX_BOTHER==0, then the coupling is a ring: nr is connected to periodize(nr+1,0,Nra->lengthra[0]) \n");}
    ra_coupling = (double *) tcalloc(Nra->lengthra[nt],sizeof(double));
    for (nr=0;nr<Nra->lengthra[nt];nr++){
      n=nget(Nra,nt,nr); 
      *(n->vpra[VARNAME_REGISTRY_snx_spike_flag]) = ((int)*(n->vpra[VARNAME_REGISTRY_snx_Vs])==maximum(1,total_Vs_states-1));
      if ((int)*(n->vpra[VARNAME_REGISTRY_snx_spike_flag])){
	if (verbose){ printf(" neuron %d fired\n",nr);}
	n->spiketime_guess = t+/* rand01 */0.5*DT; n->spiketime = n->spiketime_guess; n->spiketime_guess_flag = 0; 
	n->spikelast = n->spiketime; n->spikenext = n->spikelast;}
      if (SNXDATA_BOTHER && (int)*(n->vpra[VARNAME_REGISTRY_snx_spike_flag])){
	for (nr1=0;nr1<GLOBAL_SNXDATA->length;nr1++){
	  ra_coupling[periodize(nr1,0,GLOBAL_SNXDATA->length)] += GLOBAL_SNXDATA->connectivity[nr+nr1*GLOBAL_SNXDATA->length];}}
      else if (!SNXDATA_BOTHER && (int)*(n->vpra[VARNAME_REGISTRY_snx_spike_flag])){
	ra_coupling[periodize(nr+1,0,Nra->lengthra[nt])] += coupling_strength*((int)*(n->vpra[VARNAME_REGISTRY_snx_spike_flag]));}}
    if (verbose){ raprintf(ra_coupling,"double",1,Nra->lengthra[nt],"coupling: ");}
    ra = (double *) tcalloc(total_states,sizeof(double));
    nt=0; for (nr=0;nr<Nra->lengthra[nt];nr++){
      n=nget(Nra,nt,nr); 
      ns1 = (int)*(n->vpra[VARNAME_REGISTRY_snx_Vs]); ns1 = periodize(ns1,0,total_Vs_states);
      snx_statematrix_1(ns1,total_Vs_states,(SNXDATA_BOTHER ? GLOBAL_SNXDATA->input[nr] : input)+ra_coupling[nr],ra,NULL);
      if (verbose){ printf(" %% neuron (%d,%d) input %f+%f ",nt,nr,input,ra_coupling[nr]); raprintf(ra,"double",1,total_states,"ra: ");}
      temp=R01GET(&(n->spikeinput_rseed));
      if (verbose){ printf(" %% temp %f",temp);}
      ns=0;
      while (temp>ra[ns] && ns<(total_states-1)){ temp -= ra[ns]; ns += 1;}
      if (verbose){ printf(" goes from ns1 %d --> ns %d",ns1,ns);}
      ns1 = ns;
      if (verbose){ printf(" ns1 %d\n",ns1);}
      *(n->vpra[VARNAME_REGISTRY_snx_Vs]) = ns1;}
    tfree(ra);ra=NULL;
    tfree(ra_coupling);ra_coupling=NULL;
    break;
  case 2: 
    if (verbose){ printf(" %% version %d conductance based integrate-and-fire with two-way ring coupling\n",version_flag);}
    total_Vs_states=GLOBAL_SNX_NSTATES_[VARNAME_REGISTRY_snx_Vs]; 
    total_nAch_states=GLOBAL_SNX_NSTATES_[VARNAME_REGISTRY_snx_nAch]; 
    total_states=total_Vs_states*total_nAch_states;
    coupling_strength = CS__[TYPENAME_REGISTRY_snx_PN+TYPENAME_REGISTRY_snx_PN*GLOBAL_NTYPES+GLOBAL_INDEXING_REFILE_sra[VARNAME_REGISTRY_snx_nAch]*GLOBAL_NTYPES*GLOBAL_NTYPES];
    input = CS_ORN_[TYPENAME_REGISTRY_snx_PN];
    nt=0;
    ra_coupling = (double *) tcalloc(Nra->lengthra[nt],sizeof(double));
    for (nt=0;nt<Nra->ntypes;nt++){ for (nr=0;nr<Nra->lengthra[nt];nr++){
      n=nget(Nra,nt,nr); 
      *(n->vpra[VARNAME_REGISTRY_snx_spike_flag]) = ((int)*(n->vpra[VARNAME_REGISTRY_snx_Vs])==maximum(1,total_Vs_states-1));
      if ((int)*(n->vpra[VARNAME_REGISTRY_snx_spike_flag])){
	n->spiketime_guess = t+/* rand01 */0.5*DT; n->spiketime = n->spiketime_guess; n->spiketime_guess_flag = 0; 
	n->spikelast = n->spiketime; n->spikenext = n->spikelast;}
      ra_coupling[periodize(nr+1,0,Nra->lengthra[nt])] += coupling_strength*((int)*(n->vpra[VARNAME_REGISTRY_snx_spike_flag]));
      ra_coupling[periodize(nr-1,0,Nra->lengthra[nt])] += coupling_strength*((int)*(n->vpra[VARNAME_REGISTRY_snx_spike_flag]));}}
    ra = (double *) tcalloc(total_states,sizeof(double));
    for (nt=0;nt<Nra->ntypes;nt++){ for (nr=0;nr<Nra->lengthra[nt];nr++){
      n=nget(Nra,nt,nr); 
      ns1 = (int)*(n->vpra[VARNAME_REGISTRY_snx_Vs]); ns1 = periodize(ns1,0,total_Vs_states);
      ns2 = (int)*(n->vpra[VARNAME_REGISTRY_snx_nAch]); ns2 = periodize(ns2,0,total_nAch_states);
      snx_statematrix_0(ns1,total_Vs_states,ns2,total_nAch_states,input+ra_coupling[nr],ra);
      if (verbose){ printf(" %% neuron %d,%d ",nt,nr); raprintf(ra,"double",1,total_states,"ra: ");}
      temp=R01GET(&(n->spikeinput_rseed));
      if (verbose){ printf(" %% temp %f",temp);}
      ns=0;
      while (temp>ra[ns] && ns<(total_states-1)){ temp -= ra[ns]; ns += 1;}
      if (verbose){ printf(" ns %d",ns);}
      ns1 = ns%total_Vs_states; ns2 = ns/total_Vs_states;
      if (verbose){ printf(" ns1 %d ns2 %d\n",ns1,ns2);}
      *(n->vpra[VARNAME_REGISTRY_snx_Vs]) = ns1;
      *(n->vpra[VARNAME_REGISTRY_snx_nAch]) = ns2;}}
    tfree(ra);ra=NULL;
    tfree(ra_coupling);ra_coupling=NULL;
    break;
  default: break;}
  if (verbose){ printf(" %% [finishing snx_evolve]\n");}
}

void specialized_inverse_if(double M00,double R0,double *V0)
{
  /* specialized 1 flop inverse (row operations) for the matrix
    X    0 
    0  M00 
    with RHS [R0]
  */
  *V0 = R0/M00;
}

void nllistevolve_cif_helper(struct llist *L,struct neuron *s,double t_spike,int *iteration)
{
  int verbose=0;
  int nr=0;
  char lvlchar[512];
  double *sra = NULL;
  struct litem *l=NULL;
  struct neuron *n=NULL;
  double Vt=VOLTAGE_THRESHOLD_S,Vr=VOLTAGE_[VARNAME_REGISTRY_s_LEAK_S];
  if (verbose){ sprintf(lvlchar," %%"); for (nr=0;nr<*iteration;nr++){sprintf(lvlchar,"%s %%",lvlchar);}}
  if (verbose){ printf(" %% %s [entering nllistevolve_cif_helper] s(%d,%d) fired, *iteration %d\n",lvlchar,s->type,s->index,*iteration);}
  if (verbose){ printf(" %% spiketime_llist:\n"); l=L->first; while(l!=NULL){ n=(struct neuron *)l->item; printf("n (%d,%d), spikeinput_time %f\n",n->type,n->index,n->spikeinput_time); l=l->child;}}
  sra = (double *) tcalloc(GLOBAL_INDEXING_sra_LENGTH,sizeof(double));
  l=L->first;
  while (l!=NULL){
    n = (struct neuron *)l->item;
    if (verbose){ printf(" %% %s currently considering n(%d,%d) Vs %f\n",lvlchar,n->type,n->index,n->vra[VARNAME_REGISTRY_Vs]);}
    slink(s,n,sra,NULL); 
    if (n->spikenext<=t_spike){ if (verbose){ printf(" %% %s planning on adding %f\n",lvlchar,sra[0]);} sra_add(n,sra);}
    else /* if (n->spikenext>t_spike) */{ if (verbose){ printf(" %% %s not adding %f, since n is refracted\n",lvlchar,sra[0]);}}
    if (n->vra[VARNAME_REGISTRY_Vs]>Vt){
      n->vra[VARNAME_REGISTRY_Vs] -= (Vt-Vr);
      n->vra[VARNAME_REGISTRY_spike_flag] += 1;
      n->spiketime_guess = t_spike;
      n->spiketime = t_spike;
      n->spiketime_guess_flag = 0;
      n->spikenext = n->spiketime + TAU_REF;
      *iteration+=1; nllistevolve_cif_helper(L,n,n->spiketime,iteration);}
    l=l->child;}
  if (verbose){ printf(" %% %s [exiting nllistevolve_cif_helper]\n",lvlchar);}
  tfree(sra);sra=NULL;
}

void nllistevolve_cif_all2all_helper(struct llist *L,double t_spike,int *nspikes)
{
  /* assumes that *L is sorted with respect to Vs, and that only sra[0] matters */
  int verbose=0;
  double *sra = NULL,*sra_sum=NULL;
  struct litem *l=NULL;
  struct neuron *n=NULL;
  int continue_flag=0;
  double Vs=0,Vt=VOLTAGE_THRESHOLD_S,Vr=VOLTAGE_[VARNAME_REGISTRY_s_LEAK_S];
  if (verbose){ printf(" %% [entering nllistevolve_cif_all2all_helper] time %f, *nspikes %d\n",t_spike,*nspikes);}
  if (verbose){ printf(" %% Vs_llist:\n"); l=L->first; while(l!=NULL){ n=(struct neuron *)l->item; printf("n (%d,%d), Vs %f\n",n->type,n->index,n->vra[VARNAME_REGISTRY_Vs]); l=l->child;}}
  sra = (double *) tcalloc(GLOBAL_INDEXING_sra_LENGTH,sizeof(double));
  if (L->length>1){ slink((struct neuron *)L->first->item,(struct neuron *)L->first->child->item,sra,NULL);}
  sra_sum = (double *) tcalloc(GLOBAL_INDEXING_sra_LENGTH,sizeof(double)); sra_sum[0]=sra[0];
  if (verbose){ printf(" %% jump in Vs %f\n",sra[0]);}
  l=L->last;continue_flag=1;
  while (l!=NULL && continue_flag){
    n = (struct neuron *)l->item;
    Vs = n->vra[VARNAME_REGISTRY_Vs];
    if (verbose){ printf(" %% currently considering n(%d,%d) Vs %f\n",n->type,n->index,n->vra[VARNAME_REGISTRY_Vs]);}
    if (n->spikenext<=t_spike && Vs+sra_sum[0]>=Vt){
      if (verbose){ printf(" %% n fires...");}
      n->vra[VARNAME_REGISTRY_Vs] -= (Vt-Vr);
      n->vra[VARNAME_REGISTRY_spike_flag] += 1;
      n->spiketime_guess = t_spike;
      n->spiketime = t_spike;
      n->spiketime_guess_flag = 0;
      n->spikenext = n->spiketime + TAU_REF;
      if (verbose){ printf(" increasing jump in Vs from %f to %f\n",sra_sum[0],sra_sum[0]+sra[0]);}
      sra_sum[0] += sra[0];
      *nspikes+=1;
      continue_flag=1;}
    else if (n->spikenext>t_spike && Vs+sra_sum[0]>=Vt){
      if (verbose){ printf(" %% n would fire, but is refracted, skipping...\n");}
      continue_flag=1;}
    else if (Vs+sra_sum[0]<Vt){ if (verbose){ printf(" %% n won't fire, stopping\n");} continue_flag=0;}
    if (continue_flag){ l=l->parent;}}
  if (continue_flag==1){ if (verbose){ printf(" %% every neuron fired!\n");}}
  else /* if (continue_flag==0) */{
    while (l!=NULL){
      n = (struct neuron *)l->item;
      Vs = n->vra[VARNAME_REGISTRY_Vs];
      if (verbose){ printf(" %% currently considering n(%d,%d) Vs %f\n",n->type,n->index,n->vra[VARNAME_REGISTRY_Vs]);}
      if (n->spikenext<=t_spike && Vs+sra_sum[0]>=Vt){ printf(" %% error! this neuron shouldn't be firing\n");}
      else if (n->spikenext>t_spike && Vs+sra_sum[0]>=Vt){
	if (verbose){ printf(" %% n would fire, but is refracted, skipping...\n");}}
      else if (Vs+sra_sum[0]<Vt){ if (verbose){ printf(" %% as expected, n doesn't fire\n");} sra_add(n,sra_sum);}
      l=l->parent;}}
  if (verbose){ printf(" %% [exiting nllistevolve_cif_all2all_helper]\n");}
  tfree(sra);sra=NULL; tfree(sra_sum);sra_sum=NULL;
}

void nllistevolve_cif(struct llist *L,struct llist *L2,double t_start,double dt_start)
{
  /* current based integrate-and-fire */
  int verbose=0;
  double t_current=t_start;
  double dt_current=0;
  struct litem *l=NULL,*l2=NULL;
  struct neuron *n=NULL,*n2=NULL;
  double *sra=NULL;
  double Vs=0,Vt=VOLTAGE_THRESHOLD_S,Vr=VOLTAGE_[VARNAME_REGISTRY_s_LEAK_S];
  int input_flag=0,nr=0,nspikes=0;
  if (verbose){ printf(" %% [entering nllistevolve_cif], t_start %f, dt_start %f\n",t_start,dt_start);}
  sra = (double *) tcalloc(GLOBAL_INDEXING_sra_LENGTH,sizeof(double));
  llistsort(L->first,L->last,L->length,&spikeinput_time_compare);
  if (verbose>2){ printf(" %% spiketime_llist:\n"); l=L->first; while(l!=NULL){ n=(struct neuron *)l->item; printf("n (%d,%d), spikeinput_time %f\n",n->type,n->index,n->spikeinput_time); l=l->child;}}
  l=L->first; nr=0;
  while (l!=NULL){
    n = (struct neuron *)l->item;
    if (verbose){ printf(" %% at neuron(%d,%d), Vs %f, spikeinput_time %f, %d out of %d\n",n->type,n->index,n->vra[VARNAME_REGISTRY_Vs],n->spikeinput_time,nr,L->length);}
    input_flag=0;
    dt_current = maximum(0,t_start+dt_start-t_current);
    if (dt_current>n->spikeinput_time-t_current){ 
      dt_current = maximum(0,n->spikeinput_time-t_current);
      if (dt_current>0){ input_flag=1;}}
    if (verbose){ printf(" %% setting t_current %f, dt_current %f, input_flag=%d\n",t_current,dt_current,input_flag);}
    l2 = L->first;
    while (l2!=NULL){
      n2 = (struct neuron *)l2->item;
      n2->vra[VARNAME_REGISTRY_spike_flag] *= exp(-dt_current/0.25);
      Vs = n2->vra[VARNAME_REGISTRY_Vs];
      n2->vra[VARNAME_REGISTRY_Vs] = Vr + (Vs-Vr)*exp(-dt_current/TAU_[VARNAME_REGISTRY_Vs]) + (1-exp(-dt_current/TAU_[VARNAME_REGISTRY_Vs]))*TAU_[VARNAME_REGISTRY_Vs]*CURRENT_INJECTION_S;
      n2->vra[VARNAME_REGISTRY_Vd] = n2->vra[VARNAME_REGISTRY_Vs];
      l2 = l2->child;}
    if (input_flag){
      ilink(n,sra); 
      if (n->spikenext<=t_current+dt_current){ if (verbose){ printf(" %% planning on adding %f\n",sra[0]);} sra_add(n,sra);}
      else /* if (n->spikenext>t_current+dt_current) */{ if (verbose){ printf(" %% not adding %f, n refracted\n",sra[0]);}}}
    if (n->vra[VARNAME_REGISTRY_Vs]>Vt){
      if (verbose){ printf(" %% n (%d,%d) fired at time %f\n",n->type,n->index,t_current+dt_current);}
      n->vra[VARNAME_REGISTRY_Vs] -= (Vt-Vr);
      n->vra[VARNAME_REGISTRY_spike_flag] += 1;
      n->spiketime_guess = t_current+dt_current;
      n->spiketime = t_current+dt_current;
      n->spiketime_guess_flag = 0;
      n->spikenext = n->spiketime + TAU_REF;
/*       nspikes=1; llistsort(L2->first,L2->last,L2->length,&vra_Vs_compare);  */
/*       nllistevolve_cif_all2all_helper(L2,n->spiketime,&nspikes); */
      nspikes=1; nllistevolve_cif_helper(L,n,n->spiketime,&nspikes);
      if (nspikes==L->length && ISI_BOTHER){ isiupdate(GLOBAL_ISI,t_current+dt_current);}}
    t_current += dt_current;
    l=l->child; nr+=1;}
  tfree(sra);sra=NULL;
}

void vraevolve_if(struct neuron *n,double t_start,double dt_start)
{
  /* implicit euler */
  int verbose=0;//(n->index==0&&n->type==0);
  double temp1=0,temp2=0,temp3=0,temp7=0;
  double temp_s_A=0;
  double temp_s_N=0;
  double temp_s_NMDA=0;
  double temp_s_G=0;
  double temp_s_ORN=0;
  double temp_Vs=0;
  double temp_s_A_start=0;
  double temp_s_N_start=0;
  double temp_s_NMDA_start=0;
  double temp_s_G_start=0;
  double temp_s_ORN_start=0;
  double temp_I_syn=0;
  double temp_I_syn_Vs=0;
  double temp_I_leak_S=0;
  double temp_I_leak_S_Vs=0;
  double temp_rhs_Vs=0;
  double temp_rhs_Vs_Vs=0;
  double temp_Vs_start=0;
  double temp_Vs_end=0;
  double dt=dt_start,t=t_start;
  int microstep=0,error_flag=0;
  int iteration=0,iteration_max=4;
  double update_Vs=0;
  double norm_new=0,norm_old=0,norm_threshold=0.0000001;
  int nv=0;
  int non_refracted_flag=1;
  if (verbose>1){ printf(" %% [entering vraevolve_if]: \n"); raprintf(n->vra,"double",1,n->nvars," start ");}
  if (n->spikenext>=t_start+dt_start){
    if (verbose){ printf(" %% fully refracted, evolving %0.16f --> %0.16f\n",t_start,t_start+dt_start);}
    temp_s_A_start = n->vra[VARNAME_REGISTRY_s_A];
    temp_s_N_start = n->vra[VARNAME_REGISTRY_s_N];
    temp_s_NMDA_start = n->vra[VARNAME_REGISTRY_s_NMDA];
    temp_s_G_start = n->vra[VARNAME_REGISTRY_s_G];
    temp_s_ORN_start = n->vra[VARNAME_REGISTRY_s_ORN];
    temp_Vs_start=n->vra[VARNAME_REGISTRY_Vs];
    temp_s_A = temp_s_A_start * exp(-dt_start/TAU_[VARNAME_REGISTRY_s_A]);
    temp_s_N = temp_s_N_start * exp(-dt_start/TAU_[VARNAME_REGISTRY_s_N]);
    temp_s_NMDA = exp(-dt_start/TAU_[VARNAME_REGISTRY_s_NMDA])*(temp_s_NMDA_start + temp_s_N_start*(-1 + exp(dt_start/(1.0/TAU_[VARNAME_REGISTRY_s_NMDA]-1.0/TAU_[VARNAME_REGISTRY_s_N])))/(1.0/TAU_[VARNAME_REGISTRY_s_NMDA]-1.0/TAU_[VARNAME_REGISTRY_s_N]));
    temp_s_G = temp_s_G_start * exp(-dt_start/TAU_[VARNAME_REGISTRY_s_G]);
    temp_s_ORN = temp_s_ORN_start * exp(-dt_start/TAU_[VARNAME_REGISTRY_s_ORN]);
    temp_Vs = VOLTAGE_[VARNAME_REGISTRY_s_LEAK_S];
    temp_s_A_start = temp_s_A;
    temp_s_N_start = temp_s_N;
    temp_s_NMDA_start = temp_s_NMDA;
    temp_s_G_start = temp_s_G;
    temp_s_ORN_start = temp_s_ORN;
    temp_Vs_start = temp_Vs;
    t_start=t_start+dt_start;
    dt_start=0;
    temp_Vs_end = temp_Vs_start; temp_Vs_start = n->vra[VARNAME_REGISTRY_Vs];
    n->vra[VARNAME_REGISTRY_spike_flag] *= exp(-dt_start/(TAU_REF/16));
    n->vra[VARNAME_REGISTRY_s_A] = temp_s_A_start;
    n->vra[VARNAME_REGISTRY_s_N] = temp_s_N_start;
    n->vra[VARNAME_REGISTRY_s_NMDA] = temp_s_NMDA_start;
    n->vra[VARNAME_REGISTRY_s_G] = temp_s_G_start;
    n->vra[VARNAME_REGISTRY_s_ORN] = temp_s_ORN_start;
    n->vra[VARNAME_REGISTRY_Vs] = temp_Vs_end;
    non_refracted_flag=0;}  
  else if (n->spikenext>=t_start && n->spikenext < t_start+dt_start){
    if (verbose){ printf(" %% partially refracted, evolving %0.16f --> %0.16f\n",t_start,n->spikenext);}
    temp_s_A_start = n->vra[VARNAME_REGISTRY_s_A];
    temp_s_N_start = n->vra[VARNAME_REGISTRY_s_N];
    temp_s_NMDA_start = n->vra[VARNAME_REGISTRY_s_NMDA];
    temp_s_G_start = n->vra[VARNAME_REGISTRY_s_G];
    temp_s_ORN_start = n->vra[VARNAME_REGISTRY_s_ORN];
    temp_Vs_start=n->vra[VARNAME_REGISTRY_Vs];
    temp_s_A = temp_s_A_start * exp(-(n->spikenext-t_start)/TAU_[VARNAME_REGISTRY_s_A]);
    temp_s_N = temp_s_N_start * exp(-(n->spikenext-t_start)/TAU_[VARNAME_REGISTRY_s_N]);
    temp_s_NMDA = exp(-(n->spikenext-t_start)/TAU_[VARNAME_REGISTRY_s_NMDA])*(temp_s_NMDA_start + temp_s_N_start*(-1 + exp((n->spikenext-t_start)/(1.0/TAU_[VARNAME_REGISTRY_s_NMDA]-1.0/TAU_[VARNAME_REGISTRY_s_N])))/(1.0/TAU_[VARNAME_REGISTRY_s_NMDA]-1.0/TAU_[VARNAME_REGISTRY_s_N]));
    temp_s_G = temp_s_G_start * exp(-(n->spikenext-t_start)/TAU_[VARNAME_REGISTRY_s_G]);
    temp_s_ORN = temp_s_ORN_start * exp(-(n->spikenext-t_start)/TAU_[VARNAME_REGISTRY_s_ORN]);
    temp_Vs = VOLTAGE_[VARNAME_REGISTRY_s_LEAK_S];
    temp_s_A_start = temp_s_A;
    temp_s_N_start = temp_s_N;
    temp_s_NMDA_start = temp_s_NMDA;
    temp_s_G_start = temp_s_G;
    temp_s_ORN_start = temp_s_ORN;
    temp_Vs_start = temp_Vs;
    dt_start = (t_start+dt_start)-n->spikenext;
    t_start = n->spikenext;
    temp_Vs_end = temp_Vs_start; temp_Vs_start = n->vra[VARNAME_REGISTRY_Vs];
    n->vra[VARNAME_REGISTRY_spike_flag] *= exp(-dt_start/(TAU_REF/16));
    n->vra[VARNAME_REGISTRY_s_A] = temp_s_A_start;
    n->vra[VARNAME_REGISTRY_s_N] = temp_s_N_start;
    n->vra[VARNAME_REGISTRY_s_NMDA] = temp_s_NMDA_start;
    n->vra[VARNAME_REGISTRY_s_G] = temp_s_G_start;
    n->vra[VARNAME_REGISTRY_s_ORN] = temp_s_ORN_start;
    n->vra[VARNAME_REGISTRY_Vs] = temp_Vs_end;
    non_refracted_flag=1;}
  if (non_refracted_flag){
    if (verbose>1){ printf(" %% non refracted, evolving %0.16f --> %0.16f\n",t_start,t_start+dt_start);}
    do{ /* outer loop, choose n->microstep */
      error_flag=0;dt=dt_start/n->microstep;t=t_start;
      if (verbose>1){ printf(" %% initiating pass -- n->microstep %d dt %f\n",n->microstep,dt);}
      temp_s_A_start = n->vra[VARNAME_REGISTRY_s_A];
      temp_s_N_start = n->vra[VARNAME_REGISTRY_s_N];
      temp_s_NMDA_start = n->vra[VARNAME_REGISTRY_s_NMDA];
      temp_s_G_start = n->vra[VARNAME_REGISTRY_s_G];
      temp_s_ORN_start = n->vra[VARNAME_REGISTRY_s_ORN];
      temp_Vs_start=n->vra[VARNAME_REGISTRY_Vs];
      microstep=0;
      do{ /* each microstep */
	if (verbose>1){ printf(" %% %% microstep %d\n",microstep);}
	norm_old=16.0;norm_new=16.0;
	temp_s_A = temp_s_A_start * exp(-dt/TAU_[VARNAME_REGISTRY_s_A]);
	temp_s_N = temp_s_N_start * exp(-dt/TAU_[VARNAME_REGISTRY_s_N]);
	temp_s_NMDA = exp(-dt/TAU_[VARNAME_REGISTRY_s_NMDA])*(temp_s_NMDA_start + temp_s_N_start*(-1 + exp(dt*(1.0/TAU_[VARNAME_REGISTRY_s_NMDA]-1.0/TAU_[VARNAME_REGISTRY_s_N])))/(1.0/TAU_[VARNAME_REGISTRY_s_NMDA]-1.0/TAU_[VARNAME_REGISTRY_s_N]));
	temp_s_G = temp_s_G_start * exp(-dt/TAU_[VARNAME_REGISTRY_s_G]);
	temp_s_ORN = temp_s_ORN_start * exp(-dt/TAU_[VARNAME_REGISTRY_s_ORN]);
	temp_Vs = temp_Vs_start;
	iteration=0;
	do{ /* inner loop, linearly implicit euler */
	  if (verbose>2){ printf(" %% %% %% iteration %d\n",iteration);}
	  if (verbose>3){ printf(" %% %% %% %% calculating I_syn: ");}
	  temp1 = exp(-0.062*(temp_Vs+VOLTAGE_[VARNAME_REGISTRY_s_LEAK_D]-60)); temp2 = temp1*-0.062; temp1+=1; temp3 = -temp2/pow(temp1,2);
	  temp_I_syn = temp_s_ORN*(temp_Vs - VOLTAGE_[VARNAME_REGISTRY_s_ORN]) + temp_s_A*(temp_Vs - VOLTAGE_[VARNAME_REGISTRY_s_A]) + temp_s_NMDA/temp1*(temp_Vs - VOLTAGE_[VARNAME_REGISTRY_s_NMDA]) + temp_s_G*(temp_Vs - VOLTAGE_[VARNAME_REGISTRY_s_G]);
	  temp_I_syn_Vs = temp_s_ORN + temp_s_A + temp_s_NMDA*(1.0/temp1 + temp3*(temp_Vs-VOLTAGE_[VARNAME_REGISTRY_s_NMDA])) + temp_s_G;
	  if (verbose>3){ printf(" I_syn %f I_syn_Vs %f\n",temp_I_syn,temp_I_syn_Vs);}
	  if (verbose>3){ printf(" %% %% %% %% calculating I_leak_S: ");}
	  temp_I_leak_S_Vs = CONDUCTANCE_[VARNAME_REGISTRY_Vs];
	  temp_I_leak_S = temp_I_leak_S_Vs*(temp_Vs - VOLTAGE_[VARNAME_REGISTRY_s_LEAK_S]);
	  if (verbose>3){ printf(" I_leak_S %f, I_leak_S_Vs %f\n",temp_I_leak_S,temp_I_leak_S_Vs);}
	  if (verbose>2){ printf(" %% %% %% %% calculating: implicit euler update\n");}
	  temp_rhs_Vs = -(temp_I_leak_S+temp_I_syn+CURRENT_INJECTION_S)/TAU_[VARNAME_REGISTRY_Vs];
	  temp_rhs_Vs_Vs = -(temp_I_leak_S_Vs+temp_I_syn_Vs)/TAU_[VARNAME_REGISTRY_Vs];
	  specialized_inverse_if(1-dt*temp_rhs_Vs_Vs,temp_Vs_start - temp_Vs + dt*temp_rhs_Vs,&update_Vs);
	  norm_old = pow(update_Vs,2);
	  temp_Vs += update_Vs;
	  temp7=norm_old;norm_old=norm_new;norm_new=temp7;
	  iteration += 1;
	  if (verbose>2){ printf(" %% %% %% norm_new %0.16f norm_old %0.16f\n",norm_new,norm_old);}}
	while (norm_new>norm_threshold && norm_new<norm_old && iteration<iteration_max);
	if (verbose>1){ printf(" %% %% iteration=%d,norm_new %0.16f\n",iteration,norm_new);}
	temp_s_A_start = temp_s_A;
	temp_s_N_start = temp_s_N;
	temp_s_NMDA_start = temp_s_NMDA;
	temp_s_G_start = temp_s_G;
	temp_s_ORN_start = temp_s_ORN;
	temp_Vs_start = temp_Vs;
	microstep += 1; t+=dt;
	if (norm_new>norm_old || iteration==iteration_max || norm_new>norm_threshold){
	  if (verbose>2){ printf(" %% %% convergence poor, setting error_flag=1\n");}
	  error_flag=1;}
	else{
	  if (verbose>2){ printf(" %% %% convergence good, setting error_flag=0\n");}
	  error_flag=0;}}
      while (microstep<n->microstep && !error_flag);
      if (error_flag){ n->microstep*=2;} else{ n->microstep=maximum(1,n->microstep/2);}
      if (verbose>2){ printf(" %% total dt %f error_flag=%d, n->microstep set to %d\n",t-t_start,error_flag,n->microstep);}}
    while(error_flag);
    temp_Vs_end = temp_Vs_start; temp_Vs_start = n->vra[VARNAME_REGISTRY_Vs];
    if (verbose>2){ printf(" %% determining spiking properties\n");}
    if (n->spikenext<t_start+dt_start){
      temp1=-1; 
      if (temp_Vs_end>=VOLTAGE_THRESHOLD_S && temp_Vs_start<VOLTAGE_THRESHOLD_S){
	temp1 = linerootfinder(temp_Vs_start,temp_Vs_end,VOLTAGE_THRESHOLD_S,dt_start);}
      if (temp1>=0){
	n->spiketime_guess = t_start+temp1;
	if (n->spiketime_guess<=n->spikenext){
	  n->spiketime_guess=-2*dt_start; n->spiketime_guess_flag=0;}
	else /* if (n->spiketime_guess>n->spikenext) */{
	  n->vra[VARNAME_REGISTRY_spike_flag] += 1;
	  n->spiketime_guess_flag=1;
	  temp_Vs_end=VOLTAGE_[VARNAME_REGISTRY_s_LEAK_S];}}}
    n->vra[VARNAME_REGISTRY_spike_flag] *= exp(-dt_start/(TAU_REF/8));
    n->vra[VARNAME_REGISTRY_s_A] = temp_s_A_start;
    n->vra[VARNAME_REGISTRY_s_N] = temp_s_N_start;
    n->vra[VARNAME_REGISTRY_s_NMDA] = temp_s_NMDA_start;
    n->vra[VARNAME_REGISTRY_s_G] = temp_s_G_start;
    n->vra[VARNAME_REGISTRY_s_ORN] = temp_s_ORN_start;
    n->vra[VARNAME_REGISTRY_Vs] = temp_Vs_end;
    n->vra[VARNAME_REGISTRY_Vd] = n->vra[VARNAME_REGISTRY_Vs];}
  if (verbose>1){ printf(" %% exiting vraevolve: \n"); raprintf(n->vra,"double",1,n->nvars," end ");}
  if (verbose>1){ printf(" %% vra:\n"); for (nv=0;nv<n->nvars;nv++){ printf("%s %f \t",GLOBAL_VARNAMES[nv],n->vra[nv]);} printf("\n");}
}

void jacobian_to_lyapunov_morrislecar(double M00,double M01,double M02,double M10,double M11,double M20,double M22,double *output)
{
  /* return maximum real part of spectrum of
    X    0    1    2
    0  M00  M01  M02
    1  M10  M11   .
    2  M20   .   M22
  */
  int verbose=0;
  integer three=3,nine=9;
  doublecomplex A[9];
  doublecomplex B[3];
  integer INFO=0;
  doublecomplex WORK[9];
  doublereal RWORK[6];
  char en='N';
  int nr=0,nc=0;
  A[0 + 0*three].r = M00; A[0 + 0*three].i = 0.0;
  A[1 + 0*three].r = M10; A[1 + 0*three].i = 0.0;
  A[2 + 0*three].r = M20; A[2 + 0*three].i = 0.0;
  A[0 + 1*three].r = M01; A[0 + 1*three].i = 0.0;
  A[1 + 1*three].r = M11; A[1 + 1*three].i = 0.0;
  A[2 + 1*three].r = 0.0; A[2 + 1*three].i = 0.0;
  A[0 + 2*three].r = M02; A[0 + 2*three].i = 0.0;
  A[1 + 2*three].r = 0.0; A[1 + 2*three].i = 0.0;
  A[2 + 2*three].r = M22; A[2 + 2*three].i = 0.0;
  if (verbose){ printf("A=:\n");for (nr=0;nr<three;nr++){ for (nc=0;nc<three;nc++){ printf("%0.3f+%0.3fi\t",A[nr+nc*three].r,A[nr+nc*three].i);} printf("\n");}}
  zgeev_(&en,&en,&three,A,&three,B,NULL,&three,NULL,&three,WORK,&nine,RWORK,&INFO);
  if (INFO!=0){ printf("warning! info=%d in jacobian_to_lyapunov_morrislecar\n",(int)INFO);}
  if (verbose){ printf("B=:\n");for (nr=0;nr<three;nr++){ printf("%0.3f+%0.3fi\n",B[nr].r,B[nr].i);}}
  if (output!=NULL){
    *output=B[0].r;
    for (nr=0;nr<three;nr++){ *output=maximum(*output,B[nr].r);}
    if (verbose){ printf("biggest %f\n",*output);}}
}

void specialized_inverse_morrislecar(double M00,double M01,double M02,double M03,double M10,double M11,double M20,double M22,double M30,double M33,double R0,double R1,double R2,double R3,double *V0,double *V1,double *V2,double *V3)
{
  /* specialized 19 flop inverse (row operations) for the matrix
    X    0    1    2    3
    0  M00  M01  M02  M03
    1  M10  M11   .    .
    2  M20   .   M22   . 
    2  M30   .    .   M33
    with RHS [R0 R1 R2 R3]
  */
  double temp=0;
  temp = M03/M33; M00 -= temp*M30; R0 -= temp*R3; /* M03=0; */
  temp = M02/M22; M00 -= temp*M20; R0 -= temp*R2; /* M02=0; */
  temp = M01/M11; M00 -= temp*M10; R0 -= temp*R1; /* M01=0; */
  /* now matrix should look like
    X    0    1    2    3
    0  M00   .    .    .
    1  M10  M11   .    .
    2  M20   .   M22   . 
    2  M30   .    .   M33
  */
  temp = M10/M00; R1 -= temp*R0; /* M10=0; */
  temp = M20/M00; R2 -= temp*R0; /* M20=0; */
  temp = M30/M00; R3 -= temp*R0; /* M30=0; */
  /* now matrix should be diagonal */
  *V0 = R0/M00;
  *V1 = R1/M11;
  *V2 = R2/M22;
  *V3 = R3/M33;
}

void vraevolve_morrislecar(struct neuron *n,double t_start,double dt_start)
{
  /* linearly implicit euler */
  /* uses m_K as activation variable for K-reversal */
  int verbose=0;
  double morris_lecar_V1=-1.2;
  double morris_lecar_V2=18;
  double morris_lecar_V3=12;
  double morris_lecar_V4=17.4;
  double morris_lecar_phi=0.23;
  double morris_lecar_Ca0=1;
  double morris_lecar_mu=0.025;
  double temp1=0,temp2=0,temp3=0,temp7=0,temp_lyapunov=0;
  double alpha1=0,alpha1p=0,beta1=0,beta1p=0;
  double temp_s_A=0;
  double temp_s_N=0;
  double temp_s_NMDA=0;
  double temp_s_G=0;
  double temp_s_ORN=0;
  double temp_Vs=0;
  double temp_m_K=0;
  double temp_Vd=0;
  double temp_Ca=0;
  double temp_m_Ca_slow=0;
  double temp_s_K=0;
  double temp_s_Ca=0;
  double temp_s_KCa=0;
  double temp_s_Ca_slow=0;
  double temp_I_syn=0;
  double temp_I_syn_Vs=0;
  double temp_I_leak_S=0;
  double temp_I_leak_S_Vs=0;
  double temp_I_K=0;
  double temp_I_K_Vs=0;
  double temp_I_K_m_K=0;
  double temp_rhs_m_K=0;
  double temp_rhs_m_K_Vs=0;
  double temp_rhs_m_K_m_K=0;
  double temp_I_Ca=0;
  double temp_I_Ca_Vs=0;
  double temp_I_KCa=0;
  double temp_I_KCa_Vs=0;
  double temp_I_KCa_Ca=0;
  double temp_rhs_Ca=0;
  double temp_rhs_Ca_Vs=0;
  double temp_rhs_Ca_Ca=0;
  double temp_I_Ca_slow=0;
  double temp_I_Ca_slow_Vs=0;
  double temp_I_Ca_slow_m_Ca_slow=0;
  double temp_rhs_m_Ca_slow=0;
  double temp_rhs_m_Ca_slow_Vs=0;
  double temp_rhs_m_Ca_slow_m_Ca_slow=0;
  double temp_rhs_Vs=0;
  double temp_rhs_Vs_Vs=0;
  double temp_rhs_Vs_m_K=0;
  double temp_rhs_Vs_Ca=0;
  double temp_rhs_Vs_m_Ca_slow=0;
  double temp_s_A_start=0;
  double temp_s_N_start=0;
  double temp_s_NMDA_start=0;
  double temp_s_G_start=0;
  double temp_s_ORN_start=0;
  double temp_Vs_start=0;
  double temp_m_K_start=0;
  double temp_Vd_start=0;
  double temp_Ca_start=0;
  double temp_m_Ca_slow_start=0;
  double temp_Vs_end=0;
  double temp_Vd_end=0;
  double dt=dt_start,t=t_start;
  int microstep=0,error_flag=0;
  int iteration=0,iteration_max=4;
  double update_Vs=0;
  double update_m_K=0;
  double update_Ca=0;
  double update_m_Ca_slow=0;
  double norm_new=0,norm_old=0,norm_threshold=0.0000001;
  int nv=0;
  if (verbose>1){ printf(" %% [entering vraevolve_morrisslecar]: \n"); raprintf(n->vra,"double",1,n->nvars," start ");}
  do{ /* outer loop, choose n->microstep */
    error_flag=0;dt=dt_start/n->microstep;t=t_start;
    if (verbose>1){ printf(" %% initiating pass -- n->microstep %d dt %f\n",n->microstep,dt);}
    temp_s_A_start = n->vra[VARNAME_REGISTRY_s_A];
    temp_s_N_start = n->vra[VARNAME_REGISTRY_s_N];
    temp_s_NMDA_start = n->vra[VARNAME_REGISTRY_s_NMDA];
    temp_s_G_start = n->vra[VARNAME_REGISTRY_s_G];
    temp_s_ORN_start = n->vra[VARNAME_REGISTRY_s_ORN];
    temp_Vs_start=n->vra[VARNAME_REGISTRY_Vs];
    temp_m_K_start=n->vra[VARNAME_REGISTRY_m_K];
    temp_Vd_start=n->vra[VARNAME_REGISTRY_Vd];
    temp_Ca_start=n->vra[VARNAME_REGISTRY_Ca];
    temp_m_Ca_slow_start=n->vra[VARNAME_REGISTRY_m_Ca_slow];
    microstep=0;
    do{ /* each microstep */
      if (verbose>1){ printf(" %% %% microstep %d\n",microstep);}
      norm_old=16.0;norm_new=16.0;
      temp_s_A = temp_s_A_start * exp(-dt/TAU_[VARNAME_REGISTRY_s_A]);
      temp_s_N = temp_s_N_start * exp(-dt/TAU_[VARNAME_REGISTRY_s_N]);
      temp_s_NMDA = exp(-dt/TAU_[VARNAME_REGISTRY_s_NMDA])*(temp_s_NMDA_start + temp_s_N_start*(-1 + exp(dt*(1.0/TAU_[VARNAME_REGISTRY_s_NMDA]-1.0/TAU_[VARNAME_REGISTRY_s_N])))/(1.0/TAU_[VARNAME_REGISTRY_s_NMDA]-1.0/TAU_[VARNAME_REGISTRY_s_N]));
      temp_s_G = temp_s_G_start * exp(-dt/TAU_[VARNAME_REGISTRY_s_G]);
      temp_s_ORN = temp_s_ORN_start * exp(-dt/TAU_[VARNAME_REGISTRY_s_ORN]);
      temp_m_K = temp_m_K_start;
      temp_Ca = temp_Ca_start;
      temp_m_Ca_slow = temp_m_Ca_slow_start;
      temp_Vs = temp_Vs_start;
      temp_Vd = temp_Vd_start;
      iteration=0;
      do{ /* inner loop, linearly implicit euler */
	if (verbose>2){ printf(" %% %% %% iteration %d\n",iteration);}
	if (verbose>3){ printf(" %% %% %% %% calculating I_syn: ");}
	temp1 = exp(-0.062*(temp_Vs+VOLTAGE_[VARNAME_REGISTRY_s_LEAK_D]-60)); temp2 = temp1*-0.062; temp1+=1; temp3 = -temp2/pow(temp1,2);
	temp_I_syn = temp_s_ORN*(temp_Vs - VOLTAGE_[VARNAME_REGISTRY_s_ORN]) + temp_s_A*(temp_Vs - VOLTAGE_[VARNAME_REGISTRY_s_A]) + temp_s_NMDA/temp1*(temp_Vs - VOLTAGE_[VARNAME_REGISTRY_s_NMDA]) + temp_s_G*(temp_Vs - VOLTAGE_[VARNAME_REGISTRY_s_G]);
	temp_I_syn_Vs = temp_s_ORN + temp_s_A + temp_s_NMDA*(1.0/temp1 + temp3*(temp_Vs-VOLTAGE_[VARNAME_REGISTRY_s_NMDA])) + temp_s_G;
	if (verbose>3){ printf(" I_syn %f\n",temp_I_syn);}
	if (verbose>3){ printf(" %% %% %% %% calculating I_leak_S and I_leak_D: ");}
	temp_I_leak_S_Vs = CONDUCTANCE_[VARNAME_REGISTRY_Vs];
	temp_I_leak_S = temp_I_leak_S_Vs*(temp_Vs - VOLTAGE_[VARNAME_REGISTRY_s_LEAK_S]);
	if (verbose>3){ printf(" I_leak_S %f, I_leak_S_Vs %f\n",temp_I_leak_S,temp_I_leak_S_Vs);}
	if (verbose>3){ printf(" %% %% %% %% calculating I_K: ");}
	temp1 = (temp_Vs - VOLTAGE_[VARNAME_REGISTRY_s_K]);
	temp_I_K_m_K = CONDUCTANCE_[VARNAME_REGISTRY_s_K]*temp1;
	temp_I_K_Vs = CONDUCTANCE_[VARNAME_REGISTRY_s_K]*temp_m_K;
	temp_I_K = CONDUCTANCE_[VARNAME_REGISTRY_s_K]*temp_m_K*temp1;
	temp_s_K = temp_I_K_Vs;
	if (verbose>3){ printf(" I_K %f, I_K_Vs %f I_K_m_K %f\n",temp_I_K,temp_I_K_Vs,temp_I_K_m_K);}
	if (verbose>3){ printf(" %% %% %% %% calcuating rhs_m_K: ");}
	alpha1 = 0.5*(1+tanh((temp_Vs-morris_lecar_V3)/morris_lecar_V4)); 
	alpha1p = 0.5/pow(cosh((temp_Vs-morris_lecar_V3)/morris_lecar_V4),2)/morris_lecar_V4;
	beta1 = cosh((temp_Vs-morris_lecar_V3)/2/morris_lecar_V4); 
	beta1p = sinh((temp_Vs-morris_lecar_V3)/2/morris_lecar_V4)/2/morris_lecar_V4;
	temp_rhs_m_K = morris_lecar_phi*(alpha1 - temp_m_K)*beta1;
	temp_rhs_m_K_Vs = morris_lecar_phi*(alpha1p*beta1+alpha1*beta1p - temp_m_K*beta1p);
	temp_rhs_m_K_m_K = morris_lecar_phi*beta1;
	if (verbose>3){ printf(" rhs_m_K %f rhs_m_K_Vs %f rhs_m_K_m_K %f\n",temp_rhs_m_K,temp_rhs_m_K_Vs,temp_rhs_m_K_m_K);}
	if (verbose>3){ printf(" %% %% %% %% calculating I_Ca: ");}
	alpha1 = 0.5*(1+tanh((temp_Vs-morris_lecar_V1)/morris_lecar_V2)); 
	alpha1p = 0.5/pow(cosh((temp_Vs-morris_lecar_V1)/morris_lecar_V2),2)/morris_lecar_V2;
	temp1 = (temp_Vs - VOLTAGE_[VARNAME_REGISTRY_s_Ca]);
	temp_I_Ca = CONDUCTANCE_[VARNAME_REGISTRY_s_Ca]*alpha1*temp1;
	temp_I_Ca_Vs = CONDUCTANCE_[VARNAME_REGISTRY_s_Ca]*(alpha1p*temp1 + alpha1);
	temp_s_Ca = temp_I_Ca/temp1;
	if (verbose>3){ printf(" I_Ca %f, I_Ca_Vs %f\n",temp_I_Ca,temp_I_Ca_Vs);}
	if (verbose>3){ printf(" %% %% %% %% calculating I_KCa: ");}
	temp1 = (temp_Vs - VOLTAGE_[VARNAME_REGISTRY_s_KCa]);
	temp2 = temp_Ca/(temp_Ca+morris_lecar_Ca0); temp3 = morris_lecar_Ca0/pow(temp_Ca+morris_lecar_Ca0,2); 
	temp_I_KCa = CONDUCTANCE_[VARNAME_REGISTRY_s_KCa]*temp2*temp1;
	temp_I_KCa_Vs = CONDUCTANCE_[VARNAME_REGISTRY_s_KCa]*temp2;
	temp_I_KCa_Ca = CONDUCTANCE_[VARNAME_REGISTRY_s_KCa]*temp3*temp1;
	temp_s_KCa = temp_I_KCa/temp1;
	if (verbose>3){ printf(" I_KCa %f, I_KCa_Vs %f, I_KCa_Ca %f\n",temp_I_KCa,temp_I_KCa_Vs,temp_I_KCa_Ca);}
	if (verbose>3){ printf(" %% %% %% %% calculating rhs_Ca: ");}
	temp_rhs_Ca = (-morris_lecar_mu*temp_I_Ca - temp_Ca)/TAU_[VARNAME_REGISTRY_Ca];
	temp_rhs_Ca_Ca = -1/TAU_[VARNAME_REGISTRY_Ca];
	temp_rhs_Ca_Vs = -morris_lecar_mu*temp_I_Ca_Vs/TAU_[VARNAME_REGISTRY_Ca];
	if (verbose>3){ printf(" rhs_Ca %f, rhs_Ca_Vs %f, rhs_Ca_Ca %f\n",temp_rhs_Ca,temp_rhs_Ca_Vs,temp_rhs_Ca_Ca);}
	if (verbose>3){ printf(" %% %% %% %% calculating I_Ca_slow: ");}
	temp1 = (temp_Vs - VOLTAGE_[VARNAME_REGISTRY_s_Ca_slow]);
	temp_I_Ca_slow = CONDUCTANCE_[VARNAME_REGISTRY_s_Ca_slow]*temp_m_Ca_slow*temp1;
	temp_I_Ca_slow_Vs = CONDUCTANCE_[VARNAME_REGISTRY_s_Ca_slow]*temp_m_Ca_slow;
	temp_I_Ca_slow_m_Ca_slow = CONDUCTANCE_[VARNAME_REGISTRY_s_Ca_slow]*temp1;
	temp_s_Ca_slow = temp_I_Ca_slow/temp1;
	if (verbose>3){ printf(" I_Ca_slow %f, I_Ca_slow_Vs %f, I_Ca_slow_m_Ca_slow %f\n",temp_I_Ca_slow,temp_I_Ca_slow_Vs,temp_I_Ca_slow_m_Ca_slow);}
	if (verbose>3){ printf(" %% %% %% %% calculating rhs_m_Ca_slow: ");}
	temp1 = 0.5*(1+tanh((temp_Vs-12)/24)); temp2 = 0.5/pow(cosh((temp_Vs-12)/24),2)/24;
	temp_rhs_m_Ca_slow = (temp1 - temp_m_Ca_slow)/TAU_[VARNAME_REGISTRY_m_Ca_slow];
	temp_rhs_m_Ca_slow_Vs = temp2/TAU_[VARNAME_REGISTRY_m_Ca_slow];
	temp_rhs_m_Ca_slow_m_Ca_slow = -1.0/TAU_[VARNAME_REGISTRY_m_Ca_slow];
	if (verbose>3){ printf(" rhs_m_Ca_slow %f m_Ca_slow_Vs %f m_Ca_slow_m_Ca_slow %f\n",temp_rhs_m_Ca_slow,temp_rhs_m_Ca_slow_Vs,temp_rhs_m_Ca_slow_m_Ca_slow);}
	if (verbose>2){ printf(" %% %% %% %% calculating: implicit euler update\n");}
	temp_rhs_Vs = -(temp_I_leak_S+temp_I_Ca+temp_I_K+temp_I_KCa+temp_I_Ca_slow+temp_I_syn+CURRENT_INJECTION_S)/TAU_[VARNAME_REGISTRY_Vs];
	temp_rhs_Vs_Vs = -(temp_I_leak_S_Vs+temp_I_Ca_Vs+temp_I_K_Vs+temp_I_KCa_Vs+temp_I_Ca_slow_Vs+temp_I_syn_Vs)/TAU_[VARNAME_REGISTRY_Vs];
	temp_rhs_Vs_m_K = -(temp_I_K_m_K)/TAU_[VARNAME_REGISTRY_Vs];
	temp_rhs_Vs_Ca = -(temp_I_KCa_Ca)/TAU_[VARNAME_REGISTRY_Vs];
	temp_rhs_Vs_m_Ca_slow = -(temp_I_Ca_slow_m_Ca_slow)/TAU_[VARNAME_REGISTRY_Vs];
	specialized_inverse_morrislecar(1-dt*temp_rhs_Vs_Vs,-dt*temp_rhs_Vs_m_K,-dt*temp_rhs_Vs_Ca,-dt*temp_rhs_Vs_m_Ca_slow,-dt*temp_rhs_m_K_Vs,1-dt*temp_rhs_m_K_m_K,-dt*temp_rhs_Ca_Vs,1-dt*temp_rhs_Ca_Ca,-dt*temp_rhs_m_Ca_slow_Vs,1-dt*temp_rhs_m_Ca_slow_m_Ca_slow,temp_Vs_start - temp_Vs + dt*temp_rhs_Vs,temp_m_K_start - temp_m_K + dt*temp_rhs_m_K,temp_Ca_start - temp_Ca + dt*temp_rhs_Ca,temp_m_Ca_slow_start - temp_m_Ca_slow + dt*temp_rhs_m_Ca_slow,&update_Vs,&update_m_K,&update_Ca,&update_m_Ca_slow);
	switch (LYAPUNOV_BOTHER){
	case 0: break;
	case 1: jacobian_to_lyapunov_morrislecar(temp_rhs_Vs_Vs,temp_rhs_Vs_m_K,temp_rhs_Vs_Ca,temp_rhs_m_K_Vs,temp_rhs_m_K_m_K,temp_rhs_Ca_Vs,temp_rhs_Ca_Ca,&temp_lyapunov); break;
	default: break;}
	norm_old = pow(update_Vs,2)+pow(update_m_K,2)+pow(update_Ca,2);
	temp_Vs += update_Vs;
	temp_m_K += update_m_K;
	temp_Ca += update_Ca;
	temp_m_Ca_slow += update_m_Ca_slow;
	temp7=norm_old;norm_old=norm_new;norm_new=temp7;
	iteration += 1;
	if (verbose>2){ printf(" %% %% %% norm_new %0.16f norm_old %0.16f\n",norm_new,norm_old);}}
      while (norm_new>norm_threshold && norm_new<norm_old && iteration<iteration_max);
      if (verbose>1){ printf(" %% %% iteration=%d,norm_new %0.16f\n",iteration,norm_new);}
      temp_s_A_start = temp_s_A;
      temp_s_N_start = temp_s_N;
      temp_s_NMDA_start = temp_s_NMDA;
      temp_s_G_start = temp_s_G;
      temp_s_ORN_start = temp_s_ORN;
      temp_m_K_start = maximum(0,minimum(1,temp_m_K));
      temp_Ca_start = maximum(0,temp_Ca);
      temp_m_Ca_slow_start = maximum(0,minimum(1,temp_m_Ca_slow));
      temp_Vs_start = temp_Vs;
      temp_Vd_start = temp_Vd;
      microstep += 1; t+=dt;
      if (norm_new>norm_old || iteration==iteration_max || norm_new>norm_threshold){
	if (verbose>2){ printf(" %% %% convergence poor, setting error_flag=1\n");}
	error_flag=1;}
      else{
	if (verbose>2){ printf(" %% %% convergence good, setting error_flag=0\n");}
	error_flag=0;}}
    while (microstep<n->microstep && !error_flag);
    if (error_flag){ n->microstep*=2;} else{ n->microstep=maximum(1,n->microstep/2);}
    if (verbose>2){ printf(" %% total dt %f error_flag=%d, n->microstep set to %d\n",t-t_start,error_flag,n->microstep);}}
  while(error_flag);
  temp_Vs_end = temp_Vs_start; temp_Vs_start = n->vra[VARNAME_REGISTRY_Vs];
  temp_Vd_end = temp_Vd_start; temp_Vd_start = n->vra[VARNAME_REGISTRY_Vd];
  if (verbose>2){ printf(" %% determining spiking properties\n");}
  if (n->spikenext<t_start+dt_start){
    temp1=-1;
    if (temp_Vs_end>=VOLTAGE_THRESHOLD_S && temp_Vs_start<VOLTAGE_THRESHOLD_S){
      temp1 = linerootfinder(temp_Vs_start,temp_Vs_end,VOLTAGE_THRESHOLD_S,dt_start);}
    if (temp1>=0){
      n->spiketime_guess = t_start+temp1;
      if (n->spiketime_guess<=n->spikenext){
	n->spiketime_guess=-2*dt_start; n->spiketime_guess_flag=0;}
      else /* if (n->spiketime_guess>n->spikenext) */{
	n->vra[VARNAME_REGISTRY_spike_flag] += 1;
	n->spiketime_guess_flag=1;}}}
  n->vra[VARNAME_REGISTRY_spike_flag] *= exp(-dt_start/(TAU_REF/16));
  n->vra[VARNAME_REGISTRY_s_A] = temp_s_A_start;
  n->vra[VARNAME_REGISTRY_s_N] = temp_s_N_start;
  n->vra[VARNAME_REGISTRY_s_NMDA] = temp_s_NMDA_start;
  n->vra[VARNAME_REGISTRY_s_G] = temp_s_G_start;
  n->vra[VARNAME_REGISTRY_s_ORN] = temp_s_ORN_start;
  n->vra[VARNAME_REGISTRY_m_K] = temp_m_K_start;
  n->vra[VARNAME_REGISTRY_Ca] = temp_Ca_start;
  n->vra[VARNAME_REGISTRY_m_Ca_slow] = temp_m_Ca_slow_start;
  n->vra[VARNAME_REGISTRY_Vs] = temp_Vs_end;
  n->vra[VARNAME_REGISTRY_Vd] = n->vra[VARNAME_REGISTRY_Vs];
  n->vra[VARNAME_REGISTRY_s_K] = temp_s_K;
  n->vra[VARNAME_REGISTRY_s_Ca] = temp_s_Ca;
  n->vra[VARNAME_REGISTRY_s_KCa] = temp_s_KCa;
  n->vra[VARNAME_REGISTRY_s_Ca_slow] = temp_s_Ca_slow;
  n->vra[VARNAME_REGISTRY_lyapunov] = temp_lyapunov;
  if (verbose>1){ printf(" %% exiting vraevolve: \n"); raprintf(n->vra,"double",1,n->nvars," end ");}
  if (verbose>0){ printf(" %% vra:\n"); for (nv=0;nv<n->nvars;nv++){ printf("%s %f \t",GLOBAL_VARNAMES[nv],n->vra[nv]);} printf("\n");}
}

void jacobian_to_lyapunov_adi(double M00,double M01,double M02,double M03,double M10,double M11,double M20,double M22,double M30,double M33,double M34,double M43,double M44,double *output)
{
  /* return maximum real part of spectrum of
    X    0    1    2    3    4
    0  M00  M01  M02  M03   .
    1  M10  M11   .    .    .
    2  M20   .   M22   .    .
    3  M30   .    .   M33  M34
    4   .    .    .   M43  M44
  */
  int verbose=0;
  integer five=5,twentyfive=25;
  doublecomplex A[25];
  doublecomplex B[5];
  integer INFO=0;
  doublecomplex WORK[25];
  doublereal RWORK[10];
  char en='N';
  int nr=0,nc=0;
  A[0 + 0*five].r = M00; A[0 + 0*five].i = 0.0;
  A[1 + 0*five].r = M10; A[1 + 0*five].i = 0.0;
  A[2 + 0*five].r = M20; A[2 + 0*five].i = 0.0;
  A[3 + 0*five].r = M30; A[3 + 0*five].i = 0.0;
  A[4 + 0*five].r = 0.0; A[4 + 0*five].i = 0.0;
  A[0 + 1*five].r = M01; A[0 + 1*five].i = 0.0;
  A[1 + 1*five].r = M11; A[1 + 1*five].i = 0.0;
  A[2 + 1*five].r = 0.0; A[2 + 1*five].i = 0.0;
  A[3 + 1*five].r = 0.0; A[3 + 1*five].i = 0.0;
  A[4 + 1*five].r = 0.0; A[4 + 1*five].i = 0.0;
  A[0 + 2*five].r = M02; A[0 + 2*five].i = 0.0;
  A[1 + 2*five].r = 0.0; A[1 + 2*five].i = 0.0;
  A[2 + 2*five].r = M22; A[2 + 2*five].i = 0.0;
  A[3 + 2*five].r = 0.0; A[3 + 2*five].i = 0.0;
  A[4 + 2*five].r = 0.0; A[4 + 2*five].i = 0.0;
  A[0 + 3*five].r = M03; A[0 + 3*five].i = 0.0;
  A[1 + 3*five].r = 0.0; A[1 + 3*five].i = 0.0;
  A[2 + 3*five].r = 0.0; A[2 + 3*five].i = 0.0;
  A[3 + 3*five].r = M33; A[3 + 3*five].i = 0.0;
  A[4 + 3*five].r = M43; A[4 + 3*five].i = 0.0;
  A[0 + 4*five].r = 0.0; A[0 + 4*five].i = 0.0;
  A[1 + 4*five].r = 0.0; A[1 + 4*five].i = 0.0;
  A[2 + 4*five].r = 0.0; A[2 + 4*five].i = 0.0;
  A[3 + 4*five].r = M34; A[3 + 4*five].i = 0.0;
  A[4 + 4*five].r = M44; A[4 + 4*five].i = 0.0;
  if (verbose){ printf("A=:\n");for (nr=0;nr<five;nr++){ for (nc=0;nc<five;nc++){ printf("%0.3f+%0.3fi\t",A[nr+nc*five].r,A[nr+nc*five].i);} printf("\n");}}
  zgeev_(&en,&en,&five,A,&five,B,NULL,&five,NULL,&five,WORK,&twentyfive,RWORK,&INFO);
  if (INFO!=0){ printf("warning! info=%d in jacobian_to_lyapunov_adi\n",(int)INFO);}
  if (verbose){ printf("B=:\n");for (nr=0;nr<five;nr++){ printf("%0.3f+%0.3fi\n",B[nr].r,B[nr].i);}}
  if (output!=NULL){
    *output=B[0].r;
    for (nr=0;nr<five;nr++){ *output=maximum(*output,B[nr].r);}
    if (verbose){ printf("biggest %f\n",*output);}}
}

void specialized_inverse_adi(double M00,double M01,double M02,double M03,double M10,double M11,double M20,double M22,double M30,double M33,double M34,double M43,double M44,double R0,double R1,double R2,double R3,double R4,double *V0,double *V1,double *V2,double *V3,double *V4)
{
  /* specialized 37 flop inverse (row operations) for the matrix
    X    0    1    2    3    4
    0  M00  M01  M02  M03   .
    1  M10  M11   .    .    .
    2  M20   .   M22   .    .
    3  M30   .    .   M33  M34
    4   .    .    .   M43  M44
    with RHS [R0 R1 R2 R3 R4]
  */
  double temp=0;
  temp = M34/M44; M33 -= temp*M43; R3 -= temp*R4; /* M34=0; */
  temp = M03/M33; M00 -= temp*M30; R0 -= temp*R3; /* M03=0; */
  temp = M02/M22; M00 -= temp*M20; R0 -= temp*R2; /* M02=0; */
  temp = M01/M11; M00 -= temp*M10; R0 -= temp*R1; /* M01=0; */
  /* now matrix should look like
    X    0    1    2    3    4
    0  M00   .    .    .    .
    1  M10  M11   .    .    .
    2  M20   .   M22   .    .
    3  M30   .    .   M33   .
    4   .    .    .   M43  M44
  */
  temp = M10/M00; R1 -= temp*R0; /* M10=0; */
  temp = M20/M00; R2 -= temp*R0; /* M20=0; */
  temp = M30/M00; R3 -= temp*R0; /* M30=0; */
  temp = M43/M33; R4 -= temp*R3; /* M43=0; */
  /* now matrix should be diagonal */
  *V0 = R0/M00;
  *V1 = R1/M11;
  *V2 = R2/M22;
  *V3 = R3/M33;
  *V4 = R4/M44;
}

void vraevolve_adi(struct neuron *n,double t_start,double dt_start)
{
  /* linearly implicit euler */
  int verbose=0;
  int Ca2KCa_power=0;
  double temp1=0,temp2=0,temp3=0,temp5=0,temp6=0,temp7=0,temp_lyapunov=0;
  double alpha1=0,alpha1p=0,beta1=0,beta1p=0;
  double temp_s_A=0;
  double temp_s_N=0;
  double temp_s_NMDA=0;
  double temp_s_G=0;
  double temp_s_ORN=0;
  double temp_Vs=0;
  double temp_h_Na=0;
  double temp_m_K=0;
  double temp_Vd=0;
  double temp_Ca=0;
  double temp_s_Na=0;
  double temp_s_K=0;
  double temp_s_Ca=0;
  double temp_s_KCa=0;
  double temp_I_syn=0;
  double temp_I_syn_Vd=0;
  double temp_I_leak_S=0;
  double temp_I_leak_S_Vs=0;
  double temp_I_leak_D=0;
  double temp_I_leak_D_Vd=0;
  double temp_I_Na=0;
  double temp_I_Na_Vs=0;
  double temp_I_Na_h_Na=0;
  double temp_rhs_h_Na=0;
  double temp_rhs_h_Na_Vs=0;
  double temp_rhs_h_Na_h_Na=0;
  double temp_I_K=0;
  double temp_I_K_Vs=0;
  double temp_I_K_m_K=0;
  double temp_rhs_m_K=0;
  double temp_rhs_m_K_Vs=0;
  double temp_rhs_m_K_m_K=0;
  double temp_I_SD=0;
  double temp_I_SD_Vs=0;
  double temp_I_SD_Vd=0;
  double temp_I_DS=0;
  double temp_I_DS_Vs=0;
  double temp_I_DS_Vd=0;
  double temp_I_Ca=0;
  double temp_I_Ca_Vd=0;
  double temp_I_KCa=0;
  double temp_I_KCa_Vd=0;
  double temp_I_KCa_Ca=0;
  double temp_rhs_Ca=0;
  double temp_rhs_Ca_Vd=0;
  double temp_rhs_Ca_Ca=0;
  double temp_rhs_Vs=0;
  double temp_rhs_Vs_Vs=0;
  double temp_rhs_Vs_h_Na=0;
  double temp_rhs_Vs_m_K=0;
  double temp_rhs_Vs_Vd=0;
  double temp_rhs_Vd=0;
  double temp_rhs_Vd_Vs=0;
  double temp_rhs_Vd_Vd=0;
  double temp_rhs_Vd_Ca=0;
  double temp_s_A_start=0;
  double temp_s_N_start=0;
  double temp_s_NMDA_start=0;
  double temp_s_G_start=0;
  double temp_s_ORN_start=0;
  double temp_Vs_start=0;
  double temp_h_Na_start=0;
  double temp_m_K_start=0;
  double temp_Vd_start=0;
  double temp_Ca_start=0;
  double temp_Vs_end=0;
  double temp_Vd_end=0;
  double dt=dt_start,t=t_start;
  int microstep=0,error_flag=0;
  int iteration=0,iteration_max=4;
  double update_Vs=0;
  double update_h_Na=0;
  double update_m_K=0;
  double update_Vd=0;
  double update_Ca=0;
  double norm_new=0,norm_old=0,norm_threshold=0.0000001;
  int nv=0;
  if (verbose>1){ printf(" %% [entering vraevolve_adi]: \n"); raprintf(n->vra,"double",1,n->nvars," start ");}
  do{ /* outer loop, choose n->microstep */
    error_flag=0;dt=dt_start/n->microstep;t=t_start;
    if (verbose>1){ printf(" %% initiating pass -- n->microstep %d dt %f\n",n->microstep,dt);}
    temp_s_A_start = n->vra[VARNAME_REGISTRY_s_A];
    temp_s_N_start = n->vra[VARNAME_REGISTRY_s_N];
    temp_s_NMDA_start = n->vra[VARNAME_REGISTRY_s_NMDA];
    temp_s_G_start = n->vra[VARNAME_REGISTRY_s_G];
    temp_s_ORN_start = n->vra[VARNAME_REGISTRY_s_ORN];
    temp_Vs_start=n->vra[VARNAME_REGISTRY_Vs];
    temp_h_Na_start=n->vra[VARNAME_REGISTRY_h_Na];
    temp_m_K_start=n->vra[VARNAME_REGISTRY_m_K];
    temp_Vd_start=n->vra[VARNAME_REGISTRY_Vd];
    temp_Ca_start=n->vra[VARNAME_REGISTRY_Ca];
    microstep=0;
    do{ /* each microstep */
      if (verbose>1){ printf(" %% %% microstep %d\n",microstep);}
      norm_old=16.0;norm_new=16.0;
      temp_s_A = temp_s_A_start * exp(-dt/TAU_[VARNAME_REGISTRY_s_A]);
      temp_s_N = temp_s_N_start * exp(-dt/TAU_[VARNAME_REGISTRY_s_N]);
      temp_s_NMDA = exp(-dt/TAU_[VARNAME_REGISTRY_s_NMDA])*(temp_s_NMDA_start + temp_s_N_start*(-1 + exp(dt*(1.0/TAU_[VARNAME_REGISTRY_s_NMDA]-1.0/TAU_[VARNAME_REGISTRY_s_N])))/(1.0/TAU_[VARNAME_REGISTRY_s_NMDA]-1.0/TAU_[VARNAME_REGISTRY_s_N]));
      temp_s_G = temp_s_G_start * exp(-dt/TAU_[VARNAME_REGISTRY_s_G]);
      temp_s_ORN = temp_s_ORN_start * exp(-dt/TAU_[VARNAME_REGISTRY_s_ORN]);
      temp_h_Na = temp_h_Na_start;
      temp_m_K = temp_m_K_start;
      temp_Ca = temp_Ca_start;
      temp_Vs = temp_Vs_start;
      temp_Vd = temp_Vd_start;
      iteration=0;
      do{ /* inner loop, linearly implicit euler */
	if (verbose>2){ printf(" %% %% %% iteration %d\n",iteration);}
	if (verbose>3){ printf(" %% %% %% %% calculating I_syn: ");}
	temp1 = exp(-0.062*(temp_Vd+VOLTAGE_[VARNAME_REGISTRY_s_LEAK_D]-60)); temp2 = temp1*-0.062; temp1+=1; temp3 = -temp2/pow(temp1,2);
	temp_I_syn = temp_s_ORN*(temp_Vd - VOLTAGE_[VARNAME_REGISTRY_s_ORN]) + temp_s_A*(temp_Vd - VOLTAGE_[VARNAME_REGISTRY_s_A]) + temp_s_NMDA/temp1*(temp_Vd - VOLTAGE_[VARNAME_REGISTRY_s_NMDA]) + temp_s_G*(temp_Vd - VOLTAGE_[VARNAME_REGISTRY_s_G]);
	temp_I_syn_Vd = temp_s_ORN + temp_s_A + temp_s_NMDA*(1.0/temp1 + temp3*(temp_Vd-VOLTAGE_[VARNAME_REGISTRY_s_NMDA])) + temp_s_G;
	if (verbose>3){ printf(" I_syn %f\n",temp_I_syn);}
	if (verbose>3){ printf(" %% %% %% %% calculating I_leak_S and I_leak_D: ");}
	temp_I_leak_S_Vs = CONDUCTANCE_[VARNAME_REGISTRY_Vs];
	temp_I_leak_S = temp_I_leak_S_Vs*(temp_Vs - VOLTAGE_[VARNAME_REGISTRY_s_LEAK_S]);
	temp_I_leak_D_Vd = CONDUCTANCE_[VARNAME_REGISTRY_Vd];
	temp_I_leak_D = temp_I_leak_D_Vd*(temp_Vd - VOLTAGE_[VARNAME_REGISTRY_s_LEAK_D]);
	if (verbose>3){ printf(" I_leak_S %f, I_leak_S_Vs %f I_leak_D %f I_leak_D_Vd %f\n",temp_I_leak_S,temp_I_leak_S_Vs,temp_I_leak_D,temp_I_leak_D_Vd);}
	if (verbose>3){ printf(" %% %% %% %% calculating I_Na: ");}
	temp1 = (temp_Vs+33)/10; temp2 = -exp(-temp1); temp3 = -temp2/10; temp2+=1;
	alpha1 = temp1/temp2; alpha1p = (0.1*temp2-temp1*temp3)/pow(temp2,2);
	temp1 = (temp_Vs+53.7)/12; temp2 = 4*exp(-temp1); temp3 = -temp2/12;
	beta1 = temp2; beta1p = temp3;
	temp5 = alpha1/(alpha1+beta1); temp6 = (alpha1p*(alpha1+beta1) - alpha1*(alpha1p+beta1p))/pow(alpha1+beta1,2);
	temp1 = (temp_Vs - VOLTAGE_[VARNAME_REGISTRY_s_Na]);
	temp2 = pow(temp5,3)*temp1;
	temp3 = 3*pow(temp5,2)*temp6*temp1 + pow(temp5,3);
	temp_I_Na_h_Na = CONDUCTANCE_[VARNAME_REGISTRY_s_Na]*temp2;
	temp_I_Na = temp_I_Na_h_Na*temp_h_Na;
	temp_I_Na_Vs = CONDUCTANCE_[VARNAME_REGISTRY_s_Na]*temp_h_Na*temp3;
	temp_s_Na = CONDUCTANCE_[VARNAME_REGISTRY_s_Na]*pow(temp5,3)*temp_h_Na;
	if (verbose>3){ printf(" m_Na_inf %f, h_Na %f, I_Na %f, I_Na_Vs %f I_Na_h_Na %f\n",temp5,temp_h_Na,temp_I_Na,temp_I_Na_Vs,temp_I_Na_h_Na);}
	if (verbose>3){ printf(" %% %% %% %% calculating rhs_h_Na: ");}
	temp1 = (temp_Vs+50)/10; temp2 = 0.07*exp(-temp1); temp3 = -temp2/10;
	alpha1 = temp2; alpha1p = temp3;
	temp1 = (temp_Vs+20)/10; temp2 = exp(-temp1); temp3 = -temp2/10; temp2+=1;
	beta1 = 1/temp2; beta1p = -temp3/pow(temp2,2);
	temp_rhs_h_Na_h_Na = -4*(alpha1+beta1);
	temp_rhs_h_Na = 4*alpha1 + temp_h_Na*temp_rhs_h_Na_h_Na;
	temp_rhs_h_Na_Vs = 4*(alpha1p - temp_h_Na*(alpha1p + beta1p));
	if (verbose>3){ printf(" rhs_h_Na %f rhs_h_Na_Vs %f rhs_h_Na_h_Na %f\n",temp_rhs_h_Na,temp_rhs_h_Na_Vs,temp_rhs_h_Na_h_Na);}
	if (verbose>3){ printf(" %% %% %% %% calculating I_K: ");}
	temp1 = (temp_Vs - VOLTAGE_[VARNAME_REGISTRY_s_K]);
	temp_I_K_m_K = CONDUCTANCE_[VARNAME_REGISTRY_s_K]*4*pow(temp_m_K,3)*temp1;
	temp_I_K_Vs = CONDUCTANCE_[VARNAME_REGISTRY_s_K]*pow(temp_m_K,4);
	temp_I_K = temp_I_K_Vs*temp1;
	temp_s_K = temp_I_K_Vs;
	if (verbose>3){ printf(" I_K %f, I_K_Vs %f I_K_m_K %f\n",temp_I_K,temp_I_K_Vs,temp_I_K_m_K);}
	if (verbose>3){ printf(" %% %% %% %% calcuating rhs_m_K: ");}
	temp1 = (temp_Vs+34)/10; temp2 = -exp(-temp1); temp3 = -temp2/10; temp2+=1;
	alpha1 = temp1/temp2/10; alpha1p = (0.1*temp2-temp1*temp3)/pow(temp2,2)/10;
	temp1 = (temp_Vs+44)/25; temp2 = 0.125*exp(-temp1); temp3 = -temp2/25;
	beta1 = temp2; beta1p = temp3;
	temp_rhs_m_K_m_K = -4*(alpha1+beta1);
	temp_rhs_m_K = 4*alpha1 + temp_m_K*temp_rhs_m_K_m_K;
	temp_rhs_m_K_Vs = 4*(alpha1p - temp_m_K*(alpha1p + beta1p));
	if (verbose>3){ printf(" rhs_m_K %f rhs_m_K_Vs %f rhs_m_K_m_K %f\n",temp_rhs_m_K,temp_rhs_m_K_Vs,temp_rhs_m_K_m_K);}
	if (verbose>3){ printf(" %% %% %% %% calculating I_SD and I_DS: ");}
	temp1 = (temp_Vs - temp_Vd);
	temp_I_SD = CONDUCTANCE_SD*temp1;
	temp_I_SD_Vs = CONDUCTANCE_SD;
	temp_I_SD_Vd = -CONDUCTANCE_SD;
	temp_I_DS = -CONDUCTANCE_DS*temp1;
	temp_I_DS_Vs = -CONDUCTANCE_DS;
	temp_I_DS_Vd = CONDUCTANCE_DS;
	if (verbose>3){ printf(" I_SD %f I_SD_Vs %f I_SD_Vd %f I_DS %f I_DS_Vs %f I_DS_Vd %f\n",temp_I_SD,temp_I_SD_Vs,temp_I_SD_Vd,temp_I_DS,temp_I_DS_Vs,temp_I_DS_Vd);}
	if (verbose>3){ printf(" %% %% %% %% calculating I_Ca: ");}
	temp1 = (temp_Vd+20)/9; temp2 = exp(-temp1); temp3 = -temp2/9; temp2+=1;
	alpha1 = 1/temp2; alpha1p = -temp3/pow(temp2,2);
	temp1 = (temp_Vd - VOLTAGE_[VARNAME_REGISTRY_s_Ca]);
	temp_I_Ca = CONDUCTANCE_[VARNAME_REGISTRY_s_Ca]*pow(alpha1,2)*temp1;
	temp_I_Ca_Vd = CONDUCTANCE_[VARNAME_REGISTRY_s_Ca]*(2*alpha1*alpha1p*temp1 + pow(alpha1,2));
	temp_s_Ca = temp_I_Ca/temp1;
	if (verbose>3){ printf(" I_Ca %f, I_Ca_Vd %f\n",temp_I_Ca,temp_I_Ca_Vd);}
	if (verbose>3){ printf(" %% %% %% %% calculating I_KCa: ");}
	temp1 = (temp_Vd - VOLTAGE_[VARNAME_REGISTRY_s_KCa]);
	Ca2KCa_power = 8;
/* 	temp2 = pow(temp_Ca/(temp_Ca+30),Ca2KCa_power); */
/* 	temp3 = Ca2KCa_power*pow(temp_Ca/(temp_Ca+30),Ca2KCa_power-1)*30/pow(temp_Ca+30,2); */
	temp2 = pow(temp_Ca/0.01,Ca2KCa_power);
	temp3 = Ca2KCa_power*pow(temp_Ca/0.05,Ca2KCa_power-1)/0.05;
	temp_I_KCa = CONDUCTANCE_[VARNAME_REGISTRY_s_KCa]*temp2*temp1;
	temp_I_KCa_Vd = CONDUCTANCE_[VARNAME_REGISTRY_s_KCa]*temp2;
	temp_I_KCa_Ca = CONDUCTANCE_[VARNAME_REGISTRY_s_KCa]*temp3*temp1;
	temp_s_KCa = temp_I_KCa/temp1;
	if (verbose>3){ printf(" I_KCa %f, I_KCa_Vd %f, I_KCa_Ca %f\n",temp_I_KCa,temp_I_KCa_Vd,temp_I_KCa_Ca);}
	if (verbose>3){ printf(" %% %% %% %% calculating rhs_Ca: ");}
	temp_rhs_Ca = -0.00175*temp_I_Ca - temp_Ca/TAU_[VARNAME_REGISTRY_Ca];
	temp_rhs_Ca_Ca = -1/TAU_[VARNAME_REGISTRY_Ca];
	temp_rhs_Ca_Vd = -0.00175*temp_I_Ca_Vd;
	if (verbose>3){ printf(" rhs_Ca %f, rhs_Ca_Vd %f, rhs_Ca_Ca %f\n",temp_rhs_Ca,temp_rhs_Ca_Vd,temp_rhs_Ca_Ca);}
	if (verbose>2){ printf(" %% %% %% %% calculating: implicit euler update\n");}
	temp_rhs_Vs = -(temp_I_leak_S+temp_I_Na+temp_I_K+temp_I_SD+CURRENT_INJECTION_S)/TAU_[VARNAME_REGISTRY_Vs];
	temp_rhs_Vs_Vs = -(temp_I_leak_S_Vs+temp_I_Na_Vs+temp_I_K_Vs+temp_I_SD_Vs)/TAU_[VARNAME_REGISTRY_Vs];
	temp_rhs_Vs_h_Na = -(temp_I_Na_h_Na)/TAU_[VARNAME_REGISTRY_Vs];
	temp_rhs_Vs_m_K = -(temp_I_K_m_K)/TAU_[VARNAME_REGISTRY_Vs];
	temp_rhs_Vs_Vd = -(temp_I_SD_Vd)/TAU_[VARNAME_REGISTRY_Vs];
	temp_rhs_Vd = -(temp_I_leak_D+temp_I_Ca+temp_I_KCa+temp_I_DS+temp_I_syn+CURRENT_INJECTION_D)/TAU_[VARNAME_REGISTRY_Vd];
	temp_rhs_Vd_Vd = -(temp_I_leak_D_Vd+temp_I_Ca_Vd+temp_I_KCa_Vd+temp_I_DS_Vd+temp_I_syn_Vd)/TAU_[VARNAME_REGISTRY_Vd];
	temp_rhs_Vd_Vs = -(temp_I_DS_Vs)/TAU_[VARNAME_REGISTRY_Vd];
	temp_rhs_Vd_Ca = -(temp_I_KCa_Ca)/TAU_[VARNAME_REGISTRY_Vd];
	specialized_inverse_adi(1-dt*temp_rhs_Vs_Vs,-dt*temp_rhs_Vs_h_Na,-dt*temp_rhs_Vs_m_K,-dt*temp_rhs_Vs_Vd,-dt*temp_rhs_h_Na_Vs,1-dt*temp_rhs_h_Na_h_Na,-dt*temp_rhs_m_K_Vs,1-dt*temp_rhs_m_K_m_K,-dt*temp_rhs_Vd_Vs,1-dt*temp_rhs_Vd_Vd,-dt*temp_rhs_Vd_Ca,-dt*temp_rhs_Ca_Vd,1-dt*temp_rhs_Ca_Ca,temp_Vs_start - temp_Vs + dt*temp_rhs_Vs,temp_h_Na_start - temp_h_Na + dt*temp_rhs_h_Na,temp_m_K_start - temp_m_K + dt*temp_rhs_m_K,temp_Vd_start - temp_Vd + dt*temp_rhs_Vd,temp_Ca_start - temp_Ca + dt*temp_rhs_Ca,&update_Vs,&update_h_Na,&update_m_K,&update_Vd,&update_Ca);
	switch (LYAPUNOV_BOTHER){
	case 0: break;
	case 1: jacobian_to_lyapunov_adi(temp_rhs_Vs_Vs,temp_rhs_Vs_h_Na,temp_rhs_Vs_m_K,temp_rhs_Vs_Vd,temp_rhs_h_Na_Vs,temp_rhs_h_Na_h_Na,temp_rhs_m_K_Vs,temp_rhs_m_K_m_K,temp_rhs_Vd_Vs,temp_rhs_Vd_Vd,temp_rhs_Vd_Ca,temp_rhs_Ca_Vd,temp_rhs_Ca_Ca,&temp_lyapunov); break;
	case 2: jacobian_to_lyapunov_morrislecar(temp_rhs_Vs_Vs,temp_rhs_Vs_h_Na,temp_rhs_Vs_m_K,temp_rhs_h_Na_Vs,temp_rhs_h_Na_h_Na,temp_rhs_m_K_Vs,temp_rhs_m_K_m_K,&temp_lyapunov); break;
	default: break;}
	norm_old = pow(update_Vs,2)+pow(update_h_Na,2)+pow(update_m_K,2)+pow(update_Vd,2)+pow(update_Ca,2);
	temp_Vs += update_Vs;
	temp_h_Na += update_h_Na;
	temp_m_K += update_m_K;
	temp_Vd += update_Vd;
	temp_Ca += update_Ca;
	temp7=norm_old;norm_old=norm_new;norm_new=temp7;
	iteration += 1;
	if (verbose>2){ printf(" %% %% %% norm_new %0.16f norm_old %0.16f\n",norm_new,norm_old);}}
      while (norm_new>norm_threshold && norm_new<norm_old && iteration<iteration_max);
      if (verbose>1){ printf(" %% %% iteration=%d,norm_new %0.16f\n",iteration,norm_new);}
      temp_s_A_start = temp_s_A;
      temp_s_N_start = temp_s_N;
      temp_s_NMDA_start = temp_s_NMDA;
      temp_s_G_start = temp_s_G;
      temp_s_ORN_start = temp_s_ORN;
      temp_h_Na_start = maximum(0,minimum(1,temp_h_Na));
      temp_m_K_start = maximum(0,minimum(1,temp_m_K));
      temp_Ca_start = temp_Ca;
      temp_Vs_start = temp_Vs;
      temp_Vd_start = temp_Vd;
      microstep += 1; t+=dt;
      if (norm_new>norm_old || iteration==iteration_max || norm_new>norm_threshold){
	if (verbose>2){ printf(" %% %% convergence poor, setting error_flag=1\n");}
	error_flag=1;}
      else{
	if (verbose>2){ printf(" %% %% convergence good, setting error_flag=0\n");}
	error_flag=0;}}
    while (microstep<n->microstep && !error_flag);
    if (error_flag){ n->microstep*=2;} else{ n->microstep=maximum(1,n->microstep/2);}
    if (verbose>2){ printf(" %% total dt %f error_flag=%d, n->microstep set to %d\n",t-t_start,error_flag,n->microstep);}}
  while(error_flag);
  temp_Vs_end = temp_Vs_start; temp_Vs_start = n->vra[VARNAME_REGISTRY_Vs];
  temp_Vd_end = temp_Vd_start; temp_Vd_start = n->vra[VARNAME_REGISTRY_Vd];
  if (verbose>2){ printf(" %% determining spiking properties\n");}
  if (n->spikenext<t_start+dt_start){
    temp1=-1;
    if (temp_Vs_end>=VOLTAGE_THRESHOLD_S && temp_Vs_start<VOLTAGE_THRESHOLD_S){
      temp1 = linerootfinder(temp_Vs_start,temp_Vs_end,VOLTAGE_THRESHOLD_S,dt_start);}
    if (temp_Vd_end>=VOLTAGE_THRESHOLD_D && temp_Vd_start<VOLTAGE_THRESHOLD_D){
      temp1 = linerootfinder(temp_Vd_start,temp_Vd_end,VOLTAGE_THRESHOLD_D,dt_start);}
    if (temp1>=0){
      n->spiketime_guess = t_start+temp1;
      if (n->spiketime_guess<=n->spikenext){
	n->spiketime_guess=-2*dt_start; n->spiketime_guess_flag=0;}
      else /* if (n->spiketime_guess>n->spikenext) */{
	n->vra[VARNAME_REGISTRY_spike_flag] += 1;
	n->spiketime_guess_flag=1;}}}
  n->vra[VARNAME_REGISTRY_spike_flag] *= exp(-dt_start/(TAU_REF/16));
  n->vra[VARNAME_REGISTRY_s_A] = temp_s_A_start;
  n->vra[VARNAME_REGISTRY_s_N] = temp_s_N_start;
  n->vra[VARNAME_REGISTRY_s_NMDA] = temp_s_NMDA_start;
  n->vra[VARNAME_REGISTRY_s_G] = temp_s_G_start;
  n->vra[VARNAME_REGISTRY_s_ORN] = temp_s_ORN_start;
  n->vra[VARNAME_REGISTRY_h_Na] = temp_h_Na_start;
  n->vra[VARNAME_REGISTRY_m_K] = temp_m_K_start;
  n->vra[VARNAME_REGISTRY_Ca] = temp_Ca_start;
  n->vra[VARNAME_REGISTRY_Vs] = temp_Vs_end;
  n->vra[VARNAME_REGISTRY_Vd] = temp_Vd_end;
  n->vra[VARNAME_REGISTRY_s_Na] = temp_s_Na;
  n->vra[VARNAME_REGISTRY_s_K] = temp_s_K;
  n->vra[VARNAME_REGISTRY_s_Ca] = temp_s_Ca;
  n->vra[VARNAME_REGISTRY_s_KCa] = temp_s_KCa;
  n->vra[VARNAME_REGISTRY_lyapunov] = temp_lyapunov;
  if (verbose>1){ printf(" %% exiting vraevolve: \n"); raprintf(n->vra,"double",1,n->nvars," end ");}
  if (verbose>0){ printf(" %% vra:\n"); for (nv=0;nv<n->nvars;nv++){ printf("%s %f \t",GLOBAL_VARNAMES[nv],n->vra[nv]);} printf("\n");}
}

void jacobian_to_lyapunov_yisun_2(double M00,double M01,double M02,double M03,double M10,double M11,double M20,double M22,double M30,double M33,double *output)
{
  /* return maximum real part of spectrum of
    X    0    1    2    3
    0  M00  M01  M02  M03
    1  M10  M11   .    .
    2  M20   .   M22   .
    3  M30   .    .   M33
  */
  int verbose=0;
  integer four=4,sixteen=16;
  doublecomplex A[16];
  doublecomplex B[4];
  integer INFO=0;
  doublecomplex WORK[16];
  doublereal RWORK[8];
  char en='N';
  int nr=0,nc=0;
  A[0 + 0*four].r = M00; A[0 + 0*four].i = 0.0;
  A[1 + 0*four].r = M10; A[1 + 0*four].i = 0.0;
  A[2 + 0*four].r = M20; A[2 + 0*four].i = 0.0;
  A[3 + 0*four].r = M30; A[3 + 0*four].i = 0.0;
  A[0 + 1*four].r = M01; A[0 + 1*four].i = 0.0;
  A[1 + 1*four].r = M11; A[1 + 1*four].i = 0.0;
  A[2 + 1*four].r = 0.0; A[2 + 1*four].i = 0.0;
  A[3 + 1*four].r = 0.0; A[3 + 1*four].i = 0.0;
  A[0 + 2*four].r = M02; A[0 + 2*four].i = 0.0;
  A[1 + 2*four].r = 0.0; A[1 + 2*four].i = 0.0;
  A[2 + 2*four].r = M22; A[2 + 2*four].i = 0.0;
  A[3 + 2*four].r = 0.0; A[3 + 2*four].i = 0.0;
  A[0 + 3*four].r = M03; A[0 + 3*four].i = 0.0;
  A[1 + 3*four].r = 0.0; A[1 + 3*four].i = 0.0;
  A[2 + 3*four].r = 0.0; A[2 + 3*four].i = 0.0;
  A[3 + 3*four].r = M33; A[3 + 3*four].i = 0.0;
  if (verbose){ printf("A=:\n");for (nr=0;nr<four;nr++){ for (nc=0;nc<four;nc++){ printf("%0.3f+%0.3fi\t",A[nr+nc*four].r,A[nr+nc*four].i);} printf("\n");}}
  zgeev_(&en,&en,&four,A,&four,B,NULL,&four,NULL,&four,WORK,&sixteen,RWORK,&INFO);
  if (INFO!=0){ printf("warning! info=%d in jacobian_to_lyapunov_yisun\n",(int)INFO);}
  if (verbose){ printf("B=:\n");for (nr=0;nr<four;nr++){ printf("%0.3f+%0.3fi\n",B[nr].r,B[nr].i);}}
  if (output!=NULL){
    *output=B[0].r;
    for (nr=0;nr<four;nr++){ *output=maximum(*output,B[nr].r);}
    if (verbose){ printf("biggest %f\n",*output);}}
}

void jacobian_to_lyapunov_yisun_1(double M00,double M01,double M02,double M03,double M10,double M11,double M20,double M22,double M30,double M33,double *output)
{
  /* return maximum real part of spectrum of
    X    0    1    2    3
    0  M00  M01  M02  M03
    1  M10  M11   .    .
    2  M20   .   M22   .
    3  M30   .    .   M33
  */
  int verbose=1;
  int length=5;
  double *pr=NULL,*pi=NULL,*rr=NULL,*ri=NULL;
  pr = (double *)tcalloc(length,sizeof(double));
  pi = (double *)tcalloc(length,sizeof(double));
  rr = (double *)tcalloc(length,sizeof(double));
  ri = (double *)tcalloc(length,sizeof(double));
  pr[0] = 1;
  pr[1] = -(M00+M11+M22+M33);
  pr[2] = M00*M11 + M00*M22 + M00*M33 + M11*M22 + M11*M33 + M22*M33;
  pr[3] = -(M00*M11*M22 + M00*M11*M33 + M00*M22*M33 + M11*M22*M33);
  pr[4] = M00*M11*M22*M33;
  pr[2] -= 1*M10*M01;
  pr[3] -= -(M22+M33)*M10*M01;
  pr[4] -= M22*M33*M10*M01;
  pr[2] -= 1*M20*M02;
  pr[3] -= -(M11+M33)*M20*M02;
  pr[4] -= M11*M33*M20*M02;
  pr[2] -= 1*M30*M03;
  pr[3] -= -(M11+M22)*M30*M03;
  pr[4] -= M11*M22*M30*M03;
  if (durand_kerner_rootfind(length,pr,pi,rr,ri) && verbose){
    printf(" warning! roots don't converge\n");
    raprintf(pr,"double",1,length,"p");
    raprintf(rr,"double",1,length-1,"rr");
    raprintf(ri,"double",1,length-1,"ri");}
  stats("double",rr,length-1,output,NULL,NULL,NULL);
  if (pr!=NULL){ tfree(pr);pr=NULL;}
  if (pi!=NULL){ tfree(pi);pi=NULL;}
  if (rr!=NULL){ tfree(rr);rr=NULL;}
  if (ri!=NULL){ tfree(ri);ri=NULL;}
}

void specialized_inverse_yisun(double M00,double M01,double M02,double M03,double M10,double M11,double M20,double M22,double M30,double M33,double R0,double R1,double R2,double R3,double *V0,double *V1,double *V2,double *V3)
{
  /* specialized 28 flop inverse (row operations) for the matrix
    X    0    1    2    3 
    0  M00  M01  M02  M03 
    1  M10  M11   .    .  
    2  M20   .   M22   .  
    3  M30   .    .   M33 
    with RHS [R0 R1 R2 R3]
  */
  double temp=0;
  temp = M03/M33; M00 -= temp*M30; R0 -= temp*R3; /* M03=0; */
  temp = M02/M22; M00 -= temp*M20; R0 -= temp*R2; /* M02=0; */
  temp = M01/M11; M00 -= temp*M10; R0 -= temp*R1; /* M01=0; */
  /* now matrix should look like
    X    0    1    2    3 
    0  M00   .    .    .  
    1  M10  M11   .    .  
    2  M20   .   M22   .  
    3  M30   .    .   M33 
  */
  temp = M10/M00; R1 -= temp*R0; /* M10=0; */
  temp = M20/M00; R2 -= temp*R0; /* M20=0; */
  temp = M30/M00; R3 -= temp*R0; /* M30=0; */
  /* now matrix should be diagonal */
  *V0 = R0/M00;
  *V1 = R1/M11;
  *V2 = R2/M22;
  *V3 = R3/M33;
}

void vraevolve_yisun(struct neuron *n,double t_start,double dt_start)
{
  /* linearly implicit euler, can use hhlib */
  int verbose=0;
  double temp1=0,temp2=0,temp3=0,temp7=0,temp_lyapunov=0;
  double alpha1=0,alpha1p=0,beta1=0,beta1p=0;
  double temp_s_A=0;
  double temp_s_N=0;
  double temp_s_NMDA=0;
  double temp_s_G=0;
  double temp_s_ORN=0;
  double temp_Vs=0;
  double temp_m_Na=0;
  double temp_h_Na=0;
  double temp_m_K=0;
  double temp_s_Na=0;
  double temp_s_K=0;
  double temp_I_syn=0;
  double temp_I_syn_Vs=0;
  double temp_I_leak_S=0;
  double temp_I_leak_S_Vs=0;
  double temp_I_Na=0;
  double temp_I_Na_Vs=0;
  double temp_I_Na_m_Na=0;
  double temp_I_Na_h_Na=0;
  double temp_rhs_m_Na=0;
  double temp_rhs_m_Na_Vs=0;
  double temp_rhs_m_Na_m_Na=0;
  double temp_rhs_h_Na=0;
  double temp_rhs_h_Na_Vs=0;
  double temp_rhs_h_Na_h_Na=0;
  double temp_I_K=0;
  double temp_I_K_Vs=0;
  double temp_I_K_m_K=0;
  double temp_rhs_m_K=0;
  double temp_rhs_m_K_Vs=0;
  double temp_rhs_m_K_m_K=0;
  double temp_rhs_Vs=0;
  double temp_rhs_Vs_Vs=0;
  double temp_rhs_Vs_m_Na=0;
  double temp_rhs_Vs_h_Na=0;
  double temp_rhs_Vs_m_K=0;
  double temp_s_A_start=0;
  double temp_s_N_start=0;
  double temp_s_NMDA_start=0;
  double temp_s_G_start=0;
  double temp_s_ORN_start=0;
  double temp_Vs_start=0;
  double temp_m_Na_start=0;
  double temp_h_Na_start=0;
  double temp_m_K_start=0;
  double temp_Vs_end=0;
  double dt=dt_start,t=t_start;
  int microstep=0,error_flag=0;
  int iteration=0,iteration_max=4;
  double update_Vs=0;
  double update_m_Na=0;
  double update_h_Na=0;
  double update_m_K=0;
  double norm_new=0,norm_old=0,norm_threshold=0.0000001;
  int nv=0,tab=0;
  double *vra=NULL;
  int non_refracted_flag=1;
  struct hhlib *hl=GLOBAL_HHLIB;
  double refracted_time=0;
  if (verbose>1){ printf(" %% [entering vraevolve_yisun]: \n"); raprintf(n->vra,"double",1,n->nvars," start ");}
  if (hl!=NULL){ refracted_time = hl->trigger_rarara[n->type][n->index][1]+hl->tau_skip;}
  if (hl!=NULL && hl->trigger_rarara[n->type][n->index][0]==2 && refracted_time>=t_start+dt_start){
    if (verbose){ printf(" %% fully refracted, evolving %0.16f --> %0.16f\n",t_start,t_start+dt_start);}
    temp_s_A_start = n->vra[VARNAME_REGISTRY_s_A];
    temp_s_N_start = n->vra[VARNAME_REGISTRY_s_N];
    temp_s_NMDA_start = n->vra[VARNAME_REGISTRY_s_NMDA];
    temp_s_G_start = n->vra[VARNAME_REGISTRY_s_G];
    temp_s_ORN_start = n->vra[VARNAME_REGISTRY_s_ORN];
    temp_Vs_start=n->vra[VARNAME_REGISTRY_Vs];
    temp_s_A = temp_s_A_start * exp(-dt_start/TAU_[VARNAME_REGISTRY_s_A]);
    temp_s_N = temp_s_N_start * exp(-dt_start/TAU_[VARNAME_REGISTRY_s_N]);
    temp_s_NMDA = exp(-dt_start/TAU_[VARNAME_REGISTRY_s_NMDA])*(temp_s_NMDA_start + temp_s_N_start*(-1 + exp(dt_start/(1.0/TAU_[VARNAME_REGISTRY_s_NMDA]-1.0/TAU_[VARNAME_REGISTRY_s_N])))/(1.0/TAU_[VARNAME_REGISTRY_s_NMDA]-1.0/TAU_[VARNAME_REGISTRY_s_N]));
    temp_s_G = temp_s_G_start * exp(-dt_start/TAU_[VARNAME_REGISTRY_s_G]);
    temp_s_ORN = temp_s_ORN_start * exp(-dt_start/TAU_[VARNAME_REGISTRY_s_ORN]);
    temp_Vs = VOLTAGE_[VARNAME_REGISTRY_s_LEAK_S];
    temp_s_A_start = temp_s_A;
    temp_s_N_start = temp_s_N;
    temp_s_NMDA_start = temp_s_NMDA;
    temp_s_G_start = temp_s_G;
    temp_s_ORN_start = temp_s_ORN;
    temp_Vs_start = temp_Vs;
    t_start=t_start+dt_start;
    dt_start=0;
    temp_Vs_end = temp_Vs_start; temp_Vs_start = n->vra[VARNAME_REGISTRY_Vs];
    n->vra[VARNAME_REGISTRY_spike_flag] *= exp(-dt_start/(TAU_REF/16));
    n->vra[VARNAME_REGISTRY_s_A] = temp_s_A_start;
    n->vra[VARNAME_REGISTRY_s_N] = temp_s_N_start;
    n->vra[VARNAME_REGISTRY_s_NMDA] = temp_s_NMDA_start;
    n->vra[VARNAME_REGISTRY_s_G] = temp_s_G_start;
    n->vra[VARNAME_REGISTRY_s_ORN] = temp_s_ORN_start;
    n->vra[VARNAME_REGISTRY_m_Na] = 0;
    n->vra[VARNAME_REGISTRY_h_Na] = 0;
    n->vra[VARNAME_REGISTRY_m_K] = 0;
    n->vra[VARNAME_REGISTRY_Vs] = temp_Vs_end;
    n->vra[VARNAME_REGISTRY_Vd] = temp_Vs_end;
    n->vra[VARNAME_REGISTRY_s_Na] = 0;
    n->vra[VARNAME_REGISTRY_s_K] = 0;
    n->vra[VARNAME_REGISTRY_lyapunov] = 0;
    non_refracted_flag=0;}
  else if (hl!=NULL && hl->trigger_rarara[n->type][n->index][0]==2 && refracted_time>=t_start && refracted_time < t_start+dt_start){
    if (verbose){ printf(" %% partially refracted, evolving %0.16f --> %0.16f\n",t_start,n->spikenext);}
    temp_s_A_start = n->vra[VARNAME_REGISTRY_s_A];
    temp_s_N_start = n->vra[VARNAME_REGISTRY_s_N];
    temp_s_NMDA_start = n->vra[VARNAME_REGISTRY_s_NMDA];
    temp_s_G_start = n->vra[VARNAME_REGISTRY_s_G];
    temp_s_ORN_start = n->vra[VARNAME_REGISTRY_s_ORN];
    temp_Vs_start=n->vra[VARNAME_REGISTRY_Vs];
    temp_s_A = temp_s_A_start * exp(-(refracted_time-t_start)/TAU_[VARNAME_REGISTRY_s_A]);
    temp_s_N = temp_s_N_start * exp(-(refracted_time-t_start)/TAU_[VARNAME_REGISTRY_s_N]);
    temp_s_NMDA = exp(-(refracted_time-t_start)/TAU_[VARNAME_REGISTRY_s_NMDA])*(temp_s_NMDA_start + temp_s_N_start*(-1 + exp((refracted_time-t_start)/(1.0/TAU_[VARNAME_REGISTRY_s_NMDA]-1.0/TAU_[VARNAME_REGISTRY_s_N])))/(1.0/TAU_[VARNAME_REGISTRY_s_NMDA]-1.0/TAU_[VARNAME_REGISTRY_s_N]));
    temp_s_G = temp_s_G_start * exp(-(refracted_time-t_start)/TAU_[VARNAME_REGISTRY_s_G]);
    temp_s_ORN = temp_s_ORN_start * exp(-(refracted_time-t_start)/TAU_[VARNAME_REGISTRY_s_ORN]);
    temp_Vs = VOLTAGE_[VARNAME_REGISTRY_s_LEAK_S];
    temp_s_A_start = temp_s_A;
    temp_s_N_start = temp_s_N;
    temp_s_NMDA_start = temp_s_NMDA;
    temp_s_G_start = temp_s_G;
    temp_s_ORN_start = temp_s_ORN;
    temp_Vs_start = temp_Vs;
    temp_Vs_end = temp_Vs_start; temp_Vs_start = n->vra[VARNAME_REGISTRY_Vs];
    n->vra[VARNAME_REGISTRY_spike_flag] *= exp(-(refracted_time-t_start)/(TAU_REF/16));
    n->vra[VARNAME_REGISTRY_s_A] = temp_s_A_start;
    n->vra[VARNAME_REGISTRY_s_N] = temp_s_N_start;
    n->vra[VARNAME_REGISTRY_s_NMDA] = temp_s_NMDA_start;
    n->vra[VARNAME_REGISTRY_s_G] = temp_s_G_start;
    n->vra[VARNAME_REGISTRY_s_ORN] = temp_s_ORN_start;
    n->vra[VARNAME_REGISTRY_m_Na] = 0;
    n->vra[VARNAME_REGISTRY_h_Na] = 0;
    n->vra[VARNAME_REGISTRY_m_K] = 0;
    n->vra[VARNAME_REGISTRY_Vs] = temp_Vs_end;
    n->vra[VARNAME_REGISTRY_Vd] = temp_Vs_end;
    n->vra[VARNAME_REGISTRY_s_Na] = 0;
    n->vra[VARNAME_REGISTRY_s_K] = 0;
    n->vra[VARNAME_REGISTRY_lyapunov] = 0;
    if (hl!=NULL){ 
      vra = (double *) tcalloc(hl->hlhra[n->type]->nvars,sizeof(double));
      hhlib_histmean(hl->hlhra[n->type],&(hl->trigger_rarara[n->type][n->index][2]),vra);
      if (verbose>0){ 
	printf(" %% library for neuron (%d,%d) at time (%f,%f,%f), replacing:\n",n->type,n->index,t_start,refracted_time,t_start+dt_start);
	for (nv=0;nv<hl->indexing_nvar_length;nv++){ 
	  tab = hl->indexing_nvar_checkout[nv];
	  printf(" %% %% %s = %f with %f\n",GLOBAL_VARNAMES[tab],n->vra[tab],vra[nv]);}}
      for (nv=0;nv<hl->indexing_nvar_length;nv++){ n->vra[hl->indexing_nvar_checkout[nv]] = vra[nv];}
      tfree(vra);vra=NULL;}
    dt_start = (t_start+dt_start)-refracted_time;
    t_start = refracted_time;
    non_refracted_flag=1;}
  if (non_refracted_flag){
    if (verbose>1){ printf(" %% non refracted, evolving %0.16f --> %0.16f\n",t_start,t_start+dt_start);}
    do{ /* outer loop, choose n->microstep */
      error_flag=0;dt=dt_start/n->microstep;t=t_start;
      if (verbose>1){ printf(" %% initiating pass -- n->microstep %d dt %f\n",n->microstep,dt);}
      temp_s_A_start = n->vra[VARNAME_REGISTRY_s_A];
      temp_s_N_start = n->vra[VARNAME_REGISTRY_s_N];
      temp_s_NMDA_start = n->vra[VARNAME_REGISTRY_s_NMDA];
      temp_s_G_start = n->vra[VARNAME_REGISTRY_s_G];
      temp_s_ORN_start = n->vra[VARNAME_REGISTRY_s_ORN];
      temp_Vs_start=n->vra[VARNAME_REGISTRY_Vs];
      temp_m_Na_start=n->vra[VARNAME_REGISTRY_m_Na];
      temp_h_Na_start=n->vra[VARNAME_REGISTRY_h_Na];
      temp_m_K_start=n->vra[VARNAME_REGISTRY_m_K];
      microstep=0;
      do{ /* each microstep */
	if (verbose>1){ printf(" %% %% microstep %d\n",microstep);}
	norm_old=16.0;norm_new=16.0;
	temp_s_A = temp_s_A_start * exp(-dt/TAU_[VARNAME_REGISTRY_s_A]);
	temp_s_N = temp_s_N_start * exp(-dt/TAU_[VARNAME_REGISTRY_s_N]);
	temp_s_NMDA = exp(-dt/TAU_[VARNAME_REGISTRY_s_NMDA])*(temp_s_NMDA_start + temp_s_N_start*(-1 + exp(dt*(1.0/TAU_[VARNAME_REGISTRY_s_NMDA]-1.0/TAU_[VARNAME_REGISTRY_s_N])))/(1.0/TAU_[VARNAME_REGISTRY_s_NMDA]-1.0/TAU_[VARNAME_REGISTRY_s_N]));
	temp_s_G = temp_s_G_start * exp(-dt/TAU_[VARNAME_REGISTRY_s_G]);
	temp_s_ORN = temp_s_ORN_start * exp(-dt/TAU_[VARNAME_REGISTRY_s_ORN]);
	temp_m_Na = temp_m_Na_start;
	temp_h_Na = temp_h_Na_start;
	temp_m_K = temp_m_K_start;
	temp_Vs = temp_Vs_start;
	iteration=0;
	do{ /* inner loop, linearly implicit euler */
	  if (verbose>2){ printf(" %% %% %% iteration %d\n",iteration);}
	  if (verbose>3){ printf(" %% %% %% %% calculating I_syn: ");}
	  temp_I_syn = temp_s_ORN*(temp_Vs - VOLTAGE_[VARNAME_REGISTRY_s_ORN]) + temp_s_A*(temp_Vs - VOLTAGE_[VARNAME_REGISTRY_s_A]) + temp_s_NMDA*(temp_Vs - VOLTAGE_[VARNAME_REGISTRY_s_NMDA]) + temp_s_G*(temp_Vs - VOLTAGE_[VARNAME_REGISTRY_s_G]);
	  temp_I_syn_Vs = temp_s_ORN + temp_s_A + temp_s_NMDA + temp_s_G;
	  if (verbose>3){ printf(" I_syn %f\n",temp_I_syn);}
	  if (verbose>3){ printf(" %% %% %% %% calculating I_leak_S : ");}
	  temp_I_leak_S_Vs = CONDUCTANCE_[VARNAME_REGISTRY_Vs];
	  temp_I_leak_S = temp_I_leak_S_Vs*(temp_Vs - VOLTAGE_[VARNAME_REGISTRY_s_LEAK_S]);
	  if (verbose>3){ printf(" I_leak_S %f, I_leak_S_Vs %f \n",temp_I_leak_S,temp_I_leak_S_Vs);}
	  if (verbose>3){ printf(" %% %% %% %% calculating I_Na: ");}
	  temp1 = (temp_Vs - VOLTAGE_[VARNAME_REGISTRY_s_Na]);
	  temp_I_Na_h_Na = CONDUCTANCE_[VARNAME_REGISTRY_s_Na]*pow(temp_m_Na,3)*temp1;
	  temp_I_Na_m_Na = CONDUCTANCE_[VARNAME_REGISTRY_s_Na]*3*pow(temp_m_Na,2)*temp_h_Na*temp1;
	  temp_I_Na = temp_I_Na_h_Na*temp_h_Na;
	  temp_I_Na_Vs = CONDUCTANCE_[VARNAME_REGISTRY_s_Na]*pow(temp_m_Na,3)*temp_h_Na;
	  temp_s_Na = temp_I_Na_Vs;
	  if (verbose>3){ printf(" m_Na %f, h_Na %f, I_Na %f, I_Na_Vs %f I_Na_m_Na %f I_Na_h_Na %f\n",temp_m_Na,temp_h_Na,temp_I_Na,temp_I_Na_Vs,temp_I_Na_m_Na,temp_I_Na_h_Na);}
	  if (verbose>3){ printf(" %% %% %% %% calculating rhs_m_Na: ");}
	  temp1 = (temp_Vs+40)/10; temp2 = -exp(-temp1); temp3 = -temp2/10; temp2+=1;
	  alpha1 = temp1/temp2; alpha1p = (0.1*temp2-temp1*temp3)/pow(temp2,2);
	  temp1 = (temp_Vs+65)/18; temp2 = 4*exp(-temp1); temp3 = -temp2/18;
	  beta1 = temp2; beta1p = temp3;
	  temp_rhs_m_Na_m_Na = -1*(alpha1+beta1);
	  temp_rhs_m_Na = 1*alpha1 + temp_m_Na*temp_rhs_m_Na_m_Na;
	  temp_rhs_m_Na_Vs = 1*(alpha1p - temp_m_Na*(alpha1p + beta1p));
	  if (verbose>3){ printf(" rhs_m_Na %f rhs_m_Na_Vs %f rhs_m_Na_m_Na %f\n",temp_rhs_m_Na,temp_rhs_m_Na_Vs,temp_rhs_m_Na_m_Na);}
	  if (verbose>3){ printf(" %% %% %% %% calculating rhs_h_Na: ");}
	  temp1 = (temp_Vs+65)/20; temp2 = 0.07*exp(-temp1); temp3 = -temp2/20;
	  alpha1 = temp2; alpha1p = temp3;
	  temp1 = (temp_Vs+35)/10; temp2 = exp(-temp1); temp3 = -temp2/10; temp2+=1;
	  beta1 = 1/temp2; beta1p = -temp3/pow(temp2,2);
	  temp_rhs_h_Na_h_Na = -1*(alpha1+beta1);
	  temp_rhs_h_Na = 1*alpha1 + temp_h_Na*temp_rhs_h_Na_h_Na;
	  temp_rhs_h_Na_Vs = 1*(alpha1p - temp_h_Na*(alpha1p + beta1p));
	  if (verbose>3){ printf(" rhs_h_Na %f rhs_h_Na_Vs %f rhs_h_Na_h_Na %f\n",temp_rhs_h_Na,temp_rhs_h_Na_Vs,temp_rhs_h_Na_h_Na);}
	  if (verbose>3){ printf(" %% %% %% %% calculating I_K: ");}
	  temp1 = (temp_Vs - VOLTAGE_[VARNAME_REGISTRY_s_K]);
	  temp_I_K_m_K = CONDUCTANCE_[VARNAME_REGISTRY_s_K]*4*pow(temp_m_K,3)*temp1;
	  temp_I_K_Vs = CONDUCTANCE_[VARNAME_REGISTRY_s_K]*pow(temp_m_K,4);
	  temp_I_K = temp_I_K_Vs*temp1;
	  temp_s_K = temp_I_K_Vs;
	  if (verbose>3){ printf(" I_K %f, I_K_Vs %f I_K_m_K %f\n",temp_I_K,temp_I_K_Vs,temp_I_K_m_K);}
	  if (verbose>3){ printf(" %% %% %% %% calcuating rhs_m_K: ");}
	  temp1 = (temp_Vs+55)/10; temp2 = -exp(-temp1); temp3 = -temp2/10; temp2+=1;
	  alpha1 = temp1/temp2/10; alpha1p = (0.1*temp2-temp1*temp3)/pow(temp2,2)/10;
	  temp1 = (temp_Vs+65)/80; temp2 = 0.125*exp(-temp1); temp3 = -temp2/80;
	  beta1 = temp2; beta1p = temp3;
	  temp_rhs_m_K_m_K = -1*(alpha1+beta1);
	  temp_rhs_m_K = 1*alpha1 + temp_m_K*temp_rhs_m_K_m_K;
	  temp_rhs_m_K_Vs = 1*(alpha1p - temp_m_K*(alpha1p + beta1p));
	  if (verbose>3){ printf(" rhs_m_K %f rhs_m_K_Vs %f rhs_m_K_m_K %f\n",temp_rhs_m_K,temp_rhs_m_K_Vs,temp_rhs_m_K_m_K);}
	  if (verbose>2){ printf(" %% %% %% %% calculating: implicit euler update\n");}
	  temp_rhs_Vs = -(temp_I_leak_S+temp_I_Na+temp_I_K+temp_I_syn+CURRENT_INJECTION_S)/TAU_[VARNAME_REGISTRY_Vs];
	  temp_rhs_Vs_Vs = -(temp_I_leak_S_Vs+temp_I_Na_Vs+temp_I_K_Vs+temp_I_syn_Vs)/TAU_[VARNAME_REGISTRY_Vs];
	  temp_rhs_Vs_m_Na = -(temp_I_Na_m_Na)/TAU_[VARNAME_REGISTRY_Vs];
	  temp_rhs_Vs_h_Na = -(temp_I_Na_h_Na)/TAU_[VARNAME_REGISTRY_Vs];
	  temp_rhs_Vs_m_K = -(temp_I_K_m_K)/TAU_[VARNAME_REGISTRY_Vs];
	  specialized_inverse_yisun(1-dt*temp_rhs_Vs_Vs,-dt*temp_rhs_Vs_m_Na,-dt*temp_rhs_Vs_h_Na,-dt*temp_rhs_Vs_m_K,-dt*temp_rhs_m_Na_Vs,1-dt*temp_rhs_m_Na_m_Na,-dt*temp_rhs_h_Na_Vs,1-dt*temp_rhs_h_Na_h_Na,-dt*temp_rhs_m_K_Vs,1-dt*temp_rhs_m_K_m_K,temp_Vs_start - temp_Vs + dt*temp_rhs_Vs,temp_m_Na_start - temp_m_Na + dt*temp_rhs_m_Na,temp_h_Na_start - temp_h_Na + dt*temp_rhs_h_Na,temp_m_K_start - temp_m_K + dt*temp_rhs_m_K,&update_Vs,&update_m_Na,&update_h_Na,&update_m_K);
	  norm_old = pow(update_Vs,2)+pow(update_m_Na,2)+pow(update_h_Na,2)+pow(update_m_K,2);
	  temp_Vs += update_Vs;
	  temp_m_Na += update_m_Na;
	  temp_h_Na += update_h_Na;
	  temp_m_K += update_m_K;
	  temp7=norm_old;norm_old=norm_new;norm_new=temp7;
	  iteration += 1;
	  if (verbose>2){ printf(" %% %% %% norm_new %0.16f norm_old %0.16f\n",norm_new,norm_old);}}
	while (norm_new>norm_threshold && norm_new<norm_old && iteration<iteration_max);
	if (verbose>1){ printf(" %% %% iteration=%d,norm_new %0.16f\n",iteration,norm_new);}
	temp_s_A_start = temp_s_A;
	temp_s_N_start = temp_s_N;
	temp_s_NMDA_start = temp_s_NMDA;
	temp_s_G_start = temp_s_G;
	temp_s_ORN_start = temp_s_ORN;
	temp_m_Na_start = maximum(0,minimum(1,temp_m_Na));
	temp_h_Na_start = maximum(0,minimum(1,temp_h_Na));
	temp_m_K_start = maximum(0,minimum(1,temp_m_K));
	temp_Vs_start = temp_Vs;
	microstep += 1; t+=dt;
	if (norm_new>norm_old || iteration==iteration_max || norm_new>norm_threshold){
	  if (verbose>2){ printf(" %% %% convergence poor, setting error_flag=1\n");}
	  error_flag=1;}
	else{
	  if (verbose>2){ printf(" %% %% convergence good, setting error_flag=0\n");}
	  error_flag=0;}}
      while (microstep<n->microstep && !error_flag);
      if (error_flag){ n->microstep*=2;} else{ n->microstep=maximum(1,n->microstep/2);}
      if (verbose>2){ printf(" %% total dt %f error_flag=%d, n->microstep set to %d\n",t-t_start,error_flag,n->microstep);}}
    while(error_flag);
    temp_Vs_end = temp_Vs_start; temp_Vs_start = n->vra[VARNAME_REGISTRY_Vs];
    if (verbose>2){ printf(" %% determining spiking properties\n");}
    if (n->spikenext<t_start+dt_start){
      temp1=-1;
      if (temp_Vs_end>=VOLTAGE_THRESHOLD_S && temp_Vs_start<VOLTAGE_THRESHOLD_S){
	temp1 = linerootfinder(temp_Vs_start,temp_Vs_end,VOLTAGE_THRESHOLD_S,dt_start);}
      if (temp1>=0){
	n->spiketime_guess = t_start+temp1;
	if (n->spiketime_guess<=n->spikenext){
	  n->spiketime_guess=-2*dt_start; n->spiketime_guess_flag=0;}
	else /* if (n->spiketime_guess>n->spikenext) */{
	  n->vra[VARNAME_REGISTRY_spike_flag] += 1;
	  n->spiketime_guess_flag=1;}}}
    switch (LYAPUNOV_BOTHER){ /* use stale values from final microstep */
    case 0: break;
    case 1: jacobian_to_lyapunov_yisun_2(temp_rhs_Vs_Vs,temp_rhs_Vs_m_Na,temp_rhs_Vs_h_Na,temp_rhs_Vs_m_K,temp_rhs_m_Na_Vs,temp_rhs_m_Na_m_Na,temp_rhs_h_Na_Vs,temp_rhs_h_Na_h_Na,temp_rhs_m_K_Vs,temp_rhs_m_K_m_K,&temp_lyapunov); break;
    default: break;}
    n->vra[VARNAME_REGISTRY_spike_flag] *= exp(-dt_start/(TAU_REF/16));
    n->vra[VARNAME_REGISTRY_s_A] = temp_s_A_start;
    n->vra[VARNAME_REGISTRY_s_N] = temp_s_N_start;
    n->vra[VARNAME_REGISTRY_s_NMDA] = temp_s_NMDA_start;
    n->vra[VARNAME_REGISTRY_s_G] = temp_s_G_start;
    n->vra[VARNAME_REGISTRY_s_ORN] = temp_s_ORN_start;
    n->vra[VARNAME_REGISTRY_m_Na] = temp_m_Na_start;
    n->vra[VARNAME_REGISTRY_h_Na] = temp_h_Na_start;
    n->vra[VARNAME_REGISTRY_m_K] = temp_m_K_start;
    n->vra[VARNAME_REGISTRY_Vs] = temp_Vs_end;
    n->vra[VARNAME_REGISTRY_Vd] = temp_Vs_end;
    n->vra[VARNAME_REGISTRY_s_Na] = temp_s_Na;
    n->vra[VARNAME_REGISTRY_s_K] = temp_s_K;
    n->vra[VARNAME_REGISTRY_lyapunov] = temp_lyapunov;}
  if (verbose>1){ printf(" %% exiting vraevolve: \n"); raprintf(n->vra,"double",1,n->nvars," end ");}
  if (verbose>0){ printf(" %% vra:\n"); for (nv=0;nv<n->nvars;nv++){ printf("%s %f \t",GLOBAL_VARNAMES[nv],n->vra[nv]);} printf("\n");}
}

void gizmointegrateeiforif(struct neuron *n,double previous_spiketime,double next_spiketime,double *sra)
{
 /* this carries out integration over spike-free time interval d */
  int verbose=0;//(n->index==0);
  double dt=0;
  if (verbose){ 
    printf("%% [entering gizmointegrateeiforif] with n (%d,%d):",n->type,n->index); 
    raprintf(n->vra,"double",1,n->nvars," "); 
    raprintf(sra,"double",1,GLOBAL_INDEXING_sra_LENGTH," adding");}
  dt = next_spiketime - previous_spiketime;
  switch (GLOBAL_NEURON_MODEL){
  case 0: vraevolve_if(n,previous_spiketime,dt); break;
  case 1: vraevolve_adi(n,previous_spiketime,dt); break;
  case 2: vraevolve_yisun(n,previous_spiketime,dt); break;
  case 3: vraevolve_morrislecar(n,previous_spiketime,dt); break;
  case 4: printf(" %% warning! cannot use typical framework for current_based integrate and fire!\n"); break;
  default: vraevolve_adi(n,previous_spiketime,dt); break;}
  sra_add(n,sra);
  n->vra[VARNAME_REGISTRY_inputrate] = n->inputrate;
  if (verbose){ 
    printf("%% ending with n (%d,%d) --> spike %d:",n->type,n->index,n->spiketime_guess_flag); 
    raprintf(n->vra,"double",1,n->nvars," ");}
}

void clumpcorrect(struct llist *L1,struct llist *L2,double t,double DT)
{
  /* We speed up integration by only accounting for spikes which actually connect *L1 and *L2 */
  int verbose=0;
  int finished_with_input=0,number_of_spikes_left=0,perform_integration=0;
  int nv=0;
  double previous_integration_time=0,next_integration_time=0;
  struct neuron *n=NULL,*s=NULL;
  struct litem *lneuron=NULL,*lspike=NULL;
  double epsilon = 0.0000000001;
  double *sra=NULL;
  if (verbose){ printf(" %% [entering clumpcorrect]\n");}
  lneuron = L1->first;
  while (lneuron != NULL){
    n = (struct neuron *) lneuron->item;
    sra = (double *) tcalloc(GLOBAL_INDEXING_sra_LENGTH,sizeof(double));
    for (nv=0;nv<n->nvars;nv++){ n->vra[nv] = *(n->vpra[nv]);}
    /* verbose=(n->index==0); */
    if (verbose){ printf("neuron %d time %f, dt %f:",n->index,t,DT); raprintf(n->vra,"double",1,n->nvars," ");}
    n->spiketime_guess = n->spiketime;
    n->spiketime_guess_flag = 0;
    finished_with_input = !(n->spikeinput_flag && n->spikeinput_time >= t && n->spikeinput_time <= t+DT);
    number_of_spikes_left = L2->length;
    previous_integration_time = t;
    next_integration_time = t;
    perform_integration=0;
    lspike = L2->first;
    while (!finished_with_input || number_of_spikes_left>0 || previous_integration_time < t+DT-epsilon){
      if (verbose){ printf(" fwi=%d,nos=%d,pit=%f<t+DT=%f\n",finished_with_input,number_of_spikes_left,previous_integration_time,t+DT);}
      if (lspike==NULL){
	if (verbose){ printf("  lspike==NULL\n");}
	number_of_spikes_left=0;
	if (!finished_with_input && n->spikeinput_time >= previous_integration_time && n->spikeinput_time <= t+DT){
	  if (verbose){ printf("   but still have input at pit=%f,%f,%f=t+DT\n",previous_integration_time,n->spikeinput_time,t+DT);}
	  finished_with_input=1;
	  next_integration_time = n->spikeinput_time;
	  perform_integration = ilink(n,sra);}
	else /* (if finished_with_input || n->spikeinput_time < previous_integration_time || n->spikeinput_time > t+DT) */{
	  if (verbose){ printf("   but no more input\n");}
	  finished_with_input=1;
	  ilink(NULL,sra);
	  next_integration_time = t+DT;
	  perform_integration=1;}}
      else if (lspike!=NULL){
	if (verbose){ printf("  lspike exists\n");}
	s = (struct neuron *) lspike->item;
	if (!finished_with_input && n->spikeinput_time >= previous_integration_time && n->spikeinput_time <= minimum(t+DT,s->spiketime)){
	  if (verbose){ printf("   but still have input at pit=%f,%f,%f=s->spk\n",previous_integration_time,n->spikeinput_time,s->spiketime);}
	  finished_with_input=1;
	  next_integration_time = n->spikeinput_time;
	  perform_integration = ilink(n,sra);}
	else /* (if finished_with_input || n->spikeinput_time < previous_integration_time || n->spikeinput_time > minimum(t+DT,s->spiketime)) */{
	  if (verbose){ printf("   no input before spiketime\n");}
	  if (s->spiketime >= previous_integration_time && s->spiketime <= t+DT){
	    if (verbose){ printf("    pit=%f < spk=%f < t+DT=%f\n",previous_integration_time,s->spiketime,t+DT);}
	    next_integration_time = s->spiketime;
	    perform_integration = slink(s,n,sra,NULL);
	    number_of_spikes_left -= 1;
	    lspike = lspike->child;}
	  else /* (if s->spiketime < previous_integration_time || s->spiketime > t+DT) */{
	    if (verbose){ printf("    spiketime %f out of bounds pit=%f,t+DT=%f\n",s->spiketime,previous_integration_time,t+DT);}
	    number_of_spikes_left=0;
	    if (!finished_with_input && n->spikeinput_time >= previous_integration_time && n->spikeinput_time <= t+DT){
	      if (verbose){ printf("     not yet finished with input pit=%f,%f,%f=t+DT\n",previous_integration_time,n->spikeinput_time,t+DT);}
	      finished_with_input=1;
	      next_integration_time = n->spikeinput_time;
	      perform_integration = ilink(n,sra);}
	    else /* (if finished_with_input || n->spikeinput_time < previous_integration_time || n->spikeinput_time > t+DT) */{
	      if (verbose){ printf("     no more input, we are exhausted\n");}
	      finished_with_input=1;
	      ilink(NULL,sra);
	      next_integration_time = t+DT;
	      perform_integration=1;}}}}
      if (perform_integration==1){ 
	if (verbose){ printf("      -integrating %f %f\n",previous_integration_time,next_integration_time);}
	gizmointegrateeiforif(n,previous_integration_time,next_integration_time,sra);
	previous_integration_time = next_integration_time;}}
    if (verbose){ printf(" moving on to another neuron \n\n");}
    tfree(sra);sra=NULL;
    lneuron = lneuron->child;}
  lneuron = L1->first;
  while (lneuron != NULL){
    n = (struct neuron*) lneuron->item;
    if (n->spiketime_guess_flag){ n->spiketime = n->spiketime_guess;}
    else{ n->spiketime = t+2*DT;}
    lneuron = lneuron->child;}
}

void spikecorrect(struct llist *LS,double t,double DT)
{
  llistsort(LS->first,LS->last,LS->length,&spiketime_compare);
  clumpcorrect(LS,LS,t,DT);     
  llistsort(LS->first,LS->last,LS->length,&spiketime_compare);
  clumpcorrect(LS,LS,t,DT);     
  llistsort(LS->first,LS->last,LS->length,&spiketime_compare);
}

void spikeconduct(struct llist *LS,double t,double DT,int *totalspikes,int *totalspikes_by_type)
{
  /* presumes totalspikes and totalspikes_by_type[GLOBAL_NTYPES] are initialized to zero */
  struct neuron *n=NULL;
  struct litem *ln=NULL;
  int nv=0;
  double vtemp=0;
  ln=LS->first;
  while (ln!=NULL){
    n=(struct neuron *)ln->item;
    for (nv=0;nv<n->nvars;nv++){ vtemp = *(n->vpra[nv]); *(n->vpra[nv]) = n->vra[nv]; n->vra[nv] = vtemp;}
    if (n->spiketime >= t && n->spiketime <= t+DT){
      if (totalspikes!=NULL){ *totalspikes += 1;}
      if (totalspikes_by_type!=NULL){ totalspikes_by_type[n->type] += 1;}
      if (n->type==TYPENAME_REGISTRY_PN){ n->spikelast = n->spiketime; n->spikenext = n->spiketime + TAU_REF;}
      if (n->type==TYPENAME_REGISTRY_LN){ n->spikelast = n->spiketime; n->spikenext = n->spiketime + TAU_REF;}}
    ln=ln->child;}
}

double linerootfinder(double VI,double VF,double VT,double DT)
{
  double r=-0;
  if (VI<=VT && VF>=VT){ r=DT*(VT-VI)/(VF-VI);} else{ r=2*DT;} return r;
}

int gizmospiked(struct llist *LG,double t,double DT)
{
  int verbose=0;
  struct neuron *n=NULL;
  struct litem *ln=NULL;
  int spikes_found=0;
  ln=LG->first; spikes_found=0;
  while (ln!=NULL){ 
    n=(struct neuron *)ln->item; 
    if (n->spiketime >= t && n->spiketime <= t+DT){ 
      if (verbose){ printf(" %% neuron (%d,%d) fired at %0.16e<=%0.16e<=%0.16e\n",n->type,n->index,t,n->spiketime,t+DT);}
      spikes_found += 1;}
    ln=ln->child;}
  return spikes_found;
}

void gizmoconduct(struct neuronarray *Nra,struct llist *LS,struct llist *LG,double t,double DT)
{
  int verbose=1;
  struct neuron *n=NULL,*n2=NULL;
  struct litem *ln=NULL;
  int nt=0,nr=0,nv=0;
  double *sra=NULL,vtemp=0;
  ln=LG->first;
  while (ln!=NULL){
    n=(struct neuron *)ln->item;
    for (nv=0;nv<n->nvars;nv++){ vtemp = *(n->vpra[nv]); *(n->vpra[nv]) = n->vra[nv]; n->vra[nv] = vtemp;}
    if (n->spiketime >= t && n->spiketime <= t+DT){ 
      if (verbose && DT>2*GLOBAL_DTmin){ 
	printf(" %% warning! neuron (%d,%d) fired at %f<=%f<=%f, (DT=%0.16e)\n",n->type,n->index,t,n->spiketime,t+DT,DT);}
      for (nt=0;nt<Nra->ntypes;nt++){ for (nr=0;nr<Nra->lengthra[nt];nr++){
	n2 = nget(Nra,nt,nr);
	sra = (double *) tcalloc(GLOBAL_INDEXING_sra_LENGTH,sizeof(double));
	slink(n,n2,sra,NULL);
	for (nv=0;nv<GLOBAL_INDEXING_sra_LENGTH;nv++){ *(n2->vpra[GLOBAL_INDEXING_CHECKOUT_sra[nv]]) += sra[nv];}
	tfree(sra);sra=NULL;}}
      if (n->type==TYPENAME_REGISTRY_PN){ 
	n->spikelast = n->spiketime;
	n->spikenext = n->spiketime + TAU_REF;}
      if (n->type==TYPENAME_REGISTRY_LN){
	n->spikelast = n->spiketime;
	n->spikenext = n->spiketime + TAU_REF;}}
    ln=ln->child;}
}

/* compte functions */

/* void specialized_inverse_compte(double M00,double M01,double M02,double M03,double M04,double M05,double M06,double M10,double M11,double M20,double M22,double M30,double M33,double M40,double M44,double M50,double M51,double M55,double M56,double M60,double M66,double M67,double M76,double M77,double R0,double R1,double R2,double R3,double R4,double R5,double R6,double R7,double *V0,double *V1,double *V2,double *V3,double *V4,double *V5,double *V6,double *V7) */
/* { */
/*   /\* specialized 72 flop inverse (row operations) for the matrix */
/*     X    0    1    2    3    4    5    6    7 */
/*     0  M00  M01  M02  M03  M04  M05  M06   . */
/*     1  M10  M11   .    .    .    .    .    . */
/*     2  M20   .   M22   .    .    .    .    . */
/*     3  M30   .    .   M33   .    .    .    .  */
/*     4  M40   .    .    .   M44   .    .    . */
/*     5  M50  M51   .    .    .   M55  M56   . */
/*     6  M60   .    .    .    .    .   M66  M67 */
/*     7   .    .    .    .    .    .   M76  M77 */
/*     with RHS [R1 R2 R3 R4 R5 R6 R7] */
/*   *\/ */
/*   double temp=0; */
/*   temp = M67/M77; M66 -= temp*M76; R6 -= temp*R7; /\* M67=0; *\/ */
/*   temp = M56/M66; M50 -= temp*M60; R5 -= temp*R6; /\* M56=0; *\/ */
/*   temp = M06/M66; M00 -= temp*M60; R0 -= temp*R6; /\* M06=0; *\/ */
/*   temp = M51/M11; M50 -= temp*M10; R5 -= temp*R1; /\* M51=0; *\/ */
/*   temp = M05/M55; M00 -= temp*M50; R0 -= temp*R5; /\* M05=0; *\/ */
/*   temp = M04/M44; M00 -= temp*M40; R0 -= temp*R4; /\* M04=0; *\/ */
/*   temp = M03/M33; M00 -= temp*M30; R0 -= temp*R3; /\* M03=0; *\/ */
/*   temp = M02/M22; M00 -= temp*M20; R0 -= temp*R2; /\* M02=0; *\/ */
/*   temp = M01/M11; M00 -= temp*M10; R0 -= temp*R1; /\* M01=0; *\/ */
/*   /\* now matrix should look like  */
/*     X    0    1    2    3    4    5    6    7 */
/*     0  M00   .    .    .    .    .    .    . */
/*     1  M10  M11   .    .    .    .    .    . */
/*     2  M20   .   M22   .    .    .    .    . */
/*     3  M30   .    .   M33   .    .    .    .  */
/*     4  M40   .    .    .   M44   .    .    . */
/*     5  M50   .    .    .    .   M55   .    . */
/*     6  M60   .    .    .    .    .   M66   .  */
/*     7   .    .    .    .    .    .   M76  M77 */
/*   *\/ */
/*   temp = M10/M00; R1 -= temp*R0; /\* M10=0; *\/ */
/*   temp = M20/M00; R2 -= temp*R0; /\* M20=0; *\/ */
/*   temp = M30/M00; R3 -= temp*R0; /\* M30=0; *\/ */
/*   temp = M40/M00; R4 -= temp*R0; /\* M40=0; *\/ */
/*   temp = M50/M00; R5 -= temp*R0; /\* M50=0; *\/ */
/*   temp = M60/M00; R6 -= temp*R0; /\* M60=0; *\/ */
/*   temp = M76/M66; R7 -= temp*R6; /\* M76=0; *\/ */
/*   /\* now matrix should be diagonal *\/ */
/*   *V0 = R0/M00; */
/*   *V1 = R1/M11; */
/*   *V2 = R2/M22; */
/*   *V3 = R3/M33; */
/*   *V4 = R4/M44; */
/*   *V5 = R5/M55; */
/*   *V6 = R6/M66; */
/*   *V7 = R7/M77; */
/* } */

/* void vraevolve_compte(struct neuron *n,double t_start,double dt_start) */
/* { */
/*   /\* linearly implicit euler *\/ */
/*   int verbose=0; */
/*   double temp1=0,temp2=0,temp3=0,temp4=0,temp5=0,temp6=0,temp7=0; */
/*   double alpha1=0,alpha1p=0,beta1=0,beta1p=0; */
/*   double temp_s_A=0;  */
/*   double temp_s_N=0;  */
/*   double temp_s_NMDA=0;  */
/*   double temp_s_G=0;  */
/*   double temp_s_ORN=0; */
/*   double temp_Vs=0;  */
/*   double temp_h_Na=0;  */
/*   double temp_m_K=0;  */
/*   double temp_h_KA=0; */
/*   double temp_m_KS=0; */
/*   double temp_Na=0; */
/*   double temp_Vd=0; */
/*   double temp_Ca=0; */
/*   double temp_s_Na=0; */
/*   double temp_s_K=0; */
/*   double temp_s_KA=0; */
/*   double temp_s_KS=0; */
/*   double temp_s_KNa=0; */
/*   double temp_s_Ca=0; */
/*   double temp_s_NaP=0; */
/*   double temp_s_AR=0; */
/*   double temp_s_KCa=0; */
/*   double temp_I_syn=0; */
/*   double temp_I_syn_Vd=0; */
/*   double temp_I_leak_S=0; */
/*   double temp_I_leak_S_Vs=0; */
/*   double temp_I_leak_D=0; */
/*   double temp_I_leak_D_Vd=0; */
/*   double temp_I_Na=0; */
/*   double temp_I_Na_Vs=0; */
/*   double temp_I_Na_h_Na=0; */
/*   double temp_rhs_h_Na=0; */
/*   double temp_rhs_h_Na_Vs=0; */
/*   double temp_rhs_h_Na_h_Na=0; */
/*   double temp_I_K=0; */
/*   double temp_I_K_Vs=0; */
/*   double temp_I_K_m_K=0; */
/*   double temp_rhs_m_K=0; */
/*   double temp_rhs_m_K_Vs=0; */
/*   double temp_rhs_m_K_m_K=0; */
/*   double temp_I_KA=0; */
/*   double temp_I_KA_Vs=0; */
/*   double temp_I_KA_h_KA=0; */
/*   double temp_rhs_h_KA=0; */
/*   double temp_rhs_h_KA_Vs=0; */
/*   double temp_rhs_h_KA_h_KA=0; */
/*   double temp_I_KS=0; */
/*   double temp_I_KS_Vs=0; */
/*   double temp_I_KS_m_KS=0; */
/*   double temp_rhs_m_KS=0; */
/*   double temp_rhs_m_KS_Vs=0; */
/*   double temp_rhs_m_KS_m_KS=0; */
/*   double temp_I_KNa=0; */
/*   double temp_I_KNa_Vs=0; */
/*   double temp_I_KNa_Na=0; */
/*   double temp_I_SD=0; */
/*   double temp_I_SD_Vs=0; */
/*   double temp_I_SD_Vd=0; */
/*   double temp_I_DS=0; */
/*   double temp_I_DS_Vs=0; */
/*   double temp_I_DS_Vd=0; */
/*   double temp_I_Ca=0; */
/*   double temp_I_Ca_Vd=0; */
/*   double temp_I_NaP=0; */
/*   double temp_I_NaP_Vd=0; */
/*   double temp_I_AR=0; */
/*   double temp_I_AR_Vd=0; */
/*   double temp_I_KCa=0; */
/*   double temp_I_KCa_Vd=0; */
/*   double temp_I_KCa_Ca=0; */
/*   double temp_rhs_Na=0; */
/*   double temp_rhs_Na_Vs=0; */
/*   double temp_rhs_Na_Vd=0; */
/*   double temp_rhs_Na_h_Na=0; */
/*   double temp_rhs_Na_Na=0; */
/*   double temp_rhs_Ca=0; */
/*   double temp_rhs_Ca_Vd=0; */
/*   double temp_rhs_Ca_Ca=0; */
/*   double temp_rhs_Vs=0; */
/*   double temp_rhs_Vs_Vs=0; */
/*   double temp_rhs_Vs_h_Na=0; */
/*   double temp_rhs_Vs_m_K=0; */
/*   double temp_rhs_Vs_h_KA=0; */
/*   double temp_rhs_Vs_m_KS=0; */
/*   double temp_rhs_Vs_Na=0; */
/*   double temp_rhs_Vs_Vd=0; */
/*   double temp_rhs_Vd=0; */
/*   double temp_rhs_Vd_Vs=0; */
/*   double temp_rhs_Vd_Vd=0; */
/*   double temp_rhs_Vd_Ca=0; */
/*   double temp_s_A_start=0; */
/*   double temp_s_N_start=0; */
/*   double temp_s_NMDA_start=0; */
/*   double temp_s_G_start=0; */
/*   double temp_s_ORN_start=0; */
/*   double temp_Vs_start=0; */
/*   double temp_h_Na_start=0; */
/*   double temp_m_K_start=0; */
/*   double temp_h_KA_start=0; */
/*   double temp_m_KS_start=0; */
/*   double temp_Na_start=0; */
/*   double temp_Vd_start=0; */
/*   double temp_Ca_start=0; */
/*   double temp_Vs_end=0; */
/*   double temp_Vd_end=0; */
/*   double dt=dt_start,t=t_start; */
/*   int microstep=0,error_flag=0; */
/*   int iteration=0,iteration_max=4; */
/*   double update_Vs=0; */
/*   double update_h_Na=0; */
/*   double update_m_K=0; */
/*   double update_h_KA=0; */
/*   double update_m_KS=0; */
/*   double update_Na=0; */
/*   double update_Vd=0; */
/*   double update_Ca=0; */
/*   double norm_new=0,norm_old=0,norm_threshold=0.0000001; */
/*   int nv=0; */
/*   if (verbose>1){ printf(" %% [entering vraevolve_compte]: \n"); raprintf(n->vra,"double",1,n->nvars," start ");} */
/*   do{ /\* outer loop, choose n->microstep *\/ */
/*     error_flag=0;dt=dt_start/n->microstep;t=t_start; */
/*     if (verbose>1){ printf(" %% initiating pass -- n->microstep %d dt %f\n",n->microstep,dt);} */
/*     temp_s_A_start = n->vra[VARNAME_REGISTRY_s_A]; */
/*     temp_s_N_start = n->vra[VARNAME_REGISTRY_s_N]; */
/*     temp_s_NMDA_start = n->vra[VARNAME_REGISTRY_s_NMDA]; */
/*     temp_s_G_start = n->vra[VARNAME_REGISTRY_s_G]; */
/*     temp_s_ORN_start = n->vra[VARNAME_REGISTRY_s_ORN]; */
/*     temp_Vs_start=n->vra[VARNAME_REGISTRY_Vs]; */
/*     temp_h_Na_start=n->vra[VARNAME_REGISTRY_h_Na]; */
/*     temp_m_K_start=n->vra[VARNAME_REGISTRY_m_K]; */
/*     temp_h_KA_start=n->vra[VARNAME_REGISTRY_h_KA]; */
/*     temp_m_KS_start=n->vra[VARNAME_REGISTRY_m_KS]; */
/*     temp_Na_start=n->vra[VARNAME_REGISTRY_Na]; */
/*     temp_Vd_start=n->vra[VARNAME_REGISTRY_Vd]; */
/*     temp_Ca_start=n->vra[VARNAME_REGISTRY_Ca]; */
/*     microstep=0; */
/*     do{ /\* each microstep *\/ */
/*       if (verbose>1){ printf(" %% %% microstep %d\n",microstep);} */
/*       norm_old=16.0;norm_new=16.0; */
/*       temp_s_A = temp_s_A_start * exp(-dt/TAU_[VARNAME_REGISTRY_s_A]); */
/*       temp_s_N = temp_s_N_start * exp(-dt/TAU_[VARNAME_REGISTRY_s_N]); */
/*       temp_s_NMDA = exp(-dt/TAU_[VARNAME_REGISTRY_s_NMDA])*(temp_s_NMDA_start + temp_s_N_start*(-1 + exp(dt/(1.0/TAU_[VARNAME_REGISTRY_s_NMDA]-1.0/TAU_[VARNAME_REGISTRY_s_N])))/(1.0/TAU_[VARNAME_REGISTRY_s_NMDA]-1.0/TAU_[VARNAME_REGISTRY_s_N])); */
/*       temp_s_G = temp_s_G_start * exp(-dt/TAU_[VARNAME_REGISTRY_s_G]); */
/*       temp_s_ORN = temp_s_ORN_start * exp(-dt/TAU_[VARNAME_REGISTRY_s_ORN]); */
/*       temp_h_Na = temp_h_Na_start; */
/*       temp_m_K = temp_m_K_start; */
/*       temp_h_KA = temp_h_KA_start; */
/*       temp_m_KS = temp_m_KS_start; */
/*       temp_Na = temp_Na_start; */
/*       temp_Ca = temp_Ca_start; */
/*       temp_Vs = temp_Vs_start; */
/*       temp_Vd = temp_Vd_start; */
/*       iteration=0; */
/*       do{ /\* inner loop, linearly implicit euler *\/ */
/* 	if (verbose>2){ printf(" %% %% %% iteration %d\n",iteration);} */
/* 	if (verbose>3){ printf(" %% %% %% %% calculating I_syn: ");} */
/* 	temp1 = exp(-0.062*(temp_Vd+VOLTAGE_[VARNAME_REGISTRY_s_LEAK_D]-60)); temp2 = temp1*-0.062; temp1+=1; temp3 = -temp2/pow(temp1,2); */
/* 	temp_I_syn = temp_s_ORN*(temp_Vd - VOLTAGE_[VARNAME_REGISTRY_s_ORN]) + temp_s_A*(temp_Vd - VOLTAGE_[VARNAME_REGISTRY_s_A]) + temp_s_NMDA/temp1*(temp_Vd - VOLTAGE_[VARNAME_REGISTRY_s_NMDA]) + temp_s_G*(temp_Vd - VOLTAGE_[VARNAME_REGISTRY_s_G]); */
/* 	temp_I_syn_Vd = temp_s_ORN + temp_s_A + temp_s_NMDA*(1.0/temp1 + temp3*(temp_Vd-VOLTAGE_[VARNAME_REGISTRY_s_NMDA])) + temp_s_G; */
/* 	if (verbose>3){ printf(" I_syn %f\n",temp_I_syn);} */
/* 	if (verbose>3){ printf(" %% %% %% %% calculating I_leak_S and I_leak_D: ");} */
/* 	temp_I_leak_S_Vs = CONDUCTANCE_[VARNAME_REGISTRY_Vs]; */
/* 	temp_I_leak_S = temp_I_leak_S_Vs*(temp_Vs - VOLTAGE_[VARNAME_REGISTRY_s_LEAK_S]); */
/* 	temp_I_leak_D_Vd = CONDUCTANCE_[VARNAME_REGISTRY_Vd]; */
/* 	temp_I_leak_D = temp_I_leak_D_Vd*(temp_Vd - VOLTAGE_[VARNAME_REGISTRY_s_LEAK_D]); */
/* 	if (verbose>3){ printf(" I_leak_S %f, I_leak_S_Vs %f I_leak_D %f I_leak_D_Vd %f\n",temp_I_leak_S,temp_I_leak_S_Vs,temp_I_leak_D,temp_I_leak_D_Vd);} */
/* 	if (verbose>3){ printf(" %% %% %% %% calculating I_Na: ");} */
/* 	temp1 = (temp_Vs+33)/10; temp2 = -exp(-temp1); temp3 = -temp2/10; temp2+=1; */
/* 	alpha1 = temp1/temp2; alpha1p = (0.1*temp2-temp1*temp3)/pow(temp2,2); */
/* 	temp1 = (temp_Vs+53.7)/12; temp2 = 4*exp(-temp1); temp3 = -temp2/12; */
/* 	beta1 = temp2; beta1p = temp3; */
/* 	temp5 = alpha1/(alpha1+beta1); temp6 = (alpha1p*(alpha1+beta1) - alpha1*(alpha1p+beta1p))/pow(alpha1+beta1,2); */
/* 	temp1 = (temp_Vs - VOLTAGE_[VARNAME_REGISTRY_s_Na]); */
/* 	temp2 = pow(temp5,3)*temp1; */
/* 	temp3 = 3*pow(temp5,2)*temp6*temp1 + pow(temp5,3); */
/* 	temp_I_Na_h_Na = CONDUCTANCE_[VARNAME_REGISTRY_s_Na]*temp2; */
/* 	temp_I_Na = temp_I_Na_h_Na*temp_h_Na; */
/* 	temp_I_Na_Vs = CONDUCTANCE_[VARNAME_REGISTRY_s_Na]*temp_h_Na*temp3; */
/* 	temp_s_Na = CONDUCTANCE_[VARNAME_REGISTRY_s_Na]*pow(temp5,3)*temp_h_Na; */
/* 	if (verbose>3){ printf(" m_Na_inf %f, h_Na %f, I_Na %f, I_Na_Vs %f I_Na_h_Na %f\n",temp5,temp_h_Na,temp_I_Na,temp_I_Na_Vs,temp_I_Na_h_Na);} */
/* 	if (verbose>3){ printf(" %% %% %% %% calculating rhs_h_Na: ");} */
/* 	temp1 = (temp_Vs+50)/10; temp2 = 0.07*exp(-temp1); temp3 = -temp2/10; */
/* 	alpha1 = temp2; alpha1p = temp3; */
/* 	temp1 = (temp_Vs+20)/10; temp2 = exp(-temp1); temp3 = -temp2/10; temp2+=1; */
/* 	beta1 = 1/temp2; beta1p = -temp3/pow(temp2,2); */
/* 	temp_rhs_h_Na_h_Na = -4*(alpha1+beta1); */
/* 	temp_rhs_h_Na = 4*alpha1 + temp_h_Na*temp_rhs_h_Na_h_Na; */
/* 	temp_rhs_h_Na_Vs = 4*(alpha1p - temp_h_Na*(alpha1p + beta1p)); */
/* 	if (verbose>3){ printf(" rhs_h_Na %f rhs_h_Na_Vs %f rhs_h_Na_h_Na %f\n",temp_rhs_h_Na,temp_rhs_h_Na_Vs,temp_rhs_h_Na_h_Na);}   */
/* 	if (verbose>3){ printf(" %% %% %% %% calculating I_K: ");} */
/* 	temp1 = (temp_Vs - VOLTAGE_[VARNAME_REGISTRY_s_K]); */
/* 	temp_I_K_m_K = CONDUCTANCE_[VARNAME_REGISTRY_s_K]*4*pow(temp_m_K,3)*temp1; */
/* 	temp_I_K_Vs = CONDUCTANCE_[VARNAME_REGISTRY_s_K]*pow(temp_m_K,4); */
/* 	temp_I_K = temp_I_K_Vs*temp1; */
/* 	temp_s_K = temp_I_K_Vs; */
/* 	if (verbose>3){ printf(" I_K %f, I_K_Vs %f I_K_m_K %f\n",temp_I_K,temp_I_K_Vs,temp_I_K_m_K);} */
/* 	if (verbose>3){ printf(" %% %% %% %% calcuating rhs_m_K: ");} */
/* 	temp1 = (temp_Vs+34)/10; temp2 = -exp(-temp1); temp3 = -temp2/10; temp2+=1; */
/* 	alpha1 = temp1/temp2/10; alpha1p = (0.1*temp2-temp1*temp3)/pow(temp2,2)/10; */
/* 	temp1 = (temp_Vs+44)/25; temp2 = 0.125*exp(-temp1); temp3 = -temp2/25; */
/* 	beta1 = temp2; beta1p = temp3; */
/* 	temp_rhs_m_K_m_K = -4*(alpha1+beta1); */
/* 	temp_rhs_m_K = 4*alpha1 + temp_m_K*temp_rhs_m_K_m_K; */
/* 	temp_rhs_m_K_Vs = 4*(alpha1p - temp_m_K*(alpha1p + beta1p)); */
/* 	if (verbose>3){ printf(" rhs_m_K %f rhs_m_K_Vs %f rhs_m_K_m_K %f\n",temp_rhs_m_K,temp_rhs_m_K_Vs,temp_rhs_m_K_m_K);} */
/* 	if (verbose>3){ printf(" %% %% %% %% calculating I_KA: ");} */
/* 	temp1 = (temp_Vs+50)/20; temp2 = exp(-temp1); temp3 = -temp2/20; temp2+=1; */
/* 	temp4 = 1/temp2; temp5 = -temp3/pow(temp2,2); */
/* 	temp1 = (temp_Vs - VOLTAGE_[VARNAME_REGISTRY_s_KA]); */
/* 	temp_I_KA_h_KA = CONDUCTANCE_[VARNAME_REGISTRY_s_KA]*pow(temp4,3)*temp1; */
/* 	temp_I_KA_Vs = CONDUCTANCE_[VARNAME_REGISTRY_s_KA]*temp_h_KA*(3*pow(temp4,2)*temp5*temp1+pow(temp4,3)); */
/* 	temp_I_KA = temp_I_KA_h_KA*temp_h_KA; */
/* 	temp_s_KA = temp_I_KA/temp1;; */
/* 	if (verbose>3){ printf(" I_KA %f, I_KA_Vs %f I_KA_h_KA %f\n",temp_I_KA,temp_I_KA_Vs,temp_I_KA_h_KA);} */
/* 	if (verbose>3){ printf(" %% %% %% %% calcuating rhs_h_KA: ");} */
/* 	temp1 = (temp_Vs+80)/6; temp2 = exp(temp1); temp3 = temp2/6; temp2+=1; */
/* 	alpha1 = 1/temp2; alpha1p = -temp3/pow(temp2,2); */
/* 	temp_rhs_h_KA = (alpha1-temp_h_KA)/TAU_[VARNAME_REGISTRY_h_KA]; */
/* 	temp_rhs_h_KA_h_KA = -1/TAU_[VARNAME_REGISTRY_h_KA]; */
/* 	temp_rhs_h_KA_Vs = alpha1p/TAU_[VARNAME_REGISTRY_h_KA]; */
/* 	if (verbose>3){ printf(" rhs_h_KA %f rhs_h_KA_Vs %f rhs_h_KA_h_KA %f\n",temp_rhs_h_KA,temp_rhs_h_KA_Vs,temp_rhs_h_KA_h_KA);} */
/* 	if (verbose>3){ printf(" %% %% %% %% calculating I_KS: ");} */
/* 	temp1 = (temp_Vs - VOLTAGE_[VARNAME_REGISTRY_s_KS]); */
/* 	temp_I_KS_m_KS = CONDUCTANCE_[VARNAME_REGISTRY_s_KS]*temp1; */
/* 	temp_I_KS = temp_m_KS*temp_I_KS_m_KS; */
/* 	temp_I_KS_Vs = CONDUCTANCE_[VARNAME_REGISTRY_s_KS]*temp_m_KS; */
/* 	if (verbose>3){ printf(" I_KS_m_KS %f I_KS %f I_KS_Vs %f\n",temp_I_KS_m_KS,temp_I_KS,temp_I_KS_Vs);} */
/* 	if (verbose>3){ printf(" %% %% %% %% calculating rhs_m_KS: ");} */
/* 	temp1 = (temp_Vs+34)/6.5; temp2 = exp(-temp1); temp3 = -temp2/6.5; temp2+=1; */
/* 	alpha1 = 1/temp2; alpha1p = -temp3/pow(temp2,2); */
/* 	temp1 = (temp_Vs+55)/30; temp2 = 4/cosh(temp1); temp3 = -temp2*tanh(temp1)/30; */
/* 	beta1 = temp2; beta1p = temp3;  */
/* 	temp_rhs_m_KS = (alpha1 - temp_m_KS)/beta1; */
/* 	temp_rhs_m_KS_m_KS = -1/beta1; */
/* 	temp_rhs_m_KS_Vs = (alpha1p*beta1 - (alpha1-temp_m_KS)*beta1p)/pow(beta1,2); */
/* 	if (verbose>3){ printf(" rhs_m_KS %f rhs_m_KS_m_KS %f rhs_m_KS_Vs %f\n",temp_rhs_m_KS,temp_rhs_m_KS_m_KS,temp_rhs_m_KS_Vs);} */
/* 	if (verbose>3){ printf(" %% %% %% %% calculating I_KNa: ");} */
/* 	temp1 = (temp_Vs - VOLTAGE_[VARNAME_REGISTRY_s_KNa]); */
/* 	temp2 = pow(maximum(0,temp_Na),3.5); temp3 = temp2 + 360570; temp4 = 3.5*pow(maximum(0,temp_Na),2.5); */
/* 	alpha1 = 0.37*temp2/temp3; alpha1p = 0.37*(temp4*temp3-temp2*temp4)/pow(temp3,2); */
/* 	temp_I_KNa_Vs = CONDUCTANCE_[VARNAME_REGISTRY_s_KNa]*alpha1; */
/* 	temp_I_KNa_Na = CONDUCTANCE_[VARNAME_REGISTRY_s_KNa]*alpha1p*temp1; */
/* 	temp_I_KNa = CONDUCTANCE_[VARNAME_REGISTRY_s_KNa]*alpha1p*temp1; */
/* 	if (verbose>3){ printf(" I_KNa_Na %f I_KNa %f I_KNa_Vs %f\n",temp_I_KNa_Na,temp_I_KNa,temp_I_KNa_Vs);} */
/* 	if (verbose>3){ printf(" %% %% %% %% calculating I_SD and I_DS: ");} */
/* 	temp1 = (temp_Vs - temp_Vd); */
/* 	temp_I_SD = CONDUCTANCE_SD*temp1; */
/* 	temp_I_SD_Vs = CONDUCTANCE_SD; */
/* 	temp_I_SD_Vd = -CONDUCTANCE_SD; */
/* 	temp_I_DS = -CONDUCTANCE_DS*temp1; */
/* 	temp_I_DS_Vs = -CONDUCTANCE_DS; */
/* 	temp_I_DS_Vd = CONDUCTANCE_DS; */
/* 	if (verbose>3){ printf(" I_SD %f I_SD_Vs %f I_SD_Vd %f I_DS %f I_DS_Vs %f I_DS_Vd %f\n",temp_I_SD,temp_I_SD_Vs,temp_I_SD_Vd,temp_I_DS,temp_I_DS_Vs,temp_I_DS_Vd);}   */
/* 	if (verbose>3){ printf(" %% %% %% %% calculating I_Ca: ");} */
/* 	temp1 = (temp_Vd+20)/9; temp2 = exp(-temp1); temp3 = -temp2/9; temp2+=1; */
/* 	alpha1 = 1/temp2; alpha1p = -temp3/pow(temp2,2); */
/* 	temp1 = (temp_Vd - VOLTAGE_[VARNAME_REGISTRY_s_Ca]); */
/* 	temp_I_Ca = CONDUCTANCE_[VARNAME_REGISTRY_s_Ca]*pow(alpha1,2)*temp1; */
/* 	temp_I_Ca_Vd = CONDUCTANCE_[VARNAME_REGISTRY_s_Ca]*(2*alpha1*alpha1p*temp1 + pow(alpha1,2)); */
/* 	temp_s_Ca = temp_I_Ca/temp1; */
/* 	if (verbose>3){ printf(" I_Ca %f, I_Ca_Vd %f\n",temp_I_Ca,temp_I_Ca_Vd);} */
/* 	if (verbose>3){ printf(" %% %% %% %% calculating I_NaP: ");} */
/* 	temp1 = (temp_Vd+55.7)/7.7; temp2 = exp(-temp1); temp3 = -temp2/7.7; temp2+=1; */
/* 	alpha1 = 1/temp2; alpha1p = -temp3/pow(temp2,2); */
/* 	temp1 = (temp_Vd - VOLTAGE_[VARNAME_REGISTRY_s_NaP]); */
/* 	temp_I_NaP = CONDUCTANCE_[VARNAME_REGISTRY_s_NaP]*pow(alpha1,3)*temp1; */
/* 	temp_I_NaP_Vd = CONDUCTANCE_[VARNAME_REGISTRY_s_NaP]*(3*pow(alpha1,2)*alpha1p*temp1 + pow(alpha1,3)); */
/* 	temp_s_NaP = temp_I_NaP/temp1; */
/* 	if (verbose>3){ printf(" I_NaP %f, I_NaP_Vd %f\n",temp_I_NaP,temp_I_NaP_Vd);} */
/* 	if (verbose>3){ printf(" %% %% %% %% calculating I_AR: ");} */
/* 	temp1 = (temp_Vd+75)/4; temp2 = exp(temp1); temp3 = temp2/4; temp2+=1; */
/* 	alpha1 = 1/temp2; alpha1p = -temp3/pow(temp2,2); */
/* 	temp1 = (temp_Vd - VOLTAGE_[VARNAME_REGISTRY_s_AR]); */
/* 	temp_I_AR = CONDUCTANCE_[VARNAME_REGISTRY_s_AR]*alpha1*temp1; */
/* 	temp_I_AR_Vd = CONDUCTANCE_[VARNAME_REGISTRY_s_AR]*(alpha1p*temp1 + alpha1); */
/* 	temp_s_AR = temp_I_AR/temp1; */
/* 	if (verbose>3){ printf(" I_AR %f, I_AR_Vd %f\n",temp_I_AR,temp_I_AR_Vd);} */
/* 	if (verbose>3){ printf(" %% %% %% %% calculating I_KCa: ");} */
/* 	temp1 = (temp_Vd - VOLTAGE_[VARNAME_REGISTRY_s_KCa]); */
/* 	temp2 = temp_Ca/(temp_Ca+30); temp3 = 30/pow(temp_Ca+30,2); */
/* 	temp_I_KCa = CONDUCTANCE_[VARNAME_REGISTRY_s_KCa]*temp2*temp1; */
/* 	temp_I_KCa_Vd = CONDUCTANCE_[VARNAME_REGISTRY_s_KCa]*temp2; */
/* 	temp_I_KCa_Ca = CONDUCTANCE_[VARNAME_REGISTRY_s_KCa]*temp3*temp1; */
/* 	temp_s_KCa = temp_I_KCa/temp1; */
/* 	if (verbose>3){ printf(" I_KCa %f, I_KCa_Vd %f, I_KCa_Ca %f\n",temp_I_KCa,temp_I_KCa_Vd,temp_I_KCa_Ca);} */
/* 	if (verbose>3){ printf(" %% %% %% %% calculating rhs_Ca: ");} */
/* 	temp_rhs_Ca = -0.00175*temp_I_Ca - temp_Ca/TAU_[VARNAME_REGISTRY_Ca]; */
/* 	temp_rhs_Ca_Ca = -1/TAU_[VARNAME_REGISTRY_Ca]; */
/* 	temp_rhs_Ca_Vd = -0.00175*temp_I_Ca_Vd; */
/* 	if (verbose>3){ printf(" rhs_Ca %f, rhs_Ca_Vd %f, rhs_Ca_Ca %f\n",temp_rhs_Ca,temp_rhs_Ca_Vd,temp_rhs_Ca_Ca);} */
/* 	if (verbose>3){ printf(" %% %% %% %% calculating rhs_Na: ");} */
/* 	temp1 = pow(temp_Na,3); temp2 = temp1+pow(15,3); temp3 = 3*pow(temp_Na,2); */
/* 	alpha1 = temp1/temp2; alpha1p = (temp3*temp2-temp1*temp3)/pow(temp2,2); */
/* 	temp_rhs_Na = -0.0015*temp_I_Na + -0.0035*temp_I_NaP - (alpha1 - 0.2026)/TAU_[VARNAME_REGISTRY_Na]; */
/* 	temp_rhs_Na_Na = -alpha1p/TAU_[VARNAME_REGISTRY_Na]; */
/* 	temp_rhs_Na_Vs = -0.0015*temp_I_Na_Vs; */
/* 	temp_rhs_Na_h_Na = -0.0015*temp_I_Na_h_Na; */
/* 	temp_rhs_Na_Vd = -0.0035*temp_I_NaP_Vd; */
/* 	if (verbose>3){ printf(" rhs_Na %f, rhs_Na_Vd %f, rhs_Na_Na %f, rhs_Na_h_Na %f rhs_Na_Vs %f\n",temp_rhs_Na,temp_rhs_Na_Vd,temp_rhs_Na_Na,temp_rhs_Na_h_Na,temp_rhs_Na_Vs);} */
/* 	if (verbose>2){ printf(" %% %% %% %% calculating: implicit euler update\n");} */
/* 	temp_rhs_Vs = -(temp_I_leak_S+temp_I_Na+temp_I_K+temp_I_KA+temp_I_KS+temp_I_KNa+temp_I_SD+CURRENT_INJECTION_S)/TAU_[VARNAME_REGISTRY_Vs]; */
/* 	temp_rhs_Vs_Vs = -(temp_I_leak_S_Vs+temp_I_Na_Vs+temp_I_K_Vs+temp_I_KA_Vs+temp_I_KS_Vs+temp_I_KNa_Vs+temp_I_SD_Vs)/TAU_[VARNAME_REGISTRY_Vs]; */
/* 	temp_rhs_Vs_h_Na = -(temp_I_Na_h_Na)/TAU_[VARNAME_REGISTRY_Vs]; */
/* 	temp_rhs_Vs_m_K = -(temp_I_K_m_K)/TAU_[VARNAME_REGISTRY_Vs]; */
/* 	temp_rhs_Vs_h_KA = -(temp_I_KA_h_KA)/TAU_[VARNAME_REGISTRY_Vs]; */
/* 	temp_rhs_Vs_m_KS = -(temp_I_KS_m_KS)/TAU_[VARNAME_REGISTRY_Vs]; */
/* 	temp_rhs_Vs_Na = -(temp_I_KNa_Na)/TAU_[VARNAME_REGISTRY_Vs]; */
/* 	temp_rhs_Vs_Vd = -(temp_I_SD_Vd)/TAU_[VARNAME_REGISTRY_Vs]; */
/* 	temp_rhs_Vd = -(temp_I_leak_D+temp_I_Ca+temp_I_NaP+temp_I_AR+temp_I_KCa+temp_I_DS+temp_I_syn+CURRENT_INJECTION_D)/TAU_[VARNAME_REGISTRY_Vd]; */
/* 	temp_rhs_Vd_Vd = -(temp_I_leak_D_Vd+temp_I_Ca_Vd+temp_I_NaP_Vd+temp_I_AR_Vd+temp_I_KCa_Vd+temp_I_DS_Vd+temp_I_syn_Vd)/TAU_[VARNAME_REGISTRY_Vd]; */
/* 	temp_rhs_Vd_Vs = -(temp_I_DS_Vs)/TAU_[VARNAME_REGISTRY_Vd]; */
/* 	temp_rhs_Vd_Ca = -(temp_I_KCa_Ca)/TAU_[VARNAME_REGISTRY_Vd]; */
/* 	specialized_inverse_compte(1-dt*temp_rhs_Vs_Vs,-dt*temp_rhs_Vs_h_Na,-dt*temp_rhs_Vs_m_K,-dt*temp_rhs_Vs_h_KA,-dt*temp_rhs_Vs_m_KS,-dt*temp_rhs_Vs_Na,-dt*temp_rhs_Vs_Vd,-dt*temp_rhs_h_Na_Vs,1-dt*temp_rhs_h_Na_h_Na,-dt*temp_rhs_m_K_Vs,1-dt*temp_rhs_m_K_m_K,-dt*temp_rhs_h_KA_Vs,1-dt*temp_rhs_h_KA_h_KA,-dt*temp_rhs_m_KS_Vs,1-dt*temp_rhs_m_KS_m_KS,-dt*temp_rhs_Na_Vs,-dt*temp_rhs_Na_h_Na,1-dt*temp_rhs_Na_Na,-dt*temp_rhs_Na_Vd,-dt*temp_rhs_Vd_Vs,1-dt*temp_rhs_Vd_Vd,-dt*temp_rhs_Vd_Ca,-dt*temp_rhs_Ca_Vd,1-dt*temp_rhs_Ca_Ca,temp_Vs_start - temp_Vs + dt*temp_rhs_Vs,temp_h_Na_start - temp_h_Na + dt*temp_rhs_h_Na,temp_m_K_start - temp_m_K + dt*temp_rhs_m_K,temp_h_KA_start - temp_h_KA + dt*temp_rhs_h_KA,temp_m_KS_start - temp_m_KS + dt*temp_rhs_m_KS,temp_Na_start - temp_Na + dt*temp_rhs_Na,temp_Vd_start - temp_Vd + dt*temp_rhs_Vd,temp_Ca_start - temp_Ca + dt*temp_rhs_Ca,&update_Vs,&update_h_Na,&update_m_K,&update_h_KA,&update_m_KS,&update_Na,&update_Vd,&update_Ca); */
/* 	norm_old = pow(update_Vs,2)+pow(update_h_Na,2)+pow(update_m_K,2)+pow(update_h_KA,2)+pow(update_m_KS,2)+pow(update_Na,2)+pow(update_Vd,2)+pow(update_Ca,2); */
/* 	temp_Vs += update_Vs; */
/* 	temp_h_Na += update_h_Na; */
/* 	temp_m_K += update_m_K; */
/* 	temp_h_KA += update_h_KA; */
/* 	temp_m_KS += update_m_KS; */
/* 	temp_Na += update_Na; */
/* 	temp_Vd += update_Vd; */
/* 	temp_Ca += update_Ca; */
/* 	temp7=norm_old;norm_old=norm_new;norm_new=temp7; */
/* 	iteration += 1; */
/* 	if (verbose>2){ printf(" %% %% %% norm_new %0.16f norm_old %0.16f\n",norm_new,norm_old);}} */
/*       while (norm_new>norm_threshold && norm_new<norm_old && iteration<iteration_max); */
/*       if (verbose>1){ printf(" %% %% iteration=%d,norm_new %0.16f\n",iteration,norm_new);} */
/*       temp_s_A_start = temp_s_A; */
/*       temp_s_N_start = temp_s_N; */
/*       temp_s_NMDA_start = temp_s_NMDA; */
/*       temp_s_G_start = temp_s_G; */
/*       temp_s_ORN_start = temp_s_ORN; */
/*       temp_h_Na_start = maximum(0,minimum(1,temp_h_Na)); */
/*       temp_m_K_start = maximum(0,minimum(1,temp_m_K)); */
/*       temp_h_KA_start = maximum(0,minimum(1,temp_h_KA)); */
/*       temp_m_KS_start = maximum(0,minimum(1,temp_m_KS)); */
/*       temp_Na_start = temp_Na; */
/*       temp_Ca_start = temp_Ca; */
/*       temp_Vs_start = temp_Vs; */
/*       temp_Vd_start = temp_Vd; */
/*       microstep += 1; t+=dt; */
/*       if (norm_new>norm_old || iteration==iteration_max || norm_new>norm_threshold){  */
/* 	if (verbose>2){ printf(" %% %% convergence poor, setting error_flag=1\n");}  */
/* 	error_flag=1;}  */
/*       else{  */
/* 	if (verbose>2){ printf(" %% %% convergence good, setting error_flag=0\n");}  */
/* 	error_flag=0;}} */
/*     while (microstep<n->microstep && !error_flag);       */
/*     if (error_flag){ n->microstep*=2;} else{ n->microstep=maximum(1,n->microstep/2);} */
/*     if (verbose>2){ printf(" %% total dt %f error_flag=%d, n->microstep set to %d\n",t-t_start,error_flag,n->microstep);}} */
/*   while(error_flag); */
/*   temp_Vs_end = temp_Vs_start; temp_Vs_start = n->vra[VARNAME_REGISTRY_Vs]; */
/*   temp_Vd_end = temp_Vd_start; temp_Vd_start = n->vra[VARNAME_REGISTRY_Vd]; */
/*   if (verbose>2){ printf(" %% determining spiking properties\n");} */
/*   if (n->spikenext<t_start+dt_start){ */
/*     temp1=-1; */
/*     if (temp_Vs_end>=VOLTAGE_THRESHOLD_S && temp_Vs_start<VOLTAGE_THRESHOLD_S){  */
/*       temp1 = linerootfinder(temp_Vs_start,temp_Vs_end,VOLTAGE_THRESHOLD_S,dt_start);} */
/*     if (temp_Vd_end>=VOLTAGE_THRESHOLD_D && temp_Vd_start<VOLTAGE_THRESHOLD_D){  */
/*       temp1 = linerootfinder(temp_Vd_start,temp_Vd_end,VOLTAGE_THRESHOLD_D,dt_start);} */
/*     if (temp1>=0){ */
/*       n->spiketime_guess = t_start+temp1; */
/*       if (n->spiketime_guess<=n->spikenext){ */
/* 	n->spiketime_guess=-2*dt_start; n->spiketime_guess_flag=0;} */
/*       else /\* if (n->spiketime_guess>n->spikenext) *\/{ 	 */
/* 	n->vra[VARNAME_REGISTRY_spike_flag] += 1; */
/* 	n->spiketime_guess_flag=1;}}} */
/*   n->vra[VARNAME_REGISTRY_spike_flag] *= exp(-dt_start/(TAU_REF/16)); */
/*   n->vra[VARNAME_REGISTRY_s_A] = temp_s_A_start; */
/*   n->vra[VARNAME_REGISTRY_s_N] = temp_s_N_start; */
/*   n->vra[VARNAME_REGISTRY_s_NMDA] = temp_s_NMDA_start; */
/*   n->vra[VARNAME_REGISTRY_s_G] = temp_s_G_start; */
/*   n->vra[VARNAME_REGISTRY_s_ORN] = temp_s_ORN_start; */
/*   n->vra[VARNAME_REGISTRY_h_Na] = temp_h_Na_start; */
/*   n->vra[VARNAME_REGISTRY_m_K] = temp_m_K_start; */
/*   n->vra[VARNAME_REGISTRY_h_KA] = temp_h_KA_start; */
/*   n->vra[VARNAME_REGISTRY_m_KS] = temp_m_KS_start; */
/*   n->vra[VARNAME_REGISTRY_Na] = temp_Na_start; */
/*   n->vra[VARNAME_REGISTRY_Ca] = temp_Ca_start; */
/*   n->vra[VARNAME_REGISTRY_Vs] = temp_Vs_end; */
/*   n->vra[VARNAME_REGISTRY_Vd] = temp_Vd_end;   */
/*   n->vra[VARNAME_REGISTRY_s_Na] = temp_s_Na; */
/*   n->vra[VARNAME_REGISTRY_s_K] = temp_s_K; */
/*   n->vra[VARNAME_REGISTRY_s_KA] = temp_s_KA; */
/*   n->vra[VARNAME_REGISTRY_s_KS] = temp_s_KS; */
/*   n->vra[VARNAME_REGISTRY_s_KNa] = temp_s_KNa; */
/*   n->vra[VARNAME_REGISTRY_s_Ca] = temp_s_Ca; */
/*   n->vra[VARNAME_REGISTRY_s_NaP] = temp_s_NaP; */
/*   n->vra[VARNAME_REGISTRY_s_AR] = temp_s_AR; */
/*   n->vra[VARNAME_REGISTRY_s_KCa] = temp_s_KCa; */
/*   if (verbose>1){ printf(" %% exiting vraevolve: \n"); raprintf(n->vra,"double",1,n->nvars," end ");} */
/*   if (verbose>0){ printf(" %% vra:\n"); for (nv=0;nv<n->nvars;nv++){ printf("%s %f \t",GLOBAL_VARNAMES[nv],n->vra[nv]);} printf("\n");} */
/* }   */

/* pinsky functions */

/* void specialized_inverse_pinsky(double M00,double M01,double M02,double M03,double M10,double M11,double M20,double M22,double M30,double M33,double M34,double M35,double M36,double M37,double M43,double M44,double M53,double M55,double M63,double M64,double M66,double M76,double M77,double R0,double R1,double R2,double R3,double R4,double R5,double R6,double R7,double *V0,double *V1,double *V2,double *V3,double *V4,double *V5,double *V6,double *V7) */
/* { */
/*   /\* specialized 69 flop inverse (row operations) for the matrix */
/*     X  0    1    2    3    4    5    6    7 */
/*     0  M00  M01  M02  M03   */
/*     1  M10  M11 */
/*     2  M20       M22 */
/*     3  M30            M33  M34  M35  M36  M37 */
/*     4                 M43  M44 */
/*     5                 M53       M55 */
/*     6                 M63  M64       M66 */
/*     7                                M76  M77 */
/*     with RHS [R1 R2 R3 R4 R5 R6 R7] */
/*   *\/ */
/*   double temp=0; */
/*   temp = M37/M77; M36 -= temp*M76; R3 -= temp*R7; /\* M37=0; *\/ */
/*   temp = M36/M66; M34 -= temp*M64; M33 -= temp*M63; R3 -=temp*R6; /\* M36=0; *\/ */
/*   temp = M35/M55; M33 -= temp*M53; R3 -= temp*R5; /\* M35=0; *\/ */
/*   temp = M34/M44; M33 -= temp*M43; R3 -= temp*R4; /\* M34=0; *\/ */
/*   temp = M03/M33; M00 -= temp*M30; R0 -= temp*R3; /\* M03=0; *\/ */
/*   temp = M02/M22; M00 -= temp*M20; R0 -= temp*R2; /\* M02=0; *\/ */
/*   temp = M01/M11; M00 -= temp*M10; R0 -= temp*R1; /\* M01=0; *\/ */
/*   /\* now matrix should look like  */
/*     X  0    1    2    3    4    5    6    7 */
/*     0  M00   */
/*     1  M10  M11 */
/*     2  M20       M22 */
/*     3  M30            M33   */
/*     4                 M43  M44 */
/*     5                 M53       M55 */
/*     6                 M63  M64       M66 */
/*     7                                M76  M77 */
/*   *\/ */
/*   temp = M10/M00; R1 -= temp*R0; /\* M10=0; *\/ */
/*   temp = M20/M00; R2 -= temp*R0; /\* M20=0; *\/ */
/*   temp = M30/M00; R3 -= temp*R0; /\* M30=0; *\/ */
/*   temp = M43/M33; R4 -= temp*R3; /\* M43=0; *\/ */
/*   temp = M53/M33; R5 -= temp*R3; /\* M53=0; *\/ */
/*   temp = M63/M33; R6 -= temp*R3; /\* M63=0; *\/ */
/*   temp = M64/M44; R6 -= temp*R4; /\* M64=0; *\/ */
/*   temp = M76/M66; R7 -= temp*R6; /\* M76=0; *\/ */
/*   /\* now matrix should be diagonal *\/ */
/*   *V0 = R0/M00; */
/*   *V1 = R1/M11; */
/*   *V2 = R2/M22; */
/*   *V3 = R3/M33; */
/*   *V4 = R4/M44; */
/*   *V5 = R5/M55; */
/*   *V6 = R6/M66; */
/*   *V7 = R7/M77; */
/* } */

/* void vraevolve_pinsky(struct neuron *n,double t_start,double dt_start) */
/* { */
/*   /\* linearly implicit euler *\/ */
/*   int verbose=0; */
/*   double temp1=0,temp2=0,temp3=0,temp4=0,temp5=0,temp6=0,temp7=0,temp8=0; */
/*   double alpha1=0,alpha1p=0,beta1=0,beta1p=0; */
/*   double alpha2=0,alpha2p=0,beta2=0,beta2p=0; */
/*   double alpha3=0,alpha3p=0,beta3=0,beta3p=0; */
/*   double temp_s_A=0;  */
/*   double temp_s_N=0;  */
/*   double temp_s_NMDA=0;  */
/*   double temp_s_G=0;  */
/*   double temp_s_ORN=0;  */
/*   double temp_h_Na=0;  */
/*   double temp_m_K=0;  */
/*   double temp_m_KAHP=0;  */
/*   double temp_m_Ca=0;  */
/*   double temp_m_KCa=0;  */
/*   double temp_Ca=0;  */
/*   double temp_Vs=0;  */
/*   double temp_Vd=0; */
/*   double temp_s_Na=0; */
/*   double temp_s_K=0; */
/*   double temp_s_Ca=0; */
/*   double temp_s_KCa=0; */
/*   double temp_s_KAHP=0; */
/*   double temp_I_syn=0; */
/*   double temp_I_syn_Vd=0; */
/*   double temp_I_leak_S=0; */
/*   double temp_I_leak_S_Vs=0; */
/*   double temp_I_leak_D=0; */
/*   double temp_I_leak_D_Vd=0; */
/*   double temp_I_Na=0; */
/*   double temp_I_Na_Vs=0; */
/*   double temp_I_Na_h_Na=0; */
/*   double temp_rhs_h_Na=0; */
/*   double temp_rhs_h_Na_Vs=0; */
/*   double temp_rhs_h_Na_h_Na=0; */
/*   double temp_I_K=0; */
/*   double temp_I_K_Vs=0; */
/*   double temp_I_K_m_K=0; */
/*   double temp_rhs_m_K=0; */
/*   double temp_rhs_m_K_Vs=0; */
/*   double temp_rhs_m_K_m_K=0; */
/*   double temp_I_SD=0; */
/*   double temp_I_SD_Vs=0; */
/*   double temp_I_SD_Vd=0; */
/*   double temp_I_DS=0; */
/*   double temp_I_DS_Vs=0; */
/*   double temp_I_DS_Vd=0; */
/*   double temp_I_Ca=0; */
/*   double temp_I_Ca_Vd=0; */
/*   double temp_I_Ca_m_Ca=0; */
/*   double temp_rhs_m_Ca=0; */
/*   double temp_rhs_m_Ca_Vd=0; */
/*   double temp_rhs_m_Ca_m_Ca=0; */
/*   double temp_rhs_Ca=0; */
/*   double temp_rhs_Ca_Vd=0; */
/*   double temp_rhs_Ca_m_Ca=0; */
/*   double temp_rhs_Ca_Ca=0; */
/*   double temp_I_KCa=0; */
/*   double temp_I_KCa_Vd=0; */
/*   double temp_I_KCa_m_KCa=0; */
/*   double temp_I_KCa_Ca=0; */
/*   double temp_rhs_m_KCa=0; */
/*   double temp_rhs_m_KCa_Vd=0; */
/*   double temp_rhs_m_KCa_m_KCa=0; */
/*   double temp_I_KAHP=0; */
/*   double temp_I_KAHP_Vd=0; */
/*   double temp_I_KAHP_m_KAHP=0; */
/*   double temp_rhs_m_KAHP=0; */
/*   double temp_rhs_m_KAHP_Ca=0; */
/*   double temp_rhs_m_KAHP_m_KAHP=0; */
/*   double temp_rhs_Vs=0; */
/*   double temp_rhs_Vs_Vs=0; */
/*   double temp_rhs_Vs_h_Na=0; */
/*   double temp_rhs_Vs_m_K=0; */
/*   double temp_rhs_Vs_Vd=0; */
/*   double temp_rhs_Vd=0; */
/*   double temp_rhs_Vd_Vs=0; */
/*   double temp_rhs_Vd_Vd=0; */
/*   double temp_rhs_Vd_m_Ca=0; */
/*   double temp_rhs_Vd_m_KCa=0; */
/*   double temp_rhs_Vd_Ca=0; */
/*   double temp_rhs_Vd_m_KAHP=0; */
/*   double temp_s_A_start=0; */
/*   double temp_s_N_start=0; */
/*   double temp_s_NMDA_start=0; */
/*   double temp_s_G_start=0; */
/*   double temp_s_ORN_start=0; */
/*   double temp_Vs_start=0; */
/*   double temp_h_Na_start=0; */
/*   double temp_m_K_start=0; */
/*   double temp_Vd_start=0; */
/*   double temp_m_Ca_start=0; */
/*   double temp_m_KCa_start=0; */
/*   double temp_Ca_start=0; */
/*   double temp_m_KAHP_start=0; */
/*   double temp_Vs_end=0; */
/*   double temp_Vd_end=0; */
/*   double dt=dt_start,t=t_start; */
/*   int microstep=0,error_flag=0; */
/*   int iteration=0,iteration_max=4; */
/*   double update_Vs=0; */
/*   double update_h_Na=0; */
/*   double update_m_K=0; */
/*   double update_Vd=0; */
/*   double update_m_Ca=0; */
/*   double update_m_KCa=0; */
/*   double update_Ca=0; */
/*   double update_m_KAHP=0; */
/*   double norm_new=0,norm_old=0,norm_threshold=0.0000001; */
/*   int nv=0; */
/*   if (verbose>1){ printf(" %% [entering vraevolve_pinsky]: \n"); raprintf(n->vra,"double",1,n->nvars," start ");} */
/*   do{ /\* outer loop, choose n->microstep *\/ */
/*     error_flag=0;dt=dt_start/n->microstep;t=t_start; */
/*     if (verbose>1){ printf(" %% initiating pass -- n->microstep %d dt %f\n",n->microstep,dt);} */
/*     temp_s_A_start = n->vra[VARNAME_REGISTRY_s_A]; */
/*     temp_s_N_start = n->vra[VARNAME_REGISTRY_s_N]; */
/*     temp_s_NMDA_start = n->vra[VARNAME_REGISTRY_s_NMDA]; */
/*     temp_s_G_start = n->vra[VARNAME_REGISTRY_s_G]; */
/*     temp_s_ORN_start = n->vra[VARNAME_REGISTRY_s_ORN]; */
/*     temp_Vs_start=n->vra[VARNAME_REGISTRY_Vs]; */
/*     temp_h_Na_start=n->vra[VARNAME_REGISTRY_h_Na]; */
/*     temp_m_K_start=n->vra[VARNAME_REGISTRY_m_K]; */
/*     temp_Vd_start=n->vra[VARNAME_REGISTRY_Vd]; */
/*     temp_m_Ca_start=n->vra[VARNAME_REGISTRY_m_Ca]; */
/*     temp_m_KCa_start=n->vra[VARNAME_REGISTRY_m_KCa]; */
/*     temp_Ca_start=n->vra[VARNAME_REGISTRY_Ca]; */
/*     temp_m_KAHP_start=n->vra[VARNAME_REGISTRY_m_KAHP]; */
/*     microstep=0; */
/*     do{ /\* each microstep *\/ */
/*       if (verbose>1){ printf(" %% %% microstep %d\n",microstep);} */
/*       norm_old=16.0;norm_new=16.0; */
/*       temp_s_A = temp_s_A_start * exp(-dt/TAU_[VARNAME_REGISTRY_s_A]); */
/*       temp_s_N = temp_s_N_start * exp(-dt/TAU_[VARNAME_REGISTRY_s_N]); */
/*       temp_s_NMDA = exp(-dt/TAU_[VARNAME_REGISTRY_s_NMDA])*(temp_s_NMDA_start + temp_s_N_start*(-1 + exp(dt/(1.0/TAU_[VARNAME_REGISTRY_s_NMDA]-1.0/TAU_[VARNAME_REGISTRY_s_N])))/(1.0/TAU_[VARNAME_REGISTRY_s_NMDA]-1.0/TAU_[VARNAME_REGISTRY_s_N])); */
/*       temp_s_G = temp_s_G_start * exp(-dt/TAU_[VARNAME_REGISTRY_s_G]); */
/*       temp_s_ORN = temp_s_ORN_start * exp(-dt/TAU_[VARNAME_REGISTRY_s_ORN]); */
/*       temp_h_Na = temp_h_Na_start; */
/*       temp_m_K = temp_m_K_start; */
/*       temp_m_KAHP = temp_m_KAHP_start; */
/*       temp_m_Ca = temp_m_Ca_start; */
/*       temp_m_KCa = temp_m_KCa_start; */
/*       temp_Ca = temp_Ca_start; */
/*       temp_Vs = temp_Vs_start; */
/*       temp_Vd = temp_Vd_start; */
/*       iteration=0; */
/*       do{ /\* inner loop, linearly implicit euler *\/ */
/* 	if (verbose>2){ printf(" %% %% %% iteration %d\n",iteration);} */
/* 	if (verbose>3){ printf(" %% %% %% %% calculating I_syn: ");} */
/* 	temp1 = exp(-0.062*(temp_Vd-60)); temp2 = temp1*-0.062; temp1+=1; temp3=-temp2/pow(temp1,2); */
/* 	temp_I_syn = temp_s_ORN*(temp_Vd - VOLTAGE_[VARNAME_REGISTRY_s_ORN]) + temp_s_A*(temp_Vd - VOLTAGE_[VARNAME_REGISTRY_s_A]) + temp_s_NMDA/temp1*(temp_Vd - VOLTAGE_[VARNAME_REGISTRY_s_NMDA]) + temp_s_G*(temp_Vd - VOLTAGE_[VARNAME_REGISTRY_s_G]); */
/* 	temp_I_syn_Vd = temp_s_ORN + temp_s_A + temp_s_NMDA*(1.0/temp1 + temp3*(temp_Vd-VOLTAGE_[VARNAME_REGISTRY_s_NMDA])) + temp_s_G; */
/* 	if (verbose>3){ printf(" I_syn %f\n",temp_I_syn);} */
/* 	if (verbose>3){ printf(" %% %% %% %% calculating I_leak_S and I_leak_D: ");} */
/* 	temp_I_leak_S_Vs = CONDUCTANCE_[VARNAME_REGISTRY_Vs]; */
/* 	temp_I_leak_S = temp_I_leak_S_Vs*(temp_Vs - VOLTAGE_[VARNAME_REGISTRY_s_LEAK_S]); */
/* 	temp_I_leak_D_Vd = CONDUCTANCE_[VARNAME_REGISTRY_Vd]; */
/* 	temp_I_leak_D = temp_I_leak_D_Vd*(temp_Vd - VOLTAGE_[VARNAME_REGISTRY_s_LEAK_D]); */
/* 	if (verbose>3){ printf(" I_leak_S %f, I_leak_S_Vs %f I_leak_D %f I_leak_D_Vd %f\n",temp_I_leak_S,temp_I_leak_S_Vs,temp_I_leak_D,temp_I_leak_D_Vd);} */
/* 	if (verbose>3){ printf(" %% %% %% %% calculating I_Na: ");} */
/* 	temp1 = 13.1-temp_Vs; temp2 = exp(temp1/4); temp3 = -temp2/4; temp2-=1; */
/* 	alpha1 = 0.32*temp1/temp2; alpha1p = 0.32*(-temp2-temp1*temp3)/pow(temp2,2); */
/* 	temp1 = temp_Vs-40.1; temp2 = exp(temp1/5); temp3 = temp2/5; temp2-=1; */
/* 	beta1 = 0.28*temp1/temp2; beta1p = 0.28*(temp2-temp1*temp3)/pow(temp2,2); */
/* 	temp5 = alpha1/(alpha1+beta1); temp6 = (alpha1p*(alpha1+beta1) - alpha1*(alpha1p+beta1p))/pow(alpha1+beta1,2); */
/* 	temp1 = (temp_Vs - VOLTAGE_[VARNAME_REGISTRY_s_Na]); */
/* 	temp2 = pow(temp5,2)*temp1; */
/* 	temp3 = 2*temp5*temp6*temp1 + pow(temp5,2); */
/* 	temp_I_Na_h_Na = CONDUCTANCE_[VARNAME_REGISTRY_s_Na]*temp2; */
/* 	temp_I_Na = temp_I_Na_h_Na*temp_h_Na; */
/* 	temp_I_Na_Vs = CONDUCTANCE_[VARNAME_REGISTRY_s_Na]*temp_h_Na*temp3; */
/* 	temp_s_Na = CONDUCTANCE_[VARNAME_REGISTRY_s_Na]*pow(temp5,2)*temp_h_Na; */
/* 	if (verbose>3){ printf(" m_Na_inf %f, h_Na %f, I_Na %f, I_Na_Vs %f I_Na_h_Na %f\n",temp5,temp_h_Na,temp_I_Na,temp_I_Na_Vs,temp_I_Na_h_Na);} */
/* 	if (verbose>3){ printf(" %% %% %% %% calculating rhs_h_Na: ");} */
/* 	temp1 = 17-temp_Vs; temp2 = 0.128*exp(temp1/18); temp3 = -temp2/18; */
/* 	alpha2 = temp2; alpha2p = temp3; */
/* 	temp1 = 40-temp_Vs; temp2 = exp(temp1/5); temp3 = -temp2/5; temp2+=1; */
/* 	beta2 = 4/temp2; beta2p = -4*temp3/pow(temp2,2); */
/* 	temp_rhs_h_Na_h_Na = -(alpha2+beta2); */
/* 	temp_rhs_h_Na = alpha2 + temp_h_Na*temp_rhs_h_Na_h_Na; */
/* 	temp_rhs_h_Na_Vs = alpha2p - temp_h_Na*(alpha2p + beta2p); */
/* 	if (verbose>3){ printf(" rhs_h_Na %f rhs_h_Na_Vs %f rhs_h_Na_h_Na %f\n",temp_rhs_h_Na,temp_rhs_h_Na_Vs,temp_rhs_h_Na_h_Na);}   */
/* 	if (verbose>3){ printf(" %% %% %% %% calculating I_K: ");} */
/* 	temp1 = (temp_Vs - VOLTAGE_[VARNAME_REGISTRY_s_K]); */
/* 	temp_I_K_m_K = CONDUCTANCE_[VARNAME_REGISTRY_s_K]*temp1; */
/* 	temp_I_K_Vs = CONDUCTANCE_[VARNAME_REGISTRY_s_K]*temp_m_K; */
/* 	temp_I_K = temp_I_K_Vs*temp1; */
/* 	temp_s_K = temp_I_K_Vs; */
/* 	if (verbose>3){ printf(" I_K %f, I_K_Vs %f I_K_m_K %f\n",temp_I_K,temp_I_K_Vs,temp_I_K_m_K);} */
/* 	if (verbose>3){ printf(" %% %% %% %% calcuating rhs_m_K: ");} */
/* 	temp1 = 35.1-temp_Vs; temp2 = exp(temp1/5); temp3 = -temp2/5; temp2-=1; */
/* 	alpha1 = 0.016*temp1/temp2; alpha1p = 0.16*(-temp2-temp1*temp3)/pow(temp2,2); */
/* 	temp1 = 0.5-0.025*temp_Vs; temp2 = exp(temp1); temp3 = -0.025*temp2; temp2*=0.25; */
/* 	beta1 = temp2; beta1p = 0.25*temp3; */
/* 	temp_rhs_m_K_m_K = -(alpha1+beta1); */
/* 	temp_rhs_m_K = alpha1 + temp_m_K*temp_rhs_m_K_m_K; */
/* 	temp_rhs_m_K_Vs = alpha1p - temp_m_K*(alpha1p + beta1p); */
/* 	if (verbose>3){ printf(" rhs_m_K %f rhs_m_K_Vs %f rhs_m_K_m_K %f\n",temp_rhs_m_K,temp_rhs_m_K_Vs,temp_rhs_m_K_m_K);} */
/* 	if (verbose>3){ printf(" %% %% %% %% calculating I_SD and I_DS: ");} */
/* 	temp1 = (temp_Vs - temp_Vd); */
/* 	temp_I_SD = CONDUCTANCE_SD*temp1; */
/* 	temp_I_SD_Vs = CONDUCTANCE_SD; */
/* 	temp_I_SD_Vd = -CONDUCTANCE_SD; */
/* 	temp_I_DS = -CONDUCTANCE_DS*temp1; */
/* 	temp_I_DS_Vs = -CONDUCTANCE_DS; */
/* 	temp_I_DS_Vd = CONDUCTANCE_DS; */
/* 	if (verbose>3){ printf(" I_SD %f I_SD_Vs %f I_SD_Vd %f I_DS %f I_DS_Vs %f I_DS_Vd %f\n",temp_I_SD,temp_I_SD_Vs,temp_I_SD_Vd,temp_I_DS,temp_I_DS_Vs,temp_I_DS_Vd);}   */
/* 	if (verbose>3){ printf(" %% %% %% %% calculating I_Ca: ");} */
/* 	temp1 = (temp_Vd - VOLTAGE_[VARNAME_REGISTRY_s_Ca]); */
/* 	temp_I_Ca_m_Ca = CONDUCTANCE_[VARNAME_REGISTRY_s_Ca]*2*temp_m_Ca*temp1; */
/* 	temp_I_Ca_Vd = CONDUCTANCE_[VARNAME_REGISTRY_s_Ca]*pow(temp_m_Ca,2); */
/* 	temp_I_Ca = temp_I_Ca_Vd*temp1; */
/* 	temp_s_Ca = temp_I_Ca_Vd; */
/* 	if (verbose>3){ printf(" I_Ca %f, I_Ca_Vd %f I_Ca_m_Ca %f\n",temp_I_Ca,temp_I_Ca_Vd,temp_I_Ca_m_Ca);} */
/* 	if (verbose>3){ printf(" %% %% %% %% calculating rhs_m_Ca: ");} */
/* 	temp1 = temp_Vd - 65; temp2 = exp(-0.072*temp1); temp3 = -0.072*temp2; temp2+=1; */
/* 	alpha1 = 1.6/temp2; alpha1p = -1.6*temp3/pow(temp2,2); */
/* 	temp1 = temp_Vd - 51.1; temp2 = exp(temp1/5); temp3 = temp2/5; temp2-=1; */
/* 	beta1 = 0.02*temp1/temp2; beta1p = 0.02*(temp2-temp1*temp3)/pow(temp2,2); */
/* 	temp_rhs_m_Ca_m_Ca = -(alpha1+beta1); */
/* 	temp_rhs_m_Ca = alpha1 + temp_m_Ca*temp_rhs_m_Ca_m_Ca; */
/* 	temp_rhs_m_Ca_Vd = alpha1p - temp_m_Ca*(alpha1p + beta1p); */
/* 	if (verbose>3){ printf(" rhs_m_Ca %f rhs_m_Ca_Vd %f rhs_m_Ca_m_Ca %f\n",temp_rhs_m_Ca,temp_rhs_m_Ca_Vd,temp_rhs_m_Ca_m_Ca);} */
/* 	if (verbose>3){ printf(" %% %% %% %% calculating rhs_Ca: ");} */
/* 	temp_rhs_Ca = -0.13*temp_I_Ca - temp_Ca/TAU_[VARNAME_REGISTRY_Ca]; */
/* 	temp_rhs_Ca_Vd = -0.13*temp_I_Ca_Vd; */
/* 	temp_rhs_Ca_m_Ca = -0.13*temp_I_Ca_m_Ca; */
/* 	temp_rhs_Ca_Ca = -1.0/TAU_[VARNAME_REGISTRY_Ca]; */
/* 	if (verbose>3){ printf(" rhs_Ca %f rhs_Ca_Vd %f rhs_Ca_m_Ca %f rhs_Ca_Ca %f\n",temp_rhs_Ca,temp_rhs_Ca_Vd,temp_rhs_Ca_m_Ca,temp_rhs_Ca_Ca);} */
/* 	if (verbose>3){ printf(" %% %% %% %% calculating I_KCa: ");} */
/* 	temp1 = (temp_Vd - VOLTAGE_[VARNAME_REGISTRY_s_KCa]); */
/* 	hump(250,1,10,temp_Ca,&temp4,&temp8); */
/* 	temp2 = 0.004*temp_Ca*(1-temp4)+1*temp4; */
/* 	temp3 = 0.004*((1-temp4)-temp_Ca*temp8)+1*temp8; */
/* 	temp_I_KCa_m_KCa = CONDUCTANCE_[VARNAME_REGISTRY_s_KCa]*temp2*temp1; */
/* 	temp_I_KCa_Vd = CONDUCTANCE_[VARNAME_REGISTRY_s_KCa]*temp2*temp_m_KCa; */
/* 	temp_I_KCa = temp_I_KCa_Vd*temp1; */
/* 	temp_I_KCa_Ca = CONDUCTANCE_[VARNAME_REGISTRY_s_KCa]*temp_m_KCa*temp3*temp1; */
/* 	temp_s_KCa = temp_I_KCa_Vd; */
/* 	if (verbose>3){ printf(" I_KCa %f, I_KCa_Vd %f I_KCa_m_KCa %f I_KCa_Ca %f\n",temp_I_KCa,temp_I_KCa_Vd,temp_I_KCa_m_KCa,temp_I_KCa_Ca);} */
/* 	if (verbose>3){ printf(" %% %% %% %% calculating rhs_m_KCa: ");} */
/* 	temp1 = (1.0/11-1.0/27)*temp_Vd + (6.5/27 - 10.0/11); temp2 = exp(temp1)/18.975; temp3 = temp2*(1.0/11-1.0/27); */
/* 	alpha1 = temp2; alpha1p = temp3; */
/* 	alpha2 = 2.0*exp((6.5-temp_Vd)/27); alpha2p = -alpha2/27; */
/* 	beta1 = alpha2 - alpha1; beta1p = alpha2p - alpha1p; */
/* 	beta2 = 0; beta2p = 0; */
/* 	hump(50,1,5,temp_Vd,&temp1,&temp2); */
/* 	alpha3 = alpha1*(1-temp1)+alpha2*(temp1); alpha3p = alpha1p*(1-temp1) + -alpha1*temp2 + alpha2p*temp1 + alpha2*temp2; */
/* 	beta3 = beta1*(1-temp1)+beta2*(temp1); beta3p = beta1p*(1-temp1) + -beta1*temp2 + beta2p*temp1 + beta2*temp2; */
/* 	/\* fix later *\/	alpha3 *= STD_VIEW; alpha3p *= STD_VIEW; beta3 *= STD_VIEW; beta3p *= STD_VIEW; */
/* 	temp_rhs_m_KCa_m_KCa = -(alpha3+beta3); */
/* 	temp_rhs_m_KCa = alpha3 + temp_m_KCa*temp_rhs_m_KCa_m_KCa; */
/* 	temp_rhs_m_KCa_Vd = alpha3p - temp_m_KCa*(alpha3p + beta3p); */
/* 	if (verbose>3){ printf(" rhs_m_KCa %f rhs_m_KCa_Vd %f rhs_m_KCa_m_KCa %f\n",temp_rhs_m_KCa,temp_rhs_m_KCa_Vd,temp_rhs_m_KCa_m_KCa);} */
/* 	if (verbose>3){ printf(" %% %% %% %% calculating I_KAHP:");} */
/* 	temp1 = (temp_Vd - VOLTAGE_[VARNAME_REGISTRY_s_KAHP]); */
/* 	temp_I_KAHP_m_KAHP = CONDUCTANCE_[VARNAME_REGISTRY_s_KAHP]*temp1; */
/* 	temp_I_KAHP_Vd = CONDUCTANCE_[VARNAME_REGISTRY_s_KAHP]*temp_m_KAHP; */
/* 	temp_I_KAHP = temp_I_KAHP_Vd*temp1; */
/* 	temp_s_KAHP = temp_I_KAHP_Vd; */
/* 	if (verbose>3){ printf(" I_KAHP %f, I_KAHP_Vd %f I_KAHP_m_KAHP %f\n",temp_I_KAHP,temp_I_KAHP_Vd,temp_I_KAHP_m_KAHP);} */
/* 	if (verbose>3){ printf(" %% %% %% %% calculating rhs_m_KAHP:");} */
/* 	hump(500,1,10,temp_Ca,&temp4,&temp8); */
/* 	temp2 = 0.00002*temp_Ca*(1-temp4)+0.01*temp4; */
/* 	temp3 = 0.00002*((1-temp4)-temp_Ca*temp8)+0.01*temp8; */
/* 	alpha1 = temp2; alpha1p = temp3; */
/* 	beta1 = 0.001; beta1p = 0; */
/* 	temp_rhs_m_KAHP_m_KAHP = -(alpha1+beta1); */
/* 	temp_rhs_m_KAHP = alpha1 + temp_m_KAHP*temp_rhs_m_KAHP_m_KAHP; */
/* 	temp_rhs_m_KAHP_Ca = alpha1p - temp_m_KAHP*(alpha1p + beta1p); */
/* 	if (verbose>3){ printf(" rhs_m_KAHP %f rhs_m_KAHP_Ca %f rhs_m_KAHP_m_KAHP %f\n",temp_rhs_m_KAHP,temp_rhs_m_KAHP_Ca,temp_rhs_m_KAHP_m_KAHP);} */
/* 	if (verbose>2){ printf(" %% %% %% %% calculating: implicit euler update\n");} */
/* 	temp_rhs_Vs = -(temp_I_leak_S+temp_I_Na+temp_I_K+temp_I_SD+CURRENT_INJECTION_S)/TAU_[VARNAME_REGISTRY_Vs]; */
/* 	temp_rhs_Vs_Vs = -(temp_I_leak_S_Vs+temp_I_Na_Vs+temp_I_K_Vs+temp_I_SD_Vs)/TAU_[VARNAME_REGISTRY_Vs]; */
/* 	temp_rhs_Vs_h_Na = -(temp_I_Na_h_Na)/TAU_[VARNAME_REGISTRY_Vs]; */
/* 	temp_rhs_Vs_m_K = -(temp_I_K_m_K)/TAU_[VARNAME_REGISTRY_Vs]; */
/* 	temp_rhs_Vs_Vd = -(temp_I_SD_Vd)/TAU_[VARNAME_REGISTRY_Vs]; */
/* 	temp_rhs_Vd = -(temp_I_leak_D+temp_I_Ca+temp_I_KCa+temp_I_KAHP+temp_I_DS+temp_I_syn+CURRENT_INJECTION_D)/TAU_[VARNAME_REGISTRY_Vd]; */
/* 	temp_rhs_Vd_Vs = -(temp_I_DS_Vs)/TAU_[VARNAME_REGISTRY_Vd]; */
/* 	temp_rhs_Vd_Vd = -(temp_I_leak_D_Vd+temp_I_Ca_Vd+temp_I_KCa_Vd+temp_I_KAHP_Vd+temp_I_DS_Vd+temp_I_syn_Vd)/TAU_[VARNAME_REGISTRY_Vd]; */
/* 	temp_rhs_Vd_m_Ca = -(temp_I_Ca_m_Ca)/TAU_[VARNAME_REGISTRY_Vd]; */
/* 	temp_rhs_Vd_m_KCa = -(temp_I_KCa_m_KCa)/TAU_[VARNAME_REGISTRY_Vd]; */
/* 	temp_rhs_Vd_Ca = -(temp_I_KCa_Ca)/TAU_[VARNAME_REGISTRY_Vd]; */
/* 	temp_rhs_Vd_m_KAHP = -(temp_I_KAHP_m_KAHP)/TAU_[VARNAME_REGISTRY_Vd]; */
/* 	specialized_inverse_pinsky(1-dt*temp_rhs_Vs_Vs,-dt*temp_rhs_Vs_h_Na,-dt*temp_rhs_Vs_m_K,-dt*temp_rhs_Vs_Vd,-dt*temp_rhs_h_Na_Vs,1-dt*temp_rhs_h_Na_h_Na,-dt*temp_rhs_m_K_Vs,1-dt*temp_rhs_m_K_m_K,-dt*temp_rhs_Vd_Vs,1-dt*temp_rhs_Vd_Vd,-dt*temp_rhs_Vd_m_Ca,-dt*temp_rhs_Vd_m_KCa,-dt*temp_rhs_Vd_Ca,-dt*temp_rhs_Vd_m_KAHP,-dt*temp_rhs_m_Ca_Vd,1-dt*temp_rhs_m_Ca_m_Ca,-dt*temp_rhs_m_KCa_Vd,1-dt*temp_rhs_m_KCa_m_KCa,-dt*temp_rhs_Ca_Vd,-dt*temp_rhs_Ca_m_Ca,1-dt*temp_rhs_Ca_Ca,-dt*temp_rhs_m_KAHP_Ca,1-dt*temp_rhs_m_KAHP_m_KAHP,temp_Vs_start - temp_Vs + dt*temp_rhs_Vs,temp_h_Na_start - temp_h_Na + dt*temp_rhs_h_Na,temp_m_K_start - temp_m_K + dt*temp_rhs_m_K, temp_Vd_start - temp_Vd + dt*temp_rhs_Vd,temp_m_Ca_start - temp_m_Ca + dt*temp_rhs_m_Ca, temp_m_KCa_start - temp_m_KCa + dt*temp_rhs_m_KCa,temp_Ca_start - temp_Ca + dt*temp_rhs_Ca,temp_m_KAHP_start - temp_m_KAHP + dt*temp_rhs_m_KAHP,&update_Vs,&update_h_Na,&update_m_K,&update_Vd,&update_m_Ca,&update_m_KCa,&update_Ca,&update_m_KAHP); */
/* 	norm_old = pow(update_Vs,2)+pow(update_h_Na,2)+pow(update_m_K,2)+pow(update_Vd,2)+pow(update_m_Ca,2)+pow(update_m_KCa,2)+pow(update_Ca,2)+pow(update_m_KAHP,2); */
/* 	temp_Vs += update_Vs; */
/* 	temp_h_Na += update_h_Na; */
/* 	temp_m_K += update_m_K; */
/* 	temp_Vd += update_Vd; */
/* 	temp_m_Ca += update_m_Ca; */
/* 	temp_m_KCa += update_m_KCa; */
/* 	temp_Ca += update_Ca; */
/* 	temp_m_KAHP += update_m_KAHP; */
/* 	temp7=norm_old;norm_old=norm_new;norm_new=temp7; */
/* 	iteration += 1; */
/* 	if (verbose>2){ printf(" %% %% %% norm_new %0.16f norm_old %0.16f\n",norm_new,norm_old);}} */
/*       while (norm_new>norm_threshold && norm_new<norm_old && iteration<iteration_max); */
/*       if (verbose>1){ printf(" %% %% iteration=%d,norm_new %0.16f\n",iteration,norm_new);} */
/*       temp_s_A_start = temp_s_A; */
/*       temp_s_N_start = temp_s_N; */
/*       temp_s_NMDA_start = temp_s_NMDA; */
/*       temp_s_G_start = temp_s_G; */
/*       temp_s_ORN_start = temp_s_ORN; */
/*       temp_h_Na_start = maximum(0,minimum(1,temp_h_Na)); */
/*       temp_m_K_start = maximum(0,minimum(1,temp_m_K)); */
/*       temp_m_KAHP_start = maximum(0,minimum(1,temp_m_KAHP)); */
/*       temp_m_Ca_start = maximum(0,minimum(1,temp_m_Ca)); */
/*       temp_m_KCa_start = maximum(0,minimum(1,temp_m_KCa)); */
/*       temp_Ca_start = temp_Ca; */
/*       temp_Vs_start = temp_Vs; */
/*       temp_Vd_start = temp_Vd; */
/*       microstep += 1; t+=dt; */
/*       if (norm_new>norm_old || iteration==iteration_max || norm_new>norm_threshold){  */
/* 	if (verbose>2){ printf(" %% %% convergence poor, setting error_flag=1\n");}  */
/* 	error_flag=1;}  */
/*       else{  */
/* 	if (verbose>2){ printf(" %% %% convergence good, setting error_flag=0\n");}  */
/* 	error_flag=0;}} */
/*     while (microstep<n->microstep && !error_flag);       */
/*     if (error_flag){ n->microstep*=2;} else{ n->microstep=maximum(1,n->microstep/2);} */
/*     if (verbose>2){ printf(" %% total dt %f error_flag=%d, n->microstep set to %d\n",t-t_start,error_flag,n->microstep);}} */
/*   while(error_flag); */
/*   temp_Vs_end = temp_Vs_start; temp_Vs_start = n->vra[VARNAME_REGISTRY_Vs]; */
/*   temp_Vd_end = temp_Vd_start; temp_Vd_start = n->vra[VARNAME_REGISTRY_Vd]; */
/*   if (verbose>2){ printf(" %% determining spiking properties\n");} */
/*   if (n->spikenext<t_start+dt_start){ */
/*     temp1=-1; */
/*     if (temp_Vs_end>=VOLTAGE_THRESHOLD_S && temp_Vs_start<VOLTAGE_THRESHOLD_S){  */
/*       temp1 = linerootfinder(temp_Vs_start,temp_Vs_end,VOLTAGE_THRESHOLD_S,dt_start);} */
/*     if (temp_Vd_end>=VOLTAGE_THRESHOLD_D && temp_Vd_start<VOLTAGE_THRESHOLD_D){  */
/*       temp1 = linerootfinder(temp_Vd_start,temp_Vd_end,VOLTAGE_THRESHOLD_D,dt_start);} */
/*     if (temp1>=0){ */
/*       n->spiketime_guess = t_start+temp1; */
/*       if (n->spiketime_guess<=n->spikenext){ */
/* 	n->spiketime_guess=-2*dt_start; n->spiketime_guess_flag=0;} */
/*       else /\* if (n->spiketime_guess>n->spikenext) *\/{ 	 */
/* 	n->vra[VARNAME_REGISTRY_spike_flag] += 1; */
/* 	n->spiketime_guess_flag=1;}}} */
/*   n->vra[VARNAME_REGISTRY_spike_flag] *= exp(-dt_start/(TAU_REF/16)); */
/*   n->vra[VARNAME_REGISTRY_s_A] = temp_s_A_start; */
/*   n->vra[VARNAME_REGISTRY_s_N] = temp_s_N_start; */
/*   n->vra[VARNAME_REGISTRY_s_NMDA] = temp_s_NMDA_start; */
/*   n->vra[VARNAME_REGISTRY_s_G] = temp_s_G_start; */
/*   n->vra[VARNAME_REGISTRY_s_ORN] = temp_s_ORN_start; */
/*   n->vra[VARNAME_REGISTRY_h_Na] = temp_h_Na_start; */
/*   n->vra[VARNAME_REGISTRY_m_K] = temp_m_K_start; */
/*   n->vra[VARNAME_REGISTRY_m_KAHP] = temp_m_KAHP_start; */
/*   n->vra[VARNAME_REGISTRY_m_Ca] = temp_m_Ca_start; */
/*   n->vra[VARNAME_REGISTRY_m_KCa] = temp_m_KCa_start; */
/*   n->vra[VARNAME_REGISTRY_Ca] = temp_Ca_start; */
/*   n->vra[VARNAME_REGISTRY_Vs] = temp_Vs_end; */
/*   n->vra[VARNAME_REGISTRY_Vd] = temp_Vd_end;   */
/*   n->vra[VARNAME_REGISTRY_s_Na] = temp_s_Na; */
/*   n->vra[VARNAME_REGISTRY_s_K] = temp_s_K; */
/*   n->vra[VARNAME_REGISTRY_s_Ca] = temp_s_Ca; */
/*   n->vra[VARNAME_REGISTRY_s_KCa] = temp_s_KCa; */
/*   n->vra[VARNAME_REGISTRY_s_KAHP] = temp_s_KAHP; */
/*   if (verbose>1){ printf(" %% exiting vraevolve: \n"); raprintf(n->vra,"double",1,n->nvars," end ");} */
/*   if (verbose>0){ printf(" %% vra:\n"); for (nv=0;nv<n->nvars;nv++){ printf("%s %f \t",GLOBAL_VARNAMES[nv],n->vra[nv]);} printf("\n");} */
/* }   */

/* void jacobian_to_lyapunov_adi_bullshit(double M00,double M01,double M02,double M03,double M10,double M11,double M20,double M22,double M30,double M33,double M34,double M43,double M44,double *output) */
/* { */
/*   /\* return maximum real part of spectrum of */
/*     X    0    1    2    3    4 */
/*     0  M00  M01  M02  M03   . */
/*     1  M10  M11   .    .    . */
/*     2  M20   .   M22   .    . */
/*     3  M30   .    .   M33  M34 */
/*     4   .    .    .   M43  M44 */
/*   *\/ */
/*   int verbose=0; */
/*   integer five=5,twentyfive=25; */
/*   doublecomplex A[25]; */
/*   doublecomplex B[5]; */
/*   integer INFO=0; */
/*   doublecomplex WORK[25]; */
/*   doublereal RWORK[10]; */
/*   char en='N'; */
/*   int nr=0,nc=0; */
/*   A[0 + 0*five].r = M00; A[0 + 0*five].i = 0.0; */
/*   A[1 + 0*five].r = M10; A[1 + 0*five].i = 0.0; */
/*   A[2 + 0*five].r = M20; A[2 + 0*five].i = 0.0; */
/*   A[3 + 0*five].r = M30; A[3 + 0*five].i = 0.0; */
/*   A[4 + 0*five].r = 0.0; A[4 + 0*five].i = 0.0; */
/*   A[0 + 1*five].r = M01; A[0 + 1*five].i = 0.0; */
/*   A[1 + 1*five].r = M11; A[1 + 1*five].i = 0.0; */
/*   A[2 + 1*five].r = 0.0; A[2 + 1*five].i = 0.0; */
/*   A[3 + 1*five].r = 0.0; A[3 + 1*five].i = 0.0; */
/*   A[4 + 1*five].r = 0.0; A[4 + 1*five].i = 0.0; */
/*   A[0 + 2*five].r = M02; A[0 + 2*five].i = 0.0; */
/*   A[1 + 2*five].r = 0.0; A[1 + 2*five].i = 0.0; */
/*   A[2 + 2*five].r = M22; A[2 + 2*five].i = 0.0; */
/*   A[3 + 2*five].r = 0.0; A[3 + 2*five].i = 0.0; */
/*   A[4 + 2*five].r = 0.0; A[4 + 2*five].i = 0.0; */
/*   A[0 + 3*five].r = M03; A[0 + 3*five].i = 0.0; */
/*   A[1 + 3*five].r = 0.0; A[1 + 3*five].i = 0.0; */
/*   A[2 + 3*five].r = 0.0; A[2 + 3*five].i = 0.0; */
/*   A[3 + 3*five].r = M33; A[3 + 3*five].i = 0.0; */
/*   A[4 + 3*five].r = M43; A[4 + 3*five].i = 0.0; */
/*   A[0 + 4*five].r = 0.0; A[0 + 4*five].i = 0.0; */
/*   A[1 + 4*five].r = 0.0; A[1 + 4*five].i = 0.0; */
/*   A[2 + 4*five].r = 0.0; A[2 + 4*five].i = 0.0; */
/*   A[3 + 4*five].r = M34; A[3 + 4*five].i = 0.0; */
/*   A[4 + 4*five].r = M44; A[4 + 4*five].i = 0.0; */
/*   if (verbose){ printf("A=:\n");for (nr=0;nr<five;nr++){ for (nc=0;nc<five;nc++){ printf("%0.3f+%0.3fi\t",A[nr+nc*five].r,A[nr+nc*five].i);} printf("\n");}} */
/*   zgeev_(&en,&en,&five,A,&five,B,NULL,&five,NULL,&five,WORK,&twentyfive,RWORK,&INFO); */
/*   if (INFO!=0){ printf("warning! info=%d in jacobian_to_lyapunov_adi_bullshit\n",(int)INFO);} */
/*   if (verbose){ printf("B=:\n");for (nr=0;nr<five;nr++){ printf("%0.3f+%0.3fi\n",B[nr].r,B[nr].i);}} */
/*   if (output!=NULL){ */
/*     *output=B[0].r; */
/*     for (nr=0;nr<five;nr++){ *output=maximum(*output,B[nr].r);} */
/*     if (verbose){ printf("biggest %f\n",*output);}} */
/* } */

/* void specialized_inverse_adi_bullshit(double M00,double M01,double M02,double M03,double M10,double M11,double M20,double M22,double M30,double M33,double M34,double M44,double M45,double M53,double M55,double R0,double R1,double R2,double R3,double R4,double R5,double *V0,double *V1,double *V2,double *V3,double *V4,double *V5) */
/* { */
/*   /\* specialized 51 flop inverse (row operations) for the matrix */
/*     X    0    1    2    3    4    5 */
/*     0  M00  M01  M02  M03   .    . */
/*     1  M10  M11   .    .    .    . */
/*     2  M20   .   M22   .    .    . */
/*     3  M30   .    .   M33  M34   . */
/*     4   .    .    .    .   M44  M45 */
/*     5   .    .    .   M53   .   M55 */
/*     with RHS [R0 R1 R2 R3 R4 R5] */
/*   *\/ */
/*   double temp=0; */
/*   double M43=0; */
/*   temp = M45/M55; M43 -= temp*M53; R4 -= temp*R5; /\* M45=0; *\/ */
/*   temp = M34/M44; M33 -= temp*M43; R3 -= temp*R4; /\* M34=0; *\/ */
/*   temp = M34/M44; M33 -= temp*M43; R3 -= temp*R4; /\* M34=0; *\/ */
/*   temp = M03/M33; M00 -= temp*M30; R0 -= temp*R3; /\* M03=0; *\/ */
/*   temp = M02/M22; M00 -= temp*M20; R0 -= temp*R2; /\* M02=0; *\/ */
/*   temp = M01/M11; M00 -= temp*M10; R0 -= temp*R1; /\* M01=0; *\/ */
/*   /\* now matrix should look like */
/*     X    0    1    2    3    4    5 */
/*     0  M00   .    .    .    .    . */
/*     1  M10  M11   .    .    .    . */
/*     2  M20   .   M22   .    .    . */
/*     3  M30   .    .   M33   .    . */
/*     4   .    .    .   M43  M44   . */
/*     5   .    .    .   M53   .   M55 */
/*   *\/ */
/*   temp = M10/M00; R1 -= temp*R0; /\* M10=0; *\/ */
/*   temp = M20/M00; R2 -= temp*R0; /\* M20=0; *\/ */
/*   temp = M30/M00; R3 -= temp*R0; /\* M30=0; *\/ */
/*   temp = M43/M33; R4 -= temp*R3; /\* M43=0; *\/ */
/*   temp = M53/M33; R5 -= temp*R3; /\* M53=0; *\/ */
/*   /\* now matrix should be diagonal *\/ */
/*   *V0 = R0/M00; */
/*   *V1 = R1/M11; */
/*   *V2 = R2/M22; */
/*   *V3 = R3/M33; */
/*   *V4 = R4/M44; */
/*   *V5 = R5/M55; */
/* } */

/* void vraevolve_adi_bullshit(struct neuron *n,double t_start,double dt_start) */
/* { */
/*   /\* linearly implicit euler *\/ */
/*   int verbose=0; */
/*   double temp1=0,temp2=0,temp3=0,temp5=0,temp6=0,temp7=0,temp_lyapunov=0; */
/*   double alpha1=0,alpha1p=0,beta1=0,beta1p=0; */
/*   double temp_s_A=0; */
/*   double temp_s_N=0; */
/*   double temp_s_NMDA=0; */
/*   double temp_s_G=0; */
/*   double temp_s_ORN=0; */
/*   double temp_Vs=0; */
/*   double temp_h_Na=0; */
/*   double temp_m_K=0; */
/*   double temp_Vd=0; */
/*   double temp_Ca=0; */
/*   double temp_m_Ca_slow=0; */
/*   double temp_s_Na=0; */
/*   double temp_s_K=0; */
/*   double temp_s_Ca=0; */
/*   double temp_s_KCa=0; */
/*   double temp_I_syn=0; */
/*   double temp_I_syn_Vd=0; */
/*   double temp_I_leak_S=0; */
/*   double temp_I_leak_S_Vs=0; */
/*   double temp_I_leak_D=0; */
/*   double temp_I_leak_D_Vd=0; */
/*   double temp_I_Na=0; */
/*   double temp_I_Na_Vs=0; */
/*   double temp_I_Na_h_Na=0; */
/*   double temp_rhs_h_Na=0; */
/*   double temp_rhs_h_Na_Vs=0; */
/*   double temp_rhs_h_Na_h_Na=0; */
/*   double temp_I_K=0; */
/*   double temp_I_K_Vs=0; */
/*   double temp_I_K_m_K=0; */
/*   double temp_rhs_m_K=0; */
/*   double temp_rhs_m_K_Vs=0; */
/*   double temp_rhs_m_K_m_K=0; */
/*   double temp_I_SD=0; */
/*   double temp_I_SD_Vs=0; */
/*   double temp_I_SD_Vd=0; */
/*   double temp_I_DS=0; */
/*   double temp_I_DS_Vs=0; */
/*   double temp_I_DS_Vd=0; */
/*   double temp_I_Ca=0; */
/*   double temp_I_Ca_Vd=0; */
/*   double temp_I_KCa=0; */
/*   double temp_I_KCa_Vd=0; */
/*   double temp_I_KCa_m_Ca_slow=0; */
/*   double temp_rhs_m_Ca_slow=0; */
/*   double temp_rhs_m_Ca_slow_Ca=0; */
/*   double temp_rhs_m_Ca_slow_m_Ca_slow=0; */
/*   double temp_rhs_Ca=0; */
/*   double temp_rhs_Ca_Vd=0; */
/*   double temp_rhs_Ca_Ca=0; */
/*   double temp_rhs_Vs=0; */
/*   double temp_rhs_Vs_Vs=0; */
/*   double temp_rhs_Vs_h_Na=0; */
/*   double temp_rhs_Vs_m_K=0; */
/*   double temp_rhs_Vs_Vd=0; */
/*   double temp_rhs_Vd=0; */
/*   double temp_rhs_Vd_Vs=0; */
/*   double temp_rhs_Vd_Vd=0; */
/*   double temp_rhs_Vd_m_Ca_slow=0; */
/*   double temp_s_A_start=0; */
/*   double temp_s_N_start=0; */
/*   double temp_s_NMDA_start=0; */
/*   double temp_s_G_start=0; */
/*   double temp_s_ORN_start=0; */
/*   double temp_Vs_start=0; */
/*   double temp_h_Na_start=0; */
/*   double temp_m_K_start=0; */
/*   double temp_Vd_start=0; */
/*   double temp_Ca_start=0; */
/*   double temp_m_Ca_slow_start=0; */
/*   double temp_Vs_end=0; */
/*   double temp_Vd_end=0; */
/*   double dt=dt_start,t=t_start; */
/*   int microstep=0,error_flag=0; */
/*   int iteration=0,iteration_max=4; */
/*   double update_Vs=0; */
/*   double update_h_Na=0; */
/*   double update_m_K=0; */
/*   double update_Vd=0; */
/*   double update_Ca=0; */
/*   double update_m_Ca_slow=0; */
/*   double norm_new=0,norm_old=0,norm_threshold=0.0000001; */
/*   int nv=0; */
/*   if (verbose>1){ printf(" %% [entering vraevolve_adi_bullshit]: \n"); raprintf(n->vra,"double",1,n->nvars," start ");} */
/*   do{ /\* outer loop, choose n->microstep *\/ */
/*     error_flag=0;dt=dt_start/n->microstep;t=t_start; */
/*     if (verbose>1){ printf(" %% initiating pass -- n->microstep %d dt %f\n",n->microstep,dt);} */
/*     temp_s_A_start = n->vra[VARNAME_REGISTRY_s_A]; */
/*     temp_s_N_start = n->vra[VARNAME_REGISTRY_s_N]; */
/*     temp_s_NMDA_start = n->vra[VARNAME_REGISTRY_s_NMDA]; */
/*     temp_s_G_start = n->vra[VARNAME_REGISTRY_s_G]; */
/*     temp_s_ORN_start = n->vra[VARNAME_REGISTRY_s_ORN]; */
/*     temp_Vs_start=n->vra[VARNAME_REGISTRY_Vs]; */
/*     temp_h_Na_start=n->vra[VARNAME_REGISTRY_h_Na]; */
/*     temp_m_K_start=n->vra[VARNAME_REGISTRY_m_K]; */
/*     temp_Vd_start=n->vra[VARNAME_REGISTRY_Vd]; */
/*     temp_Ca_start=n->vra[VARNAME_REGISTRY_Ca]; */
/*     temp_m_Ca_slow_start=n->vra[VARNAME_REGISTRY_m_Ca_slow]; */
/*     microstep=0; */
/*     do{ /\* each microstep *\/ */
/*       if (verbose>1){ printf(" %% %% microstep %d\n",microstep);} */
/*       norm_old=16.0;norm_new=16.0; */
/*       temp_s_A = temp_s_A_start * exp(-dt/TAU_[VARNAME_REGISTRY_s_A]); */
/*       temp_s_N = temp_s_N_start * exp(-dt/TAU_[VARNAME_REGISTRY_s_N]); */
/*       temp_s_NMDA = exp(-dt/TAU_[VARNAME_REGISTRY_s_NMDA])*(temp_s_NMDA_start + temp_s_N_start*(-1 + exp(dt*(1.0/TAU_[VARNAME_REGISTRY_s_NMDA]-1.0/TAU_[VARNAME_REGISTRY_s_N])))/(1.0/TAU_[VARNAME_REGISTRY_s_NMDA]-1.0/TAU_[VARNAME_REGISTRY_s_N])); */
/*       temp_s_G = temp_s_G_start * exp(-dt/TAU_[VARNAME_REGISTRY_s_G]); */
/*       temp_s_ORN = temp_s_ORN_start * exp(-dt/TAU_[VARNAME_REGISTRY_s_ORN]); */
/*       temp_h_Na = temp_h_Na_start; */
/*       temp_m_K = temp_m_K_start; */
/*       temp_Ca = temp_Ca_start; */
/*       temp_m_Ca_slow = temp_m_Ca_slow_start; */
/*       temp_Vs = temp_Vs_start; */
/*       temp_Vd = temp_Vd_start; */
/*       iteration=0; */
/*       do{ /\* inner loop, linearly implicit euler *\/ */
/* 	if (verbose>2){ printf(" %% %% %% iteration %d\n",iteration);} */
/* 	if (verbose>3){ printf(" %% %% %% %% calculating I_syn: ");} */
/* 	temp1 = exp(-0.062*(temp_Vd+VOLTAGE_[VARNAME_REGISTRY_s_LEAK_D]-60)); temp2 = temp1*-0.062; temp1+=1; temp3 = -temp2/pow(temp1,2); */
/* 	temp_I_syn = temp_s_ORN*(temp_Vd - VOLTAGE_[VARNAME_REGISTRY_s_ORN]) + temp_s_A*(temp_Vd - VOLTAGE_[VARNAME_REGISTRY_s_A]) + temp_s_NMDA/temp1*(temp_Vd - VOLTAGE_[VARNAME_REGISTRY_s_NMDA]) + temp_s_G*(temp_Vd - VOLTAGE_[VARNAME_REGISTRY_s_G]); */
/* 	temp_I_syn_Vd = temp_s_ORN + temp_s_A + temp_s_NMDA*(1.0/temp1 + temp3*(temp_Vd-VOLTAGE_[VARNAME_REGISTRY_s_NMDA])) + temp_s_G; */
/* 	if (verbose>3){ printf(" I_syn %f\n",temp_I_syn);} */
/* 	if (verbose>3){ printf(" %% %% %% %% calculating I_leak_S and I_leak_D: ");} */
/* 	temp_I_leak_S_Vs = CONDUCTANCE_[VARNAME_REGISTRY_Vs]; */
/* 	temp_I_leak_S = temp_I_leak_S_Vs*(temp_Vs - VOLTAGE_[VARNAME_REGISTRY_s_LEAK_S]); */
/* 	temp_I_leak_D_Vd = CONDUCTANCE_[VARNAME_REGISTRY_Vd]; */
/* 	temp_I_leak_D = temp_I_leak_D_Vd*(temp_Vd - VOLTAGE_[VARNAME_REGISTRY_s_LEAK_D]); */
/* 	if (verbose>3){ printf(" I_leak_S %f, I_leak_S_Vs %f I_leak_D %f I_leak_D_Vd %f\n",temp_I_leak_S,temp_I_leak_S_Vs,temp_I_leak_D,temp_I_leak_D_Vd);} */
/* 	if (verbose>3){ printf(" %% %% %% %% calculating I_Na: ");} */
/* 	temp1 = (temp_Vs+33)/10; temp2 = -exp(-temp1); temp3 = -temp2/10; temp2+=1; */
/* 	alpha1 = temp1/temp2; alpha1p = (0.1*temp2-temp1*temp3)/pow(temp2,2); */
/* 	temp1 = (temp_Vs+53.7)/12; temp2 = 4*exp(-temp1); temp3 = -temp2/12; */
/* 	beta1 = temp2; beta1p = temp3; */
/* 	temp5 = alpha1/(alpha1+beta1); temp6 = (alpha1p*(alpha1+beta1) - alpha1*(alpha1p+beta1p))/pow(alpha1+beta1,2); */
/* 	temp1 = (temp_Vs - VOLTAGE_[VARNAME_REGISTRY_s_Na]); */
/* 	temp2 = pow(temp5,3)*temp1; */
/* 	temp3 = 3*pow(temp5,2)*temp6*temp1 + pow(temp5,3); */
/* 	temp_I_Na_h_Na = CONDUCTANCE_[VARNAME_REGISTRY_s_Na]*temp2; */
/* 	temp_I_Na = temp_I_Na_h_Na*temp_h_Na; */
/* 	temp_I_Na_Vs = CONDUCTANCE_[VARNAME_REGISTRY_s_Na]*temp_h_Na*temp3; */
/* 	temp_s_Na = CONDUCTANCE_[VARNAME_REGISTRY_s_Na]*pow(temp5,3)*temp_h_Na; */
/* 	if (verbose>3){ printf(" m_Na_inf %f, h_Na %f, I_Na %f, I_Na_Vs %f I_Na_h_Na %f\n",temp5,temp_h_Na,temp_I_Na,temp_I_Na_Vs,temp_I_Na_h_Na);} */
/* 	if (verbose>3){ printf(" %% %% %% %% calculating rhs_h_Na: ");} */
/* 	temp1 = (temp_Vs+50)/10; temp2 = 0.07*exp(-temp1); temp3 = -temp2/10; */
/* 	alpha1 = temp2; alpha1p = temp3; */
/* 	temp1 = (temp_Vs+20)/10; temp2 = exp(-temp1); temp3 = -temp2/10; temp2+=1; */
/* 	beta1 = 1/temp2; beta1p = -temp3/pow(temp2,2); */
/* 	temp_rhs_h_Na_h_Na = -4*(alpha1+beta1); */
/* 	temp_rhs_h_Na = 4*alpha1 + temp_h_Na*temp_rhs_h_Na_h_Na; */
/* 	temp_rhs_h_Na_Vs = 4*(alpha1p - temp_h_Na*(alpha1p + beta1p)); */
/* 	if (verbose>3){ printf(" rhs_h_Na %f rhs_h_Na_Vs %f rhs_h_Na_h_Na %f\n",temp_rhs_h_Na,temp_rhs_h_Na_Vs,temp_rhs_h_Na_h_Na);} */
/* 	if (verbose>3){ printf(" %% %% %% %% calculating I_K: ");} */
/* 	temp1 = (temp_Vs - VOLTAGE_[VARNAME_REGISTRY_s_K]); */
/* 	temp_I_K_m_K = CONDUCTANCE_[VARNAME_REGISTRY_s_K]*4*pow(temp_m_K,3)*temp1; */
/* 	temp_I_K_Vs = CONDUCTANCE_[VARNAME_REGISTRY_s_K]*pow(temp_m_K,4); */
/* 	temp_I_K = temp_I_K_Vs*temp1; */
/* 	temp_s_K = temp_I_K_Vs; */
/* 	if (verbose>3){ printf(" I_K %f, I_K_Vs %f I_K_m_K %f\n",temp_I_K,temp_I_K_Vs,temp_I_K_m_K);} */
/* 	if (verbose>3){ printf(" %% %% %% %% calcuating rhs_m_K: ");} */
/* 	temp1 = (temp_Vs+34)/10; temp2 = -exp(-temp1); temp3 = -temp2/10; temp2+=1; */
/* 	alpha1 = temp1/temp2/10; alpha1p = (0.1*temp2-temp1*temp3)/pow(temp2,2)/10; */
/* 	temp1 = (temp_Vs+44)/25; temp2 = 0.125*exp(-temp1); temp3 = -temp2/25; */
/* 	beta1 = temp2; beta1p = temp3; */
/* 	temp_rhs_m_K_m_K = -4*(alpha1+beta1); */
/* 	temp_rhs_m_K = 4*alpha1 + temp_m_K*temp_rhs_m_K_m_K; */
/* 	temp_rhs_m_K_Vs = 4*(alpha1p - temp_m_K*(alpha1p + beta1p)); */
/* 	if (verbose>3){ printf(" rhs_m_K %f rhs_m_K_Vs %f rhs_m_K_m_K %f\n",temp_rhs_m_K,temp_rhs_m_K_Vs,temp_rhs_m_K_m_K);} */
/* 	if (verbose>3){ printf(" %% %% %% %% calculating I_SD and I_DS: ");} */
/* 	temp1 = (temp_Vs - temp_Vd); */
/* 	temp_I_SD = CONDUCTANCE_SD*temp1; */
/* 	temp_I_SD_Vs = CONDUCTANCE_SD; */
/* 	temp_I_SD_Vd = -CONDUCTANCE_SD; */
/* 	temp_I_DS = -CONDUCTANCE_DS*temp1; */
/* 	temp_I_DS_Vs = -CONDUCTANCE_DS; */
/* 	temp_I_DS_Vd = CONDUCTANCE_DS; */
/* 	if (verbose>3){ printf(" I_SD %f I_SD_Vs %f I_SD_Vd %f I_DS %f I_DS_Vs %f I_DS_Vd %f\n",temp_I_SD,temp_I_SD_Vs,temp_I_SD_Vd,temp_I_DS,temp_I_DS_Vs,temp_I_DS_Vd);} */
/* 	if (verbose>3){ printf(" %% %% %% %% calculating I_Ca: ");} */
/* 	temp1 = (temp_Vd+20)/9; temp2 = exp(-temp1); temp3 = -temp2/9; temp2+=1; */
/* 	alpha1 = 1/temp2; alpha1p = -temp3/pow(temp2,2); */
/* 	temp1 = (temp_Vd - VOLTAGE_[VARNAME_REGISTRY_s_Ca]); */
/* 	temp_I_Ca = CONDUCTANCE_[VARNAME_REGISTRY_s_Ca]*pow(alpha1,2)*temp1; */
/* 	temp_I_Ca_Vd = CONDUCTANCE_[VARNAME_REGISTRY_s_Ca]*(2*alpha1*alpha1p*temp1 + pow(alpha1,2)); */
/* 	temp_s_Ca = temp_I_Ca/temp1; */
/* 	if (verbose>3){ printf(" I_Ca %f, I_Ca_Vd %f\n",temp_I_Ca,temp_I_Ca_Vd);} */
/* 	if (verbose>3){ printf(" %% %% %% %% calculating I_KCa: ");} */
/* 	temp1 = (temp_Vd - VOLTAGE_[VARNAME_REGISTRY_s_KCa]); */
/* 	temp_I_KCa = CONDUCTANCE_[VARNAME_REGISTRY_s_KCa]*temp_m_Ca_slow*temp1; */
/* 	temp_I_KCa_Vd = CONDUCTANCE_[VARNAME_REGISTRY_s_KCa]*temp_m_Ca_slow; */
/* 	temp_I_KCa_m_Ca_slow = CONDUCTANCE_[VARNAME_REGISTRY_s_KCa]*temp1; */
/* 	temp_s_KCa = temp_I_KCa/temp1; */
/* 	if (verbose>3){ printf(" I_KCa %f, I_KCa_Vd %f, I_KCa_m_Ca_slow %f\n",temp_I_KCa,temp_I_KCa_Vd,temp_I_KCa_m_Ca_slow);} */
/* 	if (verbose>3){ printf(" %% %% %% %% calcuating rhs_m_Ca_slow: ");} */
/* 	hump(3.00,1,0.1,temp_Ca,&temp1,&temp2);  */
/* 	temp_rhs_m_Ca_slow = temp1*(1-temp_m_Ca_slow) - temp_m_Ca_slow/TAU_[VARNAME_REGISTRY_m_Ca_slow]; */
/* 	temp_rhs_m_Ca_slow_m_Ca_slow = -temp1 -1.0/TAU_[VARNAME_REGISTRY_m_Ca_slow]; */
/* 	temp_rhs_m_Ca_slow_Ca = temp2*(1-temp_m_Ca_slow); */
/* 	if (verbose>3){ printf(" rhs_m_Ca_slow %f m_Ca_slow_m_Ca_slow %f m_Ca_slow_Ca %f\n",temp_rhs_m_Ca_slow,temp_rhs_m_Ca_slow_m_Ca_slow,temp_rhs_m_Ca_slow_Ca);}	 */
/* 	if (verbose>3){ printf(" %% %% %% %% calculating rhs_Ca: ");} */
/* 	temp_rhs_Ca = -0.00175*temp_I_Ca - temp_Ca/TAU_[VARNAME_REGISTRY_Ca]; */
/* 	temp_rhs_Ca_Ca = -1/TAU_[VARNAME_REGISTRY_Ca]; */
/* 	temp_rhs_Ca_Vd = -0.00175*temp_I_Ca_Vd; */
/* 	if (verbose>3){ printf(" rhs_Ca %f, rhs_Ca_Vd %f, rhs_Ca_Ca %f\n",temp_rhs_Ca,temp_rhs_Ca_Vd,temp_rhs_Ca_Ca);} */
/* 	if (verbose>2){ printf(" %% %% %% %% calculating: implicit euler update\n");} */
/* 	temp_rhs_Vs = -(temp_I_leak_S+temp_I_Na+temp_I_K+temp_I_SD+CURRENT_INJECTION_S)/TAU_[VARNAME_REGISTRY_Vs]; */
/* 	temp_rhs_Vs_Vs = -(temp_I_leak_S_Vs+temp_I_Na_Vs+temp_I_K_Vs+temp_I_SD_Vs)/TAU_[VARNAME_REGISTRY_Vs]; */
/* 	temp_rhs_Vs_h_Na = -(temp_I_Na_h_Na)/TAU_[VARNAME_REGISTRY_Vs]; */
/* 	temp_rhs_Vs_m_K = -(temp_I_K_m_K)/TAU_[VARNAME_REGISTRY_Vs]; */
/* 	temp_rhs_Vs_Vd = -(temp_I_SD_Vd)/TAU_[VARNAME_REGISTRY_Vs]; */
/* 	temp_rhs_Vd = -(temp_I_leak_D+temp_I_Ca+temp_I_KCa+temp_I_DS+temp_I_syn+CURRENT_INJECTION_D)/TAU_[VARNAME_REGISTRY_Vd]; */
/* 	temp_rhs_Vd_Vd = -(temp_I_leak_D_Vd+temp_I_Ca_Vd+temp_I_KCa_Vd+temp_I_DS_Vd+temp_I_syn_Vd)/TAU_[VARNAME_REGISTRY_Vd]; */
/* 	temp_rhs_Vd_Vs = -(temp_I_DS_Vs)/TAU_[VARNAME_REGISTRY_Vd]; */
/* 	temp_rhs_Vd_m_Ca_slow = -(temp_I_KCa_m_Ca_slow)/TAU_[VARNAME_REGISTRY_Vd]; */
/* 	specialized_inverse_adi_bullshit(1-dt*temp_rhs_Vs_Vs,-dt*temp_rhs_Vs_h_Na,-dt*temp_rhs_Vs_m_K,-dt*temp_rhs_Vs_Vd,-dt*temp_rhs_h_Na_Vs,1-dt*temp_rhs_h_Na_h_Na,-dt*temp_rhs_m_K_Vs,1-dt*temp_rhs_m_K_m_K,-dt*temp_rhs_Vd_Vs,1-dt*temp_rhs_Vd_Vd,-dt*temp_rhs_Vd_m_Ca_slow,1-dt*temp_rhs_m_Ca_slow_m_Ca_slow,-dt*temp_rhs_m_Ca_slow_Ca,-dt*temp_rhs_Ca_Vd,1-dt*temp_rhs_Ca_Ca,temp_Vs_start - temp_Vs + dt*temp_rhs_Vs,temp_h_Na_start - temp_h_Na + dt*temp_rhs_h_Na,temp_m_K_start - temp_m_K + dt*temp_rhs_m_K,temp_Vd_start - temp_Vd + dt*temp_rhs_Vd,temp_m_Ca_slow_start - temp_m_Ca_slow + dt*temp_rhs_m_Ca_slow,temp_Ca_start - temp_Ca + dt*temp_rhs_Ca,&update_Vs,&update_h_Na,&update_m_K,&update_Vd,&update_m_Ca_slow,&update_Ca); */
/* 	switch (LYAPUNOV_BOTHER){ */
/* 	case 0: break; */
/* 	case 1: printf(" %% warning! jacobian outdated in vraevolve_adi_bullshit!\n"); */
/* 	case 2: jacobian_to_lyapunov_morrislecar(temp_rhs_Vs_Vs,temp_rhs_Vs_h_Na,temp_rhs_Vs_m_K,temp_rhs_h_Na_Vs,temp_rhs_h_Na_h_Na,temp_rhs_m_K_Vs,temp_rhs_m_K_m_K,&temp_lyapunov); break; */
/* 	default: break;} */
/* 	norm_old = pow(update_Vs,2)+pow(update_h_Na,2)+pow(update_m_K,2)+pow(update_Vd,2)+pow(update_m_Ca_slow,2)+pow(update_Ca,2); */
/* 	temp_Vs += update_Vs; */
/* 	temp_h_Na += update_h_Na; */
/* 	temp_m_K += update_m_K; */
/* 	temp_Vd += update_Vd; */
/* 	temp_m_Ca_slow += update_m_Ca_slow; */
/* 	temp_Ca += update_Ca; */
/* 	temp7=norm_old;norm_old=norm_new;norm_new=temp7; */
/* 	iteration += 1; */
/* 	if (verbose>2){ printf(" %% %% %% norm_new %0.16f norm_old %0.16f\n",norm_new,norm_old);}} */
/*       while (norm_new>norm_threshold && norm_new<norm_old && iteration<iteration_max); */
/*       if (verbose>1){ printf(" %% %% iteration=%d,norm_new %0.16f\n",iteration,norm_new);} */
/*       temp_s_A_start = temp_s_A; */
/*       temp_s_N_start = temp_s_N; */
/*       temp_s_NMDA_start = temp_s_NMDA; */
/*       temp_s_G_start = temp_s_G; */
/*       temp_s_ORN_start = temp_s_ORN; */
/*       temp_h_Na_start = maximum(0,minimum(1,temp_h_Na)); */
/*       temp_m_K_start = maximum(0,minimum(1,temp_m_K)); */
/*       temp_m_Ca_slow_start = maximum(0,minimum(1,temp_m_Ca_slow)); */
/*       temp_Ca_start = temp_Ca; */
/*       temp_Vs_start = temp_Vs; */
/*       temp_Vd_start = temp_Vd; */
/*       microstep += 1; t+=dt; */
/*       if (norm_new>norm_old || iteration==iteration_max || norm_new>norm_threshold){ */
/* 	if (verbose>2){ printf(" %% %% convergence poor, setting error_flag=1\n");} */
/* 	error_flag=1;} */
/*       else{ */
/* 	if (verbose>2){ printf(" %% %% convergence good, setting error_flag=0\n");} */
/* 	error_flag=0;}} */
/*     while (microstep<n->microstep && !error_flag); */
/*     if (error_flag){ n->microstep*=2;} else{ n->microstep=maximum(1,n->microstep/2);} */
/*     if (verbose>2){ printf(" %% total dt %f error_flag=%d, n->microstep set to %d\n",t-t_start,error_flag,n->microstep);}} */
/*   while(error_flag); */
/*   temp_Vs_end = temp_Vs_start; temp_Vs_start = n->vra[VARNAME_REGISTRY_Vs]; */
/*   temp_Vd_end = temp_Vd_start; temp_Vd_start = n->vra[VARNAME_REGISTRY_Vd]; */
/*   if (verbose>2){ printf(" %% determining spiking properties\n");} */
/*   if (n->spikenext<t_start+dt_start){ */
/*     temp1=-1; */
/*     if (temp_Vs_end>=VOLTAGE_THRESHOLD_S && temp_Vs_start<VOLTAGE_THRESHOLD_S){ */
/*       temp1 = linerootfinder(temp_Vs_start,temp_Vs_end,VOLTAGE_THRESHOLD_S,dt_start);} */
/*     if (temp_Vd_end>=VOLTAGE_THRESHOLD_D && temp_Vd_start<VOLTAGE_THRESHOLD_D){ */
/*       temp1 = linerootfinder(temp_Vd_start,temp_Vd_end,VOLTAGE_THRESHOLD_D,dt_start);} */
/*     if (temp1>=0){ */
/*       n->spiketime_guess = t_start+temp1; */
/*       if (n->spiketime_guess<=n->spikenext){ */
/* 	n->spiketime_guess=-2*dt_start; n->spiketime_guess_flag=0;} */
/*       else /\* if (n->spiketime_guess>n->spikenext) *\/{ */
/* 	n->vra[VARNAME_REGISTRY_spike_flag] += 1; */
/* 	n->spiketime_guess_flag=1;}}} */
/*   n->vra[VARNAME_REGISTRY_spike_flag] *= exp(-dt_start/(TAU_REF/16)); */
/*   n->vra[VARNAME_REGISTRY_s_A] = temp_s_A_start; */
/*   n->vra[VARNAME_REGISTRY_s_N] = temp_s_N_start; */
/*   n->vra[VARNAME_REGISTRY_s_NMDA] = temp_s_NMDA_start; */
/*   n->vra[VARNAME_REGISTRY_s_G] = temp_s_G_start; */
/*   n->vra[VARNAME_REGISTRY_s_ORN] = temp_s_ORN_start; */
/*   n->vra[VARNAME_REGISTRY_h_Na] = temp_h_Na_start; */
/*   n->vra[VARNAME_REGISTRY_m_K] = temp_m_K_start; */
/*   n->vra[VARNAME_REGISTRY_m_Ca_slow] = temp_m_Ca_slow_start; */
/*   n->vra[VARNAME_REGISTRY_Ca] = temp_Ca_start; */
/*   n->vra[VARNAME_REGISTRY_Vs] = temp_Vs_end; */
/*   n->vra[VARNAME_REGISTRY_Vd] = temp_Vd_end; */
/*   n->vra[VARNAME_REGISTRY_s_Na] = temp_s_Na; */
/*   n->vra[VARNAME_REGISTRY_s_K] = temp_s_K; */
/*   n->vra[VARNAME_REGISTRY_s_Ca] = temp_s_Ca; */
/*   n->vra[VARNAME_REGISTRY_s_KCa] = temp_s_KCa; */
/*   n->vra[VARNAME_REGISTRY_lyapunov] = temp_lyapunov; */
/*   if (verbose>1){ printf(" %% exiting vraevolve: \n"); raprintf(n->vra,"double",1,n->nvars," end ");} */
/*   if (verbose>0){ printf(" %% vra:\n"); for (nv=0;nv<n->nvars;nv++){ printf("%s %f \t",GLOBAL_VARNAMES[nv],n->vra[nv]);} printf("\n");} */
/* } */

