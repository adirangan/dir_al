/* /\* here are cortex-dependent llists functions *\/ */
/* void lrupdate(struct neuronarray *,double,double); */
/* int spiketime_compare(void *,void *); */
/* int lgnangle_compare(void *,void *); */
/* void llistsort(struct litem *,struct litem *,int,(int);(*compare);(void *,void *);); */
/* void llistunmake(struct llist **,struct llist **,struct llist *); */
/* void llistremake(struct llist ***,struct llist ***,struct llist **); */

/* /\* Here are the neuron/neuronarray functions *\/ */
/* void makeshuffle(int *,int,int,int); */
/* fftw_complex * makelrkernel(int); */
/* void neuronsetang(struct neuronarray *); */
/* void neuronmake(struct neuronarray *,int,int); */
/* void neurontfree(struct neuron *); */
/* struct neuron * nget(struct neuronarray *,int,int); */
/* void nset(struct neuronarray *,int,int,void *); */
/* struct neuronarray * neuronarraymake(int,int); */
/* void neuronarraytfree(struct neuronarray *); */
/* void setblocksize(double);; */
/* void setinputrate(struct neuronarray *,double); */
/* void spikeinput(struct neuronarray *,double,double); */
/* int spikescan(struct neuronarray *,struct llist **,struct llist **,double,double); */
/* double spikeguesseif(struct neuron *,double,double); */
/* double spikeguess2(struct neuron *,double,double); */
/* double spikeguess1(struct neuron *,double,double); */
/* void spikesort(struct llist **,struct llist *); */
/* void llprintf(struct llist *); */
/* void llraprintf(struct llist **,int,int); */
/* void lreconnect(struct llist *,struct llist *); */
/* int distance(struct llist *,struct llist *); */
/* int ilink(struct neuron *,double *,double *,double *); */
/* int slink(struct neuron *,struct neuron *,double *,double *,double *); */
/* double eifrhs(double,double,double,double,double); */
/* void slavecalc(double,double,double,double *,double *,double *); */
/* void gizmointegrate2(struct neuron *,double,double,double,double,double); */
/* void gizmointegrateeif(struct neuron *,double,double,double,double,double); */
/* void clumpcorrect(struct llist *,struct llist *,double,double); */
/* void spikecorrect(struct llist *,double,double); */
/* void spikeconduct(struct llist *,double,double,int *,int *,int *,int *); */
/* double quadrootfinder(double,double,double,double,double); */
/* double cuberootfinder(double,double,double,double,double,double); */
/* void gizmoconduct(struct llist **,struct llist **,double,double); */
/* void gizmospiked(struct neuron *,double,double); */

/* here are cortex-dependent llists functions */

void lrupdate(struct neuronarray *Nra,double t,double DT)
{
  int nr=0,nc=0,nl=0,tab=0;
  int bign=NPIEROWS*NPIECOLS*PIE_ROW_DIA*PIE_COL_DIA;
  struct neuron *n=NULL;
  double temp=0,failure=0;
  for (nl=0;nl<NLRKERNELS;nl++){ for (nr=0;nr<bign;nr++){ 
    GLOBAL_LRSWAP_X[nl][nr][0]=0; GLOBAL_LRSWAP_X[nl][nr][1]=0;}}
  for (nr=0;nr<Nra->rows;nr++){ for (nc=0;nc<Nra->cols;nc++){
    n = nget(Nra,nr,nc);
    tab = n->pierow+n->piecol*NPIEROWS+n->ang*NPIEROWS*NPIECOLS;
    if (n->spiketime >= t && n->spiketime <=t+DT){
      switch((int)*(n->t2s)){
      case 0: /* IC */ GLOBAL_LRSWAP_X[n->lr_kernel_type][tab][0] += CS_ICLR; break;
      case 1: /* IS */ GLOBAL_LRSWAP_X[n->lr_kernel_type][tab][0] += CS_ISLR; break;
      case 2: /* EC */ GLOBAL_LRSWAP_X[n->lr_kernel_type][tab][0] += CS_ECLR; break;
      case 3: /* ES */ GLOBAL_LRSWAP_X[n->lr_kernel_type][tab][0] += CS_ESLR; break;
      default: printf("error! wrong (int)*(n->t2s)=%d in lrupdate\n",(int)*(n->t2s));}}}}
  for (nl=0;nl<NLRKERNELS;nl++){
    fftw_execute_dft(GLOBAL_FFTWPLAN_SWAP_FORWARD[nl],GLOBAL_LRSWAP_X[nl],GLOBAL_LRSWAP_K[nl]);
    for (nr=0;nr<bign;nr++){
      GLOBAL_LRSWAP_X[nl][nr][0] = GLOBAL_LRKERNEL_K[nl][nr][0]*GLOBAL_LRSWAP_K[nl][nr][0] - GLOBAL_LRKERNEL_K[nl][nr][1]*GLOBAL_LRSWAP_K[nl][nr][1];
      GLOBAL_LRSWAP_X[nl][nr][1] = GLOBAL_LRKERNEL_K[nl][nr][0]*GLOBAL_LRSWAP_K[nl][nr][1] + GLOBAL_LRKERNEL_K[nl][nr][1]*GLOBAL_LRSWAP_K[nl][nr][0];}
    fftw_execute_dft(GLOBAL_FFTWPLAN_SWAP_BACKWARD[nl],GLOBAL_LRSWAP_X[nl],GLOBAL_LRSWAP_K[nl]);}
  for (nr=0;nr<Nra->rows;nr++){ for (nc=0;nc<Nra->cols;nc++){
    n = nget(Nra,nr,nc);
    tab = n->ang + n->pierow*PIE_ROW_DIA*PIE_COL_DIA + n->piecol*PIE_ROW_DIA*PIE_COL_DIA*NPIEROWS;
    temp=0;
    for (nl=0;nl<NLRKERNELS;nl++){ temp += GLOBAL_LRSWAP_K[nl][n->pierow+n->piecol*NPIEROWS+GLOBAL_SHUFFLES[tab]*NPIEROWS*NPIECOLS][0];}
    temp /= (double)bign;
    failure = rand01<LR_RATE/LR_RATE;
    switch((int)*(n->t2s)){
    case 0: /* IC */ *(n->sA) += LR_TO_AMPA*failure*temp*CS_LRIC; *(n->sN) += (1-LR_TO_AMPA)*failure*temp*CS_LRIC; break;
    case 1: /* IS */ *(n->sA) += LR_TO_AMPA*failure*temp*CS_LRIS; *(n->sN) += (1-LR_TO_AMPA)*failure*temp*CS_LRIS; break;
    case 2: /* EC */ *(n->sA) += LR_TO_AMPA*failure*temp*CS_LREC; *(n->sN) += (1-LR_TO_AMPA)*failure*temp*CS_LREC; break;
    case 3: /* ES */ *(n->sA) += LR_TO_AMPA*failure*temp*CS_LRES; *(n->sN) += (1-LR_TO_AMPA)*failure*temp*CS_LRES; break;
    default: printf("error! wrong (int)*(n->t2s)=%d in lrupdate\n",(int)*(n->t2s));}}}
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

int lgnangle_compare(void *item1,void *item2)
{
  /* if n1->lgnangle > n2->lgnangle */
  int comparison=0;
  struct neuron *n1=(struct neuron *) item1;
  struct neuron *n2=(struct neuron *) item2;
  if (n1->lgnangle < n2->lgnangle){ comparison = -1;}
  else if (n1->lgnangle > n2->lgnangle){ comparison = +1;}
  else if (n1->lgnangle == n2->lgnangle && n1 < n2 ){ comparison = -1;}
  else if (n1->lgnangle == n2->lgnangle && n1 > n2 ){ comparison = +1;}
  return comparison;
}

void llistunmake(struct llist **Lra1, struct llist **Lra2, struct llist *LL)
{
  /* this tfrees (GLOBAL_NBLOCKS_TALL by GLOBAL_NBLOCKS_WIDE) element llist arrays Lra1 and Lra2 as well as llist llist LL */
  struct litem *l=NULL;
  int ar=0,ac=0;
  struct llist *L=NULL;
  if (Lra1!=NULL){
    for (ar=0;ar<GLOBAL_NBLOCKS_TALL;ar++){ for (ac=0;ac<GLOBAL_NBLOCKS_WIDE;ac++){
      L = Lra1[ar+ac*GLOBAL_NBLOCKS_TALL]; if (L!=NULL){ llisttfree(L); L=NULL;}}}
    tfree(Lra1); Lra1=NULL;}
  if (Lra2!=NULL){
    for (ar=0;ar<GLOBAL_NBLOCKS_TALL;ar++){ for (ac=0;ac<GLOBAL_NBLOCKS_WIDE;ac++){
      L = Lra2[ar+ac*GLOBAL_NBLOCKS_TALL]; if (L!=NULL){ llisttfree(L); L=NULL;}}}
    tfree(Lra2); Lra2=NULL;}
  if (LL!=NULL){
    l = LL->first;
    while (l!=NULL){ L = (struct llist *) l->item; llisttfree(L); L=NULL; l = l->child;}
    llisttfree(LL); LL=NULL;}
}

void llistremake(struct llist ***Lra1, struct llist ***Lra2, struct llist **LL)
{
  /* this remakes llist arrays Lra1 and Lra2 and llist LL, assuming that they are all NULL */
  int ar=0,ac=0;
  *Lra1 = (struct llist **) tmalloc(sizeof(struct llist *)*GLOBAL_NBLOCKS_TALL*GLOBAL_NBLOCKS_WIDE);
  *Lra2 = (struct llist **) tmalloc(sizeof(struct llist *)*GLOBAL_NBLOCKS_TALL*GLOBAL_NBLOCKS_WIDE);
  for (ar=0;ar<GLOBAL_NBLOCKS_TALL;ar++){ for (ac=0;ac<GLOBAL_NBLOCKS_WIDE;ac++){
    (*Lra1)[ar + ac*GLOBAL_NBLOCKS_TALL] = llistmake();
    (*Lra2)[ar + ac*GLOBAL_NBLOCKS_TALL] = llistmake();}}
  *LL = llistmake();
}

/* Here are the neuron/neuronarray functions */

void makeshuffle(int *ra,int length,int subd,int rseed)
{
  /* this makes a shuffled index array of length where each index is periodically within subd of its original position */
  /* the output ra must be pre-allocated with length */
  int nr=0,nr2=0,d=0,temp=0;
  int np=0,npasses=1;
  srand(rseed);
  for (nr=0;nr<length;nr++){ ra[nr]=nr;}
  for (np=0;np<npasses;np++){
  for (nr=0;nr<length;nr++){
    d = rand()%(1+2*subd)-subd;
    nr2 = periodize(nr+d,0,length);
    if (fabs(periodize(ra[nr]-nr2,-length/2,length/2))<=subd && fabs(periodize(ra[nr2]-nr,-length/2,length/2))<=subd){
      temp=ra[nr]; ra[nr]=ra[nr2]; ra[nr2]=temp;}}}
}

fftw_complex * makelrkernel(int rseed)
{
  /* this makes a random lrkernel */
  /* note that this is in row,col,ang space, which is 3 dimensional */
  int nr=0,nc=0,na=0,tab=0;
  double rdist=0,adist=0,weight=0;
  double wcutoff=0.0625;
  fftw_complex *lrk=NULL;
  int total_link=0;
  srand(rseed);
  lrk = (fftw_complex *) tcalloc(NPIEROWS*NPIECOLS*PIE_ROW_DIA*PIE_COL_DIA,sizeof(fftw_complex));
  total_link=0;
  for (nr=0;nr<NPIEROWS;nr++){ for (nc=0;nc<NPIECOLS;nc++){ for (na=0;na<PIE_ROW_DIA*PIE_COL_DIA;na++){
    tab = nr+nc*NPIEROWS+na*NPIEROWS*NPIECOLS;
    rdist = exp(-(pow(periodize(nr,-NPIEROWS/2,NPIEROWS/2),2)+pow(periodize(nc,-NPIECOLS/2,NPIECOLS/2),2))/2.0/LR_DIST_DIA_INPIES);
    adist = exp(-pow(periodize(na,-(PIE_ROW_DIA*PIE_COL_DIA)/2,(PIE_ROW_DIA*PIE_COL_DIA)/2)*PI/PIE_ROW_DIA/PIE_COL_DIA,2)/2.0/LR_ANGLE_DIA);
    weight = rdist*adist/* *(nr!=0 || nc!=0) */; /* no pie to pie autapses */
    if (weight>wcutoff){ lrk[tab][0] = weight; total_link++;}}}}
  for (nr=0;nr<NPIEROWS;nr++){ for (nc=0;nc<NPIECOLS;nc++){ for (na=0;na<PIE_ROW_DIA*PIE_COL_DIA;na++){
    tab = nr+nc*NPIEROWS+na*NPIEROWS*NPIECOLS;
    if (lrk[tab][0]>wcutoff){ lrk[tab][0] *= (rand01<LR_PCONNECT)/LR_PCONNECT;} /* sparsity of a sort */}}}
  return lrk;
}

void neuronsetang(struct neuronarray *Nra)
{
  /* this provides the ang index, the ordering of neurons within a pinwheel based on lgnangle */
  int ar=0,ac=0,nr2=0,nc2=0,nr=0,nc=0;
  struct neuron *n=NULL;
  struct llist *Ln=NULL;
  struct litem *l=NULL;
  int ang=0;
  for (ar=0;ar<NPIEROWS;ar++){ for (ac=0;ac<NPIECOLS;ac++){
    if (Ln!=NULL){ llisttfree(Ln); Ln=NULL;} Ln=llistmake();
    for (nr2=0;nr2<PIE_ROW_DIA;nr2++){ for (nc2=0;nc2<PIE_COL_DIA;nc2++){
      nr = ar*PIE_ROW_DIA+nr2;
      nc = ac*PIE_COL_DIA+nc2;
      n = nget(Nra,nr,nc);
      litemadd(Ln,n);}}
    llistsort(Ln->first,Ln->last,Ln->length,&lgnangle_compare);
    ang=0; l=Ln->first; while(l!=NULL){ n=(struct neuron *)l->item; n->ang=ang; ang++; l=l->child;}}}
  if (Ln!=NULL){ llisttfree(Ln); Ln=NULL;}
}

void neuronmake(struct neuronarray *Nra,int row,int col)
{
  struct neuron *n=NULL;
  int tab=0;
  double r=0,angle=0;
  double radmax = sqrt(pow(PIE_ROW_DIA/2.0,2)+pow(PIE_COL_DIA/2.0,2))+1;
  double *lvalue=NULL;
  n = (struct neuron *) tmalloc(sizeof(struct neuron));
  n->row=row;
  n->col=col;
  tab = n->row + n->col*Nra->rows;
  n->pierow = n->row/PIE_ROW_DIA;
  n->piecol = n->col/PIE_COL_DIA;
  n->inpierow = -(double)(n->row%PIE_ROW_DIA) + (double)PIE_ROW_DIA/2.0 - 1.0;
  n->inpiecol = (double)(n->col%PIE_COL_DIA - PIE_COL_DIA/2);
  n->rad = minimum(sqrt(pow(n->inpierow,2)+pow(n->inpiecol,2))/radmax,1);
  n->ang = 0;
  n->V = &(Nra->Vra[tab]);
  n->sA = &(Nra->sAra[tab]);
  n->sN = &(Nra->sNra[tab]);
  n->sG = &(Nra->sGra[tab]);
  n->VS = &(Nra->VSra[tab]);
  n->sV1=0;
  n->sV2=0;
  n->sA1=0;
  n->sA2=0;
  n->sN1=0;
  n->sN2=0;
  n->sG1=0;
  n->sG2=0;
  n->g=0;
  n->spiketime_guess=0;
  n->spiketime_guess_flag=0;
  n->inputrate=0;
  n->spikeinput_flag=0;
  n->spikeinput_time=0;
  n->spikeinput_multiplicity=0;
  n->spikeinput_rseed= GLOBAL_SPIKEINPUT_RSEED + n->row + n->col*Nra->rows + GLOBAL_RECORD_NUMBER*Nra->rows*Nra->cols;
  n->spikelast=0;
  n->spiketime=0;
  n->spikenext=0;
  n->sparse_out=0;
  n->sparse_in=0;
  n->homogenizing_knob=maximum(0,minimum(1,STARTING_HOMOGENIZATION));
  n->warning=0;
  n->t2s = &(Nra->t2sra[tab]);
  r = rand01;
  if (r<P_ES/(P_ES+P_EC+P_IS+P_IC)){ /* excitatory simple */ n->type=+1;n->sox=+1;lvalue=(double *)(n->t2s);*lvalue=3;}
  else if (r<(P_ES+P_EC)/(P_ES+P_EC+P_IS+P_IC)){ /* excitatory complex */ n->type=+1;n->sox=-1;lvalue=(double *)(n->t2s);*lvalue=2;}
  else if (r<(P_ES+P_EC+P_IS)/(P_ES+P_EC+P_IS+P_IC)){ /* inhibitory simple */ n->type=-1;n->sox=+1;lvalue=(double *)(n->t2s);*lvalue=1;}
  else /*if (r<(P_ES+P_EC+P_IS+P_IC)/(P_ES+P_EC+P_IS+P_IC))*/{ /* inhibitory complex */ n->type=-1;n->sox=-1;lvalue=(double *)(n->t2s);*lvalue=0;}
  if (n->pierow%2==0 && n->piecol%2==0){
    angle = atan2(n->inpierow,n->inpiecol)/2.0;
    n->lgnangle = angle;}
  else if (n->pierow%2==1 && n->piecol%2==0){
    angle = atan2(-n->inpierow,n->inpiecol)/2.0;
    n->lgnangle = angle;}
  else if (n->pierow%2==0 && n->piecol%2==1){
    angle = atan2(n->inpierow,-n->inpiecol)/2.0;
    n->lgnangle = angle;}
  else if (n->pierow%2==1 && n->piecol%2==1){
    angle = atan2(-n->inpierow,-n->inpiecol)/2.0;
    n->lgnangle = angle;}
  if (!finite(angle)){ angle=0;}
  if (angle > PI/2.0 || angle < -PI/2.0){ printf(" %% error, bad angle\n");}
  /* The above sets the lgnangle to be tiled continuously across a square lattice.
     But if we want to have slight columnar skips for iso-orientation symmetry, we can add: */
  n->lgnangle = periodize(n->lgnangle + -PI/2 + (double)(n->col/PIE_COL_DIA)/(double)NPIECOLS*PI,-PI/2,PI/2);
  n->lgnangle_drift = LGNANGLE_DRIFT*2*(rand01 - 0.5);
  n->lgndriftx = GABOR_DRIFT*2*PI*(rand01 - 0.5);
  n->lgndrifty = GABOR_DRIFT*2*PI*(rand01 - 0.5);
  n->lgnphase = rand01*2*PI;
  *(n->V) = VOLTAGE_REST + rand01*(VOLTAGE_THRESHOLD - VOLTAGE_REST);
  *(n->sG) = ((int)*(n->t2s)>1) ? 3.0*rand01*(CS_ICEC_G+CS_ICES_G+CS_ISEC_G+CS_ISES_G)/4.0 : 0;
  n->spikeinput_time=-2*TAU_REF;
  n->spikelast=-2*TAU_REF;
  n->spiketime=-2*TAU_REF;
  n->spikenext=-TAU_REF;
  n->sparse_out=rand();
  n->sparse_in=rand();
  n->sparse_link=NULL;
  n->lr_kernel_type=rand()%NLRKERNELS;
  n->sog=0;
  nset(Nra,n->row,n->col,n);
}

void neurontfree(struct neuron *n)
{
  if (n->sparse_link!=NULL){ llitemtfree(n->sparse_link,NULL); n->sparse_link=NULL;}
  tfree(n);
  n = NULL;
}

struct neuron * nget(struct neuronarray *Nra,int row,int col)
{
  if (row<0 || row>=Nra->rows || col<0 || col>=Nra->cols){ 
    printf(" %% Warning! attempting to retrieve row=%d,col=%d out of (%d,%d)\n",row,col,Nra->rows,Nra->cols); exit(EXIT_FAILURE);}
  assert(row>=0 && row<Nra->rows && col>=0 && col<Nra->cols);
  return Nra->N[row + col*Nra->rows];
}

void nset(struct neuronarray *Nra,int row,int col,void *n)
{
  if (n==NULL){ Nra->N[row + col*Nra->rows] = NULL;}
  else if (n!=NULL){ Nra->N[row + col*Nra->rows] = (struct neuron *) n;}
}

struct neuronarray * neuronarraymake(int rows,int cols)
{
  int verbose=0;
  int nr=0,nc=0,na=0,nr2=0,nc2=0,nr3=0,nc3=0,ntemp=0,nmax=0,nmin=0;
  struct neuronarray *Nra=NULL;
  struct neuron *n=NULL,*n2=NULL;
  int sparse_row_dia=SPARSE_ROW_DIA;
  int sparse_col_dia=SPARSE_COL_DIA;
  int tab=0,db=0;
  double *lvalue=NULL;
  Nra = (struct neuronarray *) tmalloc(sizeof(struct neuronarray));
  Nra->rows = rows;
  Nra->cols = cols;
  Nra->Vra = (double *) tcalloc(Nra->rows*Nra->cols,sizeof(double));
  Nra->sAra = (double *) tcalloc(Nra->rows*Nra->cols,sizeof(double));
  Nra->sNra = (double *) tcalloc(Nra->rows*Nra->cols,sizeof(double));
  Nra->sGra = (double *) tcalloc(Nra->rows*Nra->cols,sizeof(double));
  Nra->VSra = (double *) tcalloc(Nra->rows*Nra->cols,sizeof(double));
  Nra->t2sra = (double *) tcalloc(Nra->rows*Nra->cols,sizeof(double));
  Nra->N = (struct neuron **) tmalloc(sizeof(struct neuron *)*rows*cols);
  for (nr=0;nr<rows;nr++){ for (nc=0;nc<cols;nc++){ neuronmake(Nra,nr,nc);}}
  neuronsetang(Nra);
  Nra->t2stotal = (int *) tcalloc(4,sizeof(int));
  Nra->natotal = (int *) tcalloc(NSLICES,sizeof(int));
  for (nr=0;nr<rows;nr++){ for (nc=0;nc<cols;nc++){
    n = nget(Nra,nr,nc);
    na = periodize((int)floor(NSLICES*(n->lgnangle+PI/2)/PI),0,NSLICES);
    Nra->t2stotal[(int)*(n->t2s)] += 1;
    Nra->natotal[na] += 1;}}
  Nra->mES = 0; Nra->mEC = 0; Nra->mIS = 0; Nra->mIC = 0;
  if (BIG_SYSTEM_FLAG<0){
    for (nr=0;nr<rows;nr++){ for (nc=0;nc<cols;nc++){
      n=nget(Nra,nr,nc);
      n->sparse_link = llitemmake(); n->sparse_link->item=n;
      if (SPARSE_CONNECTION_TYPE==-1){ /* ring architecture, with rewiring if desired (P_REWIRE) */
	if (SPARSE_T2S_ORDERING>=0){ /* only ES and IS allowed */
	  db = (int)floor(pow(2,n->row+n->col*Nra->rows));
	  lvalue=(double *)(n->t2s); *lvalue = 1+2*((SPARSE_T2S_ORDERING%(2*db))/db);
	  n->type = ((int)(*lvalue))>=2 ? +1 : -1 ; n->sox = ((int)(*lvalue))%2 ? +1 : -1;}
	for (nr2=-sparse_row_dia/2;nr2<=sparse_row_dia/2;nr2++){ for (nc2=-sparse_col_dia/2;nc2<=sparse_col_dia/2;nc2++){
	  if (nr2!=0 || nc2!=0){
	    nr3 = periodize(nr+nr2,0,SYSTEM_ROW_SIZE);
	    nc3 = periodize(nc+nc2,0,SYSTEM_COL_SIZE);
	    if (rand01<P_REWIRE){
	      if (verbose){ printf("rewiring %d,%d,  %d,%d .... ",nr,nc,nr3,nc3);}
	      ntemp = nr + SYSTEM_ROW_SIZE/2 + (rand()%(SYSTEM_ROW_SIZE-sparse_row_dia-2) - SYSTEM_ROW_SIZE/2 + sparse_row_dia/2 + 1);
	      nmin = 0;
	      nmax = SYSTEM_ROW_SIZE;
	      periodify("int",&ntemp,&nmin,&nmax,&nr3);
	      ntemp = nc + SYSTEM_COL_SIZE/2 + (rand()%(SYSTEM_COL_SIZE-sparse_col_dia-2) - SYSTEM_COL_SIZE/2 + sparse_col_dia/2 + 1);
	      nmin = 0;
	      nmax = SYSTEM_COL_SIZE;
	      periodify("int",&ntemp,&nmin,&nmax,&nc3);
	      if (verbose){ printf(" to %d,%d\n",nr3,nc3);}}
	    n2 = nget(Nra,nr3,nc3);
	    if (llitemaddorfind(0,n->sparse_link,n2,&void_compare)==NULL){ llitemaddorfind(1,n->sparse_link,n2,&void_compare);}}}}}
      else if (SPARSE_CONNECTION_TYPE==-2){ /* collaborative effect, changes Nra->t2sra, intended for AUTAPSES_OFF, SYSTEM_SIZE>2 */
	tab = n->row + n->col*Nra->rows;
	if (tab==0){ /* simple cell */ 
	  *(n->t2s)=3; n->type = ((int)(*(n->t2s)))>=2 ? +1 : -1 ; n->sox = ((int)(*(n->t2s)))%2 ? +1 : -1;}
	else if (tab==1){ /* complex inhibitory cell "suppressor" */ 
	  *(n->t2s)=0; n->type = ((int)(*(n->t2s)))>=2 ? +1 : -1 ; n->sox = ((int)(*(n->t2s)))%2 ? +1 : -1;
	  n2=Nra->N[0];if(llitemaddorfind(0,n->sparse_link,n2,&void_compare)==NULL){llitemaddorfind(1,n->sparse_link,n2,&void_compare);}}
	else if (tab>1){ /* complex excitatory engine */
	  *(n->t2s)=2; n->type = ((int)(*(n->t2s)))>=2 ? +1 : -1 ; n->sox = ((int)(*(n->t2s)))%2 ? +1 : -1;
	  n2=Nra->N[0];if(llitemaddorfind(0,n->sparse_link,n2,&void_compare)==NULL){llitemaddorfind(1,n->sparse_link,n2,&void_compare);}}}
      else if (SPARSE_CONNECTION_TYPE==-3){ /* engine-killer effect, changes Nra->t2sra, intended for AUTAPSES_OFF, SYSTEM_SIZE>4 */
	tab = n->row + n->col*Nra->rows;
	if (tab==0){ /* simple cell */ 
	  *(n->t2s)=3; n->type = ((int)(*(n->t2s)))>=2 ? +1 : -1 ; n->sox = ((int)(*(n->t2s)))%2 ? +1 : -1;}
	else if (tab==1){ /* complex inhibitory cell "feeder and linker" */
	  *(n->t2s)=0; n->type = ((int)(*(n->t2s)))>=2 ? +1 : -1 ; n->sox = ((int)(*(n->t2s)))%2 ? +1 : -1;
	  n2=Nra->N[0];if(llitemaddorfind(0,n->sparse_link,n2,&void_compare)==NULL){llitemaddorfind(1,n->sparse_link,n2,&void_compare);}
	  n2=Nra->N[3];if(llitemaddorfind(0,n->sparse_link,n2,&void_compare)==NULL){llitemaddorfind(1,n->sparse_link,n2,&void_compare);}}
	else if (tab==2){ /* complex inhibitory cell "linker" */
	  *(n->t2s)=0; n->type = ((int)(*(n->t2s)))>=2 ? +1 : -1 ; n->sox = ((int)(*(n->t2s)))%2 ? +1 : -1;
	  n2=Nra->N[1];if(llitemaddorfind(0,n->sparse_link,n2,&void_compare)==NULL){llitemaddorfind(1,n->sparse_link,n2,&void_compare);}}
	else if (tab==3){ /* complex excitatory cell "feeder" */
	  *(n->t2s)=2; n->type = ((int)(*(n->t2s)))>=2 ? +1 : -1 ; n->sox = ((int)(*(n->t2s)))%2 ? +1 : -1;
	  n2=Nra->N[0];if(llitemaddorfind(0,n->sparse_link,n2,&void_compare)==NULL){llitemaddorfind(1,n->sparse_link,n2,&void_compare);}
	  n2=Nra->N[2];if(llitemaddorfind(0,n->sparse_link,n2,&void_compare)==NULL){llitemaddorfind(1,n->sparse_link,n2,&void_compare);}}
	else if (tab>3){ /* complex excitatory engine */
	  *(n->t2s)=2; n->type = ((int)(*(n->t2s)))>=2 ? +1 : -1 ; n->sox = ((int)(*(n->t2s)))%2 ? +1 : -1;
	  n2=Nra->N[1];if(llitemaddorfind(0,n->sparse_link,n2,&void_compare)==NULL){llitemaddorfind(1,n->sparse_link,n2,&void_compare);}
	  n2=Nra->N[3];if(llitemaddorfind(0,n->sparse_link,n2,&void_compare)==NULL){llitemaddorfind(1,n->sparse_link,n2,&void_compare);}}}
      else if (SPARSE_CONNECTION_TYPE==-4){ /* symmetric engine effect, changes Nra->t2sra, intended for AUTAPSES_OFF, SYSTEM_SIZE>7 */
	tab = n->row + n->col*Nra->rows;
	if (tab==0){ /* simple cell */ 
	  *(n->t2s)=3; n->type = ((int)(*(n->t2s)))>=2 ? +1 : -1 ; n->sox = ((int)(*(n->t2s)))%2 ? +1 : -1;}
	else if (tab==1){ /* complex inhibitory cell "feeder" */
	  *(n->t2s)=0; n->type = ((int)(*(n->t2s)))>=2 ? +1 : -1 ; n->sox = ((int)(*(n->t2s)))%2 ? +1 : -1;
	  n2=Nra->N[0];if(llitemaddorfind(0,n->sparse_link,n2,&void_compare)==NULL){llitemaddorfind(1,n->sparse_link,n2,&void_compare);}}
	else if (tab==2){ /* complex excitatory cell "feeder" */
	  *(n->t2s)=2; n->type = ((int)(*(n->t2s)))>=2 ? +1 : -1 ; n->sox = ((int)(*(n->t2s)))%2 ? +1 : -1;
	  n2=Nra->N[0];if(llitemaddorfind(0,n->sparse_link,n2,&void_compare)==NULL){llitemaddorfind(1,n->sparse_link,n2,&void_compare);}}
	else if (tab==3){ /* complex excitatory cell "feeder" */
	  *(n->t2s)=2; n->type = ((int)(*(n->t2s)))>=2 ? +1 : -1 ; n->sox = ((int)(*(n->t2s)))%2 ? +1 : -1;
	  n2=Nra->N[1];if(llitemaddorfind(0,n->sparse_link,n2,&void_compare)==NULL){llitemaddorfind(1,n->sparse_link,n2,&void_compare);}
	  n2=Nra->N[6];if(llitemaddorfind(0,n->sparse_link,n2,&void_compare)==NULL){llitemaddorfind(1,n->sparse_link,n2,&void_compare);}}
	else if (tab==4){ /* complex excitatory cell "feeder" */
	  *(n->t2s)=2; n->type = ((int)(*(n->t2s)))>=2 ? +1 : -1 ; n->sox = ((int)(*(n->t2s)))%2 ? +1 : -1;
	  n2=Nra->N[2];if(llitemaddorfind(0,n->sparse_link,n2,&void_compare)==NULL){llitemaddorfind(1,n->sparse_link,n2,&void_compare);}
	  n2=Nra->N[5];if(llitemaddorfind(0,n->sparse_link,n2,&void_compare)==NULL){llitemaddorfind(1,n->sparse_link,n2,&void_compare);}}
	else if (tab==5){ /* complex inhibitory cell "linker" */
	  *(n->t2s)=0; n->type = ((int)(*(n->t2s)))>=2 ? +1 : -1 ; n->sox = ((int)(*(n->t2s)))%2 ? +1 : -1;
	  n2=Nra->N[3];if(llitemaddorfind(0,n->sparse_link,n2,&void_compare)==NULL){llitemaddorfind(1,n->sparse_link,n2,&void_compare);}}
	else if (tab==6){ /* complex inhibitory cell "linker" */
	  *(n->t2s)=0; n->type = ((int)(*(n->t2s)))>=2 ? +1 : -1 ; n->sox = ((int)(*(n->t2s)))%2 ? +1 : -1;
	  n2=Nra->N[4];if(llitemaddorfind(0,n->sparse_link,n2,&void_compare)==NULL){llitemaddorfind(1,n->sparse_link,n2,&void_compare);}}
	else if (tab>6){ /* complex excitatory engine */
	  *(n->t2s)=2; n->type = ((int)(*(n->t2s)))>=2 ? +1 : -1 ; n->sox = ((int)(*(n->t2s)))%2 ? +1 : -1;
	  n2=Nra->N[3];if(llitemaddorfind(0,n->sparse_link,n2,&void_compare)==NULL){llitemaddorfind(1,n->sparse_link,n2,&void_compare);}
	  n2=Nra->N[4];if(llitemaddorfind(0,n->sparse_link,n2,&void_compare)==NULL){llitemaddorfind(1,n->sparse_link,n2,&void_compare);}}}
      else if (SPARSE_CONNECTION_TYPE==-5){ /* discriminable? changes Nra->t2sra, intended for AUTAPSES_OFF==0, SYSTEM_SIZE==8 */
	tab = n->row + n->col*Nra->rows;
	if (tab==0){ /* simple excitatory */
	  *(n->t2s)=3; n->type = ((int)(*(n->t2s)))>=2 ? +1 : -1 ; n->sox = ((int)(*(n->t2s)))%2 ? +1 : -1;
	  n2=Nra->N[1];if(llitemaddorfind(0,n->sparse_link,n2,&void_compare)==NULL){llitemaddorfind(1,n->sparse_link,n2,&void_compare);}
	  n2=Nra->N[3];if(llitemaddorfind(0,n->sparse_link,n2,&void_compare)==NULL){llitemaddorfind(1,n->sparse_link,n2,&void_compare);}
	  n2=Nra->N[6];if(llitemaddorfind(0,n->sparse_link,n2,&void_compare)==NULL){llitemaddorfind(1,n->sparse_link,n2,&void_compare);}
	  n2=Nra->N[7];if(llitemaddorfind(0,n->sparse_link,n2,&void_compare)==NULL){llitemaddorfind(1,n->sparse_link,n2,&void_compare);}}
	else if (tab==1){ /* simple excitatory */
	  *(n->t2s)=3; n->type = ((int)(*(n->t2s)))>=2 ? +1 : -1 ; n->sox = ((int)(*(n->t2s)))%2 ? +1 : -1;
	  n2=Nra->N[0];if(llitemaddorfind(0,n->sparse_link,n2,&void_compare)==NULL){llitemaddorfind(1,n->sparse_link,n2,&void_compare);}
	  n2=Nra->N[3];if(llitemaddorfind(0,n->sparse_link,n2,&void_compare)==NULL){llitemaddorfind(1,n->sparse_link,n2,&void_compare);}
	  n2=Nra->N[4];if(llitemaddorfind(0,n->sparse_link,n2,&void_compare)==NULL){llitemaddorfind(1,n->sparse_link,n2,&void_compare);}
	  n2=Nra->N[5];if(llitemaddorfind(0,n->sparse_link,n2,&void_compare)==NULL){llitemaddorfind(1,n->sparse_link,n2,&void_compare);}}
	else if (tab==2){ /* simple excitatory */
	  *(n->t2s)=3; n->type = ((int)(*(n->t2s)))>=2 ? +1 : -1 ; n->sox = ((int)(*(n->t2s)))%2 ? +1 : -1;
	  n2=Nra->N[6];if(llitemaddorfind(0,n->sparse_link,n2,&void_compare)==NULL){llitemaddorfind(1,n->sparse_link,n2,&void_compare);}}
	else if (tab==3){ /* simple excitatory */
	  *(n->t2s)=3; n->type = ((int)(*(n->t2s)))>=2 ? +1 : -1 ; n->sox = ((int)(*(n->t2s)))%2 ? +1 : -1;
	  n2=Nra->N[2];if(llitemaddorfind(0,n->sparse_link,n2,&void_compare)==NULL){llitemaddorfind(1,n->sparse_link,n2,&void_compare);}
	  n2=Nra->N[4];if(llitemaddorfind(0,n->sparse_link,n2,&void_compare)==NULL){llitemaddorfind(1,n->sparse_link,n2,&void_compare);}
	  n2=Nra->N[7];if(llitemaddorfind(0,n->sparse_link,n2,&void_compare)==NULL){llitemaddorfind(1,n->sparse_link,n2,&void_compare);}}
	else if (tab==4){ /* simple excitatory */
	  *(n->t2s)=3; n->type = ((int)(*(n->t2s)))>=2 ? +1 : -1 ; n->sox = ((int)(*(n->t2s)))%2 ? +1 : -1;
	  n2=Nra->N[0];if(llitemaddorfind(0,n->sparse_link,n2,&void_compare)==NULL){llitemaddorfind(1,n->sparse_link,n2,&void_compare);}
	  n2=Nra->N[3];if(llitemaddorfind(0,n->sparse_link,n2,&void_compare)==NULL){llitemaddorfind(1,n->sparse_link,n2,&void_compare);}
	  n2=Nra->N[7];if(llitemaddorfind(0,n->sparse_link,n2,&void_compare)==NULL){llitemaddorfind(1,n->sparse_link,n2,&void_compare);}}
	else if (tab==5){ /* simple inhibitory */
	  *(n->t2s)=1; n->type = ((int)(*(n->t2s)))>=2 ? +1 : -1 ; n->sox = ((int)(*(n->t2s)))%2 ? +1 : -1;
	  n2=Nra->N[1];if(llitemaddorfind(0,n->sparse_link,n2,&void_compare)==NULL){llitemaddorfind(1,n->sparse_link,n2,&void_compare);}
	  n2=Nra->N[4];if(llitemaddorfind(0,n->sparse_link,n2,&void_compare)==NULL){llitemaddorfind(1,n->sparse_link,n2,&void_compare);}
	  n2=Nra->N[7];if(llitemaddorfind(0,n->sparse_link,n2,&void_compare)==NULL){llitemaddorfind(1,n->sparse_link,n2,&void_compare);}}
	else if (tab==6){ /* simple inhibitory */
	  *(n->t2s)=1; n->type = ((int)(*(n->t2s)))>=2 ? +1 : -1 ; n->sox = ((int)(*(n->t2s)))%2 ? +1 : -1;
	  n2=Nra->N[0];if(llitemaddorfind(0,n->sparse_link,n2,&void_compare)==NULL){llitemaddorfind(1,n->sparse_link,n2,&void_compare);}
	  n2=Nra->N[1];if(llitemaddorfind(0,n->sparse_link,n2,&void_compare)==NULL){llitemaddorfind(1,n->sparse_link,n2,&void_compare);}
	  n2=Nra->N[2];if(llitemaddorfind(0,n->sparse_link,n2,&void_compare)==NULL){llitemaddorfind(1,n->sparse_link,n2,&void_compare);}}
	else if (tab==7){ /* simple excitatory */
	  *(n->t2s)=3; n->type = ((int)(*(n->t2s)))>=2 ? +1 : -1 ; n->sox = ((int)(*(n->t2s)))%2 ? +1 : -1;
	  n2=Nra->N[0];if(llitemaddorfind(0,n->sparse_link,n2,&void_compare)==NULL){llitemaddorfind(1,n->sparse_link,n2,&void_compare);}
	  n2=Nra->N[1];if(llitemaddorfind(0,n->sparse_link,n2,&void_compare)==NULL){llitemaddorfind(1,n->sparse_link,n2,&void_compare);}
	  n2=Nra->N[3];if(llitemaddorfind(0,n->sparse_link,n2,&void_compare)==NULL){llitemaddorfind(1,n->sparse_link,n2,&void_compare);}
	  n2=Nra->N[5];if(llitemaddorfind(0,n->sparse_link,n2,&void_compare)==NULL){llitemaddorfind(1,n->sparse_link,n2,&void_compare);}}}
      else if (SPARSE_CONNECTION_TYPE==-6){ /* differentiator? changes Nra->t2sra, intended for AUTAPSES_OFF==0, SYSTEM_SIZE==8 */
	tab = n->row + n->col*Nra->rows;
	if (tab==0){ /* complex inhibitory */
	  *(n->t2s)=0; n->type = ((int)(*(n->t2s)))>=2 ? +1 : -1 ; n->sox = ((int)(*(n->t2s)))%2 ? +1 : -1;
	  n2=Nra->N[1];if(llitemaddorfind(0,n->sparse_link,n2,&void_compare)==NULL){llitemaddorfind(1,n->sparse_link,n2,&void_compare);}}
	else if (tab==1){ /* complex inhibitory */
	  *(n->t2s)=0; n->type = ((int)(*(n->t2s)))>=2 ? +1 : -1 ; n->sox = ((int)(*(n->t2s)))%2 ? +1 : -1;
	  n2=Nra->N[2];if(llitemaddorfind(0,n->sparse_link,n2,&void_compare)==NULL){llitemaddorfind(1,n->sparse_link,n2,&void_compare);}}
	else if (tab==2){ /* simple inhibitory */
	  *(n->t2s)=1; n->type = ((int)(*(n->t2s)))>=2 ? +1 : -1 ; n->sox = ((int)(*(n->t2s)))%2 ? +1 : -1;
	  n2=Nra->N[0];if(llitemaddorfind(0,n->sparse_link,n2,&void_compare)==NULL){llitemaddorfind(1,n->sparse_link,n2,&void_compare);}}
	else if (tab==3){ /* simple excitatory */
	  *(n->t2s)=3; n->type = ((int)(*(n->t2s)))>=2 ? +1 : -1 ; n->sox = ((int)(*(n->t2s)))%2 ? +1 : -1;
	  n2=Nra->N[0];if(llitemaddorfind(0,n->sparse_link,n2,&void_compare)==NULL){llitemaddorfind(1,n->sparse_link,n2,&void_compare);}}
	else if (tab==4){ /* simple excitatory */
	  *(n->t2s)=3; n->type = ((int)(*(n->t2s)))>=2 ? +1 : -1 ; n->sox = ((int)(*(n->t2s)))%2 ? +1 : -1;
	  n2=Nra->N[1];if(llitemaddorfind(0,n->sparse_link,n2,&void_compare)==NULL){llitemaddorfind(1,n->sparse_link,n2,&void_compare);}}
	else if (tab==5){ /* simple excitatory */
	  *(n->t2s)=3; n->type = ((int)(*(n->t2s)))>=2 ? +1 : -1 ; n->sox = ((int)(*(n->t2s)))%2 ? +1 : -1;
	  n2=Nra->N[2];if(llitemaddorfind(0,n->sparse_link,n2,&void_compare)==NULL){llitemaddorfind(1,n->sparse_link,n2,&void_compare);}}
	else if (tab==6){ /* simple excitatory */
	  *(n->t2s)=3; n->type = ((int)(*(n->t2s)))>=2 ? +1 : -1 ; n->sox = ((int)(*(n->t2s)))%2 ? +1 : -1;
	  n2=Nra->N[0];if(llitemaddorfind(0,n->sparse_link,n2,&void_compare)==NULL){llitemaddorfind(1,n->sparse_link,n2,&void_compare);}
	  n2=Nra->N[1];if(llitemaddorfind(0,n->sparse_link,n2,&void_compare)==NULL){llitemaddorfind(1,n->sparse_link,n2,&void_compare);}
	  n2=Nra->N[2];if(llitemaddorfind(0,n->sparse_link,n2,&void_compare)==NULL){llitemaddorfind(1,n->sparse_link,n2,&void_compare);}}
	else if (tab==7){ /* simple excitatory */
	  *(n->t2s)=3; n->type = ((int)(*(n->t2s)))>=2 ? +1 : -1 ; n->sox = ((int)(*(n->t2s)))%2 ? +1 : -1;
	  n2=Nra->N[0];if(llitemaddorfind(0,n->sparse_link,n2,&void_compare)==NULL){llitemaddorfind(1,n->sparse_link,n2,&void_compare);}
	  n2=Nra->N[1];if(llitemaddorfind(0,n->sparse_link,n2,&void_compare)==NULL){llitemaddorfind(1,n->sparse_link,n2,&void_compare);}
	  n2=Nra->N[2];if(llitemaddorfind(0,n->sparse_link,n2,&void_compare)==NULL){llitemaddorfind(1,n->sparse_link,n2,&void_compare);}}}
      llitembalance(n->sparse_link); n->sparse_link=llitemclimb(n->sparse_link); llitemcheck(0,n->sparse_link,&void_compare);}}}
  return Nra;
}

void neuronarraytfree(struct neuronarray *Nra)
{
  int nr=0,nc=0;
  for (nr=0;nr<Nra->rows;nr++){ for (nc=0;nc<Nra->cols;nc++){
    neurontfree(nget(Nra,nr,nc));
    nset(Nra,nr,nc,NULL);}}
  tfree(Nra->N);Nra->N=NULL;
  tfree(Nra->t2sra);Nra->t2sra=NULL;
  tfree(Nra);Nra=NULL;    
}

void setblocksize(double DT)
{ 
  GLOBAL_BLOCK_ROW_DIA = GLOBAL_BLOCK_ROW_DIA_MIN;
  GLOBAL_BLOCK_COL_DIA = GLOBAL_BLOCK_COL_DIA_MIN;
  GLOBAL_NBLOCKS_TALL=SYSTEM_ROW_SIZE/GLOBAL_BLOCK_ROW_DIA;
  GLOBAL_NBLOCKS_WIDE=SYSTEM_COL_SIZE/GLOBAL_BLOCK_COL_DIA;
}

void setinputrate(struct neuronarray *Nra, double t)
{
  /* sets input spiketimes */
  int nr=0,nc=0,nr2=0,nc2=0,na=0,np=0;
  double rateplus=0;
  struct lgn *p = NULL;
  struct neuron *n=NULL;
  int nang = PIE_ROW_DIA*PIE_COL_DIA;
  int local_odor=0;
  if (LGN_BOTHER){
    p = GLOBAL_LGN;
    for (nr=0;nr<Nra->rows;nr++){ for (nc=0;nc<Nra->cols;nc++){
      n = nget(Nra,nr,nc);
      if (n->sox == +1){ 
	na = periodize((int)floor(NSLICES*(n->lgnangle+PI/2)/PI),0,NSLICES);
	np = periodize((int)floor(NPHASES*n->lgnphase/PI),0,2*NPHASES);
	nr2 = (int)floor((double)nr/(double)Nra->rows*p->rows);
	nc2 = (int)floor((double)nc/(double)Nra->cols*p->cols);
	/* fix later */
	if (LGN_DUMB==0){ 
	  rateplus = p->anglerate[na+(np%NPHASES)*NSLICES][4+nr2*5+nc2*5*p->rows];
	  rateplus *= pow(-1,np/NPHASES);}
	else if (LGN_DUMB==1){
	  rateplus = p->peakrate*p->anglein[na+(np%NPHASES)*NSLICES][nr2+nc2*p->rows];}
	else if (LGN_DUMB==2){
	  rateplus = p->peakrate*exp(-pow(periodize(INPUT_SPACEANGLE/PI-(double)n->ang/nang,0,1)*4,2));}
	else if (LGN_DUMB==3){ /* global pulse */
	  rateplus = p->peakrate*cos(t*2*PI*GRATING_PULSE/1024.0);}
	else if (LGN_DUMB==4){ 
	  rateplus = 0;}
	rateplus *= INPUT_CONTRAST;
	n->inputrate = LGN_CELLS_PER_V1_CELL*maximum(0,LGN_BACKRATE + rateplus);}
      else if (n->sox == -1){
	n->inputrate = OTHERLAYER_BACKRATE + OTHERLAYER_INPUTRATE*INPUT_CONTRAST;}}}}
  else{ 
    for (nr=0;nr<Nra->rows;nr++){ for (nc=0;nc<Nra->cols;nc++){
      n = nget(Nra,nr,nc);
      if (n->sox==+1){ 
	if (LGN_DUMB>=0){
	  switch (LGN_DUMB){
	  case 0: /* first harmonic */
	    rateplus = 0; break;
	  case 1: /* global pulse */ 
	    rateplus = cos(t*2*PI*GRATING_PULSE/1024.0); break;
	  case 2: /* phase change */ 
	    rateplus = cos(t*2*PI*GRATING_PULSE/1024.0) + cos(t*2*PI*2*GRATING_PULSE/1024.0 + INPUT_SPACEPHASE); break;
	  case 3: /* phase change 2 */
	    rateplus = cos(t*2*PI*GRATING_PULSE/1024.0) + cos(t*2*PI*2*GRATING_PULSE/1024.0 + INPUT_SPACEPHASE + INPUT_SPACEANGLE_BACON); break;
	  case 4: /* square wave */ 
	    rateplus = cos(t*2*PI*GRATING_PULSE/1024.0)>0?+1:-1; break;
	  default: rateplus=0; break;}}
	else /* if (LGN_DUMB<0) */{
	  local_odor = (GLOBAL_ODOR%(long int)powl(GLOBAL_ODOR_BASE,n->row+Nra->rows*n->col+1));
	  local_odor /= powl(GLOBAL_ODOR_BASE,n->row+Nra->rows*n->col);
	  rateplus = (double)periodize(local_odor,0,GLOBAL_ODOR_BASE)/(double)(GLOBAL_ODOR_BASE-1);}
	rateplus *= INPUT_CONTRAST;
	n->inputrate = LGN_CELLS_PER_V1_CELL*maximum(0,LGN_BACKRATE+rateplus);}
      else if (n->sox==-1){ n->inputrate = OTHERLAYER_BACKRATE + OTHERLAYER_INPUTRATE*INPUT_CONTRAST;}}}}
}

void spikeinput(struct neuronarray *Nra, double t, double DT)
{
  /* sets input spiketimes */
  int verbose=0;
  int nr=0,nc=0;
  struct neuron *n=NULL;
  double d=0,threshold_time=10;
  if (LGN_DETERMINISM==0){
    for (nr=0;nr<Nra->rows;nr++){ for (nc=0;nc<Nra->cols;nc++){
      n = nget(Nra,nr,nc);
      if (n->inputrate>0){ 
	if (n->spikeinput_time < t){
	  if (nr==0 && nc==0 && verbose){ printf("t=%0.2f DT=%0.2f n(%d,%d)->spikeinput_time=%0.2f, inputrate=%0.2f\n",t,DT,nr,nc,n->spikeinput_time,n->inputrate);}
	  n->spikeinput_multiplicity=0; d=0;
	  while (n->spikeinput_time < t){ 
	    d = RISIGET(&(n->spikeinput_rseed),n->inputrate); 
	    if (d<threshold_time){ n->spikeinput_time +=d; n->spikeinput_multiplicity += 1;}
	    else{ n->spikeinput_time += 1;}
	    if (nr==0 && nc==0 && verbose){ printf("\tinpttm=%0.2f, mltplcty=%d\n",n->spikeinput_time,n->spikeinput_multiplicity);}}}
	if (n->spikeinput_time >= t){
	  if (n->spikeinput_time <= t+DT && n->spikeinput_multiplicity>0){ n->spikeinput_flag=1;}
	  else{ n->spikeinput_flag=0;}}}
      else{ n->spikeinput_flag=0;}}}}
  else if (LGN_DETERMINISM==1){
    for (nr=0;nr<Nra->rows;nr++){ for (nc=0;nc<Nra->cols;nc++){
      n = nget(Nra,nr,nc);
      if (nr==0 && nc==0 && verbose){ printf("t=%0.2f DT=%0.2f n(%d,%d)->spikeinput_time=%0.2f\n",t,DT,nr,nc,n->spikeinput_time);}
      if (n->spikeinput_time >= t){
	if (nr==0 && nc==0 && verbose){ printf("\tspikeinput_time >= t... ");}
	if (n->spikeinput_time <= t+DT){ 
	  if (nr==0 && nc==0 && verbose){ printf("but spikeinput_time <= t+DT, spikeinput_flag set to 1\n");}
	  n->spikeinput_flag=1;} 
	else{ 
	  if (nr==0 && nc==0 && verbose){ printf("but spikeinput_time > t+DT, spikeinput_flag set to 0\n");}
	  n->spikeinput_flag=0;}}
      else if (n->spikeinput_time < t){ 
	if (nr==0 && nc==0 && verbose){ printf("\tspikeinput_time < t... ");}
	do{
	  if (nr==0 && nc==0 && verbose){ printf("changing from %0.2f ",n->spikeinput_time);}
	  n->spikeinput_time = ceil(n->spikeinput_time)+ceil(1.0/(n->inputrate>0 ? n->inputrate : 1.0))+(double)(0.5+n->row+n->col*Nra->rows)/(double)(Nra->rows*Nra->cols+1);
	  if (nr==0 && nc==0 && verbose){ printf("to %0.2f\n",n->spikeinput_time);}} 
	while (n->spikeinput_time<t); 
	n->spikeinput_multiplicity = 1;
	if (n->spikeinput_time>=t && n->spikeinput_time <= t+DT){ n->spikeinput_flag=1;} else{ n->spikeinput_flag=0;}}}}}
}

int spikescan(struct neuronarray *Nra, struct llist **Sra, struct llist **Gra, double t, double DT)
{
  int nr=0,nc=0,lr=0,lc=0,pspikes=0;
  double spiketime_guess=0;
  struct neuron *n=NULL;
  for (nr=0;nr<Nra->rows;nr++){ for (nc=0;nc<Nra->cols;nc++){
    lr = nr/GLOBAL_BLOCK_ROW_DIA;
    lc = nc/GLOBAL_BLOCK_COL_DIA;
    n = nget(Nra,nr,nc);
    if (EIFvsIF>0){
      spiketime_guess = spikeguesseif(n,t,DT);}
    else{
      spiketime_guess = spikeguess2(n,t,DT);}
    if (spiketime_guess >= t && spiketime_guess <= t+DT){ n->sog=1; litemadd(Sra[lr + lc*GLOBAL_NBLOCKS_TALL],n); pspikes++;}
    else if (spiketime_guess < t || spiketime_guess > t+DT){ n->sog=-1; litemadd(Gra[lr + lc*GLOBAL_NBLOCKS_TALL],n);}}}
  return pspikes;
}

double spikeguesseif(struct neuron *n,double t,double DT)
{
  double overestimate=2.0,underestimate=1.0;
  double sA=*(n->sA)*overestimate,sN=*(n->sN)*overestimate,sG=*(n->sG)*underestimate,V=*(n->V);
  double gt=0,vsgt=0,V2=0,k1=0;
  double d=0;
  double af=0,nf=0,gf=0;
  if (n->spikenext >= t+DT){ n->spiketime = t+2*DT;}
  else{ 
    d = t+DT-maximum(t,n->spikenext);
    sA *= exp(-0*d/TAU_AMPA); sN *= exp(-0*d/TAU_NMDA); sG *= exp(-1*d/TAU_GABA);
    if (n->spikeinput_flag){ sA += LGN_STRENGTH;}
    af = 1.0/P_AMPA;
    nf = 1.0/P_NMDA;
    gf = 1.0/P_GABA;
    switch ((int)*(n->t2s)){
    case 0: /* inhibitory complex */ sA += (CS_ESIC_A + CS_ECIC_A)/2.0; sN += (CS_ESIC_N + CS_ECIC_N); sG += 0; break;
    case 1: /* inhibitory simple  */ sA += (CS_ESIS_A + CS_ECIS_A)/2.0; sN += (CS_ESIS_N + CS_ECIS_N); sG += 0; break;
    case 2: /* excitatory complex */ sA += (CS_ESEC_A + CS_ECEC_A)/2.0; sN += (CS_ESEC_N + CS_ECEC_N); sG += 0; break;
    case 3: /* excitatory simple  */ sA += (CS_ESES_A + CS_ECES_A)/2.0; sN += (CS_ESES_N + CS_ECES_N); sG += 0; break;
    default: printf(" %% error, wrong type %d in spikeguesseif\n",(int)*(n->t2s)); break;}
    k1 = eifrhs(V,sA,sN,sG,0);
    gt = CONDUCTANCE_LK+sA+sN+sG;
    vsgt = (CONDUCTANCE_LK*VOLTAGE_REST + (sA+sN)*VOLTAGE_EX + sG*VOLTAGE_IN);
    V2 = (V + d*(vsgt+CONDUCTANCE_LK*VOLTAGE_DELTAT*exp((V-VOLTAGE_TAKEOFF)/VOLTAGE_DELTAT)))/(1+d*gt);
    if (V2>VOLTAGE_THRESHOLD_EIF){ n->spiketime = maximum(t,n->spikenext) + quadrootfinder(V,k1,V2,VOLTAGE_THRESHOLD_EIF,d);}
    else{ n->spiketime = t+2*DT;}}
  return n->spiketime;
}

/* double spikeguesseif(struct neuron *n,double t,double DT) */
/* { */
/*   /\* cheat and use clumpcorrect *\/ */
/*   struct llist *L1=llistmake(),*L2=llistmake(); */
/*   litemadd(L1,n); */
/*   clumpcorrect(L1,L2,t,DT); */
/*   n->spiketime = n->spiketime_guess_flag ? n->spiketime_guess : t+2*DT; */
/*   llisttfree(L1);llisttfree(L2); */
/*   //if (n->spiketime_guess_flag && n->spiketime_guess >=t && n->spiketime_guess <= t+DT){ printf("spike at %f\n",n->spiketime_guess);} */
/*   return n->spiketime; */
/* } */

double spikeguess2(struct neuron *n, double t, double DT)
{
  /* uses second order method to approximate spike time */
  double d=0,sA1=0,sN1=0,sG1=0,sA2=0,sN2=0,sG2=0,V1=0,V2=0,k1=0,eA=0,eN=0,eG=0,if1=0,if2=0,if3=0,vs1=0,vsp1=0,vs2=0,vsp2=0,vs3=0,vsp3=0,g=0;
  if (HANSHELLEY_FLAG){
    /* do nothing */
    n->spiketime = t+2*DT;
    return n->spiketime;}
  else /* if (!HANSHELLEY_FLAG) */{
    sA1=*(n->sA);sN1=*(n->sN);sG1=*(n->sG);V1=*(n->V);V2=*(n->V); if (n->spikeinput_flag){ sA1 += LGN_STRENGTH;}
    if (n->spikenext > t && *(n->V)!=VOLTAGE_RESET){ printf(" %% neuron (%d,%d) fired, but has V %f > VR in spikeguess2\n",n->row,n->col,*(n->V));}
    if (n->spikenext>t+DT){
      d = DT;
      V2 = VOLTAGE_RESET;}
    else /* if (n->spikenext<=t+DT) */{
      if (n->spikenext >= t && n->spikenext < t+DT){
	d = n->spikenext - t;
	sA1 *= exp(-d/TAU_AMPA);
	sN1 *= exp(-d/TAU_NMDA);
	sG1 *= exp(-d/TAU_GABA);}
      d = t+DT-maximum(t,n->spikenext);
      k1 = (CONDUCTANCE_LK+sA1+sN1+sG1)*V2 + (VOLTAGE_REST*CONDUCTANCE_LK + VOLTAGE_EX*(n->sA1+n->sN1) + VOLTAGE_IN*(n->sG1));
      eA = exp(-d/2/TAU_AMPA); eN = exp(-d/2/TAU_NMDA); eG = exp(-d/2/TAU_GABA);
      if1 = 0*CONDUCTANCE_LK+sA2+sN2+sG2;
      slavecalc(sA1,sN1,sG1,NULL,&vs1,&vsp1);
      sA2 = TAU_AMPA*sA1*(1-eA);
      sN2 = TAU_NMDA*sN1*(1-eN);
      sG2 = TAU_GABA*sG1*(1-eG);
      sA1 *= eA;
      sN1 *= eN;
      sG1 *= eG;
      if2 = d/2*CONDUCTANCE_LK+sA2+sN2+sG2;
      slavecalc(sA1,sN1,sG1,NULL,&vs2,&vsp2);
      sA2 += TAU_AMPA*sA1*(1-eA);
      sN2 += TAU_NMDA*sN1*(1-eN);
      sG2 += TAU_GABA*sG1*(1-eG);
      sA1 *= eA;
      sN1 *= eN;
      sG1 *= eG;
      if3 = d*CONDUCTANCE_LK+sA2+sN2+sG2;
      slavecalc(sA1,sN1,sG1,NULL,&vs3,&vsp3);
      g = d/6.0*(exp(if1-if3)*vsp1 + 4*exp(if2-if3)*vsp2 + exp(if3-if3)*vsp3);
      V2 = ((V1-vs1)*exp(-if3) - g) + vs3;}
    if (V2>VOLTAGE_THRESHOLD){
      n->spiketime = maximum(t,n->spikenext) + quadrootfinder(V1,k1,V2,VOLTAGE_THRESHOLD,d);}
    else{
      n->spiketime = t+2*DT;}
    return n->spiketime;}
}

double spikeguess1(struct neuron *n, double t, double DT)
{
  double V=0,Vs=0,sA=0,sN=0,sG=0,g=0,d=0,dt2=0;
  double cs = n->type>0 ? maximum(CS_ESES_A,CS_ESIS_A) : maximum(CS_ESEC_A,CS_ESIC_A);
  double spiketime_guess=0;
  if (n->spikenext > t && *(n->V) != VOLTAGE_RESET){ printf("%% neuron (%d,%d) fired but has V %f > VR in spikeguess1\n",n->row,n->col,*(n->V));}
  if (n->spikenext > t+DT){
    n->sV1=VOLTAGE_RESET;
    n->spiketime = t+2*DT;
    spiketime_guess = t+2*DT;}
  else if (n->spikenext >= t && n->spikenext <= t+DT){
    V = VOLTAGE_RESET;
    d = n->spikenext - t;
    dt2 = t + DT - n->spikenext;
    sA = *(n->sA)*exp(-0.0*d/TAU_AMPA) + SPIKETOL/5*cs/P_AMPA*exp(-1.0/AMPA_DIA)/PI/AMPA_DIA;
    sN = *(n->sN)*exp(-0.5*d/TAU_NMDA);
    sG = *(n->sG)*exp(-1.0*d/TAU_GABA);
    g = CONDUCTANCE_LK + sA + sN + sG;
    Vs = (7*sA + 7*sN - sG)/g*2.0/3.0;
    n->sV1 = V*exp(-dt2*g) + dt2*exp(-dt2/2*g)*Vs*g;
    spiketime_guess = t + dt2*(VOLTAGE_THRESHOLD-V)/(n->sV1-V);
    if (spiketime_guess < n->spikenext){ spiketime_guess = t+2*DT;}
    sA = *(n->sA)*exp(-0.0*d/TAU_AMPA);
    sN = *(n->sN)*exp(-0.5*d/TAU_NMDA);
    sG = *(n->sG)*exp(-1.0*d/TAU_GABA);
    g = CONDUCTANCE_LK + sA + sN + sG;
    Vs = (7*sA + 7*sN - sG)/g*2.0/3.0;
    n->sV1 = V*exp(-dt2*g) + dt2*exp(-dt2/2*g)*Vs*g;
    n->spiketime = t + dt2*(VOLTAGE_THRESHOLD-V)/(n->sV1-V);
    if (n->spiketime < n->spikenext){ n->spiketime = t+2*DT;}}
  else if (n->spikenext < t){
    V = *(n->V);
    sA = *(n->sA)*exp(-0.0*DT/TAU_AMPA) + SPIKETOL/5*cs/P_AMPA*exp(-1.0/AMPA_DIA)/PI/AMPA_DIA;
    sN = *(n->sN)*exp(-0.5*DT/TAU_NMDA);
    sG = *(n->sG)*exp(-1.0*DT/TAU_GABA);
    g = CONDUCTANCE_LK + sA + sN + sG;
    Vs = (7*sA + 7*sN - sG)/g*2.0/3.0;
    /* n->sV1 = V*exp(-DT*g) + DT*exp(-DT/2*g)*Vs*g; */ n->sV1 = V + DT*(g*(Vs-V));
    spiketime_guess = t + DT*(VOLTAGE_THRESHOLD-V)/(n->sV1-V);
    if (spiketime_guess < t){ spiketime_guess = t+2*DT;}
    sA = *(n->sA)*exp(-0.0*DT/TAU_AMPA);
    sN = *(n->sN)*exp(-0.5*DT/TAU_NMDA);
    sG = *(n->sG)*exp(-1.0*DT/TAU_GABA);
    g = CONDUCTANCE_LK + sA + sN + sG;
    Vs = (7*sA + 7*sN - sG)/g*2.0/3.0;
    n->sV1 = V*exp(-DT*g) + DT*exp(-DT/2*g)*Vs*g;
    n->spiketime = t + DT*(VOLTAGE_THRESHOLD-V)/(n->sV1-V);
    if (n->spiketime < t){ n->spiketime = t+2*DT;}}
  return spiketime_guess;
}

void spikesort_old(struct llist **Sra,struct llist *LL)
{
  int ar=0,ac=0,ar2=0,ac2=0,du=0,dl=0;
  struct litem *l=NULL;
  struct llist **vra=NULL,*vc=NULL,*vu=NULL,*vl=NULL;
  struct llist *Lc=NULL,*Lu=NULL,*Ll=NULL;
  int gnt=GLOBAL_NBLOCKS_TALL,gnw=GLOBAL_NBLOCKS_WIDE;
  vra = (struct llist **) tmalloc(sizeof(struct llist *)*gnt*gnw);
  for (ar=0;ar<gnt;ar++){ for (ac=0;ac<gnw;ac++){ vra[ar+ac*gnt]=NULL;}}
  for (ar=0;ar<gnt;ar++){ for (ac=0;ac<gnw;ac++){
    Lc = Sra[ar + ac*gnt]; 
    if (Lc->length > 0){
      vc = llistcopy(Lc); vra[ar + ac*gnt] = vc;  litemadd(LL,vc);
      Lu = NULL; vu = NULL; if (ar > 0){ Lu = Sra[(ar-1) + ac*gnt]; vu = vra[(ar-1) + ac*gnt];}
      Ll = NULL; vl = NULL; if (ac > 0){ Ll = Sra[ar + (ac-1)*gnt]; vl = vra[ar + (ac-1)*gnt];}
      du = distance(Lc,Lu);
      dl = distance(Lc,Ll);
      /* printf(" arbor(%d,%d), du=%d,dl=%d\n",ar,ac,du,dl); */
      /* llraprintf(vra,gnt,gnw);llprintf(LL); */
      if (dl <= GLOBAL_BLOCK_COL_DIA){
	/* printf("reconnecting left\n"); */
	if (vl!=vc){ lreconnect(vl,vc); litemminus(LL,vc);}
	vra[ar + ac*gnt] = vl;}
      vc = (struct llist *) vra[ar + ac*gnt];
      if (vc != vu && du <= GLOBAL_BLOCK_ROW_DIA){
	/* printf("reconnecting up\n"); */
	if (vu!=vc){ 
	  lreconnect(vu,vc); litemminus(LL,vc);
	  for (ar2=0;ar2<=ar;ar2++){ for (ac2=0;ac2<gnw;ac2++){ 
	    if (vra[ar2+ac2*gnt]==vc){ vra[ar2+ac2*gnt]=vu;}}}}}}
    else if (Lc->length == 0){
      /* do nothing? */}}}
  /* now match topmost and bottommost blocks */
  if (gnt>1){
    for (ac=0;ac<gnw;ac++){
      Lc = Sra[0 + ac*gnt]; vc = vra[0 + ac*gnt];
      Ll = Sra[gnt-1 + ac*gnt]; vl = vra[gnt-1 + ac*gnt];
      if (vc!=NULL && vl!=NULL && Lc->length > 0 && Ll->length > 0 && vc!=vl){
	dl = distance(Lc,Ll);
	if (dl <= GLOBAL_BLOCK_ROW_DIA){
	  /* printf("reconnecting across ns boundary\n"); */
	  lreconnect(vc,vl); litemminus(LL,vl);
	  for (ar2=0;ar2<gnt;ar2++){ for (ac2=0;ac2<gnw;ac2++){
	    if (vra[ar2+ac2*gnt]==vl){ vra[ar2+ac2*gnt]=vc;}}}}}}}
  /* now match leftmost and rightmost blocks */
  if (gnw>1){
    for (ar=0;ar<gnt;ar++){
      Lc = Sra[ar + 0*gnt]; vc = vra[ar + 0*gnt];
      Lu = Sra[ar + (gnw-1)*gnt]; vu = vra[ar + (gnw-1)*gnt];
      if (vc!=NULL && vu!=NULL && Lc->length > 0 && Lu->length > 0 && vc!=vu){
	du = distance(Lc,Lu);
	if (du <= GLOBAL_BLOCK_ROW_DIA){
	  /* printf("reconnecting across ew boundary\n"); */
	  lreconnect(vc,vu); litemminus(LL,vu);
	  for (ar2=0;ar2<gnt;ar2++){ for (ac2=0;ac2<gnw;ac2++){
	    if (vra[ar2+ac2*gnt]==vu){ vra[ar2+ac2*gnt]=vc;}}}}}}}
  /* llraprintf(vra,gnt,gnw); */
  tfree(vra);
  l = LL->first;
  while (l!=NULL){ 
    Lc = (struct llist *) l->item; 
    llistsort(Lc->first,Lc->last,Lc->length,&spiketime_compare);
    l = l->child;}
}

void spikesort(struct llist **Sra,struct llist *LL)
{
  int verbose=0;
  int ar=0,ac=0,tab=0,dww=0,dnw=0,dnn=0,dne=0;
  struct litem *l=NULL;
  struct llist **vra=NULL,*vc=NULL,*vww=NULL,*vnw=NULL,*vnn=NULL,*vne=NULL;
  struct llist *Lc=NULL,*Lww=NULL,*Lnw=NULL,*Lnn=NULL,*Lne=NULL;
  int gnt=GLOBAL_NBLOCKS_TALL,gnw=GLOBAL_NBLOCKS_WIDE;
  int *ira=NULL,verbose_connectivity=0;
  if (verbose){ printf(" %% [entering spikesort]\n");}
  vra = (struct llist **) tmalloc(sizeof(struct llist *)*gnt*gnw);
  if (verbose_connectivity){ 
    ira = (int *) tcalloc(gnw*gnt*9,sizeof(int)); 
    for (ar=0;ar<gnt;ar++){ for (ac=0;ac<gnw;ac++){ ira[(3*ar+1)+(3*ac+1)*3*gnt]=1;}}}
  for (ar=0;ar<gnt;ar++){ for (ac=0;ac<gnw;ac++){ 
    tab=ar+ac*gnt; vra[tab] = (Sra[tab]!=NULL && Sra[tab]->length>0) ?llistcopy(Sra[tab]) : NULL;}}
  for (ar=0;ar<gnt;ar++){ for (ac=0;ac<gnw;ac++){
    if (verbose>3){ printf(" %% currently at ar %d ac %d\n",ar,ac);}
    Lc = Sra[ar + ac*gnt]; 
    if (Lc->length > 0){
      vc = (struct llist *) vra[ar+ac*gnt];
      if (verbose_connectivity){ ira[(3*ar+1)+(3*ac+1)*3*gnt]=2;}
      tab=periodize(ar+0,0,gnt)+periodize(ac-1,0,gnw)*gnt;assert(tab>=0&&tab<gnw*gnt); Lww=Sra[tab]; vww=vra[tab]; dww = distance(Lc,Lww); 
      if (Lww!=NULL && dww <= GLOBAL_BLOCK_COL_DIA){
	if (verbose_connectivity){ ira[(3*ar+1)+(3*ac+0)*3*gnt]=2; ira[(3*periodize(ar+0,0,gnt)+1)+(3*periodize(ac-1,0,gnw)+2)*3*gnt]=2;}
	if (vww!=vc){ if (verbose>2){ printf(" %% ww\n");} lreconnect(vww,vc); spikesort_tagchange(vra,gnt,gnw,ar,ac,vc,vww); vc=vww;}}
      tab=periodize(ar-1,0,gnt)+periodize(ac-1,0,gnw)*gnt;assert(tab>=0&&tab<gnw*gnt); Lnw=Sra[tab]; vnw=vra[tab]; dnw = distance(Lc,Lnw);
      if (Lnw!=NULL && dnw <= maximum(GLOBAL_BLOCK_COL_DIA,GLOBAL_BLOCK_ROW_DIA)){
	if (verbose_connectivity){ ira[(3*ar+0)+(3*ac+0)*3*gnt]=2; ira[(3*periodize(ar-1,0,gnt)+2)+(3*periodize(ac-1,0,gnw)+2)*3*gnt]=2;}
	if (vnw!=vc){ if (verbose>2){ printf(" %% nw\n");} lreconnect(vnw,vc); spikesort_tagchange(vra,gnt,gnw,ar,ac,vc,vnw); vc=vnw;}}
      tab=periodize(ar-1,0,gnt)+periodize(ac+0,0,gnw)*gnt;assert(tab>=0&&tab<gnw*gnt); Lnn=Sra[tab]; vnn=vra[tab]; dnn = distance(Lc,Lnn); 
      if (Lnn!=NULL && dnn <= GLOBAL_BLOCK_ROW_DIA){
	if (verbose_connectivity){ ira[(3*ar+0)+(3*ac+1)*3*gnt]=2; ira[(3*periodize(ar-1,0,gnt)+2)+(3*periodize(ac+0,0,gnw)+1)*3*gnt]=2;}
	if (vnn!=vc){ if (verbose>2){ printf(" %% nn\n");} lreconnect(vnn,vc); spikesort_tagchange(vra,gnt,gnw,ar,ac,vc,vnn); vc=vnn;}}
      tab=periodize(ar-1,0,gnt)+periodize(ac+1,0,gnw)*gnt;assert(tab>=0&&tab<gnw*gnt); Lne=Sra[tab]; vne=vra[tab]; dne = distance(Lc,Lne);
      if (Lne!=NULL && dne <= maximum(GLOBAL_BLOCK_COL_DIA,GLOBAL_BLOCK_ROW_DIA)){
	if (verbose_connectivity){ ira[(3*ar+0)+(3*ac+2)*3*gnt]=2; ira[(3*periodize(ar-1,0,gnt)+2)+(3*periodize(ac+1,0,gnw)+0)*3*gnt]=2;}
	if (vne!=vc){ if (verbose>2){ printf(" %% ne\n");} lreconnect(vne,vc); spikesort_tagchange(vra,gnt,gnw,ar,ac,vc,vne); vc=vne;}}
      if (verbose>2){ printf(" %% finished with arbor(%d,%d), vc->length %d.\n",ar,ac,vra[ar+ac*gnt]->length);}}
    else /* if (Lc->length == 0) */{ /* do nothing */}}}
  if (verbose_connectivity){ printf(" %% cluster connectivity: \n"); for (ar=0;ar<gnt*3;ar++){ for (ac=0;ac<gnw*3;ac++){ printf("%s",ira[ar+ac*3*gnt]==0?" ":ira[ar+ac*3*gnt]==1?".":"X");} printf("\n");}}
  if (verbose){ llraprintf(vra,gnt,gnw);}
  for (ar=0;ar<gnt;ar++){ for (ac=0;ac<gnw;ac++){ 
    vc = vra[ar+ac*gnt];
    if (vc!=NULL){
      if (verbose>2){ printf(" %% adding llist %.3d->%d to LL\n",1+(int)vc%999,vc->length);}
      litemadd(LL,vc); spikesort_tagchange(vra,gnt,gnw,ar,ac,vc,NULL);}}}
  tfree(vra); vra=NULL;
  if (verbose){ llprintf(LL);}
  l = LL->first;
  while (l!=NULL){ 
    Lc = (struct llist *) l->item; 
    if (verbose>2){ printf(" %% sorting llist %.3d->%d\n",1+(int)Lc%999,Lc->length);}
    llistsort(Lc->first,Lc->last,Lc->length,&spiketime_compare);
    l = l->child;}
  if (verbose){ printf(" %% [finished with spikesort]\n");}
  if (ira!=NULL){ tfree(ira);ira=NULL;}
}

void spikesort_tagchange(struct llist **vra,int rows,int cols,int nr,int nc,struct llist *v_old,struct llist *v_new)
{
  if (v_old==v_new){ /* do nothing */}
  else /* if (v_old!=v_new) */{ 
    if (vra[nr+nc*rows]==v_old){
      vra[nr+nc*rows]=v_new;
      spikesort_tagchange(vra,rows,cols,periodize(nr-1,0,rows),periodize(nc-1,0,cols),v_old,v_new);
      spikesort_tagchange(vra,rows,cols,periodize(nr-1,0,rows),periodize(nc+0,0,cols),v_old,v_new);
      spikesort_tagchange(vra,rows,cols,periodize(nr-1,0,rows),periodize(nc+1,0,cols),v_old,v_new);
      spikesort_tagchange(vra,rows,cols,periodize(nr+0,0,rows),periodize(nc-1,0,cols),v_old,v_new);
      spikesort_tagchange(vra,rows,cols,periodize(nr+0,0,rows),periodize(nc+1,0,cols),v_old,v_new);
      spikesort_tagchange(vra,rows,cols,periodize(nr+1,0,rows),periodize(nc-1,0,cols),v_old,v_new);
      spikesort_tagchange(vra,rows,cols,periodize(nr+1,0,rows),periodize(nc+0,0,cols),v_old,v_new);
      spikesort_tagchange(vra,rows,cols,periodize(nr+1,0,rows),periodize(nc+1,0,cols),v_old,v_new);}}
}

void llprintf(struct llist *LL)
{
  struct litem *l=NULL,*l2=NULL;
  struct llist *L=NULL;
  struct neuron *n=NULL;
  if (LL==NULL){ printf("NULL list\n");}
  else if (LL->length == 0){ printf("LL=%d and length=%d\n",(int)LL,LL->length);}
  else if (LL->length >= 1){
    printf("LL=%d of length %d\n",(int)LL,LL->length);
    l = LL->first;
    while (l!=NULL){
      L = (struct llist *) l->item;
      printf("%.3d->%d:",1+(int)L%999,L->length);
      l2 = L->first; while (l2!=NULL){ n=(struct neuron *)l2->item; printf(" n(%d,%d)",n->row,n->col); l2=l2->child;} printf("\n");
      l = l->child;}}
}

void llraprintf(struct llist **vra, int rows, int cols)
{
  int nr=0,nc=0,d=0;
  printf("vra %d\n",(int)vra);
  for (nr=0;nr<rows;nr++){
    for (nc=0;nc<cols;nc++){
      d = (int)vra[nr + nc*rows];
      printf("%.3d->%.3d,",(d==0 ? 0 : 1+d%999),(vra[nr+nc*rows]==NULL ? 0 : vra[nr+nc*rows]->length));}
    printf("\n");}
}

void lreconnect(struct llist *L1, struct llist *L2)
{
  /*   printf("lreconnect L1=%d->%d,L2=%d->%d\n",1+(int)L1%999,L1->length,1+(int)L2%999,L2->length); */
  if (L1->length > 0 && L2->length > 0){
    L1->last->child = L2->first;
    L2->first->parent = L1->last;
    L1->last = L2->last;
    L1->length = L1->length + L2->length;}
  else if (L1->length > 0 && L2->length == 0){
    /* do nothing */;}
  else if (L1->length == 0 && L2->length > 0){
    L1->first = L2->first;
    L1->last = L2->last;
    L1->length = L2->length;}
  else if (L1->length == 0 && L2->length == 0){
    /* do nothing */;}
  tfree(L2);
}

int distance(struct llist *Lc, struct llist *Lp)
{
  double D=0,dr=0,dc=0,d=0;
  struct litem *lc=NULL,*lp=NULL;
  struct neuron *nc=NULL,*np=NULL;
  D = (double) 2*maximum(GLOBAL_BLOCK_ROW_DIA,GLOBAL_BLOCK_COL_DIA)+1;
  if (Lc==NULL || Lp==NULL){ return (int) D;}
  else if (Lc!=NULL && Lp!=NULL){
    lc = Lc->first;
    while (lc!=NULL){
      lp = Lp->first;
      while (lp!=NULL){
	nc = (struct neuron *) lc->item;
	np = (struct neuron *) lp->item;
	dr = fabs(nc->row-np->row); if (dr > GLOBAL_Nra->rows/2.0){ dr -= GLOBAL_Nra->rows;};
	dc = fabs(nc->col-np->col); if (dc > GLOBAL_Nra->cols/2.0){ dc -= GLOBAL_Nra->cols;}
	d = sqrt(dr*dr+dc*dc);
	if (d<D){ D=d;}
	lp = lp->child;}
      lc = lc->child;}
    return (int) ceil(D);}
  return (int) D;
}

int ilink(struct neuron *n,double *A,double *N,double *G)
{
  /* Here we assume that simple-complex differentiation happens at setinputrate */
  if (n->sox==-1){ *A=n->spikeinput_multiplicity*OTHERLAYER_STRENGTH;*N=0;*G=0; return 1;}
  else /* if (n->sox==+1) */{ *A=n->spikeinput_multiplicity*LGN_STRENGTH;*N=0;*G=0; return 1;}
}

int slink(struct neuron *s,struct neuron *n,double *A,double *N,double *G,double *Amax,double *Nmax,double *Gmax)
{
  /* returns synaptic links for s and n */
  int connected=0;
  double distr=0,distc=0,dist=0,gf=0,af=0,nf=0,newnum=0;
  *A=0;*N=0;*G=0;connected=0;
  if (AUTAPSES_OFF && s==n){ /* do nothing */ }
  else{
    if (s->homogenizing_knob==1 && s->sparse_link==NULL){
      switch ((int)*(s->t2s) * 4 + (int)*(n->t2s) * 1){
      case 0: /* inhibitory complex -> inhibitory complex */
	if (periodize((s->sparse_out-n->sparse_in)%SPARSEHIT_ICIC,0,SPARSEHIT_ICIC)<SPARSE_ICIC){
	  distr = fabs(s->row - n->row); if (distr > GLOBAL_Nra->rows/2){ distr -= GLOBAL_Nra->rows;}
	  distc = fabs(s->col - n->col); if (distc > GLOBAL_Nra->cols/2){ distc -= GLOBAL_Nra->cols;}
	  dist = distr*distr + distc*distc;
	  gf = (rand01 < P_GABA)*exp(-dist/(double)GABA_DIA)/P_GABA/PI/(double)GABA_DIA;
	  *A = 0; *N = 0; *G = CS_ICIC_G*gf; connected=1;}
	if (Amax!=NULL){ *Amax=0;} if (Nmax!=NULL){ *Nmax=0;} if (Gmax!=NULL){ *Gmax=CS_ICIC_G/P_GABA/PI/(double)GABA_DIA;} break;
      case 1: /* inhibitory complex -> inhibitory simple */
	if (periodize((s->sparse_out-n->sparse_in)%SPARSEHIT_ICIS,0,SPARSEHIT_ICIS)<SPARSE_ICIS){
	  distr = fabs(s->row - n->row); if (distr > GLOBAL_Nra->rows/2){ distr -= GLOBAL_Nra->rows;}
	  distc = fabs(s->col - n->col); if (distc > GLOBAL_Nra->cols/2){ distc -= GLOBAL_Nra->cols;}
	  dist = distr*distr + distc*distc;
	  gf = (rand01 < P_GABA)*exp(-dist/(double)GABA_DIA)/P_GABA/PI/(double)GABA_DIA;
	  *A = 0; *N = 0; *G = CS_ICIS_G*gf; connected=1;}
	if (Amax!=NULL){ *Amax=0;} if (Nmax!=NULL){ *Nmax=0;} if (Gmax!=NULL){ *Gmax=CS_ICIS_G/P_GABA/PI/(double)GABA_DIA;} break;
      case 2: /* inhibitory complex -> excitatory complex */
	if (periodize((s->sparse_out-n->sparse_in)%SPARSEHIT_ICEC,0,SPARSEHIT_ICEC)<SPARSE_ICEC){
	  distr = fabs(s->row - n->row); if (distr > GLOBAL_Nra->rows/2){ distr -= GLOBAL_Nra->rows;}
	  distc = fabs(s->col - n->col); if (distc > GLOBAL_Nra->cols/2){ distc -= GLOBAL_Nra->cols;}
	  dist = distr*distr + distc*distc;
	  gf = (rand01 < P_GABA)*exp(-dist/(double)GABA_DIA)/P_GABA/PI/(double)GABA_DIA;
	  *A = 0; *N = 0; *G = CS_ICEC_G*gf; connected=1;}
	if (Amax!=NULL){ *Amax=0;} if (Nmax!=NULL){ *Nmax=0;} if (Gmax!=NULL){ *Gmax=CS_ICEC_G/P_GABA/PI/(double)GABA_DIA;} break;
      case 3: /* inhibitory complex -> excitatory simple */
	if (periodize((s->sparse_out-n->sparse_in)%SPARSEHIT_ICES,0,SPARSEHIT_ICES)<SPARSE_ICES){
	  distr = fabs(s->row - n->row); if (distr > GLOBAL_Nra->rows/2){ distr -= GLOBAL_Nra->rows;}
	  distc = fabs(s->col - n->col); if (distc > GLOBAL_Nra->cols/2){ distc -= GLOBAL_Nra->cols;}
	  dist = distr*distr + distc*distc;
	  gf = (rand01 < P_GABA)*exp(-dist/(double)GABA_DIA)/P_GABA/PI/(double)GABA_DIA;
	  *A = 0; *N = 0; *G = CS_ICES_G*gf; connected=1;}
	if (Amax!=NULL){ *Amax=0;} if (Nmax!=NULL){ *Nmax=0;} if (Gmax!=NULL){ *Gmax=CS_ICES_G/P_GABA/PI/(double)GABA_DIA;} break;
      case 4: /* inhibitory simple -> inhibitory complex */
	if (periodize((s->sparse_out-n->sparse_in)%SPARSEHIT_ISIC,0,SPARSEHIT_ISIC)<SPARSE_ISIC){
	  distr = fabs(s->row - n->row); if (distr > GLOBAL_Nra->rows/2){ distr -= GLOBAL_Nra->rows;}
	  distc = fabs(s->col - n->col); if (distc > GLOBAL_Nra->cols/2){ distc -= GLOBAL_Nra->cols;}
	  dist = distr*distr + distc*distc;
	  gf = (rand01 < P_GABA)*exp(-dist/(double)GABA_DIA)/P_GABA/PI/(double)GABA_DIA;
	  *A = 0; *N = 0; *G = CS_ISIC_G*gf; connected=1;}
	if (Amax!=NULL){ *Amax=0;} if (Nmax!=NULL){ *Nmax=0;} if (Gmax!=NULL){ *Gmax=CS_ISIC_G/P_GABA/PI/(double)GABA_DIA;} break;
      case 5: /* inhibitory simple -> inhibitory simple */
	if (periodize((s->sparse_out-n->sparse_in)%SPARSEHIT_ISIS,0,SPARSEHIT_ISIS)<SPARSE_ISIS){
	  distr = fabs(s->row - n->row); if (distr > GLOBAL_Nra->rows/2){ distr -= GLOBAL_Nra->rows;}
	  distc = fabs(s->col - n->col); if (distc > GLOBAL_Nra->cols/2){ distc -= GLOBAL_Nra->cols;}
	  dist = distr*distr + distc*distc;
	  gf = (rand01 < P_GABA)*exp(-dist/(double)GABA_DIA)/P_GABA/PI/(double)GABA_DIA;
	  *A = 0; *N = 0; *G = CS_ISIS_G*gf; connected=1;}
	if (Amax!=NULL){ *Amax=0;} if (Nmax!=NULL){ *Nmax=0;} if (Gmax!=NULL){ *Gmax=CS_ISIS_G/P_GABA/PI/(double)GABA_DIA;} break;
      case 6: /* inhibitory simple -> excitatory complex */
	if (periodize((s->sparse_out-n->sparse_in)%SPARSEHIT_ISEC,0,SPARSEHIT_ISEC)<SPARSE_ISEC){
	  distr = fabs(s->row - n->row); if (distr > GLOBAL_Nra->rows/2){ distr -= GLOBAL_Nra->rows;}
	  distc = fabs(s->col - n->col); if (distc > GLOBAL_Nra->cols/2){ distc -= GLOBAL_Nra->cols;}
	  dist = distr*distr + distc*distc;
	  gf = (rand01 < P_GABA)*exp(-dist/(double)GABA_DIA)/P_GABA/PI/(double)GABA_DIA;
	  *A = 0; *N = 0; *G = CS_ISEC_G*gf; connected=1;}
	if (Amax!=NULL){ *Amax=0;} if (Nmax!=NULL){ *Nmax=0;} if (Gmax!=NULL){ *Gmax=CS_ISEC_G/P_GABA/PI/(double)GABA_DIA;} break;
      case 7: /* inhibitory simple -> excitatory simple */
	if (periodize((s->sparse_out-n->sparse_in)%SPARSEHIT_ISES,0,SPARSEHIT_ISES)<SPARSE_ISES){
	  distr = fabs(s->row - n->row); if (distr > GLOBAL_Nra->rows/2){ distr -= GLOBAL_Nra->rows;}
	  distc = fabs(s->col - n->col); if (distc > GLOBAL_Nra->cols/2){ distc -= GLOBAL_Nra->cols;}
	  dist = distr*distr + distc*distc;
	  gf = (rand01 < P_GABA)*exp(-dist/(double)GABA_DIA)/P_GABA/PI/(double)GABA_DIA;
	  *A = 0; *N = 0; *G = CS_ISES_G*gf; connected=1;}
	if (Amax!=NULL){ *Amax=0;} if (Nmax!=NULL){ *Nmax=0;} if (Gmax!=NULL){ *Gmax=CS_ISES_G/P_GABA/PI/(double)GABA_DIA;} break;
      case 8: /* excitatory complex -> inhibitory complex */
	if (periodize((s->sparse_out-n->sparse_in)%SPARSEHIT_ECIC,0,SPARSEHIT_ECIC)<SPARSE_ECIC){
	  distr = fabs(s->row - n->row); if (distr > GLOBAL_Nra->rows/2){ distr -= GLOBAL_Nra->rows;}
	  distc = fabs(s->col - n->col); if (distc > GLOBAL_Nra->cols/2){ distc -= GLOBAL_Nra->cols;}
	  dist = distr*distr + distc*distc;
	  af = (rand01 < P_AMPA)*exp(-dist/(double)AMPA_DIA)/P_AMPA/PI/(double)AMPA_DIA;
	  nf = (rand01 < P_NMDA)*exp(-dist/(double)NMDA_DIA)/P_NMDA/PI/(double)NMDA_DIA;
	  *A = CS_ECIC_A*af; *N = CS_ECIC_N*nf; *G = 0; connected=1;}
	if (Amax!=NULL){ *Amax=CS_ECIC_A/P_AMPA/PI/(double)AMPA_DIA;} if (Nmax!=NULL){ *Nmax=CS_ECIC_N/P_NMDA/PI/(double)NMDA_DIA;} if (Gmax!=NULL){ *Gmax=0;} break;
      case 9: /* excitatory complex -> inhibitory simple */
	if (periodize((s->sparse_out-n->sparse_in)%SPARSEHIT_ECIS,0,SPARSEHIT_ECIS)<SPARSE_ECIS){
	  distr = fabs(s->row - n->row); if (distr > GLOBAL_Nra->rows/2){ distr -= GLOBAL_Nra->rows;}
	  distc = fabs(s->col - n->col); if (distc > GLOBAL_Nra->cols/2){ distc -= GLOBAL_Nra->cols;}
	  dist = distr*distr + distc*distc;
	  af = (rand01 < P_AMPA)*exp(-dist/(double)AMPA_DIA)/P_AMPA/PI/(double)AMPA_DIA;
	  nf = (rand01 < P_NMDA)*exp(-dist/(double)NMDA_DIA)/P_NMDA/PI/(double)NMDA_DIA;
	  *A = CS_ECIS_A*af; *N = CS_ECIS_N*nf; *G = 0; connected=1;}
	if (Amax!=NULL){ *Amax=CS_ECIS_A/P_AMPA/PI/(double)AMPA_DIA;} if (Nmax!=NULL){ *Nmax=CS_ECIS_N/P_NMDA/PI/(double)NMDA_DIA;} if (Gmax!=NULL){ *Gmax=0;} break;
      case 10: /* excitatory complex -> excitatory complex */
	if (periodize((s->sparse_out-n->sparse_in)%SPARSEHIT_ECEC,0,SPARSEHIT_ECEC)<SPARSE_ECEC){
	  distr = fabs(s->row - n->row); if (distr > GLOBAL_Nra->rows/2){ distr -= GLOBAL_Nra->rows;}
	  distc = fabs(s->col - n->col); if (distc > GLOBAL_Nra->cols/2){ distc -= GLOBAL_Nra->cols;}
	  dist = distr*distr + distc*distc;
	  af = (rand01 < P_AMPA)*exp(-dist/(double)AMPA_DIA)/P_AMPA/PI/(double)AMPA_DIA;
	  nf = (rand01 < P_NMDA)*exp(-dist/(double)NMDA_DIA)/P_NMDA/PI/(double)NMDA_DIA;
	  *A = CS_ECEC_A*af; *N = CS_ECEC_N*nf; *G = 0; connected=1;}
	if (Amax!=NULL){ *Amax=CS_ECEC_A/P_AMPA/PI/(double)AMPA_DIA;} if (Nmax!=NULL){ *Nmax=CS_ECEC_N/P_NMDA/PI/(double)NMDA_DIA;} if (Gmax!=NULL){ *Gmax=0;} break;
      case 11: /* excitatory complex -> excitatory simple */
	if (periodize((s->sparse_out-n->sparse_in)%SPARSEHIT_ECES,0,SPARSEHIT_ECES)<SPARSE_ECES){
	  distr = fabs(s->row - n->row); if (distr > GLOBAL_Nra->rows/2){ distr -= GLOBAL_Nra->rows;}
	  distc = fabs(s->col - n->col); if (distc > GLOBAL_Nra->cols/2){ distc -= GLOBAL_Nra->cols;}
	  dist = distr*distr + distc*distc;
	  af = (rand01 < P_AMPA)*exp(-dist/(double)AMPA_DIA)/P_AMPA/PI/(double)AMPA_DIA;
	  nf = (rand01 < P_NMDA)*exp(-dist/(double)NMDA_DIA)/P_NMDA/PI/(double)NMDA_DIA;
	  *A = CS_ECES_A*af; *N = CS_ECES_N*nf; *G = 0; connected=1;}
	if (Amax!=NULL){ *Amax=CS_ECES_A/P_AMPA/PI/(double)AMPA_DIA;} if (Nmax!=NULL){ *Nmax=CS_ECES_N/P_NMDA/PI/(double)NMDA_DIA;} if (Gmax!=NULL){ *Gmax=0;} break;
      case 12: /* excitatory simple -> inhibitory complex */
	if (periodize((s->sparse_out-n->sparse_in)%SPARSEHIT_ESIC,0,SPARSEHIT_ESIC)<SPARSE_ESIC){
	  distr = fabs(s->row - n->row); if (distr > GLOBAL_Nra->rows/2){ distr -= GLOBAL_Nra->rows;}
	  distc = fabs(s->col - n->col); if (distc > GLOBAL_Nra->cols/2){ distc -= GLOBAL_Nra->cols;}
	  dist = distr*distr + distc*distc;
	  af = (rand01 < P_AMPA)*exp(-dist/(double)AMPA_DIA)/P_AMPA/PI/(double)AMPA_DIA;
	  nf = (rand01 < P_NMDA)*exp(-dist/(double)NMDA_DIA)/P_NMDA/PI/(double)NMDA_DIA;
	  *A = CS_ESIC_A*af; *N = CS_ESIC_N*nf; *G = 0; connected=1;}
	if (Amax!=NULL){ *Amax=CS_ESIC_A/P_AMPA/PI/(double)AMPA_DIA;} if (Nmax!=NULL){ *Nmax=CS_ESIC_N/P_NMDA/PI/(double)NMDA_DIA;} if (Gmax!=NULL){ *Gmax=0;} break;
      case 13: /* excitatory simple -> inhibitory simple */
	if (periodize((s->sparse_out-n->sparse_in)%SPARSEHIT_ESIS,0,SPARSEHIT_ESIS)<SPARSE_ESIS){
	  distr = fabs(s->row - n->row); if (distr > GLOBAL_Nra->rows/2){ distr -= GLOBAL_Nra->rows;}
	  distc = fabs(s->col - n->col); if (distc > GLOBAL_Nra->cols/2){ distc -= GLOBAL_Nra->cols;}
	  dist = distr*distr + distc*distc;
	  af = (rand01 < P_AMPA)*exp(-dist/(double)AMPA_DIA)/P_AMPA/PI/(double)AMPA_DIA;
	  nf = (rand01 < P_NMDA)*exp(-dist/(double)NMDA_DIA)/P_NMDA/PI/(double)NMDA_DIA;
	  *A = CS_ESIS_A*af; *N = CS_ESIS_N*nf; *G = 0; connected=1;}
	if (Amax!=NULL){ *Amax=CS_ESIS_A/P_AMPA/PI/(double)AMPA_DIA;} if (Nmax!=NULL){ *Nmax=CS_ESIS_N/P_NMDA/PI/(double)NMDA_DIA;} if (Gmax!=NULL){ *Gmax=0;} break;
      case 14: /* excitatory simple -> excitatory complex */
	if (periodize((s->sparse_out-n->sparse_in)%SPARSEHIT_ESEC,0,SPARSEHIT_ESEC)<SPARSE_ESEC){
	  distr = fabs(s->row - n->row); if (distr > GLOBAL_Nra->rows/2){ distr -= GLOBAL_Nra->rows;}
	  distc = fabs(s->col - n->col); if (distc > GLOBAL_Nra->cols/2){ distc -= GLOBAL_Nra->cols;}
	  dist = distr*distr + distc*distc;
	  af = (rand01 < P_AMPA)*exp(-dist/(double)AMPA_DIA)/P_AMPA/PI/(double)AMPA_DIA;
	  nf = (rand01 < P_NMDA)*exp(-dist/(double)NMDA_DIA)/P_NMDA/PI/(double)NMDA_DIA;
	  *A = CS_ESEC_A*af; *N = CS_ESEC_N*nf; *G = 0; connected=1;}
	if (Amax!=NULL){ *Amax=CS_ESEC_A/P_AMPA/PI/(double)AMPA_DIA;} if (Nmax!=NULL){ *Nmax=CS_ESEC_N/P_NMDA/PI/(double)NMDA_DIA;} if (Gmax!=NULL){ *Gmax=0;} break;
      case 15: /* excitatory simple -> excitatory simple */
	if (periodize((s->sparse_out-n->sparse_in)%SPARSEHIT_ESES,0,SPARSEHIT_ESES)<SPARSE_ESES){
	  distr = fabs(s->row - n->row); if (distr > GLOBAL_Nra->rows/2){ distr -= GLOBAL_Nra->rows;}
	  distc = fabs(s->col - n->col); if (distc > GLOBAL_Nra->cols/2){ distc -= GLOBAL_Nra->cols;}
	  dist = distr*distr + distc*distc;
	  af = (rand01 < P_AMPA)*exp(-dist/(double)AMPA_DIA)/P_AMPA/PI/(double)AMPA_DIA;
	  nf = (rand01 < P_NMDA)*exp(-dist/(double)NMDA_DIA)/P_NMDA/PI/(double)NMDA_DIA;
	  *A = CS_ESES_A*af; *N = CS_ESES_N*nf; *G = 0; connected=1;}
	if (Amax!=NULL){ *Amax=CS_ESES_A/P_AMPA/PI/(double)AMPA_DIA;} if (Nmax!=NULL){ *Nmax=CS_ESES_N/P_NMDA/PI/(double)NMDA_DIA;} if (Gmax!=NULL){ *Gmax=0;} break;
      default: printf(" %% error, wrong types (%d,%d)->(%d,%d) in slink\n",s->type,s->sox,n->type,n->sox); break;}}
    else if (s->homogenizing_knob<1 && s->sparse_link==NULL){
      switch ((int)*(s->t2s) * 4 + (int)*(n->t2s) * 1){
      case 0: /* inhibitory complex -> inhibitory complex */
	newnum = minimum(SPARSEHIT_ICIC,(double)SPARSE_ICIC/s->homogenizing_knob);
	if (rand01 < (double)SPARSE_ICIC/newnum){
	  if (periodize((s->sparse_out-n->sparse_in)%SPARSEHIT_ICIC,0,SPARSEHIT_ICIC)<newnum){
	    distr = fabs(s->row - n->row); if (distr > GLOBAL_Nra->rows/2){ distr -= GLOBAL_Nra->rows;}
	    distc = fabs(s->col - n->col); if (distc > GLOBAL_Nra->cols/2){ distc -= GLOBAL_Nra->cols;}
	    dist = distr*distr + distc*distc;
	    gf = (rand01 < P_GABA)*exp(-dist/(double)GABA_DIA)/P_GABA/PI/(double)GABA_DIA;
	    *A = 0; *N = 0; *G = CS_ICIC_G*gf; connected=1;}}
	if (Amax!=NULL){ *Amax=0;} if (Nmax!=NULL){ *Nmax=0;} if (Gmax!=NULL){ *Gmax=CS_ICIC_G/P_GABA/PI/(double)GABA_DIA;} break;
      case 1: /* inhibitory complex -> inhibitory simple */
	newnum = minimum(SPARSEHIT_ICIS,(double)SPARSE_ICIS/s->homogenizing_knob);
	if (rand01 < (double)SPARSE_ICIS/newnum){
	  if (periodize((s->sparse_out-n->sparse_in)%SPARSEHIT_ICIS,0,SPARSEHIT_ICIS)<newnum){
	    distr = fabs(s->row - n->row); if (distr > GLOBAL_Nra->rows/2){ distr -= GLOBAL_Nra->rows;}
	    distc = fabs(s->col - n->col); if (distc > GLOBAL_Nra->cols/2){ distc -= GLOBAL_Nra->cols;}
	    dist = distr*distr + distc*distc;
	    gf = (rand01 < P_GABA)*exp(-dist/(double)GABA_DIA)/P_GABA/PI/(double)GABA_DIA;
	    *A = 0; *N = 0; *G = CS_ICIS_G*gf; connected=1;}}
	if (Amax!=NULL){ *Amax=0;} if (Nmax!=NULL){ *Nmax=0;} if (Gmax!=NULL){ *Gmax=CS_ICIS_G/P_GABA/PI/(double)GABA_DIA;} break;
      case 2: /* inhibitory complex -> excitatory complex */
	newnum = minimum(SPARSEHIT_ICEC,(double)SPARSE_ICEC/s->homogenizing_knob);
	if (rand01 < (double)SPARSE_ICEC/newnum){
	  if (periodize((s->sparse_out-n->sparse_in)%SPARSEHIT_ICEC,0,SPARSEHIT_ICEC)<newnum){
	    distr = fabs(s->row - n->row); if (distr > GLOBAL_Nra->rows/2){ distr -= GLOBAL_Nra->rows;}
	    distc = fabs(s->col - n->col); if (distc > GLOBAL_Nra->cols/2){ distc -= GLOBAL_Nra->cols;}
	    dist = distr*distr + distc*distc;
	    gf = (rand01 < P_GABA)*exp(-dist/(double)GABA_DIA)/P_GABA/PI/(double)GABA_DIA;
	    *A = 0; *N = 0; *G = CS_ICEC_G*gf; connected=1;}}
	if (Amax!=NULL){ *Amax=0;} if (Nmax!=NULL){ *Nmax=0;} if (Gmax!=NULL){ *Gmax=CS_ICEC_G/P_GABA/PI/(double)GABA_DIA;} break;
      case 3: /* inhibitory complex -> excitatory simple */
	newnum = minimum(SPARSEHIT_ICES,(double)SPARSE_ICES/s->homogenizing_knob);
	if (rand01 < (double)SPARSE_ICES/newnum){
	  if (periodize((s->sparse_out-n->sparse_in)%SPARSEHIT_ICES,0,SPARSEHIT_ICES)<newnum){
	    distr = fabs(s->row - n->row); if (distr > GLOBAL_Nra->rows/2){ distr -= GLOBAL_Nra->rows;}
	    distc = fabs(s->col - n->col); if (distc > GLOBAL_Nra->cols/2){ distc -= GLOBAL_Nra->cols;}
	    dist = distr*distr + distc*distc;
	    gf = (rand01 < P_GABA)*exp(-dist/(double)GABA_DIA)/P_GABA/PI/(double)GABA_DIA;
	    *A = 0; *N = 0; *G = CS_ICES_G*gf; connected=1;}}
	if (Amax!=NULL){ *Amax=0;} if (Nmax!=NULL){ *Nmax=0;} if (Gmax!=NULL){ *Gmax=CS_ICES_G/P_GABA/PI/(double)GABA_DIA;} break;
      case 4: /* inhibitory simple -> inhibitory complex */
	newnum = minimum(SPARSEHIT_ISIC,(double)SPARSE_ISIC/s->homogenizing_knob);
	if (rand01 < (double)SPARSE_ISIC/newnum){
	  if (periodize((s->sparse_out-n->sparse_in)%SPARSEHIT_ISIC,0,SPARSEHIT_ISIC)<newnum){
	    distr = fabs(s->row - n->row); if (distr > GLOBAL_Nra->rows/2){ distr -= GLOBAL_Nra->rows;}
	    distc = fabs(s->col - n->col); if (distc > GLOBAL_Nra->cols/2){ distc -= GLOBAL_Nra->cols;}
	    dist = distr*distr + distc*distc;
	    gf = (rand01 < P_GABA)*exp(-dist/(double)GABA_DIA)/P_GABA/PI/(double)GABA_DIA;
	    *A = 0; *N = 0; *G = CS_ISIC_G*gf; connected=1;}}
	if (Amax!=NULL){ *Amax=0;} if (Nmax!=NULL){ *Nmax=0;} if (Gmax!=NULL){ *Gmax=CS_ISIC_G/P_GABA/PI/(double)GABA_DIA;} break;
      case 5: /* inhibitory simple -> inhibitory simple */
	newnum = minimum(SPARSEHIT_ISIS,(double)SPARSE_ISIS/s->homogenizing_knob);
	if (rand01 < (double)SPARSE_ISIS/newnum){
	  if (periodize((s->sparse_out-n->sparse_in)%SPARSEHIT_ISIS,0,SPARSEHIT_ISIS)<newnum){
	    distr = fabs(s->row - n->row); if (distr > GLOBAL_Nra->rows/2){ distr -= GLOBAL_Nra->rows;}
	    distc = fabs(s->col - n->col); if (distc > GLOBAL_Nra->cols/2){ distc -= GLOBAL_Nra->cols;}
	    dist = distr*distr + distc*distc;
	    gf = (rand01 < P_GABA)*exp(-dist/(double)GABA_DIA)/P_GABA/PI/(double)GABA_DIA;
	    *A = 0; *N = 0; *G = CS_ISIS_G*gf; connected=1;}}
	if (Amax!=NULL){ *Amax=0;} if (Nmax!=NULL){ *Nmax=0;} if (Gmax!=NULL){ *Gmax=CS_ISIS_G/P_GABA/PI/(double)GABA_DIA;} break;
      case 6: /* inhibitory simple -> excitatory complex */
	newnum = minimum(SPARSEHIT_ISEC,(double)SPARSE_ISEC/s->homogenizing_knob);
	if (rand01 < (double)SPARSE_ISEC/newnum){
	  if (periodize((s->sparse_out-n->sparse_in)%SPARSEHIT_ISEC,0,SPARSEHIT_ISEC)<newnum){
	    distr = fabs(s->row - n->row); if (distr > GLOBAL_Nra->rows/2){ distr -= GLOBAL_Nra->rows;}
	    distc = fabs(s->col - n->col); if (distc > GLOBAL_Nra->cols/2){ distc -= GLOBAL_Nra->cols;}
	    dist = distr*distr + distc*distc;
	    gf = (rand01 < P_GABA)*exp(-dist/(double)GABA_DIA)/P_GABA/PI/(double)GABA_DIA;
	    *A = 0; *N = 0; *G = CS_ISEC_G*gf; connected=1;}}
	if (Amax!=NULL){ *Amax=0;} if (Nmax!=NULL){ *Nmax=0;} if (Gmax!=NULL){ *Gmax=CS_ISEC_G/P_GABA/PI/(double)GABA_DIA;} break;
      case 7: /* inhibitory simple -> excitatory simple */
	newnum = minimum(SPARSEHIT_ISES,(double)SPARSE_ISES/s->homogenizing_knob);
	if (rand01 < (double)SPARSE_ISES/newnum){
	  if (periodize((s->sparse_out-n->sparse_in)%SPARSEHIT_ISES,0,SPARSEHIT_ISES)<newnum){
	    distr = fabs(s->row - n->row); if (distr > GLOBAL_Nra->rows/2){ distr -= GLOBAL_Nra->rows;}
	    distc = fabs(s->col - n->col); if (distc > GLOBAL_Nra->cols/2){ distc -= GLOBAL_Nra->cols;}
	    dist = distr*distr + distc*distc;
	    gf = (rand01 < P_GABA)*exp(-dist/(double)GABA_DIA)/P_GABA/PI/(double)GABA_DIA;
	    *A = 0; *N = 0; *G = CS_ISES_G*gf; connected=1;}}
	if (Amax!=NULL){ *Amax=0;} if (Nmax!=NULL){ *Nmax=0;} if (Gmax!=NULL){ *Gmax=CS_ISES_G/P_GABA/PI/(double)GABA_DIA;} break;
      case 8: /* excitatory complex -> inhibitory complex */
	newnum = minimum(SPARSEHIT_ECIC,(double)SPARSE_ECIC/s->homogenizing_knob);
	if (rand01 < (double)SPARSE_ECIC/newnum){
	  if (periodize((s->sparse_out-n->sparse_in)%SPARSEHIT_ECIC,0,SPARSEHIT_ECIC)<newnum){
	    distr = fabs(s->row - n->row); if (distr > GLOBAL_Nra->rows/2){ distr -= GLOBAL_Nra->rows;}
	    distc = fabs(s->col - n->col); if (distc > GLOBAL_Nra->cols/2){ distc -= GLOBAL_Nra->cols;}
	    dist = distr*distr + distc*distc;
	    af = (rand01 < P_AMPA)*exp(-dist/(double)AMPA_DIA)/P_AMPA/PI/(double)AMPA_DIA;
	    nf = (rand01 < P_NMDA)*exp(-dist/(double)NMDA_DIA)/P_NMDA/PI/(double)NMDA_DIA;
	    *A = CS_ECIC_A*af; *N = CS_ECIC_N*nf; *G = 0; connected=1;}}
	if (Amax!=NULL){ *Amax=CS_ECIC_A/P_AMPA/PI/(double)AMPA_DIA;} if (Nmax!=NULL){ *Nmax=CS_ECIC_N/P_NMDA/PI/(double)NMDA_DIA;} if (Gmax!=NULL){ *Gmax=0;} break;
      case 9: /* excitatory complex -> inhibitory simple */
	newnum = minimum(SPARSEHIT_ECIS,(double)SPARSE_ECIS/s->homogenizing_knob);
	if (rand01 < (double)SPARSE_ECIS/newnum){
	  if (periodize((s->sparse_out-n->sparse_in)%SPARSEHIT_ECIS,0,SPARSEHIT_ECIS)<newnum){
	    distr = fabs(s->row - n->row); if (distr > GLOBAL_Nra->rows/2){ distr -= GLOBAL_Nra->rows;}
	    distc = fabs(s->col - n->col); if (distc > GLOBAL_Nra->cols/2){ distc -= GLOBAL_Nra->cols;}
	    dist = distr*distr + distc*distc;
	    af = (rand01 < P_AMPA)*exp(-dist/(double)AMPA_DIA)/P_AMPA/PI/(double)AMPA_DIA;
	    nf = (rand01 < P_NMDA)*exp(-dist/(double)NMDA_DIA)/P_NMDA/PI/(double)NMDA_DIA;
	    *A = CS_ECIS_A*af; *N = CS_ECIS_N*nf; *G = 0; connected=1;}}
	if (Amax!=NULL){ *Amax=CS_ECIS_A/P_AMPA/PI/(double)AMPA_DIA;} if (Nmax!=NULL){ *Nmax=CS_ECIS_N/P_NMDA/PI/(double)NMDA_DIA;} if (Gmax!=NULL){ *Gmax=0;} break;
      case 10: /* excitatory complex -> excitatory complex */
	newnum = minimum(SPARSEHIT_ECEC,(double)SPARSE_ECEC/s->homogenizing_knob);
	if (rand01 < (double)SPARSE_ECEC/newnum){
	  if (periodize((s->sparse_out-n->sparse_in)%SPARSEHIT_ECEC,0,SPARSEHIT_ECEC)<newnum){
	    distr = fabs(s->row - n->row); if (distr > GLOBAL_Nra->rows/2){ distr -= GLOBAL_Nra->rows;}
	    distc = fabs(s->col - n->col); if (distc > GLOBAL_Nra->cols/2){ distc -= GLOBAL_Nra->cols;}
	    dist = distr*distr + distc*distc;
	    af = (rand01 < P_AMPA)*exp(-dist/(double)AMPA_DIA)/P_AMPA/PI/(double)AMPA_DIA;
	    nf = (rand01 < P_NMDA)*exp(-dist/(double)NMDA_DIA)/P_NMDA/PI/(double)NMDA_DIA;
	    *A = CS_ECEC_A*af; *N = CS_ECEC_N*nf; *G = 0; connected=1;}}
	if (Amax!=NULL){ *Amax=CS_ECEC_A/P_AMPA/PI/(double)AMPA_DIA;} if (Nmax!=NULL){ *Nmax=CS_ECEC_N/P_NMDA/PI/(double)NMDA_DIA;} if (Gmax!=NULL){ *Gmax=0;} break;
      case 11: /* excitatory complex -> excitatory simple */
	newnum = minimum(SPARSEHIT_ECES,(double)SPARSE_ECES/s->homogenizing_knob);
	if (rand01 < (double)SPARSE_ECES/newnum){
	  if (periodize((s->sparse_out-n->sparse_in)%SPARSEHIT_ECES,0,SPARSEHIT_ECES)<newnum){
	    distr = fabs(s->row - n->row); if (distr > GLOBAL_Nra->rows/2){ distr -= GLOBAL_Nra->rows;}
	    distc = fabs(s->col - n->col); if (distc > GLOBAL_Nra->cols/2){ distc -= GLOBAL_Nra->cols;}
	    dist = distr*distr + distc*distc;
	    af = (rand01 < P_AMPA)*exp(-dist/(double)AMPA_DIA)/P_AMPA/PI/(double)AMPA_DIA;
	    nf = (rand01 < P_NMDA)*exp(-dist/(double)NMDA_DIA)/P_NMDA/PI/(double)NMDA_DIA;
	    *A = CS_ECES_A*af; *N = CS_ECES_N*nf; *G = 0; connected=1;}}
	if (Amax!=NULL){ *Amax=CS_ECES_A/P_AMPA/PI/(double)AMPA_DIA;} if (Nmax!=NULL){ *Nmax=CS_ECES_N/P_NMDA/PI/(double)NMDA_DIA;} if (Gmax!=NULL){ *Gmax=0;} break;
      case 12: /* excitatory simple -> inhibitory complex */
	newnum = minimum(SPARSEHIT_ESIC,(double)SPARSE_ESIC/s->homogenizing_knob);
	if (rand01 < (double)SPARSE_ESIC/newnum){
	  if (periodize((s->sparse_out-n->sparse_in)%SPARSEHIT_ESIC,0,SPARSEHIT_ESIC)<newnum){
	    distr = fabs(s->row - n->row); if (distr > GLOBAL_Nra->rows/2){ distr -= GLOBAL_Nra->rows;}
	    distc = fabs(s->col - n->col); if (distc > GLOBAL_Nra->cols/2){ distc -= GLOBAL_Nra->cols;}
	    dist = distr*distr + distc*distc;
	    af = (rand01 < P_AMPA)*exp(-dist/(double)AMPA_DIA)/P_AMPA/PI/(double)AMPA_DIA;
	    nf = (rand01 < P_NMDA)*exp(-dist/(double)NMDA_DIA)/P_NMDA/PI/(double)NMDA_DIA;
	    *A = CS_ESIC_A*af; *N = CS_ESIC_N*nf; *G = 0; connected=1;}}
	if (Amax!=NULL){ *Amax=CS_ESIC_A/P_AMPA/PI/(double)AMPA_DIA;} if (Nmax!=NULL){ *Nmax=CS_ESIC_N/P_NMDA/PI/(double)NMDA_DIA;} if (Gmax!=NULL){ *Gmax=0;} break;
      case 13: /* excitatory simple -> inhibitory simple */
	newnum = minimum(SPARSEHIT_ESIS,(double)SPARSE_ESIS/s->homogenizing_knob);
	if (rand01 < (double)SPARSE_ESIS/newnum){
	  if (periodize((s->sparse_out-n->sparse_in)%SPARSEHIT_ESIS,0,SPARSEHIT_ESIS)<newnum){
	    distr = fabs(s->row - n->row); if (distr > GLOBAL_Nra->rows/2){ distr -= GLOBAL_Nra->rows;}
	    distc = fabs(s->col - n->col); if (distc > GLOBAL_Nra->cols/2){ distc -= GLOBAL_Nra->cols;}
	    dist = distr*distr + distc*distc;
	    af = (rand01 < P_AMPA)*exp(-dist/(double)AMPA_DIA)/P_AMPA/PI/(double)AMPA_DIA;
	    nf = (rand01 < P_NMDA)*exp(-dist/(double)NMDA_DIA)/P_NMDA/PI/(double)NMDA_DIA;
	    *A = CS_ESIS_A*af; *N = CS_ESIS_N*nf; *G = 0; connected=1;}}
	if (Amax!=NULL){ *Amax=CS_ESIS_A/P_AMPA/PI/(double)AMPA_DIA;} if (Nmax!=NULL){ *Nmax=CS_ESIS_N/P_NMDA/PI/(double)NMDA_DIA;} if (Gmax!=NULL){ *Gmax=0;} break;
      case 14: /* excitatory simple -> excitatory complex */
	newnum = minimum(SPARSEHIT_ESEC,(double)SPARSE_ESEC/s->homogenizing_knob);
	if (rand01 < (double)SPARSE_ESEC/newnum){
	  if (periodize((s->sparse_out-n->sparse_in)%SPARSEHIT_ESEC,0,SPARSEHIT_ESEC)<newnum){
	    distr = fabs(s->row - n->row); if (distr > GLOBAL_Nra->rows/2){ distr -= GLOBAL_Nra->rows;}
	    distc = fabs(s->col - n->col); if (distc > GLOBAL_Nra->cols/2){ distc -= GLOBAL_Nra->cols;}
	    dist = distr*distr + distc*distc;
	    af = (rand01 < P_AMPA)*exp(-dist/(double)AMPA_DIA)/P_AMPA/PI/(double)AMPA_DIA;
	    nf = (rand01 < P_NMDA)*exp(-dist/(double)NMDA_DIA)/P_NMDA/PI/(double)NMDA_DIA;
	    *A = CS_ESEC_A*af; *N = CS_ESEC_N*nf; *G = 0; connected=1;}}
	if (Amax!=NULL){ *Amax=CS_ESEC_A/P_AMPA/PI/(double)AMPA_DIA;} if (Nmax!=NULL){ *Nmax=CS_ESEC_N/P_NMDA/PI/(double)NMDA_DIA;} if (Gmax!=NULL){ *Gmax=0;} break;
      case 15: /* excitatory simple -> excitatory simple */
	newnum = minimum(SPARSEHIT_ESES,(double)SPARSE_ESES/s->homogenizing_knob);
	if (rand01 < (double)SPARSE_ESES/newnum){
	  if (periodize((s->sparse_out-n->sparse_in)%SPARSEHIT_ESES,0,SPARSEHIT_ESES)<newnum){
	    distr = fabs(s->row - n->row); if (distr > GLOBAL_Nra->rows/2){ distr -= GLOBAL_Nra->rows;}
	    distc = fabs(s->col - n->col); if (distc > GLOBAL_Nra->cols/2){ distc -= GLOBAL_Nra->cols;}
	    dist = distr*distr + distc*distc;
	    af = (rand01 < P_AMPA)*exp(-dist/(double)AMPA_DIA)/P_AMPA/PI/(double)AMPA_DIA;
	    nf = (rand01 < P_NMDA)*exp(-dist/(double)NMDA_DIA)/P_NMDA/PI/(double)NMDA_DIA;
	    *A = CS_ESES_A*af; *N = CS_ESES_N*nf; *G = 0; connected=1;}}
	if (Amax!=NULL){ *Amax=CS_ESES_A/P_AMPA/PI/(double)AMPA_DIA;} if (Nmax!=NULL){ *Nmax=CS_ESES_N/P_NMDA/PI/(double)NMDA_DIA;} if (Gmax!=NULL){ *Gmax=0;} break;
      default: printf(" %% error, wrong types (%d,%d)->(%d,%d) in slink\n",s->type,s->sox,n->type,n->sox); break;}}
    else if (s->homogenizing_knob==1 && s->sparse_link!=NULL){ 
      if (llitemaddorfind(0,s->sparse_link,(void *)n,&void_compare)!=NULL){
	connected=1;
	switch ((int)*(s->t2s) * 4 + (int)*(n->t2s) * 1){
	case 0: /* inhibitory complex -> inhibitory complex */  *A = 0; *N = 0; *G = CS_ICIC_G*(rand01<P_GABA)/P_GABA; 
	  if (Amax!=NULL){ *Amax=0;} if (Nmax!=NULL){ *Nmax=0;} if (Gmax!=NULL){ *Gmax=CS_ICIC_G/P_GABA;} break;
	case 1: /* inhibitory complex -> inhibitory simple */   *A = 0; *N = 0; *G = CS_ICIS_G*(rand01<P_GABA)/P_GABA; 
	  if (Amax!=NULL){ *Amax=0;} if (Nmax!=NULL){ *Nmax=0;} if (Gmax!=NULL){ *Gmax=CS_ICIS_G/P_GABA;} break;
	case 2: /* inhibitory complex -> excitatory complex */  *A = 0; *N = 0; *G = CS_ICEC_G*(rand01<P_GABA)/P_GABA; 
	  if (Amax!=NULL){ *Amax=0;} if (Nmax!=NULL){ *Nmax=0;} if (Gmax!=NULL){ *Gmax=CS_ICEC_G/P_GABA;} break;
	case 3: /* inhibitory complex -> excitatory simple */   *A = 0; *N = 0; *G = CS_ICES_G*(rand01<P_GABA)/P_GABA; 
	  if (Amax!=NULL){ *Amax=0;} if (Nmax!=NULL){ *Nmax=0;} if (Gmax!=NULL){ *Gmax=CS_ICES_G/P_GABA;} break;
	case 4: /* inhibitory simple -> inhibitory complex */   *A = 0; *N = 0; *G = CS_ISIC_G*(rand01<P_GABA)/P_GABA; 
	  if (Amax!=NULL){ *Amax=0;} if (Nmax!=NULL){ *Nmax=0;} if (Gmax!=NULL){ *Gmax=CS_ISIC_G/P_GABA;} break;
	case 5: /* inhibitory simple -> inhibitory simple */    *A = 0; *N = 0; *G = CS_ISIS_G*(rand01<P_GABA)/P_GABA; 
	  if (Amax!=NULL){ *Amax=0;} if (Nmax!=NULL){ *Nmax=0;} if (Gmax!=NULL){ *Gmax=CS_ISIS_G/P_GABA;} break;
	case 6: /* inhibitory simple -> excitatory complex */   *A = 0; *N = 0; *G = CS_ISEC_G*(rand01<P_GABA)/P_GABA; 
	  if (Amax!=NULL){ *Amax=0;} if (Nmax!=NULL){ *Nmax=0;} if (Gmax!=NULL){ *Gmax=CS_ISEC_G/P_GABA;} break;
	case 7: /* inhibitory simple -> excitatory simple */    *A = 0; *N = 0; *G = CS_ISES_G*(rand01<P_GABA)/P_GABA; 
	  if (Amax!=NULL){ *Amax=0;} if (Nmax!=NULL){ *Nmax=0;} if (Gmax!=NULL){ *Gmax=CS_ISES_G/P_GABA;} break;
	case 8: /* excitatory complex -> inhibitory complex */  
	  *A = CS_ECIC_A*(rand01<P_AMPA)/P_AMPA; *N = CS_ECIC_N*(rand01<P_NMDA)/P_NMDA; *G = 0; 
	  if (Amax!=NULL){ *Amax=CS_ECIC_A/P_AMPA;} if (Nmax!=NULL){ *Nmax=CS_ECIC_N/P_NMDA;} if (Gmax!=NULL){ *Gmax=0;} break;
	case 9: /* excitatory complex -> inhibitory simple */   
	  *A = CS_ECIS_A*(rand01<P_AMPA)/P_AMPA; *N = CS_ECIS_N*(rand01<P_NMDA)/P_NMDA; *G = 0; 
	  if (Amax!=NULL){ *Amax=CS_ECIS_A/P_AMPA;} if (Nmax!=NULL){ *Nmax=CS_ECIS_N/P_NMDA;} if (Gmax!=NULL){ *Gmax=0;} break;
	case 10: /* excitatory complex -> excitatory complex */ 
	  *A = CS_ECEC_A*(rand01<P_AMPA)/P_AMPA; *N = CS_ECEC_N*(rand01<P_NMDA)/P_NMDA; *G = 0; 
	  if (Amax!=NULL){ *Amax=CS_ECEC_A/P_AMPA;} if (Nmax!=NULL){ *Nmax=CS_ECEC_N/P_NMDA;} if (Gmax!=NULL){ *Gmax=0;} break;
	case 11: /* excitatory complex -> excitatory simple */  
	  *A = CS_ECES_A*(rand01<P_AMPA)/P_AMPA; *N = CS_ECES_N*(rand01<P_NMDA)/P_NMDA; *G = 0; 
	  if (Amax!=NULL){ *Amax=CS_ECES_A/P_AMPA;} if (Nmax!=NULL){ *Nmax=CS_ECES_N/P_NMDA;} if (Gmax!=NULL){ *Gmax=0;} break;
	case 12: /* excitatory simple -> inhibitory complex */  
	  *A = CS_ESIC_A*(rand01<P_AMPA)/P_AMPA; *N = CS_ESIC_N*(rand01<P_NMDA)/P_NMDA; *G = 0; 
	  if (Amax!=NULL){ *Amax=CS_ESIC_A/P_AMPA;} if (Nmax!=NULL){ *Nmax=CS_ESIC_N/P_NMDA;} if (Gmax!=NULL){ *Gmax=0;} break;
	case 13: /* excitatory simple -> inhibitory simple */   
	  *A = CS_ESIS_A*(rand01<P_AMPA)/P_AMPA; *N = CS_ESIS_N*(rand01<P_NMDA)/P_NMDA; *G = 0; 
	  if (Amax!=NULL){ *Amax=CS_ESIS_A/P_AMPA;} if (Nmax!=NULL){ *Nmax=CS_ESIS_N/P_NMDA;} if (Gmax!=NULL){ *Gmax=0;} break;
	case 14: /* excitatory simple -> excitatory complex */  
	  *A = CS_ESEC_A*(rand01<P_AMPA)/P_AMPA; *N = CS_ESEC_N*(rand01<P_NMDA)/P_NMDA; *G = 0; 
	  if (Amax!=NULL){ *Amax=CS_ESEC_A/P_AMPA;} if (Nmax!=NULL){ *Nmax=CS_ESEC_N/P_NMDA;} if (Gmax!=NULL){ *Gmax=0;} break;
	case 15: /* excitatory simple -> excitatory simple */   
	  *A = CS_ESES_A*(rand01<P_AMPA)/P_AMPA; *N = CS_ESES_N*(rand01<P_NMDA)/P_NMDA; *G = 0; 
	  if (Amax!=NULL){ *Amax=CS_ESES_A/P_AMPA;} if (Nmax!=NULL){ *Nmax=CS_ESES_N/P_NMDA;} if (Gmax!=NULL){ *Gmax=0;} break;
	default: printf(" %% error, wrong types (%d,%d)->(%d,%d) in slink\n",s->type,s->sox,n->type,n->sox); break;}}}
    else if (s->homogenizing_knob<1 && s->sparse_link!=NULL){ 
      /* p = (double)llitemlength(s->sparse_link)/(double)SYSTEM_ROW_SIZE/(double)SYSTEM_COL_SIZE)
	 s is connected with probability 
	 p*(1-s->homogenizing_knob)
	 to any neuron not within s->sparse_link, and with probability
	 p + (1-p)*s->homogenizing_knob
	 to any neuron within s->sparse_link */
      /* this way we have l*(l/s+(1-l/s)*h) + (s-l)*(l/s*(1-h)) = ll/s + lh - llh/s + l(1-h) - ll/s(1-h) = l average connections */
      newnum = minimum(1,(double)llitemlength(s->sparse_link)/(double)SYSTEM_ROW_SIZE*(double)SYSTEM_COL_SIZE);
      if (llitemaddorfind(0,s->sparse_link,(void *)n,&void_compare)!=NULL){
	connected = rand01 < (newnum + (1-newnum)*s->homogenizing_knob);}
      else{ connected = rand01 < newnum*(1-s->homogenizing_knob);}
      if (connected){
	switch ((int)*(s->t2s) * 4 + (int)*(n->t2s) * 1){
	case 0: /* inhibitory complex -> inhibitory complex */  *A = 0; *N = 0; *G = CS_ICIC_G*(rand01<P_GABA)/P_GABA; 
	  if (Amax!=NULL){ *Amax=0;} if (Nmax!=NULL){ *Nmax=0;} if (Gmax!=NULL){ *Gmax=CS_ICIC_G/P_GABA;} break;
	case 1: /* inhibitory complex -> inhibitory simple */   *A = 0; *N = 0; *G = CS_ICIS_G*(rand01<P_GABA)/P_GABA; 
	  if (Amax!=NULL){ *Amax=0;} if (Nmax!=NULL){ *Nmax=0;} if (Gmax!=NULL){ *Gmax=CS_ICIS_G/P_GABA;} break;
	case 2: /* inhibitory complex -> excitatory complex */  *A = 0; *N = 0; *G = CS_ICEC_G*(rand01<P_GABA)/P_GABA; 
	  if (Amax!=NULL){ *Amax=0;} if (Nmax!=NULL){ *Nmax=0;} if (Gmax!=NULL){ *Gmax=CS_ICEC_G/P_GABA;} break;
	case 3: /* inhibitory complex -> excitatory simple */   *A = 0; *N = 0; *G = CS_ICES_G*(rand01<P_GABA)/P_GABA; 
	  if (Amax!=NULL){ *Amax=0;} if (Nmax!=NULL){ *Nmax=0;} if (Gmax!=NULL){ *Gmax=CS_ICES_G/P_GABA;} break;
	case 4: /* inhibitory simple -> inhibitory complex */   *A = 0; *N = 0; *G = CS_ISIC_G*(rand01<P_GABA)/P_GABA; 
	  if (Amax!=NULL){ *Amax=0;} if (Nmax!=NULL){ *Nmax=0;} if (Gmax!=NULL){ *Gmax=CS_ISIC_G/P_GABA;} break;
	case 5: /* inhibitory simple -> inhibitory simple */    *A = 0; *N = 0; *G = CS_ISIS_G*(rand01<P_GABA)/P_GABA; 
	  if (Amax!=NULL){ *Amax=0;} if (Nmax!=NULL){ *Nmax=0;} if (Gmax!=NULL){ *Gmax=CS_ISIS_G/P_GABA;} break;
	case 6: /* inhibitory simple -> excitatory complex */   *A = 0; *N = 0; *G = CS_ISEC_G*(rand01<P_GABA)/P_GABA; 
	  if (Amax!=NULL){ *Amax=0;} if (Nmax!=NULL){ *Nmax=0;} if (Gmax!=NULL){ *Gmax=CS_ISEC_G/P_GABA;} break;
	case 7: /* inhibitory simple -> excitatory simple */    *A = 0; *N = 0; *G = CS_ISES_G*(rand01<P_GABA)/P_GABA; 
	  if (Amax!=NULL){ *Amax=0;} if (Nmax!=NULL){ *Nmax=0;} if (Gmax!=NULL){ *Gmax=CS_ISES_G/P_GABA;} break;
	case 8: /* excitatory complex -> inhibitory complex */  
	  *A = CS_ECIC_A*(rand01<P_AMPA)/P_AMPA; *N = CS_ECIC_N*(rand01<P_NMDA)/P_NMDA; *G = 0; 
	  if (Amax!=NULL){ *Amax=CS_ECIC_A/P_AMPA;} if (Nmax!=NULL){ *Nmax=CS_ECIC_N/P_NMDA;} if (Gmax!=NULL){ *Gmax=0;} break;
	case 9: /* excitatory complex -> inhibitory simple */   
	  *A = CS_ECIS_A*(rand01<P_AMPA)/P_AMPA; *N = CS_ECIS_N*(rand01<P_NMDA)/P_NMDA; *G = 0; 
	  if (Amax!=NULL){ *Amax=CS_ECIS_A/P_AMPA;} if (Nmax!=NULL){ *Nmax=CS_ECIS_N/P_NMDA;} if (Gmax!=NULL){ *Gmax=0;} break;
	case 10: /* excitatory complex -> excitatory complex */ 
	  *A = CS_ECEC_A*(rand01<P_AMPA)/P_AMPA; *N = CS_ECEC_N*(rand01<P_NMDA)/P_NMDA; *G = 0; 
	  if (Amax!=NULL){ *Amax=CS_ECEC_A/P_AMPA;} if (Nmax!=NULL){ *Nmax=CS_ECEC_N/P_NMDA;} if (Gmax!=NULL){ *Gmax=0;} break;
	case 11: /* excitatory complex -> excitatory simple */  
	  *A = CS_ECES_A*(rand01<P_AMPA)/P_AMPA; *N = CS_ECES_N*(rand01<P_NMDA)/P_NMDA; *G = 0; 
	  if (Amax!=NULL){ *Amax=CS_ECES_A/P_AMPA;} if (Nmax!=NULL){ *Nmax=CS_ECES_N/P_NMDA;} if (Gmax!=NULL){ *Gmax=0;} break;
	case 12: /* excitatory simple -> inhibitory complex */  
	  *A = CS_ESIC_A*(rand01<P_AMPA)/P_AMPA; *N = CS_ESIC_N*(rand01<P_NMDA)/P_NMDA; *G = 0; 
	  if (Amax!=NULL){ *Amax=CS_ESIC_A/P_AMPA;} if (Nmax!=NULL){ *Nmax=CS_ESIC_N/P_NMDA;} if (Gmax!=NULL){ *Gmax=0;} break;
	case 13: /* excitatory simple -> inhibitory simple */   
	  *A = CS_ESIS_A*(rand01<P_AMPA)/P_AMPA; *N = CS_ESIS_N*(rand01<P_NMDA)/P_NMDA; *G = 0; 
	  if (Amax!=NULL){ *Amax=CS_ESIS_A/P_AMPA;} if (Nmax!=NULL){ *Nmax=CS_ESIS_N/P_NMDA;} if (Gmax!=NULL){ *Gmax=0;} break;
	case 14: /* excitatory simple -> excitatory complex */  
	  *A = CS_ESEC_A*(rand01<P_AMPA)/P_AMPA; *N = CS_ESEC_N*(rand01<P_NMDA)/P_NMDA; *G = 0; 
	  if (Amax!=NULL){ *Amax=CS_ESEC_A/P_AMPA;} if (Nmax!=NULL){ *Nmax=CS_ESEC_N/P_NMDA;} if (Gmax!=NULL){ *Gmax=0;} break;
	case 15: /* excitatory simple -> excitatory simple */   
	  *A = CS_ESES_A*(rand01<P_AMPA)/P_AMPA; *N = CS_ESES_N*(rand01<P_NMDA)/P_NMDA; *G = 0; 
	  if (Amax!=NULL){ *Amax=CS_ESES_A/P_AMPA;} if (Nmax!=NULL){ *Nmax=CS_ESES_N/P_NMDA;} if (Gmax!=NULL){ *Gmax=0;} break;
	default: printf(" %% error, wrong types (%d,%d)->(%d,%d) in slink\n",s->type,s->sox,n->type,n->sox); break;}}
      else /* if (!connected) */{ 
	switch ((int)*(s->t2s) * 4 + (int)*(n->t2s) * 1){
	case 0: /* inhibitory complex -> inhibitory complex */  *A = 0; *N = 0; *G = 0;
	  if (Amax!=NULL){ *Amax=0;} if (Nmax!=NULL){ *Nmax=0;} if (Gmax!=NULL){ *Gmax=CS_ICIC_G/P_GABA;} break;
	case 1: /* inhibitory complex -> inhibitory simple */   *A = 0; *N = 0; *G = 0;
	  if (Amax!=NULL){ *Amax=0;} if (Nmax!=NULL){ *Nmax=0;} if (Gmax!=NULL){ *Gmax=CS_ICIS_G/P_GABA;} break;
	case 2: /* inhibitory complex -> excitatory complex */  *A = 0; *N = 0; *G = 0;
	  if (Amax!=NULL){ *Amax=0;} if (Nmax!=NULL){ *Nmax=0;} if (Gmax!=NULL){ *Gmax=CS_ICEC_G/P_GABA;} break;
	case 3: /* inhibitory complex -> excitatory simple */   *A = 0; *N = 0; *G = 0;
	  if (Amax!=NULL){ *Amax=0;} if (Nmax!=NULL){ *Nmax=0;} if (Gmax!=NULL){ *Gmax=CS_ICES_G/P_GABA;} break;
	case 4: /* inhibitory simple -> inhibitory complex */   *A = 0; *N = 0; *G = 0;
	  if (Amax!=NULL){ *Amax=0;} if (Nmax!=NULL){ *Nmax=0;} if (Gmax!=NULL){ *Gmax=CS_ISIC_G/P_GABA;} break;
	case 5: /* inhibitory simple -> inhibitory simple */    *A = 0; *N = 0; *G = 0;
	  if (Amax!=NULL){ *Amax=0;} if (Nmax!=NULL){ *Nmax=0;} if (Gmax!=NULL){ *Gmax=CS_ISIS_G/P_GABA;} break;
	case 6: /* inhibitory simple -> excitatory complex */   *A = 0; *N = 0; *G = 0;
	  if (Amax!=NULL){ *Amax=0;} if (Nmax!=NULL){ *Nmax=0;} if (Gmax!=NULL){ *Gmax=CS_ISEC_G/P_GABA;} break;
	case 7: /* inhibitory simple -> excitatory simple */    *A = 0; *N = 0; *G = 0;
	  if (Amax!=NULL){ *Amax=0;} if (Nmax!=NULL){ *Nmax=0;} if (Gmax!=NULL){ *Gmax=CS_ISES_G/P_GABA;} break;
	case 8: /* excitatory complex -> inhibitory complex */  *A = 0; *N = 0; *G = 0;
	  if (Amax!=NULL){ *Amax=CS_ECIC_A/P_AMPA;} if (Nmax!=NULL){ *Nmax=CS_ECIC_N/P_NMDA;} if (Gmax!=NULL){ *Gmax=0;} break;
	case 9: /* excitatory complex -> inhibitory simple */   *A = 0; *N = 0; *G = 0;
	  if (Amax!=NULL){ *Amax=CS_ECIS_A/P_AMPA;} if (Nmax!=NULL){ *Nmax=CS_ECIS_N/P_NMDA;} if (Gmax!=NULL){ *Gmax=0;} break;
	case 10: /* excitatory complex -> excitatory complex */ *A = 0; *N = 0; *G = 0;
	  if (Amax!=NULL){ *Amax=CS_ECEC_A/P_AMPA;} if (Nmax!=NULL){ *Nmax=CS_ECEC_N/P_NMDA;} if (Gmax!=NULL){ *Gmax=0;} break;
	case 11: /* excitatory complex -> excitatory simple */  *A = 0; *N = 0; *G = 0;
	  if (Amax!=NULL){ *Amax=CS_ECES_A/P_AMPA;} if (Nmax!=NULL){ *Nmax=CS_ECES_N/P_NMDA;} if (Gmax!=NULL){ *Gmax=0;} break;
	case 12: /* excitatory simple -> inhibitory complex */  *A = 0; *N = 0; *G = 0;
	  if (Amax!=NULL){ *Amax=CS_ESIC_A/P_AMPA;} if (Nmax!=NULL){ *Nmax=CS_ESIC_N/P_NMDA;} if (Gmax!=NULL){ *Gmax=0;} break;
	case 13: /* excitatory simple -> inhibitory simple */   *A = 0; *N = 0; *G = 0;
	  if (Amax!=NULL){ *Amax=CS_ESIS_A/P_AMPA;} if (Nmax!=NULL){ *Nmax=CS_ESIS_N/P_NMDA;} if (Gmax!=NULL){ *Gmax=0;} break;
	case 14: /* excitatory simple -> excitatory complex */  *A = 0; *N = 0; *G = 0;
	  if (Amax!=NULL){ *Amax=CS_ESEC_A/P_AMPA;} if (Nmax!=NULL){ *Nmax=CS_ESEC_N/P_NMDA;} if (Gmax!=NULL){ *Gmax=0;} break;
	case 15: /* excitatory simple -> excitatory simple */   *A = 0; *N = 0; *G = 0;
	  if (Amax!=NULL){ *Amax=CS_ESES_A/P_AMPA;} if (Nmax!=NULL){ *Nmax=CS_ESES_N/P_NMDA;} if (Gmax!=NULL){ *Gmax=0;} break;
	default: printf(" %% error, wrong types (%d,%d)->(%d,%d) in slink\n",s->type,s->sox,n->type,n->sox); break;}}}}
  return connected;
}

double ifrhs(double V,double sA,double sN,double sG,double dt)
{ 
  return (-CONDUCTANCE_LK*(V-VOLTAGE_REST) - (sA*exp(-dt/TAU_AMPA)+sN*exp(-dt/TAU_NMDA))*(V-VOLTAGE_EX) - (sG*exp(-dt/TAU_GABA))*(V-VOLTAGE_IN));
}

double eifrhs(double V,double sA,double sN,double sG,double dt)
{ 
  return (-CONDUCTANCE_LK*(V-VOLTAGE_REST) - (sA*exp(-dt/TAU_AMPA)+sN*exp(-dt/TAU_NMDA))*(V-VOLTAGE_EX) - (sG*exp(-dt/TAU_GABA))*(V-VOLTAGE_IN) + CONDUCTANCE_LK*VOLTAGE_DELTAT*exp((V-VOLTAGE_TAKEOFF)/VOLTAGE_DELTAT));
}

void slavecalc(double sa,double sn,double sg,double *gs,double *vs,double *vsprime)
{
  /* given conductances sa,sn,sg that decay exponentially, we compute the 
     slaving conductance gs
     slaving voltage vs
     slaving voltage prime vsprime */
  double a=0,b=0,c=0,d=0;
  a = CONDUCTANCE_LK + sa + sn + sg;
  b = -sa/TAU_AMPA - sn/TAU_NMDA - sg/TAU_GABA;
  c = CONDUCTANCE_LK*VOLTAGE_REST + (sa+sn)*VOLTAGE_EX + (sg)*VOLTAGE_IN;
  d = -(sa/TAU_AMPA+sn/TAU_NMDA)*VOLTAGE_EX - (sg/TAU_GABA)*VOLTAGE_IN;
  if (gs!=NULL){ *gs = a;}
  if (vs!=NULL){ *vs = c/a;}
  if (vsprime!=NULL){ *vsprime = (a*d - c*b)/(a*a);}
}

void gizmointegrate2(struct neuron *n,double previous_spiketime,double next_spiketime,double cA,double cN,double cG)
{
  /* this carries out the stable 2nd (not 3rd) order integration over spike-free time interval d, 
     assuming that the n->?1 hold initial conditions and the n->?2 are available for storage */
  int verbose=0;//(n->row*n->col==1);
  double d=0,eA=0,eN=0,eG=0,k1=0,if1=0,if2=0,if3=0,vs1=0,vs2=0,vs3=0,vsp1=0,vsp2=0,vsp3=0;
  double y1=0,v2=0,y2=0;
  if (verbose){ printf("start at *(n->sA)=%f,*(n->sN)=%f,*(n->sG)=%f,*(n->V)=%f\n",n->sA1,n->sN1,n->sG1,n->sV2);}
  if (HANSHELLEY_FLAG){
    /* perform rk instead */
    if (n->spikenext >= next_spiketime){
      d = next_spiketime - previous_spiketime;
      n->sA1 *= (1-d/TAU_AMPA); n->sA1 += cA;
      n->sN1 *= (1-d/TAU_NMDA); n->sN1 += cN;
      n->sG1 *= (1-d/TAU_GABA); n->sG1 += cG;
      n->sV1 = VOLTAGE_RESET;
      n->sV2 = VOLTAGE_RESET;}
    else{
      if (n->spikenext >= previous_spiketime && n->spikenext < next_spiketime){
	d = n->spikenext - previous_spiketime;
	n->sA1 *= (1-d/TAU_AMPA);
	n->sN1 *= (1-d/TAU_NMDA);
	n->sG1 *= (1-d/TAU_GABA);}
      n->sA2 = 0; n->sN2 = 0; n->sG2 = 0;
      d = next_spiketime - maximum(previous_spiketime,n->spikenext);
      if (verbose){ printf("    HANSHELLEY integrating with d=%f\n",d);}
      if (verbose){ printf("    we start with %f,%f,%f,%f,%f,%f\n",n->sA1,n->sA2,n->sN1,n->sN2,n->sG1,n->sG2);}
      k1 = (CONDUCTANCE_LK+n->sA1+n->sN1+n->sG1)*n->sV2 + (VOLTAGE_REST*CONDUCTANCE_LK + VOLTAGE_EX*(n->sA1+n->sN1) + VOLTAGE_IN*(n->sG1));
      eA = (1-d/TAU_AMPA); eN = (1-d/TAU_NMDA); eG = (1-d/TAU_GABA);
      n->sV1 = n->sV2;
      y1 = ifrhs(n->sV1,n->sA1,n->sN1,n->sG1,0);
      v2 = n->sV1 + 2.0*d/3.0*y1;
      y2 = ifrhs(v2,n->sA1,n->sN1,n->sG1,2.0*d/3.0);
      n->sV2 = n->sV1 + d/4.0*(y1+3*y2);
      n->sA1 *= eA; n->sN1 *= eN; n->sG1 *= eG;
      n->sA1 += cA; n->sN1 += cN; n->sG1 += cG;
      if (n->sV1 < VOLTAGE_THRESHOLD && n->sV2 >= VOLTAGE_THRESHOLD && !n->spiketime_guess_flag){ 
	n->spiketime_guess = maximum(previous_spiketime,n->spikenext) + quadrootfinder(n->sV1,k1,n->sV2,VOLTAGE_THRESHOLD,d);
	n->spiketime_guess_flag = 1;}}
    if (verbose){ printf("  now *(n->sA)=%f,*(n->sN)=%f,*(n->sG)=%f,*(n->V)=%f\n",n->sA1,n->sN1,n->sG1,n->sV2);}}
  else /* if (!HANSHELLEY_FLAG) */{
    if (n->spikenext >= next_spiketime){
      d = next_spiketime - previous_spiketime;
      n->sA1 *= exp(-d/TAU_AMPA); n->sA1 += cA;
      n->sN1 *= exp(-d/TAU_NMDA); n->sN1 += cN;
      n->sG1 *= exp(-d/TAU_GABA); n->sG1 += cG;
      n->sV1 = VOLTAGE_RESET;
      n->sV2 = VOLTAGE_RESET;}
    else{
      if (n->spikenext >= previous_spiketime && n->spikenext < next_spiketime){
	d = n->spikenext - previous_spiketime;
	n->sA1 *= exp(-d/TAU_AMPA);
	n->sN1 *= exp(-d/TAU_NMDA);
	n->sG1 *= exp(-d/TAU_GABA);}
      n->sA2 = 0; n->sN2 = 0; n->sG2 = 0;
      d = next_spiketime - maximum(previous_spiketime,n->spikenext);
      if (verbose){ printf("    integrating with d=%f\n",d);}
      if (verbose){ printf("    we start with %f,%f,%f,%f,%f,%f\n",n->sA1,n->sA2,n->sN1,n->sN2,n->sG1,n->sG2);}
      if (verbose){ printf("    we have %f,%f,%f,%f,%f,%f,%f,%f,%f\n",if1,if2,if3,vs1,vs2,vs3,vsp1,vsp2,vsp3);}
      k1 = (CONDUCTANCE_LK+n->sA1+n->sN1+n->sG1)*n->sV2 + (VOLTAGE_REST*CONDUCTANCE_LK + VOLTAGE_EX*(n->sA1+n->sN1) + VOLTAGE_IN*(n->sG1));
      eA = exp(-d/2/TAU_AMPA); eN = exp(-d/2/TAU_NMDA); eG = exp(-d/2/TAU_GABA);
      if1 = 0*CONDUCTANCE_LK+n->sA2+n->sN2+n->sG2;
      slavecalc(n->sA1,n->sN1,n->sG1,NULL,&vs1,&vsp1);
      n->sA2 = TAU_AMPA*n->sA1*(1 - eA);
      n->sN2 = TAU_NMDA*n->sN1*(1 - eN);
      n->sG2 = TAU_GABA*n->sG1*(1 - eG);
      n->sA1 *= eA;
      n->sN1 *= eN;
      n->sG1 *= eG;
      if2 = d/2*CONDUCTANCE_LK+n->sA2+n->sN2+n->sG2;
      slavecalc(n->sA1,n->sN1,n->sG1,NULL,&vs2,&vsp2);
      n->sA2 += TAU_AMPA*n->sA1*(1 - eA);
      n->sN2 += TAU_NMDA*n->sN1*(1 - eN);
      n->sG2 += TAU_GABA*n->sG1*(1 - eG);
      n->sA1 *= eA;
      n->sN1 *= eN;
      n->sG1 *= eG;
      if3 = d*CONDUCTANCE_LK+n->sA2+n->sN2+n->sG2;
      slavecalc(n->sA1,n->sN1,n->sG1,NULL,&vs3,&vsp3);
      n->g = d/6.0*(exp(if1-if3)*vsp1 + 4*exp(if2-if3)*vsp2 + exp(if3-if3)*vsp3);
      n->sA1 += cA;
      n->sN1 += cN;
      n->sG1 += cG;
      n->sV1 = n->sV2;
      n->sV2 = ((n->sV1-vs1)*exp(-if3) - n->g) + vs3;
      if (verbose){ printf("    now we have %f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",if1,if2,n->sA2,if3,vs1,vs2,vs3,vsp1,vsp2,vsp3);}
      if (n->sV1 < VOLTAGE_THRESHOLD && n->sV2 >= VOLTAGE_THRESHOLD && !n->spiketime_guess_flag){ 
	n->spiketime_guess = maximum(previous_spiketime,n->spikenext) + quadrootfinder(n->sV1,k1,n->sV2,VOLTAGE_THRESHOLD,d);
	n->spiketime_guess_flag = 1;}}
    if (verbose){ printf("  now *(n->sA)=%f,*(n->sN)=%f,*(n->sG)=%f,*(n->V)=%f\n",n->sA1,n->sN1,n->sG1,n->sV2);}}
}

void gizmointegrateeif(struct neuron *n,double previous_spiketime,double next_spiketime,double cA,double cN,double cG)
{
 /* this carries out the stable 2nd (not 3rd) order integration over spike-free time interval d */
  int verbose=0;//(n->row*n->col==1);
  double d=0,eA=0,eN=0,eG=0,k1=0,if1=0,if2=0,if3=0,vs1=0,vs2=0,vs3=0,vsp1=0,vsp2=0,vsp3=0;
  double Wnew=0,deltaw=0;
  int iteration=0,iteration_max=10;
  double epsilon_small=0.01,epsilon_big=10000000;
  if (verbose){ printf("%% Starting with n->sA1=%f,n->sN1=%f,n->sG1=%f,n->sV2=%f\n",n->sA1,n->sN1,n->sG1,n->sV2);}
  if (n->spikenext >= next_spiketime){
    d = next_spiketime - previous_spiketime;
    n->sA1 *= exp(-d/TAU_AMPA); n->sA1 += cA;
    n->sN1 *= exp(-d/TAU_NMDA); n->sN1 += cN;
    n->sG1 *= exp(-d/TAU_GABA); n->sG1 += cG;
    n->sV1 = VOLTAGE_RESET;
    n->sV2 = VOLTAGE_RESET;}
  else{
    if (n->spikenext >= previous_spiketime && n->spikenext < next_spiketime){
      d = n->spikenext - previous_spiketime;
      n->sA1 *= exp(-d/TAU_AMPA);
      n->sN1 *= exp(-d/TAU_NMDA);
      n->sG1 *= exp(-d/TAU_GABA);}
    d = next_spiketime - maximum(previous_spiketime,n->spikenext);
    eA = exp(-d/TAU_AMPA); eN = exp(-d/TAU_NMDA); eG = exp(-d/TAU_GABA);
    n->sV1 = n->sV2;
    k1 = eifrhs(n->sV1,n->sA1,n->sN1,n->sG1,0);
    d = next_spiketime - maximum(previous_spiketime,n->spikenext);
    n->sA2 = 0; n->sN2 = 0; n->sG2 = 0;
    eA = exp(-d/2/TAU_AMPA); eN = exp(-d/2/TAU_NMDA); eG = exp(-d/2/TAU_GABA);
    if1 = 0*CONDUCTANCE_LK+n->sA2+n->sN2+n->sG2;
    slavecalc(n->sA1,n->sN1,n->sG1,NULL,&vs1,&vsp1);
    n->sA2 += TAU_AMPA*n->sA1*(1 - eA); n->sN2 += TAU_NMDA*n->sN1*(1 - eN); n->sG2 += TAU_GABA*n->sG1*(1 - eG);
    n->sA1 *= eA; n->sN1 *= eN; n->sG1 *= eG;
    if2 = d/2*CONDUCTANCE_LK+n->sA2+n->sN2+n->sG2;
    slavecalc(n->sA1,n->sN1,n->sG1,NULL,&vs2,&vsp2);
    n->sA2 += TAU_AMPA*n->sA1*(1 - eA); n->sN2 += TAU_NMDA*n->sN1*(1 - eN); n->sG2 += TAU_GABA*n->sG1*(1 - eG);
    n->sA1 *= eA; n->sN1 *= eN; n->sG1 *= eG;
    if3 = d*CONDUCTANCE_LK+n->sA2+n->sN2+n->sG2;
    slavecalc(n->sA1,n->sN1,n->sG1,NULL,&vs3,&vsp3);
    n->g = exp(-if3)*(n->sV1-vs1) - d/6.0*(exp(if1-if3)*vsp1 + 4*exp(if2-if3)*vsp2 + exp(if3-if3)*vsp3);
    n->sA1 += cA; n->sN1 += cN; n->sG1 += cG;
    Wnew = n->g + d*CONDUCTANCE_LK*VOLTAGE_DELTAT*exp((n->g+vs3-VOLTAGE_TAKEOFF)/VOLTAGE_DELTAT); iteration = 0;
    do{
      deltaw = exp((Wnew + vs3 - VOLTAGE_TAKEOFF)/VOLTAGE_DELTAT);
      deltaw = (n->g + CONDUCTANCE_LK*VOLTAGE_DELTAT*d*deltaw - Wnew)/(1 - d*CONDUCTANCE_LK*deltaw);
      if (verbose){ printf("%% %% iteration with Wnew=%0.3f, deltaw=%0.3f\n",Wnew,deltaw);}
      Wnew += deltaw;
      deltaw = fabs(deltaw/(Wnew==0? 1 : Wnew));
      iteration++;}
    while (finite(deltaw) && fabs(deltaw) > epsilon_small && fabs(deltaw) < epsilon_big && iteration < iteration_max);
    if (finite(deltaw) && fabs(deltaw) <= epsilon_small){
      n->sV2 = Wnew + vs3;}
    else if (!finite(deltaw)){
      if (verbose){ printf("%% warning, deltaw not finite in gizmointegrateeif\n");} 
      n->sV2 = n->g + vs3 + d*CONDUCTANCE_LK*VOLTAGE_DELTAT*exp((n->g+vs3-VOLTAGE_TAKEOFF)/VOLTAGE_DELTAT);}
    else if (fabs(deltaw) >= epsilon_big || iteration >= iteration_max){
      if (verbose){ printf("%% warning, deltaw=%0.3f after %d iterations too big in gizmointegrateeif\n",deltaw,iteration);} 
      n->sV2 = n->g + vs3 + d*CONDUCTANCE_LK*VOLTAGE_DELTAT*exp((n->g+vs3-VOLTAGE_TAKEOFF)/VOLTAGE_DELTAT);}}
  if (finite(n->sV2) && n->sV1 < VOLTAGE_THRESHOLD_EIF && n->sV2 >= VOLTAGE_THRESHOLD_EIF){
    if (!n->spiketime_guess_flag){ 
      n->spiketime_guess = maximum(previous_spiketime,n->spikenext) + quadrootfinder(n->sV1,k1,n->sV2,VOLTAGE_THRESHOLD_EIF,d);
      n->spiketime_guess_flag = 1;}
    n->sV2 = VOLTAGE_RESET;}
  else if (!finite(n->sV2)){
    if (!n->spiketime_guess_flag){ 
      n->spiketime_guess = maximum(previous_spiketime,n->spikenext) + quadrootfinder(n->sV1,k1,n->sV1 + 200*(VOLTAGE_THRESHOLD_EIF-VOLTAGE_REST),VOLTAGE_THRESHOLD_EIF,d);
      n->spiketime_guess_flag = 1;}
    n->sV2 = VOLTAGE_RESET;}
  if (verbose){ printf("%% Ending with n->sA1=%f,n->sN1=%f,n->sG1=%f,n->sV2=%f\n",n->sA1,n->sN1,n->sG1,n->sV2);}
}

/* void gizmointegrateeif(struct neuron *n,double previous_spiketime,double next_spiketime,double cA,double cN,double cG) */
/* { */
/*   /\* this carries out the stable 2nd (not 3rd) order integration over spike-free time interval d *\/ */
/*   double d=0,y1=0,v2=0,y2=0,eA=0,eN=0,eG=0; */
/*   //if (verbose){ printf("start at *(n->sA)=%f,*(n->sN)=%f,*(n->sG)=%f,*(n->V)=%f\n",n->sA1,n->sN1,n->sG1,n->sV2);} */
/*   if (n->spikenext >= next_spiketime){ */
/*     d = next_spiketime - previous_spiketime; */
/*     n->sA1 *= exp(-d/TAU_AMPA); n->sA1 += cA; */
/*     n->sN1 *= exp(-d/TAU_NMDA); n->sN1 += cN; */
/*     n->sG1 *= exp(-d/TAU_GABA); n->sG1 += cG; */
/*     n->sV1 = VOLTAGE_RESET; */
/*     n->sV2 = VOLTAGE_RESET;} */
/*   else{ */
/*     if (n->spikenext >= previous_spiketime && n->spikenext < next_spiketime){ */
/*       d = n->spikenext - previous_spiketime; */
/*       n->sA1 *= exp(-d/TAU_AMPA); */
/*       n->sN1 *= exp(-d/TAU_NMDA); */
/*       n->sG1 *= exp(-d/TAU_GABA);} */
/*     d = next_spiketime - maximum(previous_spiketime,n->spikenext); */
/*     eA = exp(-d/TAU_AMPA); eN = exp(-d/TAU_NMDA); eG = exp(-d/TAU_GABA); */
/*     n->sV1 = n->sV2; */
/*     y1 = eifrhs(n->sV1,n->sA1,n->sN1,n->sG1,0); */
/*     v2 = n->sV1 + 2.0*d/3.0*y1; */
/*     y2 = eifrhs(v2,n->sA1,n->sN1,n->sG1,2.0*d/3.0); */
/*     n->sV2 = n->sV1 + d/4.0*(y1+3*y2); */
/*     n->sA1 *= eA; n->sN1 *= eN; n->sG1 *= eG; */
/*     n->sA1 += cA; n->sN1 += cN; n->sG1 += cG; */
/*     if (finite(n->sV2) && n->sV1 < VOLTAGE_THRESHOLD_EIF && n->sV2 >= VOLTAGE_THRESHOLD_EIF){ */
/*       if (!n->spiketime_guess_flag){  */
/* 	n->spiketime_guess = maximum(previous_spiketime,n->spikenext) + quadrootfinder(n->sV1,y1,n->sV2,VOLTAGE_THRESHOLD_EIF,d); */
/* 	//assert(n->spiketime_guess >= previous_spiketime && n->spiketime_guess <= next_spiketime); */
/* 	n->spiketime_guess_flag = 1;} */
/*       n->sV2 = VOLTAGE_RESET;} */
/*     else if (!finite(n->sV2)){ */
/*       if (!n->spiketime_guess_flag){  */
/* 	n->spiketime_guess = maximum(previous_spiketime,n->spikenext) + quadrootfinder(n->sV1,y1,n->sV1 + 200*(VOLTAGE_THRESHOLD_EIF-VOLTAGE_REST),VOLTAGE_THRESHOLD_EIF,d); */
/* 	//assert(n->spiketime_guess >= previous_spiketime && n->spiketime_guess <= next_spiketime); */
/* 	n->spiketime_guess_flag = 1;} */
/*       n->sV2 = VOLTAGE_RESET;}} */
/*   //if (n->row*n->col==1){ printf("V=%f..%f, spike=%f _ final\n",n->sV1,n->sV2,n->spiketime);} */
/* } */

void clumpcorrect(struct llist *L1,struct llist *L2,double t,double DT)
{
  /* We speed up integration by only accounting for spikes which actually connect *L1 and *L2 */
  //int verbose=0;
  int finished_with_input=0,number_of_spikes_left=0,perform_integration=0;
  double previous_integration_time=0,next_integration_time=0;
  double cA=0,cN=0,cG=0;
  struct neuron *n=NULL,*s=NULL;
  struct litem *lneuron=NULL,*lspike=NULL;
  double epsilon = 0.0000000001;
  lneuron = L1->first;
  while (lneuron != NULL){
    n = (struct neuron *) lneuron->item;
    //verbose = (n->row*n->col==1);
    //if (verbose){ printf("neuron (%d,%d) at time %f, dt %f\n",n->row,n->col,t,DT);}
    n->sA1 = *(n->sA);
    n->sN1 = *(n->sN);
    n->sG1 = *(n->sG);
    n->sA2 = 0;
    n->sN2 = 0;
    n->sG2 = 0;
    n->g = 0;
    n->sV1 = *(n->V);
    n->sV2 = *(n->V);
    n->spiketime_guess = n->spiketime;
    n->spiketime_guess_flag = 0;
    finished_with_input = !(n->spikeinput_flag && n->spikeinput_time >= t && n->spikeinput_time <= t+DT);
    number_of_spikes_left = L2->length;
    previous_integration_time = t;
    next_integration_time = t;
    perform_integration=0;
    lspike = L2->first;
    while (!finished_with_input || number_of_spikes_left>0 || previous_integration_time < t+DT-epsilon){
      //if (verbose){ printf(" fwi=%d,nos=%d,pit=%f<t+DT=%f\n",finished_with_input,number_of_spikes_left,previous_integration_time,t+DT);}
      if (lspike==NULL){
	//if (verbose){ printf("  lspike==NULL\n");}
	number_of_spikes_left=0;
	if (!finished_with_input && n->spikeinput_time >= previous_integration_time && n->spikeinput_time <= t+DT){
	  //if (verbose){ printf("   but still have input at pit=%f,%f,%f=t+DT\n",previous_integration_time,n->spikeinput_time,t+DT);}
	  finished_with_input=1;
	  next_integration_time = n->spikeinput_time;
	  perform_integration = ilink(n,&cA,&cN,&cG);}
	else /* (if finished_with_input || n->spikeinput_time < previous_integration_time || n->spikeinput_time > t+DT) */{
	  //if (verbose){ printf("   but no more input\n");}
	  finished_with_input=1;
	  cA=0;cN=0;cG=0;
	  next_integration_time = t+DT;
	  perform_integration=1;}}
      else if (lspike!=NULL){
	//if (verbose){ printf("  lspike exists\n");}
	s = (struct neuron *) lspike->item;
	if (!finished_with_input && n->spikeinput_time >= previous_integration_time && n->spikeinput_time <= minimum(t+DT,s->spiketime)){
	  //if (verbose){ printf("   but still have input at pit=%f,%f,%f=s->spk\n",previous_integration_time,n->spikeinput_time,s->spiketime);}
	  finished_with_input=1;
	  next_integration_time = n->spikeinput_time;
	  perform_integration = ilink(n,&cA,&cN,&cG);}
	else /* (if finished_with_input || n->spikeinput_time < previous_integration_time || n->spikeinput_time > minimum(t+DT,s->spiketime)) */{
	  //if (verbose){ printf("   no input before spiketime\n");}
	  if (s->spiketime >= previous_integration_time && s->spiketime <= t+DT){
	    //if (verbose){ printf("    pit=%f < spk=%f < t+DT=%f\n",previous_integration_time,s->spiketime,t+DT);}
	    next_integration_time = s->spiketime;
	    perform_integration = slink(s,n,&cA,&cN,&cG,NULL,NULL,NULL);
	    number_of_spikes_left -= 1;
	    lspike = lspike->child;}
	  else /* (if s->spiketime < previous_integration_time || s->spiketime > t+DT) */{
	    //if (verbose){ printf("    spiketime %f out of bounds pit=%f,t+DT=%f\n",s->spiketime,previous_integration_time,t+DT);}
	    number_of_spikes_left=0;
	    if (!finished_with_input && n->spikeinput_time >= previous_integration_time && n->spikeinput_time <= t+DT){
	      //if (verbose){ printf("     not yet finished with input pit=%f,%f,%f=t+DT\n",previous_integration_time,n->spikeinput_time,t+DT);}
	      finished_with_input=1;
	      next_integration_time = n->spikeinput_time;
	      perform_integration = ilink(n,&cA,&cN,&cG);}
	    else /* (if finished_with_input || n->spikeinput_time < previous_integration_time || n->spikeinput_time > t+DT) */{
	      //if (verbose){ printf("     no more input, we are exhausted\n");}
	      finished_with_input=1;
	      cA=0;cN=0;cG=0;
	      next_integration_time = t+DT;
	      perform_integration=1;}}}}
      if (perform_integration==1){ 
	//if (verbose){ printf("      -integrating %f %f\n",previous_integration_time,next_integration_time);}
	if (EIFvsIF>0){ gizmointegrateeif(n,previous_integration_time,next_integration_time,cA,cN,cG);}
	else{ gizmointegrate2(n,previous_integration_time,next_integration_time,cA,cN,cG);}
	previous_integration_time = next_integration_time;}}
    //if (verbose){ printf(" moving on to another neuron \n\n");}
    lneuron = lneuron->child;}
  lneuron = L1->first;
  while (lneuron != NULL){
    n = (struct neuron*) lneuron->item;
    assert(!(n->sV2 >= ((EIFvsIF>0) ? VOLTAGE_THRESHOLD_EIF : VOLTAGE_THRESHOLD) && (n->spiketime_guess < t || n->spiketime_guess > t+DT || n->spiketime_guess_flag==0)));
    if (n->spiketime_guess_flag){ n->spiketime = n->spiketime_guess;}
    else{ n->spiketime = t+2*DT;}
    lneuron = lneuron->child;}
}

void spikecorrect(struct llist *LL, double t, double DT)
{
  struct llist *Spikes=NULL;
  struct litem *Lc=NULL;
  Lc = LL->first;
  while (Lc!=NULL){
    Spikes = (struct llist *) Lc->item;
    clumpcorrect(Spikes,Spikes,t,DT);     
    llistsort(Spikes->first,Spikes->last,Spikes->length,&spiketime_compare);
    clumpcorrect(Spikes,Spikes,t,DT);
    /* set spiketimes */
    llistsort(Spikes->first,Spikes->last,Spikes->length,&spiketime_compare);
    Lc = Lc->child;}
}

void spikeconduct(struct llist *LL, double t, double DT,int *ES_nspk,int *IS_nspk,int *EC_nspk,int *IC_nspk)
{
  struct neuron *n=NULL;
  struct llist *Spikes=NULL;
  struct litem *Lc=NULL,*lneuron=NULL;
  *ES_nspk=0,*IS_nspk=0,*EC_nspk=0,*IC_nspk=0;
  Lc = LL->first;
  while (Lc!=NULL){
    Spikes = (struct llist *) Lc->item;
    lneuron = Spikes->first;
    /* evolve neurons */
    while (lneuron != NULL){
      n = (struct neuron *) lneuron->item;
      *(n->sA) = n->sA1;
      *(n->sN) = n->sN1;
      *(n->sG) = n->sG1;
      if (n->spiketime >= t && n->spiketime <= t+DT){ n->spikenext = n->spiketime + TAU_REF; n->sV2=VOLTAGE_RESET;}
      *(n->V) = n->sV2; 
      slavecalc(*(n->sA),*(n->sN),*(n->sG),NULL,n->VS,NULL);
      if (*(n->V)>= ((EIFvsIF>0) ? VOLTAGE_THRESHOLD_EIF : VOLTAGE_THRESHOLD)){ 
	n->warning = 2; printf(" %% warning! neuron (%d,%d) with v=%f in spikeconduct\n",n->row,n->col,*(n->V));}
      lneuron = lneuron->child;}
    /* update blocks and other data structures */
    lneuron = Spikes->first;
    while (lneuron != NULL){
      n = (struct neuron *) lneuron->item; 
      if (n->spiketime >= t && n->spiketime <= t+DT){ 
	n->warning = -2;
	n->spikelast=n->spiketime;
	switch ((int)*(n->t2s)){
	case 0: /* IC */ *IC_nspk += 1; break;
	case 1: /* IS */ *IS_nspk += 1; break;
	case 2: /* EC */ *EC_nspk += 1; break;
	case 3: /* ES */ *ES_nspk += 1; break;
	default: break;}}
      else{ n->warning = -1;}
      lneuron = lneuron->child;}
    Lc = Lc->child;}
}

double quadrootfinder(double VI, double k1, double VF, double VT, double DT)
{
  int i=0,imax=10;
  double a=0,b=0,c=0,e=0.00000000001;
  double r=0,dr=0;
  if (VI<= VT && VF >= VT){
    c = VI;
    b = k1;
    a = (VF - c - DT*b)/pow(DT,2);
    r = DT*(VT-VI)/(VF-VI);
    do{
      dr = -((a*r+b)*r+c-VT)/(2.0*a*r+b);
      r = r+dr;
      i++;}
    while (dr/r > e && i<imax);
    if (i>=imax || r<=0 || r>=DT){ 
      r = DT*(VT-VI)/(VF-VI);}}
  else{ r = 2*DT;}
  return r;
}

double cuberootfinder(double VI, double k1, double VF, double k3, double VT, double DT)
{
  double a=0,b=0,c=0,d=0;
  double r=0,dr=1.0;
  double e=0.00000000001;
  int i=0,imax=10;
  d = VI;
  c = k1;
  b = -(k3 - c)/DT + 3*(VF - d - c*DT)/pow(DT,2);
  a = (k3 - c)/pow(DT,2) - 2*(VF - d - c*DT)/pow(DT,3);
  r = quadrootfinder(VI,k1,VF,VT,DT);  
  if (r < 0 || r > DT){ r = DT*(VT-VI)/(VF-VI);}
  if (r < 0 || r > DT){ r = DT/2;}
  do{
    dr = -(((a*r+b)*r+c)*r+d-VT)/((3*a*r+2*b)*r+c);
    r = r+dr;
    i++;}
  while (dr/r > e && i < imax);
  if (i>=imax){ /* printf(" %% warning, dr/r=%f and i==imax=%d in cuberootfinder\n",dr/r,i); */ r = 2*DT;}
  if (r < 0){ r = 2*DT;}
  if (r > DT && VI < VT && VF > VT){ r = quadrootfinder(VI,k1,VF,VT,DT);}
  return r;
}

void gizmoconduct(struct llist **Sra , struct llist **Gra, double t ,double DT)
{
  int ar=0,ac=0,ar1=0,ar2=0,ar3=0,ac1=0,ac2=0,ac3=0;
  struct neuron *n=NULL;
  struct litem *lneuron=NULL;
  struct llist *Spikes=NULL,*S1=NULL,*S2=NULL,*S3=NULL,*S4=NULL,*S5=NULL,*S6=NULL,*S7=NULL,*S8=NULL,*S9=NULL,*Spikeadd=NULL;
  struct llist *Gizmos=NULL;
  int gnt=GLOBAL_NBLOCKS_TALL,gnw=GLOBAL_NBLOCKS_WIDE;
  for (ar=0;ar<gnt;ar++){ for (ac=0;ac<gnw;ac++){
    /* we examine nearby arbors 
       147
       258
       369
    */
    ar2 = ar; ar1 = periodize(ar2-1,0,gnt); ar3 = periodize(ar2+1,0,gnt);
    ac2 = ac; ac1 = periodize(ac2-1,0,gnw); ac3 = periodize(ac2+1,0,gnw);
    S1 = Sra[ar1 + ac1*gnt];
    S2 = Sra[ar2 + ac1*gnt];
    S3 = Sra[ar3 + ac1*gnt];
    S4 = Sra[ar1 + ac2*gnt];
    S5 = Sra[ar2 + ac2*gnt];
    S6 = Sra[ar3 + ac2*gnt];
    S7 = Sra[ar1 + ac3*gnt];
    S8 = Sra[ar2 + ac3*gnt];
    S9 = Sra[ar3 + ac3*gnt];
    Spikes = llistmake();
    if (S5!=NULL){ if (S5->length > 0){ Spikeadd = llistcopy(S5); lreconnect(Spikes,Spikeadd);}}
    if (S1!=NULL){ if (S1!=S5 && S1->length > 0){ Spikeadd = llistcopy(S1); lreconnect(Spikes,Spikeadd);}}
    if (S2!=NULL){ if (S2!=S5 && S2!=S1 && S2->length > 0){ Spikeadd = llistcopy(S2); lreconnect(Spikes,Spikeadd);}}
    if (S3!=NULL){ if (S3!=S5 && S3!=S1 && S3!=S2 && S3->length > 0){ Spikeadd = llistcopy(S3); lreconnect(Spikes,Spikeadd);}}
    if (S4!=NULL){ if (S4!=S5 && S4!=S1 && S4!=S2 && S4!=S3 && S4->length > 0){ Spikeadd = llistcopy(S4); lreconnect(Spikes,Spikeadd);}}
    if (S6!=NULL){ if (S6!=S5 && S6!=S1 && S6!=S2 && S6!=S3 && S6!=S4 && S6->length > 0){ Spikeadd = llistcopy(S6); lreconnect(Spikes,Spikeadd);}}
    if (S7!=NULL){ if (S7!=S5 && S7!=S1 && S7!=S2 && S7!=S3 && S7!=S4 && S7!=S6 && S7->length > 0){ Spikeadd = llistcopy(S7); lreconnect(Spikes,Spikeadd);}}
    if (S8!=NULL){ if (S8!=S5 && S8!=S1 && S8!=S2 && S8!=S3 && S8!=S4 && S8!=S6 && S8!=S7 && S8->length > 0){ Spikeadd = llistcopy(S8); lreconnect(Spikes,Spikeadd);}}
    if (S9!=NULL){ if (S9!=S5 && S9!=S1 && S9!=S2 && S9!=S3 && S9!=S4 && S9!=S6 && S9!=S7 && S9!=S8 && S9->length > 0){ Spikeadd = llistcopy(S9); lreconnect(Spikes,Spikeadd);}}
    llistsort(Spikes->first,Spikes->last,Spikes->length,&spiketime_compare);
    Gizmos = Gra[ar + ac*gnt];
    clumpcorrect(Gizmos,Spikes,t,DT);
    lneuron = Gizmos->first;
    while (lneuron!=NULL){
      n = (struct neuron *) lneuron->item;
      *(n->sA) = n->sA1;
      *(n->sN) = n->sN1;
      *(n->sG) = n->sG1;
      *(n->V) = n->sV2;
      slavecalc(*(n->sA),*(n->sN),*(n->sG),NULL,n->VS,NULL);
      lneuron = lneuron->child;}
    lneuron = Gizmos->first;
    while (lneuron!=NULL){
      n = (struct neuron *) lneuron->item;
      if (n->spiketime >= t && n->spiketime <= t+DT){ gizmospiked(n,t,DT); n->warning=1;} else{ n->warning=0;}
      if (*(n->V)>=((EIFvsIF>0) ? VOLTAGE_THRESHOLD_EIF : VOLTAGE_THRESHOLD)){ n->warning=2; printf(" %% warning! neuron (%d,%d) with v=%f in gizmoconduct\n",n->row,n->col,*(n->V));}
      lneuron = lneuron->child;}
    llisttfree(Spikes);}}
}

void gizmospiked(struct neuron *n,double t,double DT)
{
  /* special treatment for surprising spikes */
  int verbose=0;
  int ar=0,ac=0,ar2=0,ac2=0;
  double cA=0,cN=0,cG=0,d=0;
  struct neuron *n2=NULL;
  int rowbndry = maximum(minimum(GLOBAL_BLOCK_ROW_DIA,SYSTEM_ROW_SIZE/2),1);
  int colbndry = maximum(minimum(GLOBAL_BLOCK_COL_DIA,SYSTEM_COL_SIZE/2),1);
  if (verbose){ printf(" %% entering gizmospiked with rowbndry=%d,colbndry=%d\n",rowbndry,colbndry);}
  assert(n->spiketime >= t && n->spiketime <= t+DT);
  d = t+DT-n->spiketime;
  for (ar=1-rowbndry;ar<rowbndry;ar++){
    for (ac=1-colbndry;ac<colbndry;ac++){
      ar2 = periodize(ar+n->row,0,GLOBAL_Nra->rows);
      ac2 = periodize(ac+n->col,0,GLOBAL_Nra->cols);
      if (verbose){ printf(" %% trying to access relative (%d,%d) absolute (%d,%d)\n",ar,ac,ar2,ac2);}
      n2 = nget(GLOBAL_Nra,ar2,ac2);
      slink(n,n2,&cA,&cN,&cG,NULL,NULL,NULL);
      *(n2->sA) += exp(-d/TAU_AMPA)*cA;
      *(n2->sN) += exp(-d/TAU_NMDA)*cN;
      *(n2->sG) += exp(-d/TAU_GABA)*cG;}}
  n->spikenext = n->spiketime + TAU_REF;
  *(n->V) = VOLTAGE_RESET;
  slavecalc(*(n->sA),*(n->sN),*(n->sG),NULL,n->VS,NULL);
  n->spikelast=n->spiketime;
  if (verbose){ printf(" %% leaving gizmospiked!\n");}
}
