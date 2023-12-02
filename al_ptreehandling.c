
struct pnode * pnodemake(struct pnode *parent,struct region *region,double weight,double relevance)
{
  struct pnode *p=NULL;
  p = (struct pnode *) tmalloc(sizeof(struct pnode));
  p->parent = parent;
  p->region = region;
  p->weight = weight;
  p->relevance = relevance;
  p->childllitem = llitemmake();
  p->broodsize = 1;
  p->temp = NULL;
  return p;
}

void pnodetfree(void *vp)
{
  /* recursively frees pnode *p and its children */
  int verbose=0;
  struct pnode *p = (struct pnode *) vp;
  if (verbose){ 
    printf(" trying to free pnode %d with parent %d, region %d, and childllitem...\n",(int)p,(int)p->parent,p->region->label);
    llitemprintf(p->childllitem,&pnodeprintf_label);}
  llitemtfree(p->childllitem,&pnodetfree); 
  p->childllitem=NULL;
  tfree(p);p=NULL;
}

struct ptree * ptreemake(int nregions,double event_within,int event_threshold,int region_type,int nlegs,int legtime)
{
  /* update to use region data structure */
  struct ptree *p=NULL;
  int nr=0;
  p = (struct ptree *) tmalloc(sizeof(struct ptree));
  p->nregions = nregions;
  p->regionra = (struct region **) tcalloc(p->nregions,sizeof(struct region*));
  regionramake(p,event_within,event_threshold,region_type);
  p->nlegs=nlegs;
  p->legtime=legtime;
  p->length = p->nlegs*p->legtime+1;
  p->eventra = (struct llitem **) tmalloc(p->length*sizeof(struct llitem *));
  for (nr=0;nr<p->length;nr++){ p->eventra[nr]=llitemmake();}
  p->tab=0;
  p->update_every=1;
  p->update_last=GLOBAL_TI;
  p->total_time=0;
  p->gate_flag=1;
  ptreemakepospre(p);
  p->wh = histmake(1,+1,-1);
  p->rh = histmake(1,+1,-1);
  return p;
}

void ptreemakepospre(struct ptree *p)
{
  /* assumes p->pretree and p->postree are both NULL */
  int nr=0;
  p->pretree = llitemmake();
  for (nr=0;nr<p->nregions;nr++){
    if (llitemaddorfind(0,p->pretree,p->regionra[nr],&region2region_compare_label)!=NULL){ printf(" %% duplicate region label in ptreemakepospre\n");}
    else{ llitemaddorfind(1,p->pretree,pnodemake(NULL,p->regionra[nr],0,0),&pnode2pnode_compare_label);}}
  llitembalance(p->pretree); p->pretree = llitemclimb(p->pretree);
  p->postree = llitemmake();
  for (nr=0;nr<p->nregions;nr++){
    if (llitemaddorfind(0,p->postree,p->regionra[nr],&region2region_compare_label)!=NULL){ printf(" %% duplicate region label in ptreemakepospre\n");}
    else{ llitemaddorfind(1,p->postree,pnodemake(NULL,p->regionra[nr],0,0),&pnode2pnode_compare_label);}}
  llitembalance(p->postree); p->postree = llitemclimb(p->postree);
}

void ptreereset(struct ptree *p)
{
  /* destroys postree and pretree, but retains regionra and eventra */
  p->total_time=0;
  llitemtfree(p->pretree,&pnodetfree); p->pretree=NULL;
  llitemtfree(p->postree,&pnodetfree); p->postree=NULL;
  ptreemakepospre(p);
}

void ptreetfree(struct ptree *p)
{
  int nr=0;
  if (p->pretree!=NULL){ llitemtfree(p->pretree,&pnodetfree); p->pretree=NULL;}
  if (p->postree!=NULL){ llitemtfree(p->postree,&pnodetfree); p->postree=NULL;}
  if (p->eventra!=NULL){ 
    for (nr=0;nr<p->length;nr++){ if (p->eventra[nr]!=NULL){ llitemtfree(p->eventra[nr],NULL); p->eventra[nr]=NULL;}}
    tfree(p->eventra); p->eventra=NULL;}
  if (p->regionra!=NULL){ regionratfree(p); p->regionra=NULL;}
  if (p->wh!=NULL){ histtfree(p->wh); p->wh=NULL;}
  if (p->rh!=NULL){ histtfree(p->rh); p->rh=NULL;}
  if (p!=NULL){ tfree(p); p=NULL;}
}

int region2region_compare_last_event(void *v1,void *v2)
{
  struct region *r1 = (struct region *)v1;
  struct region *r2 = (struct region *)v2;
  return (r1->last_event > r2->last_event ? 1 : r1->last_event < r2->last_event ? -1 : 0);
}

int region2region_compare_label(void *v1,void *v2)
{
  struct region *r1 = (struct region *)v1;
  struct region *r2 = (struct region *)v2;
  return (r1->label > r2->label ? 1 : r1->label < r2->label ? -1 : 0);
}

int region2pnode_compare_label(void *v1,void *v2)
{
  struct region *r1 = (struct region *)v1;
  struct pnode *p2 = (struct pnode *)v2;
  return (r1->label > p2->region->label ? 1 : r1->label < p2->region->label ? -1 : 0);
}

int pnode2pnode_compare_label(void *v1,void *v2)
{
  struct pnode *p1 = (struct pnode *)v1;
  struct pnode *p2 = (struct pnode *)v2;
  return (p1->region->label > p2->region->label ? 1 : p1->region->label < p2->region->label ? -1 : 0);
}

int pnode2pnode_compare_relevance(void *v1,void *v2)
{
  struct pnode *p1 = (struct pnode *)v1;
  struct pnode *p2 = (struct pnode *)v2;
  return (p1->relevance > p2->relevance ? 1 : p1->relevance < p2->relevance ? -1 : 0);
}

int pnode2pnode_compare_weight(void *v1,void *v2)
{
  struct pnode *p1 = (struct pnode *)v1;
  struct pnode *p2 = (struct pnode *)v2;
  return (p1->weight > p2->weight ? 1 : p1->weight < p2->weight ? -1 : 0);
}

void ptreeupdate_helper(struct ptree *p,double t,double DT,double weight,int verbose)
{
  /* call if t>= p->update-last + p->update_every */
  int nt=0,tab_back=0,exit_flag=0,tab=0;
  struct llist *LL=NULL,*L2=NULL;
  struct litem *l=NULL,*l2=NULL;
  struct llist *l0=NULL;
  struct region *r=NULL;
  if (verbose){ printf(" %% [calling ptreeupdate_helper] \n");}
  assert(t >= p->update_last + p->update_every);
  if (llitemlength(p->eventra[p->tab])>0){
    if (verbose){ printf(" %% have %d events, continuing\n",llitemlength(p->eventra[p->tab]));}
    if (verbose>1){ printf(" %% creating LL\n");}
    LL = llistmake();
    if (verbose>1){ printf(" %% creating l0\n");}
    l0 = llistmake();
    if (verbose>1){ printf(" %% growing l0 for tab %d\n",p->tab);}
    llistgrowllitem(l0,p->eventra[p->tab]);
    if (verbose>1){ printf(" %% pruning and sorting l0 \n");}
    llistprune(l0); llistsort(l0->first,l0->last,l0->length,&region2region_compare_last_event);
    if (verbose>1){ printf(" %% adding l0 to LL\n");}
    litemadd(LL,l0);
    if (verbose>1){ printf(" %% now adding further to LL\n");}
    tab_back=1;exit_flag=(LL->length==p->nlegs+1);
    while (!exit_flag){
      l0 = llistmake();
      for (nt=0;nt<p->legtime;nt++){
	tab = periodize(p->tab-tab_back-nt,0,p->length);
	if (verbose>5){ printf(" %% growing l0 for tab %d\n",tab);}
	llistgrowllitem(l0,p->eventra[tab]);}
      if (l0->length>0){
	if (verbose>5){ printf(" %% pruning and sorting l0\n");}
	llistprune(l0); llistsort(l0->first,l0->last,l0->length,&region2region_compare_last_event);
	if (verbose>5){ printf(" now adding further to LL\n");}
	litemadd(LL,l0);
	if (verbose>1){ printf(" %% old tab_back=%d . . .",tab_back);}
	tab_back += p->legtime;
	if (verbose>1){ printf(" new tab_back=%d\n",tab_back);}
	exit_flag = LL->length==p->nlegs+1;
	if (verbose>1){ printf(" exit_flag set to %d\n",exit_flag);}}
      else /* if (l0->length==0) */{
	if (verbose>5){ printf(" %% no events between tabs %d and %d\n",p->tab-tab_back,p->tab-tab_back-nt);}
	llisttfree(l0); l0=NULL;
	exit_flag = 1;}}
    if (verbose>5){
      printf(" %% at this point LL is composed of %d region llists:\n",LL->length);
      l = LL->first;
      while (l!=NULL){
	L2 = (struct llist *)l->item;
	printf(" %% %% llist %d of length %d\n",(int)L2,L2->length);
	l2 = L2->first;
	while (l2!=NULL){
	  r = (struct region *) l2->item;
	  printf(" %% %% %% region %d with label %d and last_event %0.3f\n",(int)r,r->label,r->last_event);
	  l2=l2->child;}
	l=l->child;}}
    if (verbose>6){ printf(" %% about to start pstrengthen_starter, ptree is:\n"); pnodeprintf(NULL,p->postree,0,0);}
    pstrengthen_starter(p,LL,weight);
    if (verbose>6){ printf(" %% just finished pstrengthen_starter, ptree is:\n"); pnodeprintf(NULL,p->postree,0,0);}}
  else{
    if (verbose){ printf(" %% have %d events, not updating\n",llitemlength(p->eventra[p->tab]));}}
  p->tab += 1; p->tab = periodize(p->tab,0,p->length);
  llitemtfree(p->eventra[p->tab],NULL); p->eventra[p->tab] = llitemmake();
  p->update_last = maximum(t,p->update_last + p->update_every);
}

void ptreeupdate(struct ptree *p,double t,double DT,double weight)
{
  /* updates ptree *p */
  int verbose=0;
  int nr=0;
  if (verbose){ printf(" %% [starting ptreeupdate] with p=%d, t=%0.3f, DT=%0.3f, weight=%0.3f\n",(int)p,t,DT,weight);}
  for (nr=0;nr<p->nregions;nr++){
    if (p->gate_flag && region_has_event(p->regionra[nr],t,DT) /* determine activity and set last_event_flag */){ 
      if (llitemaddorfind(0,p->eventra[p->tab],p->regionra[nr],&region2region_compare_label)!=NULL){ 
	if (verbose>1){ printf(" %% region %d active, but already in llitem\n",nr);}}
      else{
	if (verbose>1){ printf(" %% region %d active, adding to llitem\n",nr);}
	llitemaddorfind(1,p->eventra[p->tab],p->regionra[nr],&region2region_compare_label);}}
    else if (verbose>1){ printf(" %% region %d inactive\n",nr);}}
  p->total_time += DT;
  if (t >= p->update_last + p->update_every){ ptreeupdate_helper(p,t,DT,weight,verbose);}
}

void pstrengthen_starter(struct ptree *p,struct llist *LL,double weight)
{
  /* starting recursive dropdown, and might as well free LL while we are at it */
  int verbose=0;
  struct litem *l=NULL,*l2=NULL;
  struct region *r=NULL;
  struct llist *L=NULL;
  if (verbose){ 
    printf(" %% [starting pstrengthen_starter] p=%d, LL->length=%d\n",(int)p,LL->length);
    printf(" %% LL is composed of llists of regions:\n");
    l = LL->first;
    while (l!=NULL){
      L = (struct llist *) l->item;
      printf(" %% %% length %d llist\n",L->length);
      l2 = L->first;
      while (l2!=NULL){
	r = (struct region *) l2->item;
	printf(" %% %% %% region %d\n",r->label);
	l2=l2->child;}
      l=l->child;}}
  while (LL->length>=1){
    pstrengthen_helper(1,p,NULL,p->postree,LL,LL->length,weight);
    pstrengthen_helper(0,p,NULL,p->pretree,LL,LL->length,weight);
    L = (struct llist *) LL->last->item;
    if (verbose){ printf(" %% trying to free final llist L=%d\n",(int)L);}
    litemminus(LL,L);
    llisttfree(L);L=NULL;
    if (verbose){ printf(" %% now LL->length=%d\n",LL->length);}}
  if (verbose){ printf(" %% now we have exhausted all the llists in LL (now of length %d)\n",LL->length);}
  llisttfree(LL); LL=NULL;
}

void pstrengthen_helper(int postorpre,struct ptree *p,struct pnode *parent,struct llitem *childllitem,struct llist *LL,int length,double weight)
{
  /* given a llitem of pnode2pnode_compare_label sorted pnodes, and a length-long llist *LL of region2region_compare_last_event
     sorted llists of regional events */
  int verbose=0;
  struct llist *L=NULL;
  struct litem *l=NULL,*l2=NULL;
  struct region *r=NULL,*r2=NULL;
  struct llitem *l0=NULL,*l02=NULL;
  struct pnode *pn=NULL,*pn2=NULL;
  int depth=0;
  if (verbose){ printf(" %% [starting pstrengthen_helper] with postorpre=%d, parent %d, childllitem %d of length %d, LL %d of length %d, length %d, weight %0.3f\n",postorpre,(int)parent,(int)childllitem,llitemlength(childllitem),(int)LL,LL->length,length,weight);}
  assert(LL->length>=length);
  if (parent!=NULL){ assert(parent->childllitem==childllitem);}
  if (postorpre==1){ /* postree - read length from front */
    l=LL->first; depth=1; while(l!=NULL && depth<length){ depth+=1; l=l->child;} L = (struct llist *) l->item;}
  else /* if (postorpre==0) */ { /* pretree - read length from back */
    l=LL->last; depth=1; while(l!=NULL && depth<length){ depth+=1; l=l->parent;} L = (struct llist *) l->item;}
  if (verbose){ printf(" %% read out llist L=%d of length %d\n",(int)L,L->length);}
  if (length==1){
    if (parent==NULL && LL->length==1){ /* singleton llist, must loop through causal events */
      if (postorpre==0){ /* pretree, starting at the end and heading back */
	if (verbose){ printf(" %% starting at the end and heading back\n");}
	l=L->last;
	while (l!=NULL){
	  r = (struct region *) l->item;
	  if (verbose>5){ printf(" %% at region %d\n",r->label);}
	  if ((l0=llitemaddorfind(0,childllitem,r,region2pnode_compare_label))!=NULL){
	    pn = (struct pnode *) l0->item;
	    pn->weight += weight;
	    if (verbose){ printf(" %% now at region %d out of %d\n",r->label,L->length);}
	    if (p->nlegs>0){
	      l2=l->parent;
	      while (l2!=NULL){
		r2 = (struct region *) l2->item;
		if ((l02=llitemaddorfind(0,pn->childllitem,r2,region2pnode_compare_label))!=NULL){
		  pn2 = (struct pnode *) l02->item;
		  pn2->weight += weight;}
		else{
		  pn2 = pnodemake(pn,r2,weight,0);
		  llitemaddorfind(1,pn->childllitem,pn2,pnode2pnode_compare_label);}
		l2 = l2->parent;}}}
	  else{ printf(" %% something's fishy in pstrengthen_helper\n");}
	  l = l->parent;}}
      else /* if (postorpre==1) */{ /* postree, starting at the beginning and heading forward */
	if (verbose){ printf(" %% starting at the beginning and heading forward\n");}
	l=L->first;
	while (l!=NULL){
	  r = (struct region *) l->item;
	  if (verbose>5){ printf(" %% at region %d\n",r->label);}
	  if ((l0=llitemaddorfind(0,childllitem,r,region2pnode_compare_label))!=NULL){
	    pn = (struct pnode *) l0->item;
	    pn->weight += weight;
	    if (verbose){ printf(" %% now at region %d out of %d\n",r->label,L->length);}
	    if (p->nlegs>0){
	      l2=l->child;
	      while (l2!=NULL){
		r2 = (struct region *) l2->item;
		if ((l02=llitemaddorfind(0,pn->childllitem,r2,region2pnode_compare_label))!=NULL){
		  pn2 = (struct pnode *) l02->item;
		  pn2->weight += weight;}
		else{
		  pn2 = pnodemake(pn,r2,weight,0);
		  llitemaddorfind(1,pn->childllitem,pn2,pnode2pnode_compare_label);}
		l2 = l2->child;}}}
	  else{ printf(" %% something's fishy in pstrengthen_helper\n");}
	  l = l->child;}}}
    else /* if (parent!=NULL || LL->length>1) */{ /* end of long llist, simply update */
      l=L->first;
      while (l!=NULL){
	r = (struct region *) l->item;
	if ((l0=llitemaddorfind(0,childllitem,r,region2pnode_compare_label))!=NULL){ /* already exists */
	  pn = (struct pnode *) l0->item;
	  pn->weight += weight;}
	else{ /* must make */
	  pn = pnodemake(parent,r,weight,0);
	  llitemaddorfind(1,childllitem,pn,pnode2pnode_compare_label);}
	l=l->child;}
      /* this is the most suspicious step */
      if (parent!=NULL && llitemlength(parent->childllitem)>2*parent->broodsize){
	llitembalance(parent->childllitem); parent->childllitem = llitemclimb(parent->childllitem); parent->broodsize *= 2;}}}
  else if (length>1){
    l=L->first;
    while (l!=NULL){
      r = (struct region *) l->item;
      if ((l0=llitemaddorfind(0,childllitem,r,region2pnode_compare_label))!=NULL){ /* descend */
	pn = (struct pnode *) l0->item;
	pstrengthen_helper(postorpre,p,pn,pn->childllitem,LL,length-1,weight);}
      else{ /* do nothing */ }
      l=l->child;}}
}

void pnodeprintf(struct pnode *parent,struct llitem *l0,int maxlevel,int level)
{
  /* recursively prints out the pnodes in llitem *l0 */
  char text[128];
  int nl=0;
  struct llist *L=llistmake();
  struct litem *l=NULL;
  struct pnode *pn2=NULL;
  if (l0!=NULL){ llistgrowllitem(L,l0);}
  if (level==0){ sprintf(text,"|__");}
  else{ sprintf(text,"|  "); for (nl=0;nl<level-1;nl++){ sprintf(text,"%s   ",text);} sprintf(text,"%s|__",text);}
  if (parent!=NULL){ assert(parent->childllitem==llitemclimb(l0));}
  //if (L->length>0){ printf(" %%%s %d->childllist of length %d contains:\n",text,(int)parent,L->length);}
  l = L->first;
  while (l!=NULL){
    pn2 = (struct pnode *) l->item;
    printf(" %%%s pnode %d->parent %d, region %d->label %d, weight %0.1f relevance %0.1f",text,(int)pn2,(int)pn2->parent,(int)pn2->region,(int)pn2->region->label,pn2->weight,pn2->relevance);
    if (llitemlength(pn2->childllitem)>0 && (maxlevel==-1 || level<maxlevel)){
      printf(", now descending into %d children...\n",llitemlength(pn2->childllitem));
      pnodeprintf(pn2,pn2->childllitem,maxlevel,level+1);}
    else{ printf("\n");}
    l=l->child;}
  llisttfree(L);L=NULL;
}

void pnodeprune_starter(int posorpre,struct ptree *p,struct pnode *parent,struct llitem *childllitem,struct llist *L)
{
  /* starts pruning process by descending p[re,os]tree and passing leaves (allong with [de,as]cendent llist *L) to pnodeprune_helper */
  int firstentry=0;
  struct pnode *pn=NULL;
  if (parent!=NULL && childllitem->parent==NULL){ /* first descent into parent->childllitem */ 
    assert(parent->childllitem==childllitem); firstentry=1;}
  if (firstentry){ litemadd(L,parent->region);}
  if (firstentry){ pnodeprune_helper(posorpre,p,parent,L);}
  if (childllitem->kidl!=NULL){ pnodeprune_starter(posorpre,p,parent,childllitem->kidl,L);}
  if (childllitem->item!=NULL){ pn=(struct pnode *)childllitem->item; pnodeprune_starter(posorpre,p,pn,pn->childllitem,L);}
  if (childllitem->kidr!=NULL){ pnodeprune_starter(posorpre,p,parent,childllitem->kidr,L);}
  if (firstentry){ llistkilllast(L);}
}

void pnodeprune_helper(int posorpre,struct ptree *p,struct pnode *pn,struct llist *L)
{
  /* given a llist of regions *L and a leaf node *pn (of depth L->length) */
  int verbose=0;
  struct llitem *childllitem1=NULL,*childllitem2=NULL,*childllitem3=NULL;
  struct litem *l=NULL;
  struct llitem *l0=NULL;
  struct region *r=NULL;
  struct pnode *pn2=NULL;
  int depth=0;
  double weight_expected=0,weight_first=0,weight_second=0,weight_third=0;
  double *ra1=NULL,*ra2=NULL,*ra3=NULL;
  double temp=0,relevance=0;
  if (verbose){
    if (posorpre){ childllitem1=p->postree;} else{ childllitem1=p->pretree;}
    l=L->first; depth=0; 
    while (l!=NULL){
      r = (struct region *) l->item;
      if ((l0=llitemaddorfind(0,childllitem1,r,region2pnode_compare_label))!=NULL){
	pn2 = (struct pnode *) l0->item;
	childllitem1 = pn2->childllitem;
	l=l->child;
	depth++;}
      else{ l=NULL;}}
    if (depth==L->length && pn2==pn){ printf(" %% passed correct node && lineage to pnodeprune_helper\n");}}
  if (L->length==1){
    weight_expected = pn->weight;
    pn->relevance = pnode_shear(pn->weight,weight_expected);}
  else if (L->length==2){ /* single leg */
    if (posorpre){ childllitem1=p->postree;} else{ childllitem1=p->pretree;}
    if (verbose){ printf(" %% %% only one leg\n");}
    r = (struct region *) L->first->item; weight_first=0;
    if ((l0=llitemaddorfind(0,childllitem1,r,region2pnode_compare_label))!=NULL){
      pn2 = (struct pnode *) l0->item;
      weight_first = pn2->weight;}
    r = (struct region *) L->last->item; weight_second=0;
    if ((l0=llitemaddorfind(0,childllitem1,r,region2pnode_compare_label))!=NULL){
      pn2 = (struct pnode *) l0->item;
      weight_second = pn2->weight;}
    weight_expected = weight_first*p->legtime*weight_second/p->total_time;
    pn->relevance = pnode_shear(pn->weight,weight_expected);}
  else if (L->length>2){ /* multiple legs */
    if (posorpre){ childllitem1=p->postree; childllitem2=p->pretree;} 
    else{ childllitem1=p->pretree; childllitem2=p->postree;}
    childllitem3=childllitem1;
    ra1 = (double *) tcalloc(L->length,sizeof(double));
    ra2 = (double *) tcalloc(L->length,sizeof(double));
    ra3 = (double *) tcalloc(L->length,sizeof(double));
    l=L->first; depth=0; 
    while (l!=NULL){
      r = (struct region *) l->item;
      if ((l0=llitemaddorfind(0,childllitem1,r,region2pnode_compare_label))!=NULL){
	pn2 = (struct pnode *) l0->item;
	ra1[depth] = pn2->weight;
	childllitem1 = pn2->childllitem;}
      else{ childllitem1=NULL;}
      l=l->child; depth++;}
    l=L->last; depth=0;
    while (l!=NULL){
      r = (struct region *) l->item;
      if ((l0=llitemaddorfind(0,childllitem2,r,region2pnode_compare_label))!=NULL){
	pn2 = (struct pnode *) l0->item;
	ra2[depth] = pn2->weight;
	childllitem2 = pn2->childllitem;}
      else{ childllitem2=NULL;}
      l=l->parent; depth++;}
    l=L->first; depth=0;
    while (l!=NULL){
      r = (struct region *) l->item;
      if ((l0=llitemaddorfind(0,childllitem3,r,region2pnode_compare_label))!=NULL){
	pn2 = (struct pnode *) l0->item;
	ra3[depth] = pn2->weight;}
      l=l->child;
      depth++;}
    relevance=0;
    for (depth=1;depth<L->length-1;depth++){
      weight_first = ra1[depth];
      weight_second = ra2[L->length-1-depth];
      weight_third = ra3[depth];
      weight_expected = weight_first*weight_second/weight_third;
      temp = pnode_shear(pn->weight,weight_expected);
      if (fabs(temp)>fabs(relevance)){ relevance=temp;}}
    pn->relevance = relevance;
    tfree(ra1);tfree(ra2);tfree(ra3);}
}

double pnode_shear(double pn_weight,double weight_expected)
{
  /* some sort of shear for pruning trees */
  double weight_tol=10;
  double relevance=0;
  if (pn_weight<weight_tol && weight_expected<weight_tol){ relevance=0;}
  else{ relevance = log(pn_weight/weight_expected);}
  if (!finite(relevance)){ relevance=0;}
  //if (!finite(relevance)){ if (pn_weight==0){ relevance=-10;} else if (weight_expected==0){ relevance=10;} else{ relevance=0;}}
  return relevance;
}

void pnodeZZZ_starter(int posorpre,struct ptree *p,struct pnode *parent,struct llitem *childllitem,double threshold,int depth)
{
  /* starts ZZZing process by descending p[re,os]tree and passing leaves (allong with [de,as]cendent llist *L) to pnodeZZZ_helper */
  int firstentry=0;
  struct pnode *pn=NULL;
  if (parent!=NULL && childllitem->parent==NULL){ /* first descent into parent->childllitem */ 
    assert(parent->childllitem==childllitem); firstentry=1;}
  if (firstentry && depth<p->nlegs){ pnodeZZZ_helper(posorpre,p,parent,threshold);}
  if (childllitem->kidl!=NULL){ pnodeZZZ_starter(posorpre,p,parent,childllitem->kidl,threshold,depth+1);}
  if (childllitem->item!=NULL){ pn=(struct pnode *)childllitem->item; pnodeZZZ_starter(posorpre,p,pn,pn->childllitem,threshold,depth+1);}
  if (childllitem->kidr!=NULL){ pnodeZZZ_starter(posorpre,p,parent,childllitem->kidr,threshold,depth+1);}
}

void pnodeZZZ_helper(int posorpre,struct ptree *p,struct pnode *pn,double threshold)
{
  /* checks pnode *pn for missing nodes given reference ptree *p */
  struct litem *l=NULL;
  struct llitem *l0=NULL;
  struct pnode *pn2=NULL;
  double expectance=0;
  struct llist *L=llistmake();
  if (posorpre){ llistgrowllitem(L,p->postree);} else{ llistgrowllitem(L,p->pretree);}
  l=L->first;
  while (l!=NULL){
    pn2 = (struct pnode *) l->item;
    if ((l0=llitemaddorfind(0,pn->childllitem,pn2->region,&region2pnode_compare_label))==NULL){
      expectance = pn->weight*p->legtime*pn2->weight/p->total_time;
      if (abs(expectance)>threshold){
	llitemaddorfind(1,pn->childllitem,pnodemake(pn,pn2->region,0,pnode_shear(0,expectance)),&pnode2pnode_compare_label);}}
    l=l->child;}
  llisttfree(L);L=NULL;
}

void pnodebalance_starter(struct pnode *parent,struct llitem *childllitem)
{
  int firstentry=0;
  struct pnode *pn=NULL;
  if (parent!=NULL && childllitem->parent==NULL){ /* first descent into parent->childllitem */ 
    assert(parent->childllitem==childllitem); firstentry=1;}
  if (childllitem->kidl!=NULL){ pnodebalance_starter(parent,childllitem->kidl);}
  if (childllitem->kidr!=NULL){ pnodebalance_starter(parent,childllitem->kidr);}
  if (childllitem->item!=NULL){ pn=(struct pnode *)childllitem->item; pnodebalance_starter(pn,pn->childllitem);}
  if (firstentry){ 
    if (parent!=NULL && llitemlength(parent->childllitem)>2*parent->broodsize){ 
      llitembalance(parent->childllitem); parent->childllitem = llitemclimb(parent->childllitem); 
      parent->broodsize = maximum(parent->broodsize*2,llitemlength(parent->childllitem));}}
}

void ptreerate(struct ptree *p)
{
  int verbose=0;
  struct llist *L=NULL;
  struct pnode *pn=NULL;
  if (verbose){ printf(" %% starting ptreerate, calling pnodeprune_starter\n");}
  L=llistmake(); pnodeprune_starter(1,p,NULL,p->postree,L); llisttfree(L); L=NULL;
  L=llistmake(); pnodeprune_starter(0,p,NULL,p->pretree,L); llisttfree(L); L=NULL;
  if (GLOBAL_PTREE_ZZZ){
    if (verbose){ printf(" %% finished calling pnodeprune_starter, calling pnodeZZZ_starter\n");}
    pnodeZZZ_starter(1,p,NULL,p->postree,1,0);
    pnodeZZZ_starter(0,p,NULL,p->pretree,1,0);
    if (verbose){ printf(" %% finished calling pnodeZZZ_starter, calling pnodebalance_starter\n");}}
  pnodebalance_starter(NULL,p->postree);
  pnodebalance_starter(NULL,p->pretree);
  L=llistmake(); pnode2llist_starter(NULL,p->postree,1,0,L); llistsort(L->first,L->last,L->length,pnode2pnode_compare_weight);
  pn = (struct pnode *) L->last->item;
  histtfree(p->wh); p->wh=histmake(128,pn->weight,0);
  llisttfree(L); L=NULL;
  histtfree(p->rh); p->rh=histmake(128,+4,-4);
  pnode2hist_starter(NULL,p->postree,-1,0,p->wh,p->rh);
}

struct llist * ptreextract_starter(struct ptree *p,double minworth,int nregions)
{
  /* extracts nregions most relevant regions from rated ptree *p, 
     simply by scanning through all legs and determining which 
     regions have the greatest worth = (relevance>threshold)*weight
     last_event is used to hold worth */
  int verbose=0;
  struct region **regionra=NULL;
  int nr=0,depth=0;
  struct llist *Lr=NULL,*Ll=NULL;
  struct litem *l=NULL;
  struct region *r=NULL;
  double *label=NULL;
  if (verbose){ printf(" %% [entering ptreextract_starter] with minworth %0.2f and nregions %d\n",minworth,nregions);}
  if (verbose){ printf(" %% making temporary regionra\n");}
  regionra = (struct region **) tcalloc(p->nregions,sizeof(struct region *));
  for (nr=0;nr<p->nregions;nr++){
    regionra[nr] = (struct region *) tcalloc(1,sizeof(struct region));
    regionra[nr]->label = p->regionra[nr]->label;
    regionra[nr]->neuronllist = NULL;}
  if (verbose){ printf(" %% calling ptreextract_helper\n");}
  ptreextract_helper(NULL,p->postree,regionra,minworth);
  if (verbose){ printf(" %% making and sorting regionllist\n");}
  Lr = llistmake(); for (nr=0;nr<p->nregions;nr++){ litemadd(Lr,regionra[nr]);} tfree(regionra); regionra=NULL;
  llistsort(Lr->first,Lr->last,Lr->length,&region2region_compare_last_event);
  if (verbose){ printf(" %% making labellist\n");}
  Ll = llistmake(); depth=0;
  l=Lr->last; 
  while (l!=NULL && depth<nregions){
    r = (struct region *) l->item;
    if (verbose){ printf(" %% %% adding label %d\n",r->label);}
    label = (double *) tcalloc(1,sizeof(double)); *label = (double)r->label; litemadd(Ll,label);
    l=l->parent; depth+=1;}
  if (verbose){ printf(" %% freeing temporary regions\n");}
  l=Lr->first;
  while (l!=NULL){
    r = (struct region *) l->item;
    tfree(r); r=NULL; l->item=NULL;
    l=l->child;}
  llisttfree(Lr);
  llistsort(Ll->first,Ll->last,Ll->length,&double_compare);
  return Ll;
}

void ptreextract_helper(struct pnode *parent,struct llitem *l0,struct region **regionra,double minworth)
{
  int firstentry=0;
  struct pnode *pn=NULL;
  if (parent!=NULL && l0->parent==NULL){ /* first descent into parent->childllitem */ 
    assert(parent->childllitem==l0); firstentry=1;}
  if (firstentry){ regionra[parent->region->label]->last_event += pnode_worth(parent->weight,parent->relevance,minworth);}
  if (l0->kidl!=NULL){ ptreextract_helper(parent,l0->kidl,regionra,minworth);}
  if (l0->item!=NULL){ pn=(struct pnode *)l0->item; ptreextract_helper(pn,pn->childllitem,regionra,minworth);}
  if (l0->kidr!=NULL){ ptreextract_helper(parent,l0->kidr,regionra,minworth);}
}

double pnode_worth(double weight,double relevance,double minworth){ return weight*(fabs(relevance)+minworth);}

void ptree2jpg_starter(struct ptree *p,char *filename_given,double wmax,double wmin,double rmax,double rmin)
{
  /* assuming filename_given starts with "./" */
  int remove_flag=0,jpg_flag=0,eps_flag=(ON_MY_COMPUTER ? 1 : 0);
  int coloringtype=7;
  char filename_w[256],filename2_w[512],filename_r[256],filename2_r[512],filename_wh[512],filename_rh[512];
  char filename_ws[256],filename2_ws[512],filename_rs[256],filename2_rs[512];
  FILE *fpw=NULL,*fpr=NULL,*fpws=NULL,*fprs=NULL;
  struct llist *L=NULL;
  char command[1024];
  double wmax2=0,wmin2=0,wmean=0,wstd=0,rmax2=0,rmin2=0,rmean=0,rstd=0;
  pnodestats_starter(NULL,p->postree,3,-1,0,NULL,&wmax2,&wmin2,&wmean,&wstd,&rmax2,&rmin2,&rmean,&rstd);
  if (wmax<wmin){ wmax = wmean + STD_VIEW*wstd; wmin = wmean - STD_VIEW*wstd;}
  else if (wmax==wmin){ wmax = wmax2; wmin = wmin2;}
  if (rmax<rmin){ rmax = rmean + STD_VIEW*rstd; rmin = rmean - STD_VIEW*rstd;}
  else if (rmax==rmin){ rmax = rmax2; rmin = rmin2;}
  if (filename_given==NULL){ sprintf(filename_w,"./ptree_%srecord_w",GLOBAL_STRING_2);}
  else{ sprintf(filename_w,"%s_w",filename_given);}
  sprintf(filename2_w,"%s.fig",filename_w);
  if ((fpw=fopen(filename2_w,"w"))==NULL){ printf("warning, cannot create %s in ptree2jpg_starter",filename2_w); fpw=stdout;}
  if (filename_given==NULL){ sprintf(filename_r,"./ptree_%srecord_r",GLOBAL_STRING_2);}
  else{ sprintf(filename_r,"%s_r",filename_given);}
  sprintf(filename2_r,"%s.fig",filename_r);
  if ((fpr=fopen(filename2_r,"w"))==NULL){ printf("warning, cannot create %s in ptree2jpg_starter",filename2_r); fpr=stdout;}
  if (filename_given==NULL){ sprintf(filename_ws,"./ptree_%srecord_ws",GLOBAL_STRING_2);}
  else{ sprintf(filename_ws,"%s_ws",filename_given);}
  sprintf(filename2_ws,"%s.fig",filename_ws);
  if ((fpws=fopen(filename2_ws,"w"))==NULL){ printf("warning, cannot create %s in ptree2jpg_starter",filename2_ws); fpws=stdout;}
  if (filename_given==NULL){ sprintf(filename_rs,"./ptree_%srecord_rs",GLOBAL_STRING_2);}
  else{ sprintf(filename_rs,"%s_rs",filename_given);}
  sprintf(filename2_rs,"%s.fig",filename_rs);
  if ((fprs=fopen(filename2_rs,"w"))==NULL){ printf("warning, cannot create %s in ptree2jpg_starter",filename2_rs); fprs=stdout;}
  if (filename_given==NULL){ sprintf(filename_wh,"./ptree_%srecord_wh",GLOBAL_STRING_2);}
  else{ sprintf(filename_wh,"%s_wh",filename_given);}
  if (filename_given==NULL){ sprintf(filename_rh,"./ptree_%srecord_rh",GLOBAL_STRING_2);}
  else{ sprintf(filename_rh,"%s_rh",filename_given);}
  histdump(p->wh,0,filename_wh," ",0); histdump(p->wh,1,filename_wh,"wh",0);
  histdump(p->rh,0,filename_rh," ",0); histdump(p->rh,1,filename_rh,"rh",0);
  fprintf(fpw,"%s",FIG_PREAMBLE); fprintf(fpr,"%s",FIG_PREAMBLE);
  fprintf(fpws,"%s",FIG_PREAMBLE); fprintf(fprs,"%s",FIG_PREAMBLE);
  switch (coloringtype){
  case 5: fprintf(fpw,"%s",FIG_PREAMBLE_COLOR_5); fprintf(fpr,"%s",FIG_PREAMBLE_COLOR_5); break;
  case 7: fprintf(fpw,"%s",FIG_PREAMBLE_COLOR_7); fprintf(fpr,"%s",FIG_PREAMBLE_COLOR_7); break;
  default: break;}
  fprintf(fpws,"%s",FIG_PREAMBLE_COLOR_7); fprintf(fprs,"%s",FIG_PREAMBLE_COLOR_7);
  ptree2jpg_helper(+1,p->nregions,wmax,wmin,rmax,rmin,NULL,NULL,fpw,fpr,NULL);
  ptree2jpg_helper(-1,p->nregions,wmax,wmin,rmax,rmin,NULL,NULL,fpw,fpr,NULL);
  L=llistmake();
  ptree2jpg_helper(+1,p->nregions,wmax,wmin,rmax,rmin,NULL,p->postree,fpw,fpr,L);
  llisttfree(L); L=NULL;
  L=llistmake();
  ptree2jpg_helper(-1,p->nregions,wmax,wmin,rmax,rmin,NULL,p->pretree,fpw,fpr,L);
  llisttfree(L); L=NULL;
  L=llistmake();
  ptree2star_helper(+1,p->nlegs,p->nregions,wmax2,1,rmax,rmin,NULL,p->postree,fpws,fprs,L);
  ptree2star_helper(-1,p->nlegs,p->nregions,wmax2,1,rmax,rmin,NULL,p->pretree,fpws,fprs,L);
  llisttfree(L); L=NULL;
  if (fpw!=stdout){ fclose(fpw);} if (fpr!=stdout){ fclose(fpr);}
  if (fpws!=stdout){ fclose(fpws);} if (fprs!=stdout){ fclose(fprs);}
  if (jpg_flag){ sprintf(command,"fig2dev -Ljpeg -q %d %s.fig %s.jpg;",/*quality*/5,filename_w,filename_w); system(command);}
  if (jpg_flag){ sprintf(command,"fig2dev -Ljpeg -q %d %s.fig %s.jpg;",/*quality*/5,filename_r,filename_r); system(command);}
  if (jpg_flag){ sprintf(command,"fig2dev -Ljpeg -q %d %s.fig %s.jpg;",/*quality*/5,filename_ws,filename_ws); system(command);}
  if (jpg_flag){ sprintf(command,"fig2dev -Ljpeg -q %d %s.fig %s.jpg;",/*quality*/5,filename_rs,filename_rs); system(command);}
  if (eps_flag){ sprintf(command,"fig2dev -Leps %s.fig %s.eps;",filename_w,filename_w); system(command);}
  if (eps_flag){ sprintf(command,"fig2dev -Leps %s.fig %s.eps;",filename_r,filename_r); system(command);}
  if (eps_flag){ sprintf(command,"fig2dev -Leps %s.fig %s.eps;",filename_ws,filename_ws); system(command);}
  if (eps_flag){ sprintf(command,"fig2dev -Leps %s.fig %s.eps;",filename_rs,filename_rs); system(command);}
  if (remove_flag){ sprintf(command,"rm %s.fig;",filename_w); system(command);}
  if (remove_flag){ sprintf(command,"rm %s.fig;",filename_r); system(command);}
  if (remove_flag){ sprintf(command,"rm %s.fig;",filename_ws); system(command);}
  if (remove_flag){ sprintf(command,"rm %s.fig;",filename_rs); system(command);}
}

void ptree2jpg_helper(int posorpre,int nregions,double wmax,double wmin,double rmax,double rmin,struct pnode *parent,struct llitem *l0,FILE *fpw,FILE *fpr,struct llist *L)
{
  /* runs through llitem *l0 of pnode items, the relevance is dumped up to depth 3
     any pnode order works, but we do reachfirst.
  */
  int firstentry=0;
  int nr=0,maxdia=10000;
  struct pnode *pn=NULL;
  double rcolor=0,gcolor=0,bcolor=0;
  int colorcode=0;
  struct region *r1=NULL,*r2=NULL,*r3=NULL,*r4=NULL;
  double rord = 1.0/(double)nregions,side=0,xside=0,yside=0;
  double xord=0,yord=0;
  if (l0==NULL){ /* draw skeleton */
    xord = 0.5 + 0.1*posorpre;
    yord = 1.25*(double)nregions/2.0 + 1*(double)(nregions-1.0)/2.0*rord;
    xside = rord; yside = 1.0;
    if (posorpre==+1){
      fprintf(fpw,"4 0 0 %d 0 %d %d 0.0000 4 %d %d %d %d %0.5f\\001\n",/*depth*/1,/*font*/12,/*point*/72,/*textheight*/(int)(24*13.5),/*textwidth*/(int)(24*9*18),/*xpos*/(int)floor(maxdia*(xord-2.5*rord)),/*ypos*/(int)floor(maxdia*(yord-(nregions+1)*rord/2)),wmax);
      fprintf(fpw,"4 0 0 %d 0 %d %d 0.0000 4 %d %d %d %d %0.5f\\001\n",/*depth*/1,/*font*/12,/*point*/72,/*textheight*/(int)(24*13.5),/*textwidth*/(int)(24*9*18),/*xpos*/(int)floor(maxdia*(xord-2.5*rord)),/*ypos*/(int)floor(maxdia*(yord+(nregions+1.5)*rord/2)),wmin);
      fprintf(fpr,"4 0 0 %d 0 %d %d 0.0000 4 %d %d %d %d %0.5f\\001\n",/*depth*/1,/*font*/12,/*point*/72,/*textheight*/(int)(24*13.5),/*textwidth*/(int)(24*9*18),/*xpos*/(int)floor(maxdia*(xord-2.5*rord)),/*ypos*/(int)floor(maxdia*(yord-(nregions+1)*rord/2)),rmax);
      fprintf(fpr,"4 0 0 %d 0 %d %d 0.0000 4 %d %d %d %d %0.5f\\001\n",/*depth*/1,/*font*/12,/*point*/72,/*textheight*/(int)(24*13.5),/*textwidth*/(int)(24*9*18),/*xpos*/(int)floor(maxdia*(xord-2.5*rord)),/*ypos*/(int)floor(maxdia*(yord+(nregions+1.5)*rord/2)),rmin);}
    fprintf(fpw,"2 2 0 4 0 7 900 0 -1 0.000 0 0 -1 0 0 5\n");
    fprintf(fpw,"\t %d %d %d %d %d %d %d %d %d %d\n",(int)floor(maxdia*(xord+xside/2)),(int)floor(maxdia*(yord+yside/2)),(int)floor(maxdia*(xord+xside/2)),(int)floor(maxdia*(yord-yside/2)),(int)floor(maxdia*(xord-xside/2)),(int)floor(maxdia*(yord-yside/2)),(int)floor(maxdia*(xord-xside/2)),(int)floor(maxdia*(yord+yside/2)),(int)floor(maxdia*(xord+xside/2)),(int)floor(maxdia*(yord+yside/2)));
    fprintf(fpr,"2 2 0 4 0 7 900 0 -1 0.000 0 0 -1 0 0 5\n");
    fprintf(fpr,"\t %d %d %d %d %d %d %d %d %d %d\n",(int)floor(maxdia*(xord+xside/2)),(int)floor(maxdia*(yord+yside/2)),(int)floor(maxdia*(xord+xside/2)),(int)floor(maxdia*(yord-yside/2)),(int)floor(maxdia*(xord-xside/2)),(int)floor(maxdia*(yord-yside/2)),(int)floor(maxdia*(xord-xside/2)),(int)floor(maxdia*(yord+yside/2)),(int)floor(maxdia*(xord+xside/2)),(int)floor(maxdia*(yord+yside/2)));
    for (nr=0;nr<nregions;nr++){
      xord = 0.5 + 0.75*posorpre;
      yord = 1.25*(double)nr + 1*(double)(nregions-1.0)/2.0*rord;
      xside = rord; yside = 1.0;
      fprintf(fpw,"2 2 0 4 0 7 900 0 -1 0.000 0 0 -1 0 0 5\n");
      fprintf(fpw,"\t %d %d %d %d %d %d %d %d %d %d\n",(int)floor(maxdia*(xord+xside/2)),(int)floor(maxdia*(yord+yside/2)),(int)floor(maxdia*(xord+xside/2)),(int)floor(maxdia*(yord-yside/2)),(int)floor(maxdia*(xord-xside/2)),(int)floor(maxdia*(yord-yside/2)),(int)floor(maxdia*(xord-xside/2)),(int)floor(maxdia*(yord+yside/2)),(int)floor(maxdia*(xord+xside/2)),(int)floor(maxdia*(yord+yside/2)));
      fprintf(fpr,"2 2 0 4 0 7 900 0 -1 0.000 0 0 -1 0 0 5\n");
      fprintf(fpr,"\t %d %d %d %d %d %d %d %d %d %d\n",(int)floor(maxdia*(xord+xside/2)),(int)floor(maxdia*(yord+yside/2)),(int)floor(maxdia*(xord+xside/2)),(int)floor(maxdia*(yord-yside/2)),(int)floor(maxdia*(xord-xside/2)),(int)floor(maxdia*(yord-yside/2)),(int)floor(maxdia*(xord-xside/2)),(int)floor(maxdia*(yord+yside/2)),(int)floor(maxdia*(xord+xside/2)),(int)floor(maxdia*(yord+yside/2)));}
    for (nr=0;nr<nregions;nr++){
      xord = 1.75*posorpre + (double)(nregions-1.0)/2.0*rord;
      yord = 1.25*(double)nr + 1*(double)(nregions-1.0)/2.0*rord;
      xside = 1.0; yside = 1.0;
      fprintf(fpw,"2 2 0 4 0 7 900 0 -1 0.000 0 0 -1 0 0 5\n");
      fprintf(fpw,"\t %d %d %d %d %d %d %d %d %d %d\n",(int)floor(maxdia*(xord+xside/2)),(int)floor(maxdia*(yord+yside/2)),(int)floor(maxdia*(xord+xside/2)),(int)floor(maxdia*(yord-yside/2)),(int)floor(maxdia*(xord-xside/2)),(int)floor(maxdia*(yord-yside/2)),(int)floor(maxdia*(xord-xside/2)),(int)floor(maxdia*(yord+yside/2)),(int)floor(maxdia*(xord+xside/2)),(int)floor(maxdia*(yord+yside/2)));
      fprintf(fpr,"2 2 0 4 0 7 900 0 -1 0.000 0 0 -1 0 0 5\n");
      fprintf(fpr,"\t %d %d %d %d %d %d %d %d %d %d\n",(int)floor(maxdia*(xord+xside/2)),(int)floor(maxdia*(yord+yside/2)),(int)floor(maxdia*(xord+xside/2)),(int)floor(maxdia*(yord-yside/2)),(int)floor(maxdia*(xord-xside/2)),(int)floor(maxdia*(yord-yside/2)),(int)floor(maxdia*(xord-xside/2)),(int)floor(maxdia*(yord+yside/2)),(int)floor(maxdia*(xord+xside/2)),(int)floor(maxdia*(yord+yside/2)));}
    for (nr=0;nr<nregions;nr++){
      xord = 3.0*posorpre + (double)(nregions-1.0)/2.0*rord;
      yord = 1.25*(double)nr + 1*(double)(nregions-1.0)/2.0*rord;
      xside = 1.0; yside = 1.0;
      fprintf(fpw,"2 2 0 4 0 7 900 0 -1 0.000 0 0 -1 0 0 5\n");
      fprintf(fpw,"\t %d %d %d %d %d %d %d %d %d %d\n",(int)floor(maxdia*(xord+xside/2)),(int)floor(maxdia*(yord+yside/2)),(int)floor(maxdia*(xord+xside/2)),(int)floor(maxdia*(yord-yside/2)),(int)floor(maxdia*(xord-xside/2)),(int)floor(maxdia*(yord-yside/2)),(int)floor(maxdia*(xord-xside/2)),(int)floor(maxdia*(yord+yside/2)),(int)floor(maxdia*(xord+xside/2)),(int)floor(maxdia*(yord+yside/2)));
      fprintf(fpr,"2 2 0 4 0 7 900 0 -1 0.000 0 0 -1 0 0 5\n");
      fprintf(fpr,"\t %d %d %d %d %d %d %d %d %d %d\n",(int)floor(maxdia*(xord+xside/2)),(int)floor(maxdia*(yord+yside/2)),(int)floor(maxdia*(xord+xside/2)),(int)floor(maxdia*(yord-yside/2)),(int)floor(maxdia*(xord-xside/2)),(int)floor(maxdia*(yord-yside/2)),(int)floor(maxdia*(xord-xside/2)),(int)floor(maxdia*(yord+yside/2)),(int)floor(maxdia*(xord+xside/2)),(int)floor(maxdia*(yord+yside/2)));}}
  else /* if (l0!=NULL) */{ /* fill out ptree */
    if (parent!=NULL && l0->parent==NULL){ /* first descent into parent->childllitem */ 
      assert(parent->childllitem==l0); firstentry=1;}
    if (firstentry){ litemadd(L,parent->region);}
    if (firstentry){ 
      switch(L->length){
      case 0: break;
      case 1: 
	r1 = (struct region *) L->last->item;
	colorscale(0,parent->relevance,rmax,rmin,&rcolor,&gcolor,&bcolor);
	colorcode = crop((int)floor(512*rcolor),0,511);
	xord = 0.5 + 0.1*posorpre;
	yord = 1.25*(double)nregions/2.0 + 1*(double)r1->label*rord;
	side = rord;
	fprintf(fpr,"2 2 0 0 -1 %d %d 0 %d 0.000 0 0 -1 0 0 5\n",/*colorcode*/colorcode+32,/*depth*/1+r1->label,/*fill*/20);
	fprintf(fpr,"\t %d %d %d %d %d %d %d %d %d %d\n",(int)floor(maxdia*(xord+side/2)),(int)floor(maxdia*(yord+side/2)),(int)floor(maxdia*(xord+side/2)),(int)floor(maxdia*(yord-side/2)),(int)floor(maxdia*(xord-side/2)),(int)floor(maxdia*(yord-side/2)),(int)floor(maxdia*(xord-side/2)),(int)floor(maxdia*(yord+side/2)),(int)floor(maxdia*(xord+side/2)),(int)floor(maxdia*(yord+side/2)));
	colorscale(0,parent->weight,wmax,wmin,&rcolor,&gcolor,&bcolor);
	colorcode = crop((int)floor(512*rcolor),0,511);
	fprintf(fpw,"2 2 0 0 -1 %d %d 0 %d 0.000 0 0 -1 0 0 5\n",/*colorcode*/colorcode+32,/*depth*/1+r1->label,/*fill*/20);
	fprintf(fpw,"\t %d %d %d %d %d %d %d %d %d %d\n",(int)floor(maxdia*(xord+side/2)),(int)floor(maxdia*(yord+side/2)),(int)floor(maxdia*(xord+side/2)),(int)floor(maxdia*(yord-side/2)),(int)floor(maxdia*(xord-side/2)),(int)floor(maxdia*(yord-side/2)),(int)floor(maxdia*(xord-side/2)),(int)floor(maxdia*(yord+side/2)),(int)floor(maxdia*(xord+side/2)),(int)floor(maxdia*(yord+side/2)));
	break;
      case 2:
	r1 = (struct region *) L->last->parent->item;
	r2 = (struct region *) L->last->item;
	colorscale(0,parent->relevance,rmax,rmin,&rcolor,&gcolor,&bcolor);
	colorcode = crop((int)floor(512*rcolor),0,511);
	xord = 0.5 + 0.75*posorpre;
	yord = 1.25*(double)r1->label + 1*(double)r2->label*rord;
	side = rord;
	fprintf(fpr,"2 2 0 0 -1 %d %d 0 %d 0.000 0 0 -1 0 0 5\n",/*colorcode*/colorcode+32,/*depth*/1+r2->label,/*fill*/20);
	fprintf(fpr,"\t %d %d %d %d %d %d %d %d %d %d\n",(int)floor(maxdia*(xord+side/2)),(int)floor(maxdia*(yord+side/2)),(int)floor(maxdia*(xord+side/2)),(int)floor(maxdia*(yord-side/2)),(int)floor(maxdia*(xord-side/2)),(int)floor(maxdia*(yord-side/2)),(int)floor(maxdia*(xord-side/2)),(int)floor(maxdia*(yord+side/2)),(int)floor(maxdia*(xord+side/2)),(int)floor(maxdia*(yord+side/2)));
	colorscale(0,parent->weight,wmax,wmin,&rcolor,&gcolor,&bcolor);
	colorcode = crop((int)floor(512*rcolor),0,511);
	fprintf(fpw,"2 2 0 0 -1 %d %d 0 %d 0.000 0 0 -1 0 0 5\n",/*colorcode*/colorcode+32,/*depth*/1+r2->label,/*fill*/20);
	fprintf(fpw,"\t %d %d %d %d %d %d %d %d %d %d\n",(int)floor(maxdia*(xord+side/2)),(int)floor(maxdia*(yord+side/2)),(int)floor(maxdia*(xord+side/2)),(int)floor(maxdia*(yord-side/2)),(int)floor(maxdia*(xord-side/2)),(int)floor(maxdia*(yord-side/2)),(int)floor(maxdia*(xord-side/2)),(int)floor(maxdia*(yord+side/2)),(int)floor(maxdia*(xord+side/2)),(int)floor(maxdia*(yord+side/2)));
	break;
      case 3:
	r1 = (struct region *) L->last->parent->parent->item;
	r2 = (struct region *) L->last->parent->item;
	r3 = (struct region *) L->last->item;
	colorscale(0,parent->relevance,rmax,rmin,&rcolor,&gcolor,&bcolor);
	colorcode = crop((int)floor(512*rcolor),0,511);
	xord = 1.75*posorpre + (double)r2->label*rord;
	yord = 1.25*(double)r1->label + 1*(double)r3->label*rord;
	side = rord;
	fprintf(fpr,"2 2 0 0 -1 %d %d 0 %d 0.000 0 0 -1 0 0 5\n",/*colorcode*/colorcode+32,/*depth*/1+r3->label,/*fill*/20);
	fprintf(fpr,"\t %d %d %d %d %d %d %d %d %d %d\n",(int)floor(maxdia*(xord+side/2)),(int)floor(maxdia*(yord+side/2)),(int)floor(maxdia*(xord+side/2)),(int)floor(maxdia*(yord-side/2)),(int)floor(maxdia*(xord-side/2)),(int)floor(maxdia*(yord-side/2)),(int)floor(maxdia*(xord-side/2)),(int)floor(maxdia*(yord+side/2)),(int)floor(maxdia*(xord+side/2)),(int)floor(maxdia*(yord+side/2)));
	colorscale(0,parent->weight,wmax,wmin,&rcolor,&gcolor,&bcolor);
	colorcode = crop((int)floor(512*rcolor),0,511);
	fprintf(fpw,"2 2 0 0 -1 %d %d 0 %d 0.000 0 0 -1 0 0 5\n",/*colorcode*/colorcode+32,/*depth*/1+r3->label,/*fill*/20);
	fprintf(fpw,"\t %d %d %d %d %d %d %d %d %d %d\n",(int)floor(maxdia*(xord+side/2)),(int)floor(maxdia*(yord+side/2)),(int)floor(maxdia*(xord+side/2)),(int)floor(maxdia*(yord-side/2)),(int)floor(maxdia*(xord-side/2)),(int)floor(maxdia*(yord-side/2)),(int)floor(maxdia*(xord-side/2)),(int)floor(maxdia*(yord+side/2)),(int)floor(maxdia*(xord+side/2)),(int)floor(maxdia*(yord+side/2)));
	break;
      case 4:
	r1 = (struct region *) L->last->parent->parent->parent->item;
	r2 = (struct region *) L->last->parent->parent->item;
	r3 = (struct region *) L->last->parent->item;
	r4 = (struct region *) L->last->item;
	colorscale(0,parent->relevance,rmax,rmin,&rcolor,&gcolor,&bcolor);
	colorcode = crop((int)floor(512*rcolor),0,511);
	xord = 3.0*posorpre + (double)r2->label*rord;
	yord = 1.25*(double)r1->label + 1*(double)r3->label*rord;
	side = (double)(1 + r4->label)*rord*rord;
	fprintf(fpr,"2 2 0 0 -1 %d %d 0 %d 0.000 0 0 -1 0 0 5\n",/*colorcode*/colorcode+32,/*depth*/1+r4->label,/*fill*/20);
	fprintf(fpr,"\t %d %d %d %d %d %d %d %d %d %d\n",(int)floor(maxdia*(xord+side/2)),(int)floor(maxdia*(yord+side/2)),(int)floor(maxdia*(xord+side/2)),(int)floor(maxdia*(yord-side/2)),(int)floor(maxdia*(xord-side/2)),(int)floor(maxdia*(yord-side/2)),(int)floor(maxdia*(xord-side/2)),(int)floor(maxdia*(yord+side/2)),(int)floor(maxdia*(xord+side/2)),(int)floor(maxdia*(yord+side/2)));
	colorscale(0,parent->weight,wmax,wmin,&rcolor,&gcolor,&bcolor);
	colorcode = crop((int)floor(512*rcolor),0,511);
	fprintf(fpw,"2 2 0 0 -1 %d %d 0 %d 0.000 0 0 -1 0 0 5\n",/*colorcode*/colorcode+32,/*depth*/1+r4->label,/*fill*/20);
	fprintf(fpw,"\t %d %d %d %d %d %d %d %d %d %d\n",(int)floor(maxdia*(xord+side/2)),(int)floor(maxdia*(yord+side/2)),(int)floor(maxdia*(xord+side/2)),(int)floor(maxdia*(yord-side/2)),(int)floor(maxdia*(xord-side/2)),(int)floor(maxdia*(yord-side/2)),(int)floor(maxdia*(xord-side/2)),(int)floor(maxdia*(yord+side/2)),(int)floor(maxdia*(xord+side/2)),(int)floor(maxdia*(yord+side/2)));
	break;      
      default: /* do nothing */ break;}}
    if (l0->kidl!=NULL){ ptree2jpg_helper(posorpre,nregions,wmax,wmin,rmax,rmin,parent,l0->kidl,fpw,fpr,L);}
    if (l0->item!=NULL){ pn=(struct pnode *)l0->item; ptree2jpg_helper(posorpre,nregions,wmax,wmin,rmax,rmin,pn,pn->childllitem,fpw,fpr,L);}
    if (l0->kidr!=NULL){ ptree2jpg_helper(posorpre,nregions,wmax,wmin,rmax,rmin,parent,l0->kidr,fpw,fpr,L);}
    if (firstentry){ llistkilllast(L);}}
}

void ptree2star_helper(int posorpre,int nlegs,int nregions,double wmax,double wmin,double rmax,double rmin,struct pnode *parent,struct llitem *l0,FILE *fpws,FILE *fprs,struct llist *L)
{
  /* runs through llitem *l0 of pnode items, the relevance is dumped up to depth 3
     any pnode order works, but we do reachfirst.
  */
  int firstentry=0;
  int nr=0,maxdia=10000;
  struct pnode *pn=NULL;
  struct litem *l=NULL;
  double rcolor=0,gcolor=0,bcolor=0;
  int colorcode=0;
  struct region *rg=NULL;
  int nl=0;
  double xord=0,yord=0,xord2=0,yord2=0,rord=0;
  double vx=0,vy=0,vrad=0,epsilon=0;
  if (parent!=NULL && l0->parent==NULL){ /* first descent into parent->childllitem */ 
    assert(parent->childllitem==l0); firstentry=1;}
  if (firstentry){ litemadd(L,parent->region);}
  if (firstentry){ 
    l=L->last;nl=0;xord=0,yord=0,vx=0;vy=0;vrad=PI/(double)nregions;vrad=0.5*sin(vrad)/(1+sin(vrad));if (nregions>3){ epsilon=0.5;}
    while (l!=NULL){
      rg=(struct region *)l->item;
      vx = cos(2*PI*(rg->label+0.5)/(double)nregions)*(0.5-vrad);
      vy = sin(2*PI*(rg->label+0.5)/(double)nregions)*(0.5-vrad);
      xord += pow(3+epsilon,nl)*vx; yord += pow(3+epsilon,nl)*vy;
      nl+=1;
      l=l->parent;}
    xord /= pow(3+epsilon,nlegs); yord /= pow(3+epsilon,nlegs); rord = vrad*pow(3+epsilon,-nlegs);
    xord2 = xord + 0.5*posorpre; yord2 = yord;
    colorscale(0,parent->relevance,rmax,rmin,&rcolor,&gcolor,&bcolor);colorcode=crop((int)floor(512*rcolor),0,511);
    fprintf(fprs,"2 1 0 0 0 %d %d 0 %d 0.000 0 0 -1 0 0 %d\n",/*colorcode*/colorcode+32,/*depth*/1+L->length,/*fill*/20,/*npoints*/nregions+1); fprintf(fprs,"\t"); for (nr=0;nr<=nregions;nr++){ fprintf(fprs,"%d %d ",(int)floor(maxdia*(xord2+rord*cos(2*PI*((double)nr+0.5)/(double)nregions))),(int)floor(maxdia*(yord2+rord*sin(2*PI*((double)nr+0.5)/(double)nregions))));} fprintf(fprs,"\n");
    colorscale(0,log(maximum(1,parent->weight)),log(wmax),log(wmin),&rcolor,&gcolor,&bcolor);colorcode=crop((int)floor(512*rcolor),0,511);
    fprintf(fpws,"2 1 0 0 0 %d %d 0 %d 0.000 0 0 -1 0 0 %d\n",/*colorcode*/colorcode+32,/*depth*/1+L->length,/*fill*/20,/*npoints*/nregions+1); fprintf(fpws,"\t"); for (nr=0;nr<=nregions;nr++){ fprintf(fpws,"%d %d ",(int)floor(maxdia*(xord2+rord*cos(2*PI*((double)nr+0.5)/(double)nregions))),(int)floor(maxdia*(yord2+rord*sin(2*PI*((double)nr+0.5)/(double)nregions))));} fprintf(fpws,"\n");}
  if (l0->kidl!=NULL){ ptree2star_helper(posorpre,nlegs,nregions,wmax,wmin,rmax,rmin,parent,l0->kidl,fpws,fprs,L);}
  if (l0->item!=NULL){ pn=(struct pnode *)l0->item; ptree2star_helper(posorpre,nlegs,nregions,wmax,wmin,rmax,rmin,pn,pn->childllitem,fpws,fprs,L);}
  if (l0->kidr!=NULL){ ptree2star_helper(posorpre,nlegs,nregions,wmax,wmin,rmax,rmin,parent,l0->kidr,fpws,fprs,L);}
  if (firstentry){ llistkilllast(L);}
}

void ptreedump_helper_ascii(struct ptree *p,char *filename,int level)
{
  /* assumes filename starts with "./"
     if level==-1 it plots all levels 
     each level is dumped as 1-dimensional ascii array ordered with standard multi-index convention */
  int verbose=1;
  int nl=0,nr=0,ni=0;
  char filename_extend[1024];
  FILE *fp=NULL;
  struct llist *L=NULL;
  int *ira=NULL,*ira_d=NULL;
  double temp=0;
  if (verbose){ printf(" %% [entering ptreedump_helper_ascii] level %d\n",level);}
  if (level<0 || level>p->nlegs){ for (nl=0;nl<=p->nlegs;nl++){ ptreedump_helper_ascii(p,filename,nl);}}
  else if (level>=0 && level<=p->nlegs){
    nl=level;
    sprintf(filename_extend,"%s_level_%d",filename,nl);
    if ((fp=fopen(filename_extend,"w"))==NULL){ printf("warning, cannot create %s in ptreedump_helper_ascii",filename_extend); fp=stdout;}
    ira = (int *) tcalloc(nl+1,sizeof(int));
    ira_d = (int *) tcalloc(nl+1,sizeof(int)); for (ni=0;ni<nl+1;ni++){ ira_d[ni] = p->nregions;}
    for (nr=0;nr<(int)pow(p->nregions,nl+1);nr++){
      indextract(nr,nl+1,ira_d,ira);
      if (verbose){ printf(" %% region multi-index %d corresponding to: ",nr); raprintf(ira,"int",1,nl+1,"ira: ");}
      L=llistmake(); for (ni=0;ni<nl+1;ni++){ litemadd(L,p->regionra[ira[ni]]);}
      temp=0; pnode_dig_and_evaluate(1,L->first,p,NULL,p->postree,&pnode_evaluate_weight,&temp);
      if (verbose){ printf(" %% found weight %d corresponding to rate %f\n",(int)temp,temp/p->total_time);}
      fprintf(fp,"%0.16f\n",temp/p->total_time);
      llisttfree(L);L=NULL;}
    tfree(ira_d);ira_d=NULL;
    tfree(ira);ira=NULL;
    if (fp!=stdout){ fclose(fp);fp=NULL;}}
}

void ptreedump_starter(struct ptree *p,char *fgvn,int dump_type,int fullname_flag,double wmax,double wmin,double rmax,double rmin,int keepzeroweights)
{
  /* dumps ptree *p as a sparse array, first dumping
     p->nregions
     p->nlegs
     p->legtime
     p->length 
     p->update_every
     p->total_time
     assumes fgvn does NOT start with "./"
     unless fullname_flag==1
  */
  char filename[512];
  FILE *fp=NULL;
  struct llist *L=NULL;
  int bitbybit=0,continue_flag=0;
  if (fullname_flag==1){
    sprintf(filename,"%s",fgvn);}
  else{
    if (GLOBAL_PTREE_BITBYBIT==0){
      if (fgvn==NULL){ sprintf(filename,"./ptree_%srecord_time%d",GLOBAL_STRING_2,(int)floor(GLOBAL_time));}
      else /* if (fgvn!=NULL) */{ sprintf(filename,"./ptree_%s_%srecord_time%d",fgvn,GLOBAL_STRING_2,(int)floor(GLOBAL_time));}}
    else /* if (GLOBAL_PTREE_BITBYBIT==1) */{
      bitbybit=GLOBAL_PTREE_BITBYBIT_CURRENT_NUMBER;continue_flag=0;
      do{
	bitbybit+=1;
	if (fgvn==NULL){ sprintf(filename,"./ptree_%srecord_%d",GLOBAL_STRING_2,bitbybit);}
	else /* if (fgvn!=NULL) */{ sprintf(filename,"./ptree_%s_%srecord_%d",fgvn,GLOBAL_STRING_2,bitbybit);}
	continue_flag=checktofind(filename);}
      while (continue_flag);}}
  if (dump_type==0){
    sprintf(filename,"%s.txt",filename);
    ptreedump_helper_ascii(p,filename,-1);}
  if (dump_type==2){
    if ((fp=fopen(filename,"w"))==NULL){ printf("warning, cannot create %s in ptreedump_starter",filename); fp=stdout;}
    fwrite(&(p->nregions),sizeof(int),1,fp);
    fwrite(&(p->nlegs),sizeof(int),1,fp);
    fwrite(&(p->legtime),sizeof(int),1,fp);
    fwrite(&(p->length),sizeof(int),1,fp);
    fwrite(&(p->update_every),sizeof(double),1,fp);
    fwrite(&(p->total_time),sizeof(double),1,fp);
    L=llistmake();
    ptreedump_helper(1,NULL,p->postree,fp,L,keepzeroweights);
    llisttfree(L); L=NULL;
    L=llistmake();
    ptreedump_helper(0,NULL,p->pretree,fp,L,keepzeroweights);
    llisttfree(L); L=NULL;
    if (fp!=stdout){ fclose(fp);}
    if (fullname_flag!=0 || GLOBAL_PTREE_BITBYBIT==0){ ptree2jpg_starter(p,filename,wmax,wmin,rmax,rmin);}}
}

void ptreedump_helper(int posorpre,struct pnode *parent,struct llitem *l0,FILE *fp,struct llist *L,int keepzeroweights)
{
  /* runs through llitem *l0 of pnode items, the format is:
     parent->weight,parent->relevance,label_0,label_1,...,label_n,posorpre ? -1 : -2
     any pnode order works, but we do reachfirst.
  */
  int verbose=0;
  int firstentry=0;
  int minusone=-1,minustwo=-2;
  struct pnode *pn=NULL;
  struct litem *l=NULL;
  struct region *r=NULL;
  if (parent!=NULL && l0->parent==NULL){ /* first descent into parent->childllitem */ 
    assert(parent->childllitem==l0); firstentry=1;}
  if (firstentry){ litemadd(L,parent->region);}
  if (firstentry){ 
    if (verbose){ printf(" %% firstentry in ptreedump_helper, weight %f relevance %f\n",parent->weight,parent->relevance);}
    fwrite(&(parent->weight),sizeof(double),1,fp);
    fwrite(&(parent->relevance),sizeof(double),1,fp);
    l = L->first;
    while (l!=NULL){
      r = (struct region *) l->item;
      if (verbose){ printf(" %% %% region label %d\n",r->label);}
      fwrite(&(r->label),sizeof(int),1,fp);
      l=l->child;}
    if (posorpre){ fwrite(&minusone,sizeof(int),1,fp);}
    else{ fwrite(&minustwo,sizeof(int),1,fp);}}
  if (l0->kidl!=NULL){ ptreedump_helper(posorpre,parent,l0->kidl,fp,L,keepzeroweights);}
  if (l0->item!=NULL){ pn=(struct pnode *)l0->item; if (keepzeroweights || (!keepzeroweights && pn->weight!=0)){ ptreedump_helper(posorpre,pn,pn->childllitem,fp,L,keepzeroweights);}}
  if (l0->kidr!=NULL){ ptreedump_helper(posorpre,parent,l0->kidr,fp,L,keepzeroweights);}
  if (firstentry){ llistkilllast(L);}
}

struct ptree * ptreadback(char *filename,int region_type)
{
  int verbose=0;
  char filename2[256];
  FILE *fp=NULL;
  struct ptree *p=NULL;
  int nr=0;
  struct pnode *pn=NULL;
  struct llitem *l0=NULL,*l1=NULL;
  int readout=0;
  double weight=0,relevance=0;
  struct llist *L=NULL;
  struct litem *l=NULL;
  int *label=NULL;
  int event_threshold = GLOBAL_PTREE_EVENT_THRESHOLD;
  double event_within = GLOBAL_PTREE_EVENT_WITHIN;
  int exit_flag=0,depth=0;
  if (verbose){ printf(" %% \n");}
  if (verbose){ printf(" %% [entering ptreadback] with filename %s\n",filename==NULL ? "<null>" : filename);}
  if (filename==NULL){ sprintf(filename2,"./ptree_%srecord",GLOBAL_STRING_2);} else{ sprintf(filename2,"%s",filename);}
  if (verbose){ printf(" %% decided on filename %s\n",filename2);}
  if ((fp=fopen(filename2,"r"))==NULL){ printf("warning, cannot read %s in ptreadback\n",filename2);}
  else{
    p = (struct ptree *) tmalloc(sizeof(struct ptree));
    fread(&(p->nregions),sizeof(int),1,fp);
    fread(&(p->nlegs),sizeof(int),1,fp);
    fread(&(p->legtime),sizeof(int),1,fp);
    fread(&(p->length),sizeof(int),1,fp);
    fread(&(p->update_every),sizeof(double),1,fp);
    fread(&(p->total_time),sizeof(double),1,fp);
    if (verbose){ printf(" %% nregions %d, nlegs %d, legtime %d, length %d, update_every %f, total_time %f\n",p->nregions,p->nlegs,p->legtime,p->length,p->update_every,p->total_time);}
    p->regionra = (struct region **) tcalloc(p->nregions,sizeof(struct region*));
    if (verbose){ printf(" %% now making regions\n");}
    regionramake(p,event_within,event_threshold,region_type);
    assert(p->length>p->nlegs*p->legtime);
    p->eventra = (struct llitem **) tmalloc(p->length*sizeof(struct llitem *)); for (nr=0;nr<p->length;nr++){ p->eventra[nr]=llitemmake();}
    p->tab=0;
    p->update_last=GLOBAL_TI;
    if (verbose){ printf(" %% now making pretree and postree\n");}
    p->pretree=llitemmake();
    p->postree=llitemmake();
    p->wh=histmake(1,+1,-1);
    p->rh=histmake(1,+1,-1);
    exit_flag=0;
    while(exit_flag!=3){
      if ((readout = fread(&weight,sizeof(double),1,fp))==1 && (readout = fread(&relevance,sizeof(double),1,fp))==1){
	if (verbose){ printf(" %% read weight %f, relevance %f\n",weight,relevance);}
	L=llistmake();exit_flag=0;
	while (!exit_flag){
	  label = (int *) tmalloc(sizeof(int));
	  readout = fread(label,sizeof(int),1,fp);
	  if (readout!=1){ tfree(label); label=NULL; exit_flag=3;}
	  else{ 
	    if (*label>=0 && *label<p->nregions){ litemadd(L,label); exit_flag=0;}
	    else if (*label==-1){ tfree(label); label=NULL; exit_flag=1;}
	    else if (*label==-2){ tfree(label); label=NULL; exit_flag=2;}
	    else{ printf("warning, funny terminator %d in ptreadback\n",*label);}}}
	if (L->length>0 && (exit_flag==1 || exit_flag==2)){ /* postree or pretree */ 
	  l=L->first;depth=1; if (exit_flag==1){ l0=p->postree;pn=NULL;} else{ l0=p->pretree;pn=NULL;}
	  while(l!=NULL && depth<L->length){ /* traverse through nonlast elements of llist *L */
	    label = (int *) l->item;
	    if ((l1=llitemaddorfind(0,l0,p->regionra[*label],&region2pnode_compare_label))==NULL){
	      l1 = llitemaddorfind(1,l0,pnodemake(pn,p->regionra[*label],0,0),&pnode2pnode_compare_label);}
	    pn = (struct pnode *) l1->item;
	    l0 = pn->childllitem;
	    l=l->child;
	    depth+=1;}
	  if (depth==L->length){ /* at end of llist *L, must add pnode */
	    label = (int *) l->item;
	    if ((l1=llitemaddorfind(0,l0,p->regionra[*label],&region2pnode_compare_label))==NULL){
	      l1 = llitemaddorfind(1,l0,pnodemake(pn,p->regionra[*label],weight,relevance),&pnode2pnode_compare_label);}
	    else{ pn = (struct pnode *)l1->item; pn->weight += weight; pn->relevance=relevance;}}}
	llisttfree3(L);L=NULL;}
      else{ if (verbose){ printf(" %% couldn't read weight and relevance\n");} exit_flag=3;}}
    if (fp!=stdout){ fclose(fp);}
    if (verbose){ printf(" %% now balancing pretree and postree\n");}
    llitembalance(p->pretree); p->pretree = llitemclimb(p->pretree);
    llitembalance(p->postree); p->postree = llitemclimb(p->postree);}
  return p;
}

void ptreeplusequals_starter(struct ptree *p0,struct ptree *p1)
{
  ptreeplusequals_helper(p0,NULL,p0->postree,p1,NULL,p1->postree);
  ptreeplusequals_helper(p0,NULL,p0->pretree,p1,NULL,p1->pretree);
  pnodebalance_starter(NULL,p0->postree);
  pnodebalance_starter(NULL,p0->pretree);
  p0->total_time += p1->total_time;
}

void ptreeplusequals_helper(struct ptree *p0,struct pnode *parent0,struct llitem *l0,struct ptree *p1,struct pnode *parent1,struct llitem *l1)
{
  /* *p0 += *p1 */
  int firstentry=0;
  struct pnode *pn0=NULL,*pn1=NULL;
  struct llitem *l2=NULL;
  if (parent0!=NULL){ assert(parent0->childllitem==l0);}
  if (parent1!=NULL && l1->parent==NULL){ /* first descent into parent1->childllitem */ 
    assert(parent1->childllitem==l1); firstentry=1;}
  if (l1->kidl!=NULL){ ptreeplusequals_helper(p0,parent0,l0,p1,parent1,l1->kidl);}
  if (l1->kidr!=NULL){ ptreeplusequals_helper(p0,parent0,l0,p1,parent1,l1->kidr);}
  if (l1->item!=NULL){ 
    pn1 = (struct pnode *)l1->item; 
    if ((l2 = llitemaddorfind(0,l0,p0->regionra[pn1->region->label],&region2pnode_compare_label))==NULL){
      l2 = llitemaddorfind(1,l0,pnodemake(parent0,p0->regionra[pn1->region->label],0,0),&pnode2pnode_compare_label);}
    pn0 = (struct pnode *)l2->item;
    ptreeplusequals_helper(p0,pn0,pn0->childllitem,p1,pn1,pn1->childllitem);}
  if (firstentry){ 
    parent0->weight += parent1->weight;
    parent0->relevance = parent1->relevance;}
}

struct ptree * ptreesetequalto_starter(struct ptree *p1)
{
  struct ptree *p=NULL;
  int nr=0;
  p = (struct ptree *) tmalloc(sizeof(struct ptree));
  p->nregions = p1->nregions;
  p->regionra = (struct region **) tcalloc(p->nregions,sizeof(struct region*));
  for (nr=0;nr<p->nregions;nr++){ 
    p->regionra[nr] = (struct region *) tmalloc(sizeof(struct region)); 
    p->regionra[nr]->label = nr; /* critical to have label nr match location in regionra */
    p->regionra[nr]->last_event = p1->regionra[nr]->last_event;
    p->regionra[nr]->event_within = p1->regionra[nr]->event_within;
    p->regionra[nr]->event_threshold = p1->regionra[nr]->event_threshold;
    p->regionra[nr]->neuronllist = llistcopy(p1->regionra[nr]->neuronllist);
    p->regionra[nr]->pn = p1->regionra[nr]->pn;}
  p->nlegs=p1->nlegs;
  p->legtime=p1->legtime;
  p->length = p1->length;
  p->eventra = (struct llitem **) tmalloc(p->length*sizeof(struct llitem *));
  for (nr=0;nr<p->length;nr++){ 
    p->eventra[nr]=llitemmake(); 
    llitemgrowllitem(p->eventra[nr],p1->eventra[nr],&region2region_compare_label); 
    llitembalance(p->eventra[nr]);
    p->eventra[nr]=llitemclimb(p->eventra[nr]);}
  p->tab=p1->tab;
  p->update_every=p1->update_every;
  p->update_last=p1->update_last;
  p->total_time=p1->total_time;
  p->gate_flag=p1->gate_flag;
  ptreemakepospre(p);
  p->wh = histmake(1,+1,-1);
  p->rh = histmake(1,+1,-1);
  ptreesetequalto_helper(p,NULL,p->postree,p1,NULL,p1->postree);
  ptreesetequalto_helper(p,NULL,p->pretree,p1,NULL,p1->pretree);
  pnodebalance_starter(NULL,p->postree);
  pnodebalance_starter(NULL,p->pretree);
  return p;
}

void ptreesetequalto_helper(struct ptree *p0,struct pnode *parent0,struct llitem *l0,struct ptree *p1,struct pnode *parent1,struct llitem *l1)
{
  /* *p0 = *p1 */
  int firstentry=0;
  struct pnode *pn0=NULL,*pn1=NULL;
  struct llitem *l2=NULL;
  if (parent0!=NULL){ assert(parent0->childllitem==l0);}
  if (parent1!=NULL && l1->parent==NULL){ /* first descent into parent1->childllitem */ 
    assert(parent1->childllitem==l1); firstentry=1;}
  if (l1->kidl!=NULL){ ptreesetequalto_helper(p0,parent0,l0,p1,parent1,l1->kidl);}
  if (l1->kidr!=NULL){ ptreesetequalto_helper(p0,parent0,l0,p1,parent1,l1->kidr);}
  if (l1->item!=NULL){ 
    pn1 = (struct pnode *)l1->item; 
    if ((l2 = llitemaddorfind(0,l0,p0->regionra[pn1->region->label],&region2pnode_compare_label))==NULL){
      l2 = llitemaddorfind(1,l0,pnodemake(parent0,p0->regionra[pn1->region->label],0,0),&pnode2pnode_compare_label);}
    pn0 = (struct pnode *)l2->item;
    ptreesetequalto_helper(p0,pn0,pn0->childllitem,p1,pn1,pn1->childllitem);}
  if (firstentry){ 
    parent0->weight = parent1->weight;
    parent0->relevance = parent1->relevance;}
}

struct ptree * ptreesubtptree_starter(struct ptree *px,struct ptree *py,int ratematch)
{
  struct llist *pnLx=NULL,*pnLy=NULL;
  struct litem *lx=NULL,*ly=NULL;
  struct pnode *pnx=NULL,*pny=NULL;
  double xy=0,yy=0,ynormalizer=0;
  struct ptree *pz=NULL;
  if (ratematch){
    pnLx=llistmake(); pnLy=llistmake();
    pnode2llist_starter(NULL,px->postree,1,0,pnLx); pnode2llist_starter(NULL,py->postree,1,0,pnLy);
    llistsort(pnLx->first,pnLx->last,pnLx->length,&pnode2pnode_compare_label);
    llistsort(pnLy->first,pnLy->last,pnLy->length,&pnode2pnode_compare_label);
    lx=pnLx->first; ly=pnLy->first; xy=0; yy=0;
    while (lx!=NULL && ly!=NULL){
      pnx=(struct pnode *)lx->item; pny=(struct pnode *)ly->item;
      xy += pnx->weight*pny->weight; yy += pny->weight*pny->weight;      
      lx=lx->child; ly=ly->child;}
    llisttfree(pnLx); pnLx=NULL; llisttfree(pnLy); pnLy=NULL;
    ynormalizer=xy/yy; if (!finite(ynormalizer)){ ynormalizer=1;}}
  else /* if (!ratematch) */{ ynormalizer=1;}
  pz = ptreesetequalto_starter(px);
  ptreesubtptree_helper(pz,NULL,pz->postree,py,NULL,py->postree,ynormalizer);
  ptreesubtptree_helper(pz,NULL,pz->pretree,py,NULL,py->pretree,ynormalizer);
  return pz;
}

void ptreesubtptree_helper(struct ptree *px,struct pnode *parentx,struct llitem *lx,struct ptree *py,struct pnode *parenty,struct llitem *ly,double ynormalizer)
{
  /* *px -= *py*ynormalizer */
  int firstentry=0;
  struct pnode *pnx=NULL,*pny=NULL;
  struct llitem *l2=NULL;
  if (parentx!=NULL){ assert(parentx->childllitem==lx);}
  if (parenty!=NULL && ly->parent==NULL){ /* first descent into parenty->childllitem */ 
    assert(parenty->childllitem==ly); firstentry=1;}
  if (ly->kidl!=NULL){ ptreesubtptree_helper(px,parentx,lx,py,parenty,ly->kidl,ynormalizer);}
  if (ly->kidr!=NULL){ ptreesubtptree_helper(px,parentx,lx,py,parenty,ly->kidr,ynormalizer);}
  if (ly->item!=NULL){ 
    pny = (struct pnode *)ly->item; 
    if ((l2 = llitemaddorfind(0,lx,px->regionra[pny->region->label],&region2pnode_compare_label))==NULL){
      l2 = llitemaddorfind(1,lx,pnodemake(parentx,px->regionra[pny->region->label],0,0),&pnode2pnode_compare_label);}
    pnx = (struct pnode *)l2->item;
    ptreesubtptree_helper(px,pnx,pnx->childllitem,py,pny,pny->childllitem,ynormalizer);}
  if (firstentry){ 
    parentx->weight -= parenty->weight*ynormalizer;
    parentx->relevance = 0;}
}

void pnodefrob_starter(struct pnode *parent,struct llitem *l0,int maxlevel,int level,double *wnorm,double *rnorm)
{
  /* returns frobenius norm */
  int firstentry=0;
  struct pnode *pn0=NULL;
  if (parent!=NULL && l0->parent==NULL){ /* first descent into parent->childllitem */ 
    assert(parent->childllitem==l0); firstentry=1;}
  if (l0->kidl!=NULL){ pnodefrob_starter(parent,l0->kidl,maxlevel,level,wnorm,rnorm);}
  if (l0->kidr!=NULL){ pnodefrob_starter(parent,l0->kidr,maxlevel,level,wnorm,rnorm);}
  if (l0->item!=NULL){ 
    if (maxlevel==-1 || level<maxlevel){ 
      pn0 = (struct pnode *)l0->item; 
      pnodefrob_starter(pn0,pn0->childllitem,maxlevel,level+1,wnorm,rnorm);}}
  if (firstentry){ 
    if (wnorm!=NULL){ *wnorm += parent->weight*parent->weight;}
    if (rnorm!=NULL){ *rnorm += parent->relevance*parent->relevance;}}
}

void ptreex2(struct ptree *p0,struct pnode *parent0,struct llitem *l0,struct ptree *p1,struct pnode *parent1,struct llitem *l1,int retain_self,int maxlevel,int level,double *output_x2,double *output_x2p,double *output_fabs,double *output_fabs2)
{
  /* computes x2 estimate with p0 expected and p1 sample 
     returns estimators for chi_squared, forced_finite_chi_squared, fabsenius, squared_fabsenius */
  int verbose=0,nr=0;
  char lvlchar[64];
  int firstentry=0;
  struct pnode *pn0=NULL,*pn1=NULL;
  struct llitem *l2=NULL;
  double weight_expected=0,weight_sample=0,weight_minimum=0;
  if (verbose){ sprintf(lvlchar," %%"); for (nr=0;nr<level;nr++){sprintf(lvlchar,"%s %%",lvlchar);}}
  if (verbose>1){ 
    printf(" %s [entering ptreex2]\n",lvlchar);
    printf(" %s p0 %d parent0 %d->%d l0 %d \n",lvlchar,(int)p0,(int)parent0,parent0==NULL ? 0 : (int)(parent0->childllitem),(int)l0);
    printf(" %s p1 %d parent1 %d->%d l1 %d\n",lvlchar,(int)p1,(int)parent1,parent1==NULL ? 0 : (int)(parent1->childllitem),(int)l1);
    printf(" %s retain_self %d maxlevel %d level %d\n",lvlchar,retain_self,maxlevel,level);
    if (output_x2!=NULL){ printf(" %s output_x2 %f\n",lvlchar,*output_x2);}
    if (output_x2p!=NULL){ printf(" %s output_x2p %f\n",lvlchar,*output_x2p);}
    if (output_fabs!=NULL){ printf(" %s output_fabs %f\n",lvlchar,*output_fabs);}
    if (output_fabs2!=NULL){ printf(" %s output_fabs2 %f\n",lvlchar,*output_fabs2);}}
  if (parent1!=NULL){ assert(parent1->childllitem==l1);}
  if (parent0!=NULL && l0->parent==NULL){ /* first descent into parent0->childllitem */ 
    if (verbose>1){ printf(" %s first descent into parent0->childllitem, setting firstentry to 1\n",lvlchar);}
    assert(parent0->childllitem==l0); firstentry=1;}
  if (l0->kidl!=NULL){ 
    if (verbose>1){ printf(" %s l0->kidl found, moving left to new llitem l0->kidl=%d\n",lvlchar,(int)(l0->kidl));}
    ptreex2(p0,parent0,l0->kidl,p1,parent1,l1,retain_self,maxlevel,level,output_x2,output_x2p,output_fabs,output_fabs2);}
  if (l0->kidr!=NULL){ 
    if (verbose>1){ printf(" %s l0->kidr found, moving right to new llitem l0->kidr=%d\n",lvlchar,(int)(l0->kidr));}
    ptreex2(p0,parent0,l0->kidr,p1,parent1,l1,retain_self,maxlevel,level,output_x2,output_x2p,output_fabs,output_fabs2);}
  if (l0->item!=NULL){ 
    if (verbose>1){ printf(" %s l0->item %d exists\n",lvlchar,(int)(struct pnode *)(l0->item));}
    pn0 = (struct pnode *)l0->item; 
    if (maxlevel==-1 || level<maxlevel){
      if (verbose>1){ printf(" %s level %d less than maxlevel %d, continuing...\n",lvlchar,level,maxlevel);}
      if (l1!=NULL && (l2=llitemaddorfind(0,l1,p1->regionra[pn0->region->label],&region2pnode_compare_label))!=NULL){
	pn1 = (struct pnode *)l2->item;
	if (verbose>1){ printf(" %s pn0 label %d matches pn1 label %d, descending...\n",lvlchar,pn0->region->label,pn1->region->label);}
	ptreex2(p0,pn0,pn0->childllitem,p1,pn1,pn1->childllitem,retain_self,maxlevel,level+1,output_x2,output_x2p,output_fabs,output_fabs2);}
      else{
	pn1 = NULL;
	if (verbose>1){ printf(" %s pn0 label %d not found within l1, descending...\n",lvlchar,pn0->region->label);}
	ptreex2(p0,pn0,pn0->childllitem,p1,NULL,NULL,retain_self,maxlevel,level+1,output_x2,output_x2p,output_fabs,output_fabs2);}}}
  if (firstentry){ 
    if (!retain_self){ 
      weight_expected = (parent0->weight - (parent1==NULL ? 0 : parent1->weight))*p1->total_time/(p0->total_time-p1->total_time); 
      weight_minimum = 1.0*p1->total_time/(p0->total_time-p1->total_time);}
    else /* if (retain_self) */{
      weight_expected = parent0->weight*p1->total_time/p0->total_time; weight_minimum = 1.0*p1->total_time/p0->total_time;}
    if (verbose>0){ printf(" %s first entry, and weight_expected = %f...",lvlchar,weight_expected);}
    if (parent1!=NULL){ 
      if (verbose>0){ printf(" parent1 exists, so using sample weight %f\n",parent1->weight);}
      weight_sample = parent1->weight;}
    else{       
      if (verbose>0){ printf(" parent1 does not exist, so using sample weight 0\n");}
      weight_sample = 0;}
    if (output_x2!=NULL){ *output_x2 += pow(weight_sample - weight_expected,2)/pow((weight_expected > 0 ? weight_expected : 1),2);}
    if (output_x2p!=NULL){ *output_x2p += pow(weight_sample - weight_expected,2)/maximum(weight_minimum,weight_expected);}
    if (output_fabs!=NULL){ *output_fabs += fabs(weight_sample - weight_expected);}
    if (output_fabs2!=NULL){ *output_fabs2 += pow(weight_sample - weight_expected,2);}}
}

void ptreerelent_breadth(struct ptree *p0,struct pnode *parent0,struct llitem *l0,struct ptree *p1,struct pnode *parent1,struct llitem *l1,int retain_self,int maxlevel,int minlevel,int level,double *output,double *rnormalizer,double *onormalizer)
{
  /* computes "breadthwise" relative entropy with p0 expected and p1 observed 
     this assumes that rnormalizer and onormalizer are both maxlevel long, 
     and contain the mean weights at each level */
  int verbose=0,nr=0;
  char lvlchar[64];
  int firstentry=0;
  struct pnode *pn0=NULL,*pn1=NULL;
  struct llitem *l2=NULL;
  double weight_expected=0,weight_sample=0;
  if (verbose){ sprintf(lvlchar," %%"); for (nr=0;nr<level;nr++){sprintf(lvlchar,"%s %%",lvlchar);}}
  if (verbose>1){ 
    printf(" %s [entering ptreerelent_breadth]\n",lvlchar);
    printf(" %s p0 %d parent0 %d->%d l0 %d \n",lvlchar,(int)p0,(int)parent0,parent0==NULL ? 0 : (int)(parent0->childllitem),(int)l0);
    printf(" %s p1 %d parent1 %d->%d l1 %d\n",lvlchar,(int)p1,(int)parent1,parent1==NULL ? 0 : (int)(parent1->childllitem),(int)l1);
    printf(" %s retain_self %d maxlevel %d minlevel %d level %d\n",lvlchar,retain_self,maxlevel,minlevel,level);
    if (output!=NULL){ printf(" %s output %f\n",lvlchar,*output);}
    printf(" rnormalizer:\n"); raprintf(rnormalizer,"double",1,maxlevel,lvlchar);
    printf(" onormalizer:\n"); raprintf(onormalizer,"double",1,maxlevel,lvlchar);}
  if (parent1!=NULL){ assert(parent1->childllitem==l1);}
  if (parent0!=NULL && l0->parent==NULL){ /* first descent into parent0->childllitem */ 
    if (verbose>1){ printf(" %s first descent into parent0->childllitem, setting firstentry to 1\n",lvlchar);}
    assert(parent0->childllitem==l0); firstentry=1;}
  if (l0->kidl!=NULL){ 
    if (verbose>1){ printf(" %s l0->kidl found, moving left to new llitem l0->kidl=%d\n",lvlchar,(int)(l0->kidl));}
    ptreerelent_breadth(p0,parent0,l0->kidl,p1,parent1,l1,retain_self,maxlevel,minlevel,level,output,rnormalizer,onormalizer);}
  if (l0->kidr!=NULL){ 
    if (verbose>1){ printf(" %s l0->kidr found, moving right to new llitem l0->kidr=%d\n",lvlchar,(int)(l0->kidr));}
    ptreerelent_breadth(p0,parent0,l0->kidr,p1,parent1,l1,retain_self,maxlevel,minlevel,level,output,rnormalizer,onormalizer);}
  if (l0->item!=NULL){ 
    if (verbose>1){ printf(" %s l0->item %d exists\n",lvlchar,(int)(struct pnode *)(l0->item));}
    pn0 = (struct pnode *)l0->item; 
    if (maxlevel==-1 || level<maxlevel){
      if (verbose>1){ printf(" %s level %d less than maxlevel %d, continuing...\n",lvlchar,level,maxlevel);}
      if (l1!=NULL && (l2=llitemaddorfind(0,l1,p1->regionra[pn0->region->label],&region2pnode_compare_label))!=NULL){
	pn1 = (struct pnode *)l2->item;
	if (verbose>1){ printf(" %s pn0 label %d matches pn1 label %d, descending...\n",lvlchar,pn0->region->label,pn1->region->label);}
	ptreerelent_breadth(p0,pn0,pn0->childllitem,p1,pn1,pn1->childllitem,retain_self,maxlevel,minlevel,level+1,output,rnormalizer,onormalizer);}
      else{
	pn1 = NULL;
	if (verbose>1){ printf(" %s pn0 label %d not found within l1, descending...\n",lvlchar,pn0->region->label);}
	ptreerelent_breadth(p0,pn0,pn0->childllitem,p1,NULL,NULL,retain_self,maxlevel,minlevel,level+1,output,rnormalizer,onormalizer);}}}
  if (firstentry){ if (minlevel==-1 || level>=minlevel){
    if (!retain_self){ 
      weight_expected = maximum(1,(parent0->weight - (parent1==NULL ? 0 : parent1->weight)))/maximum(1,rnormalizer[level-1]-onormalizer[level-1]);}
    else /* if (retain_self) */{
      weight_expected = maximum(1,parent0->weight)/maximum(1,rnormalizer[level-1]);}
    if (verbose>0){ printf(" %s first entry, and weight_expected = %f...",lvlchar,weight_expected);}
    if (parent1!=NULL){ 
      weight_sample = parent1->weight/maximum(1,onormalizer[level-1]);
      if (verbose>0){ printf(" parent1 exists, so using sample weight %f\n",weight_sample);}}
    else{       
      if (verbose>0){ printf(" parent1 does not exist, so using sample weight 0\n");}
      weight_sample = 0;}
    if (output!=NULL){ *output += weight_sample>0 ? weight_sample*log(weight_sample/weight_expected) : 0;}}}
}

void ptreerelent_depth(struct ptree *p0,struct pnode *parent0,struct llitem *l0,struct ptree *p1,struct pnode *parent1,struct llitem *l1,int retain_self,int maxlevel,int level,double *output)
{
  /* computes "depthwise" relative entropy with p0 expected and p1 observed */
  int verbose=0,nr=0;
  char lvlchar[64];
  int firstentry=0;
  struct pnode *pn0=NULL,*pn1=NULL;
  struct llitem *l2=NULL;
  double weight_expected_prior=0,weight_expected=0,weight_observed_prior=0,weight_observed=0;
  double p_expected_prior=0,p_expected=0,p_observed_prior=0,p_observed=0;
  int maxdepth = (maxlevel==-1 ? p0->nlegs : minimum(p0->nlegs,maxlevel));
  if (verbose){ sprintf(lvlchar," %%"); for (nr=0;nr<level;nr++){sprintf(lvlchar,"%s %%",lvlchar);}}
  if (verbose>1){ 
    printf(" %s [entering ptreerelent_depth]\n",lvlchar);
    printf(" %s p0 %d parent0 %d->%d l0 %d \n",lvlchar,(int)p0,(int)parent0,parent0==NULL ? 0 : (int)(parent0->childllitem),(int)l0);
    printf(" %s p1 %d parent1 %d->%d l1 %d\n",lvlchar,(int)p1,(int)parent1,parent1==NULL ? 0 : (int)(parent1->childllitem),(int)l1);
    printf(" %s retain_self %d maxlevel %d level %d (maxdepth %d)\n",lvlchar,retain_self,maxlevel,level,maxdepth);
    if (output!=NULL){ printf(" %s output %f\n",lvlchar,*output);}}
  if (parent1!=NULL){ assert(parent1->childllitem==l1);}
  if (parent0!=NULL && l0->parent==NULL){ /* first descent into parent0->childllitem */ 
    if (verbose>1){ printf(" %s first descent into parent0->childllitem, setting firstentry to 1\n",lvlchar);}
    assert(parent0->childllitem==l0); firstentry=1;}
  if (l0->kidl!=NULL){ 
    if (verbose>1){ printf(" %s l0->kidl found, moving left to new llitem l0->kidl=%d\n",lvlchar,(int)(l0->kidl));}
    ptreerelent_depth(p0,parent0,l0->kidl,p1,parent1,l1,retain_self,maxlevel,level,output);}
  if (l0->kidr!=NULL){ 
    if (verbose>1){ printf(" %s l0->kidr found, moving right to new llitem l0->kidr=%d\n",lvlchar,(int)(l0->kidr));}
    ptreerelent_depth(p0,parent0,l0->kidr,p1,parent1,l1,retain_self,maxlevel,level,output);}
  if (l0->item!=NULL){ 
    if (verbose>1){ printf(" %s l0->item %d exists\n",lvlchar,(int)(struct pnode *)(l0->item));}
    pn0 = (struct pnode *)l0->item; 
    if (maxlevel==-1 || level<maxlevel){
      if (verbose>1){ printf(" %s level %d less than maxlevel %d, continuing...\n",lvlchar,level,maxlevel);}
      if (l1!=NULL && (l2=llitemaddorfind(0,l1,p1->regionra[pn0->region->label],&region2pnode_compare_label))!=NULL){
	pn1 = (struct pnode *)l2->item;
	if (verbose>1){ printf(" %s pn0 label %d matches pn1 label %d, descending...\n",lvlchar,pn0->region->label,pn1->region->label);}
	ptreerelent_depth(p0,pn0,pn0->childllitem,p1,pn1,pn1->childllitem,retain_self,maxlevel,level+1,output);}
      else{
	pn1 = NULL;
	if (verbose>1){ printf(" %s pn0 label %d not found within l1, descending...\n",lvlchar,pn0->region->label);}
	ptreerelent_depth(p0,pn0,pn0->childllitem,p1,NULL,NULL,retain_self,maxlevel,level+1,output);}}}
  if (firstentry){ 
    if (parent0->parent==NULL){ 
      assert(level==1);
      weight_expected_prior = p0->total_time;
      if (verbose>1){ printf(" %s %d==1 level event, no parent, must be a firing rate, weight_expected_prior=%f\n",lvlchar,level,weight_expected_prior);}}
    else{
      weight_expected_prior = parent0->parent->weight;
      if (verbose>1){ printf(" %s %d>1 level event, parent exists, weight_expected_prior=%f\n",lvlchar,level,weight_expected_prior);}}
    if (level==1){
      weight_observed_prior=p1->total_time;
      if (verbose>1){ printf(" %s level 1, must use weight_observed_prior=%f\n",lvlchar,weight_observed_prior);}}
    else{
      if (parent1!=NULL){
	assert(parent1->parent!=NULL);
	weight_observed_prior=parent1->parent->weight;
	if (verbose>1){ printf(" %s level %d, parent exists, weight_observed_prior=%f\n",lvlchar,level,weight_observed_prior);}}
      else{ 
	weight_observed_prior=0;
	if (verbose>1){ printf(" %s level %d, parent doesn't exist, weight_observed_prior=%f\n",lvlchar,level,weight_observed_prior);}}}
    weight_expected = parent0->weight;
    if (verbose>1){ printf(" %s first entry, weight_expected=%f...",lvlchar,weight_expected);}
    if (parent1!=NULL){ 
      weight_observed = parent1->weight;
      if (verbose>1){ printf(" %s parent1 exists, so weight_observed %f\n",lvlchar,weight_observed);}}
    else{       
      weight_observed = 0;
      if (verbose>1){ printf(" %s parent1 does not exist, so weight_observed %f\n",lvlchar,weight_observed);}}
    if (!retain_self){
      p_expected = minimum(1,maximum(0,(weight_expected-weight_observed)/(p0->total_time-p1->total_time)));
      p_expected_prior = minimum(1,maximum(0,((weight_expected_prior-weight_observed_prior)-(weight_expected-weight_observed))/(p0->total_time-p1->total_time)));
      p_observed = minimum(1,maximum(0,weight_observed/p1->total_time));
      p_observed_prior = minimum(1,maximum(0,(weight_observed_prior-weight_observed)/p1->total_time));}
    else /* if (retain_self) */{
      p_expected = minimum(1,maximum(0,weight_expected/p0->total_time));
      p_expected_prior = minimum(1,maximum(0,(weight_expected_prior-weight_expected)/p0->total_time));
      p_observed = minimum(1,maximum(0,weight_observed/p1->total_time));
      p_observed_prior = minimum(1,maximum(0,(weight_observed_prior-weight_observed)/p1->total_time));}
    if (verbose>1){ printf(" %s probabilities e %f ep %f o %f op %f,\n",lvlchar,p_expected,p_expected_prior,p_observed,p_observed_prior);}
    if (p_observed_prior>0 && p_expected_prior>0){ if (output!=NULL){ 
      *output += pow(p0->nregions,maxdepth-level)*p_observed_prior*log(p_observed_prior/p_expected_prior);}}
    if (level==maxdepth){ if (p_observed>0 && p_expected>0){ if (output!=NULL){ 
      *output += p_observed*log(p_observed/p_expected);}}}
    else if (level<maxdepth){ if (p_observed>0 && p_expected>0){ if (output!=NULL){ 
      *output += pow(p0->nregions,maxdepth-(level+1))*(p0->nregions-llitemlength(l1))*p_observed*log(p_observed/p_expected);}}}}
}

void pnodeclear_starter(struct pnode *parent,struct llitem *l0)
{
  int firstentry=0;
  struct pnode *pn=NULL;
  if (parent!=NULL && l0->parent==NULL){ /* first descent into parent->childllitem */ 
    assert(parent->childllitem==l0); firstentry=1;}
  if (l0->kidl!=NULL){ pnodeclear_starter(parent,l0->kidl);}
  if (l0->kidr!=NULL){ pnodeclear_starter(parent,l0->kidr);}
  if (l0->item!=NULL){ pn=(struct pnode *)l0->item; pnodeclear_starter(pn,pn->childllitem);}
  if (firstentry){ if (parent!=NULL){ parent->weight=0; parent->relevance=0;}}
}

void pnodehist_starter(struct pnode *parent,struct llitem *l0,int maxlevel,int level,struct hist *histw,struct hist *histr)
{
  /* puts weights and relevances into histw and histr */
  int firstentry=0;
  struct pnode *pn=NULL;
  if (parent!=NULL && l0->parent==NULL){ /* first descent into parent->childllitem */ 
    assert(parent->childllitem==l0); firstentry=1;}
  if (l0->kidl!=NULL){ pnodehist_starter(parent,l0->kidl,maxlevel,level,histw,histr);}
  if (l0->kidr!=NULL){ pnodehist_starter(parent,l0->kidr,maxlevel,level,histw,histr);}
  if (l0->item!=NULL){
    if (maxlevel==-1 || level<maxlevel){ 
      pn=(struct pnode *)l0->item; 
      pnodehist_starter(pn,pn->childllitem,maxlevel,level+1,histw,histr);}}
  if (firstentry){ if (parent!=NULL){ 
    if (histw!=NULL){ histadd(histw,parent->weight,1);} if (histr!=NULL){ histadd(histr,parent->relevance,1);}}}
}

void pnodestats_starter(struct pnode *parent,struct llitem *l0,int maxlevel,int minlevel,int level,int *nelements,double *wmax,double *wmin,double *wmean,double *wstdev,double *rmax,double *rmin,double *rmean,double *rstdev)
{
  /* start out with nelements==NULL, stdev only accumulates if mean!=NULL */
  int firstentry=0,veryfirstentry=0,madenelements_flag=0;
  struct pnode *pn=NULL;
  if (parent==NULL && l0->parent==NULL){ /* very first descent into ptree */
    if (nelements==NULL){ madenelements_flag=1; nelements=(int *)tmalloc(sizeof(int));} else{ madenelements_flag=0;} *nelements=0;
    if (wmin!=NULL){ *wmin=0;} if (wmax!=NULL){ *wmax=0;} if (wmean!=NULL){ *wmean=0; if (wstdev!=NULL){ *wstdev=0;}}
    if (rmin!=NULL){ *rmin=0;} if (rmax!=NULL){ *rmax=0;} if (rmean!=NULL){ *rmean=0; if (rstdev!=NULL){ *rstdev=0;}}
    veryfirstentry=1;}
  if (parent!=NULL && l0->parent==NULL){ /* first descent into parent->childllitem */ 
    assert(parent->childllitem==l0); assert(nelements!=NULL); firstentry=1;}
  if (l0->kidl!=NULL){ pnodestats_starter(parent,l0->kidl,maxlevel,minlevel,level,nelements,wmax,wmin,wmean,wstdev,rmax,rmin,rmean,rstdev);}
  if (l0->kidr!=NULL){ pnodestats_starter(parent,l0->kidr,maxlevel,minlevel,level,nelements,wmax,wmin,wmean,wstdev,rmax,rmin,rmean,rstdev);}
  if (l0->item!=NULL){
    if (maxlevel==-1 || level<maxlevel){ 
      pn=(struct pnode *)l0->item; 
      pnodestats_starter(pn,pn->childllitem,maxlevel,minlevel,level+1,nelements,wmax,wmin,wmean,wstdev,rmax,rmin,rmean,rstdev);}}
  if (firstentry){ if (parent!=NULL){ if (nelements!=NULL){ if (minlevel==-1 || level>=minlevel){
    if (*nelements==0){ if (wmax!=NULL){ *wmax=parent->weight;} if (wmin!=NULL){ *wmin=parent->weight;}}
    else if (*nelements>0){ if (wmax!=NULL){ *wmax=maximum(*wmax,parent->weight);} if (wmin!=NULL){ *wmin=minimum(*wmin,parent->weight);}}
    if (wmean!=NULL){ *wmean += parent->weight; if (wstdev!=NULL){ *wstdev += pow(parent->weight,2);}}
    if (*nelements==0){ if (rmax!=NULL){ *rmax=parent->relevance;} if (rmin!=NULL){ *rmin=parent->relevance;}}
    else if (*nelements>0){ if (rmax!=NULL){ *rmax=maximum(*rmax,parent->relevance);} if (rmin!=NULL){ *rmin=minimum(*rmin,parent->relevance);}}
    if (rmean!=NULL){ *rmean += parent->relevance; if (rstdev!=NULL){ *rstdev += pow(parent->relevance,2);}}
    *nelements += 1;}}}}
  if (veryfirstentry){ if (nelements!=NULL){ 
    if (wmean!=NULL){ *wmean /= maximum(1.0,(double)(*nelements)); if (wstdev!=NULL){ *wstdev /= maximum(1.0,(double)(*nelements)); *wstdev -= pow(*wmean,2);}}
    if (rmean!=NULL){ *rmean /= maximum(1.0,(double)(*nelements)); if (rstdev!=NULL){ *rstdev /= maximum(1.0,(double)(*nelements)); *rstdev -= pow(*rmean,2); *rstdev = sqrt(*rstdev);}}
    if (madenelements_flag==1){ tfree(nelements);nelements=NULL;}}}
}

void pnode2llist_starter(struct pnode *parent,struct llitem *l0,int maxlevel,int level,struct llist *L)
{
  /* start out with empty llist *L */
  int firstentry=0;
  struct pnode *pn=NULL;
  if (parent!=NULL && l0->parent==NULL){ /* first descent into parent->childllitem */ 
    assert(parent->childllitem==l0); firstentry=1;}
  if (l0->kidl!=NULL){ pnode2llist_starter(parent,l0->kidl,maxlevel,level,L);}
  if (l0->kidr!=NULL){ pnode2llist_starter(parent,l0->kidr,maxlevel,level,L);}
  if (l0->item!=NULL){
    if (maxlevel==-1 || level<maxlevel){ 
      pn=(struct pnode *)l0->item; 
      pnode2llist_starter(pn,pn->childllitem,maxlevel,level+1,L);}}
  if (firstentry){ if (parent!=NULL){ litemadd(L,parent);}}
}

void pnode2hist_starter(struct pnode *parent,struct llitem *l0,int maxlevel,int level,struct hist *wh,struct hist *rh)
{
  /* adds weights to *wh and relevances to *rh */
  int firstentry=0;
  struct pnode *pn=NULL;
  if (parent!=NULL && l0->parent==NULL){ /* first descent into parent->childllitem */ 
    assert(parent->childllitem==l0); firstentry=1;}
  if (l0->kidl!=NULL){ pnode2hist_starter(parent,l0->kidl,maxlevel,level,wh,rh);}
  if (l0->kidr!=NULL){ pnode2hist_starter(parent,l0->kidr,maxlevel,level,wh,rh);}
  if (l0->item!=NULL){
    if (maxlevel==-1 || level<maxlevel){ 
      pn=(struct pnode *)l0->item; 
      pnode2hist_starter(pn,pn->childllitem,maxlevel,level+1,wh,rh);}}
  if (firstentry){ if (parent!=NULL){ 
    if (wh!=NULL){ histadd(wh,parent->weight,1);} 
    if (rh!=NULL && parent->relevance!=0){ histadd(rh,parent->relevance,1);}}}
}

void pnodetimesd_starter(struct pnode *parent,struct llitem *l0,double weightx,double relevancex)
{
  int firstentry=0;
  struct pnode *pn=NULL;
  if (parent!=NULL && l0->parent==NULL){ /* first descent into parent->childllitem */ 
    assert(parent->childllitem==l0); firstentry=1;}
  if (l0->kidl!=NULL){ pnodetimesd_starter(parent,l0->kidl,weightx,relevancex);}
  if (l0->kidr!=NULL){ pnodetimesd_starter(parent,l0->kidr,weightx,relevancex);}
  if (l0->item!=NULL){ pn=(struct pnode *)l0->item; pnodetimesd_starter(pn,pn->childllitem,weightx,relevancex);}
  if (firstentry){ if (parent!=NULL){ parent->weight*=weightx; parent->relevance*=relevancex;}}
}

void pnodefabswr_starter(struct pnode *parent,struct llitem *l0,double weightf,double relevancef)
{
  int firstentry=0;
  struct pnode *pn=NULL;
  if (parent!=NULL && l0->parent==NULL){ /* first descent into parent->childllitem */ 
    assert(parent->childllitem==l0); firstentry=1;}
  if (l0->kidl!=NULL){ pnodefabswr_starter(parent,l0->kidl,weightf,relevancef);}
  if (l0->kidr!=NULL){ pnodefabswr_starter(parent,l0->kidr,weightf,relevancef);}
  if (l0->item!=NULL){ pn=(struct pnode *)l0->item; pnodefabswr_starter(pn,pn->childllitem,weightf,relevancef);}
  if (firstentry){ if (parent!=NULL){ 
    if (weightf){ parent->weight = fabs(parent->weight);}
    if (relevancef){ parent->relevance = fabs(parent->relevance);}}}
}

double ptreex2_match(int nE,int nS,double *x1ra,int *minindex)
{
  /* checks to see if x1ra[nS] is the minimum value
     if x1ra[nS] matches nmin of the values, we pass back 1.0/nmin; */
  int ne=0,nmin=0;
  double min=0,wtol=0.000000001,wr=0;
  stats("double",x1ra,nE,NULL,&min,NULL,NULL);
  nmin=0; for (ne=0;ne<nE;ne++){ if (fabs(x1ra[ne]-min)<wtol){ nmin += 1; if (minindex!=NULL){ *minindex=ne;}}}
  if (fabs(x1ra[nS]-min)<wtol){ wr = 1.0/(double)nmin; /* (nmin==1); */} else{ wr=0;}
  if (minindex!=NULL && nmin>1){ *minindex=-1;}
  return wr;
}

double ptreex1_smatch_old(int nE,int nS,double *we,double wa,int *minindex)
{
  /* additive matching criterion for nE expected weights *we and actual weight wa, given that the actual index should be nS 
     if wa matches nmin of the maxima of *we, we pass back (nmin==1) */
  int ne=0,nmin=0;
  double *wf=NULL,min=0,wtol=0.000000001,wr=0;
  wf = (double *) tcalloc(nE,sizeof(double));
  for (ne=0;ne<nE;ne++){ wf[ne] = fabs(we[ne]-wa);}
  stats("double",wf,nE,NULL,&min,NULL,NULL);
  nmin=0; for (ne=0;ne<nE;ne++){ if (fabs(wf[ne]-min)<wtol){ nmin += 1; if (minindex!=NULL){ *minindex=ne;}}}
  if (fabs(wf[nS]-min)<wtol){ wr=(nmin==1);} else{ wr=0;}
  tfree(wf); wf=NULL;
  if (minindex!=NULL && nmin>1){ *minindex=-1;}
  return wr;
}

double ptreex1_smatch(int nE,int nS,double *we,double wa,int *minindex)
{
  /* additive matching criterion for nE expected weights *we and actual weight wa, given that the actual index should be nS 
     if wa matches nmin of the maxima of *we, we pass back 1.0/nmin */
  int ne=0,nmin=0;
  double *wf=NULL,min=0,wtol=0.0000001,wr=0;
  wf = (double *) tcalloc(nE,sizeof(double));
  for (ne=0;ne<nE;ne++){ wf[ne] = fabs(we[ne]-wa);}
  stats("double",wf,nE,NULL,&min,NULL,NULL);
  nmin=0; for (ne=0;ne<nE;ne++){ if (fabs(wf[ne]-min)<wtol){ nmin += 1; if (minindex!=NULL){ *minindex=ne;}}}
  if (fabs(wf[nS]-min)<wtol){ wr=1.0/(double)nmin;} else{ wr=0;}
  tfree(wf); wf=NULL;
  if (minindex!=NULL && nmin>1){ *minindex=-1;}
  return wr;
}

double ptreex1_xmatch(int nE,int nS,double *we,double wa,int *maxindex)
{
  /* multiplicative matching criterion for nE expected weights *we and actual weight wa, given that the actual index should be nS */
  int ne=0,nmax=0;
  double *wf=NULL,max=0,wtol=0.0000001,wr=0;
  wf = (double *) tcalloc(nE,sizeof(double));
  for (ne=0;ne<nE;ne++){ wf[ne] = ((we[ne]!=0 && wa!=0) ? minimum(we[ne]/wa,wa/we[ne]) : 0);}
  stats("double",wf,nE,&max,NULL,NULL,NULL);
  nmax=0; for (ne=0;ne<nE;ne++){ if (fabs(wf[ne]-max)<wtol){ nmax += 1; if (maxindex!=NULL){ *maxindex=ne;}}}
  if (fabs(wf[nS]-max)<wtol){ wr=1.0/(double)nmax;} else{ wr=0;}
  tfree(wf); wf=NULL;
  if (maxindex!=NULL && nmax>1){ *maxindex=-1;}
  return wr;
}

void ptreex1(int nE,struct ptree **pE,struct pnode **parentE,struct llitem **lE,int nS,struct ptree *pS,struct pnode *parentS,struct llitem *lS,struct ptree *pR,struct pnode *parentR,struct llitem *lR,int maxlevel,int level)
{
  /* computes individual x1 estimate of sample ptree (*pS) against array of nE expected ptrees (**pE),
     output stored in record ptree (*pR), 
     for which it is assumed that the structure of *pR is the union of all the **pE (or more precisely, all possible *pS) */
  int verbose=0;
  char lvlchar[64];
  int firstentry=0,ne=0;
  struct pnode *pnR=NULL,*pnS=NULL,**pnE=NULL;
  struct llitem *lS2=NULL,**lE2=NULL;
  double *weight_expected=NULL,weight_actual=0;
  if (verbose){ sprintf(lvlchar," %%"); for (ne=0;ne<level;ne++){ sprintf(lvlchar,"%s %%",lvlchar);}}
  if (verbose>1){ printf(" %s \n",lvlchar);}
  if (verbose>1){ 
    printf(" %s [entering ptreex1] with nE %d and nS %d\n",lvlchar,nE,nS);
    for (ne=0;ne<nE;ne++){ printf(" %s pE[%d] %d parentE[%d] %d->%d lE[%d] %d \n",lvlchar,ne,(int)pE[ne],ne,(int)parentE[ne],parentE[ne]==NULL ? 0 : (int)(parentE[ne]->childllitem),ne,(int)lE[ne]);}
    printf(" %s pS %d parentS %d->%d lS %d\n",lvlchar,(int)pS,(int)parentS,parentS==NULL ? 0 : (int)(parentS->childllitem),(int)lS);
    printf(" %s pR %d parentR %d->%d lR %d\n",lvlchar,(int)pR,(int)parentR,parentR==NULL ? 0 : (int)(parentR->childllitem),(int)lR);
    printf(" %s maxlevel %d level %d\n",lvlchar,maxlevel,level);}
  if (parentS!=NULL){ assert(parentS->childllitem==lS);}
  for (ne=0;ne<nE;ne++){ if (parentE[ne]!=NULL){ assert(parentE[ne]->childllitem==lE[ne]);}}
  if (parentR!=NULL && lR->parent==NULL){ /* first descent into parentR->childllitem */
    assert(parentR->childllitem==lR);
    if (verbose>1){ printf(" %s first descent into parentR->childllitem, setting firstentry to 1\n",lvlchar);}
    firstentry=1;}
  if (lR->kidl!=NULL){ 
    if (verbose>1){ printf(" %s stepping left to lR->kidl %d\n",lvlchar,(int)(lR->kidl));}
    ptreex1(nE,pE,parentE,lE,nS,pS,parentS,lS,pR,parentR,lR->kidl,maxlevel,level);}
  if (lR->kidr!=NULL){
    if (verbose>1){ printf(" %s stepping right to lR->kidr %d\n",lvlchar,(int)(lR->kidr));}
    ptreex1(nE,pE,parentE,lE,nS,pS,parentS,lS,pR,parentR,lR->kidr,maxlevel,level);}
  if (lR->item!=NULL){
    if (verbose>1){ printf(" %s lR->item %d exists, delving deeper\n",lvlchar,(int)(lR->item));}
    pnR = (struct pnode *)lR->item;
    if (level < maxlevel){
      if (verbose>1){ printf(" %s level %d less than maxlevel %d, continuing...\n",lvlchar,level,maxlevel);}
      if (lS!=NULL && (lS2 = llitemaddorfind(0,lS,pS->regionra[pnR->region->label],&region2pnode_compare_label))!=NULL){
	pnS = (struct pnode *)lS2->item;}
      else{ pnS=NULL;}
      lS2 = (pnS==NULL ? NULL : pnS->childllitem);
      if (verbose>1){ 
	if (pnS!=NULL){ printf(" %s pnS label %d matches pnR label %d, ...\n",lvlchar,pnS->region->label,pnR->region->label);}
	else{ printf(" %s pnS label doesn't match pnR label %d, ...\n",lvlchar,pnR->region->label);}}
      pnE = (struct pnode **) tcalloc(nE,sizeof(struct pnode *));
      lE2 = (struct llitem **) tcalloc(nE,sizeof(struct llitem *));
      for (ne=0;ne<nE;ne++){
	if (lE[ne]!=NULL && (lE2[ne] = llitemaddorfind(0,lE[ne],pE[ne]->regionra[pnR->region->label],&region2pnode_compare_label))!=NULL){
	  pnE[ne] = (struct pnode *)lE2[ne]->item;}
	else{ pnE[ne] = NULL;}
	lE2[ne] = (pnE[ne]==NULL ? NULL : pnE[ne]->childllitem);
	if (verbose>1){ 
	  if (pnE[ne]!=NULL){ 
	    printf(" %s pnE[%d] label %d matches pnR label %d, ...\n",lvlchar,ne,pnE[ne]->region->label,pnR->region->label);}
	  else{ printf(" %s pnE[%d] label doesn't match pnR label %d, ...\n",lvlchar,ne,pnR->region->label);}}}
      if (verbose>1){ printf(" %s descending\n",lvlchar);}
      ptreex1(nE,pE,pnE,lE2,nS,pS,pnS,lS2,pR,pnR,pnR->childllitem,maxlevel,level+1);
      tfree(pnE); pnE=NULL;
      tfree(lE2); lE2=NULL;}}
  if (firstentry){
    if (verbose>1){ printf(" %s first entry for actual index %d,\n",lvlchar,nS);}
    weight_expected = (double *) tcalloc(nE,sizeof(double));
    for (ne=0;ne<nE;ne++){ if (parentE[ne]!=NULL){ weight_expected[ne] = parentE[ne]->weight*pS->total_time/pE[ne]->total_time;}}
    if (verbose>1){ printf(" %s index expcted weights:",lvlchar); for (ne=0;ne<nE;ne++){ printf(" %f",weight_expected[ne]);} printf("\n");}
    if (parentS!=NULL){ weight_actual = parentS->weight;} else{ weight_actual=0;}
    if (verbose>1){ printf(" %s index actual weight:  %f\n",lvlchar,weight_actual);}
    if (verbose>1){ printf(" %s computing best index\n",lvlchar);}
    parentR->weight += ptreex1_smatch(nE,nS,weight_expected,weight_actual,NULL);
    parentR->relevance += ptreex1_xmatch(nE,nS,weight_expected,weight_actual,NULL);
    if (verbose>1){ printf(" %s indexing... weight set to %f and relevance set to %f\n",lvlchar,parentR->weight,parentR->relevance);}
    tfree(weight_expected);weight_expected=NULL;}
}

void ptreex1sbias(int nE,struct ptree **pE,struct pnode **parentE,struct llitem **lE,int nS,struct ptree *pS,struct pnode *parentS,struct llitem *lS,struct ptree *pR,struct pnode *parentR,struct llitem *lR,int maxlevel,int level,int nX,double *percentages,double *sbias)
{
  /* calculates the biased sum-based x1 index estimator for sample *pS with actual index nS
     in comparison to nE expected ptrees **pE
     using reference ptree *pR which has been normalized to have entries between 0 and 1 (probabilities).
     Assumes that *percentages is an ascending array of nX elements, 
     and *sbias is some pregenerated array of nbases*nX*maxelevel elements
  */
  int verbose=0;
  char lvlchar[64];
  int firstentry=0,ne=0,nx=0,nl=0;
  struct pnode *pnR=NULL,*pnS=NULL,**pnE=NULL;
  struct llitem *lS2=NULL,**lE2=NULL;
  double *weight_expected=NULL,weight_actual=0;
  int minindex=0;
  double nadd=0;
  if (verbose){ sprintf(lvlchar," %%"); for (ne=0;ne<level;ne++){ sprintf(lvlchar,"%s %%",lvlchar);}}
  if (verbose>1){ printf(" %s \n",lvlchar);}
  if (verbose>1){ 
    printf(" %s [entering ptreex1sbias] with nE %d and nS %d and nX %d\n",lvlchar,nE,nS,nX);
    for (ne=0;ne<nE;ne++){ printf(" %s pE[%d] %d parentE[%d] %d->%d lE[%d] %d \n",lvlchar,ne,(int)pE[ne],ne,(int)parentE[ne],parentE[ne]==NULL ? 0 : (int)(parentE[ne]->childllitem),ne,(int)lE[ne]);}
    printf(" %s pS %d parentS %d->%d lS %d\n",lvlchar,(int)pS,(int)parentS,parentS==NULL ? 0 : (int)(parentS->childllitem),(int)lS);
    printf(" %s pR %d parentR %d->%d lR %d\n",lvlchar,(int)pR,(int)parentR,parentR==NULL ? 0 : (int)(parentR->childllitem),(int)lR);
    printf(" %s maxlevel %d level %d\n",lvlchar,maxlevel,level);
    printf(" %s percentages:",lvlchar); for (nx=0;nx<nX;nx++){ printf(" %0.3f",percentages[nx]);} printf("\n");}
  if (parentS!=NULL){ assert(parentS->childllitem==lS);}
  for (ne=0;ne<nE;ne++){ if (parentE[ne]!=NULL){ assert(parentE[ne]->childllitem==lE[ne]);}}
  if (parentR!=NULL && lR->parent==NULL){ /* first descent into parentR->childllitem */
    assert(parentR->childllitem==lR);
    if (verbose>1){ printf(" %s first descent into parentR->childllitem, setting firstentry to 1\n",lvlchar);}
    firstentry=1;}
  if (lR->kidl!=NULL){ 
    if (verbose>1){ printf(" %s stepping left to lR->kidl %d\n",lvlchar,(int)(lR->kidl));}
    ptreex1sbias(nE,pE,parentE,lE,nS,pS,parentS,lS,pR,parentR,lR->kidl,maxlevel,level,nX,percentages,sbias);}
  if (lR->kidr!=NULL){
    if (verbose>1){ printf(" %s stepping right to lR->kidr %d\n",lvlchar,(int)(lR->kidr));}
    ptreex1sbias(nE,pE,parentE,lE,nS,pS,parentS,lS,pR,parentR,lR->kidr,maxlevel,level,nX,percentages,sbias);}
  if (lR->item!=NULL){
    if (verbose>1){ printf(" %s lR->item %d exists, delving deeper\n",lvlchar,(int)(lR->item));}
    pnR = (struct pnode *)lR->item;
    if (level < maxlevel){
      if (verbose>1){ printf(" %s level %d less than maxlevel %d, continuing...\n",lvlchar,level,maxlevel);}
      if (lS!=NULL && (lS2 = llitemaddorfind(0,lS,pS->regionra[pnR->region->label],&region2pnode_compare_label))!=NULL){
	pnS = (struct pnode *)lS2->item;}
      else{ pnS=NULL;}
      lS2 = (pnS==NULL ? NULL : pnS->childllitem);
      if (verbose>1){ 
	if (pnS!=NULL){ printf(" %s pnS label %d matches pnR label %d, ...\n",lvlchar,pnS->region->label,pnR->region->label);}
	else{ printf(" %s pnS label doesn't match pnR label %d, ...\n",lvlchar,pnR->region->label);}}
      pnE = (struct pnode **) tcalloc(nE,sizeof(struct pnode *));
      lE2 = (struct llitem **) tcalloc(nE,sizeof(struct llitem *));
      for (ne=0;ne<nE;ne++){
	if (lE[ne]!=NULL && (lE2[ne] = llitemaddorfind(0,lE[ne],pE[ne]->regionra[pnR->region->label],&region2pnode_compare_label))!=NULL){
	  pnE[ne] = (struct pnode *)lE2[ne]->item;}
	else{ pnE[ne] = NULL;}
	lE2[ne] = (pnE[ne]==NULL ? NULL : pnE[ne]->childllitem);
	if (verbose>1){ 
	  if (pnE[ne]!=NULL){ 
	    printf(" %s pnE[%d] label %d matches pnR label %d, ...\n",lvlchar,ne,pnE[ne]->region->label,pnR->region->label);}
	  else{ printf(" %s pnE[%d] label doesn't match pnR label %d, ...\n",lvlchar,ne,pnR->region->label);}}}
      if (verbose>1){ printf(" %s descending\n",lvlchar);}
      ptreex1sbias(nE,pE,pnE,lE2,nS,pS,pnS,lS2,pR,pnR,pnR->childllitem,maxlevel,level+1,nX,percentages,sbias);
      tfree(pnE);
      tfree(lE2);}}
  if (firstentry){
    if (verbose>1){ printf(" %s first entry at level %d, actual index %d, element's value is %f\n",lvlchar,level,nS,parentR->weight);}
    if (parentR->weight>=percentages[0]){ /* not bad */
      weight_expected = (double *) tcalloc(nE,sizeof(double));
      for (ne=0;ne<nE;ne++){ if (parentE[ne]!=NULL){ weight_expected[ne] = parentE[ne]->weight*pS->total_time/pE[ne]->total_time;}}
      if (parentS!=NULL){ weight_actual = parentS->weight;} else{ weight_actual=0;}
      nadd = ptreex1_smatch(nE,nS,weight_expected,weight_actual,&minindex);
      if (verbose>1){ printf(" %s expected weights:",lvlchar); for (ne=0;ne<nE;ne++){ printf(" %0.5f",weight_expected[ne]);} printf("\n"); printf(" %s actual weight %0.2f\n",lvlchar,weight_actual); printf(" %s determined minindex %d\n",lvlchar,minindex);}
      /* if we want to only consider unique minimums */
      if (minindex!=-1){ 
        if (verbose>1){ printf(" %s unique minimum found at index %d, with reference %0.4f, sbias starting at:\n",lvlchar,minindex,parentR->weight); printf(" %s\n",lvlchar); for (nl=0;nl<maxlevel;nl++){ printf(" %s level %d\n",lvlchar,nl); raprintf(&(sbias[0+0*nE+nl*nX*nE]),"double",nE,nX,lvlchar);}}
	for (nx=0;nx<nX;nx++){ if (parentR->weight>=percentages[nx]){ for (nl=level-1;nl<maxlevel;nl++){ sbias[minindex+nx*nE+nl*nE*nX] += 1;}}}
        if (verbose>1){ printf(" %s now sbias:\n",lvlchar); printf(" %s\n",lvlchar); for (nl=0;nl<maxlevel;nl++){ printf(" %s level %d\n",lvlchar,nl); raprintf(&(sbias[0+0*nE+nl*nX*nE]),"double",nE,nX,lvlchar);}}}
      /* on the other hand, if nonstrict minimums won't affect the answer if nbases==2 */
      tfree(weight_expected);}}
}

void pnodeprintf_label(void *v)
{
  /* prints pn->region->label */
  struct pnode *pn=NULL;
  if (v!=NULL){ pn=(struct pnode *)v; printf("pnode %d label %d\n",(int)pn,pn->region->label);}
  else{ printf("NULL\n");}
}

void pnodeshizuffle_starter(struct pnode *parent,struct llitem *l0,struct ptree *p)
{
  int verbose=0;
  int firstentry=0;
  struct pnode *pn=NULL;
  struct llist *L=NULL;
  struct litem *l=NULL;
  struct pnode **permutation=NULL,*pntemp=NULL;
  int npasses=5;
  int nr=0,pr=0,index1=0,index2=0;
  if (verbose){ printf(" %% \n");}
  if (verbose){ printf(" %% [entering pnodeshizuffle_starter]\n");}
  if (parent!=NULL && l0->parent==NULL){ /* first descent into parent->childllitem */ 
    assert(parent->childllitem==l0); firstentry=1;}
  if (l0->kidl!=NULL){ 
    if (verbose){ printf(" %% passing left\n");} pnodeshizuffle_starter(parent,l0->kidl,p);}
  if (l0->kidr!=NULL){ 
    if (verbose){ printf(" %% passing right\n");} pnodeshizuffle_starter(parent,l0->kidr,p);}
  if (l0->item!=NULL){ 
    if (verbose){ printf(" %% passing down\n");} pn=(struct pnode *)l0->item; pnodeshizuffle_starter(pn,pn->childllitem,p);}
  if (firstentry){ 
    if (parent!=NULL && l0!=NULL && l0->item!=NULL){ 
      pr = parent->region->label;
      if (verbose){ printf(" %% firstentry, starting with parent region %d, and childllitem:\n",pr); llitemprintf(l0,&pnodeprintf_label);}
      L=llistmake();
      llitem2llist(l0,L); 
      llistsort(L->first,L->last,L->length,&pnode2pnode_compare_label); /* should be unnecessary */
      if (verbose){ printf(" %% constructed llist:\n"); l=L->first;while(l!=NULL){ pn=(struct pnode *)l->item;printf("%d\n",pn->region->label);l=l->child;}}
      permutation = (struct pnode **) tcalloc(p->nregions,sizeof(struct pnode *));
      nr=0;l=L->first;
      while(l!=NULL && nr<p->nregions){ 
	pn=(struct pnode *)l->item; 
	if (pn->region->label==nr){ permutation[nr]=pn; nr+=1; l=l->child;}
	else{ permutation[nr]=NULL; nr+=1;}}
      llisttfree(L);L=NULL;
      if (verbose){ printf(" %% constructed permutation:\n"); for (nr=0;nr<p->nregions;nr++){ if (permutation[nr]!=NULL){ printf("pnode %d label %d\n",(int)permutation[nr],permutation[nr]->region->label);} else{ printf("NULL\n");}}}
      for (nr=0;nr<npasses*p->nregions;nr++){
	index1=periodize((int)floor(p->nregions*rand01),0,p->nregions);
	index2=periodize((int)floor(p->nregions*rand01),0,p->nregions);
	if (index1!=pr && index2!=pr){
	  pntemp = permutation[index1];
	  permutation[index1] = permutation[index2];
	  permutation[index2] = pntemp;}}
      for (nr=0;nr<p->nregions;nr++){ if (permutation[nr]!=NULL){ permutation[nr]->region = p->regionra[nr];}}
      if (verbose){ printf(" %% transposed into permutation:\n"); for (nr=0;nr<p->nregions;nr++){ if (permutation[nr]!=NULL){ printf("pnode %d label %d\n",(int)permutation[nr],permutation[nr]->region->label);} else{ printf("NULL\n");}}}
      L=llistmake();
      for (nr=0;nr<p->nregions;nr++){ if (permutation[nr]!=NULL){ litemadd(L,permutation[nr]);}}
      if (verbose){ printf(" %% transposed into llist:\n"); l=L->first;while(l!=NULL){ pn=(struct pnode *)l->item;printf("%d\n",pn->region->label);l=l->child;}}
      llitemtfree(l0,NULL);l0=llitemmake();
      llist2llitem(L,l0,&pnode2pnode_compare_label);
      llisttfree(L);L=NULL;
      tfree(permutation);permutation=NULL;
      parent->childllitem=l0;
      llitembalance(parent->childllitem);parent->childllitem=llitemclimb(parent->childllitem);l0=parent->childllitem;
      if (verbose){ printf(" %% firstentry, ending with parent region %d, and childllitem:\n",pr); llitemprintf(parent->childllitem,&pnodeprintf_label);}}}
}

void ptreex2_disthist(double *distra,double *x2ra,int nbases)
{
  /* adds (1.0/(double)nmin) to the indices of distra which correspond to the minima of x2ra */
  double max=0,min=0,wtol=0.0000001;
  int nr=0,nmin=0;
  stats("double",x2ra,nbases,&max,&min,NULL,NULL);
  for (nr=0;nr<nbases;nr++){ if (fabs(x2ra[nr]-min)<wtol){ nmin += 1;}} 
  nmin=maximum(1,nmin);
  for (nr=0;nr<nbases;nr++){ if (fabs(x2ra[nr]-min)<wtol){ distra[nr] += (1.0/(double)nmin);}}
}

double logfactorial(double n)
{
/*   return 0.5*(log(2.0*n + 1.0/3.0)+log(PI)) + n*log(n>0 ? n : 1) - n; */
  return lgamma(n+1);
}

double binomial(double n,double N,double p)
{
  /* Probability that given N choices of 2 options (option 1 has probability p), option 1 occurs n times */
  /* approximate binomial distribution is used */
  int verbose=0;
  double n2=0;
  double wtol=0.0000001;
  double b=0;
  if (verbose){ printf(" %% binomial(%f,%f,%f)",n,N,p);}
  N=maximum(0,N);n=maximum(0,n);p=minimum(1,maximum(0,p));
  n2=N-n;
  if (n>=0 && n2>=0){
    if (fabs(N-n)<wtol){ b = exp(n*log(p));}
    else if (fabs(N-n2)<wtol){ b = exp(n2*log(1-p));}
    else /* if (n<N && n2<N) */{ b = exp(logfactorial(N) - logfactorial(n) - logfactorial(n2) + n*log(p) + n2*log(1-p));}}
  else /* if (n<0 || n2<0) */{ b = 0;}
  if (verbose){ printf("=%f\n",b);}
  return b;
}

double Pmode(double n,double N,int m,double *p)
{
  /* Probability that given N choices of m options (option i has probability p[i]), the mode occurs more than n times */
  /* GLOBAL_QUADRATURE_N, GLOBAL_QUADRATURE_Z, GLOBAL_QUADRATURE_W are used as quadrature nodes */
  /* for integrals spanning more than quadthresh integers we approximate
     \int_{0}^{n}d\theta binomial(\theta,N,p[0])*Pmode(n,N-\theta,m-1,&(p[1]))
     + \int_{n}^{N}d\theta binomial(\theta,N,p[0]) */
  /* my intuition says that p should be sorted in descending order */
  int verbose=0;
  int quadthresh=100;
  int nl=0;
  double sumprob=0;
  double icenter=0,iradius=0,ix=0;
  double integral1=0,integral2=0;
  if (verbose){ printf(" %% [entering Pmode]\n %% n=%f,N=%f,m=%d",n,N,m); raprintf(p,"double",1,m," p = ");}
  if (n>=N){ if (verbose){ printf(" %% n>=N, exiting\n");} return 0;}
  else /* if (n<N) */{
    if (m<1){ printf(" %% error, m=%d in Pmode\n",m);}
    else if (m==1){ if (verbose){ printf(" %% only one choice\n");} return 1;}
    else if (m>1){ 
      if (n<(N/(double)m)){ 
	if (verbose){ printf(" %% too few choices, mode at least N/m=%f\n",N/(double)m);}
	return 1;}
      else /* if n>=(N/(double)m) */{
	if (verbose){ printf(" %% many choices, mode may not be %f or more\n",n+1);}
	stats("double",p,m,NULL,NULL,&sumprob,NULL); sumprob *= m; if (sumprob<=0){ sumprob=1;}
	if (verbose){ printf(" %% probability sums to %f\n",sumprob);}
	if (n>quadthresh){
	  integral1=0;icenter=(n+0.5+0.0)/2.0;iradius=(n+0.5-0.0)/2.0;
	  for (nl=0;nl<GLOBAL_QUADRATURE_N;nl++){
	    ix = GLOBAL_QUADRATURE_Z[nl]*iradius + icenter;
	    integral1 += GLOBAL_QUADRATURE_W[nl]*binomial(ix,N,p[0]/sumprob)*Pmode(n,N-ix,m-1,&(p[1]));}
	  integral1 *= iradius;}
	else /* if domain is small */{
	  integral1=0;
	  for (nl=0;nl<=n;nl++){
	    integral1 += binomial(nl,N,p[0]/sumprob)*Pmode(n,N-nl,m-1,&(p[1]));}}
	if ((N-n-1)>quadthresh){
	  integral2=0;icenter=(N+n+0.5)/2.0;iradius=(N-n-0.5)/2.0;
	  for (nl=0;nl<GLOBAL_QUADRATURE_N;nl++){
	    ix = GLOBAL_QUADRATURE_Z[nl]*iradius + icenter;
	    integral2 += GLOBAL_QUADRATURE_W[nl]*binomial(ix,N,p[0]/sumprob);}
	  integral2 *= iradius;}
	else /* if domain is small */{ 
	  integral2=0;
	  for (nl=n+1;nl<=N;nl++){
	    integral2 += binomial(nl,N,p[0]/sumprob);}}
	if (verbose){ printf(" %% returning %f + %f\n",integral1,integral2);}
	return integral1+integral2;}}}
  return 0;
}

double multithresh_single(double *pra,int nbases,double threshold)
{
  /* returns minimum number of independent observations which guarantees (up to threshold) 
     that the most probable observation (with probability p[0]) is the mode */  
  /* GLOBAL_QUADRATURE_N, GLOBAL_QUADRATURE_Z, GLOBAL_QUADRATURE_W are used as quadrature nodes */
  /* calculate 
     \int_{N/m}^{N}d\theta binomial(\theta,N,p[0])*(1-Pmode(\theta,N-\theta,m-1,&(p[1]))) */
  int quadthresh=100;
  int nb=0,nl=0;
  struct llist *L=NULL;
  struct litem *l=NULL;
  double *pra2=NULL,sumprob=0;
  double icenter=0,iradius=0,ix=0;
  double integral=0;
  int nobs=0,nobsmax=100;
  L=llistmake();
  for (nb=0;nb<nbases;nb++){ litemadd(L,&(pra[nb]));}
  llistsort(L->first,L->last,L->length,&double_compare);
  pra2 = (double *) tcalloc(nbases,sizeof(double));
  l=L->last; nb=0; while (l!=NULL){ pra2[nb]=*(double *)l->item; l=l->parent; nb++;}
  stats("double",pra2,nbases,NULL,NULL,&sumprob,NULL); sumprob *= nbases; if (sumprob<=0){ sumprob=1;}
  llisttfree(L);L=NULL;
  nobs=0;
  do{
    nobs += 1;
    if (nobs>quadthresh){
      integral=0;icenter=((double)nobs/(double)nbases + (double)nobs)/2.0;iradius=((double)nobs-(double)nobs/(double)nbases)/2.0;
      for (nl=0;nl<GLOBAL_QUADRATURE_N;nl++){
	ix = GLOBAL_QUADRATURE_Z[nl]*iradius + icenter;
	integral += GLOBAL_QUADRATURE_W[nl]*binomial(ix,(double)nobs,pra2[0]/sumprob)*(1-Pmode(ix,(double)nobs-ix,nbases-1,&(pra2[1])));}
      integral *= iradius;}
    else /* if domain is small */{
      integral=0;
      for (nl=(int)floor(nobs/(double)nbases);nl<=nobs;nl++){
	integral += binomial(nl,nobs,pra2[0]/sumprob)*(1-Pmode(nl,nobs-nl,nbases-1,&(pra2[1])));}}}
  while (nobs<nobsmax && integral<threshold);
  tfree(pra2);pra2=NULL;
  return nobs;
}

double multithresh(double *distra,int nbases,double threshold)
{
  /* probabilities are normalized versions of distra[observed_base + actual_base*nbases]
     returns minimum number of independent observations which guarantees (up to threshold) observing actual_base as the mode 
     or -1 if not possible */
  int verbose=0;
  double mean=0,*pra=NULL,max=0,wtol=0.0000001,temp=0;
  int nr=0,nr2=0;
  double nobs=0,nobsmax=100;
  if (verbose){ printf(" %% [entering multithresh]\n");}
  pra = (double *) tcalloc(nbases*nbases,sizeof(double));
  for (nr=0;nr<nbases;nr++){ 
    stats("double",&(distra[0+nr*nbases]),nbases,NULL,NULL,&mean,NULL); 
    for (nr2=0;nr2<nbases;nr2++){ pra[nr2+nr*nbases] = distra[nr2+nr*nbases]/(mean*(double)nbases);}}
  if (verbose){ raprintf(distra,"double",nbases,nbases," %% distra "); raprintf(pra,"double",nbases,nbases," %% pra ");}
  nobs=0;
  for (nr=0;nr<nbases;nr++){
    stats("double",&(pra[0+nr*nbases]),nbases,&max,NULL,NULL,NULL);
    if (fabs(pra[nr+nr*nbases]-max)>wtol){ 
      if (verbose){ printf(" %% actual_base %d not possible\n",nr);} return nobsmax;}
    else{ 
      temp=pra[0+nr*nbases]; pra[0+nr*nbases]=pra[nr+nr*nbases]; pra[nr+nr*nbases]=temp; 
      if (verbose){ printf(" %% shuffled pra for %d\n",nr); raprintf(&(pra[0+nr*nbases]),"double",1,nbases," %% pra ");}
      nobs += (multithresh_single(&(pra[0+nr*nbases]),nbases,threshold))/(double)nbases;}}
  tfree(pra);pra=NULL;
  if (verbose){ printf(" %% returning %f\n",nobs);}
  return nobs;
}

void llistmedian(struct llist *L1,struct llist *L2,int (*compare)(void *,void *),int sort_flag,void **iout,int *L1less,int *L2more)
{
  /* assumes that L1 and L2 are (disjoint and distinct) llists of comparable objects
     Here we return the address iout of the item (either within L1 or L2) which is the comparative median, 
     that is, which maximizes the number of elements of L1 that are < iout, and of L2 which are >= iout 
     these numbers are stored as *L1less and *L2more */
  int verbose=0;
  struct llist *L1t=NULL,*L2t=NULL,*L3t=NULL;
  struct litem *l1=NULL,*l2=NULL,*l3=NULL;
  int L3less=0,L3more=0,L3tot=0;
  if (verbose){ if (compare==&double_compare){ 
    printf(" %% double comparison of\n");llistprintf2(L1);printf(" %% and\n");llistprintf2(L2);}}
  if (sort_flag){ L1t=L1;L2t=L2;} else{ L1t=llistcopy(L1);L2t=llistcopy(L2);} L3t=llistmake();
  llistsort(L1t->first,L1t->last,L1t->length,compare);
  llistsort(L2t->first,L2t->last,L2t->length,compare);
  l1=L1t->first;l2=L2t->first;
  while (l1!=NULL || l2!=NULL){
    if (0){}
    else if (l1==NULL){ litemadd(L3t,l2->item); l2=l2->child;}
    else if (l2==NULL){ litemadd(L3t,l1->item); l1=l1->child;}
    else if (l1!=NULL && l2!=NULL){
      if ((*compare)(l1->item,l2->item)>0){ litemadd(L3t,l2->item); l2=l2->child;}
      else{ litemadd(L3t,l1->item); l1=l1->child;}}}
  if (verbose){ if (compare==&double_compare){ 
    printf(" %% sorted...\n");llistprintf2(L1t);llistprintf2(L2t);llistprintf2(L3t);}}
  if (L1t->length==0 && L2t->length==0){ /* do nothing */ *L1less=0;*L2more=0;*iout=NULL;}
  else if (L1t->length==0 && L2t->length>0){ /* easy */ *L1less=0;*L2more=L2t->length;*iout=(void *)L2t->first->item;}
  else if (L1t->length>0 && L2t->length==0){ /* easy */ *L1less=L1t->length;*L2more=0;*iout=(void *)L1t->last->item;}
  else if (L1t->length>0 && L2t->length>0){ 
    l1=L1t->first;l2=L2t->first;l3=L3t->first; *iout=l3->item;
    if ((*compare)((void *)l3->item,(void *)l1->item)==0){ *L1less=0;*L2more=L2->length;}
    else if ((*compare)((void *)l3->item,(void *)l2->item)==0){ *L1less=0;*L2more=L2->length;}
    L3less=*L1less;L3more=*L2more;L3tot=L3less+L3more;
    if (verbose){ if (compare==&double_compare){
      printf(" %% starting with l1->item %f l2->item %f l3->item %f\n",*(double *)l1->item,*(double *)l2->item,*(double *)l3->item);
      printf(" %% *iout %f L3less %d L3more %d L3tot %d\n",*(double *)*iout,L3less,L3more,L3tot);}}
    while (l3!=NULL){
      if (l1!=NULL && l3->item==l1->item){ 
	if (verbose){ printf(" %%  belongs to first llist \n");}
	l1=l1->child; l3=l3->child; L3less += 1;
	if (l3!=NULL && L3less+L3more > L3tot){ *iout=l3->item; *L1less=L3less;*L2more=L3more; L3tot=L3less+L3more;}}
      else if (l2!=NULL && l3->item==l2->item){ 
	if (verbose){ printf(" %%  belongs to second llist \n");}
	l2=l2->child; l3=l3->child; L3more -= 1;
	if (l3!=NULL && L3less+L3more > L3tot){ *iout=l3->item; *L1less=L3less;*L2more=L3more; L3tot=L3less+L3more;}}
      if (verbose){ if (compare==&double_compare){ 
	printf(" %% %% now l1->item %f l2->item %f l3->item %f\n",l1==NULL?0:*(double *)l1->item,l2==NULL?0:*(double *)l2->item,l3==NULL?0:*(double *)l3->item);
	printf(" %% %% *iout %f L3less %d L3more %d L3tot %d\n",*(double *)*iout,L3less,L3more,L3tot);}}}
    if (verbose){ if (compare==&double_compare){
      printf(" %% ending with *iout %f *L1less %d *L2more %d L3tot %d\n",*(double *)*iout,*L1less,*L2more,L3tot);}}}
  if (!sort_flag){ llisttfree(L1t);L1t=NULL;llisttfree(L2t);L2t=NULL;} llisttfree(L3t);L3t=NULL;
}

void discdata_threshold(double *discdata,int nbases,int maxlevel,int nrecords,char *outputname)
{
  /* given the indicators of discdata (sorted so that discdata[0 + level*nbases + nr*nbases*maxlevel] should be lowest) */
  int verbose=0;
  int nl=0,nb=0,nr=0;
  double val=0;
  struct llist *L=NULL,*Lsml1=NULL,*Lsml2=NULL,*Lbig1=NULL,*Lbig2=NULL;
  double wtol=0.0000001;
  int Lsmlless=0,Lsmlmore=0,Lbigless=0,Lbigmore=0;
  void *ismlout=NULL,*ibigout=NULL;
  double min1=0,min2=0,max1=0,max2=0;
  double *smlmarg=NULL,*smlpass=NULL,*bigmarg=NULL,*bigpass=NULL,*classifysml=NULL,*classifybig=NULL,*threshsml=NULL,*threshbig=NULL;
  double *outprint=NULL;
  if (maxlevel>1 && nbases>1){
    smlmarg = (double *) tcalloc(maxlevel*nrecords,sizeof(double));
    smlpass = (double *) tcalloc(maxlevel*nrecords,sizeof(double));
    bigmarg = (double *) tcalloc(maxlevel*nrecords,sizeof(double));
    bigpass = (double *) tcalloc(maxlevel*nrecords,sizeof(double));
    threshsml = (double *) tcalloc(maxlevel,sizeof(double));
    threshbig = (double *) tcalloc(maxlevel,sizeof(double));
    outprint = (double *) tcalloc(maxlevel,sizeof(double));
    for (nr=0;nr<nrecords;nr++){
      if (verbose){ 
	printf(" %% examining record %d\n",nr);
	raprintf(&(discdata[0 + 0*nbases + nr*nbases*maxlevel]),"double",nbases,maxlevel," %% ");}
      for (nl=0;nl<maxlevel;nl++){
	val = discdata[0+nl*nbases+nr*nbases*maxlevel];
	L = llistmake();
	for (nb=0;nb<nbases;nb++){ litemadd(L,&(discdata[nb+nl*nbases+nr*nbases*maxlevel]));}
	if (verbose){ printf(" %% level %d ",nl); llistprintf2(L);}
	llistsort(L->first,L->last,L->length,&double_compare);
	if (verbose){ printf(" %% sorted... "); llistprintf2(L);}
	min1=*(double *)L->first->item; min2=*(double *)L->first->child->item;
	max1=*(double *)L->last->item; max2=*(double *)L->last->parent->item;
	smlpass[nl+nr*maxlevel] = (fabs(val-min1)<wtol?+1:-1); smlmarg[nl+nr*maxlevel] = fabs((min2-min1)/min1);
	bigpass[nl+nr*maxlevel] = (fabs(val-max1)<wtol?+1:-1); bigmarg[nl+nr*maxlevel] = fabs((max2-max1)/max2);
	if (verbose){ printf(" %% min1 %f min2 %f max2 %f max1 %f, smlmarg %f smlpass %f, bigmarg %f bigpass %f\n",min1,min2,max2,max2,smlmarg[nl+nr*maxlevel],smlpass[nl+nr*maxlevel],bigmarg[nl+nr*maxlevel],bigpass[nl+nr*maxlevel]);}
	llisttfree(L);L=NULL;}}
    classifysml = (double *) tcalloc(nrecords*maxlevel,sizeof(double)); classifybig = (double *) tcalloc(nrecords*maxlevel,sizeof(double));
    for (nr=0;nr<nrecords;nr++){ 
      classifysml[nr+0*nrecords]=(smlpass[0+nr*maxlevel]>=0);classifybig[nr+0*nrecords]=(bigpass[0+nr*maxlevel]>=0);}
    for (nl=1;nl<maxlevel;nl++){
      if (verbose){ printf(" %% performing check for levels (%d,%d)\n",(nl-1),nl);}
      Lsml1=llistmake();Lsml2=llistmake(); Lbig1=llistmake();Lbig2=llistmake();
      for (nr=0;nr<nrecords;nr++){
	if (0){}
	else if (!classifysml[nr+(nl-1)*nrecords] && smlpass[nl+nr*maxlevel]>=0){ /* population 1 */ 
	  litemadd(Lsml1,&(smlmarg[(nl-1)+nr*maxlevel]));}
	else if (classifysml[nr+(nl-1)*nrecords] && smlpass[nl+nr*maxlevel]<0){ /* population 2 */ 
	  litemadd(Lsml2,&(smlmarg[(nl-1)+nr*maxlevel]));}
	if (0){}
	else if (!classifybig[nr+(nl-1)*nrecords] && bigpass[nl+nr*maxlevel]>=0){ /* population 1 */ 
	  litemadd(Lbig1,&(bigmarg[(nl-1)+nr*maxlevel]));}
	else if (classifybig[nr+(nl-1)*nrecords] && bigpass[nl+nr*maxlevel]<0){ /* population 2 */ 
	  litemadd(Lbig2,&(bigmarg[(nl-1)+nr*maxlevel]));}}
      if (verbose){ 
	printf(" %% sml bad to good "); llistprintf2(Lsml1);
	printf(" %% sml good to bad "); llistprintf2(Lsml2);
	printf(" %% big bad to good "); llistprintf2(Lbig1);
	printf(" %% big good to bad "); llistprintf2(Lbig2);}
      llistmedian(Lsml1,Lsml2,&double_compare,1,&ismlout,&Lsmlless,&Lsmlmore); threshsml[(nl-1)]=(ismlout==NULL?0:*(double *)ismlout);
      if (verbose){ printf(" %% sml divider %f, left %d, right %d\n",threshsml[(nl-1)],Lsmlless,Lsmlmore);}
      llistmedian(Lbig1,Lbig2,&double_compare,1,&ibigout,&Lbigless,&Lbigmore); threshbig[(nl-1)]=(ibigout==NULL?0:*(double *)ibigout);
      if (verbose){ printf(" %% big divider %f, left %d, right %d\n",threshsml[(nl-1)],Lbigless,Lbigmore);}
      for (nr=0;nr<nrecords;nr++){
	if (verbose){ printf(" %% classifysml[%d]=%f,smlpass[%d]=%f,smlmarg[%d]=%f,smlpass[%d]=%f,smlmarg[%d]=%f\n",nr+(nl-1)*nrecords,classifysml[nr+(nl-1)*nrecords],(nl-1)+nr*maxlevel,smlpass[(nl-1)+nr*maxlevel],(nl-1)+nr*maxlevel,smlmarg[(nl-1)+nr*maxlevel],nl+nr*maxlevel,smlpass[nl+nr*maxlevel],nl+nr*maxlevel,smlmarg[nl+nr*maxlevel]);}
	if (classifysml[nr+(nl-1)*nrecords]){
	  if (smlmarg[(nl-1)+nr*maxlevel]>=threshsml[(nl-1)]){ classifysml[nr+nl*nrecords] = 1;}
	  else /* if (smlmarg[(nl-1)+nr*maxlevel]<threshsml[(nl-1)]) */{ 
	    if (smlpass[nl+nr*maxlevel]>=0){ classifysml[nr+nl*nrecords] = 1;}
	    else /* if (smlpass[nl+nr*maxlevel]<0) */{ classifysml[nr+nl*nrecords] = 0;}}}
	else /* if (!classifysml[nr+(nl-1)*nrecords]) */{
	  if (smlmarg[(nl-1)+nr*maxlevel]>=threshsml[(nl-1)]){ classifysml[nr+nl*nrecords] = 0;}
	  else /* if (smlmarg[(nl-1)+nr*maxlevel]<threshsml[(nl-1)]) */{ 
	    if (smlpass[nl+nr*maxlevel]>=0){ classifysml[nr+nl*nrecords] = 1;}
	    else /* if (smlpass[nl+nr*maxlevel]<0) */{ classifysml[nr+nl*nrecords] = 0;}}}
	if (verbose){ printf(" %% classifysml[%d]=%f\n",nr+nl*nrecords,classifysml[nr+nl*nrecords]);}}
      for (nr=0;nr<nrecords;nr++){
	if (verbose){ printf(" %% classifybig[%d]=%f,bigpass[%d]=%f,bigmarg[%d]=%f,bigpass[%d]=%f,bigmarg[%d]=%f\n",nr+(nl-1)*nrecords,classifybig[nr+(nl-1)*nrecords],(nl-1)+nr*maxlevel,bigpass[(nl-1)+nr*maxlevel],(nl-1)+nr*maxlevel,bigmarg[(nl-1)+nr*maxlevel],nl+nr*maxlevel,bigpass[nl+nr*maxlevel],nl+nr*maxlevel,bigmarg[nl+nr*maxlevel]);}
	if (classifybig[nr+(nl-1)*nrecords]){
	  if (bigmarg[(nl-1)+nr*maxlevel]>=threshbig[(nl-1)]){ classifybig[nr+nl*nrecords] = 1;}
	  else /* if (bigmarg[(nl-1)+nr*maxlevel]<threshbig[(nl-1)]) */{ 
	    if (bigpass[nl+nr*maxlevel]>=0){ classifybig[nr+nl*nrecords] = 1;}
	    else /* if (bigpass[nl+nr*maxlevel]<0) */{ classifybig[nr+nl*nrecords] = 0;}}}
	else /* if (!classifybig[nr+(nl-1)*nrecords]) */{
	  if (bigmarg[(nl-1)+nr*maxlevel]>=threshbig[(nl-1)]){ classifybig[nr+nl*nrecords] = 0;}
	  else /* if (bigmarg[(nl-1)+nr*maxlevel]<threshbig[(nl-1)]) */{ 
	    if (bigpass[nl+nr*maxlevel]>=0){ classifybig[nr+nl*nrecords] = 1;}
	    else /* if (bigpass[nl+nr*maxlevel]<0) */{ classifybig[nr+nl*nrecords] = 0;}}}
	if (verbose){ printf(" %% classifybig[%d]=%f\n",nr+nl*nrecords,classifybig[nr+nl*nrecords]);}}
      llisttfree(Lsml1);Lsml1=NULL; llisttfree(Lsml2);Lsml2=NULL;
      llisttfree(Lbig1);Lbig1=NULL; llisttfree(Lbig2);Lbig2=NULL;}
    for (nl=0;nl<maxlevel;nl++){
      stats("double",&(classifysml[0+nl*nrecords]),nrecords,NULL,NULL,&max1,NULL);
      stats("double",&(classifybig[0+nl*nrecords]),nrecords,NULL,NULL,&max2,NULL);
      if (verbose){ printf(" %% level %d, smldisc %f, bigdisc %f\n",nl,max1,max2);}
      outprint[nl]=maximum(max1,max2);}
    ra2jpg2(outprint,"double",1,maxlevel,0,1.0,0.0,outputname);
    tfree(smlmarg);tfree(smlpass);tfree(classifysml);tfree(threshsml);
    tfree(bigmarg);tfree(bigpass);tfree(classifybig);tfree(threshbig);
    tfree(outprint);}
}

void ptreeobsdist(int nbases,struct ptree *pR,struct pnode *parentR,struct llitem *lR,int nb,struct ptree *p1,struct pnode *parent1,struct llitem *l1,int maxlevel,int level)
{
  /* creates llistra (nbases elements) at each node, corresponding to observation distribution of each event-chain */
  int verbose=1,nr=0,nb2=0;
  char lvlchar[64];
  int firstentry=0;
  struct pnode *pnR=NULL,*pn1=NULL;
  struct llitem *l2=NULL;
  struct llist **Lra=NULL;
  double wtol=0.0000001;
  double *temp=NULL;
  if (verbose){ sprintf(lvlchar," %%"); for (nr=0;nr<level;nr++){sprintf(lvlchar,"%s %%",lvlchar);}}
  if (verbose>1){ 
    printf(" %s [entering ptreeobsdist]\n",lvlchar);
    printf(" %s pR %d parentR %d->%d lR %d \n",lvlchar,(int)pR,(int)parentR,parentR==NULL ? 0 : (int)(parentR->childllitem),(int)lR);
    printf(" %s p1 %d parent1 %d->%d l1 %d\n",lvlchar,(int)p1,(int)parent1,parent1==NULL ? 0 : (int)(parent1->childllitem),(int)l1);
    printf(" %s nbases %d maxlevel %d level %d\n",lvlchar,nbases,maxlevel,level);}
  if (parentR!=NULL){ assert(parentR->childllitem==lR);}
  if (parent1!=NULL && l1->parent==NULL){ /* first descent into parent1->childllitem */ 
    if (verbose>1){ printf(" %s first descent into parent1->childllitem, setting firstentry to 1\n",lvlchar);}
    assert(parent1->childllitem==l1); firstentry=1;}
  if (l1->kidl!=NULL){ 
    if (verbose>1){ printf(" %s l1->kidl found, moving left to new llitem l1->kidl=%d\n",lvlchar,(int)(l1->kidl));}
    ptreeobsdist(nbases,pR,parentR,lR,nb,p1,parent1,l1->kidl,maxlevel,level);}
  if (l1->kidr!=NULL){ 
    if (verbose>1){ printf(" %s l1->kidr found, moving right to new llitem l1->kidr=%d\n",lvlchar,(int)(l1->kidr));}
    ptreeobsdist(nbases,pR,parentR,lR,nb,p1,parent1,l1->kidr,maxlevel,level);}
  if (l1->item!=NULL){ 
    if (verbose>1){ printf(" %s l1->item %d exists\n",lvlchar,(int)(struct pnode *)(l1->item));}
    pn1 = (struct pnode *)l1->item; 
    if (maxlevel==-1 || level<maxlevel){
      if (verbose>1){ printf(" %s level %d less than maxlevel %d, continuing...\n",lvlchar,level,maxlevel);}
      if (lR!=NULL && (l2=llitemaddorfind(0,lR,pR->regionra[pn1->region->label],&region2pnode_compare_label))!=NULL){
	pnR = (struct pnode *)l2->item;
	if (verbose>1){ printf(" %s pn1 label %d matches pnR label %d, descending...\n",lvlchar,pn1->region->label,pnR->region->label);}
	ptreeobsdist(nbases,pR,pnR,pnR->childllitem,nb,p1,pn1,pn1->childllitem,maxlevel,level+1);}
      else{
	pnR = NULL; 
	if (verbose>1){ printf(" %s pn1 label %d not found,terminating\n",lvlchar,pn1->region->label);}}}}
  if (firstentry){
    if (parent1!=NULL){ 
      assert(parentR!=NULL);
      if (parentR->temp==NULL){ 
	Lra = (struct llist **) tcalloc(nbases,sizeof(struct llist *)); parentR->temp = (void *) Lra;
	for (nb2=0;nb2<nbases;nb2++){ Lra[nb2] = llistmake();}}
      Lra = (struct llist **) parentR->temp;
      if (fabs(parent1->weight)>wtol){ 
	if (verbose>2){ printf(" %s starting with llist:",lvlchar); llistprintf2(Lra[nb]);}
	temp=(double *)tcalloc(1,sizeof(double));*temp=parent1->weight;litemadd(Lra[nb],temp);
	if (verbose>2){ printf(" %s adding %f to llist %d at pR %d:",lvlchar,parent1->weight,nb,(int)parentR); llistprintf2(Lra[nb]);}}}}
}

void obsdistdisc(int nbases,int nrecords,struct llist **Lra,double *discriminability,double *reliability)
{
  /* given llistra (nbases elements) of observation distributions, calculates discriminability (assuming good statistics for now) */
  int verbose=0;
  int nb=0,nr=0,nv=0;
  double *maxra=NULL,max=0,max2=0,mean=0;
  struct hist **histra=NULL;
  struct litem *l=NULL;
  double nright=0,nwrong=0,nweak=0;
  double wtol=0.0000001;
  int nmax=0;
  maxra = (double *) tcalloc(nbases,sizeof(double));
  for (nb=0;nb<nbases;nb++){ lliststats(Lra[nb],&(maxra[nb]),NULL,NULL,NULL);}
  stats("double",maxra,nbases,&max,NULL,NULL,NULL); tfree(maxra);maxra=NULL;
  histra = (struct hist **) tcalloc(nbases,sizeof(struct hist *));
  for (nb=0;nb<nbases;nb++){ 
    if (verbose){ printf(" %% dealing with nb %d llist:\n",nb); llistprintf2(Lra[nb]);}
    histra[nb]=histmake((int)max,max,0); 
    l=Lra[nb]->first; nr=0; while(l!=NULL){ histadd(histra[nb],*(double *)l->item,1); nr+=1; l=l->child;}
    histadd(histra[nb],0,nrecords-nr);
    if (verbose){ printf(" %% now have hist:\n"); histprintf(histra[nb]," %% ");}}
  nright=0;nwrong=0;nweak=0;
  for (nv=0;nv<(int)max;nv++){
    maxra=(double *)tcalloc(nbases,sizeof(double));
    for (nb=0;nb<nbases;nb++){ if (nv>=histra[nb]->nbins){ maxra[nb]=0;} else{ maxra[nb]=histra[nb]->data[nv];}}
    if (verbose){ printf(" %% comparing bin %d: \n",nv); raprintf(maxra,"double",1,nbases," %% ");}
    stats("double",maxra,nbases,&max2,NULL,&mean,NULL);
    if ((nbases*mean)>1.5){
      nmax=0; for (nb=0;nb<nbases;nb++){ if (fabs(maxra[nb]-max2)<wtol){ nmax +=1;}}
      for (nb=0;nb<nbases;nb++){ 
	if (fabs(maxra[nb]-max2)<wtol){ nright += 1.0/(double)nmax*maxra[nb]; nwrong += (1-1.0/(double)nmax)*maxra[nb];}
	else{ nwrong += maxra[nb];}}
      if (verbose){ printf(" %% non-singleton, nmax %d nright %f nwrong %f\n",nmax,nright,nwrong);}}
    else if ((nbases*mean)<=1.5 && fabs(nbases*mean)>wtol){ 
      nweak += 1; if (verbose){ printf(" %% singleton, ignoring, nweak %f\n",nweak);}}
    else if (fabs(nbases*mean)<wtol){ if (verbose){ printf(" %% both empty, ignoring\n");}}
    tfree(maxra);maxra=NULL;}
  if (nright+nwrong>0){ *discriminability = nright/(nright+nwrong);} else{ *discriminability = 1.0/(double)nbases;}
  *reliability = 1-nweak/maximum(1,nright+nwrong+nweak);
  if (verbose){ printf(" %% now discriminability %f reliability %f\n",*discriminability,*reliability);}
  for (nb=0;nb<nbases;nb++){ histtfree(histra[nb]);histra[nb]=NULL;} tfree(histra);histra=NULL;
}

void ptreeobsdistend(int nbases,int nrecords,struct ptree *pR,struct pnode *parentR,struct llitem *lR,int level)
{
  /* checks discriminability of each event-chain, and then tfrees llistra (nbases elements) at each node */
  int verbose=0,nr=0,nb=0;
  char lvlchar[64];
  int firstentry=0;
  struct pnode *pnR=NULL;
  struct llist **Lra=NULL;
  if (verbose){ sprintf(lvlchar," %%"); for (nr=0;nr<level;nr++){sprintf(lvlchar,"%s %%",lvlchar);}}
  if (verbose>1){ 
    printf(" %s [entering ptreeobsdistend]\n",lvlchar);
    printf(" %s pR %d parentR %d->%d lR %d \n",lvlchar,(int)pR,(int)parentR,parentR==NULL ? 0 : (int)(parentR->childllitem),(int)lR);
    printf(" %s nbases %d nrecords %d level %d\n",lvlchar,nbases,nrecords,level);}
  if (parentR!=NULL && lR->parent==NULL){ /* first descent into parentR->childllitem */ 
    if (verbose>1){ printf(" %s first descent into parentR->childllitem, setting firstentry to 1\n",lvlchar);}
    assert(parentR->childllitem==lR); firstentry=1;}
  if (lR->kidl!=NULL){ 
    if (verbose>1){ printf(" %s lR->kidl found, moving left to new llitem lR->kidl=%d\n",lvlchar,(int)(lR->kidl));}
    ptreeobsdistend(nbases,nrecords,pR,parentR,lR->kidl,level);}
  if (lR->kidr!=NULL){ 
    if (verbose>1){ printf(" %s lR->kidr found, moving right to new llitem lR->kidr=%d\n",lvlchar,(int)(lR->kidr));}
    ptreeobsdistend(nbases,nrecords,pR,parentR,lR->kidr,level);}
  if (lR->item!=NULL){ 
    if (verbose>1){ printf(" %s lR->item %d exists\n",lvlchar,(int)(struct pnode *)(lR->item));}
    pnR = (struct pnode *)lR->item; 
    ptreeobsdistend(nbases,nrecords,pR,pnR,pnR->childllitem,level+1);}
  if (firstentry){
    if (parentR!=NULL){ 
      if (parentR->temp!=NULL){ 
	Lra = (struct llist **) parentR->temp;
	obsdistdisc(nbases,nrecords,Lra,&(parentR->weight),&(parentR->relevance));
	for (nb=0;nb<nbases;nb++){ llisttfree2(Lra[nb]); Lra[nb]=NULL;} tfree(Lra);
	if (verbose>1){ printf(" %s freeing Lra at pR\n",lvlchar);}}}
    else{ if (verbose>1){ printf(" %s Lra not found, weight must be 0\n",lvlchar);} parentR->weight=0;}}
}

void ptree_trialaverage_helper(int nbases,int maxlevel,int tclump,int verbose,int X_flag,int x_flag,int s_flag,int d_flag,char **filename_base,char *outputname)
{
  /* assumes that outputname does NOT start with "./" */
  int old_region_type;
  int nb=0,nb2=0,nmin=0; 
  double ntemp=0;
  char **filename=NULL,allstring[16];
  char *gs2=GLOBAL_STRING_2;
  char outputfilename[256],outputfilename2[256];
  char outputfilename_disc_x2p[256];
  char outputfilename_disc_depthent[256];
  char outputfilename_disc_breadthent[256];
  char outputfilename_obsdist[256];
  char plotfilename2[256],plotfilename_x2p[256],plotfilename_x2pper[256];
  char plotfilename_depthent[256],plotfilename_depthentper[256];
  char plotfilename_breadthent[256],plotfilename_breadthentper[256];
  FILE **fp=NULL,*op=NULL,*op2=NULL;
  struct ptree **p0=NULL,**p1=NULL,*pR=NULL,*ptmp=NULL;
  struct pnode **pn0=NULL;
  struct llitem **l0pos=NULL,**l0pre=NULL;
  double *nright_x2=NULL,*nwrong_x2=NULL;
  double *nright_x2p=NULL,*nwrong_x2p=NULL;
  double *nright_fabs=NULL,*nwrong_fabs=NULL;
  double *nright_fabs2=NULL,*nwrong_fabs2=NULL;
  double *nright_depthent=NULL,*nwrong_depthent=NULL;
  double *nright_breadthent=NULL,*nwrong_breadthent=NULL;
  double *temp_x2=NULL,*temp_x2p=NULL,*temp_fabs=NULL,*temp_fabs2=NULL,*temp_depthent=NULL,*temp_breadthent=NULL;
  double *discra=NULL;
  double *dist_x2p=NULL;
  double *nright_x2p_perstimulus=NULL,*nwrong_x2p_perstimulus=NULL;
  double *nright_depthent_perstimulus=NULL,*nwrong_depthent_perstimulus=NULL;
  double *nright_breadthent_perstimulus=NULL,*nwrong_breadthent_perstimulus=NULL;
  int bitbybit=0,continue_flag=0,level=0,tcounter=0;
  double wmax=0,wmin=0,rmax=0,rmin=0;
  double *percentages=NULL,*sbias=NULL,maxbias=0,biastol=0.00001;
  int nx,nX=0,nR=0;
  double *nbias=NULL;
  struct hist *temphist=NULL;
  int nelements=0;
  double *referencenormalizers=NULL;
  double *observednormalizers=NULL;
  double *discdata_x2p=NULL;
  double *discdata_depthent=NULL;
  double *discdata_breadthent=NULL;
  if (verbose>1){ 
    printf(" %% [entering ptree_trialaverage_helper]\n");
    printf(" %% nbases %d maxlevel %d tclump %d verbose %d X_flag %d x_flag %d s_flag %d d_flag %d\n",nbases,maxlevel,tclump,verbose,X_flag,x_flag,s_flag,d_flag);
    for (nb=0;nb<nbases;nb++){ printf(" %% filename_base[%d]=%s\n",nb,filename_base[nb]);}
    printf(" %% outputname=%s\n",outputname);}
  if (outputname==NULL){ 
    sprintf(outputfilename,"./outprint_T%d_s%d_d%d_%srecord.m",tclump,s_flag,d_flag,gs2); 
    sprintf(outputfilename2,"./outdist_T%d_s%d_d%d_%srecord.m",tclump,s_flag,d_flag,gs2); 
    sprintf(outputfilename_disc_x2p,"./outdiscdata_x2p_T%d_s%d_d%d_%srecord",tclump,s_flag,d_flag,gs2); 
    sprintf(outputfilename_disc_depthent,"./outdiscdata_depthent_T%d_s%d_d%d_%srecord",tclump,s_flag,d_flag,gs2); 
    sprintf(outputfilename_disc_breadthent,"./outdiscdata_breadthent_T%d_s%d_d%d_%srecord",tclump,s_flag,d_flag,gs2); 
    sprintf(outputfilename_obsdist,"./outobsdist_s%d_d%d_%srecord",s_flag,d_flag,gs2);
    sprintf(plotfilename2,"./outdist_s%d_d%d_%srecord",s_flag,d_flag,gs2);
    sprintf(plotfilename_x2p,"./outprint_x2p_s%d_d%d_%srecord",s_flag,d_flag,gs2);
    sprintf(plotfilename_x2pper,"./outprintper_x2p_s%d_d%d_%srecord",s_flag,d_flag,gs2);
    sprintf(plotfilename_depthent,"./outprint_depthent_s%d_d%d_%srecord",s_flag,d_flag,gs2);
    sprintf(plotfilename_depthentper,"./outprintper_depthent_s%d_d%d_%srecord",s_flag,d_flag,gs2);
    sprintf(plotfilename_breadthent,"./outprint_breadthent_s%d_d%d_%srecord",s_flag,d_flag,gs2);
    sprintf(plotfilename_breadthentper,"./outprintper_breadthent_s%d_d%d_%srecord",s_flag,d_flag,gs2);}
  else{ 
    sprintf(outputfilename,"./outprint%s_T%d_s%d_d%d_%srecord.m",outputname,tclump,s_flag,d_flag,gs2); 
    sprintf(outputfilename2,"./outdist%s_T%d_s%d_d%d_%srecord.m",outputname,tclump,s_flag,d_flag,gs2); 
    sprintf(outputfilename_disc_x2p,"./outdiscdata_x2p_%s_T%d_s%d_d%d_%srecord",outputname,tclump,s_flag,d_flag,gs2); 
    sprintf(outputfilename_disc_depthent,"./outdiscdata_depthent_%s_T%d_s%d_d%d_%srecord",outputname,tclump,s_flag,d_flag,gs2); 
    sprintf(outputfilename_disc_breadthent,"./outdiscdata_breadthent_%s_T%d_s%d_d%d_%srecord",outputname,tclump,s_flag,d_flag,gs2); 
    sprintf(outputfilename_obsdist,"./outobsdist_%s_s%d_d%d_%srecord",outputname,s_flag,d_flag,gs2);
    sprintf(plotfilename2,"./outdist%s_s%d_d%d_%srecord",outputname,s_flag,d_flag,gs2);
    sprintf(plotfilename_x2p,"./outprint_x2p_%s_s%d_d%d_%srecord",outputname,s_flag,d_flag,gs2);
    sprintf(plotfilename_x2pper,"./outprintper_x2p_%s_s%d_d%d_%srecord",outputname,s_flag,d_flag,gs2);
    sprintf(plotfilename_depthent,"./outprint_depthent_%s_s%d_d%d_%srecord",outputname,s_flag,d_flag,gs2);
    sprintf(plotfilename_depthentper,"./outprintper_depthent_%s_s%d_d%d_%srecord",outputname,s_flag,d_flag,gs2);
    sprintf(plotfilename_breadthent,"./outprint_breadthent_%s_s%d_d%d_%srecord",outputname,s_flag,d_flag,gs2);
    sprintf(plotfilename_breadthentper,"./outprintper_breadthent_%s_s%d_d%d_%srecord",outputname,s_flag,d_flag,gs2);}
  if ((op=fopen(outputfilename,"w"))==NULL){ if (verbose>1){ printf(" %% writing to stdout\n");} op=stdout;}
  if ((op2=fopen(outputfilename2,"w"))==NULL){ if (verbose>1){ printf(" %% writing to stdout\n");} op2=stdout;}
  sprintf(allstring,"_s%d_all",s_flag);
  nright_x2 = (double *) tcalloc(maxlevel,sizeof(double));
  nwrong_x2 = (double *) tcalloc(maxlevel,sizeof(double));
  nright_x2p = (double *) tcalloc(maxlevel,sizeof(double));
  nwrong_x2p = (double *) tcalloc(maxlevel,sizeof(double));
  nright_fabs = (double *) tcalloc(maxlevel,sizeof(double));
  nwrong_fabs = (double *) tcalloc(maxlevel,sizeof(double));
  nright_fabs2 = (double *) tcalloc(maxlevel,sizeof(double));
  nwrong_fabs2 = (double *) tcalloc(maxlevel,sizeof(double));
  nright_depthent = (double *) tcalloc(maxlevel,sizeof(double));
  nwrong_depthent = (double *) tcalloc(maxlevel,sizeof(double));
  nright_breadthent = (double *) tcalloc(maxlevel,sizeof(double));
  nwrong_breadthent = (double *) tcalloc(maxlevel,sizeof(double));
  discra = (double *) tcalloc(maxlevel*nbases,sizeof(double));
  filename = (char **) tcalloc(nbases,sizeof(char *));
  for (nb=0;nb<nbases;nb++){ 
    filename[nb] = (char *) tcalloc(256,sizeof(char));}
  fp = (FILE **) tcalloc(nbases,sizeof(FILE *));
  p0 = (struct ptree **) tcalloc(nbases,sizeof(struct ptree *));
  p1 = (struct ptree **) tcalloc(nbases,sizeof(struct ptree *));
  temp_x2 = (double *) tcalloc(nbases,sizeof(double));
  temp_x2p = (double *) tcalloc(nbases,sizeof(double));
  temp_fabs = (double *) tcalloc(nbases,sizeof(double));
  temp_fabs2 = (double *) tcalloc(nbases,sizeof(double));
  temp_depthent = (double *) tcalloc(nbases,sizeof(double));
  temp_breadthent = (double *) tcalloc(nbases,sizeof(double));
  dist_x2p = (double *) tcalloc(nbases*nbases*maxlevel,sizeof(double)); /* deduced_base + actual_base*nbases + level*nbases*nbases */
  nright_x2p_perstimulus = (double *) tcalloc(maxlevel*nbases,sizeof(double));
  nwrong_x2p_perstimulus = (double *) tcalloc(maxlevel*nbases,sizeof(double));
  nright_depthent_perstimulus = (double *) tcalloc(maxlevel*nbases,sizeof(double));
  nwrong_depthent_perstimulus = (double *) tcalloc(maxlevel*nbases,sizeof(double));
  nright_breadthent_perstimulus = (double *) tcalloc(maxlevel*nbases,sizeof(double));
  nwrong_breadthent_perstimulus = (double *) tcalloc(maxlevel*nbases,sizeof(double));
  referencenormalizers = (double *) tcalloc(maxlevel*nbases,sizeof(double));
  observednormalizers = (double *) tcalloc(maxlevel,sizeof(double));
  old_region_type=GLOBAL_PTREE_REGION_TYPE; GLOBAL_PTREE_REGION_TYPE=0;
  if (verbose>1){ printf(" %% setting GLOBAL_PTREE_REGION_TYPE to %d\n",GLOBAL_PTREE_REGION_TYPE);}
  for (nb=0;nb<nbases;nb++){
    if (verbose){ printf(" %% checking for %s%s... ",filename_base[nb],allstring);}
    sprintf(filename[nb],"./%s%s",filename_base[nb],allstring);
    if (checktofind(filename[nb])){ if (verbose){ printf("found\n");} p0[nb] = ptreadback(filename[nb],0);}
    else /* if not found */{
      if (verbose){ printf("not found, making...\n");}
      bitbybit=1;
      sprintf(filename[nb],"%s_%d",filename_base[nb],bitbybit);
      p0[nb] = ptreadback(filename[nb],0); 
      if (s_flag==1){ if (verbose>1){ printf("shizuffling %d\n",s_flag);} pnodeshizuffle_starter(NULL,p0[nb]->postree,p0[nb]); pnodeshizuffle_starter(NULL,p0[nb]->pretree,p0[nb]);}
      bitbybit=2; continue_flag=1;
      do{
	sprintf(filename[nb],"./%s_%d",filename_base[nb],bitbybit);
	if (verbose>1){ printf(" %% trying %s... ",filename[nb]);}
	if (checktofind(filename[nb])){ 
	  if (verbose>1){ printf("found\n");}
	  continue_flag=1; 
	  p1[nb] = ptreadback(filename[nb],0);
	  if (s_flag==1){ if (verbose>1){ printf("shizuffling %d\n",s_flag);} pnodeshizuffle_starter(NULL,p1[nb]->postree,p1[nb]); pnodeshizuffle_starter(NULL,p1[nb]->pretree,p1[nb]);}
	  ptreeplusequals_starter(p0[nb],p1[nb]);
	  ptreetfree(p1[nb]); p1[nb]=NULL;}
	else /* if not found */{ if (verbose>1){ printf("not found, stopping\n");} continue_flag=0;}
	bitbybit+=1;}
      while (continue_flag);
      sprintf(filename[nb],"./%s%s",filename_base[nb],allstring);
      if (verbose){ printf(" %% dumping %s\n",filename[nb]);}
      ptreerate(p0[nb]); ptreedump_starter(p0[nb],filename[nb],2,1,0,0,+1,-1,0);}
    if (verbose){ printf(" %% p0[%d]->total_time %f\n",nb,p0[nb]->total_time);}
    for (level=0;level<maxlevel;level++){ nelements=0; pnodestats_starter(NULL,p0[nb]->postree,level+1,level+1,0,&nelements,NULL,NULL,&(referencenormalizers[level+nb*maxlevel]),NULL,NULL,NULL,NULL,NULL); referencenormalizers[level+nb*maxlevel]*=nelements;}
    if (verbose>1){ for (level=0;level<maxlevel;level++){ printf(" %% p0[%d] level %d normalizer = %f\n",nb,level,referencenormalizers[level+nb*maxlevel]);}}}  
  if (X_flag>=1){
    sprintf(filename[0],"./%s_%d",filename_base[0],1);
    if (verbose>1){ printf(" %% reading time from %s... ",filename[0]);}
    if (checktofind(filename[0])){ p1[0] = ptreadback(filename[0],0);} else{ printf(" %% warning, %s not found\n",filename[0]);}
    nR=0; for (nb=0;nb<nbases;nb++){ nR+=ceil(p0[nb]->total_time);} nR=ceil(2*nR/(double)(tclump*p1[0]->total_time));
    ptreetfree(p1[0]);p1[0]=NULL;
    if (verbose>1){ printf(" %% %% %% allocating discdata size %d with %d records\n",nbases*maxlevel*nR,nR);}
    discdata_x2p = (double *) tcalloc((nbases*maxlevel)*nR,sizeof(double));
    discdata_depthent = (double *) tcalloc((nbases*maxlevel)*nR,sizeof(double));
    discdata_breadthent = (double *) tcalloc((nbases*maxlevel)*nR,sizeof(double));
    nR=0;
    for (nb=0;nb<nbases;nb++){
      bitbybit=1; continue_flag=1;
      do{
	sprintf(filename[nb],"./%s_%d",filename_base[nb],bitbybit);
	if (verbose>1){ printf(" %% checking %s... ",filename[nb]);}
	if (checktofind(filename[nb])){ 
	  if (verbose>1){ printf("found, initiating tcounter\n");}
	  p1[nb] = ptreadback(filename[nb],0);
	  if (s_flag==1){ if (verbose>1){ printf("shizuffling %d\n",s_flag);} pnodeshizuffle_starter(NULL,p1[nb]->postree,p1[nb]); pnodeshizuffle_starter(NULL,p1[nb]->pretree,p1[nb]);}
	  bitbybit += 1; continue_flag=1;}
	else /* if not found */{ if (verbose>1){ printf("not found, stopping\n");} continue_flag=0;}
	tcounter=1;
	while (tcounter<tclump && continue_flag){
	  sprintf(filename[nb],"./%s_%d",filename_base[nb],bitbybit);
	  if (verbose>1){ printf(" %% checking %s... ",filename[nb]);}
	  if (checktofind(filename[nb])){ 
	    if (verbose>1){ printf("found with tcounter %d\n",tcounter);}
	    bitbybit += 1; continue_flag=1;
	    ptmp = ptreadback(filename[nb],0);
	    if (s_flag==1){ if (verbose>1){ printf("shizuffling %d\n",s_flag);} pnodeshizuffle_starter(NULL,ptmp->postree,ptmp); pnodeshizuffle_starter(NULL,ptmp->pretree,ptmp);}
	    ptreeplusequals_starter(p1[nb],ptmp);
	    ptreetfree(ptmp);ptmp=NULL;}
	  else /* if not found */{ if (verbose>1){ printf("not found, stopping\n");} continue_flag=0;}
	  tcounter += 1;}
	if (continue_flag){
	  if (verbose>1){ printf(" %% accumulated %d-%d total time %0.2f\n",bitbybit-1-tcounter+1,bitbybit-1,p1[nb]->total_time);}
	  nR += 1;
	  for (level=0;level<maxlevel;level++){ nelements=0; pnodestats_starter(NULL,p1[nb]->postree,level+1,level+1,0,&nelements,NULL,NULL,&(observednormalizers[level]),NULL,NULL,NULL,NULL,NULL); observednormalizers[level]*=nelements;}
	  if (verbose>2){ raprintf(observednormalizers,"double",1,maxlevel," %% %% observednormalizers ");}
	  for (level=1;level<=maxlevel;level++){
	    for (nb2=0;nb2<nbases;nb2++){
	      if (verbose>2){ printf(" %% %% %% nb %d level %d nb2 %d\n",nb,level,nb2);}
	      if (verbose>2){ printf(" %% %% %% calling ptreex2\n");}
	      temp_x2[nb2]=0; temp_x2p[nb2]=0; temp_fabs[nb2]=0; temp_fabs2[nb2]=0; temp_depthent[nb2]=0; temp_breadthent[nb2]=0;
	      ptreex2(p0[nb2],NULL,p0[nb2]->postree,p1[nb],NULL,p1[nb]->postree,nb==nb2 ? d_flag : 1,level,0,&(temp_x2[nb2]),&(temp_x2p[nb2]),&(temp_fabs[nb2]),&(temp_fabs2[nb2]));
	      if (verbose>2){ printf(" %% %% %% calling ptreerelent_depth\n");}
	      ptreerelent_depth(p0[nb2],NULL,p0[nb2]->postree,p1[nb],NULL,p1[nb]->postree,nb==nb2 ? d_flag : 1,level,0,&(temp_depthent[nb2]));
	      if (verbose>2){ printf(" %% %% %% calling ptreerelent_breadth\n");}
	      ptreerelent_breadth(p0[nb2],NULL,p0[nb2]->postree,p1[nb],NULL,p1[nb]->postree,nb==nb2 ? d_flag : 1,level,level,0,&(temp_breadthent[nb2]),&(referencenormalizers[0+nb2*maxlevel]),observednormalizers);}
	    if (verbose>2){ printf(" %% %% discdata tab %d+%d*nbases+%d*nbases*maxlevel\n",periodize(nb2-nb,0,nbases),level-1,nR-1);}
	    for (nb2=0;nb2<nbases;nb2++){
	      discdata_x2p[periodize(nb2-nb,0,nbases) + (level-1)*nbases + (nR-1)*nbases*maxlevel] = temp_x2p[nb2];
	      discdata_depthent[periodize(nb2-nb,0,nbases) + (level-1)*nbases + (nR-1)*nbases*maxlevel] = temp_depthent[nb2];
	      discdata_breadthent[periodize(nb2-nb,0,nbases) + (level-1)*nbases + (nR-1)*nbases*maxlevel] = temp_breadthent[nb2];}
	    ptreex2_disthist(&(dist_x2p[0 + nb*nbases + (level-1)*nbases*nbases]),temp_x2p,nbases);
	    if (verbose>1){
	      printf(" %% %% level %d, nb %d,temp_x2p[0]=%f, temp_x2p[1]=%f\n",level,nb,temp_x2p[0],temp_x2p[1]);
	      printf(" %% %% level %d, nb %d,temp_breadthent[0]=%f, temp_breadthent[1]=%f\n",level,nb,temp_breadthent[0],temp_breadthent[1]);
	      printf(" %% %% level %d, nb %d, temp_depthent[0]=%f, temp_depthent[1]=%f\n",level,nb,temp_depthent[0],temp_depthent[1]);}
	    ntemp = ptreex2_match(nbases,nb,temp_x2,NULL); nright_x2[level-1] += ntemp; nwrong_x2[level-1] += 1-ntemp;
	    ntemp = ptreex2_match(nbases,nb,temp_x2p,NULL); nright_x2p[level-1] += ntemp; nwrong_x2p[level-1] += 1-ntemp;
	    nright_x2p_perstimulus[(level-1)+nb*maxlevel] += ntemp; nwrong_x2p_perstimulus[(level-1)+nb*maxlevel] += 1-ntemp;
	    ntemp = ptreex2_match(nbases,nb,temp_fabs,NULL); nright_fabs[level-1] += ntemp; nwrong_fabs[level-1] += 1-ntemp;
	    ntemp = ptreex2_match(nbases,nb,temp_fabs2,NULL); nright_fabs2[level-1] += ntemp; nwrong_fabs2[level-1] += 1-ntemp;
	    ntemp = ptreex2_match(nbases,nb,temp_depthent,NULL); nright_depthent[level-1] += ntemp; nwrong_depthent[level-1] += 1-ntemp;
	    nright_depthent_perstimulus[(level-1)+nb*maxlevel] += ntemp; nwrong_depthent_perstimulus[(level-1)+nb*maxlevel] += 1-ntemp;
	    ntemp = ptreex2_match(nbases,nb,temp_breadthent,NULL); nright_breadthent[level-1] += ntemp; nwrong_breadthent[level-1] += 1-ntemp;
	    nright_breadthent_perstimulus[(level-1)+nb*maxlevel] += ntemp; nwrong_breadthent_perstimulus[(level-1)+nb*maxlevel] += 1-ntemp;
	    if (verbose>2){
	      printf(" %% ptreex2 for %s->postree nl %d got ",filename[nb],level);
	      for (nb2=0;nb2<nbases;nb2++){
		printf("x2 %f x2p %f fabs %f fabs2 %f ",temp_x2[nb2],temp_x2p[nb2],temp_fabs[nb2],temp_fabs2[nb2]);}
	      printf("correct index %d",nb);
	      printf("\n");}}}
	if (verbose>2){ printf(" %% %% freeing p1\n");} 
	if (p1[nb]!=NULL){ ptreetfree(p1[nb]);p1[nb]=NULL;} 
	if (verbose>2){ printf(" %% %% freed p1\n");}}
      while (continue_flag);}
    for (level=0;level<maxlevel;level++){ discra[level] = nright_x2p[level]/(nright_x2p[level]+nwrong_x2p[level]);}
    ra2jpg2(discra,"double",1,maxlevel,0,1.0,0.0,plotfilename_x2p);
    for (level=0;level<maxlevel;level++){ discra[level] = nright_depthent[level]/(nright_depthent[level]+nwrong_depthent[level]);}
    ra2jpg2(discra,"double",1,maxlevel,0,1.0,0.0,plotfilename_depthent);
    for (level=0;level<maxlevel;level++){ discra[level] = nright_breadthent[level]/(nright_breadthent[level]+nwrong_breadthent[level]);}
    ra2jpg2(discra,"double",1,maxlevel,0,1.0,0.0,plotfilename_breadthent);
    for (level=0;level<maxlevel;level++){ discra[level] = multithresh(&(dist_x2p[0 + 0*nbases + level*nbases*nbases]),nbases,0.95);}
    ra2jpg2(discra,"double",1,maxlevel,0,100,0.0,plotfilename2);
    for (nb=0;nb<nbases;nb++){ for (level=0;level<maxlevel;level++){ discra[level+nb*maxlevel]=nright_x2p_perstimulus[level+nb*maxlevel]/(nright_x2p_perstimulus[level+nb*maxlevel]+nwrong_x2p_perstimulus[level+nb*maxlevel]);}} ra2jpg2(discra,"double",maxlevel,nbases,1,1.0,0.0,plotfilename_x2pper);
    for (nb=0;nb<nbases;nb++){ for (level=0;level<maxlevel;level++){ discra[level+nb*maxlevel]=nright_depthent_perstimulus[level+nb*maxlevel]/(nright_depthent_perstimulus[level+nb*maxlevel]+nwrong_depthent_perstimulus[level+nb*maxlevel]);}} ra2jpg2(discra,"double",maxlevel,nbases,1,1.0,0.0,plotfilename_depthentper);
    for (nb=0;nb<nbases;nb++){ for (level=0;level<maxlevel;level++){ discra[level+nb*maxlevel]=nright_breadthent_perstimulus[level+nb*maxlevel]/(nright_breadthent_perstimulus[level+nb*maxlevel]+nwrong_breadthent_perstimulus[level+nb*maxlevel]);}} ra2jpg2(discra,"double",maxlevel,nbases,1,1.0,0.0,plotfilename_breadthentper);
    if (verbose){ 
      for (level=0;level<maxlevel;level++){ 
	printf(" %% level %d (%0.1f,%0.1f) or (%0.1f,%0.1f) or (%0.1f,%0.1f) or (%0.1f,%0.1f) or (%0.1f,%0.1f) or (%0.1f,%0.1f)\n",level+1,nright_x2[level],nwrong_x2[level],nright_x2p[level],nwrong_x2p[level],nright_fabs[level],nwrong_fabs[level],nright_fabs2[level],nwrong_fabs2[level],nright_depthent[level],nwrong_depthent[level],nright_breadthent[level],nwrong_breadthent[level]);}
      for (level=0;level<maxlevel;level++){ 
	fprintf(op," x2_%s_STLR(%d,%d,%d)=%0.1f; x2_%s_STLW(%d,%d,%d)=%0.1f;\n",gs2,s_flag+1,tclump,level+1,nright_x2p[level],gs2,s_flag+1,tclump,level+1,nwrong_x2p[level]);}
      fprintf(op," x2_%s_STL=x2_%s_STLR./(x2_%s_STLR + x2_%s_STLW);\n",gs2,gs2,gs2,gs2);
      fprintf(op," plot(permute(x2_%s_STL(s_to_plot,:,:),[3 2 1]),'-');\n",gs2);
      for (level=0;level<maxlevel;level++){ for (nb=0;nb<nbases;nb++){ for (nb2=0;nb2<nbases;nb2++){
	fprintf(op2,"dist_x2p_%s(%d,%d,%d)=%f;\n",gs2,nb2+1,nb+1,level+1,dist_x2p[nb2+nb*nbases+level*nbases*nbases]);}}}}
    if (verbose>1){ printf(" %% testing discdata_x2p\n");}
    discdata_threshold(discdata_x2p,nbases,maxlevel,nR,outputfilename_disc_x2p);
    if (verbose>1){ printf(" %% testing discdata_depthent\n");}
    discdata_threshold(discdata_depthent,nbases,maxlevel,nR,outputfilename_disc_depthent);
    if (verbose>1){ printf(" %% testing discdata_breadthent\n");}
    discdata_threshold(discdata_breadthent,nbases,maxlevel,nR,outputfilename_disc_breadthent);
    radump(discdata_x2p,"double",nbases*maxlevel,nR,outputfilename_disc_x2p);
    radump(discdata_depthent,"double",nbases*maxlevel,nR,outputfilename_disc_depthent);
    radump(discdata_breadthent,"double",nbases*maxlevel,nR,outputfilename_disc_breadthent);
    tfree(discdata_x2p);discdata_x2p=NULL;
    tfree(discdata_depthent);discdata_depthent=NULL;
    tfree(discdata_breadthent);discdata_breadthent=NULL;}
  if (x_flag==-1){ /* approximation of observation distribution for each event-chain */
    if (verbose){ printf(" %% rereading %s%s... ",filename_base[0],allstring);}
    sprintf(filename[0],"./%s%s",filename_base[0],allstring);
    if (checktofind(filename[0])){ if (verbose){ printf("found\n");} pR = ptreadback(filename[0],0);}
    else /* if not found */{ if (verbose){ printf("warning! couldn't find %s, aborting\n",filename[0]);} exit(EXIT_SUCCESS);}
    for (nb=1;nb<nbases;nb++){
      if (verbose){ printf(" %% adding on %s%s...\n",filename_base[nb],allstring);}
      ptreeplusequals_starter(pR,p0[nb]);}
    pnodeclear_starter(NULL,pR->postree); pnodeclear_starter(NULL,pR->pretree);
    nR=0;
    for (nb=0;nb<nbases;nb++){
      bitbybit=1; continue_flag=1;
      do{
	sprintf(filename[nb],"./%s_%d",filename_base[nb],bitbybit);
	if (verbose>1){ printf(" %% checking %s... ",filename[nb]);}
	if (checktofind(filename[nb])){ 
	  if (verbose>1){ printf("found, initiating tcounter\n");}
	  p1[nb] = ptreadback(filename[nb],0);
	  if (s_flag==1){ if (verbose>1){ printf("shizuffling %d\n",s_flag);} pnodeshizuffle_starter(NULL,p1[nb]->postree,p1[nb]); pnodeshizuffle_starter(NULL,p1[nb]->pretree,p1[nb]);}
	  bitbybit += 1; continue_flag=1;}
	else /* if not found */{ if (verbose>1){ printf("not found, stopping\n");} continue_flag=0;}
	tcounter=1;
	while (tcounter<tclump && continue_flag){
	  sprintf(filename[nb],"./%s_%d",filename_base[nb],bitbybit);
	  if (verbose>1){ printf(" %% checking %s... ",filename[nb]);}
	  if (checktofind(filename[nb])){ 
	    if (verbose>1){ printf("found with tcounter %d\n",tcounter);}
	    bitbybit += 1; continue_flag=1;
	    ptmp = ptreadback(filename[nb],0);
	    if (s_flag==1){ if (verbose>1){ printf("shizuffling %d\n",s_flag);} pnodeshizuffle_starter(NULL,ptmp->postree,ptmp); pnodeshizuffle_starter(NULL,ptmp->pretree,ptmp);}
	    ptreeplusequals_starter(p1[nb],ptmp);
	    ptreetfree(ptmp);ptmp=NULL;}
	  else /* if not found */{ if (verbose>1){ printf("not found, stopping\n");} continue_flag=0;}
	  tcounter += 1;}
	if (continue_flag){
	  if (verbose>1){ printf(" %% accumulated %d-%d total time %0.2f\n",bitbybit-1-tcounter+1,bitbybit-1,p1[nb]->total_time);}
	  nR += 1;
	  if (verbose>1){ printf(" %% calling ptreeobsdist for %s->postree\n",filename[nb]);}
	  ptreeobsdist(nbases,pR,NULL,pR->postree,nb,p1[nb],NULL,p1[nb]->postree,maxlevel,0);}
	if (p1[nb]!=NULL){ ptreetfree(p1[nb]);p1[nb]=NULL;}}
      while (continue_flag);}
    ptreeobsdistend(nbases,nR,pR,NULL,pR->postree,0);
    temphist=histmake(128,1,0);
    pnode2hist_starter(NULL,pR->postree,-1,0,temphist,NULL);
    if (verbose>1){ histprintf(temphist,"discriminability: ");}
    histdump(temphist,0,outputfilename_obsdist," ",0);
    histtfree(temphist);temphist=NULL;
    ptreetfree(pR);pR=NULL;}
  if (x_flag==1){
    if (verbose){ printf(" %% rereading %s%s... ",filename_base[0],allstring);}
    sprintf(filename[0],"./%s%s",filename_base[0],allstring);
    if (checktofind(filename[0])){ if (verbose){ printf("found\n");} pR = ptreadback(filename[0],0);}
    else /* if not found */{ if (verbose){ printf("warning! couldn't find %s, aborting\n",filename[0]);} exit(EXIT_SUCCESS);}
    for (nb=1;nb<nbases;nb++){
      if (verbose){ printf(" %% adding on %s%s...\n",filename_base[nb],allstring);}
      ptreeplusequals_starter(pR,p0[nb]);}
    pnodeclear_starter(NULL,pR->postree); pnodeclear_starter(NULL,pR->pretree);
    if (verbose>1){ printf(" %% making pn0 array and l0pos,l0pre arrays\n");}
    pn0 = (struct pnode **) tcalloc(nbases,sizeof(struct pnode *));
    l0pos = (struct llitem **) tcalloc(nbases,sizeof(struct llitem *));
    l0pre = (struct llitem **) tcalloc(nbases,sizeof(struct llitem *));
    for (nb=0;nb<nbases;nb++){ pn0[nb] = NULL; l0pos[nb] = p0[nb]->postree; l0pre[nb] = p0[nb]->pretree;}
    nR=0;
    for (nb=0;nb<nbases;nb++){
      bitbybit=1; continue_flag=1;
      do{
	sprintf(filename[nb],"./%s_%d",filename_base[nb],bitbybit);
	if (verbose>1){ printf(" %% checking %s... ",filename[nb]);}
	if (checktofind(filename[nb])){ 
	  if (verbose>1){ printf("found, initiating tcounter\n");}
	  p1[nb] = ptreadback(filename[nb],0);
	  if (s_flag==1){ if (verbose>1){ printf("shizuffling %d\n",s_flag);} pnodeshizuffle_starter(NULL,p1[nb]->postree,p1[nb]); pnodeshizuffle_starter(NULL,p1[nb]->pretree,p1[nb]);}
	  bitbybit += 1; continue_flag=1;}
	else /* if not found */{ if (verbose>1){ printf("not found, stopping\n");} continue_flag=0;}
	tcounter=1;
	while (tcounter<tclump && continue_flag){
	  sprintf(filename[nb],"./%s_%d",filename_base[nb],bitbybit);
	  if (verbose>1){ printf(" %% checking %s... ",filename[nb]);}
	  if (checktofind(filename[nb])){ 
	    if (verbose>1){ printf("found with tcounter %d\n",tcounter);}
	    bitbybit += 1; continue_flag=1;
	    ptmp = ptreadback(filename[nb],0);
	    if (s_flag==1){ if (verbose>1){ printf("shizuffling %d\n",s_flag);} pnodeshizuffle_starter(NULL,ptmp->postree,ptmp); pnodeshizuffle_starter(NULL,ptmp->pretree,ptmp);}
	    ptreeplusequals_starter(p1[nb],ptmp);
	    ptreetfree(ptmp);ptmp=NULL;}
	  else /* if not found */{ if (verbose>1){ printf("not found, stopping\n");} continue_flag=0;}
	  tcounter += 1;}
	if (continue_flag){
	  if (verbose>1){ printf(" %% accumulated %d-%d total time %0.2f\n",bitbybit-1-tcounter+1,bitbybit-1,p1[nb]->total_time);}
	  nR += 1;
	  if (verbose>1){ printf(" %% calling ptreex1 for %s->postree\n",filename[nb]);}
	  ptreex1(nbases,p0,pn0,l0pos,nb,p1[nb],NULL,p1[nb]->postree,pR,NULL,pR->postree,maxlevel,0);
	  if (verbose>1){ printf(" %% calling ptreex1 for %s->pretree\n",filename[nb]);}
	  ptreex1(nbases,p0,pn0,l0pre,nb,p1[nb],NULL,p1[nb]->pretree,pR,NULL,pR->pretree,maxlevel,0);}
	if (p1[nb]!=NULL){ ptreetfree(p1[nb]);p1[nb]=NULL;}}
      while (continue_flag);}
    pnodetimesd_starter(NULL,pR->postree,1.0/(double)nR,1.0/(double)nR);
    pnodetimesd_starter(NULL,pR->pretree,1.0/(double)nR,1.0/(double)nR);
    tfree(pn0); pn0=NULL; tfree(l0pos); l0pos=NULL; tfree(l0pre); l0pre=NULL;
    if (verbose>2){ pnodeprintf(NULL,pR->postree,maxlevel,0);}
    if (outputname==NULL){sprintf(filename[0],"./ptree_x1_t%d_%srecord%s",tclump,gs2,allstring);}
    else{ sprintf(filename[0],"./ptree_%s_x1_t%d_%srecord%s",outputname,tclump,gs2,allstring);}
    if (verbose){ printf(" %% dumping %s, with %d records found\n",filename[0],nR);}
    pnodestats_starter(NULL,pR->postree,-1,-1,0,NULL,&wmax,&wmin,NULL,NULL,&rmax,&rmin,NULL,NULL);
    ptreedump_starter(pR,filename[0],2,1,wmax,wmin,rmax,rmin,1);
    ptreetfree(pR);pR=NULL;}
  if (x_flag>=1){
    if (outputname==NULL){sprintf(filename[0],"./ptree_x1_t%d_%srecord%s",tclump,gs2,allstring);}
    else{ sprintf(filename[0],"./ptree_%s_x1_t%d_%srecord%s",outputname,tclump,gs2,allstring);}
    if (verbose){ printf(" %% checking %s...",filename[0]);}
    if (checktofind(filename[0])){
      if (verbose){ printf("found\n");}
      pR = ptreadback(filename[0],0);}
    else /* if not found */{ if (verbose){ printf("not found\n");} exit(EXIT_SUCCESS);}
    if (verbose>1){ printf(" %% making pn0 array and l0pos,l0pre arrays and percentages and nbias\n");}
    pn0 = (struct pnode **) tcalloc(nbases,sizeof(struct pnode *));
    l0pos = (struct llitem **) tcalloc(nbases,sizeof(struct llitem *));
    l0pre = (struct llitem **) tcalloc(nbases,sizeof(struct llitem *));
    for (nb=0;nb<nbases;nb++){ pn0[nb] = NULL; l0pos[nb] = p0[nb]->postree; l0pre[nb] = p0[nb]->pretree;}
    nX=11;
    percentages = (double *) tcalloc(nX,sizeof(double)); 
    pnodestats_starter(NULL,pR->postree,-1,-1,0,NULL,&(percentages[nX-1]),&(percentages[0]),NULL,NULL,NULL,NULL,NULL,NULL);
    percentages[0]=0.50;
    for (nx=1;nx<nX-1;nx++){ percentages[nx] = percentages[0] + (percentages[nX-1]-percentages[0])*(double)nx/(double)(nX-1);}
    temphist=histmake(nX,percentages[nX-1],percentages[0]);
    pnodehist_starter(NULL,pR->postree,-1,0,temphist,NULL);
    if (verbose){
      for (nx=0;nx<nX;nx++){ fprintf(op," percentages(%d)=%0.16lf;\n",nx+1,percentages[nx]);}
      for (nx=0;nx<temphist->nbins;nx++){
	fprintf(op," x1_%s_hist_T%dS%d(%d) = %0.16lf;\n",gs2,tclump,s_flag,nx+1,temphist->data[nx]);}}
    histtfree(temphist);temphist=NULL;
    nbias = (double *) tcalloc(nX*maxlevel,sizeof(double));
    nR=0;
    for (nb=0;nb<nbases;nb++){
      bitbybit=1; continue_flag=1;
      do{
	sprintf(filename[nb],"./%s_%d",filename_base[nb],bitbybit);
	if (verbose>1){ printf(" %% checking %s... ",filename[nb]);}
	if (checktofind(filename[nb])){ 
	  if (verbose>1){ printf("found, initiating tcounter\n");}
	  p1[nb] = ptreadback(filename[nb],0);
	  if (s_flag==1){ if (verbose>1){ printf("shizuffling %d\n",s_flag);} pnodeshizuffle_starter(NULL,p1[nb]->postree,p1[nb]); pnodeshizuffle_starter(NULL,p1[nb]->pretree,p1[nb]);}
	  bitbybit += 1; continue_flag=1;}
	else /* if not found */{ if (verbose>1){ printf("not found, stopping\n");} continue_flag=0;}
	tcounter=1;
	while (tcounter<tclump && continue_flag){
	  sprintf(filename[nb],"./%s_%d",filename_base[nb],bitbybit);
	  if (verbose>1){ printf(" %% checking %s... ",filename[nb]);}
	  if (checktofind(filename[nb])){ 
	    if (verbose>1){ printf("found with tcounter %d\n",tcounter);}
	    bitbybit += 1; continue_flag=1;
	    ptmp = ptreadback(filename[nb],0);
	    if (s_flag==1){ if (verbose>1){ printf("shizuffling %d\n",s_flag);} pnodeshizuffle_starter(NULL,ptmp->postree,ptmp); pnodeshizuffle_starter(NULL,ptmp->pretree,ptmp);}
	    ptreeplusequals_starter(p1[nb],ptmp);
	    ptreetfree(ptmp);ptmp=NULL;}
	  else /* if not found */{ if (verbose>1){ printf("not found, stopping\n");} continue_flag=0;}
	  tcounter += 1;}
	if (continue_flag){
	  if (verbose>1){ printf(" %% accumulated %d-%d total time %0.2f\n",bitbybit-1-tcounter+1,bitbybit-1,p1[nb]->total_time);}
	  nR += 1;
	  if (verbose>1){ printf(" %% calling ptreex1sbias for %s->postree (base number %d)\n",filename[nb],nb);}
	  if (verbose>1){ printf(" %% creating sbias\n");}
	  sbias = (double *) tcalloc(nbases*nX*maxlevel,sizeof(double));
	  ptreex1sbias(nbases,p0,pn0,l0pos,nb,p1[nb],NULL,p1[nb]->postree,pR,NULL,pR->postree,maxlevel,0,nX,percentages,sbias);
	  if (verbose>2){ for (level=0;level<maxlevel;level++){ printf(" %% sbias nl %d given by:\n",level); raprintf(&(sbias[0 + 0*nbases + level*nbases*nX]),"double",nbases,nX," %% ");}}
	  if (verbose>2){ printf(" %% now turning nbias from:\n"); raprintf(nbias,"double",nX,maxlevel," %% ");}
	  for (level=0;level<maxlevel;level++){ for (nx=0;nx<nX;nx++){
	    stats("double",&(sbias[0 + nx*nbases + level*nbases*nX]),nbases,&maxbias,NULL,NULL,NULL);
	    if (verbose>1){ printf(" %% nl %d nx %d found maxbias %f...",level,nx,maxbias);}
	    nmin=0;
	    for (nb2=0;nb2<nbases;nb2++){
	      if (maxbias>0 && fabs(sbias[nb2 + nx*nbases + level*nbases*nX]-maxbias)<biastol){ nmin+=1;}}
	    if (maxbias>0 && fabs(sbias[nb + nx*nbases + level*nbases*nX]-maxbias)<biastol){
	      if (verbose>1){ printf("at correct index %d\n",nb);}
	      nbias[nx + level*nX] += 1.0/(double)nmin;}
	    else /* if (fabs(sbias[nb+nx*nbases+level*nbases*nX]-maxbias)>=biastol) */{
	      if (verbose>1){ printf("at wrong index\n");}}}}
	  if (verbose>2){ printf(" %% into:\n"); raprintf(nbias,"double",nX,maxlevel," %% ");}
	  tfree(sbias); sbias=NULL;}
	if (p1[nb]!=NULL){ ptreetfree(p1[nb]);p1[nb]=NULL;}}
      while (continue_flag);}
    if (verbose){ 
      printf(" %% number of records = %d\n",nR);
      fprintf(op," nR_%s_T(%d) = %d;\n",gs2,tclump,nR);
      for (level=0;level<maxlevel;level++){ for (nx=0;nx<nX;nx++){ 
	printf(" %% level %d nx %d percentage %0.2f nbias %0.1f\n",level+1,nx,percentages[nx],nbias[nx+level*nX]);}}
      for (level=0;level<maxlevel;level++){ for (nx=0;nx<nX;nx++){
	fprintf(op," x1_%s_STLPB(%d,%d,%d,%d)=%0.1f;\n",gs2,s_flag+1,tclump,level+1,nx+1,nbias[nx+level*nX]);}}
      for (level=0;level<maxlevel;level++){
	fprintf(op," temp=x1_%s_STLPB(%d,%d,1:%d,:)/nR_%s_T(%d);\n",gs2,s_flag+1,tclump,level+1,gs2,tclump);
	fprintf(op," x1_%s_STL(%d,%d,%d)=max(temp(:));\n",gs2,s_flag+1,tclump,level+1);}
      fprintf(op," plot(permute(x1_%s_STL(s_to_plot,:,:),[3 2 1]),':');\n",gs2);}
    tfree(nbias);nbias=NULL;
    tfree(percentages);percentages=NULL;
    tfree(pn0); pn0=NULL; tfree(l0pos); l0pos=NULL; tfree(l0pre); l0pre=NULL;
    ptreetfree(pR);pR=NULL;}
  for (nb=0;nb<nbases;nb++){
    ptreetfree(p0[nb]);p0[nb]=NULL;
    tfree(filename[nb]);}
  tfree(p0);tfree(p1);
  tfree(filename);
  tfree(temp_x2);
  tfree(temp_x2p);
  tfree(temp_fabs);
  tfree(temp_fabs2);
  tfree(temp_depthent);
  tfree(temp_breadthent);
  tfree(nright_x2);tfree(nwrong_x2);
  tfree(nright_x2p);tfree(nwrong_x2p);
  tfree(nright_fabs);tfree(nwrong_fabs);
  tfree(nright_fabs2);tfree(nwrong_fabs2);
  tfree(nright_depthent);tfree(nwrong_depthent);
  tfree(nright_breadthent);tfree(nwrong_breadthent);
  tfree(discra);
  tfree(dist_x2p);
  tfree(nright_x2p_perstimulus);tfree(nwrong_x2p_perstimulus);
  tfree(nright_depthent_perstimulus);tfree(nwrong_depthent_perstimulus);
  tfree(nright_breadthent_perstimulus);tfree(nwrong_breadthent_perstimulus);
  tfree(referencenormalizers);
  tfree(observednormalizers);
  tfree(fp);
  GLOBAL_PTREE_REGION_TYPE=old_region_type;
  if (op!=stdout){ fclose(op); op=NULL;}
  if (op2!=stdout){ fclose(op2); op2=NULL;}
}

void ptree_trialaverage(int argc,char **argv)
{
  int verbose=2;
  int still_have_options = argc-2;
  int nbases=0,nb=0;
  char helpfile[1024];
  char **filename_base=NULL;
  int maxlevel=1,tclump=1;
  int X_flag=1,x_flag=1,s_flag=0,d_flag=0;
  char *outputname=NULL;
  if (verbose>1){ printf(" %% [entering ptree_trialaverage]\n");}
  sprintf(helpfile,"%% [-Vverboselevel] -Nnumber_of_bases [-Mmaxlevel] [-Ttclump] [-XX_flag] [-xx_flag] [-ss_flag] [-dd_flag] -Fbase_name1 -Fbase_name2 ... [-Ooutputname]\n");
  if (still_have_options==0){ printf(helpfile); exit(EXIT_SUCCESS);}
  while (still_have_options){
    switch(getopt(argc,argv,"N:n:F:f:M:m:T:t:V:v:X:x:s:d:O:o:")){
    case 'M': case 'm':
      maxlevel = maximum(1,atoi(optarg));
      if (verbose){ printf(" %% maxlevel set to %d\n",maxlevel);}
      break;
    case 'T': case 't':
      tclump = maximum(1,atoi(optarg));
      if (verbose){ printf(" %% tclump set to %d\n",tclump);}
      break;
    case 'V': case 'v':
      verbose = atoi(optarg);
      if (verbose){ printf(" %% verboselevel set to %d\n",verbose);}
      break;
    case 'N': case 'n': 
      nbases = maximum(0,atoi(optarg));
      filename_base = (char **) tcalloc(nbases,sizeof(char *));
      for (nb=0;nb<nbases;nb++){ filename_base[nb] = (char *) tcalloc(256,sizeof(char));}
      nb=0;
      if (verbose){ printf(" %% nbases read as %d\n",nbases);}
      break;
    case 'F': case 'f': 
      sprintf(filename_base[nb],"%s",optarg); 
      if (verbose){ printf(" %% filename_base %d given by %s...",nb,filename_base[nb]);}
      nb = minimum(nbases-1,nb+1);
      if (verbose){ printf(" now nb %d\n",nb);}
      break;
    case 'X':
      X_flag = atoi(optarg);
      if (verbose){ printf(" %% now X_flag=%d\n",X_flag);}
      break;
    case 'x':
      x_flag = atoi(optarg);
      if (verbose){ printf(" %% now x_flag=%d\n",x_flag);}
      break;
    case 's':
      s_flag = atoi(optarg);
      if (verbose){ printf(" %% now s_flag=%d\n",s_flag);}
      break;
    case 'd':
      d_flag = atoi(optarg);
      if (verbose){ printf(" %% now d_flag=%d\n",d_flag);}
      break;
    case 'O': case 'o': 
      if (outputname==NULL){ outputname = (char *) tcalloc(128,sizeof(char));}
      sprintf(outputname,"%s",optarg);
      if (verbose){ printf(" %% now outputname=%s\n",outputname);}
      break;
    default: printf(helpfile); exit(EXIT_SUCCESS); break;}
    still_have_options--;}
  if (verbose>1){ printf(" %% calling ptree_trialaverage_helper\n");}
  ptree_trialaverage_helper(nbases,maxlevel,tclump,verbose,X_flag,x_flag,s_flag,d_flag,filename_base,outputname);
  for (nb=0;nb<nbases;nb++){ tfree(filename_base[nb]);} tfree(filename_base);
}

void pnode2distribution_starter(struct ptree *p0,struct pnode *parent0,struct llitem *l0,struct ptree *p1,struct pnode *parent1,struct llitem *l1,int maxlevel,int level)
{
  /* adds weights of tree *p1 to llist stored using *temp in reference tree *p0 */
  int firstentry=0;
  struct pnode *pn0=NULL,*pn1=NULL;
  struct llitem *l2=NULL;
  double *temp=NULL;
  if (parent1!=NULL && l1->parent==NULL){ /* first descent into parent1->childllitem */ 
    assert(parent1->childllitem==l1); firstentry=1;}
  if (l1->kidl!=NULL){ pnode2distribution_starter(p0,parent0,l0,p1,parent1,l1->kidl,maxlevel,level);}
  if (l1->kidr!=NULL){ pnode2distribution_starter(p0,parent0,l0,p1,parent1,l1->kidr,maxlevel,level);}
  if (l1->item!=NULL){
    if (maxlevel==-1 || level<maxlevel){ 
      pn1=(struct pnode *)l1->item; 
      if ((l2=llitemaddorfind(0,l0,pn1->region,&region2pnode_compare_label))!=NULL){
	pn0=(struct pnode *)l2->item;
	pnode2distribution_starter(p0,pn0,pn0->childllitem,p1,pn1,pn1->childllitem,maxlevel,level+1);}}}
  if (firstentry){ if (parent0!=NULL && parent1!=NULL){ 
    if (parent0->temp==NULL){ parent0->temp=(void *)llistmake();}
    temp = (double *) tcalloc(1,sizeof(double)); *temp = parent1->weight; litemadd((struct llist *)parent0->temp,temp);}}
}

void pnode2distribution_helper(struct ptree *p0,struct pnode *parent0,struct llitem *l0,int maxlevel,int level,struct llist **Lra,int nm,int nR,struct hist **Phistra)
{
  /* constructs nm llists from the weight llists stored in *temp of reference tree *p0 
     Lra[moment + level*nm] stores moments obtained by various observations at level 
     in addition constructs the "average" distribution Phistra, 
     assuming that (struct llist *)(parent0->temp)->length<=nR and Phistra is maxlevel elements long */
  int firstentry=0;
  struct pnode *pn0=NULL;
  int nr=0,nl=0;
  double *mra=NULL,*temp=NULL;
  struct llist *L=NULL;
  struct litem *l=NULL;
  if (parent0!=NULL && l0->parent==NULL){ /* first descent into parent0->childllitem */ assert(parent0->childllitem==l0); firstentry=1;}
  if (l0->kidl!=NULL){ pnode2distribution_helper(p0,parent0,l0->kidl,maxlevel,level,Lra,nm,nR,Phistra);}
  if (l0->kidr!=NULL){ pnode2distribution_helper(p0,parent0,l0->kidr,maxlevel,level,Lra,nm,nR,Phistra);}
  if (l0->item!=NULL){ if (maxlevel==-1||level<maxlevel){ pn0=(struct pnode *)l0->item; pnode2distribution_helper(p0,pn0,pn0->childllitem,maxlevel,level+1,Lra,nm,nR,Phistra);}}
  if (firstentry){ if (parent0!=NULL){ if (parent0->temp!=NULL){ 
    L = (struct llist *)parent0->temp;
    mra = (double *) tcalloc(maximum(3,nm),sizeof(double));
    lliststats2(L,mra,maximum(3,nm),nR-L->length);
    for (nr=0;nr<nm;nr++){ temp = (double *) tcalloc(1,sizeof(double)); *temp = mra[nr]; litemadd(Lra[nr + level*nm],temp);}
    l=L->first;nr=0;
    while (l!=NULL && nr<nR){
      for (nl=level-1;nl<maxlevel;nl++){ histadd(Phistra[nl],((*(double *)l->item)-mra[1])/sqrt(mra[2]),1);}
      nr+=1;
      l=l->child;}
    for (nl=level-1;nl<maxlevel;nl++){ histadd(Phistra[nl],(0-mra[1])/sqrt(mra[2]),nR-nr);}
    tfree(mra);mra=NULL;}}}
}

void pnode2distribution_ender(struct ptree *p0,struct pnode *parent0,struct llitem *l0,int maxlevel,int level)
{
  /* frees the llists stored in *temp of reference tree *p0 */
  int firstentry=0;
  struct pnode *pn0=NULL;
  if (parent0!=NULL && l0->parent==NULL){ /* first descent into parent0->childllitem */  assert(parent0->childllitem==l0); firstentry=1;}
  if (l0->kidl!=NULL){ pnode2distribution_ender(p0,parent0,l0->kidl,maxlevel,level);}
  if (l0->kidr!=NULL){ pnode2distribution_ender(p0,parent0,l0->kidr,maxlevel,level);}
  if (l0->item!=NULL){
    if (maxlevel==-1||level<maxlevel){ pn0=(struct pnode *)l0->item; pnode2distribution_ender(p0,pn0,pn0->childllitem,maxlevel,level+1);}}
  if (firstentry){ if (parent0!=NULL){ if (parent0->temp!=NULL){ llisttfree2((struct llist *)parent0->temp);parent0->temp=NULL;}}}
}

void ptree_trialdistribution_helper(int nbases,int maxlevel,int tclump,int verbose,char **filename_base,char *outputname)
{
  /* assumes that outputname does NOT start with "./" */
  int old_region_type;
  int s_flag=0;
  int nmoments=4;
  int nb=0,nm=0,nl=0;
  char *filename=NULL,allstring[16];
  struct ptree **p0=NULL,**p1=NULL,*ptmp=NULL;
  int bitbybit=0,continue_flag=0,tcounter=0;
  struct llist **Lra=NULL;
  double max=0,min=0,mean=0,stdev=0;
  struct hist *temphist=NULL;
  int nR=0;
  struct hist **Phistra=NULL;
  if (verbose>1){ 
    printf(" %% [entering ptree_trialdistribution_helper]\n");
    printf(" %% nbases %d maxlevel %d tclump %d verbose %d s_flag %d\n",nbases,maxlevel,tclump,verbose,s_flag);
    for (nb=0;nb<nbases;nb++){ printf(" %% filename_base[%d]=%s\n",nb,filename_base[nb]);}
    printf(" %% outputname=%s\n",outputname);}
  sprintf(allstring,"_s%d_all",s_flag);
  filename = (char *) tcalloc(256,sizeof(char));
  p0 = (struct ptree **) tcalloc(nbases,sizeof(struct ptree *));
  p1 = (struct ptree **) tcalloc(nbases,sizeof(struct ptree *));
  old_region_type=GLOBAL_PTREE_REGION_TYPE; GLOBAL_PTREE_REGION_TYPE=0;
  if (verbose>1){ printf(" %% setting GLOBAL_PTREE_REGION_TYPE to %d\n",GLOBAL_PTREE_REGION_TYPE);}
  for (nb=0;nb<nbases;nb++){
    if (verbose){ printf(" %% checking for %s%s... ",filename_base[nb],allstring);}
    sprintf(filename,"./%s%s",filename_base[nb],allstring);
    if (checktofind(filename)){ if (verbose){ printf("found\n");} p0[nb] = ptreadback(filename,0);}
    else /* if not found */{ printf("error, not found, exiting\n"); exit(EXIT_FAILURE);}}
  for (nb=0;nb<nbases;nb++){
    bitbybit=1; continue_flag=1;nR=0;
    do{
      sprintf(filename,"./%s_%d",filename_base[nb],bitbybit); if (verbose>1){ printf(" %% checking %s... ",filename);}
      if (checktofind(filename)){ 
	if (verbose>1){ printf("found, initiating tcounter\n");} 
	p1[nb] = ptreadback(filename,0); if (s_flag==1){ if (verbose>1){ printf("shizuffling %d\n",bitbybit);} pnodeshizuffle_starter(NULL,p1[nb]->postree,p1[nb]); pnodeshizuffle_starter(NULL,p1[nb]->pretree,p1[nb]);}
	bitbybit += 1; continue_flag=1;}
      else /* if not found */{ continue_flag=0;}
      tcounter=1;
      while (tcounter<tclump && continue_flag){
	sprintf(filename,"./%s_%d",filename_base[nb],bitbybit); if (verbose>1){ printf(" %% checking %s... ",filename);}
	if (checktofind(filename)){
	  if (verbose>1){ printf("found, with tcounter %d\n",tcounter);} 
	  ptmp = ptreadback(filename,0); if (s_flag==1){ if (verbose>1){ printf("shizuffling %d\n",bitbybit);} pnodeshizuffle_starter(NULL,ptmp->postree,ptmp); pnodeshizuffle_starter(NULL,ptmp->pretree,ptmp);}
	  ptreeplusequals_starter(p1[nb],ptmp);  
	  ptreetfree(ptmp);ptmp=NULL;
	  bitbybit += 1; continue_flag=1;}
	else /* if not found */{ if (verbose>1){ printf("not found, stopping\n");} continue_flag=0;}
	tcounter += 1;}
      if (continue_flag){
	nR += 1;
	if (verbose>1){ printf(" %% accumulated %d-%d total time %0.2f\n",bitbybit-1-tcounter+1,bitbybit-1,p1[nb]->total_time);}
	pnode2distribution_starter(p0[nb],NULL,p0[nb]->postree,p1[nb],NULL,p1[nb]->postree,maxlevel,0);}
      if (p1[nb]!=NULL){ ptreetfree(p1[nb]);p1[nb]=NULL;}}
    while (continue_flag);
    if (verbose){ printf(" %% found %d total records\n",nR);}
    Lra = (struct llist **) tcalloc((maxlevel+1)*nmoments,sizeof(struct llist *));
    for (nm=0;nm<nmoments*(maxlevel+1);nm++){ Lra[nm]=llistmake();}
    Phistra=(struct hist **) tcalloc(maxlevel,sizeof(struct hist *));
    for (nl=0;nl<maxlevel;nl++){ Phistra[nl]=histmake(128,0+16*STD_VIEW,0-1*STD_VIEW);}
    pnode2distribution_helper(p0[nb],NULL,p0[nb]->postree,maxlevel,0,Lra,nmoments,nR,Phistra);
    for (nl=0;nl<maxlevel;nl++){
      sprintf(filename,"./%s_dist_Phist_level%d",filename_base[nb],nl);
      histdump(Phistra[nl],0,filename," ",0); histdump(Phistra[nl],1,filename," ",0);}
    for (nl=0;nl<maxlevel;nl++){ histtfree(Phistra[nl]);Phistra[nl]=NULL;} tfree(Phistra);Phistra=NULL;
    for (nm=0;nm<nmoments;nm++){ for (nl=0;nl<maxlevel;nl++){ 
      lliststats(Lra[nm + nl*nmoments],&max,&min,&mean,&stdev);
      temphist = histmake(128,mean+STD_VIEW*stdev,mean-STD_VIEW*stdev);
      llist2hist(Lra[nm + nl*nmoments],temphist);
      sprintf(filename,"./%s_dist_moment%d_level%d",filename_base[nb],nm,nl); 
      histdump(temphist,0,filename," ",0); histdump(temphist,1,filename," ",0);
      histtfree(temphist);temphist=NULL;}}
    pnode2distribution_ender(p0[nb],NULL,p0[nb]->postree,maxlevel,0);
    for (nm=0;nm<nmoments*(maxlevel+1);nm++){ llisttfree2(Lra[nm]);Lra[nm]=NULL;}
    tfree(Lra);Lra=NULL;}
  for (nb=0;nb<nbases;nb++){
    ptreetfree(p0[nb]);p0[nb]=NULL;}
  tfree(p0);tfree(p1);
  tfree(filename);
  GLOBAL_PTREE_REGION_TYPE=old_region_type;
}

void ptree_trialdistribution(int argc,char **argv)
{
  int verbose=2;
  int still_have_options = argc-2;
  int nbases=0,nb=0;
  char helpfile[1024];
  char **filename_base=NULL;
  int maxlevel=1,tclump=1;
  char *outputname=NULL;
  if (verbose>1){ printf(" %% [entering ptree_trialdistribution]\n");}
  sprintf(helpfile,"%% [-Vverboselevel] -Nnumber_of_bases [-Mmaxlevel] [-Ttclump] -Fbase_name1 -Fbase_name2 ... [-Ooutputname]\n");
  if (still_have_options==0){ printf(helpfile); exit(EXIT_SUCCESS);}
  while (still_have_options){
    switch(getopt(argc,argv,"N:n:F:f:M:m:T:t:V:v:X:x:s:O:o:")){
    case 'M': case 'm':
      maxlevel = maximum(1,atoi(optarg));
      if (verbose){ printf(" %% maxlevel set to %d\n",maxlevel);}
      break;
    case 'T': case 't':
      tclump = maximum(1,atoi(optarg));
      if (verbose){ printf(" %% tclump set to %d\n",tclump);}
      break;
    case 'V': case 'v':
      verbose = atoi(optarg);
      if (verbose){ printf(" %% verboselevel set to %d\n",verbose);}
      break;
    case 'N': case 'n': 
      nbases = maximum(0,atoi(optarg));
      filename_base = (char **) tcalloc(nbases,sizeof(char *));
      for (nb=0;nb<nbases;nb++){ filename_base[nb] = (char *) tcalloc(256,sizeof(char));}
      nb=0;
      if (verbose){ printf(" %% nbases read as %d\n",nbases);}
      break;
    case 'F': case 'f': 
      sprintf(filename_base[nb],"%s",optarg); 
      if (verbose){ printf(" %% filename_base %d given by %s...",nb,filename_base[nb]);}
      nb = minimum(nbases-1,nb+1);
      if (verbose){ printf(" now nb %d\n",nb);}
      break;
    case 'O': case 'o': 
      if (outputname==NULL){ outputname = (char *) tcalloc(128,sizeof(char));}
      sprintf(outputname,"%s",optarg);
      if (verbose){ printf(" %% now outputname=%s\n",outputname);}
      break;
    default: printf(helpfile); exit(EXIT_SUCCESS); break;}
    still_have_options--;}
  if (verbose>1){ printf(" %% calling ptree_trialdistribution_helper\n");}
  ptree_trialdistribution_helper(nbases,maxlevel,tclump,verbose,filename_base,outputname);
  for (nb=0;nb<nbases;nb++){ tfree(filename_base[nb]);} tfree(filename_base);
}

struct obsdisthist * obsdisthistmake(int ninputs)
{
  struct obsdisthist *odh=NULL;
  int nh=0;
  odh = (struct obsdisthist *) tmalloc(sizeof(struct obsdisthist));
  odh->length = ninputs;
  odh->hra = (struct hist **) tcalloc(odh->length,sizeof(struct hist *));
  for (nh=0;nh<odh->length;nh++){ odh->hra[nh] = histmake(0,0,0);}
  odh->nrecords=0;
  odh->value = (double *) tcalloc(odh->length*odh->length,sizeof(double));
  odh->ballot_2 = (double *) tcalloc(odh->length,sizeof(double));
  odh->ballot_3 = (double *) tcalloc(odh->length*odh->length,sizeof(double));
  return odh;
}

void obsdisthisttfree(struct obsdisthist *odh)
{
  int nh=0;
  for (nh=0;nh<odh->length;nh++){ histtfree(odh->hra[nh]); odh->hra[nh]=NULL;} tfree(odh->hra); odh->hra=NULL;
  tfree(odh->value); odh->value = NULL;
  tfree(odh->ballot_2); odh->ballot_2 = NULL;
  tfree(odh->ballot_3); odh->ballot_3 = NULL;
  tfree(odh); odh=NULL;
}

void obsdisthist2file(void *void_odh,FILE *fp)
{
  int plusone=+1,minusone=-1;
  int nh=0;
  struct obsdisthist *odh=NULL;
  if (void_odh==NULL){ fwrite(&minusone,sizeof(int),1,fp);}
  else /* if (void_odh!=NULL) */{
    fwrite(&plusone,sizeof(int),1,fp);
    odh=(struct obsdisthist *)void_odh;
    fwrite(&(odh->length),sizeof(int),1,fp);
    for (nh=0;nh<odh->length;nh++){ hist2file(odh->hra[nh],fp);}
    fwrite(&(odh->nrecords),sizeof(int),1,fp);
    fwrite(&(odh->value[0]),sizeof(double),odh->length*odh->length,fp);
    fwrite(&(odh->ballot_2[0]),sizeof(double),odh->length,fp);
    fwrite(&(odh->ballot_3[0]),sizeof(double),odh->length*odh->length,fp);}
}

void * file2obsdisthist(FILE *fp,int *no_error_passback)
{
  struct obsdisthist *odh=NULL;
  int nh=0;
  int label=0,no_error=0,no_minor_error=0,temp_no_error=0;
  no_error = 1;
  if (no_error){ no_error = (fread(&label,sizeof(int),1,fp)==(size_t)1) && !feof(fp) && !ferror(fp);}
  if (no_error && label==+1){
    odh = (struct obsdisthist *) tmalloc(sizeof(struct obsdisthist));
    if (no_error){ no_error = (fread(&(odh->length),sizeof(int),1,fp)==(size_t)1) && !feof(fp) && !ferror(fp);}
    if (odh->length>=0){
      odh->hra = (struct hist **) tcalloc(odh->length,sizeof(struct hist *));
      no_minor_error=1;
      for (nh=0;nh<odh->length;nh++){ 
	odh->hra[nh]=(struct hist *)file2hist(fp,&temp_no_error); 
	no_minor_error = no_minor_error && temp_no_error;}
      no_error = no_error && no_minor_error;
      if (no_error){ no_error = (fread(&(odh->nrecords),sizeof(int),1,fp)==(size_t)1) && !feof(fp) && !ferror(fp);}
      odh->value = (double *) tcalloc(odh->length*odh->length,sizeof(double));
      if (no_error){ no_error = (fread(&(odh->value[0]),sizeof(double),odh->length*odh->length,fp)==(size_t)(odh->length*odh->length)) && !feof(fp) && !ferror(fp);}
      odh->ballot_2 = (double *) tcalloc(odh->length,sizeof(double));
      if (no_error){ no_error = (fread(&(odh->ballot_2[0]),sizeof(double),odh->length,fp)==(size_t)(odh->length)) && !feof(fp) && !ferror(fp);}
      odh->ballot_3 = (double *) tcalloc(odh->length*odh->length,sizeof(double));
      if (no_error){ no_error = (fread(&(odh->ballot_3[0]),sizeof(double),odh->length*odh->length,fp)==(size_t)(odh->length*odh->length)) && !feof(fp) && !ferror(fp);}}
    if (!no_error){
      printf(" %% warning! improperly read obsdisthist in file2obsdisthist\n");
      obsdisthisttfree(odh);odh=NULL;}}
  if (no_error_passback!=NULL){ *no_error_passback=no_error;}
  return (void *)odh;
}

void obsdisthistupdate(struct obsdisthist *odh,int ni,double val,double ntimes)
{
  /* adds integer value val to odh->hra[ni], expanding odh->hra[input]->data as necessary */
  int verbose=0;
  int bin=0;
  int v=0,d=0,nb=0;
  struct hist *h=odh->hra[ni];
  if (verbose){ printf(" %% [entering obsdisthistupdate] ni %d v %d ntimes %d\n",ni,(int)val,(int)ntimes);}
  if (verbose){ histprintf(h," %% starting: ");}
  if (finite(ntimes) && finite(val) && ntimes>0){
    if (h->nbins>0){
      v = (int)val;
      if (v+1>(int)h->max){
	if (verbose){ printf(" %% v %d >= h->max %d ",v,(int)h->max);}
	d = v+1-(int)h->max;
	h->max = v+1;
	h->nbins += d;
	if (verbose){ printf(" now h->max %d, h->nbins %d\n",(int)h->max,h->nbins);}
	h->data = (double *) trealloc(h->data,h->nbins*sizeof(double));
	for (nb=h->nbins-d;nb<h->nbins;nb++){ h->data[nb]=0;}}
      else if (v<(int)h->min){
	if (verbose){ printf(" %% v %d < h->min %d ",v,(int)h->min);}
	d = (int)h->min-v;
	h->min = v;
	h->nbins += d;
	h->data = (double *) trealloc(h->data,h->nbins*sizeof(double));
	if (verbose){ printf(" now h->min %d, h->nbins %d\n",(int)h->min,h->nbins);}
	for (nb=h->nbins-1;nb>=d;nb--){ h->data[nb] = h->data[nb-d];} for (nb=0;nb<d;nb++){ h->data[nb]=0;}}
      bin = v-(int)h->min;
      bin = maximum(0,minimum(bin,h->nbins-1));
      h->data[bin] += ntimes;
      h->datasum += ntimes*v;
      h->ndatum += ntimes;}
    else /* if (h->nbins==0) */{ 
      h->nbins=2; 
      h->data=(double *) tcalloc(2,sizeof(double)); 
      h->data[0]=(int)ntimes; h->min=(int)val; h->max=h->min+1;
      h->datasum=(int)ntimes*(int)val; h->ndatum=ntimes;}}
  if (verbose){ histprintf(h," %% ending: ");}
  if (verbose){ printf(" %% [exiting obsdisthistupdate]\n");}
}

void obsdisthistprintf(struct obsdisthist *odh,char *prefix)
{
  int nh=0;
  for (nh=0;nh<odh->length;nh++){ printf(" %s odh->hra[%d]:\n",prefix,nh); histprintf(odh->hra[nh],prefix);}
  printf(" %s odh->value:\n",prefix);
  raprintf(odh->value,"double",odh->length,odh->length,prefix);
}

double obsdisthistappraise_helper_helper(int index,struct hist *h1,int zerofill1,struct hist *h2,int zerofill2)
{
  /* given observation distribution (zerofill1,h1) and (zerofill2,h2), return if index (of h1) is classified as h1,
     0 means no, 1 means yes, 0.5 means can't tell.
     assume that (if index>=0), h1 is nonempty */
  int verbose=0;
  int lag1=0,lead1=0,lag2=0,lead2=0,bin1=0,bin2=0,lag1d=0,lead1d=0,lag2d=0,lead2d=0,lay1d=0,lay2d=0;
  int lag1_under=0,lag2_under=0,lead1_over=0,lead2_over=0;
  double lay1=0,lay2=0;
  int lag2_split=0;
  double classify_correctly=0;
  if (verbose){ 
    printf(" %% [entering obsdisthistappraise_helper_helper], index %d, z1 %d, z2 %d\n",index,zerofill1,zerofill2);
    histprintf(h1," %% %% h1:"); histprintf(h2," %% %% h2:");}
  if (index<0){
    if (zerofill1>zerofill2){ classify_correctly=1;}
    else if (zerofill1<zerofill2){ classify_correctly=0;}
    else if (zerofill1==zerofill2){ 
      lead1=0; lead1_over=0; while (lead1<h1->nbins && (int)h1->data[lead1]<=0){ lead1+=1;} if (lead1>=h1->nbins){ lead1_over=1;}
      lead2=0; lead2_over=0; while (lead2<h2->nbins && (int)h2->data[lead2]<=0){ lead2+=1;} if (lead2>=h2->nbins){ lead2_over=1;}
      if (verbose){ printf(" %% ended with lead1=%d,lead2=%d\n",lead1,lead2);}
      if (lead1_over && lead2_over){ classify_correctly=0.5;}
      else if (lead1_over && !lead2_over){ classify_correctly=0;}
      else if (!lead1_over && lead2_over){ classify_correctly=1;}
      else if (!lead1_over && !lead2_over){
	if ((int)h1->min+lead1 < (int)h2->min+lead2){ classify_correctly=1;}
	else if ((int)h1->min+lead1 > (int)h2->min+lead2){ classify_correctly=0;}
	else if ((int)h1->min+lead1 == (int)h2->min+lead2){ 
	  if ((int)h1->data[lead1] > (int)h2->data[lead2]){ classify_correctly=1;}
	  else if ((int)h1->data[lead1] < (int)h2->data[lead2]){ classify_correctly=0;}
	  else if ((int)h1->data[lead1] == (int)h2->data[lead2]){ classify_correctly=0.5;}}}}}
  else if (index>=0){
    bin1 = maximum(0,minimum(index,h1->nbins-1));
    lag1 = bin1; lag1d=0; lag1_under=0;
    while (lag1>=0 && (int)h1->data[lag1]<=0){ lag1-=1; lag1d+=1;} if (lag1<0){ lag1_under=1; lag1d=(int)h1->min+bin1;}
    lead1 = bin1; lead1d=0; lead1_over=0;
    while (lead1<h1->nbins && (int)h1->data[lead1]<=0){ lead1+=1; lead1d-=1;} if (lead1>=h1->nbins){ lead1_over=1;}
    if (!lag1_under && !lead1_over){ 
      if (abs(lag1d)<abs(lead1d)){ lay1=(int)h1->data[lag1]; lay1d=abs(lag1d);}
      else if (abs(lag1d)>abs(lead1d)){ lay1=(int)h1->data[lead1]; lay1d=abs(lead1d);}
      else if (abs(lag1d)==abs(lead1d)){ lay1=((int)h1->data[lag1]+(int)h1->data[lead1])/2.0; lay1d=abs(lead1d);}}
    else if (lag1_under && !lead1_over){
      if (abs(lag1d)<abs(lead1d)){ 
	if (zerofill1>0){ lay1=zerofill1; lay1d=abs(lag1d);} 
	else if (zerofill1<=0){ lay1=(int)h1->data[lead1]; lay1d=abs(lead1d);}}
      else if (abs(lag1d)>abs(lead1d)){ lay1=(int)h1->data[lead1]; lay1d=abs(lead1d);}
      else if (abs(lag1d)==abs(lead1d)){ 
	if (zerofill1>0){ lay1=(zerofill1+(int)h1->data[lead1])/2.0; lay1d=abs(lead1d);}
	else if (zerofill1<=0){ lay1=(int)h1->data[lead1]; lay1d=abs(lead1d);}}}
    else if (!lag1_under && lead1_over){
      lay1=(int)h1->data[lag1]; lay1d=abs(lag1d);}
    else if (lag1_under && lead1_over){
      if (zerofill1>0){ lay1=zerofill1; lay1d=abs(lag1d);}
      else if (zerofill1<=0){ lay1=0; lay1d=abs(lag1d);}}    
    if (verbose){ printf(" %% ended with lag1=%d,lag1d=%d,lead1=%d,lead1d=%d,lay1=%0.1f,lay1d=%d\n",lag1,lag1d,lead1,lead1d,lay1,lay1d);}
    if (h2->nbins>0){
      bin2 = ((int)h1->min+bin1) - (int)h2->min;
      lag2 = maximum(0,minimum(h2->nbins-1,bin2)); lag2d=((int)h1->min+bin1) - ((int)h2->min+lag2); lag2_under=0;
      while (lag2>=0 && (int)h2->data[lag2]<=0){ lag2-=1; lag2d+=1;} 
      lag2_split=0;
      if (lag2<0 || (int)h1->min+bin1<=abs(lag2d)){ 
	if ((int)h1->min+bin1==abs(lag2d) && lag2>=0){ lag2_split=1;} 
	lag2_under=1; lag2d=(int)h1->min+bin1;}
      lead2 = maximum(0,minimum(h2->nbins-1,bin2)); lead2d=((int)h1->min+bin1) - ((int)h2->min+lead2); lead2_over=0;
      while (lead2<h2->nbins && (int)h2->data[lead2]<=0){ lead2+=1; lead2d-=1;} 
      if (lead2>=h2->nbins){ lead2_over=1;}}
    else if (h2->nbins<=0){ lag2d=(int)h1->min+bin1;lead2d=(int)h1->min+bin1;lag2_under=1;lead2_over=1;lag2_split=0;}
    if (!lag2_under && !lead2_over){ 
      if (abs(lag2d)<abs(lead2d)){ lay2=(int)h2->data[lag2]; lay2d=abs(lag2d);}
      else if (abs(lag2d)>abs(lead2d)){ lay2=(int)h2->data[lead2]; lay2d=abs(lead2d);}
      else if (abs(lag2d)==abs(lead2d)){ lay2=((int)h2->data[lag2]+(int)h2->data[lead2])/2.0; lay2d=abs(lead2d);}}
    else if (lag2_under && !lead2_over){
      if (abs(lag2d)<abs(lead2d)){ 
	if (zerofill2>0){ 
	  if (lag2_split){ lay2=(zerofill2+(int)h2->data[lag2])/2.0; lay2d=abs(lag2d);} 
	  else if (!lag2_split){ lay2=zerofill2; lay2d=abs(lag2d);}} 
	else if (zerofill2<=0){ 
	  if (lag2_split){ lay2=(int)h2->data[lag2]; lay2d=abs(lag2d);}  
	  else if (!lag2_split){ lay2=(int)h2->data[lead2]; lay2d=abs(lead2d);}}}
      else if (abs(lag2d)>abs(lead2d)){ lay2=(int)h2->data[lead2]; lay2d=abs(lead2d);}
      else if (abs(lag2d)==abs(lead2d)){ 
	if (zerofill2>0){ 
	  if (lag2_split){ lay2=(zerofill2+(int)h2->data[lead2])/2.0; lay2d=abs(lead2d);}
	  else if (!lag2_split){ lay2=(zerofill2+(int)h2->data[lead2])/2.0; lay2d=abs(lead2d);}}
	else if (zerofill2<=0){ lay2=(int)h2->data[lead2]; lay2d=abs(lead2d);}}}
    else if (!lag2_under && lead2_over){
      lay2=(int)h2->data[lag2]; lay2d=abs(lag2d);}
    else if (lag2_under && lead2_over){
      if (zerofill2>0){ lay2=zerofill2; lay2d=abs(lag2d);}
      else if (zerofill2<=0){ lay2=0; lay2d=abs(lag2d);}}    
    if (verbose){ printf(" %% ended with lag2=%d,lag2d=%d,lead2=%d,lead2d=%d,lay2=%0.1f,lay2d=%d\n",lag2,lag2d,lead2,lead2d,lay2,lay2d);}
    if (lay1d<lay2d){ classify_correctly=1;}
    else if (lay1d>lay2d){ classify_correctly=0;}
    else if (lay1d==lay2d){
      if (lay1>lay2){ classify_correctly=1;}
      else if (lay1<lay2){ classify_correctly=0;}
      else if (lay1==lay2){ classify_correctly=0.5;}}}
  if (verbose){ printf(" %% finally returning %0.1f\n",classify_correctly);}
  return classify_correctly;
}

double obsdisthistappraise_helper(struct hist *h1,int zerofill1,struct hist *h2,int zerofill2)
{
  /* given observation distribution (zerofill1,h1) and (zerofill2,h2), return the number of correctly classified h1 entries
     (self exclusion is assumed). */
  int verbose=0;
  int nb=0;
  double value=0,valtemp=0;
  if (verbose){ printf(" %% \n");}
  if (verbose){ 
    printf(" %% [entering obsdisthistappraise_helper], zerofill1=%d, zerofill2=%d\n",zerofill1,zerofill2);
    histprintf(h1," %% h1"); histprintf(h2," %% h2");}
  valtemp = obsdisthistappraise_helper_helper(-1,h1,zerofill1-1,h2,zerofill2);
  if (verbose){ printf(" %% correctly_classified==%0.1f at entry -1 (%d observations)\n",valtemp,zerofill1);}
  value += valtemp*zerofill1;
  for (nb=0;nb<h1->nbins;nb++){
    if (h1->data[nb]>=1){
      h1->data[nb]-=1; 
      valtemp = obsdisthistappraise_helper_helper(nb,h1,zerofill1,h2,zerofill2);
      h1->data[nb]+=1;
      if (verbose){ printf(" %% correctly_classified==%0.1f at entry %d (%d observations)\n",valtemp,nb,(int)h1->data[nb]);}
      value += valtemp*(int)h1->data[nb];}}
  if (verbose){ printf(" %% total of %0.1f correctly classified\n",value);}
  return value;
}      

void obsdisthistappraise(struct obsdisthist *odh,int nrecords,int self_exclude)
{
  /* runs through pairs odh->hra[nh1],odh->hra[nh2], determining value (stored at odh->value[nh1+nh2*odh->length]) */
  int verbose=0;
  int nh1=0,nh2=0,nb1=0,nb2=0,v1=0,v2=0;
  double value=0;
  struct hist *h1=NULL,*h2=NULL;
  if (verbose){ printf(" %% [entering obsdisthistappraise] nrecords %d\n",nrecords);}
  odh->nrecords=nrecords;
  for (nh1=0;nh1<odh->length;nh1++){ 
    h1 = odh->hra[nh1]; 
    for (nh2=nh1+1;nh2<odh->length;nh2++){ 
      h2 = odh->hra[nh2];
      if (verbose){ printf(" %% evaluating odh->hra[%d],odh->hra[%d]\n",nh1,nh2); histprintf(h1," %% h1 "); histprintf(h2," %% h2 ");}
      if (!self_exclude){
	if (verbose){ printf(" %% not bothering to worry about self_exclusion\n");}
	nb1=0;nb2=0; value=0;
	while (nb1<h1->nbins && nb2<h2->nbins){
	  v1 = (int)h1->min+nb1; v2 = (int)h2->min+nb2;
	  if (verbose>2){ printf(" nb1 %d v1 %d nb2 %d v2 %d...",nb1,v1,nb2,v2);}
	  if (v1<v2){ if (verbose>2){ printf(" v1<v2, value += %d\n",(int)h1->data[nb1]);} value += (int)h1->data[nb1];nb1++;}
	  else if (v1>v2){ if (verbose>2){ printf(" v1>v2, value += %d\n",(int)h2->data[nb2]);} value += (int)h2->data[nb2]; nb2++;}
	  else if (v1==v2){ 
	    if (verbose>2){ printf(" v1==v2, value += %d\n",maximum((int)h1->data[nb1],(int)h2->data[nb2]));}
	    value += maximum((int)h1->data[nb1],(int)h2->data[nb2]); nb1++; nb2++;}}
	while (nb1<h1->nbins){ if (verbose>2){ printf(" nb1 %d value += %d\n",nb1,(int)h1->data[nb1]);} 
	value += (int)h1->data[nb1]; nb1++;}
	while (nb2<h2->nbins){ if (verbose>2){ printf(" nb2 %d value += %d\n",nb2,(int)h2->data[nb2]);} 
	value += (int)h2->data[nb2]; nb2++;}
	if (verbose>2){ printf(" value += %d\n",(int)maximum(odh->nrecords-h1->ndatum,odh->nrecords-h2->ndatum));}
	value += (int)maximum(odh->nrecords-h1->ndatum,odh->nrecords-h2->ndatum);
	if (verbose){ printf(" %% greedy classification gets %0.1f right...",value);}
	value = maximum(0,(value-odh->nrecords)/(double)odh->nrecords);
	if (verbose){ printf(" value %f\n",value);}}
      else /* if (self_exclude) */{
	if (verbose){ printf(" %% busy with self_exclusion\n");}
	value = 0;
	value += obsdisthistappraise_helper(h1,odh->nrecords-h1->ndatum,h2,odh->nrecords-h2->ndatum);
	value += obsdisthistappraise_helper(h2,odh->nrecords-h2->ndatum,h1,odh->nrecords-h1->ndatum);
	if (verbose){ printf(" %% careful classification gets %0.1f right...",value);}
	value = maximum(0,(value-odh->nrecords)/(double)odh->nrecords);
	if (verbose){ printf(" value %f\n",value);}}
      odh->value[nh1 + nh2*odh->length] = value; odh->value[nh2 + nh1*odh->length] = value;}}
  if (verbose){ printf(" %% [exiting obsdisthistappraise]\n");}
}

void pnode_obs2dist_starter(int ni,int ninputs,int nrecords,struct pnode *parent,struct llitem *l0,int maxlevel,int level,int purpose_flag)
{
  /* purpose_flag equals:
     1: runs through a single observation, creating a histra at each pn->temp if it does not exist, 
     and storing the occurence number within the histra[ninput] (a hist of integers) 
     2: tfrees (struct obsdisthist *) pn->temp
     3: determines value (with appropriate self_exclusion)
     4: print (struct obsdisthist *) pn->temp
  */
  int verbose=0;
  int firstentry=0;
  struct pnode *pn=NULL;
  if (verbose && parent==NULL && l0->parent==NULL){ printf(" %% [entering pnode_obs2dist_starter]  purpose=%d\n",purpose_flag);}
  if (parent!=NULL && l0->parent==NULL){ /* first descent into parent->childllitem */ 
    assert(parent->childllitem==l0); firstentry=1;}
  if (l0->kidl!=NULL){ pnode_obs2dist_starter(ni,ninputs,nrecords,parent,l0->kidl,maxlevel,level,purpose_flag);}
  if (l0->kidr!=NULL){ pnode_obs2dist_starter(ni,ninputs,nrecords,parent,l0->kidr,maxlevel,level,purpose_flag);}
  if (l0->item!=NULL){
    if (maxlevel==-1 || level<maxlevel){ 
      pn=(struct pnode *)l0->item; 
      switch (purpose_flag){
      case 1: if (pn->weight!=0){ pnode_obs2dist_starter(ni,ninputs,nrecords,pn,pn->childllitem,maxlevel,level+1,purpose_flag);} break;
      default: pnode_obs2dist_starter(ni,ninputs,nrecords,pn,pn->childllitem,maxlevel,level+1,purpose_flag); break;}}}
  if (firstentry){ 
    if (parent!=NULL){ 
      switch (purpose_flag){
      case 1: 
	if (parent->weight!=0){
	  if (parent->temp==NULL){ if (verbose){ printf(" %% making temp\n");} parent->temp = obsdisthistmake(ninputs);}
	  obsdisthistupdate((struct obsdisthist *)parent->temp,ni,(int)parent->weight,1);}
	break;
      case 2: 
	if (parent->temp!=NULL){ obsdisthisttfree((struct obsdisthist *)parent->temp); parent->temp=NULL;}
	break;
      case 3: 
	if (parent->temp!=NULL){ obsdisthistappraise((struct obsdisthist *)parent->temp,nrecords,1);}
	break;
      case 4: 
	if (parent->temp!=NULL){ obsdisthistprintf((struct obsdisthist *)parent->temp," %%%% ");}
	break;
      default: printf(" %% warning! invalid purpose_flag=%d in pnode_obs2dist_starter\n",purpose_flag); break;}}}
  if (verbose && parent==NULL && l0->parent==NULL){ printf(" %% [exiting pnode_obs2dist_starter]\n");}
}

void obsdisthistevaluate_helper(int nbins,int max,int min,double *data,int weight,int zerofill,int self_exclude,double *value,double *distance)
{
  /* finds the value of the bin within *h closest to weight, assuming that the zero-bin has zerofill value 
     assumptions are that h is `minimal' in the sense that data[0]!=0 and data[nbins-2]!=0 */
  int verbose=0;
  int nb1=0,nb2=0;
  double d1=0,d2=0;
  if (verbose){ printf(" %% [entering obsdisthistevaluate_helper] nbins %d max %d min %d weight %d zerofill %d self_exclude %d\n",nbins,max,min,weight,zerofill,self_exclude); raprintf(data,"double",1,nbins," %% data: ");}
  if (!self_exclude){
    if (verbose){ printf(" %% no self_exclusion.\n");}
    if (nbins==0){ if (verbose){ printf(" %% %% no data, setting value=zerofill=%d\n",zerofill);} *value=zerofill;*distance=weight;}
    else /* if (nbins>0) */{
      if (verbose){ printf(" %% %% data exists, determining ideal bin number nb1=%d\n",weight-min);}
      nb1 = weight-min;
      if (nb1>nbins-2){ 
	if (verbose){ printf(" %% %% %% nb1>nbins-2, setting value=%f\n",data[nbins-2]);}
	nb1=nbins-2; d1=weight-min-nb1; *value=data[nb1]; *distance=abs(d1);}
      else if (nb1<0){ 
	if (verbose){ printf(" %% %% %% nb1<0\n");}
	if (min-weight<weight-0 || (min-weight>=weight-0 && zerofill==0)){ 
	  if (verbose){ printf(" %% %% %% %% but is closer to first entry in data, setting value=%f\n",data[0]);}
	  nb1=0;d1=weight-min; *value=data[nb1]; *distance=abs(d1);}
	else if (min-weight>weight-0 && zerofill>0){ 
	  if (verbose){ printf(" %% %% %% %% but is closer to zerofill, setting value=zerofill=%d\n",zerofill);}
	  nb1=-1;d1=weight-0; *value=zerofill; *distance=abs(d1);} 
	else if (min-weight==weight-0 && zerofill>0){ 
	  if (verbose){ printf(" %% %% %% %% but is in between zerofill and first entry, setting value=%f\n",0.5*zerofill+0.5*data[0]);}
	  nb1=-1;d1=weight-0; nb2=0;d2=weight-min;
	  *value=0.5*zerofill+0.5*data[nb2]; *distance=0.5*abs(d1)+0.5*abs(d2);}}
      else /* if (nb>=0 && nb<=nbins-2) */{
	if (verbose){ printf(" %% %% %% 0<=(nb1=%d)<=nbins-2\n",nb1);}
	nb2=nb1;d1=0;d2=0;
	while (nb1>0 && data[nb1]==0){ nb1-=1; d1+=1;}
	while (nb2<nbins-2 && data[nb2]==0){ nb2+=1; d2-=1;}
	if (verbose){ printf(" %% %% %% searched and found nb1,data[%d]=%f, nb2,data[%d]=%f\n",nb1,data[nb1],nb1,data[nb2]);}
	if (abs(d1)<abs(d2)){ *value=data[nb1]; *distance=abs(d1);}
	else if (abs(d1)>abs(d2)){ *value=data[nb2]; *distance=abs(d2);}
	else if (abs(d1)==abs(d2)){ *value=0.5*data[nb1]+0.5*data[nb2]; *distance=0.5*abs(d1)+0.5*abs(d2);}
	if (verbose){ printf(" %% %% %% set value=%f",*value);}}}}
  else /* if (self_exclude) */{
    if (verbose){ printf(" %% self exclusion required\n");}
    if (nbins==0){ 
      if (verbose){ printf(" %% %% no data, value=zerofill-1\n");} *value=zerofill-1;*distance=weight;}
    else /* if (nbins>0) */{
      if (verbose){ printf(" %% %% data exists, ideal bin number %d\n",weight-min);}
      nb1 = weight-min;
      if (nb1>nbins-2){ 
	if (verbose){ printf(" %% %% %% nb1>nbins-2, no exclusion necessary\n");}
	obsdisthistevaluate_helper(nbins,max,min,data,weight,zerofill,0,value,distance);}
      else if (nb1<0){ 
	if (verbose){ printf(" %% %% %% nb1<0\n");}
	if (min-weight>=weight-0 && zerofill>0){
	  if (verbose){ printf(" %% %% %% %% but is closer to zerofill, calling with no exclusion\n");}
	  obsdisthistevaluate_helper(nbins,max,min,data,weight,zerofill-1,0,value,distance);}
	else if (min-weight<weight-0 || (min-weight>=weight-0 && zerofill==0)){ 
	  if (verbose){ printf(" %% %% %% %% but is closer to first entry in data...");}
	  if (data[0]>=1){
	    if (verbose){ printf(" calling with no exclusion\n");}
	    data[0]-=1;
	    obsdisthistevaluate_helper(nbins,max,min,data,weight,zerofill,0,value,distance);
	    data[0]+=1;}
	  else if (data[0]<=1){
	    if (verbose){ printf(" need to exclude first entry\n");}
	    nb1=1; d1=1;
	    while (nb1<nbins-2 && data[nb1]==0){ nb1+=1;d1+=1;}
	    if (data[nb1]>0){
	      if (verbose){ printf(" found nearby entry nb1=%d, calling with no exclusion\n",nb1);}
	      obsdisthistevaluate_helper(nbins-d1,max,min+d1,&(data[nb1]),weight,zerofill,0,value,distance);}
	    else /* if (data[nb1]<=0) */{
	      if (verbose){ printf(" could not find nearby entry, calling with no exclusion\n");}
	      obsdisthistevaluate_helper(0,0,0,NULL,weight,maximum(0,zerofill-1),0,value,distance);}}}}
      else /* if (nb1>=0 && nb1<=nbins-2) */{
	if (verbose){ printf(" %% %% %% 0<=(nb1=%d)<=nbins-2\n",nb1);}
	if (nb1==0 && data[0]<=1){ 
	  if (verbose){ printf(" %% %% %% %% however, data[%d]==%f, calling with no exclusion\n",nb1,data[0]);}
	  data[nb1]-=1;
	  nb2=1; while (nb2<nbins-1 && data[nb2]<=0){ nb2+=1;} 
	  if (nb2<nbins-1){ obsdisthistevaluate_helper(nbins-nb2,max,min+nb2,&(data[nb2]),weight,zerofill,0,value,distance);}
	  else /* if (nb2>=nbins-1) */{ obsdisthistevaluate_helper(0,0,0,NULL,weight,zerofill,0,value,distance);}
	  data[nb1]+=1;}
	else if (nb1==nbins-2 && data[nb1]<=1){
	  if (verbose){ printf(" %% %% %% %% however, data[%d]==%f, calling with no exclusion\n",nb1,data[nbins-2]);}
	  data[nb1]-=1;
	  nb2 = nbins-3; while (nb2>=0 && data[nb2]<=0){ nb2-=1;}
	  if (nb2>=0){ obsdisthistevaluate_helper(nb2+2,min+nb2+2,min,&(data[0]),weight,zerofill,0,value,distance);}
	  else /* if (nb2<0) */{ obsdisthistevaluate_helper(0,0,0,NULL,weight,zerofill,0,value,distance);}
	  data[nb1]+=1;}
	else /* if ((nb1>0 && nb1<nbins-2) || data[nb1]>1) */{
	  if (verbose){ printf(" %% %% %% %% data[%d]==%f, calling with no exclusion\n",nb1,data[nb1]);}
	  data[nb1]-=1;
	  obsdisthistevaluate_helper(nbins,max,min,&(data[0]),weight,zerofill,0,value,distance);
	  data[nb1]+=1;}}}}
}

void obsdisthistevaluate(struct obsdisthist *odh,int weight,int ninput)
{
  /* given a sample weight from ninput, 
     determine (after self exclusion) which of the odh->hra[ninput],odh->hra[nh] is the most likely generating distribution.
     results are stored in odh->ballot_2, +1 means right, -1 means wrong 
     results are stored in odh->ballot_3, most likely index recorded (-1 means tie) */
  int verbose=0;
  int nh=0,nh1=0,nh2=0,tab1=0,tab2=0,test_12=0,test_23=0,test_13=0;
  double v1=0,v2=0,v3=0,d1=0,d2=0,d3=0;
  struct hist *h1=NULL,*h2=NULL,*h3=NULL;
  if (verbose){ printf(" %% [entering obsdisthistevaluate] weight %d ninput %d\n",weight,ninput);}
  if (verbose){ printf(" %% making ballot_2\n");}
  h1=odh->hra[ninput];
  for (nh=0;nh<odh->length;nh++){
    if (nh==ninput){ odh->ballot_2[nh]=0;}
    else /* if (nh!=ninput) */{
      h2 = odh->hra[nh];
      if (verbose){ printf(" %% considering nh %d:\n",nh);}
      if (verbose>2){ histprintf(h1," %% h1"); histprintf(h2," %% h2");}
      obsdisthistevaluate_helper(h1->nbins,h1->max,h1->min,h1->data,weight,odh->nrecords-h1->ndatum,1,&v1,&d1);
      obsdisthistevaluate_helper(h2->nbins,h2->max,h2->min,h2->data,weight,odh->nrecords-h2->ndatum,0,&v2,&d2);
      if (d1<d2 || (d1==d2 && v1>v2)){ odh->ballot_2[nh]=odh->value[nh+ninput*odh->length];}
      else if (d1>d2 || (d1==d2 && v1<v2)){ odh->ballot_2[nh]=-odh->value[nh+ninput*odh->length];}
      else if (d1==d2 && v1==v2){ odh->ballot_2[nh]=0;}
      if (verbose){ printf(" %% found d1 %0.1f v1 %0.1f d2 %0.1f v2 %0.1f ballot %f\n",d1,v1,d2,v2,odh->ballot_2[nh]);}}}
  if (verbose){ raprintf(odh->ballot_2,"double",1,odh->length," %% ");}
  if (verbose){ printf(" %% making ballot_3\n");}
  h1=odh->hra[ninput];
  for (nh1=0;nh1<odh->length;nh1++){
    for (nh2=0;nh2<nh1;nh2++){
      tab1 = nh1 + nh2*odh->length; tab2 = nh2 + nh1*odh->length;
      if (nh1==ninput || nh2==ninput){ odh->ballot_3[tab1]=-1; odh->ballot_3[tab2]=-1;}
      else /* if (nh1!=ninput && nh2!=ninput) */{
	h2 = odh->hra[nh1];
	h3 = odh->hra[nh2];
	if (verbose){ printf(" %% considering nh1 %d, nh2 %d\n",nh1,nh2);}
	obsdisthistevaluate_helper(h1->nbins,h1->max,h1->min,h1->data,weight,odh->nrecords-h1->ndatum,1,&v1,&d1);
	obsdisthistevaluate_helper(h2->nbins,h2->max,h2->min,h2->data,weight,odh->nrecords-h2->ndatum,0,&v2,&d2);
	obsdisthistevaluate_helper(h3->nbins,h3->max,h3->min,h3->data,weight,odh->nrecords-h3->ndatum,0,&v3,&d3);
	test_12=0;test_13=0;test_23=0;
	if (d1<d2 || (d1==d2 && v1>v2)){ test_12=1;}
	else if (d1>d2 || (d1==d2 && v1<v2)){ test_12=-1;}
	else if (d1==d2 && v1==v2){ test_12=0;}
	if (d1<d3 || (d1==d3 && v1>v3)){ test_13=1;}
	else if (d1>d3 || (d1==d3 && v1<v3)){ test_13=-1;}
	else if (d1==d3 && v1==v3){ test_13=0;}
	if (d2<d3 || (d2==d3 && v2>v3)){ test_23=1;}
	else if (d2>d3 || (d2==d3 && v2<v3)){ test_23=-1;}
	else if (d2==d3 && v2==v3){ test_23=0;}
	if (verbose){ printf(" %% found d1 %0.1f v1 %0.1f d2 %0.1f v2 %0.1f d3 %0.1f v3 %0.1f\n",d1,v1,d2,v2,d3,v3);}
	if (verbose){ printf(" %% test_12 %d test_13 %d test_23 %d\n",test_12,test_13,test_23);}
	if (test_12==1 && test_13==1){ odh->ballot_3[tab1]=ninput; odh->ballot_3[tab2]=ninput;}
	else if (test_12==-1 && test_23==1){ odh->ballot_3[tab1]=nh1; odh->ballot_3[tab2]=nh1;}
	else if (test_13==-1 && test_23==-1){ odh->ballot_3[tab1]=nh2; odh->ballot_3[tab2]=nh2;}
	else /* cycle */{ odh->ballot_3[tab1]=-1; odh->ballot_3[tab2]=-1;}
	if (verbose){ printf(" %% odh->ballot_3[%d] = odh->ballot_3[%d] = %d\n",tab1,tab2,(int)odh->ballot_3[tab1]);}}}}
  if (verbose){ printf(" %% [exiting obsdisthistevaluate]\n");}
}

void ptree_classify_helper1(struct ptree *p0,struct pnode *parent0,struct llitem *l0,struct ptree *p1,struct pnode *parent1,struct llitem *l1,int ninput,int maxlevel,int level)
{
  /* uses value within pn0->temp to assess pairwise discriminability of observation p1 of input ninput, draw a ballot */
  int verbose=0,nr=0;
  char lvlchar[64];
  int firstentry=0;
  struct pnode *pn0=NULL,*pn1=NULL;
  struct llitem *l2=NULL;
  int weight_sample=0;
  if (verbose){ sprintf(lvlchar," %%"); for (nr=0;nr<level;nr++){sprintf(lvlchar,"%s %%",lvlchar);}}
  if (verbose>1){ 
    printf(" %s [entering ptree_classify_helper1]\n",lvlchar);
    printf(" %s p0 %d parent0 %d->%d l0 %d \n",lvlchar,(int)p0,(int)parent0,parent0==NULL ? 0 : (int)(parent0->childllitem),(int)l0);
    printf(" %s p1 %d parent1 %d->%d l1 %d\n",lvlchar,(int)p1,(int)parent1,parent1==NULL ? 0 : (int)(parent1->childllitem),(int)l1);}
  if (parent1!=NULL){ assert(parent1->childllitem==l1);}
  if (parent0!=NULL && l0->parent==NULL){ /* first descent into parent0->childllitem */ 
    if (verbose>1){ printf(" %s first descent into parent0->childllitem, setting firstentry to 1\n",lvlchar);}
    assert(parent0->childllitem==l0); firstentry=1;}
  if (l0->kidl!=NULL){ 
    if (verbose>1){ printf(" %s l0->kidl found, moving left to new llitem l0->kidl=%d\n",lvlchar,(int)(l0->kidl));}
    ptree_classify_helper1(p0,parent0,l0->kidl,p1,parent1,l1,ninput,maxlevel,level);}
  if (l0->kidr!=NULL){ 
    if (verbose>1){ printf(" %s l0->kidr found, moving right to new llitem l0->kidr=%d\n",lvlchar,(int)(l0->kidr));}
    ptree_classify_helper1(p0,parent0,l0->kidr,p1,parent1,l1,ninput,maxlevel,level);}
  if (l0->item!=NULL){ 
    if (verbose>1){ printf(" %s l0->item %d exists\n",lvlchar,(int)(struct pnode *)(l0->item));}
    pn0 = (struct pnode *)l0->item; 
    if (maxlevel==-1 || level<maxlevel){
      if (verbose>1){ printf(" %s level %d less than maxlevel %d, continuing...\n",lvlchar,level,maxlevel);}
      if (l1!=NULL && (l2=llitemaddorfind(0,l1,p1->regionra[pn0->region->label],&region2pnode_compare_label))!=NULL){
	pn1 = (struct pnode *)l2->item;
	if (verbose>1){ printf(" %s pn0 label %d matches pn1 label %d, descending...\n",lvlchar,pn0->region->label,pn1->region->label);}
	ptree_classify_helper1(p0,pn0,pn0->childllitem,p1,pn1,pn1->childllitem,ninput,maxlevel,level+1);}
      else{
	pn1 = NULL;
	if (verbose>1){ printf(" %s pn0 label %d not found within l1, descending...\n",lvlchar,pn0->region->label);}
	ptree_classify_helper1(p0,pn0,pn0->childllitem,p1,NULL,NULL,ninput,maxlevel,level+1);}}}
  if (firstentry){ 
    weight_sample = (parent1==NULL ? 0 : (int)parent1->weight);
    if (verbose>0){ printf(" %s first entry, and observed weight %d\n",lvlchar,weight_sample);}
    if (parent0->temp!=NULL){ obsdisthistevaluate((struct obsdisthist *)parent0->temp,weight_sample,ninput);}}
}

void ptree_classify_helper2(struct ptree *p0,struct pnode *parent0,struct llitem *l0,int ninput,int maxlevel,int level,double *b2box,double *b3box)
{
  /* once ballot has been cast, we tally the votes
     assume *b2box is organized as [ninput2 + level*odh->length]
     assume *b3box is organized as [123 + ((level*odh->length + ninput3)*odh->length + ninput2)*3], ninput3<ninput2 */
  int verbose=0,nr=0;
  int nh=0,nl=0,nh1=0,nh2=0,tab1=0,tab2=0,tab1l=0,nb=0;
  double valmin=0;
  char lvlchar[64];
  int firstentry=0;
  struct pnode *pn=NULL;
  struct obsdisthist *odh=NULL;
  if (verbose){ sprintf(lvlchar," %%"); for (nr=0;nr<level;nr++){sprintf(lvlchar,"%s %%",lvlchar);}}
  if (verbose>1){ 
    printf(" %s [entering ptree_classify_helper2]\n",lvlchar);
    printf(" %s p0 %d parent0 %d->%d l0 %d \n",lvlchar,(int)p0,(int)parent0,parent0==NULL ? 0 : (int)(parent0->childllitem),(int)l0);}
  if (parent0!=NULL && l0->parent==NULL){ /* first descent into parent0->childllitem */ 
    assert(parent0->childllitem==l0); firstentry=1;}
  if (l0->kidl!=NULL){ ptree_classify_helper2(p0,parent0,l0->kidl,ninput,maxlevel,level,b2box,b3box);}
  if (l0->kidr!=NULL){ ptree_classify_helper2(p0,parent0,l0->kidr,ninput,maxlevel,level,b2box,b3box);}
  if (l0->item!=NULL){
    if (maxlevel==-1 || level<maxlevel){ 
      pn=(struct pnode *)l0->item; 
      ptree_classify_helper2(p0,pn,pn->childllitem,ninput,maxlevel,level+1,b2box,b3box);}}
  if (firstentry){ 
    if (parent0->temp!=NULL){
      odh = (struct obsdisthist *)parent0->temp;
      for (nh=0;nh<odh->length;nh++){ for (nl=level-1;nl<maxlevel;nl++){
	if (verbose){ printf(" %% ninput %d nh %d nl %d/%d ballot %f",ninput,nh,nl,maxlevel,odh->ballot_2[nh]);}
	b2box[nh + nl*odh->length] += odh->ballot_2[nh];
	if (verbose){ printf(" sums to %f\n",b2box[nh + nl*odh->length]);}}}
      for (nh1=0;nh1<odh->length;nh1++){ for (nh2=0;nh2<nh1;nh2++){ 
	if ((int)odh->ballot_3[nh1 + nh2*odh->length]==ninput){ 
	  tab1 = ninput + nh1*odh->length;
	  tab2 = ninput + nh2*odh->length;
	  nb = 0;
	  valmin = minimum(fabs(odh->value[tab1]),fabs(odh->value[tab2]));}
	else if ((int)odh->ballot_3[nh1 + nh2*odh->length]==nh1){
	  tab1 = ninput + nh1*odh->length;
	  tab2 = nh1 + nh2*odh->length;
	  nb = 1;
	  valmin = minimum(fabs(odh->value[tab1]),fabs(odh->value[tab2]));}
	else if ((int)odh->ballot_3[nh1 + nh2*odh->length]==nh2){
	  tab1 = ninput + nh2*odh->length;
	  tab2 = nh1 + nh2*odh->length;
	  nb = 2;
	  valmin = minimum(fabs(odh->value[tab1]),fabs(odh->value[tab2]));}
	else if ((int)odh->ballot_3[nh1 + nh2*odh->length]==-1){
	  nb = -1;
	  valmin = 0;}
	else /* warning */{ printf(" %% warning, incorrect ballot_3[%d]=%f in ptree_classify_helper2\n",nh1+nh2*odh->length,odh->ballot_3[nh1+nh2*odh->length]);}
	if (nb>=0){
	  for (nl=level-1;nl<maxlevel;nl++){
	    tab1l = ((nl*odh->length+nh2)*odh->length+nh1)*3;
	    if (verbose){ printf(" %% ninput %d nh1 %d nh2 %d nl %d/%d ballot_3 %0.1f, adding %0.1f to b3box[%d]\n",ninput,nh1,nh2,nl,maxlevel,odh->ballot_3[nh1+nh2*odh->length],valmin,nb+tab1l);}
	    b3box[nb+tab1l] += valmin;}}}}}}
}

void ptree_classify_helper3(char *filename_output,int ninputs,int maxlevel,struct ptree *p0,double *victorra_2,double *victorra_3)
{
  /* fprintfs relevant discriminability output 
   assumes victorra_2 organized as [nh1 + nh2*ninputs + level*ninputs*ninputs] 
   assumes victorra_3 organized as [nh1 + nh2*ninputs + nh2*ninputs*ninputs + level*ninputs*ninputs*ninputs] */
  FILE *fp=NULL;
  int level=0,nr1=0,nr2=0,nr3=0;
  double *tempra=NULL,*tempra2=NULL,sum=0;
  if ((fp=fopen(filename_output,"w"))==NULL){ printf("warning, cannot read %s in ptree_classify_helper3\n",filename_output); fp=stdout;}
  fprintf(fp," %% ninputs %d maxlevel %d\n",ninputs,maxlevel);
  for (level=0;level<maxlevel;level++){ 
    fprintf(fp," %% level %d victorra_2\n",level);
    for (nr1=0;nr1<ninputs;nr1++){ 
      fprintf(fp," %% ");
      for (nr2=0;nr2<ninputs;nr2++){ 
	fprintf(fp," %0.1f ",victorra_2[nr1+nr2*ninputs+level*ninputs*ninputs]);} 
      fprintf(fp,"\n");}}
  for (level=0;level<maxlevel;level++){ 
    fprintf(fp," %% level %d victorra_3\n",level);
    for (nr1=0;nr1<ninputs;nr1++){ 
      fprintf(fp," %% %% nr1 %d\n",nr1);
      for (nr2=0;nr2<ninputs;nr2++){
	fprintf(fp," %% %% ");
	for (nr3=0;nr3<ninputs;nr3++){ 
	  if (/* nr3<nr2 && nr2!=nr1 && nr3!=nr1 */1){ 
	    fprintf(fp," %0.1f ",victorra_3[nr1 + nr2*ninputs + nr3*ninputs*ninputs + level*ninputs*ninputs*ninputs]);} 
	  else /* non-cycle */ { fprintf(fp," XX ");}}
	fprintf(fp,"\n");}}}
  for (level=0;level<maxlevel;level++){
    fprintf(fp," %% level %d obsdisthist->value statistics\n",level);
    tempra = (double *) tcalloc(2,sizeof(double));
    tempra2 = (double *) tcalloc(4,sizeof(double));
    for (nr1=0;nr1<ninputs;nr1++){
      for (nr2=0;nr2<ninputs;nr2++){
	tempra[0]=nr1;tempra[1]=nr2;
	pnode_statstemp_starter(NULL,p0->postree,level+1,-1,0,NULL,&(tempra2[0]),&(tempra2[1]),&(tempra2[2]),&(tempra2[3]),&pnode_obsdisthist_value,tempra);
	fprintf(fp," %% %% inputs %d vs %d : max %f,min %f,mean %f,stdev %f\n",nr1,nr2,tempra2[0],tempra2[1],tempra2[2],tempra2[3]);}}
    tfree(tempra);tempra=NULL;
    tfree(tempra2);tempra2=NULL;}
  for (level=0;level<maxlevel;level++){ 
    fprintf(fp," %% level %d victorra_2: total\n",level);
    for (nr1=0;nr1<ninputs;nr1++){ for (nr2=0;nr2<nr1;nr2++){
      sum = victorra_2[nr1+nr2*ninputs+level*ninputs*ninputs]+victorra_2[nr2+nr1*ninputs+level*ninputs*ninputs];
      fprintf(fp," (nr1 %d,nr2 %d, %0.1f) ",nr1,nr2,sum);}}
    fprintf(fp,"\n");}
  for (level=0;level<maxlevel;level++){ 
    fprintf(fp," %% level %d victorra_3: total\n",level);
    for (nr1=0;nr1<ninputs;nr1++){ for (nr2=0;nr2<nr1;nr2++){ for (nr3=0;nr3<nr2;nr3++){
      sum = 0;
      sum += victorra_3[nr1+nr2*ninputs+nr3*ninputs*ninputs+level*ninputs*ninputs*ninputs];
      sum += victorra_3[nr2+nr1*ninputs+nr3*ninputs*ninputs+level*ninputs*ninputs*ninputs];
      sum += victorra_3[nr3+nr1*ninputs+nr2*ninputs*ninputs+level*ninputs*ninputs*ninputs];
      fprintf(fp," (nr1 %d,nr2 %d,nr3 %d, %0.1f) ",nr1,nr2,nr3,sum);}}}
    fprintf(fp,"\n");}
  if (fp!=stdout){ fclose(fp);}
}

void ptree_classify_starter(int ninputs,int maxlevel,int verbose,char **filename_base,char *outputname,int bitbybit_remove)
{
  /* assumes that outputname does NOT start with "./" 
     also assumes that the ordering of **filename_base corresponds to odh->hra ordering */
  int old_region_type=GLOBAL_PTREE_REGION_TYPE;
  char *gs2=GLOBAL_STRING_2,filename_output[256],filename[256],command[512];
  struct ptree *ptmp=NULL,*p0=GLOBAL_PTREE;
  int bitbybit=0,continue_flag=0,level=0;
  int nh=0,nR=0,nr=0,nh1=0,nh2=0,tab1=0,tab2=0;
  double *b2box=NULL,*b3box=NULL,*victorra_2=NULL,*victorra_3=NULL;
  GLOBAL_PTREE_REGION_TYPE=0;
  if (verbose>1){ 
    printf(" %% [entering ptree_classify_starter]\n");
    printf(" %% ninputs %d maxlevel %d verbose %d remove? %d\n",ninputs,maxlevel,verbose,bitbybit_remove);
    for (nh=0;nh<ninputs;nh++){ printf(" %% filename_base[%d]=%s\n",nh,filename_base[nh]);}
    printf(" %% outputname=%s\n",outputname);}
  pnodeclear_starter(NULL,p0->postree); pnodeclear_starter(NULL,p0->pretree); p0->total_time=0;
  b2box = (double *) tcalloc(ninputs*maxlevel,sizeof(double));
  b3box = (double *) tcalloc(3*ninputs*ninputs*maxlevel,sizeof(double));
  victorra_2 = (double *) tcalloc(ninputs*ninputs*maxlevel,sizeof(double));
  victorra_3 = (double *) tcalloc(ninputs*ninputs*ninputs*maxlevel,sizeof(double));
  for (nh=0;nh<ninputs;nh++){
    nR = 0;
    bitbybit=1; continue_flag=1;
    do{
      sprintf(filename,"./%s_%d",filename_base[nh],bitbybit);
      if (verbose>1){ printf(" %% checking %s... ",filename);}
      if (checktofind(filename)){ 
	if (verbose>1){ printf("found, reading... ");}
	ptmp = ptreadback(filename,0);
	if (bitbybit_remove){ 
	  if (verbose>1){ printf(" removing observation %s\n",filename);}
	  sprintf(command,"""rm"" %s ;",filename); system(command);}
	else /* if (!bitbybit_remove) */{
	  if (verbose>1){ printf(" retaining\n");}}
	bitbybit += 1; continue_flag=1;}
      else /* if not found */{ if (verbose>1){ printf("not found, stopping\n");} continue_flag=0;}
      if (continue_flag){
	nR += 1;
	ptree_classify_helper1(p0,NULL,p0->postree,ptmp,NULL,ptmp->postree,nh,maxlevel,0);
	for (nr=0;nr<ninputs*maxlevel;nr++){ b2box[nr]=0;}
	for (nr=0;nr<3*ninputs*ninputs*maxlevel;nr++){ b3box[nr]=0;}
	ptree_classify_helper2(p0,NULL,p0->postree,nh,maxlevel,0,b2box,b3box);
	for (level=0;level<maxlevel;level++){ for (nr=0;nr<ninputs;nr++){ 
	  victorra_2[nh+nr*ninputs+level*ninputs*ninputs]+=(b2box[nr+level*ninputs]>0);}}
	if (verbose>2){ 
	  for (level=0;level<maxlevel;level++){ 
	    printf(" %% level %d b2box\n",level);
	    raprintf(&(b2box[0+level*ninputs]),"double",1,ninputs," %% %% ");
	    printf(" %% level %d victorra_2\n",level);
	    raprintf(&(victorra_2[0+0*ninputs+level*ninputs*ninputs]),"double",ninputs,ninputs," %% %% ");}}
	for (nh1=0;nh1<ninputs;nh1++){ for (nh2=0;nh2<nh1;nh2++){ if (nh1!=nh && nh2!=nh){ for (level=0;level<maxlevel;level++){ 
	  tab1 = (((level*ninputs + nh2)*ninputs + nh1)*ninputs + nh);
	  tab2 = (((level*ninputs + nh2)*ninputs + nh1)*3 + 0);
	  victorra_3[tab1] += (b3box[0+tab2]>b3box[1+tab2] && b3box[0+tab2]>b3box[2+tab2]);}}}}
	if (verbose>2){
	  for (level=0;level<maxlevel;level++){ 
	    printf(" %% level %d b3box\n",level);
	    for (nr=0;nr<ninputs;nr++){
	      printf(" %% %% nh2 %d\n",nr);
	      raprintf(&(b3box[0+0+nr*3*ninputs+level*3*ninputs*ninputs]),"double",3,ninputs," %% %% %%");}
	    printf(" %% level %d victorra_3\n",level);
	    for (nr=0;nr<ninputs;nr++){
	      printf(" %% %% nh2 %d\n",nr);
	      raprintf(&(victorra_3[0+0+nr*ninputs*ninputs+level*ninputs*ninputs*ninputs]),"double",ninputs,ninputs," %% %% %%");}}}}
      if (ptmp!=NULL){ ptreetfree(ptmp);ptmp=NULL;}}
    while (continue_flag);}
  if (outputname!=NULL){ sprintf(filename_output,"ptree_classify_%s_%srecord.txt",outputname,gs2);}
  else /* if (outputname==NULL) */{ sprintf(filename_output,"ptree_classify_%srecord.txt",gs2);}
  ptree_classify_helper3(filename_output,ninputs,maxlevel,p0,victorra_2,victorra_3);
  GLOBAL_PTREE_REGION_TYPE=old_region_type;
  tfree(b2box);b2box=NULL;
  tfree(b3box);b3box=NULL;
  tfree(victorra_2);victorra_2=NULL;
  tfree(victorra_3);victorra_3=NULL;
}

void ptree_dumptemp_starter(struct ptree *p,char *fgvn,int dump_type,void (*temp2file)(void *,FILE *),int keepzeroweights)
{
  /* dumps ptree *p as a sparse array, first dumping
     p->nregions
     p->nlegs
     p->legtime
     p->length 
     p->update_every
     p->total_time
     then calling ptree_dumptemp_helper;
     assumes fgvn does NOT start with "./"
  */
  char filename[256];
  FILE *fp=NULL;
  struct llist *L=NULL;
  if (fgvn==NULL){ sprintf(filename,"./ptree_dumptemp_%srecord",GLOBAL_STRING_2);}
  else /* if (fgvn!=NULL) */{ sprintf(filename,"./ptree_%s_%srecord",fgvn,GLOBAL_STRING_2);}
  if ((fp=fopen(filename,"w"))==NULL){ printf("warning, cannot create %s in ptree_dumptemp_starter",filename); fp=stdout;}
  fwrite(&(p->nregions),sizeof(int),1,fp);
  fwrite(&(p->nlegs),sizeof(int),1,fp);
  fwrite(&(p->legtime),sizeof(int),1,fp);
  fwrite(&(p->length),sizeof(int),1,fp);
  fwrite(&(p->update_every),sizeof(double),1,fp);
  fwrite(&(p->total_time),sizeof(double),1,fp);
  L=llistmake();
  ptree_dumptemp_helper(1,NULL,p->postree,fp,L,temp2file,keepzeroweights);
  llisttfree(L); L=NULL;
  L=llistmake();
  ptree_dumptemp_helper(0,NULL,p->pretree,fp,L,temp2file,keepzeroweights);
  llisttfree(L); L=NULL;
  if (fp!=stdout){ fclose(fp);}
}

void ptree_dumptemp_helper(int posorpre,struct pnode *parent,struct llitem *l0,FILE *fp,struct llist *L,void (*temp2file)(void *,FILE *),int keepzeroweights)
{
  /* runs through llitem *l0 of pnode items, the format is:
     parent->weight,parent->relevance,parent->temp,label_0,label_1,...,label_n,posorpre ? -1 : -2
     any pnode order works, but we do reachfirst.
  */
  int firstentry=0;
  int minusone=-1,minustwo=-2;
  struct pnode *pn=NULL;
  struct litem *l=NULL;
  struct region *r=NULL;
  if (parent!=NULL && l0->parent==NULL){ /* first descent into parent->childllitem */ 
    assert(parent->childllitem==l0); firstentry=1;}
  if (firstentry){ litemadd(L,parent->region);}
  if (firstentry){ 
    fwrite(&(parent->weight),sizeof(double),1,fp);
    fwrite(&(parent->relevance),sizeof(double),1,fp);
    temp2file(parent->temp,fp);
    l = L->first;
    while (l!=NULL){
      r = (struct region *) l->item;
      fwrite(&(r->label),sizeof(int),1,fp);
      l=l->child;}
    if (posorpre){ fwrite(&minusone,sizeof(int),1,fp);}
    else{ fwrite(&minustwo,sizeof(int),1,fp);}}
  if (l0->kidl!=NULL){ ptree_dumptemp_helper(posorpre,parent,l0->kidl,fp,L,temp2file,keepzeroweights);}
  if (l0->item!=NULL){ 
    pn=(struct pnode *)l0->item; 
    if (keepzeroweights || (!keepzeroweights && pn->weight!=0)){ 
      ptree_dumptemp_helper(posorpre,pn,pn->childllitem,fp,L,temp2file,keepzeroweights);}}
  if (l0->kidr!=NULL){ ptree_dumptemp_helper(posorpre,parent,l0->kidr,fp,L,temp2file,keepzeroweights);}
  if (firstentry){ llistkilllast(L);}
}

struct ptree * ptree_readbacktemp(char *filename,void *(*file2temp)(FILE *,int *),int *no_error_passback)
{
  int verbose=0;
  FILE *fp=NULL;
  struct ptree *p=NULL;
  int nr=0;
  struct pnode *pn=NULL;
  struct llitem *l0=NULL,*l1=NULL;
  int readout=0;
  double weight=0,relevance=0;
  void *void_temp=NULL;
  struct llist *L=NULL;
  struct litem *l=NULL;
  int *label=NULL;
  int event_threshold = GLOBAL_PTREE_EVENT_THRESHOLD;
  double event_within = GLOBAL_PTREE_EVENT_WITHIN;
  int region_type = GLOBAL_PTREE_REGION_TYPE;
  int exit_flag=0,no_error=0,no_minor_error=0,depth=0;
  if (verbose){ printf(" %% \n");}
  if (verbose){ printf(" %% [entering ptree_readbacktemp] with filename %s\n",filename==NULL ? "<null>" : filename);}
  if ((fp=fopen(filename,"r"))==NULL){ printf("warning, cannot read %s in ptree_readbacktemp\n",filename);}
  else /* if (fp!=NULL) */{
    p = (struct ptree *) tcalloc(1,sizeof(struct ptree));
    no_error=1;
    if (no_error){ no_error = (fread(&(p->nregions),sizeof(int),1,fp)==(size_t)1) && !feof(fp) && !ferror(fp);}
    if (no_error){ no_error = (fread(&(p->nlegs),sizeof(int),1,fp)==(size_t)1) && !feof(fp) && !ferror(fp);}
    if (no_error){ no_error = (fread(&(p->legtime),sizeof(int),1,fp)==(size_t)1) && !feof(fp) && !ferror(fp);}
    if (no_error){ no_error = (fread(&(p->length),sizeof(int),1,fp)==(size_t)1) && !feof(fp) && !ferror(fp);}
    if (no_error){ no_error = (fread(&(p->update_every),sizeof(double),1,fp)==(size_t)1) && !feof(fp) && !ferror(fp);}
    if (no_error){ no_error = (fread(&(p->total_time),sizeof(double),1,fp)==(size_t)1) && !feof(fp) && !ferror(fp);}
    if (verbose){ printf(" %% no_error %d nregions %d, nlegs %d, legtime %d, length %d, update_every %f, total_time %f\n",no_error,p->nregions,p->nlegs,p->legtime,p->length,p->update_every,p->total_time);}
    if (no_error){
    p->regionra = (struct region **) tcalloc(p->nregions,sizeof(struct region*));
    if (verbose){ printf(" %% now making regions\n");}
    regionramake(p,event_within,event_threshold,region_type);
    assert(p->length>p->nlegs*p->legtime);
    p->eventra = (struct llitem **) tmalloc(p->length*sizeof(struct llitem *)); for (nr=0;nr<p->length;nr++){ p->eventra[nr]=llitemmake();}
    p->tab=0;
    p->update_last=GLOBAL_TI;
    if (verbose){ printf(" %% now making pretree and postree\n");}
    p->pretree=llitemmake();
    p->postree=llitemmake();
    p->wh=histmake(1,+1,-1);
    p->rh=histmake(1,+1,-1);
    exit_flag=0;
    while(no_error && exit_flag!=3){
      weight=0;relevance=0;no_minor_error=1;
      if (no_minor_error){ no_minor_error = (fread(&weight,sizeof(double),1,fp)==(size_t)1) && !feof(fp) && !ferror(fp);}
      if (no_minor_error){ no_minor_error = (fread(&relevance,sizeof(double),1,fp)==(size_t)1) && !feof(fp) && !ferror(fp);}
      if (no_minor_error){
	if (verbose>1){ printf(" %% no_minor_error %d, read weight %f, relevance %f\n",no_minor_error,weight,relevance);}
	void_temp = file2temp(fp,&no_error);
	L=llistmake();exit_flag=0;
	while (no_minor_error && no_error && !exit_flag){
	  label = (int *) tmalloc(sizeof(int));
	  readout = fread(label,sizeof(int),1,fp);
	  no_minor_error = no_minor_error && readout==1 && !feof(fp) && !ferror(fp);
	  if (readout!=1){ tfree(label); label=NULL; exit_flag=3;}
	  else{ 
	    if (*label>=0){ litemadd(L,label); exit_flag=0;}
	    else if (*label==-1){ tfree(label); label=NULL; exit_flag=1;}
	    else if (*label==-2){ tfree(label); label=NULL; exit_flag=2;}
	    else{ printf("warning, funny terminator %d in ptree_readbacktemp\n",*label);}}}
	if (no_minor_error && no_error && L->length>0 && (exit_flag==1 || exit_flag==2)){ /* postree or pretree */ 
	  l=L->first;depth=1; if (exit_flag==1){ l0=p->postree;pn=NULL;} else{ l0=p->pretree;pn=NULL;}
	  while(l!=NULL && depth<L->length){ /* traverse through nonlast elements of llist *L */
	    label = (int *) l->item;
	    if ((l1=llitemaddorfind(0,l0,p->regionra[*label],&region2pnode_compare_label))==NULL){
	      l1 = llitemaddorfind(1,l0,pnodemake(pn,p->regionra[*label],0,0),&pnode2pnode_compare_label);}
	    pn = (struct pnode *) l1->item;
	    l0 = pn->childllitem;
	    l=l->child;
	    depth+=1;}
	  if (depth==L->length){ /* at end of llist *L, must add pnode */
	    label = (int *) l->item;
	    if ((l1=llitemaddorfind(0,l0,p->regionra[*label],&region2pnode_compare_label))==NULL){
	      l1 = llitemaddorfind(1,l0,pnodemake(pn,p->regionra[*label],weight,relevance),&pnode2pnode_compare_label);}
	    else{ pn = (struct pnode *)l1->item; pn->weight += weight; pn->relevance=relevance;}
	    pn = (struct pnode *)l1->item; pn->temp = void_temp;}}
	llisttfree3(L);L=NULL;}
      else /* if (!no_minor_error) */{ if (verbose){ printf(" %% couldn't read weight and relevance\n");} exit_flag=3;}}}
    if (fp!=stdout){ fclose(fp);}
    if (no_error){
      if (verbose){ printf(" %% now balancing pretree and postree\n");}
      llitembalance(p->pretree); p->pretree = llitemclimb(p->pretree);
      llitembalance(p->postree); p->postree = llitemclimb(p->postree);}
    else /* if (!no_error) */{ 
      if (verbose){ printf(" %% major error, freeing ptree, without freeing any p->temp -- memory leak!\n");}
      ptreetfree(p);p=NULL;}}
  if (no_error_passback!=NULL){ *no_error_passback=no_error;}
  return p;
}

void hist_tofig(void *void_hist,FILE *fp,void *void_parameters,double side,double xoffset,double yoffset)
{
  /* fits histogram in a box centered on (xoffset,yoffset) */
  /* parameter list:
     double *colorcode 
     double *xmax
     double *xmin
     double *ymax
     double *ymin
     double *zerofill
     double *upordown
  */
  double bside=side/sqrt(2);
  int nb=0;
  double *parameters,colorcode=0,xmax=0,ymax=0,xmin=0,ymin=0,zerofill=0;
  double xord1=0,yord1=0,xord2=0,yord2=0;
  struct hist *h=NULL;
  if (void_hist!=NULL){ 
    h=(struct hist *)void_hist;
    if (void_parameters!=NULL){ 
      parameters=(double *)void_parameters;
      colorcode = parameters[0];
      xmax = parameters[1];
      xmin = parameters[2];
      ymax = parameters[3];
      ymin = parameters[4];
      zerofill = parameters[5];}
    else /* if (void_parameters==NULL) */{
      colorcode = 0;
      xmax = h->max; xmin = h->min;
      stats("double",h->data,h->nbins,&ymax,&ymin,NULL,NULL);
      zerofill = 0;}
    if (xmax==xmin){ xmax+=1;} if (ymax==ymin){ ymax+=1;}
    if (zerofill>0){
      yord1 = (0-ymin)/(ymax-ymin); yord2 = (minimum(zerofill,ymax)-ymin)/(ymax-ymin);
      xord1 = ((0 + 0.125/(double)maximum(1,h->nbins)*(h->max-h->min))-xmin)/(xmax-xmin);
      xord2 = ((0 + 0.875/(double)maximum(1,h->nbins)*(h->max-h->min))-xmin)/(xmax-xmin);
      fprintf(fp,"2 1 0 0 0 %d %d 0 %d 0.000 0 0 -1 0 0 %d\n",/*colorcode*/(int)colorcode+32,/*depth*/1,/*fill*/20,/*npoints*/4+1); 
      fprintf(fp,"\t"); 
      fprintf(fp,"%d %d ",(int)floor(xoffset+bside*(2*xord1-1)),(int)floor(yoffset+bside*(2*yord1-1)));
      fprintf(fp,"%d %d ",(int)floor(xoffset+bside*(2*xord1-1)),(int)floor(yoffset+bside*(2*yord2-1)));
      fprintf(fp,"%d %d ",(int)floor(xoffset+bside*(2*xord2-1)),(int)floor(yoffset+bside*(2*yord2-1)));
      fprintf(fp,"%d %d ",(int)floor(xoffset+bside*(2*xord2-1)),(int)floor(yoffset+bside*(2*yord1-1)));
      fprintf(fp,"%d %d ",(int)floor(xoffset+bside*(2*xord1-1)),(int)floor(yoffset+bside*(2*yord1-1)));
      fprintf(fp,"\n");}
    for (nb=0;nb<h->nbins;nb++){
      yord1 = (0-ymin)/(ymax-ymin); yord2 = (h->data[nb]-ymin)/(ymax-ymin);
      xord1 = ((h->min + (double)(nb+0.125)/(double)maximum(1,h->nbins)*(h->max-h->min))-xmin)/(xmax-xmin);
      xord2 = ((h->min + (double)(nb+0.875)/(double)maximum(1,h->nbins)*(h->max-h->min))-xmin)/(xmax-xmin);
      fprintf(fp,"2 1 0 0 0 %d %d 0 %d 0.000 0 0 -1 0 0 %d\n",/*colorcode*/(int)colorcode+32,/*depth*/1,/*fill*/20,/*npoints*/4+1); 
      fprintf(fp,"\t"); 
      fprintf(fp,"%d %d ",(int)floor(xoffset+bside*(2*xord1-1)),(int)floor(yoffset+bside*(2*yord1-1)));
      fprintf(fp,"%d %d ",(int)floor(xoffset+bside*(2*xord1-1)),(int)floor(yoffset+bside*(2*yord2-1)));
      fprintf(fp,"%d %d ",(int)floor(xoffset+bside*(2*xord2-1)),(int)floor(yoffset+bside*(2*yord2-1)));
      fprintf(fp,"%d %d ",(int)floor(xoffset+bside*(2*xord2-1)),(int)floor(yoffset+bside*(2*yord1-1)));
      fprintf(fp,"%d %d ",(int)floor(xoffset+bside*(2*xord1-1)),(int)floor(yoffset+bside*(2*yord1-1)));
      fprintf(fp,"\n");}}
}

void pnode_obsdisthist_tofig(void *void_pn,FILE *fp,void *void_parameters,double side,double xoffset,double yoffset)
{
  /* should fit within side of (xoffset,yoffset) */
  /* parameter list:
     double *histortotal
     double *nh1ormax
     double *nh2ormin
     double logscale_flag
     */
  double rcolor=0,gcolor=0,bcolor=0;
  struct pnode *pn=NULL;
  struct obsdisthist *odh=NULL;
  double *parameters=NULL,*parameters_hist=NULL;
  int histortotal=0,nh1=0,nh2=0;
  double max=0,min=0;
  int logscale_flag=0;
  int colorcode=0;
  if (void_pn!=NULL){ 
    pn=(struct pnode *)void_pn; 
    if (pn->temp!=NULL){
      odh=(struct obsdisthist *)pn->temp;
      if (void_parameters!=NULL){
	parameters=(double *)void_parameters;
	histortotal=(int)parameters[0];
	nh1=(int)parameters[1]; max=parameters[1];
	nh2=(int)parameters[2]; min=parameters[2];
	logscale_flag=parameters[3];}
      else /* if (void_parameters==NULL) */{
	histortotal=+1;	nh1=0; nh2=0; logscale_flag=0;}
      if (histortotal==+1){
	fprintf(fp,"4 0 0 %d 0 %d %d 0.0000 4 %d %d %d %d %0.2f\\001\n",/*depth*/98,/*font*/12,/*point*/5,/*textheight*/(int)(5*13.5),/*textwidth*/(int)(5*9*3),/*xpos*/(int)(xoffset),/*ypos*/(int)(yoffset),odh->value[nh1+nh2*odh->length]);
	parameters_hist = (double *) tcalloc(6,sizeof(double));
	colorscale(0,odh->value[nh1+nh2*odh->length],+1,0,&rcolor,&gcolor,&bcolor);
	parameters_hist[0]=crop((int)floor(512*rcolor),0,511);
	parameters_hist[1]=maximum(1,maximum(odh->hra[nh1]->max,odh->hra[nh2]->max));
	parameters_hist[2]=0;
	stats("double",odh->hra[nh1]->data,odh->hra[nh1]->nbins,&(parameters_hist[3]),NULL,NULL,NULL);
	stats("double",odh->hra[nh2]->data,odh->hra[nh2]->nbins,&(parameters_hist[4]),NULL,NULL,NULL);
	parameters_hist[3] = maximum(parameters_hist[3],parameters_hist[4]);
	parameters_hist[4]=0;
	parameters_hist[5]=odh->nrecords-odh->hra[nh1]->ndatum;
	hist_tofig(odh->hra[nh1],fp,parameters_hist,side,xoffset,yoffset);
	if (nh2!=nh1){
	  parameters_hist[5]=odh->nrecords-odh->hra[nh2]->ndatum;
	  parameters_hist[3] *= -1;
	  hist_tofig(odh->hra[nh2],fp,parameters_hist,side,xoffset,yoffset);}
	tfree(parameters_hist);parameters_hist=NULL;}
      else if (histortotal==-1){
	for (nh1=0;nh1<odh->length;nh1++){
	  fprintf(fp,"4 0 0 %d 0 %d %d 0.0000 4 %d %d %d %d %d\\001\n",/*depth*/98,/*font*/12,/*point*/5,/*textheight*/(int)(5*13.5),/*textwidth*/(int)(5*9*(int)(1+log(maximum(1,odh->hra[nh1]->datasum)))),/* xoffset */(int)floor(xoffset+side*cos(2*PI*((double)nh1+0.5+nh2*0.25)/(double)(odh->length))),/* yoffset */(int)floor(yoffset+side*sin(2*PI*((double)nh1+0.5+nh2*0.25)/(double)(odh->length))),(int)odh->hra[nh1]->datasum);
	  if (logscale_flag){
	    if (odh->hra[nh1]!=NULL && odh->hra[nh1]->nbins>0){ colorscale(0,log(maximum(1,odh->hra[nh1]->datasum)),log(maximum(1,max)),log(maximum(1,min)),&rcolor,&gcolor,&bcolor);}
	    else /* if (odh->hra[nh1]==NULL || odh->hra[nh1]->nbins==0) */{ colorscale(0,1,maximum(1,max),maximum(1,min),&rcolor,&gcolor,&bcolor);}}
	  else if (!logscale_flag){
	    if (odh->hra[nh1]!=NULL && odh->hra[nh1]->nbins>0){ colorscale(0,odh->hra[nh1]->datasum,max,min,&rcolor,&gcolor,&bcolor);}
	    else /* if (odh->hra[nh1]==NULL || odh->hra[nh1]->nbins==0) */{ colorscale(0,0,max,min,&rcolor,&gcolor,&bcolor);}}
	  colorcode=crop((int)floor(512*rcolor),0,511);
	  fprintf(fp,"2 1 0 0 0 %d %d 0 %d 0.000 0 0 -1 0 0 %d\n",/*colorcode*/colorcode+32,/*depth*/1,/*fill*/20,/*npoints*/5+2); fprintf(fp,"\t"); fprintf(fp,"%d %d ",(int)floor(xoffset),(int)floor(yoffset)); for (nh2=-2;nh2<3;nh2++){ fprintf(fp,"%d %d ",(int)floor(xoffset+side*cos(2*PI*((double)nh1+0.5+nh2*0.25)/(double)(odh->length))),(int)floor(yoffset+side*sin(2*PI*((double)nh1+0.5+nh2*0.25)/(double)(odh->length))));} fprintf(fp,"%d %d ",(int)floor(xoffset),(int)floor(yoffset)); fprintf(fp,"\n");}}
      fprintf(fp,"4 0 0 %d 0 %d %d 0.0000 4 %d %d %d %d %d\\001\n",/*depth*/99,/*font*/12,/*point*/24,/*textheight*/(int)(24*13.5),/*textwidth*/(int)(24*9*1),/*xpos*/(int)(xoffset),/*ypos*/(int)(yoffset),pn->region->label);}}
}

void pnode_tofig_helper(struct llist *L,int nregions,struct pnode *parent,struct llitem *l0,int maxlevel,int level,FILE *fp,void (*tofigtemp)(void *,FILE *,void *,double,double,double),void *tofigtemp_parameters)
{
  /* L should be NULL to start with */
  int verbose=0,nr=0;
  char lvlchar[64];
  int firstentry=0,veryfirstentry=0;
  struct litem *l=NULL;
  struct pnode *pn=NULL;
  int maxdia=100;
  struct region *rg=NULL;
  int nl=0;
  double xord=0,yord=0,vx=0,vy=0,vrad=0,vrad_next=0;
  if (verbose){ sprintf(lvlchar," %%"); for (nr=0;nr<level;nr++){sprintf(lvlchar,"%s %%",lvlchar);}}
  if (verbose){ printf(" %s [entering pnode_tofig_helper] nregions %d,maxlevel %d,level %d\n",lvlchar,nregions,maxlevel,level);}
  if (parent==NULL && l0->parent==NULL){ 
    if (verbose){ printf(" %s very first descent into pnode\n",lvlchar);} 
    veryfirstentry=1;}
  if (veryfirstentry){ 
    fprintf(fp,"%s",FIG_PREAMBLE); fprintf(fp,"%s",FIG_PREAMBLE_COLOR_7); 
    if (verbose){ printf(" %s making llist\n",lvlchar);}
    assert(L==NULL);
    L=llistmake();}
  if (parent!=NULL && l0->parent==NULL){ 
    if (verbose){ printf(" %s first descent into this node\n",lvlchar);}
    assert(parent->childllitem==l0); firstentry=1;}
  if (firstentry){ 
    if (verbose){ printf(" %s adding region %d to llist\n",lvlchar,parent->region->label);}
    litemadd(L,parent->region);}
  if (firstentry){ 
    if (verbose){ printf(" %s relevant node, stepping through llist\n",lvlchar);}
    l=L->last;nl=0;xord=0,yord=0,vx=0,vy=0;vrad=maxdia;
    while (l!=NULL){
      rg = (struct region *)l->item;
      vx = cos(2*PI*(rg->label+0.5)/(double)nregions);
      vy = sin(2*PI*(rg->label+0.5)/(double)nregions);
      vrad_next = maximum(2*vrad,vrad/sin(PI/maximum(1,nregions)));
      vrad += vrad_next;
      xord += vrad_next*vx; yord += vrad_next*vy;
      nl+=1;
      if (verbose){ printf(" %s nl %d,label %d,vx %0.2f,vy %0.2f,vrad %0.2f,vrad_next %0.2f,xord %0.2f, yord %0.2f\n",lvlchar,nl,rg->label,vx,vy,vrad,vrad_next,xord,yord);}
      l=l->parent;}
    if (verbose){ printf(" %s calling figmaker\n",lvlchar);}
    tofigtemp(parent,fp,tofigtemp_parameters,maxdia,xord,yord);}
  if (l0->kidl!=NULL){ pnode_tofig_helper(L,nregions,parent,l0->kidl,maxlevel,level,fp,tofigtemp,tofigtemp_parameters);}
  if (l0->item!=NULL){
    if (maxlevel==-1 || level<maxlevel){ 
      if (verbose){ printf(" %s level %d, maxlevel %d, digging deeper\n",lvlchar,level,maxlevel);}
      pn=(struct pnode *)l0->item; 
      pnode_tofig_helper(L,nregions,pn,pn->childllitem,maxlevel,level+1,fp,tofigtemp,tofigtemp_parameters);}}
  if (l0->kidr!=NULL){ pnode_tofig_helper(L,nregions,parent,l0->kidr,maxlevel,level,fp,tofigtemp,tofigtemp_parameters);}
  if (firstentry){ 
    if (verbose){ printf(" %s llistkilllast\n",lvlchar);}
    llistkilllast(L);}
  if (veryfirstentry){ 
    if (verbose){ printf(" %s llisttfree\n",lvlchar);}
    llisttfree(L);L=NULL;}
}

void pnode_obsdisthist_tofig_starter(struct ptree *p,char *fgvn,int ninputs)
{
  /* assumes fgvn does NOT start with "./" */
  FILE *fp=NULL;
  char filename[256];
  int nh1=0,nh2=0;
  double *parameters=NULL;
  struct llist *L=NULL;
  struct litem *l=NULL;
  struct pnode *pn=NULL;
  struct obsdisthist *odh=NULL;
  double max=0;
  if (ninputs>1){
    for (nh1=0;nh1<ninputs;nh1++){ for (nh2=nh1+1;nh2<ninputs;nh2++){
      if (fgvn!=NULL){ sprintf(filename,"./ptree_%s_i%di%d_%srecord.fig",fgvn,nh1,nh2,GLOBAL_STRING_2);}
      else /* if (fgvn==NULL) */{ sprintf(filename,"./ptree_obsdisthist_i%di%d_%srecord.fig",nh1,nh2,GLOBAL_STRING_2);}
      if ((fp=fopen(filename,"w"))==NULL){ printf("warning, cannot create %s in pnode_obsdisthist_tofig_starter\n",filename); fp=stdout;}
      parameters=(double *)tcalloc(4,sizeof(double));
      parameters[0]=+1;
      parameters[1]=nh1;
      parameters[2]=nh2;
      pnode_tofig_helper(NULL,p->nregions,NULL,p->postree,-1,0,fp,&pnode_obsdisthist_tofig,parameters);
      tfree(parameters);parameters=NULL;
      if (fp!=stdout){ fclose(fp);}}}}
  else if (ninputs==1){ 
    if (fgvn!=NULL){ sprintf(filename,"./ptree_%s_i%d_%srecord.fig",fgvn,0,GLOBAL_STRING_2);}
    else /* if (fgvn==NULL) */{ sprintf(filename,"./ptree_obsdisthist_i%d_%srecord.fig",0,GLOBAL_STRING_2);}
    if ((fp=fopen(filename,"w"))==NULL){ printf("warning, cannot create %s in pnode_obsdisthist_tofig_starter\n",filename); fp=stdout;}
    parameters=(double *)tcalloc(4,sizeof(double));
    parameters[0]=+1;
    parameters[1]=0;
    parameters[2]=0;
    pnode_tofig_helper(NULL,p->nregions,NULL,p->postree,-1,0,fp,&pnode_obsdisthist_tofig,parameters);
    tfree(parameters);parameters=NULL;
    if (fp!=stdout){ fclose(fp);}}
  L=llistmake();
  llitem2llist(p->postree,L);
  l=L->first;max=1;
  while (l!=NULL){
    pn = (struct pnode *)l->item; 
    if (pn->temp!=NULL){ 
      odh = (struct obsdisthist *) pn->temp; for (nh1=0;nh1<odh->length;nh1++){ max = maximum(max,odh->hra[nh1]->datasum);}}
    l=l->child;}
  llisttfree(L);L=NULL;
  if (fgvn!=NULL){ sprintf(filename,"./ptree_%s_sum_%srecord.fig",fgvn,GLOBAL_STRING_2);}
  else /* if (fgvn==NULL) */{ sprintf(filename,"./ptree_obsdisthist_sum_%srecord.fig",GLOBAL_STRING_2);}
  if ((fp=fopen(filename,"w"))==NULL){ printf("warning, cannot create %s in pnode_obsdisthist_tofig_starter\n",filename); fp=stdout;}
  parameters=(double *)tcalloc(4,sizeof(double));
  parameters[0]=-1;
  parameters[1]=max;
  parameters[2]=0;
  parameters[3]=0;
  pnode_tofig_helper(NULL,p->nregions,NULL,p->postree,-1,0,fp,&pnode_obsdisthist_tofig,parameters);
  tfree(parameters);parameters=NULL;
  if (fp!=stdout){ fclose(fp);}
  if (fgvn!=NULL){ sprintf(filename,"./ptree_%s_logsum_%srecord.fig",fgvn,GLOBAL_STRING_2);}
  else /* if (fgvn==NULL) */{ sprintf(filename,"./ptree_obsdisthist_logsum_%srecord.fig",GLOBAL_STRING_2);}
  if ((fp=fopen(filename,"w"))==NULL){ printf("warning, cannot create %s in pnode_obsdisthist_tofig_starter\n",filename); fp=stdout;}
  parameters=(double *)tcalloc(4,sizeof(double));
  parameters[0]=-1;
  parameters[1]=max;
  parameters[2]=0;
  parameters[3]=1;
  pnode_tofig_helper(NULL,p->nregions,NULL,p->postree,-1,0,fp,&pnode_obsdisthist_tofig,parameters);
  tfree(parameters);parameters=NULL;
  if (fp!=stdout){ fclose(fp);}
}

double pnode_obsdisthist_value(void *void_pn,void *void_parameters)
{
  /* parameter list:
     double nh1
     double nh2 */
  double output=0;
  double *parameters=NULL;
  int nh1=0,nh2=0;
  struct pnode *pn=NULL;
  struct obsdisthist *odh=NULL;
  if (void_parameters!=NULL){ parameters = (double *)void_parameters; nh1 = (int)parameters[0]; nh2 = (int)parameters[1];}
  else /* if (void_parameters==NULL) */{ nh1 = 0; nh2 = 1;}
  if (void_pn!=NULL){ 
    pn=(struct pnode *)void_pn; 
    if (pn->temp!=NULL){ 
      odh=(struct obsdisthist *)pn->temp; 
      output = odh->value[nh1+nh2*odh->length];}}
  return output;
}

void pnode_statstemp_starter(struct pnode *parent,struct llitem *l0,int maxlevel,int minlevel,int level,int *nelements,double *temp_max,double *temp_min,double *temp_mean,double *temp_stdev,double (*temp_evaluate)(void *,void *),void *temp_evaluate_parameters)
{
  /* start out with nelements==NULL, stdev only accumulates if mean!=NULL */
  int firstentry=0,veryfirstentry=0,madenelements_flag=0;
  struct pnode *pn=NULL;
  double temp=0;
  if (parent==NULL && l0->parent==NULL){ /* very first descent into ptree */
    if (nelements==NULL){ madenelements_flag=1; nelements=(int *)tmalloc(sizeof(int));} else{ madenelements_flag=0;} *nelements=0;
    if (temp_min!=NULL){ *temp_min=0;} if (temp_max!=NULL){ *temp_max=0;} if (temp_mean!=NULL){ *temp_mean=0; if (temp_stdev!=NULL){ *temp_stdev=0;}}
    veryfirstentry=1;}
  if (parent!=NULL && l0->parent==NULL){ /* first descent into parent->childllitem */ 
    assert(parent->childllitem==l0); assert(nelements!=NULL); firstentry=1;}
  if (l0->kidl!=NULL){ pnode_statstemp_starter(parent,l0->kidl,maxlevel,minlevel,level,nelements,temp_max,temp_min,temp_mean,temp_stdev,temp_evaluate,temp_evaluate_parameters);}
  if (l0->kidr!=NULL){ pnode_statstemp_starter(parent,l0->kidr,maxlevel,minlevel,level,nelements,temp_max,temp_min,temp_mean,temp_stdev,temp_evaluate,temp_evaluate_parameters);}
  if (l0->item!=NULL){
    if (maxlevel==-1 || level<maxlevel){ 
      pn=(struct pnode *)l0->item; 
      pnode_statstemp_starter(pn,pn->childllitem,maxlevel,minlevel,level+1,nelements,temp_max,temp_min,temp_mean,temp_stdev,temp_evaluate,temp_evaluate_parameters);}}
  if (firstentry){ if (parent!=NULL){ if (nelements!=NULL){ if (minlevel==-1 || level>=minlevel){
    temp = temp_evaluate(parent,temp_evaluate_parameters);
    if (*nelements==0){ 
      if (temp_max!=NULL){ *temp_max=temp;}
      if (temp_min!=NULL){ *temp_min=temp;}}
    else if (*nelements>0){
      if (temp_max!=NULL){ *temp_max=maximum(*temp_max,temp);}
      if (temp_min!=NULL){ *temp_min=minimum(*temp_min,temp);}}
    if (temp_mean!=NULL){ *temp_mean += temp; if (temp_stdev!=NULL){ *temp_stdev += pow(temp,2);}}
    *nelements += 1;}}}}
  if (veryfirstentry){ if (nelements!=NULL){ 
    if (temp_mean!=NULL){ *temp_mean /= maximum(1.0,(double)(*nelements)); if (temp_stdev!=NULL){ *temp_stdev /= maximum(1.0,(double)(*nelements)); *temp_stdev -= pow(*temp_mean,2); *temp_stdev = sqrt(*temp_stdev);}}
    if (madenelements_flag==1){ tfree(nelements);nelements=NULL;}}}
}

void pnode_normanddots(void *v1,void *v2,void *void_parameters)
{
  /* determines dot products of weights,and dot products of relevances */
  /* parameter list
     (void *) (double *) p1->weight dot p1->weight
     (void *) (double *) p2->weight dot p2->weight
     (void *) (double *) p1->weight dot p2->weight
     (void *) (double *) p1->relevance dot p1->relevance
     (void *) (double *) p2->relevance dot p2->relevance
     (void *) (double *) p1->relevance dot p2->relevance */
  int verbose=0;
  struct pnode *pn1=NULL,*pn2=NULL;
  void **vrara=NULL;
  double *w11=NULL,*w22=NULL,*w12=NULL,*r11=NULL,*r22=NULL,*r12=NULL;
  if (verbose){ printf(" %% [entering pnode_normanddots]\n");}
  if (void_parameters!=NULL){
    vrara = (void **) void_parameters;
    w11 = vrara[0]; w22 = vrara[1]; w12 = vrara[2]; r11 = vrara[3]; r22 = vrara[4]; r12 = vrara[5];
    if (v1!=NULL){
      pn1 = (struct pnode *)v1;
      if (verbose){ printf(" %% p1->weight %f relevance %f\n",pn1->weight,pn1->relevance);}
      if (w11!=NULL){ *w11 += pn1->weight*pn1->weight;}
      if (r11!=NULL){ *r11 += pn1->relevance*pn1->relevance;}}
    if (v2!=NULL){
      pn2 = (struct pnode *)v2;
      if (verbose){ printf(" %% p2->weight %f relevance %f\n",pn2->weight,pn2->relevance);}
      if (w22!=NULL){ *w22 += pn2->weight*pn2->weight;}
      if (r22!=NULL){ *r22 += pn2->relevance*pn2->relevance;}}
    if (pn1!=NULL && pn2!=NULL){
      if (w12!=NULL){ *w12 += pn1->weight*pn2->weight;}
      if (r12!=NULL){ *r12 += pn1->relevance*pn2->relevance;}}}
}

void pnode_obsdisthist2obsdisthist_appraise(void *v1,void *v2,void *void_parameters)
{
  /* determines maximum separability of odh1->hra[0] and odh2->hra[0], */
  /* parameter list
     (double *) parameters, storing [nrecords,bestsofar] */
  int verbose=0;
  int nrecords=0,h1_freeflag=0,h2_freeflag=0;
  struct pnode *pn1=NULL,*pn2=NULL;
  struct obsdisthist *odh1=NULL,*odh2=NULL;
  struct hist *h1=NULL,*h2=NULL;
  double *parameters=NULL,bestsofar=0,best1=0,best2=0;
  if (verbose){ printf(" %% [entering pnode_obsdisthist2obsdisthist_appraise]\n");}
  parameters = (double *) void_parameters;
  nrecords = (int)parameters[0]; bestsofar = (double)parameters[1];
  pn1=NULL;odh1=NULL;pn2=NULL;odh2=NULL;
  if (v1!=NULL){
    pn1=(struct pnode *)v1; if (pn1->temp!=NULL){ odh1=(struct obsdisthist *)pn1->temp;} else /* if (pn1->temp==NULL) */{ odh1=NULL;}}
  if (v2!=NULL){
    pn2=(struct pnode *)v2; if (pn2->temp!=NULL){ odh2=(struct obsdisthist *)pn2->temp;} else /* if (pn2->temp==NULL) */{ odh2=NULL;}}
  if (odh1==NULL){ h1=histmake(0,0,0);h1_freeflag=1;} else /* if (odh1!=NULL) */{ h1=odh1->hra[0]; h1_freeflag=0;}
  if (odh2==NULL){ h2=histmake(0,0,0);h2_freeflag=1;} else /* if (odh2!=NULL) */{ h2=odh2->hra[0]; h2_freeflag=0;}
  if (verbose){ 
    printf(" %% nrecords %d, bestsofar %0.1f, v1%d,v2%d\n",nrecords,bestsofar,(int)v1,(int)v2);    
    histprintf(h1," %% h1"); histprintf(h2," %% h2");}
  best1 = obsdisthistappraise_helper(h1,nrecords-h1->ndatum,h2,nrecords-h2->ndatum); 
  best2 = obsdisthistappraise_helper(h2,nrecords-h2->ndatum,h1,nrecords-h1->ndatum); 
  parameters[1] = maximum(bestsofar,best1+best2);
  if (verbose){ printf(" %% now parameters[1]=%0.1f >=  (%0.1f+%0.1f)\n",parameters[1],best1,best2);}
  if (h1_freeflag){ histtfree(h1); h1=NULL;}
  if (h2_freeflag){ histtfree(h2); h2=NULL;}
}

void pnode2pnode_void_operate(struct ptree *pA,struct pnode *parentA,struct llitem *lA,struct ptree *pB,struct pnode *parentB,struct llitem *lB,int maxlevel,int level,void (*void_operate)(void *,void *,void *),void *void_parameters,int swap_flag)
{
  int verbose=0;
  char lvlchar[64];
  int nl=0,firstentry=0;
  struct pnode *pnA=NULL,*pnB=NULL;
  struct llitem *lB2=NULL;
  if (verbose){ sprintf(lvlchar," %%"); for (nl=0;nl<level;nl++){ sprintf(lvlchar,"%s %%",lvlchar);}}
  if (verbose && parentA==NULL && lA->parent==NULL){ printf(" %s [entering pnode2pnode_void_operate]\n",lvlchar);}
  if (verbose>2){ 
    printf(" %% parentA %d, lA-> %d\n",(int)parentA,(lA!=NULL?(int)lA->item:-1)); pnodeprintf(parentA,lA,-1,0);
    printf(" %% parentB %d, lB-> %d\n",(int)parentB,(lB!=NULL?(int)lB->item:-1)); pnodeprintf(parentB,lB,-1,0);}
  if (parentA!=NULL && lA->parent==NULL){ /* first descent into parent->childllitem */ 
    if (verbose>1){ printf(" %s first descent into parentA->childllitem, setting firstentry to 1\n",lvlchar);}
    assert(parentA->childllitem==lA); firstentry=1;}
  if (lA->kidl!=NULL){ 
    if (verbose>1){ printf(" %s stepping left to lA->kidl\n",lvlchar);}
    pnode2pnode_void_operate(pA,parentA,lA->kidl,pB,parentB,lB,maxlevel,level,void_operate,void_parameters,swap_flag);}
  if (lA->kidr!=NULL){ 
    if (verbose>1){ printf(" %s stepping right to lA->kidr\n",lvlchar);}
    pnode2pnode_void_operate(pA,parentA,lA->kidr,pB,parentB,lB,maxlevel,level,void_operate,void_parameters,swap_flag);}
  if (maxlevel==-1 || level<maxlevel){ 
    if (lA->item!=NULL){
      if (verbose>1){ printf(" %s level %d less than maxlevel %d, and lA->item exists, delving deeper\n",lvlchar,level,maxlevel);}
      pnA= (struct pnode *) lA->item; 
      if (lB!=NULL && (lB2 = llitemaddorfind(0,lB,pB->regionra[pnA->region->label],&region2pnode_compare_label))!=NULL){
	pnB = (struct pnode *) lB2->item;}
      else{ pnB=NULL;}
      lB2 = (pnB==NULL ? NULL : pnB->childllitem);
      if (verbose>1){ 
	if (pnB!=NULL){ printf(" %s pnB label %d matches pnA label %d, ...\n",lvlchar,pnB->region->label,pnA->region->label);}
	else{ printf(" %s pnB label doesn't match pnA label %d, ...\n",lvlchar,pnA->region->label);}
	printf(" %s descending\n",lvlchar);}
      pnode2pnode_void_operate(pA,pnA,pnA->childllitem,pB,pnB,lB2,maxlevel,level+1,void_operate,void_parameters,swap_flag);}
    else /* if (lA->item==NULL) */{
      if (verbose>1){ printf(" %s level %d less than maxlevel %d, and lA->item doesn't exist...",lvlchar,level,maxlevel);}
      if (lB!=NULL && lB->item!=NULL){
	if (verbose>1){ printf(" but lB->item does exist, swapping.\n");}
	pnB = (struct pnode *) lB->item;
	if (pnB!=NULL && pnB->childllitem!=NULL){
	  pnode2pnode_void_operate(pB,pnB,pnB->childllitem,pA,NULL,NULL,maxlevel,level+1,void_operate,void_parameters,!swap_flag);}}
      else /* if (lB==NULL || lB->itime==NULL) */{ if (verbose>1){ printf(" lB->item doesn't exist either, stopping.\n");}}}}
  if (firstentry){
    if (verbose>1){ printf(" %s first entry for region %d,\n",lvlchar,parentA->region->label);}
    if (swap_flag){ void_operate(parentB,parentA,void_parameters);}
    else /* if (!swap_flag) */{ void_operate(parentA,parentB,void_parameters);}}
} 

void pnode_set_weight(void *void_pnode,void *void_parameters)
{
  /* parameter list
     double * temp */
  int verbose=0;
  struct pnode *pn=NULL;
  double *temp=NULL;
  if (void_parameters!=NULL){
    pn = (struct pnode *) void_pnode;
    temp = (double *) void_parameters;
    pn->weight = *temp;
    if (verbose){ printf(" %% weight %f\n",pn->weight);}}
  else /* if (void_parameters==NULL) */{ if (verbose){ printf(" %% warning, empty parameter list in pnode_set_weight\n");}}
}

void pnode_evaluate_relevance(void *void_pnode,void *void_parameters)
{
  /* parameter list
     double * temp */
  int verbose=0;
  struct pnode *pn=NULL;
  double *temp=NULL;
  if (void_parameters!=NULL){
    pn = (struct pnode *) void_pnode;
    temp = (double *) void_parameters;
    *temp = pn->relevance;
    if (verbose){ printf(" %% relevance %f\n",pn->relevance);}}
  else /* if (void_parameters==NULL) */{ if (verbose){ printf(" %% warning, empty parameter list in pnode_evaluate_relevance\n");}}
}

void pnode_evaluate_obsdisthist_mean(void *void_pnode,void *void_parameters)
{
  /* parameter list
     double nh1
     double nrecords
     double output */
  int verbose=0;
  struct pnode *pn=NULL; struct obsdisthist *odh=NULL;
  double *dra=NULL; int nh1=0; double nrecords=0; struct hist *h=NULL;
  if (void_parameters!=NULL){
    pn = (struct pnode *) void_pnode; odh = (struct obsdisthist *) pn->temp;
    dra = (double *) void_parameters; nh1 = (int) dra[0]; nrecords = (double) dra[1]; h = odh->hra[nh1];
    if (verbose){ histprintf(h,"h: ");}
    dra[2] = h->datasum/(double)maximum(1,nrecords/* and *NOT* h->ndatum, due to the unknown zerofill */);
    if (verbose){ printf(" %% mean %f\n",dra[2]);}}
  else /* if (void_parameters==NULL) */{ if (verbose){ printf(" %% warning, empty parameter list in pnode_evaluate_obsdisthist_mean\n");}} 
}

void pnode_evaluate_weight(void *void_pnode,void *void_parameters)
{
  /* parameter list
     double * output */
  int verbose=0;
  struct pnode *pn=NULL;
  double *output=NULL;
  if (void_parameters!=NULL){
    pn = (struct pnode *) void_pnode; output = (double *) void_parameters;
    *output = pn->weight; if (verbose){ printf(" %% weight %f\n",*output);}}
  else /* if (void_parameters==NULL) */{ if (verbose){ printf(" %% warning, empty parameter list in pnode_evaluate_weight\n");}}
}

void pnode_dig_and_evaluate(int llist_type,struct litem *l,struct ptree *p,struct pnode *parent,struct llitem *l0,void (*pnode_evaluate)(void *,void *),void *void_parameters)
{
  /* llist_type:
     0: pnodes
     1: regions */
  int verbose=0;
  struct region *rg=NULL;
  struct llitem *l1=NULL;
  struct pnode *pn=NULL;
  if (l!=NULL){
    if (llist_type==0){ pn = (struct pnode *) l->item; rg = pn->region;}
    else if (llist_type==1){ rg = (struct region *) l->item;}
    else{ printf("warning! improper llist_type in pnode_dig_and_evaluate\n");}
    if (verbose){ printf(" %% digging for region->%d in pnode_dig_and_evaluate... ",rg->label);}
    if (l0!=NULL && (l1 = llitemaddorfind(0,l0,rg,&region2pnode_compare_label))!=NULL){
      pn = (struct pnode *) l1->item;
      if (verbose){ printf(" found label %d which matches %d, digging deeper\n",pn->region->label,rg->label);}
      pnode_dig_and_evaluate(llist_type,l->child,p,pn,pn->childllitem,pnode_evaluate,void_parameters);}
    else /* if (l0==NULL || l1==NULL) */{
      if (verbose){ printf(" nothing found, stopping\n");}}}
  else /* if (l==NULL) */{
    if (parent!=NULL){
      if (verbose){ printf(" %% done digging, evaluating at label %d\n",parent->region->label);}
      pnode_evaluate(parent,void_parameters);}
    else /* if (parent==NULL) */{ if (verbose){ printf(" %% warning! NULL litem, but no parent in pnode_dig_and_evaluate\n");}}}
}

void pnode_obsdisthist2mtrx(void *vA,void *vB,void *void_parameters)
{
  /* it is assumed that pnA=(struct pnode *)vA is a member of pA->postree */
  /* parameter list 
     void * (struct ptree *)pA
     void * (double *)ra (of dimension pow(pA->nregions,*level))
     void * (int *)level
     void * (double *)total_time 
     void * (double *)nrecords */
  int verbose=0;
  void **vrara=NULL;
  struct ptree *pA=NULL;
  struct llist *LA=NULL;
  struct litem *l=NULL;
  struct pnode *pnA=NULL,*pnB=NULL;
  double *ra=NULL;
  int *level=NULL,tab=0,nl=0;
  double *total_time=NULL,*nrecords=NULL,weight=0;
  struct obsdisthist *odh=NULL;
  struct region *rg=NULL;
  if (verbose){ printf(" %% [entering pnode_obsdisthist2mtrx]\n");}
  if (void_parameters!=NULL){
    vrara = (void **) void_parameters;
    pA = (struct ptree *) vrara[0]; 
    ra = (double *) vrara[1]; 
    level = (int *)vrara[2]; 
    total_time = (double *)vrara[3]; 
    nrecords = (double *)vrara[4];
    if (vA!=NULL){ 
      pnA = (struct pnode *)vA; 
      LA = llistmake();
      pnB = pnA; while (pnB!=NULL){ litemadd(LA,pnB->region); pnB=pnB->parent;}
      if (LA->length==*level){
	if (verbose){ printf(" %% found valid node, with region label llist "); l=LA->last; while (l!=NULL){ rg=(struct region *)l->item; printf("-->%d",rg->label); l=l->parent;} printf("... ");}
	tab = 0;
	l=LA->last; nl=0;
	while (l!=NULL){
	  rg = (struct region *) l->item;
	  tab += rg->label*(int)pow(pA->nregions,nl);
	  l=l->parent;
	  nl+=1;}
	if (verbose){ printf("tab %d... ",tab);}
	ra[tab] = 0;
	if (pnA->temp!=NULL){
	  odh = (struct obsdisthist *) pnA->temp; 
	  weight = (double)odh->hra[0]->datasum/(double)maximum(1,*nrecords)/(*total_time);
	  if (verbose){ printf(" weight_rate %0.2f\n",weight);}
	  ra[tab] = weight;}
	else /* if (pnA->temp==NULL) */{if (verbose){ printf(" odh not found\n");}}}
      if (verbose){ printf(" %% invalid node, LA->length %d\n",LA->length);}
      llisttfree(LA);LA=NULL;}}
  if (verbose){ printf(" %% [finished pnode_obsdisthist2mtrx]\n");}
}

void pnode_obsdisthist2obsdisthist_meansubtmean(void *vA,void *vB,void *void_parameters)
{
  /* it is assumed that pnA=(struct pnode *)vA and pnB=(struct pnode *)vB are members of pA->postree and pB->postree respectively */
  /* stores the difference between the mean of pnA->odh->hra[0] and the mean of pnB->odh->hra[0] as pnD->relevance */
  /* parameter list 
     void * (struct ptree *)pA
     void * (struct ptree *)pB
     void * (struct ptree *)pD 
     void * (double *)total_time 
     void * (double *)nrecords */
  int verbose=0;
  void **vrara=NULL;
  struct ptree *pA=NULL,*pB=NULL,*pD=NULL;
  struct llist *LD=NULL;
  struct litem *l=NULL;
  struct pnode *pnA=NULL,*pnB=NULL,*pnD=NULL;
  struct llitem *lD=NULL,*lD2=NULL;
  double *total_time=NULL,*nrecords=NULL;
  double weight_A=0,weight_B=0;
  struct obsdisthist *odh=NULL;
  struct region *rg=NULL;
  if (verbose){ printf(" %% [entering pnode_obsdisthist2obsdisthist_meansubtmean]\n");}
  if (void_parameters!=NULL){
    vrara = (void **) void_parameters;
    pA = (struct ptree *) vrara[0]; pB = (struct ptree *) vrara[1]; pD = (struct ptree *) vrara[2]; 
    total_time = (double *)vrara[3]; nrecords = (double *)vrara[4];
    weight_A = 0;
    if (vA!=NULL){ 
      pnA = (struct pnode *)vA; 
      if (pnA->temp!=NULL){
	odh = (struct obsdisthist *) pnA->temp; 
	if (verbose){ obsdisthistprintf(odh,"pnA->temp: ");}
	weight_A = (double)odh->hra[0]->datasum/(double)maximum(1,*nrecords);}
      else /* if (pnA->temp==NULL) */{ if (verbose){ printf("pnA->temp==NULL\n");}}}
    weight_B = 0;
    if (vB!=NULL){ 
      pnB = (struct pnode *)vB; 
      if (pnB->temp!=NULL){
	odh = (struct obsdisthist *) pnB->temp; 
	if (verbose){ obsdisthistprintf(odh,"pnB->temp: ");}
	weight_B = (double)odh->hra[0]->datasum/(double)maximum(1,*nrecords);}
      else /* if (pnB->temp==NULL) */{ if (verbose){ printf("pnB->temp==NULL\n");}}}
    LD = llistmake(); pnD = (pnA==NULL ? pnB : pnA); while (pnD!=NULL){ if (LD->length==0){ litemadd(LD,pnD->region);} else /* if (LD->length>0) */{ l = liteminsertbefore(LD,LD->first); l->item = pnD->region;} pnD=pnD->parent;} if (verbose){ printf(" %% region list:"); l=LD->first; while (l!=NULL){ rg=(struct region *)l->item; printf(" %d",rg->label); l=l->child;} printf("\n");}
    l=LD->first; pnD = NULL; lD = pD->postree;
    while (l!=NULL){ 
      if ((lD2=llitemaddorfind(0,lD,(struct region *)l->item,&region2pnode_compare_label))==NULL){
	if (verbose){ printf(" %% region not found, making... ");}
	lD2=llitemaddorfind(1,lD,pnodemake(pnD,(struct region *)l->item,0,0),&pnode2pnode_compare_label);}
      else /* if found */{ if (verbose){ printf(" %% region found, not making... ");}}
      pnD = (struct pnode *) lD2->item; lD = pnD->childllitem;
      if (verbose){ printf("now at region %d\n",pnD->region->label);}
      l=l->child;}
    if (verbose){ printf(" %% finally at region %d, updating relevance with %f-%f=%f\n",pnD->region->label,weight_A,weight_B,weight_A - weight_B);}
    pnD->relevance = weight_A - weight_B;
    llisttfree(LD);LD=NULL;}
  if (verbose){ printf(" %% [finished pnode_obsdisthist2obsdisthist_meansubtmean]\n");}
}

void pnode_obsdisthist2reweight_obsdisthist_mean_shift(void *vA,void *vB,void *void_parameters)
{
  /* it is assumed that pnA=(struct pnode *)vA is a member of pA->postree */
  /* uses the reweight_matrix to predict a shift in the mean of pnA->odh->hra[0], and stores the result in pnW->relevance */
  /* first pass only focuses on the 1-chains (2 events per chain) */
  /* second pass only focuses on the 0-chains (1 event per chain) */
  /* parameter list 
     void * (struct ptree *)pA
     void * (double *) reweight_matrix
     void * (struct ptree *)pW 
     void * (double *)total_time
     void * (double *)nrecords
     void * (int *) secondorfirstpass */
  int verbose=0;
  void **vrara=NULL;
  struct ptree *pA=NULL,*pW=NULL;
  struct llist *LA=NULL;
  struct litem *l=NULL;
  struct pnode *pnA=NULL,*pnW=NULL;
  struct llitem *lW=NULL,*lW2=NULL;
  double *reweight_matrix=NULL,*dra=NULL,*total_time=NULL,*nrecords=NULL;
  int *secondorfirstpass=NULL;
  double weight_12A=0,weight_1A=0,weight_2A=0,weight_A=0;
  struct obsdisthist *odh=NULL;
  int nr1=0;
  if (verbose){ printf(" %% [entering pnode_obsdisthist2reweight_obsdisthist_mean_shift]\n");}
  if (void_parameters!=NULL){
    vrara = (void **) void_parameters;
    pA = (struct ptree *) vrara[0]; reweight_matrix = (double *) vrara[1]; pW = (struct ptree *) vrara[2]; 
    total_time = (double *)vrara[3]; nrecords = (double *)vrara[4]; secondorfirstpass = (int *)vrara[5];
    if (vA!=NULL){ 
      pnA = (struct pnode *)vA; 
      if (*secondorfirstpass==1 && pnA->parent==NULL){ 
	if (verbose){ printf(" %% 0-chain region %d, assume that all descendant 1-chains are processed\n",pnA->region->label);}
	weight_A=0;
	for (nr1=0;nr1<pW->nregions;nr1++){ 
	  pnW = NULL; lW = pW->postree;
	  if ((lW2=llitemaddorfind(0,lW,pW->regionra[nr1],&region2pnode_compare_label))!=NULL){
	    pnW = (struct pnode *) lW2->item;
	    if ((lW2=llitemaddorfind(0,pnW->childllitem,pW->regionra[pnA->region->label],&region2pnode_compare_label))!=NULL){
	      pnW = (struct pnode *) lW2->item;
	      if (verbose){ printf(" %% checking 2-chain %d->%d %f\n",pnW->parent->region->label,pnW->region->label,pnW->relevance);}
	      weight_A += pnW->relevance;}}}
	if (verbose){ printf(" %% accumulated %f\n",weight_A);}
	if ((lW=llitemaddorfind(0,pW->postree,pW->regionra[pnA->region->label],&region2pnode_compare_label))!=NULL){
	  pnW = (struct pnode *) lW->item;
	  pnW->relevance = weight_A;}}
      if (*secondorfirstpass==0 && pnA->parent!=NULL && pnA->parent->parent==NULL){
	if (pnA->temp!=NULL){
	  odh = (struct obsdisthist *) pnA->temp; 
	  if (verbose){ printf(" %% 1-chain, %d->%d\n",pnA->parent->region->label,pnA->region->label);}
	  if (verbose){ obsdisthistprintf(odh,"pnA->temp: ");}
	  weight_12A = (double)odh->hra[0]->datasum/(double)maximum(1,*nrecords);
	  if (pnA->parent!=NULL && pnA->parent->temp!=NULL){
	    odh = (struct obsdisthist *) pnA->parent->temp;
	    if (verbose){ obsdisthistprintf(odh,"pnA->parent->temp: ");}
	    weight_1A = (double)odh->hra[0]->datasum/(double)maximum(1,*nrecords);}
	  LA = llistmake(); litemadd(LA,pnA->region);
	  dra = (double *) tcalloc(3,sizeof(double)); dra[0]=0; dra[1]=*nrecords;
	  if (verbose){ printf(" %% pnA->temp exists, digging through pA to find label %d\n",pnA->region->label);}
	  pnode_dig_and_evaluate(1,LA->first,pA,NULL,pA->postree,&pnode_evaluate_obsdisthist_mean,dra);
	  weight_2A = dra[2];
	  tfree(dra);dra=NULL;
	  llisttfree(LA);LA=NULL;
	  weight_A = (weight_12A - weight_1A*weight_2A*pA->legtime/(double)maximum(1,*total_time));
	  weight_A *= reweight_matrix[pnA->parent->region->label + pnA->region->label*pA->nregions];
	  LA = llistmake(); litemadd(LA,pnA->parent->region); litemadd(LA,pnA->region);
	  l=LA->first; pnW = NULL; lW = pW->postree;
	  while (l!=NULL){ 
	    if ((lW2=llitemaddorfind(0,lW,(struct region *)l->item,&region2pnode_compare_label))==NULL){
	      if (verbose){ printf(" %% region not found, making... ");}
	      lW2=llitemaddorfind(1,lW,pnodemake(pnW,(struct region *)l->item,0,0),&pnode2pnode_compare_label);}
	    else /* if found */{ if (verbose){ printf(" %% region found, not making... ");}}
	    pnW = (struct pnode *) lW2->item; lW = pnW->childllitem;
	    if (verbose){ printf("now at region %d\n",pnW->region->label);}
	    l=l->child;}
	  if (verbose){ printf(" %% finally at region %d, updating relevance %0.1f\n",pnW->region->label,weight_A);}
	  pnW->relevance = weight_A;
	  llisttfree(LA);LA=NULL;}}}}
  if (verbose){ printf(" %% [finished pnode_obsdisthist2reweight_obsdisthist_mean]\n");}
}

void pnode_obsdisthist2obsdisthist_reweight(void *vA,void *vB,void *void_parameters)
{
  /* it is assumed that pnA=(struct pnode *)vA and pnB=(struct pnode *)vB are members of pA->postree and pB->postree respectively */
  /* estimates reweight_matrix (pW) which could potentially shift pA into pB (including shift of pB itself only if considerB_flag!=0) */
  /* parameter list 
     void * (struct ptree *)pA
     void * (struct ptree *)pB
     void * (struct ptree *)pW 
     void * (double *)total_time 
     void * (double *)nrecords 
     void * (int *)considerB_flag */
  int verbose=0;
  void **vrara=NULL;
  struct ptree *pA=NULL,*pB=NULL,*pW=NULL;
  struct llist *LA=NULL,*LB=NULL;
  struct litem *l=NULL;
  struct pnode *pnA=NULL,*pnB=NULL,*pnW=NULL;
  struct llitem *lW=NULL,*lW2=NULL;
  double *dra=NULL,*total_time=NULL,*nrecords=NULL;
  int *considerB_flag=NULL;
  double weight_12A=0,weight_1A=0,weight_2A=0,weight_12B=0,weight_1B=0,weight_2B=0,weight_A=0,weight_B=0;
  struct obsdisthist *odh=NULL;
  struct region *rg=NULL;
  if (verbose){ printf(" %% [entering pnode_obsdisthist2obsdisthist_reweight]\n");}
  if (void_parameters!=NULL){
    vrara = (void **) void_parameters;
    pA = (struct ptree *) vrara[0]; pB = (struct ptree *) vrara[1]; pW = (struct ptree *) vrara[2]; 
    total_time = (double *)vrara[3]; nrecords = (double *)vrara[4]; considerB_flag = (int *) vrara[5];
    if (vA!=NULL){ 
      pnA = (struct pnode *)vA; 
      if (pnA->temp!=NULL){
	odh = (struct obsdisthist *) pnA->temp; 
	if (verbose){ obsdisthistprintf(odh,"pnA->temp: ");}
	weight_12A = (double)odh->hra[0]->datasum/(double)maximum(1,*nrecords);
	if (pnA->parent!=NULL && pnA->parent->temp!=NULL){
	  odh = (struct obsdisthist *) pnA->parent->temp;
	  if (verbose){ obsdisthistprintf(odh,"pnA->parent->temp: ");}
	  weight_1A = (double)odh->hra[0]->datasum/(double)maximum(1,*nrecords);}
	LA = llistmake(); litemadd(LA,pnA);
	dra = (double *) tcalloc(3,sizeof(double)); dra[0]=0; dra[1]=*nrecords;
	if (verbose){ printf(" %% pnA->temp exists, digging through pA to find label %d\n",pnA->region->label);}
	pnode_dig_and_evaluate(0,LA->first,pA,NULL,pA->postree,&pnode_evaluate_obsdisthist_mean,dra);
	weight_2A = dra[2];
	tfree(dra);dra=NULL;
	llisttfree(LA);LA=NULL;}
      else /* if (pnA->temp==NULL) */{ if (verbose){ printf("pnA->temp==NULL\n");}}}
    if (vB!=NULL){
      pnB = (struct pnode *)vB;
      if (pnB->temp!=NULL){
	odh = (struct obsdisthist *) pnB->temp; 
	if (verbose){ obsdisthistprintf(odh,"pnB->temp: ");}
	weight_12B = (double)odh->hra[0]->datasum/(double)maximum(1,*nrecords);
	if (pnB->parent!=NULL && pnB->parent->temp!=NULL){
	  odh = (struct obsdisthist *) pnB->parent->temp;
	  if (verbose){ obsdisthistprintf(odh,"pnB->parent->temp: ");}
	  weight_1B = (double)odh->hra[0]->datasum/(double)maximum(1,*nrecords);}
	LB = llistmake(); litemadd(LB,pnB);
	dra = (double *) tcalloc(3,sizeof(double)); dra[0]=0; dra[1]=*nrecords;
	if (verbose){ printf(" %% pnB->temp exists, digging through pB to find label %d\n",pnB->region->label);}
	pnode_dig_and_evaluate(0,LB->first,pB,NULL,pB->postree,&pnode_evaluate_obsdisthist_mean,dra);
	weight_2B = dra[2];
	tfree(dra);dra=NULL;
	llisttfree(LB);LB=NULL;}
      else /* if (pnB->temp==NULL) */{ if (verbose){ printf("pnB->temp==NULL\n");}}}
    if (*considerB_flag==0){
      weight_A = (weight_12B - weight_12A)/(weight_12A - weight_1A*weight_2A*pA->legtime/maximum(1,*total_time));
      weight_B = (weight_12A - weight_12B)/(weight_12B - weight_1B*weight_2B*pB->legtime/maximum(1,*total_time));
      if (verbose){ printf(" %% considerB_flag %d, weight_A = (%f-%f)/(%f-%f*%f/%f)\n",*considerB_flag,weight_12B,weight_12A,weight_12A,weight_1A,weight_2A,*total_time/pA->legtime); printf(" %% weight_B = (%f-%f)/(%f-%f*%f/%f)\n",weight_12A,weight_12B,weight_12B,weight_1B,weight_2B,*total_time/pB->legtime);}}
    else /* if (considerB_flag!=0) */{
      weight_A = (weight_12B - weight_12A)/(weight_12A - weight_1A*weight_2A*pA->legtime/maximum(1,*total_time) - weight_12B + weight_1B*weight_2B*pB->legtime/maximum(1,*total_time));
      if (verbose){ printf(" %% considerB_flag %d, weight_A = (%f-%f)/(%f-%f*%f/%f-%f+%f*%f/%f)\n",*considerB_flag,weight_12B,weight_12A,weight_12A,weight_1A,weight_2A,*total_time/pA->legtime,weight_12B,weight_1B,weight_2B,*total_time/pB->legtime);}}
    LA = llistmake(); pnW = (pnA==NULL ? pnB : pnA); while (pnW!=NULL){ if (LA->length==0){ litemadd(LA,pnW->region);} else /* if (LA->length>0) */{ l = liteminsertbefore(LA,LA->first); l->item = pnW->region;} pnW=pnW->parent;}
    if (verbose){ printf(" %% region list:"); l=LA->first; while (l!=NULL){ rg=(struct region *)l->item; printf(" %d",rg->label); l=l->child;} printf("\n");}
    l=LA->first; pnW = NULL; lW = pW->postree;
    while (l!=NULL){ 
      if ((lW2=llitemaddorfind(0,lW,(struct region *)l->item,&region2pnode_compare_label))==NULL){
	if (verbose){ printf(" %% region not found, making... ");}
	lW2=llitemaddorfind(1,lW,pnodemake(pnW,(struct region *)l->item,0,0),&pnode2pnode_compare_label);}
      else /* if found */{ if (verbose){ printf(" %% region found, not making... ");}}
      pnW = (struct pnode *) lW2->item; lW = pnW->childllitem;
      if (verbose){ printf("now at region %d\n",pnW->region->label);}
      l=l->child;}
    if (verbose){ printf(" %% finally at region %d, updating relevance %f\n",pnW->region->label,weight_A);}
    pnW->relevance = weight_A;
    llisttfree(LA);LA=NULL;}
  if (verbose){ printf(" %% [finished pnode_obsdisthist2obsdisthist_reweight]\n");}
}

void pnode2pnode_void_hist_operate(struct llist *LpA,struct ptree *pA,struct pnode *parentA,struct llitem *lA,struct llist *LpB,struct ptree *pB,struct pnode *parentB,struct llitem *lB,int maxlevel,int level,void (*void_hist_operate)(struct llist *,void *,struct llist *,void *,void *),void *void_parameters,int swap_flag)
{
  /* warning! never tested! fix later */
  /* LpA and LpB should be NULL to start with */
  int verbose=0;
  char lvlchar[64];
  int nl=0,firstentry=0;
  struct pnode *pnA=NULL,*pnB=NULL;
  struct llitem *lA2=NULL,*lB2=NULL;
  struct litem *l=NULL;
  struct pnode *pnC=NULL;
  if (verbose){ sprintf(lvlchar," %%"); for (nl=0;nl<level;nl++){ sprintf(lvlchar,"%s %%",lvlchar);}}
  if (verbose && parentA==NULL && lA->parent==NULL){ printf(" %s [entering pnode2pnode_void_hist_operate]\n",lvlchar);}
  if (verbose>2){ 
    l=LpA->first;while (l!=NULL){ pnC=(struct pnode *)l->item;printf(" %s LpA %d_r->%d\n",lvlchar,(int)pnC,pnC->region->label);l=l->child;}
    printf(" %% parentA %d, lA-> %d\n",(int)parentA,(lA!=NULL?(int)lA->item:-1)); pnodeprintf(parentA,lA,-1,0);
    l=LpB->first;while (l!=NULL){ pnC=(struct pnode *)l->item;printf(" %s LpB %d_r->%d\n",lvlchar,(int)pnC,pnC->region->label);l=l->child;}
    printf(" %% parentB %d, lB-> %d\n",(int)parentB,(lB!=NULL?(int)lB->item:-1)); pnodeprintf(parentB,lB,-1,0);}
  if (parentA!=NULL && lA->parent==NULL){ /* first descent into parent->childllitem */ 
    if (verbose>1){ printf(" %s first descent into parentA->childllitem, setting firstentry to 1\n",lvlchar);}
    assert(parentA->childllitem==lA); firstentry=1;}
  if (lA->kidl!=NULL){ 
    if (verbose>1){ printf(" %s stepping left to lA->kidl\n",lvlchar);}
    pnode2pnode_void_hist_operate(LpA,pA,parentA,lA->kidl,LpB,pB,parentB,lB,maxlevel,level,void_hist_operate,void_parameters,swap_flag);}
  if (lA->kidr!=NULL){ 
    if (verbose>1){ printf(" %s stepping right to lA->kidr\n",lvlchar);}
    pnode2pnode_void_hist_operate(LpA,pA,parentA,lA->kidr,LpB,pB,parentB,lB,maxlevel,level,void_hist_operate,void_parameters,swap_flag);}
  if (maxlevel==-1 || level<maxlevel){ 
    if (lA->item!=NULL){
      if (verbose>1){ printf(" %s level %d less than maxlevel %d, and lA->item exists, delving deeper\n",lvlchar,level,maxlevel);}
      pnA= (struct pnode *) lA->item; lA2 = (pnA==NULL ? NULL : pnA->childllitem);
      if (lB!=NULL && (lB2 = llitemaddorfind(0,lB,pB->regionra[pnA->region->label],&region2pnode_compare_label))!=NULL){
	pnB = (struct pnode *) lB2->item;}
      else{ pnB=NULL;}
      lB2 = (pnB==NULL ? NULL : pnB->childllitem);
      if (verbose>1){ 
	if (pnB!=NULL){ printf(" %s pnB label %d matches pnA label %d, ...",lvlchar,pnB->region->label,pnA->region->label);}
	else{ printf(" %s pnB label doesn't match pnA label %d, ...",lvlchar,pnA->region->label);}
	printf(" descending\n");}
      if (LpA==NULL){ LpA = llistmake();} litemadd(LpA,pnA); if (LpB==NULL){ LpB = llistmake();} litemadd(LpB,pnB);
      pnode2pnode_void_hist_operate(LpA,pA,pnA,lA2,LpB,pB,pnB,lB2,maxlevel,level+1,void_hist_operate,void_parameters,swap_flag);}
    else /* if (lA->item==NULL) */{
      if (verbose>1){ printf(" %s level %d less than maxlevel %d, and lA->item doesn't exist...",lvlchar,level,maxlevel);}
      pnA = NULL; lA2 = (pnA==NULL ? NULL : pnA->childllitem);
      if (lB!=NULL && lB->item!=NULL){
	if (verbose>1){ printf(" but lB->item does exist, swapping.\n");}
	pnB = (struct pnode *) lB->item;
	lB2 = (pnB==NULL ? NULL : pnB->childllitem);
	if (pnB!=NULL && lB2!=NULL){
	  if (LpA==NULL){ LpA = llistmake();} litemadd(LpA,pnA); if (LpB==NULL){ LpB = llistmake();} litemadd(LpB,pnB);
	  pnode2pnode_void_hist_operate(LpB,pB,pnB,lB2,LpA,pA,pnA,lA2,maxlevel,level+1,void_hist_operate,void_parameters,!swap_flag);}}
      else /* if (lB==NULL || lB->itime==NULL) */{ if (verbose>1){ printf(" lB->item doesn't exist either, stopping.\n");}}}}
  if (firstentry){
    if (verbose>1){ printf(" %s first entry for region %d,\n",lvlchar,parentA->region->label);}
    if (swap_flag){ void_hist_operate(LpB,parentB,LpA,parentA,void_parameters);}
    else /* if (!swap_flag) */{ void_hist_operate(LpA,parentA,LpB,parentB,void_parameters);}
    if (LpA!=NULL){ llistkilllast(LpA); if (LpA->length==0){ llisttfree(LpA); LpA=NULL;}}
    if (LpB!=NULL){ llistkilllast(LpB); if (LpB->length==0){ llisttfree(LpB); LpB=NULL;}}
    if (verbose>1){ if (LpA!=NULL){ l=LpA->first;while (l!=NULL){ pnC=(struct pnode *)l->item;printf(" %s after deletion: LpA %d_r->%d\n",lvlchar,(int)pnC,pnC->region->label);l=l->child;}} if (LpB!=NULL){ l=LpB->first;while (l!=NULL){ pnC=(struct pnode *)l->item;printf(" %s after deletion: LpB %d_r->%d\n",lvlchar,(int)pnC,pnC->region->label);l=l->child;}}}}
} 
