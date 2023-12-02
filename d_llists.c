/* /\* This holds a link-listed item *\/ */
/* struct litem */
/* { */
/*   struct litem *parent; */
/*   struct litem *child; */
/*   void *item; */
/* }; */

/* /\* This holds a link-list of litems *\/ */
/* struct llist */
/* { */
/*   struct litem *first; */
/*   struct litem *last; */
/*   int length; */
/* }; */

/* /\* This holds a mem-list item *\/ */
/* struct mitem */
/* { */
/*   struct mitem *parent; */
/*   struct mitem *child; */
/*   void *item; */
/*   size_t size; */
/* }; */

/* /\* This holds a list of mitems *\/ */
/* struct mlist */
/* { */
/*   struct mitem *first; */
/*   struct mitem *last; */
/*   int length; */
/*   size_t size; */
/* }; */

/* /\* Here are memory functions *\/ */
/* void meminit(); */
/* void memlink(void *,size_t); */
/* void memdrop(void *); */
/* void memprintf(int); */
/* void memsearch(void *); */
/* void * tmalloc(size_t); */
/* void * tcalloc(size_t,size_t); */
/* void tfree(void *); */

/* /\* Here are llist functions *\/ */
/* struct litem * litemmake(); */
/* void litemtfree(struct litem * ); */
/* struct llist * llistmake(); */
/* struct llist * llistcopy(); */
/* void llisttfree(struct llist *); */
/* void llisttfree2(struct llist *); */
/* void litemadd(struct llist *,void *); */
/* void litemminus(struct llist *,void *); */
/* void llistsort(struct litem *,struct litem *,int); */
/* void llistprintf(struct llist *); */
/* void llistunmake(struct llist **,struct llist **,struct llist *); */
/* void llistremake(struct llist ***,struct llist ***,struct llist **); */

/* /\* Here are the statistics functions *\/ */
/* void lliststats(struct llist *,int,double *); */
/* double llistcorrelation(struct llist *,struct llist *); */
/* void binstats(char *,void *,int,double *,double *); */
/* void stats(char *,void *,int,double *,double *,double *,double *); */
/* double correlation(double *,double *,int); */
/* double * spacesmear(double *,int,int); */

/* Here are the RNUM functions */

double randn()
{
  /* box muller polar form */
  double u=0,v=0,s=0;
  while (s==0 || s>1){ u=2*rand01-1;v=2*rand01-1;s=u*u+v*v;}
  return u*sqrt(-2*log(s)/s);
}

long long int RGET(long long int *rseed)
{
  /* basically:
     RNEXT = (RPREV*((long long int)(pow(2,RPOW)+RADD))%((long long int) pow(2,2*RPOW-1)));
     return RNEXT */
  *rseed = (*rseed*POW2RPOWPLUSRADD%POW22RPOWMINUSONE);
  return *rseed;
}

double R01GET(long long int *rseed)
{
 /* basically:
    RNEXT = (RPREV*((long long int)(pow(2,RPOW)+RADD))%((long long int) pow(2,2*RPOW-1)));
    return = (double)RNEXT/(double)pow(2,2*RPOW-1); */
  *rseed = (*rseed*POW2RPOWPLUSRADD%POW22RPOWMINUSONE);
  return (double)*rseed/(double)POW22RPOWMINUSONE;
}

double RNGET(long long int *rseed)
{
  /* box muller polar form */
  double u=0,v=0,s=0;
  while (s==0 || s>1){ u=2*R01GET(rseed)-1;v=2*R01GET(rseed)-1;s=u*u+v*v;}
  return u*sqrt(-2*log(s)/s);
}

double RISIGET(long long int *rseed,double rate)
{
  double r = R01GET(rseed);
  if (r==0){ r=1;}
  return -log(r)/rate;
}

/* Here are memory functions */

void meminit()
{
  struct mlist *M=NULL;
  M = (struct mlist *) malloc(sizeof(struct mlist));
  M->first=NULL;
  M->last=NULL;
  M->length=0;
  M->size=0;
  MEMDEBUGLIST=M;
}

void memlink(void *p,size_t size)
{
  struct mitem *m=NULL;
  struct mlist *M=NULL;
  M = (struct mlist *) MEMDEBUGLIST;
  m = (struct mitem *) malloc(sizeof(struct mitem));
  m->parent = NULL;
  m->child = NULL;
  m->item = p;
  m->size = size;
  if (M->length==0){
    M->first = m;
    M->last = m;
    M->length = 1;
    M->size = m->size;}
  else if (M->length > 0){
    M->last->child = m;
    m->parent = M->last;
    M->last = m;
    M->length++;
    M->size += m->size;}
}

void memdrop(void *p)
{
  int found=0,unfound=1;
  struct mitem *n=NULL;
  struct mlist *M=NULL;
  M = (struct mlist *) MEMDEBUGLIST;
  if (M->length == 0){ return;}
  else if (M->length == 1){
    n = M->first;
    if (n->item == p){
      M->first = NULL;
      M->last = NULL;
      M->length = 0;
      M->size = 0;
      free(n);}
    else{ printf(" %% error, cannot free %d\n",(int)p);}}
  else if (M->length > 1){
    found = 0;
    n = M->first;
    if (!found && n->item == p){ M->first = n->child; M->first->parent = NULL; M->length--; M->size -= n->size; free(n); found=1;}
    n = M->last;
    if (!found && n->item == p){ M->last = n->parent; M->last->child = NULL; M->length--; M->size -= n->size; free(n); found=1;}
    n = M->last; unfound =1;
    while (!found && unfound && n!=NULL){
      if (n->item == p){
	if (n->parent != NULL){ n->parent->child = n->child;}
	if (n->child != NULL){ n->child->parent = n->parent;}
	M->length--;
	M->size -= n->size;
	free(n);n=NULL;
	unfound = 0;
	found = 1;}
      if (!found && unfound){ n = n->parent;}}
    if (!found && unfound){ printf(" %% error, %d not in list\n",(int) p);}}
}

void memprintf(int imax)
{
  int i=0;
  struct mitem *m=NULL;
  struct mlist *M=NULL;
  M = (struct mlist *) MEMDEBUGLIST;
  printf("%% MEM:  %d items of total size %d\n",M->length,(int) M->size);
  m = M->last;
  i=0;
  while (i<imax && m!=NULL){
    printf(" %% mitem %d parent %d child %d item %d size %d\n",(int)m,(int)m->parent,(int)m->child,(int)m->item,(int)m->size);
    m = m->parent;
    i++;}
}

void memsearch(void *p)
{
  int i=0,found=0;
  struct mitem *m=NULL;
  struct mlist *M=NULL;
  M = (struct mlist *) MEMDEBUGLIST;
  m = M->last;
  i = 0;
  while (m!=NULL){
    if (m->item == p){ printf("found %d at memlist index %d\n",(int)p,M->length-i);found++;}
    m = m->parent;
    i++;}
  printf("%d instances of %d located\n",found,(int)p);
}

void * tmalloc(size_t size)
{
  void *p;
  p = malloc(size);
  if (MEMDEBUG){ memlink(p,size);}
  return p;
}

void * tcalloc(size_t nmemb,size_t size)
{
  void *p;
  p = calloc(nmemb,size);
  if (MEMDEBUG){ memlink(p,size*nmemb);}
  return p;
}

void * trealloc(void *p,size_t size)
{
  void *p2;
  p2 = realloc(p,size);
  if (MEMDEBUG){ memdrop(p); memlink(p2,size);}
  return p2;
}

void tfree(void *p)
{
  if (MEMDEBUG){ memdrop(p);}
  if (p==NULL){ printf(" %% Warning! tfree(NULL)\n");}
  free(p);
}

/* Here are llist functions */

struct litem * litemmake()
{
  struct litem *l;
  l = (struct litem *) tmalloc(sizeof(struct litem));
  l->parent = NULL;
  l->child = NULL;
  l->item = NULL;
  return l;
}

void litemtfree(struct litem * l)
{
  if (l->parent != NULL){ l->parent->child = l->child;}
  if (l->child != NULL){ l->child->parent = l->parent;}
  tfree(l);
  l = NULL;
}

struct llist * llistmake()
{
  struct llist *L;
  L = (struct llist *) tmalloc(sizeof(struct llist));
  L->first = NULL;
  L->last = NULL;
  L->length = 0;
  return L;
}

struct llist * llistcopy(struct llist *L)
{
  struct litem *l=NULL;
  struct llist *M=NULL;
  M = llistmake();
  if (L!=NULL){ l = L->first; while (l != NULL){ litemadd(M,l->item); l = l->child;}}
  return M;
}

struct llist * llistunion(struct llist *L1,struct llist *L2)
{
  /* creates a llist by adding together *L1 and *L2 */
  struct litem *l=NULL;
  struct llist *M=llistcopy(L1);
  if (L2!=NULL){ l = L2->first; while (l!=NULL){ litemadd(M,l->item); l = l->child;}}
  return M;
}

void llistgrowllist(struct llist *L1,struct llist *L2)
{
  /* grows llist *L1 by adding elements of *L2 */
  struct litem *l=NULL;
  if (L2!=NULL){ l = L2->first; while (l!=NULL){ litemadd(L1,l->item); l = l->child;}}
}

void llistprune(struct llist *L)
{
  /* gets rid of duplicate entries in L */
  struct litem *l=NULL,*l2=NULL;
  if (L!=NULL){
    if (L->length>1){
      llistsort(L->first,L->last,L->length,&void_compare);
      l = L->first;
      while (l!=NULL){
	l2=l->child;
	if (l2!=NULL){
	  if (l2->item==l->item){
	    if (l2->child==NULL){ litemtfree(l2); L->last = l;}
	    else /* if (l2->child!=NULL) */{ litemtfree(l2);}
	    L->length -= 1;}
	  else /* if (l2->item!=l->item) */{ l = l->child;}}
	else /* if (l2==NULL) */{ l = l->child;}}}}
}

void llisttfree(struct llist *L)
{
  /* this tfrees every item of the llist */
  struct litem *l;
  if (L->length == 0){
    tfree(L);
    L=NULL;}
  else if (L->length == 1){ 
    l = L->first;
    litemtfree(l);
    l = NULL;
    tfree(L);
    L=NULL;}
  else if (L->length > 1){
   l = L->last;
   L->last = l->parent;
   litemtfree(l);
   L->length--;
   assert(L->last->child == NULL);
   llisttfree(L);}
}

void llisttfree2(struct llist *L)
{
  /* this tfrees the llist and items and (double *) elements */
  struct litem *l;
  if (L->length == 0){
    tfree(L);
    L=NULL;}
  else if (L->length == 1){ 
    l = L->first;
    tfree((double *)(l->item)); l->item=NULL;
    litemtfree(l);
    l = NULL;
    tfree(L);
    L=NULL;}
  else if (L->length > 1){
   l = L->last;
   L->last = l->parent;
   tfree((double *)(l->item)); l->item=NULL;
   litemtfree(l);
   L->length--;
   assert(L->last->child == NULL);
   llisttfree2(L);}
}

void llisttfree3(struct llist *L)
{
  /* this tfrees the llist and items and (int *) elements */
  struct litem *l;
  if (L->length == 0){
    tfree(L);
    L=NULL;}
  else if (L->length == 1){ 
    l = L->first;
    tfree((int *)(l->item)); l->item=NULL;
    litemtfree(l);
    l = NULL;
    tfree(L);
    L=NULL;}
  else if (L->length > 1){
   l = L->last;
   L->last = l->parent;
   tfree((int *)(l->item)); l->item=NULL;
   litemtfree(l);
   L->length--;
   assert(L->last->child == NULL);
   llisttfree3(L);}
}

int isin(struct llist *L,void *n)
{
  int found_flag=0;
  struct litem *l=NULL;
  l=L->first;
  while(l!=NULL){ found_flag += (l->item==n); l=l->child;} 
  return found_flag;
}

void litemadd(struct llist *L,void *n)
{
  struct litem *l=litemmake();
  l->item = n;
  if (L->length == 0){
    L->first = l;
    L->last = l;
    L->length = 1;}
  else if (L->length > 0){
    L->last->child = l;
    l->parent = L->last;
    L->last = l;
    L->length++;}
}

int litemexadd(struct llist *L,void *n)
{
  /* adds the item *n to llist *L if *n is not already contained in *L */
  struct litem *l=NULL;
  struct litem *l2=NULL;
  int found_flag=0;
  if (L->length == 0){
    l = litemmake(); l->item = n;
    L->first = l;
    L->last = l;
    L->length = 1;}
  else if (L->length > 0){
    l2 = L->first; found_flag=0;
    while (l2!=NULL){ found_flag += (l2->item==n); l2=l2->child;}
    if (found_flag){ /* do nothing */}
    else if (!found_flag){
      l = litemmake(); l->item = n;
      L->last->child = l;
      l->parent = L->last;
      L->last = l;
      L->length++;}}
  return found_flag;
}

struct litem * liteminsertbefore(struct llist *L,struct litem *l)
{
  struct litem *l2=litemmake();
  if (L->first == l){
    l2->parent = NULL;
    l2->child = l;
    l->parent = l2;
    L->first=l2;}
  else if (L->first != l){
    l2->parent = l->parent;
    l2->child = l;
    l->parent->child = l2;
    l->parent = l2;}
  L->length += 1;
  return l2;
}

struct litem * liteminsertafter(struct llist *L,struct litem *l)
{
  struct litem *l2=litemmake();
  if (L->last == l){
    l2->parent = l;
    l2->child = NULL;
    l->child = l2;
    L->last=l2;}
  else if (L->last != l){
    l2->parent = l;
    l2->child = l->child;
    l->child->parent = l2;
    l->child = l2;}
  L->length += 1;
  return l2;
}

void llistkillfirst(struct llist *L)
{
  /* simply kills the first element of a llist */
  struct litem *l=NULL;
  if (L->length == 0){ return;}
  else if (L->length==1){ litemtfree(L->first); L->first=NULL; L->last=NULL; L->length=0;}
  else if (L->length>1){ l=L->first; L->first = l->child; L->length -= 1; litemtfree(l); l=NULL;}
}

void llistkillfirst2(struct llist *L)
{
  /* kills the first element of a llist, assuming that it is a (double *) element */
  struct litem *l=NULL;
  if (L->length == 0){ return;}
  else if (L->length==1){ 
    tfree((double *)(L->first->item)); L->first->item=NULL;
    litemtfree(L->first); L->first=NULL; L->last=NULL; L->length=0;}
  else if (L->length>1){ 
    l=L->first; L->first = l->child; L->length -= 1; 
    tfree((double *)(l->item)); l->item=NULL;
    litemtfree(l); l=NULL;}
}

void llistkilllast(struct llist *L)
{
  /* simply kills the last element of a llist */
  struct litem *l=NULL;
  if (L->length == 0){ return;}
  else if (L->length==1){ litemtfree(L->last); L->first=NULL; L->last=NULL; L->length=0;}
  else if (L->length>1){ l=L->last; L->last = l->parent; L->length -= 1; litemtfree(l); l=NULL;}
}

void llistkilllast2(struct llist *L)
{
  /* kills the last element of a llist, assuming that it is a (double *) element */
  struct litem *l=NULL;
  if (L->length == 0){ return;}
  else if (L->length==1){ 
    tfree((double *)(L->first->item)); L->first->item=NULL;
    litemtfree(L->last); L->first=NULL; L->last=NULL; L->length=0;}
  else if (L->length>1){ 
    l=L->last; L->last = l->parent; L->length -= 1; 
    tfree((double *)(l->item)); l->item=NULL;
    litemtfree(l); l=NULL;}
}

void litemminus(struct llist *L,void *n)
{
  /* this assumes L is a llist of distinct objects */
  int done=0;
  struct litem *l=NULL;
  if (L->length == 0){ return;}
  else if (L->length == 1){ if (L->first->item==n){ litemtfree(L->first); L->first=NULL; L->last=NULL; L->length=0;}}
  else if (L->length > 1){
    done=0;
    l=L->first;
    if (!done && l->item==n){ L->first = l->child; L->first->parent=NULL; litemtfree(l); L->length--; done=1;}
    l=L->last;
    if (!done && l->item==n){ L->last = l->parent; L->last->child=NULL; litemtfree(l); L->length--; done=1;}
    l=L->last;
    while (!done && l!=NULL){
      if (l->item==n){
	l->parent->child = l->child;
	l->child->parent = l->parent;
	litemtfree(l);
	L->length--;
	done=1;}
      l = l->parent;}}
}

void llistprintf(struct llist *L)
{
  struct litem *l=NULL;
  if (L==NULL){ printf("NULL list\n");}
  else if (L->length == 0){ printf("L=%d with first %d, last %d and length %d\n",(int)L,(int)L->first,(int)L->last,L->length);}
  else if (L->length >= 1){
    printf("L=%d with first %d, last %d and length %d\n",(int)L,(int)L->first,(int)L->last,L->length);
    l = L->first;
    while (l != NULL){
      printf("   %d with parent %d, child %d and item %d\n",(int)l,(int)l->parent,(int)l->child,(int)l->item);
      l = l->child;}}
}

void llistprintf2(struct llist *L)
{
  struct litem *l=NULL;
  if (L==NULL){ printf("NULL list\n");}
  else if (L->length == 0){ printf(" empty list\n");}
  else if (L->length >= 1){
    l=L->first; while (l!=NULL){ printf("%f ",*(double *)l->item); l=l->child;} printf("\n");}
}

int ra2ra_generic_compare_dictionary(void *vra1,void *vra2,void *void_parameters)
{
  /* parameter list
     (int *) length
     (size_t *) varsize
     (int)(*compare)(void *,void *)
  */
  int verbose=0;
  int *length=NULL;
  size_t *varsize=NULL;
  int (*compare)(void *,void *) = NULL;
  void **vra=NULL;
  int comparison_flag=0;
  int nl=0;
  if (verbose){ printf(" %% [entering ra2ra_generic_compare_dictionary] with input addresses %d,%d\n",(int)vra1,(int)vra2);}
  if (void_parameters!=NULL){
    vra = (void **) void_parameters;
    length = vra[0];
    varsize = vra[1];
    if (verbose){ printf(" %% read length %d, size %d\n",(int)(*length),(int)(*varsize));}
    compare = vra[2];
    comparison_flag=0;
    nl=0;
    while (comparison_flag==0 && nl<*length){
      if (verbose){ printf(" %% comparing memory addresses %d+%d and %d+%d\n",(int)vra1,(int)nl*(*varsize),(int)vra2,(int)nl*(*varsize));}
      comparison_flag = (*compare)(vra1+(*varsize)*nl,vra2+(*varsize)*nl);
      if (verbose){ printf(" %% %% found comparison_flag %d at nl %d\n",comparison_flag,nl);}
      nl+=1;}}
  return comparison_flag;
}

int double_compare(void *n1,void *n2)
{
  /* if *(double *)n1 > *(double *)n2 */
  double *d1=(double *)n1,*d2=(double *)n2;
  return (*d1>*d2) ? 1 : (*d1<*d2) ? -1 : 0;
}

int int_compare(void *n1,void *n2)
{
  /* if *(int *)n1 > *(int *)n2 */
  int *i1=(int *)n1,*i2=(int *)n2;
  return (*i1>*i2) ? 1 : (*i1<*i2) ? -1 : 0;
}

int void_compare(void *n1,void *n2)
{
  /* if n1 > n2 */
  return (n1>n2) ? 1 : (n1<n2) ? -1 : 0;
}

void llistsort(struct litem *lfirst,struct litem *llast,int length,int (*compare)(void *,void *))
{
  /* this quicksorts a list of (unique) litems in order of increasing compare */
  int pivot=0,itemindex=0;
  int next=0,epl=0,epr=0;
  int comparison=0;
  struct litem *lpivot=NULL;
  struct litem *lnext=NULL;
  struct litem *lepl=NULL;
  struct litem *lepr=NULL;
  void *ipivot=NULL,*itmp=NULL,*inext=NULL;
  if (length <= 1){ return;}
  else if (length > 1){
    pivot = rand()%length;
    if (pivot <= length/2){
      lpivot = lfirst;
      for (itemindex=0;itemindex<pivot;itemindex++){ lpivot = lpivot->child;}}
    else if (pivot > length/2){ 
      lpivot = llast;
      for (itemindex=length;itemindex>pivot;itemindex--){ lpivot = lpivot->parent;}}
    ipivot = (void *) lpivot->item; /* this may change subsequently */
    next=0;epl=0;epr=length;
    lnext=lfirst;lepl=lfirst;lepr=llast;
    while (next < epr){
      comparison = 0;
      inext = (void *) lnext->item;
      comparison = (*compare)(inext,ipivot);
      if (comparison < 0){
	itmp = (void *) lepl->item;
	lepl->item = lnext->item;
	lnext->item = itmp;
	epl++;
	next++;
	lepl = lepl->child;
	lnext = lnext->child;}
      else if (comparison > 0){
	itmp = (void *) lepr->item;
	lepr->item = lnext->item;
	lnext->item = itmp;
	epr--;
	lepr = lepr->parent;}
      else if (comparison == 0){
	next ++;
	lnext = lnext->child;}}
    assert(0<=epl && epl<epr && epr<=length);
    llistsort(lfirst,lepl,epl+1,compare);
    if (epr < length){ lepr = lepr->child; llistsort(lepr,llast,length-epr,compare);}}
}

struct llist ** llistra2read(char *filename,int *nrecords)
{
  /* assume *filename dumped with radump, and starts with "./" */
  FILE *fp=NULL;
  struct llist **Lra=NULL;
  int nr=0,nl=0,length=0;
  double *temp=NULL;
  if ((fp=fopen(filename,"r"))==NULL){ printf(" %% Warning, cannot read %s in llistra2read\n",filename);}
  else{
    fread(nrecords,sizeof(int),1,fp);
    Lra = (struct llist **) tcalloc(*nrecords,sizeof(struct llist *));
    for (nr=0;nr<*nrecords;nr++){ 
      Lra[nr]=llistmake();
      fread(&length,sizeof(int),1,fp);
      for (nl=0;nl<length;nl++){
	temp = (double *) tcalloc(1,sizeof(double));
	fread(temp,sizeof(double),1,fp);
	litemadd(Lra[nr],temp);}}
    fclose(fp);fp=NULL;
    return Lra;}
  return NULL;
}

int llistra2dump(struct llist **Lra,int total_records,char *filename)
{
  /* assume filename starts with "./" 
     assumes Lra[nr] is a llist of (double *) items */
  int output_flag=1,nr=0;
  FILE *fp=NULL;
  struct litem *l=NULL;
  if ((fp=fopen(filename,"w"))==NULL){ printf(" %% Warning, cannot open %s in llistradump, writing to stdout\n",filename); fp=stdout;}
  output_flag *= fwrite(&total_records,sizeof(int),1,fp);
  for (nr=0;nr<total_records;nr++){
    output_flag *= fwrite(&(Lra[nr]->length),sizeof(int),1,fp);
    l=Lra[nr]->first;
    while (l!=NULL){
      output_flag *= fwrite((double *)l->item,sizeof(double),1,fp);
      l=l->child;}}
  if (fp!=stdout){ fclose(fp);}
  return output_flag;
}

/* Here are the llitem functions */

struct llitem * llitemmake()
{
  struct llitem *l=NULL;
  l = (struct llitem *) tmalloc(sizeof(struct llitem));
  l->parent = NULL;
  l->kidl = NULL;
  l->kidr = NULL;
  l->nleft = 0;
  l->nright = 0;
  l->item = NULL;
  return l;
}

void llitemtfree(struct llitem *l,void (*free_function)(void *))
{
  if (free_function!=NULL && l->item!=NULL){ (*free_function)(l->item); l->item=NULL;}
  if (l->kidl!=NULL){ llitemtfree(l->kidl,free_function);}
  if (l->kidr!=NULL){ llitemtfree(l->kidr,free_function);}
  tfree(l); l=NULL;
}

struct llitem * llitemaddorfind(int addorfind,struct llitem *l,void *item,int (*compare)(void *,void *))
{
  /* addorfind = 1 --> add
     addorfind = 0 --> find */
  int c_self=0,c_parent=0;
  int c_down=0;
  if (l!=NULL && l->item!=NULL){
    c_down=1;
    c_self=(*compare)(item,l->item);
    if (l->parent!=NULL){
      c_parent=(*compare)(item,l->parent->item);
      if ((l==l->parent->kidl && c_self>=0 && c_parent>=0) || (l==l->parent->kidr && c_self<=0 && c_parent<=0)){ c_down=0; return llitemaddorfind(addorfind,l->parent,item,compare);}}
    if (c_down){
      if (c_self==0){ /* selfsame item */ if (addorfind){ printf(" %% error! duplicating element in llitemaddorfind\n");} return l; }
      else if (c_self<0){ /* lesser item */
	if (addorfind){ l->nleft += 1;}
	if (l->kidl!=NULL){ return llitemaddorfind(addorfind,l->kidl,item,compare);}
	else{ if (addorfind){ l->kidl=llitemmake(); l->kidl->parent=l; l->kidl->item=item; return l->kidl;} else{ return NULL;}}}
      else if (c_self>0){ /* greater item */
	if (addorfind){ l->nright += 1;}
	if (l->kidr!=NULL){ return llitemaddorfind(addorfind,l->kidr,item,compare);}
	else{ if (addorfind){ l->kidr=llitemmake(); l->kidr->parent=l; l->kidr->item=item; return l->kidr;} else{ return NULL;}}}}}
  else if (l!=NULL && l->item==NULL){
    if (addorfind){ l->item = item; return l;}}
  else if (l==NULL){ return NULL;}
  return NULL;
}

struct llitem * llitemaddorfind_generic(int addorfind,struct llitem *l,void *item,int (*void_compare)(void *,void *,void *),void *void_parameters)
{
  /* addorfind = 1 --> add
     addorfind = 0 --> find */
  int c_self=0,c_parent=0;
  int c_down=0;
  if (l!=NULL && l->item!=NULL){
    c_down=1;
    c_self=(*void_compare)(item,l->item,void_parameters);
    if (l->parent!=NULL){
      c_parent=(*void_compare)(item,l->parent->item,void_parameters);
      if ((l==l->parent->kidl && c_self>=0 && c_parent>=0) || (l==l->parent->kidr && c_self<=0 && c_parent<=0)){ c_down=0; return llitemaddorfind_generic(addorfind,l->parent,item,void_compare,void_parameters);}}
    if (c_down){
      if (c_self==0){ /* selfsame item */ if (addorfind){ printf(" %% error! duplicating element in llitemaddorfind_generic\n");} return l; }
      else if (c_self<0){ /* lesser item */
	if (addorfind){ l->nleft += 1;}
	if (l->kidl!=NULL){ return llitemaddorfind_generic(addorfind,l->kidl,item,void_compare,void_parameters);}
	else{ if (addorfind){ l->kidl=llitemmake(); l->kidl->parent=l; l->kidl->item=item; return l->kidl;} else{ return NULL;}}}
      else if (c_self>0){ /* greater item */
	if (addorfind){ l->nright += 1;}
	if (l->kidr!=NULL){ return llitemaddorfind_generic(addorfind,l->kidr,item,void_compare,void_parameters);}
	else{ if (addorfind){ l->kidr=llitemmake(); l->kidr->parent=l; l->kidr->item=item; return l->kidr;} else{ return NULL;}}}}}
  else if (l!=NULL && l->item==NULL){
    if (addorfind){ l->item = item; return l;}}
  else if (l==NULL){ return NULL;}
  return NULL;
}

struct llitem * llitemcoup(int rightorleft,struct llitem *l)
{
  /* if rightorleft=0, places l under l->kidl and returns l->kidl
     if rightorleft=1, places l under l->kidr and returns l->kidr */
  int nadd=0;
  struct llitem *ltemp=NULL,*lkidl=l->kidl,*lkidr=l->kidr;
  if (rightorleft==0){
    assert(lkidl!=NULL);
    nadd = l->nright + 1;
    ltemp = lkidl;
    while (ltemp->kidr!=NULL){
      ltemp->nright += nadd;
      ltemp = ltemp->kidr;}
    ltemp->nright += nadd;
    ltemp->kidr = l;
    lkidl->parent = l->parent;
    l->parent = ltemp;
    l->kidl = NULL; 
    l->nleft = 0;
    if (lkidl->parent !=NULL && lkidl->parent->kidl==l){ lkidl->parent->kidl=lkidl;}
    if (lkidl->parent !=NULL && lkidl->parent->kidr==l){ lkidl->parent->kidr=lkidl;}
    return lkidl;}
  else if (rightorleft==1){
    assert(lkidr!=NULL);
    nadd = l->nleft + 1;
    ltemp = lkidr;
    while (ltemp->kidl!=NULL){
      ltemp->nleft += nadd;
      ltemp = ltemp->kidl;}
    ltemp->nleft += nadd;
    ltemp->kidl = l;
    lkidr->parent = l->parent;
    l->parent = ltemp;
    l->kidr = NULL;
    l->nright = 0;
    if (lkidr->parent !=NULL && lkidr->parent->kidl==l){ lkidr->parent->kidl=lkidr;}
    if (lkidr->parent !=NULL && lkidr->parent->kidr==l){ lkidr->parent->kidr=lkidr;}
    return lkidr;}
  return NULL;
}

void llitembalance(struct llitem *l)
{
  /* balances llitem *l */
  int enough=2;
  if (l->nleft>enough && l->nleft>enough*l->nright){ llitembalance(llitemcoup(0,l));}
  else if (l->nright>enough && l->nright>enough*l->nleft){ llitembalance(llitemcoup(1,l));}
  else{ if (l->kidl!=NULL){ llitembalance(l->kidl);} if (l->kidr!=NULL){ llitembalance(l->kidr);}}
}

struct llitem * llitemclimb(struct llitem *l){ /* returns head */ while (l->parent!=NULL){ l=l->parent;} return l;}

void llitemcheck(int verbose,struct llitem *l,int (*compare)(void *,void *))
{
  if (verbose){ printf(" %% l=%d, l->parent=%d, l->kidl=%d, l->kidr=%d, l->nleft=%d, l->nright=%d, l->item=%d\n",(int)l,(int)l->parent,(int)l->kidl,(int)l->kidr,l->nleft,l->nright,(int)l->item);}
  if (l->kidl!=NULL){ 
    assert(l->nleft == l->kidl->nleft + l->kidl->nright + 1); 
    assert((*compare)(l->kidl->item,l->item)<0); 
    llitemcheck(verbose,l->kidl,compare);}
  else{ assert(l->nleft == 0);}
  if (l->kidr!=NULL){ 
    assert(l->nright = l->kidr->nleft + l->kidr->nright + 1);
    assert((*compare)(l->item,l->kidr->item)<0); 
    llitemcheck(verbose,l->kidr,compare);}
  else{ assert(l->nright == 0);}
}

void llitemprintf(struct llitem *l,void (*vprintf)(void *))
{
  if (l->kidl!=NULL){ llitemprintf(l->kidl,vprintf);}
  if (vprintf==NULL){ printf("%d\n",(int)l->item);} else{ vprintf(l->item);}
  if (l->kidr!=NULL){ llitemprintf(l->kidr,vprintf);}
}

int llitemlength(struct llitem *l){ return (l==NULL? 0 : (l->nleft + l->nright + (l->item!=NULL)));}

void llist2llitem(struct llist *L,struct llitem *l,int (*compare)(void *,void *))
{
  /* produces llitem *l from llist *L given *compare */
  struct litem *ltemp=NULL;
  if (L->length>0){
    ltemp = L->first;
    l->item = ltemp->item;
    ltemp=ltemp->child;
    while (ltemp!=NULL){
      if (llitemaddorfind(0,l,ltemp->item,compare)==NULL){ llitemaddorfind(1,l,ltemp->item,compare);}
      ltemp=ltemp->child;}}
  /* l needs to be balanced and climbed */
}

void llitem2llist(struct llitem *l,struct llist *L)
{
  /* produces llist *L from llitem *l
     if *l is sorted, then *L should be sorted */
  if (l!=NULL){
    if (l->kidl!=NULL){ llitem2llist(l->kidl,L);}
    litemadd(L,l->item);
    if (l->kidr!=NULL){ llitem2llist(l->kidr,L);}}
}

void llitemgrowllitem(struct llitem *l1,struct llitem *l2,int (*compare)(void *,void *))
{
  /* adds *l2 to *l1, given *compare */
  if (l2->item!=NULL && llitemaddorfind(0,l1,l2->item,compare)==NULL){ llitemaddorfind(1,l1,l2->item,compare);}
  if (l2->kidl!=NULL){ llitemgrowllitem(l1,l2->kidl,compare);}
  if (l2->kidr!=NULL){ llitemgrowllitem(l1,l2->kidr,compare);}
  /* l1 needs to be balanced and climbed */
}

void llistgrowllitem(struct llist *L,struct llitem *l2)
{
  if (l2->kidl!=NULL){ llistgrowllitem(L,l2->kidl);}
  if (l2->item!=NULL){ litemadd(L,l2->item);}
  if (l2->kidr!=NULL){ llistgrowllitem(L,l2->kidr);}
}

/* Here are the histogram functions */

struct hist * histmake(int nbins,double max,double min)
{
  struct hist *h=NULL;
  h = (struct hist *) tmalloc(sizeof(struct hist));
  h->nbins = nbins;
  h->max = maximum(max,min);
  h->min = minimum(max,min);
  if (h->nbins>0){ h->data = (double *) tcalloc(h->nbins,sizeof(double));}
  else /* if (h->nbins==0) */{ h->data=NULL;}
  h->datasum = 0;
  h->ndatum = 0;
  return h;
}

void hist2file(void *void_hist,FILE *fp)
{
  int plusone=+1,minusone=-1;
  struct hist *h=NULL;
  if (void_hist==NULL){ fwrite(&minusone,sizeof(int),1,fp);}
  else /* if (void_hist!=NULL) */{
    fwrite(&plusone,sizeof(int),1,fp);
    h=(struct hist *)void_hist;
    fwrite(&(h->nbins),sizeof(int),1,fp);
    fwrite(&(h->max),sizeof(double),1,fp);
    fwrite(&(h->min),sizeof(double),1,fp);
    if (h->nbins>0){ fwrite(&(h->data[0]),sizeof(double),h->nbins,fp);}
    fwrite(&(h->datasum),sizeof(double),1,fp);
    fwrite(&(h->ndatum),sizeof(double),1,fp);}
}

void * file2hist(FILE *fp,int *no_error_passback)
{
  struct hist *h=NULL;
  int label=0,no_error=0;
  no_error=1;
  if (no_error){ no_error = (fread(&label,sizeof(int),1,fp)==(size_t)1) && !feof(fp) && !ferror(fp);}
  if (no_error && label==+1){
    h = (struct hist *) tcalloc(1,sizeof(struct hist));
    if (no_error){ no_error = (fread(&(h->nbins),sizeof(int),1,fp)==(size_t)1) && !feof(fp) && !ferror(fp);}
    if (no_error){ no_error = (fread(&(h->max),sizeof(double),1,fp)==(size_t)1) && !feof(fp) && !ferror(fp);}
    if (no_error){ no_error = (fread(&(h->min),sizeof(double),1,fp)==(size_t)1) && !feof(fp) && !ferror(fp);}
    if (no_error){
      if (h->nbins>0){
	h->data = (double *)tcalloc(h->nbins,sizeof(double));
	if (no_error){ no_error = (fread(&(h->data[0]),sizeof(double),h->nbins,fp)==(size_t)(h->nbins)) && !feof(fp) && !ferror(fp);}}
      else /* if (h->nbins==0) */{ h->data=NULL;}}
    if (no_error){ no_error = (fread(&(h->datasum),sizeof(double),1,fp)==(size_t)1) && !feof(fp) && !ferror(fp);}
    if (no_error){ no_error = (fread(&(h->ndatum),sizeof(double),1,fp)==(size_t)1) && !feof(fp) && !ferror(fp);}
    if (!no_error){
      printf(" %% warning! improperly read hist in file2hist\n");
      histtfree(h);h=NULL;}}
  if (no_error_passback!=NULL){ *no_error_passback=no_error;}
  return (void *)h;
}

void histtfree(struct hist *h){ if (h->data!=NULL){ tfree(h->data);} h->data=NULL; tfree(h); h=NULL;}

struct hist * histcopy(struct hist *h1)
{
  struct hist *h2=NULL;
  int nb=0;
  h2 = (struct hist *) tmalloc(sizeof(struct hist));
  h2->nbins = h1->nbins;
  h2->max = h1->max;
  h2->min = h1->min;
  if (h2->nbins>0){ 
    h2->data = (double *) tmalloc(h2->nbins*sizeof(double)); 
    for (nb=0;nb<h2->nbins;nb++){ h2->data[nb] = h1->data[nb];}}
  else /* if (h2->nbins==0) */{ h2->data=NULL;}
  h2->datasum = h1->datasum;
  h2->ndatum = h1->ndatum;
  return h2;
}

void histadd(struct hist *h,double val,double ntimes)
{
  int bin=0;
  if (finite(ntimes) && finite(val) && ntimes>0){
    bin = (int)floor(h->nbins*(val-h->min)/(h->max-h->min));
    bin = maximum(0,minimum(bin,h->nbins-1));
    if (h->data!=NULL){
      h->data[bin] += ntimes;
      h->datasum += ntimes*val;
      h->ndatum += ntimes;}}
}

void histprintf(struct hist *h,char *prefix)
{
  int bin=0;
  printf(" %s nbins %d, max %f, min %f, datasum %f, ndatum %f\n",prefix,h->nbins,h->max,h->min,h->datasum,h->ndatum);
  if (h->nbins<1){ printf(" %s empty!\n",prefix);}
  for (bin=0;bin<h->nbins;bin++){ printf(" %s bin %d(%f)\t%f\n",prefix,bin,h->min+bin*(h->max-h->min)/(double)(h->nbins-1),h->data[bin]);}
}

void llist2hist(struct llist *L,struct hist *h)
{
  struct litem *l=NULL;
  l=L->first;
  while (l!=NULL){
    histadd(h,*(double *)l->item,1);
    l=l->child;}
}

void histdump(struct hist *h,int dump_type,char *filename,char *prefix,int logflag)
{
  /* dumps a histogram to ("%s.m",filename)
     assumes filename starts with "./" */
  int bin=0;
  char filename2[1024];
  FILE *fp=NULL;
  if (dump_type==0){
    sprintf(filename2,"%s.m",filename);
    if (filename2==NULL || (fp=fopen(filename2,"w"))==NULL){ printf("histogram dumped to stdout\n"); fp=stdout;}
    fprintf(fp,"%% %shistogram dump of %s\n",prefix,filename);
    fprintf(fp,"%% max %0.16f min %0.16f nbins %d\n",h->max,h->min,h->nbins);
    fprintf(fp,"%shist_datasum = %0.16lf;\n",prefix,h->datasum);
    fprintf(fp,"%shist_ndatum = %0.16lf;\n",prefix,h->ndatum);
    for (bin=0;bin<h->nbins;bin++){
      fprintf(fp,"%shist_data(%d) = %0.16lf;\n",prefix,bin+1,h->data[bin]);}
    fprintf(fp,"plot(%shist_data);\n",prefix);
    if (fp!=stdout){ fclose(fp);}}
  else if (dump_type==1){ if (h->ndatum>0){ ra2jpg(h->data,"double",1,h->nbins,0,filename,logflag);}}
}

/* Here are the statistics functions */

double znorm(double *z1r,double *z1i){ return sqrt(pow(*z1r,2)+pow(*z1i,2));}
void zpz(double *z1r,double *z1i,double *z2r,double *z2i,double *z3r,double *z3i){
  double tr=0,ti=0; tr=*z1r+*z2r; ti=*z1i+*z2i; if(z3r!=NULL){*z3r=tr;}if(z3i!=NULL){*z3i=ti;}}  
void zmz(double *z1r,double *z1i,double *z2r,double *z2i,double *z3r,double *z3i){
  double tr=0,ti=0; tr=*z1r-*z2r; ti=*z1i-*z2i; if(z3r!=NULL){*z3r=tr;}if(z3i!=NULL){*z3i=ti;}}  
void zxz(double *z1r,double *z1i,double *z2r,double *z2i,double *z3r,double *z3i){
  double tr=0,ti=0; tr=*z1r**z2r-*z1i**z2i; ti=*z1r**z2i+*z1i**z2r; if(z3r!=NULL){*z3r=tr;}if(z3i!=NULL){*z3i=ti;}}
void zdz(double *z1r,double *z1i,double *z2r,double *z2i,double *z3r,double *z3i){
  double tr=0,ti=0,tn=*z2r**z2r+*z2i**z2i; if (tn==0){tn=1;}
  tr=*z1r**z2r+*z1i**z2i; ti=*z2r**z1i-*z1r**z2i; if(z3r!=NULL){*z3r=tr/tn;}if(z3i!=NULL){*z3i=ti/tn;}}
void zpeval(int length,double *coef_real,double *coef_imag,double *zr,double *zi,double *z3r,double *z3i){
  int nr=0; double tr=coef_real[0],ti=coef_imag[0];
  for (nr=1;nr<length;nr++){ zxz(&tr,&ti,zr,zi,&tr,&ti); zpz(&tr,&ti,coef_real+nr,coef_imag+nr,&tr,&ti);}
  if (z3r!=NULL){ *z3r=tr;} if (z3i!=NULL){ *z3i=ti;}}

int durand_kerner_rootfind(int length,double *coef_real,double *coef_imag,double *root_real,double *root_imag)
{
  /* simple durand_kerner rootfinding, no recursion, no analytic acceleration */
  int verbose=0;
  double *norm=NULL,ubound=0,zbound=0,relative_error_threshold=0.00001;
  double *temp=NULL,*temp2=NULL,*temp3=NULL;
  int nz=0,nz2=0;
  int iteration=0,iteration_max=2*(int)pow(length,2);
  int error_flag=0;
  norm = (double *) tcalloc(length,sizeof(double));
  temp = (double *) tcalloc(2,sizeof(double));
  temp2 = (double *) tcalloc(2,sizeof(double));
  temp3 = (double *) tcalloc(2,sizeof(double));
  for (nz=1;nz<length;nz++){ ubound = maximum(ubound,pow(znorm(coef_real+nz,coef_imag+nz)/znorm(coef_real,coef_imag),1.0/(length-nz)));}
  if (verbose>0){ printf(" starting norm %f\n",ubound);}
  for (nz=0;nz<length-1;nz++){ 
    norm[nz]=ubound;
    root_real[nz] = ubound*cos(2*PI*(double)(nz+0.5)/(double)(length-1));
    root_imag[nz] = ubound*sin(2*PI*(double)(nz+0.5)/(double)(length-1));}
  iteration=0;
  do{
    if (verbose>0){ printf(" iteration %d, ",iteration);}
    if (verbose>0){ printf("roots: "); for (nz=0;nz<length-1;nz++){ printf(" %f+%fi ",root_real[nz],root_imag[nz]);} printf("\n");}
    norm[length-1]=0;
    if (verbose>0){ printf("update: ");}
    for (nz=0;nz<length-1;nz++){
      zpeval(length,coef_real,coef_imag,root_real+nz,root_imag+nz,temp,temp+1);
      if (verbose>1){ printf("(%f+%fi)",temp[0],temp[1]);}
      temp2[0]=1;temp2[1]=0;
      for (nz2=0;nz2<length-1;nz2++){
	if (nz2!=nz){ 
	  zmz(root_real+nz,root_imag+nz,root_real+nz2,root_imag+nz2,temp3,temp3+1); 
	  zxz(temp2,temp2+1,temp3,temp3+1,temp2,temp2+1);}}
      if (verbose>1){ printf("/(%f+%fi) ",temp2[0],temp2[1]);}
      zdz(temp,temp+1,temp2,temp2+1,temp3,temp3+1);
      if (verbose>0){ printf(" %f+%fi ",temp3[0],temp3[1]);}
      norm[nz] = znorm(temp3,temp3+1)/znorm(root_real+nz,root_imag+nz);
      zmz(root_real+nz,root_imag+nz,temp3,temp3+1,root_real+nz,root_imag+nz);
      zbound=znorm(root_real+nz,root_imag+nz); if (zbound>ubound){ root_real[nz]*=ubound/zbound; root_imag[nz]*=ubound/zbound;}
      norm[length-1]=maximum(norm[length-1],norm[nz]);}
    if (verbose>0){ printf("\n");}
    iteration += 1;}
  while (iteration<iteration_max && norm[length-1]>relative_error_threshold);
  if (iteration>=iteration_max || norm[length-1]>=relative_error_threshold){
    if (verbose>-1){ printf(" warning, improper convergence in durand_kerner_rootfind, norm %f\n",norm[length-1]);} error_flag=1;}
  else /* if (iteration<iteration_max && norm[length-1]<relative_error_threshold) */{
    if (verbose>0){ printf(" successful output, iteration %d norm %f\n",iteration,norm[length-1]);} error_flag=0;}
  tfree(norm);norm=NULL;
  tfree(temp);temp=NULL;
  tfree(temp2);temp2=NULL;
  tfree(temp3);temp3=NULL;
  return error_flag;
}

int adi_round(double val)
{
  int val_low=0,val_high=0,val_round=0;
  val_low = (int)floor(val); val_high = (int)ceil(val);
  if (fabs(val-val_high)<=fabs(val-val_low)){ val_round = val_high;} else{ val_round = val_low;}
  return val_round;
}

void doubleswap(double *val1,double *val2){ double temp=*val1; *val1=*val2; *val2=temp;}

void periodify(char *type,void *number,void *min,void *max,void *output)
{
  double maxval=0,minval=0,val=0,range=0;
  if (strcmp(type,"double")==0){ 
    maxval = *(double *)max;
    minval = *(double *)min;
    val = *(double *)number;}
  else if (strcmp(type,"int")==0){ 
    maxval = *(int *)max;
    minval = *(int *)min;
    val = *(int *)number;}
  range = maxval-minval;
  while (val<minval){ val+=range;}
  while (val>=maxval){ val-=range;}
  if (strcmp(type,"double")==0){ *(double *)output = (double)val;}
  else if (strcmp(type,"int")==0){ *(int *)output = (int)val;}
}

double llistcorrelation(struct llist *L1,struct llist *L2)
{
  /* L1,L2 are link lists of double data */
  double c=0;
  double mean1=0,mean2=0,stdev1=0,stdev2=0;
  struct litem *l1=NULL,*l2=NULL;
  if (L1->length != L2->length){ printf(" %% warning, different length lists in llistcorrelation\n");}
  l1=L1->first;l2=L2->first;
  while (l1!=NULL && l2!=NULL){
    c += (*((double *) l1->item))*(*((double *) l2->item));
    l1 = l1->child;
    l2 = l2->child;}
  c /= (double)minimum(L1->length,L2->length);
  lliststats(L1,NULL,NULL,&mean1,&stdev1);
  lliststats(L2,NULL,NULL,&mean2,&stdev2);
  c -= mean1*mean2;
  c /= stdev1*stdev2;
  return c;
}

void lliststats(struct llist *L,double *max,double *min,double *mean,double *stdev)
{
  /* L is a link list of double data. */
  double min2=0,max2=0,mean2=0,stdev2=0;
  double d=0;
  struct litem *l=NULL;
  if (L->length>0){
    l=L->first; min2=*(double *)l->item; max2=*(double *)l->item;
    while (l!=NULL){
      d = *(double *)l->item;
      max2 = maximum(max2,d); min2 = minimum(min2,d);
      mean2 += d;
      l = l->child;}
    mean2 /= L->length;
    l=L->first;
    while (l!=NULL){
      d = *(double *)l->item;
      stdev2 += pow(d-mean2,2);
      l = l->child;}
    stdev2 /= L->length;
    stdev2 = sqrt(stdev2);}
  else{ /* do nothing */ }
  if (max!=NULL){ *max=max2;}
  if (min!=NULL){ *min=min2;}
  if (mean!=NULL){ *mean=mean2;}
  if (stdev!=NULL){ *stdev=stdev2;}
}

void lliststats2(struct llist *L,double *mra,int nm,int zpad)
{
  /* L is a link list of double data, zpad is the number of zeros to be "padded" into the data */
  struct litem *l=NULL;
  double mean=0;
  int nr=0;
  if (L->length>0){ l=L->first; while (l!=NULL){ mean += *(double *)l->item; l=l->child;} mean /= (double)(L->length+zpad);}
  l=L->first; for (nr=0;nr<nm;nr++){ mra[nr]=0;}
  while (l!=NULL){ for (nr=0;nr<nm;nr++){ mra[nr] += pow(*(double *)l->item - mean,nr);} l=l->child;}
  for (nr=0;nr<nm;nr++){ mra[nr] += zpad*pow(0 - mean,nr);}
  for (nr=0;nr<nm;nr++){ mra[nr] /= (double)(L->length+zpad);}
}

void raaddra(double *A,double *B,int rows,int cols)
{ 
  /* if rows are negative, we flip them, if cols are negative, we flip them */
  int nr=0,nc=0;
  if (0){}
  else if (rows>0 && cols>0){
    for (nr=0;nr<rows;nr++){ for (nc=0;nc<cols;nc++){ 
      A[nr + nc*rows] += B[nr + nc*rows];}}}
  else if (rows>0 && cols<0){
    for (nr=0;nr<rows;nr++){ for (nc=0;nc<-cols;nc++){ 
      A[nr + nc*rows] += B[nr + (-cols-nc-1)*rows];}}}
  else if (rows<0 && cols>0){
    for (nr=0;nr<-rows;nr++){ for (nc=0;nc<cols;nc++){ 
      A[nr + nc*-rows] += B[-rows-nr-1 + nc*-rows];}}}
  else if (rows<0 && cols<0){
    for (nr=0;nr<-rows;nr++){ for (nc=0;nc<-cols;nc++){ 
      A[nr + nc*-rows] += B[-rows-nr-1 + (-cols-nc-1)*-rows];}}}
}

double radotmean(double *A,double *B,int length)
{
  int i=0;
  double c=0;
  for (i=0;i<length;i++){ c+=A[i]*B[i];} 
  c /= (double)length;
  return finite(c) ? c : 0;
}

double correlation(double *A,double *B,int length)
{
  double c=0;
  double amax,amin,amean,astd;
  double bmax,bmin,bmean,bstd;
  stats("double",A,length,&amax,&amin,&amean,&astd);
  stats("double",B,length,&bmax,&bmin,&bmean,&bstd);
  c = radotmean(A,B,length);
  c -= amean*bmean;
  c /= (astd*bstd);
  if (fabs(c)>1){ printf("warning, c=%f>1 in correlation\n",c);}
  return finite(c) ? c : 0;
}

double * spacesmear(double *ra,int rows,int cols,int s)
{
  /* assume a periodic array *ra with rows*cols dimension, we smear by s*/
  int scutoff=3;
  int i=0,j=0,i2=0,j2=0,iteration=0;
  double *A=NULL,*B=NULL;
  if (s>=minimum(rows,cols)/2){ s = (int)floor(minimum(rows,cols)/2 - 1);} s = maximum(s,0);
  if (s<=scutoff){ 
    A = (double *) tcalloc(rows*cols,sizeof(double));
    for (i=0;i<rows;i++){ for (j=0;j<cols;j++){
      iteration=0;
      for (i2=-s;i2<=s;i2++){ for (j2=-s;j2<=s;j2++){
	A[i+j*rows] += ra[periodize(i+i2,0,rows) + periodize(j+j2,0,cols)*rows];
	iteration++;}}
      A[i+j*rows] /= iteration;}}}
  else /* if (s>scutoff) */{
    A = (double *) tcalloc(rows*cols,sizeof(double));
    B = (double *) tcalloc(rows*cols,sizeof(double));
    for (i=0;i<rows;i++){ for (j=0;j<cols;j++){
      B[i+j*rows] = (periodize(i+s,0,rows)<=(2*s+1) && periodize(j+s,0,cols)<=(2*s+1));}}
    fftwconvolve(NULL,NULL,ra,B,A,rows,cols);
    for (i=0;i<rows;i++){ for (j=0;j<cols;j++){
      A[i+j*rows] /= pow(s+1,2);}}
    tfree(B);}
  return A;
}

void binstats(char *s,void *ra,int length,double *mean,double *var)
{
  /* finds the stats of a histogram */
  int i=0;
  double a=0,b=0,c=0;
  if (strcmp(s,"double")==0){
    for (i=0;i<length;i++){
      a += ((double *)ra)[i];}
    for (i=0;i<length;i++){
      b += i*((double *)ra)[i]/a;}
    for (i=0;i<length;i++){
      c += pow(i-b,2)*((double *)ra)[i]/a;}
    *mean = b;
    *var = c;}
  else if (strcmp(s,"long")==0){
    for (i=0;i<length;i++){
      a += ((long *)ra)[i];}
    for (i=0;i<length;i++){
      b += i*((long *)ra)[i]/a;}
    for (i=0;i<length;i++){
      c += pow(i-b,2)*((long *)ra)[i]/a;}
    *mean = b;
    *var = c;}
  else if (strcmp(s,"int")==0){
    for (i=0;i<length;i++){
      a += ((int *)ra)[i];}
    for (i=0;i<length;i++){
      b += i*((int *)ra)[i]/a;}
    for (i=0;i<length;i++){
      c += pow(i-b,2)*((int *)ra)[i]/a;}
    *mean = b;
    *var = c;}
}

int indexpack(int nvar,int *dlength,int *indexra)
{
  /* the opposite of indextract, given indexra recompose dindex */
  int nv=0;
  int dindex=0;
  dindex=indexra[nvar-1];
  for (nv=nvar-2;nv>=0;nv--){
    dindex *= dlength[nv];
    dindex += indexra[nv];}
  return dindex;
}

void indextract(int dindex,int nvar,int *dlength,int *indexra)
{
  /* extract individual indices from a bound multi-index */
  /* for example, (127,4,[2 5 3 7]) yields [a + b*2 + c*2*5 + d*2*5*3] yields [1 + 3*2 + 0*2*5 + 4*2*5*3 ] yields [1 3 0 4] */
  /* both int *dlength and int *indexra should be initialized */
  int check_flag=1,dindex_start=dindex;
  int curvar=0;
  if (nvar > 0){
    curvar = 0;
    do{
      indexra[curvar] = dindex%dlength[curvar];
      dindex = (dindex-indexra[curvar])/dlength[curvar];
      curvar += 1;}
    while (curvar<nvar);}
  if (check_flag){ dindex = indexpack(nvar,dlength,indexra); assert(dindex==dindex_start);}  
}

void maxmindex(char *s,void *ra,int length,int **maxindex,int *nmax,int **minindex,int *nmin,double tolerance)
{
  /* finds the indices corresponding to maximum and minimum elements, assuming *maxindex and *minindex are both NULL */
  int verbose=0;
  double max=0,min=0;
  int nr=0;
  struct llist *Lmax=NULL,*Lmin=NULL;
  double *temp=NULL;
  struct litem *l=NULL;
  if (verbose){ printf(" %% [entering maxmindex] with ra: \n"); raprintf(ra,s,1,length," ");}
  stats(s,ra,length,&max,&min,NULL,NULL);
  if (verbose){ printf(" %% found max %f, min %f\n",max,min);}
  Lmax=llistmake();Lmin=llistmake();
  if (strcmp(s,"double")==0){
    for (nr=0;nr<length;nr++){
      if (fabs(((double *)ra)[nr]-max)<tolerance){ temp=(double *)tmalloc(sizeof(double));*temp=(double)nr;litemadd(Lmax,temp);}
      if (fabs(((double *)ra)[nr]-min)<tolerance){ temp=(double *)tmalloc(sizeof(double));*temp=(double)nr;litemadd(Lmin,temp);}}}
  else if (strcmp(s,"int")==0){
    for (nr=0;nr<length;nr++){
      if (fabs(((int *)ra)[nr]-max)<tolerance){ temp=(double *)tmalloc(sizeof(double));*temp=(double)nr;litemadd(Lmax,temp);}
      if (fabs(((int *)ra)[nr]-min)<tolerance){ temp=(double *)tmalloc(sizeof(double));*temp=(double)nr;litemadd(Lmin,temp);}}}
  else if (strcmp(s,"long")==0){
    for (nr=0;nr<length;nr++){
      if (fabs(((long *)ra)[nr]-max)<tolerance){ temp=(double *)tmalloc(sizeof(double));*temp=(double)nr;litemadd(Lmax,temp);}
      if (fabs(((long *)ra)[nr]-min)<tolerance){ temp=(double *)tmalloc(sizeof(double));*temp=(double)nr;litemadd(Lmin,temp);}}}
  else{ printf(" %% Warning, improper type %s in maxmindex\n",s);}
  *nmax=Lmax->length;*nmin=Lmin->length;
  *maxindex = (int *) tcalloc(*nmax,sizeof(int));
  *minindex = (int *) tcalloc(*nmin,sizeof(int));
  nr=0;l=Lmax->first;
  while (l!=NULL){
    (*maxindex)[nr] = (int)*((double *)l->item);
    nr += 1;
    l=l->child;}
  nr=0;l=Lmin->first;
  while (l!=NULL){
    (*minindex)[nr] = (int)*((double *)l->item);
    nr += 1;
    l=l->child;}
  llisttfree2(Lmax);Lmax=NULL;
  llisttfree2(Lmin);Lmin=NULL;
  if (verbose){ printf(" %% found maxindex of length %d\n",*nmax); raprintf(*maxindex,"int",1,*nmax,"max");}
  if (verbose){ printf(" %% found minindex of length %d\n",*nmin); raprintf(*minindex,"int",1,*nmin,"min");}
}

void stats(char *s,void *ra,int length,double *max,double *min,double *mean,double *stdev)
{
  /* finds the stats of an array */
  int i=0;
  double max2=0,min2=0,mean2=0,stdev2=0;
  if (ra!=NULL){
    if (strcmp(s,"double")==0){
      for (i=0;i<length;i++){
	mean2 += ((double *)ra)[i];}
      mean2 /= (double)length;
      max2=mean2;min2=mean2;
      max2=((double *)ra)[0]; min2=((double *)ra)[0];
      for (i=0;i<length;i++){
	if (((double *)ra)[i]>max2){ max2=((double *)ra)[i];}
	if (((double *)ra)[i]<min2){ min2=((double *)ra)[i];}
	stdev2 += (((double *)ra)[i]-mean2)*(((double *)ra)[i]-mean2);}
      stdev2 /= (double)length;
      stdev2 = sqrt(stdev2);}
    else if (strcmp(s,"long")==0){
      for (i=0;i<length;i++){
	mean2 += ((long *)ra)[i];}
      mean2 /= (double)length;
      max2=mean2;min2=mean2;
      max2=((long *)ra)[0]; min2=((long *)ra)[0];
      for (i=0;i<length;i++){
	if (((long *)ra)[i]>max2){ max2=((long *)ra)[i];}
	if (((long *)ra)[i]<min2){ min2=((long *)ra)[i];}
	stdev2 += (((long *)ra)[i]-mean2)*(((long *)ra)[i]-mean2);}
      stdev2 /= (double)length;
      stdev2 = sqrt(stdev2);}
    else if (strcmp(s,"int")==0){
      for (i=0;i<length;i++){
	mean2 += ((int *)ra)[i];}
      mean2 /= (double)length;
      max2=mean2;min2=mean2;
      max2=((int *)ra)[0]; min2=((int *)ra)[0];
      for (i=0;i<length;i++){
	if (((int *)ra)[i]>max2){ max2=((int *)ra)[i];}
	if (((int *)ra)[i]<min2){ min2=((int *)ra)[i];}
	stdev2 += (((int *)ra)[i]-mean2)*(((int *)ra)[i]-mean2);}
      stdev2 /= (double)length;
      stdev2 = sqrt(stdev2);}
    if (max!=NULL){*max=max2;}
    if (min!=NULL){*min=min2;}
    if (mean!=NULL){*mean=mean2;}
    if (stdev!=NULL){*stdev=stdev2;}}
}

void rgb2hsv(double r,double g,double b,double *h,double *s,double *v)
{
  /* shamelessly stolen from http://cs.haifa.ac.il/hagit/courses/ist/Lectures/Demos/ColorApplet2/t_convert.html
     r,g,b values are from 0 to 1
     h = [0,360], s = [0,1], v = [0,1]
     if s == 0, then h = -1 (undefined) */
  double min=0,max=0,delta=0;
  min = minimum(r,minimum(g,b));
  max = maximum(r,maximum(g,b));
  *v = max;
  delta = max - min;
  if (max!=0){ *s = delta/max;}
  else{ *s = 0; *h = -1; return;}
  if (r==max){ /* between yellow and magenta */ *h = (g - b)/delta;}
  else if (g==max){ /* between cyan and yellow */ *h = 2+(b-r)/delta;}
  else /* if (b==max) */{ /* between magenta and cyan */ *h = 4+(r-g)/delta;}
  *h *= 60; /* degrees */
  if (*h<0){ *h+=360;}
}

void hsv2rgb(double h,double s,double v,double *r,double *g,double *b)
{
  int i;
  double f=0,p=0,q=0,t=0;
  if (s==0){ /* greyscale */ *r=v;*g=v;*b=v; return;}
  h /= 60; i=floor(h); /* six sectors */
  f = h-i;
  p = v*(1-s);
  q = v*(1-s*f);
  t = v*(1-s*(1-f));
  switch (i){
  case 0: *r=v;*g=t;*b=p; break;
  case 1: *r=q;*g=v;*b=p; break;
  case 2: *r=p;*g=v;*b=t; break;
  case 3: *r=p;*g=q;*b=v; break;
  case 4: *r=t;*g=p;*b=v; break;
  default: /* case 5 */ *r=v;*g=p;*b=q; break;}
}

void colorscale(int coloringtype,double val,double valmax,double valmin,double *rcolor,double *gcolor,double *bcolor)
{
  double gamma=(coloringtype==5?1.0:0.5);
  int norm_flag=1;
  double norm=1;
  double v=0;
  double h=0;
  switch (coloringtype){
  case 0: /* greyscale coloring */
    if (valmax<=valmin){ valmax=valmin+1;}
    v = (val-valmin)/(valmax-valmin); v = maximum(0,minimum(1,v));
    *rcolor = v; *gcolor = v; *bcolor = v;
    break;
  case 2: /* cyclical coloring */
    if (valmax<=valmin){ valmax=valmin+1;}
    v = 3*((val-valmin)/(valmax-valmin)-0.5);
    if (v < -0.5 || v > 0.5){ *rcolor = sqrt(fabs(v)-0.5);} else{ *rcolor=0;}
    if (v > -0.5){ *gcolor = sqrt(1-fabs(v-0.5));} else{ *gcolor=0;}
    if (v < 0.5){ *bcolor = sqrt(1-fabs(v+0.5));} else{ *bcolor=0;}
    break;
  case 3: /* IC IS EC ES ALL coloring */
    switch ((int)val){
    case 0: /* IC */ *rcolor=0.25; *gcolor=0.75; *bcolor=1.00; break;
    case 1: /* IS */ *rcolor=0.00; *gcolor=0.00; *bcolor=1.00; break;
    case 2: /* EC */ *rcolor=1.00; *gcolor=0.75; *bcolor=0.25; break;
    case 3: /* ES */ *rcolor=1.00; *gcolor=0.00; *bcolor=0.00; break;
    case 4: /* ALL */ *rcolor=0.75; *gcolor=0.00; *bcolor=0.75; break;
    default: *rcolor=1; *gcolor=1; *bcolor=1; break;}
    break;
  case 4: /* coloring with gamma */
    if (valmax==valmin){ valmax++;}
    v = (val-valmin)/(valmax-valmin);
    if (v<0){ v=0;}
    if (v>1){ v=1;}
    v = pow(v,gamma);
    if (!GRAYSCALE){
      if (v<0.5){
	v = asin(2*v);
	*rcolor = 0;
	*gcolor = sin(v);
	*bcolor = cos(v);}
      else if (v>=0.5){
	v = asin(2*(v-0.5));
	*rcolor = sin(v);
	*gcolor = cos(v);
	*bcolor = 0;}}
    else if (GRAYSCALE){
      *rcolor = v;
      *gcolor = v;
      *bcolor = v;}
    break;
  case 5: /* super coloring with blue-cyan-green-yellow-red and gamma value */
    if (valmax==valmin){ valmax++;}
    v = (val-valmin)/(valmax-valmin);
    if (v<0){ v=0;}
    if (v>1){ v=1;}
    v = pow(v,gamma);
    if (!GRAYSCALE){
      if (v>=0.0 && v<0.25){
	*rcolor = 0;
	*gcolor = (v-0.0)/0.25;
	*bcolor = 1;}
      if (v>=0.25 && v<0.5){
	*rcolor = 0;
	*gcolor = 1;
	*bcolor = 1-(v-0.25)/0.25;}
      if (v>=0.5 && v<0.75){
	*rcolor = (v-0.5)/0.25;
	*gcolor = 1;
	*bcolor = 0;}
      if (v>=0.75 && v<=1){
	*rcolor = 1;
	*gcolor = 1-(v-0.75)/0.25;
	*bcolor = 0;}
      if (norm_flag){ norm = sqrt(pow(*rcolor,2)+pow(*gcolor,2)+pow(*bcolor,2)); *rcolor /= norm; *gcolor /= norm; *bcolor /= norm;}}
    else if (GRAYSCALE){
      *rcolor = v;
      *gcolor = v;
      *bcolor = v;}
    break;
  case 6: /* use hsv2rgb with h in [0,360] */
    if (valmax<=valmin){ valmax=valmin+1;}
    v = crop(val,valmin,valmax);
    h = 359.99*(v-valmin)/(valmax-valmin);
    hsv2rgb(h,1,1,rcolor,gcolor,bcolor);
    break;
  case 7: /* use hsv2rgb with h in [240,-60] */
    if (valmax<=valmin){ valmax=valmin+1;}
    v = crop(val,valmin,valmax);
    h = 300*(valmax-v)/(valmax-valmin)-60;
    hsv2rgb(h,1,1,rcolor,gcolor,bcolor);
    break;
  default:
    if (valmax==valmin){ valmax++;}
    v = (val-valmin)/(valmax-valmin);
    if (v<0){ v=0;}
    if (v>1){ v=1;}
    if (!GRAYSCALE){
      if (v<0.5){
	v = asin(2*v);
	*rcolor = 0;
	*gcolor = sin(v);
	*bcolor = cos(v);}
      else if (v>=0.5){
	v = asin(2*(v-0.5));
	*rcolor = sin(v);
	*gcolor = cos(v);
	*bcolor = 0;}}
    else if (GRAYSCALE){
      *rcolor = v;
      *gcolor = v;
      *bcolor = v;}
    break;}
}

void colorscaleinv(int coloringtype,double rcolor,double gcolor,double bcolor,double valmax,double valmin,double *val)
{
  /* reverses the color function */
  double mincolor=minimum(rcolor,minimum(gcolor,bcolor));
  double scaleval=0;
  double v=0;
  switch (coloringtype){
  case 0: /* reverses the greyscale coloring */
    v = (rcolor+gcolor+bcolor)/3;
    *val = valmin + v*(valmax-valmin);
    break;
  case 2: /* reverses standard cyclical coloring */
    rcolor = rcolor-mincolor;
    gcolor = gcolor-mincolor;
    bcolor = bcolor-mincolor;
    if (rcolor==0){ scaleval = 0.5 - pow(bcolor,2);}
    if (bcolor==0){ scaleval = 1.5 - pow(gcolor,2);}
    if (gcolor==0){ scaleval = -(pow(rcolor,2)+1/2);}
    if (val!=NULL){ *val = valmin + (scaleval+1.5)/3.0*(valmax-valmin);}
    break;
  default:
    rcolor = rcolor-mincolor;
    gcolor = gcolor-mincolor;
    bcolor = bcolor-mincolor;
    if (rcolor==0){ scaleval = gcolor/2;}
    else{ scaleval = rcolor/2+0.5;}
    if (val!=NULL){ *val = valmin + scaleval*(valmax-valmin);}
    break;}
}

void fftwconvolve(void * pforward,void * pbackwards,double *ra1,double *ra2,double *ra3,int rows,int cols)
{
  /* uses /usr/local/lib/libfftw3.a to quickly convolve two real arrays of the same size
     ra1 and ra2 are stored in column major order (row_number + col_number*number_of_rows)
     fftw takes arrays in row major order (col_number + rows_number*number_of_cols)
     the output is put into the (precreated) ra3 
     if pforward,pbackward are both NULL, plans are created and destroyed.
     otherwise pforward and pbackwards are used 
     Specifically, this implements the following operation:
     for (nr=0;nr<nrows;nr++){ for (nc=0;nc<ncols;nc++){
       ra3[nr+nc*nrows] = 0;
       for (nr2=0;nr2<nrows;nr2++){ for (nc2=0;nc2<ncols;nc2++){
         nr3 = periodize(nr2-nr,0,nrows); nc3 = periodize(nc2-nc,0,ncols);
         ra3[nr+nc*nrows] += ra1[nr2+nc2*nrows]*ra2[nr3+nc3*nrows];}}}}
     Note that the entries of ra2 must be 'flipped'.
  */
  int nr=0,nc=0,tab=0;
  int create_plans=0;
  fftw_complex *cra1=NULL,*cra2=NULL,*cra1z=NULL,*cra2z=NULL;
  fftw_plan temp_plan_forward,temp_plan_backward,*plan_forward=NULL,*plan_backward=NULL;
  if (pforward==NULL || pbackwards==NULL){ create_plans=1;} else{ create_plans=0;}
  cra1 = (fftw_complex *) fftw_malloc(rows*cols*sizeof(fftw_complex));
  cra2 = (fftw_complex *) fftw_malloc(rows*cols*sizeof(fftw_complex));
  cra1z = (fftw_complex *) fftw_malloc(rows*cols*sizeof(fftw_complex));
  cra2z = (fftw_complex *) fftw_malloc(rows*cols*sizeof(fftw_complex));
  if (create_plans){
    temp_plan_forward = fftw_plan_dft_2d(rows,cols,cra1,cra1z,-1,FFTW_ESTIMATE);
    temp_plan_backward = fftw_plan_dft_2d(rows,cols,cra1,cra2,+1,FFTW_ESTIMATE);
    plan_forward = &temp_plan_forward;
    plan_backward = &temp_plan_backward;}
  else{
    plan_forward = (fftw_plan *) pforward;
    plan_backward = (fftw_plan *) pbackwards;}
  for (nr=0;nr<rows;nr++){ for (nc=0;nc<cols;nc++){
    cra1[nc + nr*cols][0] = ra1[nr + nc*rows]; cra1[nc+nr*cols][1] = 0;
    cra2[nc + nr*cols][0] = ra2[periodize(rows-nr,0,rows) + periodize(cols-nc,0,cols)*rows]; cra2[nc+nr*cols][1] = 0;}}
  fftw_execute_dft(*plan_forward,cra1,cra1z);
  fftw_execute_dft(*plan_forward,cra2,cra2z);
  for (nr=0;nr<rows;nr++){ for (nc=0;nc<cols;nc++){
    tab = nc + nr*cols;
    cra1[tab][0] = cra1z[tab][0]*cra2z[tab][0] - cra1z[tab][1]*cra2z[tab][1];
    cra1[tab][1] = cra1z[tab][0]*cra2z[tab][1] + cra1z[tab][1]*cra2z[tab][0];}}
  fftw_execute_dft(*plan_backward,cra1,cra2);
  for (nr=0;nr<rows;nr++){ for (nc=0;nc<cols;nc++){
    ra3[nr + nc*rows] = cra2[nc + nr*cols][0]/pow(rows*cols,1);}}
  if (create_plans){
    fftw_destroy_plan(temp_plan_forward); fftw_destroy_plan(temp_plan_backward);}
  free(cra1);free(cra2);free(cra1z);free(cra2z);
}

double * ReadPNMfile(char *filename,int readra_flag,int readcomment_flag,int *colorsout,double *maxout,double *minout,int *widthout,int *heightout)
{
  /* This reads a raw-bits PNM file of width w and height h into the unsigned char * array output
     all header lines end with '\n' (10)
     comment lines start with '#' (35)
     exactly one comment line is expected
     otherwise this will not work,
     assumes filaname starts with "./" */
  int verbose=0;
  int i=0,j=0,k=0;
  unsigned char *output=NULL;
  double max=0,min=0;
  double *output2=NULL;
  unsigned char *buffer=NULL;
  FILE *fp=NULL;
  int width=0,height=0,maxval=0,colors=0;
  buffer = (unsigned char *) tmalloc(sizeof(unsigned char)); 
  if ( (fp = fopen(filename,"r"))==NULL ){ printf(" %% warning, could not open %s in ReadPNMfile\n",filename); return NULL;}
  /* first line */
  fscanf(fp,"%c",buffer); fscanf(fp,"%d",&colors);
  i=0; while (i<1){ fread(buffer,sizeof(unsigned char),1,fp); i += (int)buffer[0]==10;}
  if (verbose){ printf("reading pnm file of colors %d\n",colors);}
  colors = colors <= 5 ? 1 : 3;
  /* comment line */
  if (readcomment_flag){ fscanf(fp,"%c",buffer); fscanf(fp,"%lf",&max); fscanf(fp,"%lf",&min);}
  i=0; while (i<1){ fread(buffer,sizeof(unsigned char),1,fp); i += (int)buffer[0]==10;}
  if (verbose){ printf("reading max=%f,min=%f\n",max,min);}
  /* dimension line */
  fscanf(fp,"%d",&width); fscanf(fp,"%d",&height); i=0; while (i<1){ fread(buffer,sizeof(unsigned char),1,fp); i += (int)buffer[0]==10;}
  /* maxval line */
  fscanf(fp,"%d",&maxval); i=0; while (i<1){ fread(buffer,sizeof(unsigned char),1,fp); i += (int)buffer[0]==10;}
  if (verbose){ printf("reading width %d and height %d and maxval %d\n",width,height,maxval);}
  tfree(buffer);
  if (readra_flag){
    output = (unsigned char *) tmalloc(sizeof(unsigned char)*width*height*colors);
    fread(output,sizeof(unsigned char),width*height*colors,fp);}
  fclose(fp);
  if (readra_flag){
    if (verbose){ for (i=0;i<height;i++){ for (j=0;j<width;j++){ for (k=0;k<colors;k++){ printf("%d ",output[k+j*colors+i*width*colors]);} printf("\n");}} printf("\n");}
    output2 = (double *) tmalloc(sizeof(double)*width*height*colors);
    for (i=0;i<height;i++){ for (j=0;j<width;j++){ for (k=0;k<colors;k++){
      if (verbose){ printf("at row %d, col %d, reading %d\n",i,j,output[k+j*colors+i*width*colors]);}
      output2[k+i*colors+j*height*colors] = (double)(output[k+j*colors+i*width*colors])/(double)maxval;}}}
    if (verbose){ for (i=0;i<height;i++){ for (j=0;j<width;j++){ for (k=0;k<colors;k++){ printf("%d ",(int)output2[k+i*colors+j*height*colors]);} printf("\n");}} printf("\n");}}
  if (colorsout!=NULL){ *colorsout=colors;}
  if (maxout!=NULL){ *maxout=max;}
  if (minout!=NULL){ *minout=min;}
  if (widthout!=NULL){ *widthout=width;}
  if (heightout!=NULL){ *heightout=height;}
  if (output!=NULL){ tfree(output);}
  return output2;
}

int WritePNMfile(double *ra,int rows,int cols,double max,double min,char *filename)
{
  /* This writes a double* array of height h and width w (stored in row-major order) into a pnm file at char* filename 
     assumes filename starts with "./" */
  int nr=0,nc=0;
  double max2=0,min2=0,mean=0,stdev=0;
  double d=0;
  unsigned char ud=0;
  FILE *fp;
  if (max<=min){ stats("double",ra,rows*cols,&max2,&min2,&mean,&stdev); max2=mean+STD_VIEW*stdev; min2=mean-STD_VIEW*stdev;}
  else{ max2=max; min2=min;}
  if ((fp=fopen(filename,"w"))==NULL) { printf(" %% Warning: cannot open %s in WritePNMFile, writing to stdout\n",filename); fp=stdout;}
  //fprintf(fp,"P5\n# CREATOR: WritePNMfile()\n%d %d\n%d\n",cols,rows,UCHAR_MAX);
  fprintf(fp,"P5\n# %0.16lf %0.16lf\n%d %d\n%d\n",max2,min2,cols,rows,UCHAR_MAX);
  for (nr=0;nr<rows;nr++){ for (nc=0;nc<cols;nc++){
    d = ra[nr+nc*rows];
    ud = (unsigned char)(UCHAR_MAX*(d-min2)/(max2-min2));
    fwrite(&ud,sizeof(unsigned char),1,fp);}}
  if (fp!=stdout){ 
    fclose(fp);
    printf(" %% wrote file %s, (%d x %d pixels, %d bytes)\n",filename,rows,cols,rows*cols*sizeof(unsigned char));
    return 1;}
  else{ printf(" %% wrote file %s to stdout\n",filename); return 0;}
  return 0;
}

int WritePNMfile_color(double *ra,int rows,int cols,double max,double min,char *filename,int coloringtype)
{
  /* This writes a double* array of height h and width w (stored in row-major order) into a color pnm file at char* filename 
     assumes filename starts with "./" */
  int nr=0,nc=0;
  double max2=0,min2=0,mean=0,stdev=0;
  double d=0,rcolor=0,gcolor=0,bcolor=0;
  unsigned char *ud=NULL;
  FILE *fp;
  if (max<=min){ stats("double",ra,rows*cols,&max2,&min2,&mean,&stdev); max2=mean+STD_VIEW*stdev; min2=mean-STD_VIEW*stdev;} 
  else{ max2=max; min2=min;}
  if ((fp=fopen(filename,"w"))==NULL) { printf(" %% Warning: cannot open %s in WritePNMFile_color, writing to stdout\n",filename); fp=stdout;}
  //fprintf(fp,"P6\n# CREATOR: WritePNMfile_color()\n%d %d\n%d\n",cols,rows,UCHAR_MAX);
  fprintf(fp,"P6\n# %0.16lf %0.16lf\n%d %d\n%d\n",max2,min2,cols,rows,UCHAR_MAX);
  ud = (unsigned char *) tcalloc(3*rows*cols,sizeof(unsigned char));
  for (nr=0;nr<rows;nr++){ for (nc=0;nc<cols;nc++){
    d = ra[nr+nc*rows];
    switch(coloringtype){
    case 3: if (d>mean+(1-STD_VIEW)*stdev){ colorscale(1,d,max2,min2,&rcolor,&gcolor,&bcolor);} else{ rcolor=0,gcolor=0,bcolor=0;} break;
    default: colorscale(coloringtype,d,max2,min2,&rcolor,&gcolor,&bcolor); break;}
    ud[0+nc*3+nr*3*cols] = (unsigned char)(UCHAR_MAX*rcolor);
    ud[1+nc*3+nr*3*cols] = (unsigned char)(UCHAR_MAX*gcolor);
    ud[2+nc*3+nr*3*cols] = (unsigned char)(UCHAR_MAX*bcolor);}}
  fwrite(ud,sizeof(unsigned char),3*rows*cols,fp);
  tfree(ud);
  if (fp!=stdout){ 
    fclose(fp);
    printf(" %% wrote file %s, (%d x %d pixels, %d bytes)\n",filename,rows,cols,rows*cols*sizeof(unsigned char));
    return 1;}
  else{ printf(" %% wrote file %s to stdout\n",filename); return 0;}
  return 0;
}

void num2frame(int num,char *frame)
{
  if (num<10){ sprintf(frame,"000%d",num);}
  else if (num<100){ sprintf(frame,"00%d",num);}
  else if (num<1000){ sprintf(frame,"0%d",num);}
  else if (num<10000){ sprintf(frame,"%d",num);}
}

int RescalePNMfiles(char *filename,int framestart,int frameend)
{
  /* assuming all frames are of the form filename_framenumber.pnm, where filename starts with "./"
     assuming all frames are written by WritePNMfile() or WritePNMfile_color() */
  int verbose=0;
  char text[128],framename[16],textout[128];
  int colors=0,width=0,height=0;
  double *maxvals=NULL,*minvals=NULL,maxmean=0,minmean=0,maxnow=0,minnow=0;
  int nf=0,i=0,nframes=frameend-framestart+1;
  double *A=NULL,*B=NULL;
  double rcolor=0,gcolor=0,bcolor=0;
  int output_flag=1;
  maxvals = (double *) tcalloc(nframes,sizeof(double));
  minvals = (double *) tcalloc(nframes,sizeof(double));
  for (nf=0;nf<nframes;nf++){
    num2frame(nf+framestart,framename); sprintf(text,"%s.%s.pnm",filename,framename);
    ReadPNMfile(text,0,1,&colors,maxvals+nf,minvals+nf,&width,&height);}
  stats("double",maxvals,nframes,NULL,NULL,&maxmean,NULL);
  stats("double",minvals,nframes,NULL,NULL,&minmean,NULL);
  tfree(maxvals);tfree(minvals);
  if (verbose){ printf("%f,%f\n",maxmean,minmean);}
  for (nf=0;nf<nframes;nf++){
    num2frame(nf+framestart,framename); sprintf(text,"%s.%s.pnm",filename,framename); sprintf(textout,"%s_rs.%s.pnm",filename,framename);
    A = ReadPNMfile(text,1,1,&colors,&maxnow,&minnow,&width,&height);
    B = (double *) tcalloc(width*height,sizeof(double));
    for (i=0;i<width*height;i++){
      if (colors==1){ B[i] = minnow + A[i]*(maxnow-minnow);}
      else{ rcolor=A[0+i*colors]; gcolor=A[1+i*colors]; bcolor=A[2+i*colors]; colorscaleinv(1,rcolor,gcolor,bcolor,maxnow,minnow,B+i);}}
    if (colors==1){ output_flag *= WritePNMfile(B,height,width,maxmean,minmean,textout);}
    else{ output_flag *= WritePNMfile_color(B,height,width,maxmean,minmean,textout,1);}
    tfree(B);
    tfree(A);}
  return output_flag;
}

void rareset(void *vra,char *type,int length,void *filler)
{
  int nr=0;
  double *dra=NULL,dval=0;
  int *ira=NULL,ival=0;
  if (strcmp(type,"double")==0){
    dra = (double *) vra;
    dval = (filler!=NULL ? *(double *)filler : 0);
    for (nr=0;nr<length;nr++){ dra[nr]=dval;}}
  else if (strcmp(type,"int")==0){
    ira = (int *) vra;
    ival = (filler!=NULL ? *(int *)filler : 0);
    for (nr=0;nr<length;nr++){ ira[nr]=ival;}}
  else{ printf(" warning, poor type %s in rareset\n",type);}
}

void raprintf(void *vra,char *type,int rows,int cols,char *prefix)
{
  int nr=0,nc=0;
  double *dra=NULL;
  int *ira=NULL;
  double printftol=0.000000001;
  if (strcmp(type,"double")==0){
    dra = (double *) vra;
    for (nr=0;nr<rows;nr++){ 
      printf("%s",prefix); 
      for (nc=0;nc<cols;nc++){ 
	if (fabs(dra[nr+nc*rows]-(int)dra[nr+nc*rows])<printftol){ printf(" %d",(int)dra[nr+nc*rows]);}
	else{ printf(" %f",dra[nr+nc*rows]);}} 
      printf("\n");}}
  else if (strcmp(type,"int")==0){
    ira = (int *) vra;
    for (nr=0;nr<rows;nr++){ printf("%s",prefix); for (nc=0;nc<cols;nc++){ printf(" %d",ira[nr+nc*rows]);} printf("\n");}}
  else{ printf(" warning, poor type %s in raprintf\n",type);}
}

void ra2jpg(void *vra,char *type,int rows,int cols,int transpose_flag,char *filename_base,int logflag)
{
  /* plots rows of *vra against column index using a single axis 
     or if transpose_flag, plots cols of *vra against row index 
     filename_base should start with "./" 
  logflag=0 --> normal plot
  logflag=1 --> semilogy
  logflag=2 --> semilogx
  logflag=3 --> loglog */
  int remove_flag=0,jpeg_flag=(ON_MY_COMPUTER ? 1 : 0);
  double maxdia = 10000;
  FILE *fp=NULL;
  char filename[512];
  int rows2=0,cols2=0;
  double max=0,min=0,mean=0,stdev=0;
  double *dra=NULL;
  int *ira=NULL;
  int dori=0;
  int nr=0,nc=0,tab=0;
  int lcolor=0,lcolormin=0,lcolormax=30;
  char command[1024];
  double xpos=0,ypos=0;
  rows2 = transpose_flag ? cols : rows;
  cols2 = transpose_flag ? rows : cols;
  sprintf(filename,"%s.fig",filename_base);
  if ((fp=fopen(filename,"w"))==NULL){ printf(" %% Warning, cannot create %s in ra2jpg\n",filename); fp=stdout;}
  fprintf(fp,"%s",FIG_PREAMBLE);
  stats(type,vra,rows*cols,&max,&min,&mean,&stdev);
  if (logflag==1 || logflag==3){ min = 0.00000000001*max;}
  else{ max = mean + STD_VIEW*stdev; min = mean - STD_VIEW*stdev;}
  if (strcmp(type,"int")==0){ ira = (int *) vra; dori=0;}
  else if (strcmp(type,"double")==0){ dra = (double *) vra; dori=1;}
  else { printf(" %% Warning, improper type %s in ra2jpg\n",type);}
  for (nr=0;nr<rows2;nr++){ 
    periodify("int",&nr,&lcolormin,&lcolormax,&lcolor);
    fprintf(fp,"2 1 0 5 %d 7 %d 0 -1 0.000 0 0 -1 0 0 %d\n\t",/*color*/lcolor,/*depth*/(nr%999)+1,/*npoints*/cols2);
    for (nc=0;nc<cols2;nc++){
      tab = transpose_flag ? nc + nr*rows : nr + nc*rows;
      xpos = (logflag==2 || logflag==3) ? log((double)nc + 0.5)/log((double)cols2) : ((double)nc + 0.5)/(double)cols2;
      ypos = (logflag==1 || logflag==3) ? (log((double)fabs(maximum(dori ? dra[tab] : ira[tab],min)))-log(min))/(log(max)-log(min)) : ((double)(dori ? dra[tab] : ira[tab])-min)/(max-min);
      fprintf(fp," %d %d",(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos));}
    fprintf(fp,"\n");}
  fprintf(fp,"4 0 0 %d 0 %d %d 0.0000 4 %d %d %d %d %0.16f\\001\n",/*depth*/1,/*font*/12,/*point*/24,/*textheight*/(int)(24*13.5),/*textwidth*/(int)(24*9*18),/*xpos*/(int)(0.1*maxdia),/*ypos*/(int)maxdia-(int)(1.0*maxdia),(logflag==1 || logflag==3) ? log(max) : max);
  fprintf(fp,"4 0 0 %d 0 %d %d 0.0000 4 %d %d %d %d %0.16f\\001\n",/*depth*/1,/*font*/12,/*point*/24,/*textheight*/(int)(24*13.5),/*textwidth*/(int)(24*9*18),/*xpos*/(int)(0.1*maxdia),/*ypos*/(int)maxdia-(int)(0.0*maxdia),(logflag==1 || logflag==3) ? log(min) : min);
  if (fp!=stdout){ fclose(fp);}
  if (jpeg_flag){ sprintf(command,"fig2dev -Ljpeg -q %d %s.fig %s.jpg;",/*quality*/5,filename_base,filename_base); system(command);}
  if (remove_flag){ sprintf(command,"rm %s.fig;",filename_base); system(command);}
}

void ra2jpg2(void *vra,char *type,int rows,int cols,int transpose_flag,double max,double min,char *filename_base)
{
  /* plots rows of *vra against column index using a single axis 
     or if transpose_flag, plots cols of *vra against row index 
     filename_base should start with "./" 
     vertical axis between max and min */
  int remove_flag=0,jpeg_flag=(ON_MY_COMPUTER ? 1 : 0);
  double maxdia = 10000;
  FILE *fp=NULL;
  char filename[512];
  int rows2=0,cols2=0;
  double mean=0,stdev=0;
  double *dra=NULL;
  int *ira=NULL;
  int dori=0;
  int nr=0,nc=0,tab=0;
  int lcolor=0,lcolormin=0,lcolormax=30;
  char command[1024];
  double xpos=0,ypos=0;
  rows2 = transpose_flag ? cols : rows;
  cols2 = transpose_flag ? rows : cols;
  if (max<=min){ stats(type,vra,rows*cols,NULL,NULL,&mean,&stdev); max=mean+STD_VIEW*stdev; min=mean-STD_VIEW*stdev;}
  sprintf(filename,"%s.fig",filename_base);
  if (checktofind(filename)){ 
    if ((fp=fopen(filename,"a"))==NULL){ printf(" %% Warning, can't open %s in ra2jpg2\n",filename); fp=stdout;}}
  else{ 
    if ((fp=fopen(filename,"w"))==NULL){ printf(" %% Warning, can't create %s in ra2jpg2\n",filename); fp=stdout;} 
    fprintf(fp,"%s",FIG_PREAMBLE);
    fprintf(fp,"4 0 0 %d 0 %d %d 0.0000 4 %d %d %d %d %0.16f\\001\n",/*depth*/1,/*font*/12,/*point*/24,/*textheight*/(int)(24*13.5),/*textwidth*/(int)(24*9*18),/*xpos*/(int)(0.1*maxdia),/*ypos*/(int)maxdia-(int)(1.0*maxdia),max);
    fprintf(fp,"4 0 0 %d 0 %d %d 0.0000 4 %d %d %d %d %0.16f\\001\n",/*depth*/1,/*font*/12,/*point*/24,/*textheight*/(int)(24*13.5),/*textwidth*/(int)(24*9*18),/*xpos*/(int)(0.1*maxdia),/*ypos*/(int)maxdia-(int)(0.0*maxdia),min);}
  if (strcmp(type,"int")==0){ ira = (int *) vra; dori=0;}
  else if (strcmp(type,"double")==0){ dra = (double *) vra; dori=1;}
  else { printf(" %% Warning, improper type %s in ra2jpg\n",type);}
  for (nr=0;nr<rows2;nr++){ 
    periodify("int",&nr,&lcolormin,&lcolormax,&lcolor);
    fprintf(fp,"2 1 0 5 %d 7 %d 0 -1 0.000 0 0 -1 0 0 %d\n\t",/*color*/lcolor,/*depth*/(nr%999)+1,/*npoints*/cols2);
    for (nc=0;nc<cols2;nc++){
      tab = transpose_flag ? nc + nr*rows : nr + nc*rows;
      xpos = ((double)nc + 0.5)/(double)cols2;
      ypos = ((double)(dori ? dra[tab] : ira[tab])-min)/(max-min);
      fprintf(fp," %d %d",(int)floor(maxdia*xpos),(int)maxdia-(int)floor(maxdia*ypos));}
    fprintf(fp,"\n");}
  if (fp!=stdout){ fclose(fp);}
  if (jpeg_flag){ sprintf(command,"fig2dev -Ljpeg -q %d %s.fig %s.jpg;",/*quality*/5,filename_base,filename_base); system(command);}
  if (remove_flag){ sprintf(command,"rm %s.fig;",filename_base); system(command);}
}

int checktofind_howmany(char *fpre,int *index,char *fpos)
{
  char filename[1024];
  int continue_flag=0;
  int nfound=0;
  continue_flag=1;nfound=0;
  do{
    sprintf(filename,"%s%d%s",fpre,*index,fpos); 
    continue_flag=checktofind(filename);
    if (continue_flag){ *index += 1; nfound += 1;}}
  while (continue_flag);
  return nfound;
}

int checktofind(char *fgvn)
{
  /* assumes fgvn starts with "./" */
  int found_flag=0;
  FILE *fp;
  if ((fp=fopen(fgvn,"r"))!=NULL){ found_flag=1; fclose(fp);fp=NULL;}
  return found_flag;
}

int multidradump(double *dra,int ndim,int *ira,char *filename)
{
  /* assumes filename starts with "./" */
  int nd=0;
  int length=0;
  FILE *fp=NULL;
  int output_flag=0;
  length=1; for (nd=0;nd<ndim;nd++){ length*=ira[nd];}
  if ((fp=fopen(filename,"w"))==NULL){ printf(" %% warning, couldn't open %s in multidradump\n",filename);}
  else /* if fp!=NULL */{
    output_flag += fwrite(&ndim,sizeof(int),1,fp);
    output_flag += fwrite(ira,sizeof(int),ndim,fp);
    output_flag += fwrite(dra,sizeof(double),length,fp);
    if (fp!=NULL){ fclose(fp);fp=NULL;}}
  return output_flag;
}

int multidraread(char *filename,int **ira,double **dra)
{
  /* assumes *filename dumped with multidradump, and starts with "./" */
  FILE *fp=NULL;
  int ndim=0,nd=0,length=0;
  if ((fp=fopen(filename,"r"))==NULL){ printf(" %% warning, couldn't open %s in multidraread\n",filename);}
  else /* if fp!=NULL */{
    fread(&ndim,sizeof(int),1,fp);
    *ira = (int *) tcalloc(ndim,sizeof(int));
    fread(*ira,sizeof(int),ndim,fp);
    length=1; for (nd=0;nd<ndim;nd++){ length*=(*ira)[nd];}
    *dra = (double *) tcalloc(length,sizeof(double));
    fread(*dra,sizeof(double),length,fp);
    if (fp!=NULL){ fclose(fp);fp=NULL;}}
  return ndim;
}

void * raread(char *filename,int *type,int *rows,int *cols)
{
  /* assume *filename dumped with radump, and starts with "./" */
  FILE *fp=NULL;
  void *ra=NULL;
  if ((fp=fopen(filename,"r"))==NULL){ printf(" %% Warning, cannot read %s in raread\n",filename);}
  else{
    fread(type,sizeof(int),1,fp);
    if (*type==0){ /* int */ 
      fread(rows,sizeof(int),1,fp); fread(cols,sizeof(int),1,fp);
      ra = tcalloc(*rows**cols,sizeof(int));
      fread(ra,sizeof(int),*rows**cols,fp);}
    else if (*type==1){ /* double */
      fread(rows,sizeof(int),1,fp); fread(cols,sizeof(int),1,fp);
      ra = tcalloc(*rows**cols,sizeof(double));
      fread(ra,sizeof(double),*rows**cols,fp);}
    else{ printf(" %% Warning, improper type %d in raread\n",*type);}
    fclose(fp);
    return ra;}
  return NULL;
}

int radump(void *ra,char *type,int rows,int cols,char *filename)
{
  /* assume filename starts with "./" */
  int output_flag=1;
  FILE *fp=NULL;
  int int_type=0;
  int double_type=1;
  if ((fp=fopen(filename,"w"))==NULL){ printf(" %% Warning, cannot open %s in radump, writing to stdout\n",filename); fp=stdout;}
  if (strcmp(type,"int")==0){
    output_flag *= fwrite(&int_type,sizeof(int),1,fp);
    output_flag *= fwrite(&rows,sizeof(int),1,fp);
    output_flag *= fwrite(&cols,sizeof(int),1,fp);
    output_flag *= fwrite(ra,sizeof(int),rows*cols,fp);}
  else if (strcmp(type,"double")==0){
    output_flag *= fwrite(&double_type,sizeof(int),1,fp);
    output_flag *= fwrite(&rows,sizeof(int),1,fp);
    output_flag *= fwrite(&cols,sizeof(int),1,fp);
    output_flag *= fwrite(ra,sizeof(double),rows*cols,fp);}
  else{ printf(" %% Warning, type %s in radump\n",type);}
  if (fp!=stdout){ fclose(fp);}
  return output_flag;
}

int rarecord(int datavsstring,double *ra,int length,char *filename)
{
  /* stupid double array dump (tacks onto filename, which starts with "./") */
  int output_flag=0;
  FILE *fp=NULL;
  int na=0;
  if ((fp=fopen(filename,"a"))==NULL) { printf(" %% Warning: cannot open %s in rarecord, writing to stdout\n",filename); fp=stdout;}
  if (datavsstring>=0){ output_flag = fwrite(ra,sizeof(double),length,fp);}
  else /* if (datavsstring<0) */{ output_flag = 1; for (na=0;na<length;na++){ output_flag *= fprintf(fp,"%0.16f\n",ra[na]);}}
  if (fp!=stdout){ fclose(fp);}
  return output_flag;
}

int pnm2mpg(char *basename,int framestart,int frameend)
{
  /* assume basename starts with "./" */
  char text[128],command[256],framestartname[16],frameendname[16];
  int output_flag=1;
  FILE *fp=NULL;
  num2frame(framestart,framestartname);
  num2frame(frameend,frameendname);
  if (ON_MY_COMPUTER==1){
    sprintf(command,"nconvert -out jpeg %s*.pnm;",basename); output_flag *= system(command);
    sprintf(command,"makempeg -fs %s -fe %s -fi 0001 -base %s -ext jpg > %slog;",framestartname,frameendname,basename,basename);
    output_flag *= system(command);
    sprintf(command,"rm %s*.jpg;",basename); output_flag *= system(command);}
  else if (ON_MY_COMPUTER==-1){ /* we are on stratum */
    sprintf(text,"%s.param",basename);
    if ((fp=fopen(text,"w"))==NULL){ printf("error, cannot open %s in pnm2mpg, writing to stdout\n",text); fp=stdout;}
    fprintf(fp,"PATTERN \t IBBPBBPBBPBBPBBP\n");
    fprintf(fp,"OUTPUT \t %s.mpg\n",basename);
    fprintf(fp,"BASE_FILE_FORMAT \t PNM\n");
    fprintf(fp,"INPUT_CONVERT *\n");
    fprintf(fp,"GOP_SIZE \t 16\n");
    fprintf(fp,"SLICES_PER_FRAME \t 1\n");
    fprintf(fp,"INPUT_DIR \t ./\n");
    fprintf(fp,"INPUT\n");
    fprintf(fp,"%s.*.pnm [%s-%s]\n",basename,framestartname,frameendname);
    fprintf(fp,"END_INPUT\n");
    fprintf(fp,"PIXEL \t HALF\n");
    fprintf(fp,"RANGE \t 10\n");
    fprintf(fp,"PSEARCH_ALG \t LOGARITHMIC\n");
    fprintf(fp,"BSEARCH_ALG \t CROSS2\n");
    fprintf(fp,"IQSCALE \t 8\n");
    fprintf(fp,"PQSCALE \t 10\n");
    fprintf(fp,"BQSCALE \t 25\n");
    fprintf(fp,"REFERENCE_FRAME ORIGINAL\n");
    fprintf(fp,"BIT_RATE 1000000\n");
    fprintf(fp,"BUFFER_SIZE 327680\n");
    fprintf(fp,"FRAME_RATE 30\n");
    if (fp!=stdout){ fclose(fp);}
    if ((fp=fopen(text,"r"))==NULL){ printf("error, could not read %s in pnm2mpg, doing nothing\n",text); output_flag=0;}
    else{ fclose(fp); sprintf(command,"/home/rangan/local/bin/mpeg_encode_stratum %s;",text); output_flag *= system(command);}}
  return output_flag;
}

double ra_norm(double *ra1,int length){ int nr=0;double output=0; for (nr=0;nr<length;nr++){ output+=pow(ra1[nr],2);} return sqrt(output);}

double ra2ra_dot(double *ra1,double *ra2,int length)
{
  int i=0; double c=0; for (i=0;i<length;i++){ c+=ra1[i]*ra2[i];} return finite(c) ? c : 0;
}

double * ra2ra_plus(double *ra1,double *ra2,int length)
{
  double *ra=NULL; int nr=0;
  ra = (double *) tcalloc(length,sizeof(double));
  for (nr=0;nr<length;nr++){ ra[nr]=ra1[nr]+ra2[nr];}
  return ra;
}

double * ra2ra_minus(double *ra1,double *ra2,int length)
{
  double *ra=NULL; int nr=0;
  ra = (double *) tcalloc(length,sizeof(double));
  for (nr=0;nr<length;nr++){ ra[nr]=ra1[nr]-ra2[nr];}
  return ra;
}

double * ra2ra_matrix_multiply(double *ra1,int rows1,int cols1,int transpose_flag1,double *ra2,int rows2,int cols2,int transpose_flag2)
{
  double *ra=NULL;
  int rows1t=0,cols1t=0,rows2t=0,cols2t=0,rows=0,cols=0,tab=0,tab1=0,tab2=0,nr=0,nc=0,nl=0;
  if (transpose_flag1){ rows1t=cols1; cols1t=rows1;} else{ rows1t=rows1; cols1t=cols1;}
  if (transpose_flag2){ rows2t=cols2; cols2t=rows2;} else{ rows2t=rows2; cols2t=cols2;}
  if (cols1t!=rows2t){ printf(" %% warning! improper matrix dimensions in ra2ra_matrix_multiply\n");}
  rows=rows1t; cols=cols2t;
  ra = (double *) tcalloc(rows*cols,sizeof(double));
  for (nr=0;nr<rows;nr++){ for (nc=0;nc<cols;nc++){
    tab=nr+nc*rows; ra[tab]=0;
    for (nl=0;nl<cols1t;nl++){ 
      tab1 = transpose_flag1 ? nl+nr*cols1t : nr+nl*rows; 
      tab2 = transpose_flag2 ? nc+nl*cols : nl+nc*rows2t; 
      ra[tab] += ra1[tab1]*ra2[tab2];}}}
  return ra;
}

void raplugin(double *ra,int rows,int cols,double *ra2,int rows2,int cols2,int rowstart,int colstart)
{
  /* insert ra2 into ra at entry (rowstart,colstart) */
  int nr=0,nc=0;
  if (rowstart+rows2>rows || colstart+cols2>cols){ printf(" warning! exceeding matrix dimensions in raplugin\n");}
  for (nr=0;nr<rows2;nr++){ for (nc=0;nc<cols2;nc++){ ra[rowstart+nr+(colstart+nc)*rows]=ra2[nr+nc*rows2];}}
}

double * raplugout(double *ra,int rows,int cols,int rowstart,int colstart,int rowlength,int collength)
{
  /* if ra is rows by cols, we pull out the subarray ra[rowstart:rowstart+rowlength-1,colstart:colstart+collength-1] */
  int nr=0,nc=0;
  double *sra=NULL;
  rowlength = minimum(rows,rowlength); collength = minimum(cols,collength);
  rowstart = periodize(rowstart,0,rows); colstart = periodize(colstart,0,cols);
  sra = (double *) tcalloc(rowlength*collength,sizeof(double));
  for (nr=0;nr<rowlength;nr++){ for (nc=0;nc<collength;nc++){
    sra[nr + nc*rowlength] = ra[periodize(rowstart+nr,0,rows) + periodize(colstart+nc,0,cols)*rows];}}
  return sra;
}

void raplusequals(double *ra1,int length,double *ra2)
{
  /* adds ra2 to ra1 */ int nr=0; for (nr=0;nr<length;nr++){ ra1[nr] += ra2[nr];}
}

void ratimesequals(double *ra,int length,double multiplier)
{
  /* multiplies ra by multiplier */ int nr=0; for (nr=0;nr<length;nr++){ ra[nr] *= multiplier;}
}

void raaddequals(double *ra,int length,double adder)
{
  /* adds ra to adder */ int nr=0; for (nr=0;nr<length;nr++){ ra[nr] += adder;}
}

double * ra2power(void *pforward,double *ra,int length,double *outra,int sqrt_flag,int log_flag)
{
  /* uses /usr/local/lib/libfftw3.a to determine power of ra
     if pforward is NULL, forward plan is created and destroyed.
     otherwise pforward is used 
     if outra=NULL, creates output, otherwise, uses *outra */
  int nr=0;
  int create_plans=0;
  double *output=NULL,*outra2=NULL;
  fftw_complex *cra=NULL,*craz=NULL;
  fftw_plan temp_plan_forward,*plan_forward=NULL;
  if (pforward==NULL){ create_plans=1;} else{ create_plans=0;}
  cra = (fftw_complex *) fftw_malloc(length*sizeof(fftw_complex));
  craz = (fftw_complex *) fftw_malloc(length*sizeof(fftw_complex));
  if (create_plans){
    temp_plan_forward = fftw_plan_dft_1d(length,cra,craz,-1,FFTW_ESTIMATE);
    plan_forward = &temp_plan_forward;}
  else{
    plan_forward = (fftw_plan *) pforward;}
  for (nr=0;nr<length;nr++){
    cra[nr][0] = ra[nr]; cra[nr][1] = 0;}
  fftw_execute_dft(*plan_forward,cra,craz);
  if (outra==NULL){ output = (double *) tcalloc(length,sizeof(double)); outra2=output;}
  else /* if (outra!=NULL) */{ outra2=outra;}
  switch (sqrt_flag + 2*log_flag){
  case 1: /* sqrt */ for (nr=0;nr<length;nr++){ outra2[nr] = sqrt(pow(craz[nr][0],2)+pow(craz[nr][1],2));} break;
  case 2: /* log */ for (nr=0;nr<length;nr++){ outra2[nr] = log(pow(craz[nr][0],2)+pow(craz[nr][1],2));} break;
  case 3: /* log(sqrt()) */ for (nr=0;nr<length;nr++){ outra2[nr] = log(sqrt(pow(craz[nr][0],2)+pow(craz[nr][1],2)));} break;
  default: /* neither */ for (nr=0;nr<length;nr++){ outra2[nr] = pow(craz[nr][0],2)+pow(craz[nr][1],2);} break;}
  for (nr=0;nr<length;nr++){ if (!finite(outra2[nr])){ outra2[nr]=0;}}
  if (create_plans){
    fftw_destroy_plan(temp_plan_forward);}
  free(cra);free(craz);
  return outra2;
}

void rara2corr(void *pforward,void *pbackwards,double **rara,int N,int length,double *autocorrelation,double *crosscorrelation)
{
  /* uses /usr/local/lib/libfftw3.a to determine power of ra
     if pforward or pbackwards are NULL, plans are created and destroyed.
     otherwise plans are used */
  int nr1=0,nr2=0,nl=0;
  int create_plan_forward=0,create_plan_backwards=0;
  fftw_complex *cra1=NULL,*cra1z=NULL,*cra2=NULL,*cra2z=NULL;
  fftw_plan temp_plan_forward,*plan_forward=NULL;
  fftw_plan temp_plan_backwards,*plan_backwards=NULL;
  double **rara2=NULL;
  double mean=0,stdev=0;
  double normalizer=0;
  if (pforward==NULL){ create_plan_forward=1;} else{ create_plan_forward=0;}
  if (pbackwards==NULL){ create_plan_backwards=1;} else{ create_plan_backwards=0;}
  cra1 = (fftw_complex *) fftw_malloc(length*sizeof(fftw_complex));
  cra1z = (fftw_complex *) fftw_malloc(length*sizeof(fftw_complex));
  cra2 = (fftw_complex *) fftw_malloc(length*sizeof(fftw_complex));
  cra2z = (fftw_complex *) fftw_malloc(length*sizeof(fftw_complex));
  if (create_plan_forward){
    temp_plan_forward = fftw_plan_dft_1d(length,cra1,cra1z,-1,FFTW_ESTIMATE);
    plan_forward = &temp_plan_forward;}
  else{ plan_forward = (fftw_plan *) pforward;}
  if (create_plan_backwards){
    temp_plan_backwards = fftw_plan_dft_1d(length,cra1,cra2,+1,FFTW_ESTIMATE);
    plan_backwards = &temp_plan_backwards;}
  else{ plan_backwards = (fftw_plan *) pbackwards;}
  rara2 = (double **) tcalloc(N,sizeof(double *));
  for (nr1=0;nr1<N;nr1++){ rara2[nr1] = (double *) tcalloc(length,sizeof(double));}
  for (nr1=0;nr1<N;nr1++){
    stats("double",rara[nr1],length,NULL,NULL,&mean,&stdev);
    for (nl=0;nl<length;nl++){ rara2[nr1][nl] = (rara[nr1][nl]-mean)/stdev; if (!finite(rara2[nr1][nl])){ rara2[nr1][nl]=0;}}}
  if (autocorrelation!=NULL){
    normalizer = length*length*N;
    for (nr1=0;nr1<N;nr1++){
      for (nl=0;nl<length;nl++){
	cra1[nl][0] = rara2[nr1][nl]; cra1[nl][1] = 0;}
      fftw_execute_dft(*plan_forward,cra1,cra1z);
      for (nl=0;nl<length;nl++){
	cra1[nl][0] = cra1z[nl][0]*cra1z[nl][0] + cra1z[nl][1]*cra1z[nl][1];
	cra1[nl][1] = cra1z[nl][0]*cra1z[nl][1] - cra1z[nl][1]*cra1z[nl][0];}
      fftw_execute_dft(*plan_backwards,cra1,cra2);
      for (nl=0;nl<length;nl++){ autocorrelation[nl] += cra2[nl][0]/normalizer;}}}
  if (crosscorrelation!=NULL){
    normalizer = length*length*N*(N-1);
    for (nr1=0;nr1<N;nr1++){ for (nr2=0;nr2<N;nr2++){
      if (nr1!=nr2){
	for (nl=0;nl<length;nl++){
	  cra1[nl][0] = rara2[nr1][nl]; cra1[nl][1] = 0;
	  cra2[nl][0] = rara2[nr2][nl]; cra2[nl][1] = 0;}
	fftw_execute_dft(*plan_forward,cra1,cra1z);
	fftw_execute_dft(*plan_forward,cra2,cra2z);
	for (nl=0;nl<length;nl++){
	  cra1[nl][0] = cra1z[nl][0]*cra2z[nl][0] + cra1z[nl][1]*cra2z[nl][1];
	  cra1[nl][1] = cra1z[nl][0]*cra2z[nl][1] - cra1z[nl][1]*cra2z[nl][0];}
	fftw_execute_dft(*plan_backwards,cra1,cra2);
	for (nl=0;nl<length;nl++){ crosscorrelation[nl] += cra2[nl][0]/normalizer;}}}}}
  if (create_plan_forward){ fftw_destroy_plan(temp_plan_forward);}
  if (create_plan_backwards){ fftw_destroy_plan(temp_plan_backwards);}
  free(cra1);free(cra1z);
  free(cra2);free(cra2z);
  for (nr1=0;nr1<N;nr1++){ if (rara2[nr1]!=NULL){ tfree(rara2[nr1]);rara2[nr1]=NULL;}}
  tfree(rara2);rara2=NULL;
}

void ping(){ printf("ping\n");} void pong(){ printf("pong\n");}

