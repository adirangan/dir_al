/* file_start */

/* Here are the matrix functions */

struct matrix * mmake()
{
  struct matrix *A;
  A = (struct matrix *) tmalloc(sizeof(struct matrix));
  (*A).rows = 0;
  (*A).cols = 0;
  (*A).mtrx = NULL;
  return A;
}

void minit0(struct matrix *A,int n,int m)
{
  int index1,index2;
  if ((n==0)||(m==0)) printf(" %% error, minit0 called with (n,m)==(%d,%d)\n",n,m);
  if (A != NULL){
    mtfree(A);}
  (*A).rows = n;
  (*A).cols = m;
  (*A).mtrx = (double *) tmalloc(sizeof(double)*n*m);
  for (index1=0;index1<(*A).rows;index1++)
    for (index2=0;index2<(*A).cols;index2++)
      sentry(A,index1,index2,0);
}

void minit1(struct matrix *A,int n,int m)
{
  int index1,index2;
  if ((n==0)||(m==0)) printf(" %% error, minit1 called with (n,m)==(%d,%d)\n",n,m);
  if (A != NULL){
    mtfree(A);}
  (*A).rows = n;
  (*A).cols = m;
  (*A).mtrx = (double *) tmalloc(sizeof(double)*n*m);
  for (index1 = 0; index1 < n; index1++)
    for (index2 = 0; index2 < m; index2++)
      sentry(A,index1,index2,(index1==index2 ? 1 : 0));
}

void minitr(struct matrix *A,int n,int m)
{
  int index1,index2;
  if ((n==0)||(m==0)) printf(" %% error, minitr called with (n,m)==(%d,%d)\n",n,m);
  (*A).rows = n;
  (*A).cols = m;
  (*A).mtrx = (double *) tmalloc(sizeof(double)*n*m);
  for (index1 = 0; index1 < n; index1++)
    for (index2 = 0; index2 < m; index2++)
      sentry(A,index1,index2, ((double) (rand() - RAND_MAX/2)) / ((double) RAND_MAX) );
}

void sentry(struct matrix *A,int j,int k,double d)
{
  if ((j>=(*A).rows)||(k>=(*A).cols)||(j<0)||(k<0))
    printf(" %% error, index(%d,%d) outside matrix dimensions(%d,%d) in sentry\n",j,k,A->rows,A->cols);
  (*A).mtrx[entry(j,k,(*A).rows)] = d;
}

double gentry(struct matrix *A,int j,int k)
{
  if ((j>=(*A).rows)||(k>=(*A).cols)||(j<0)||(k<0))
    printf(" %% error, index(%d,%d) outside matrix dimensions(%d,%d) in gentry\n",j,k,A->rows,A->cols);
  return (*A).mtrx[entry(j,k,(*A).rows)];
}

void mcopy(struct matrix *A,struct matrix *B)
{
  if (B!=A){
    mtfree(B);
    (*B).rows = (*A).rows;
    (*B).cols = (*A).cols;
    (*B).mtrx = (double *) tmalloc(sizeof(double) * (*B).rows * (*B).cols);
    memcpy(B->mtrx,A->mtrx,B->rows*B->cols*sizeof(double));}
}  

void mplugin(struct matrix *A,int ma1,int ma2,int na1,int na2,struct matrix *C,int mb,int nb)
{
  /* This plugs submatrix A(ma1:ma2,na1:na2) into the submatrix with first entry C(mb,nb) */
  int index1,index2;
  struct matrix *tmp;
  if ((ma1 < 0)||(na1 < 0))
    printf(" %% error, out of bound submatrix A(%d:%d,%d:%d) when A is only %d by %d in mplugin\n",ma1,ma2,na1,na2,A->rows,A->cols);
  if ((ma2 >= (*A).rows)||(na2 >= (*A).cols))
    printf(" %% error, out of bound submatrix A(%d:%d,%d:%d) when A is only %d by %d in mplugin\n",ma1,ma2,na1,na2,A->rows,A->cols);
  if ((mb < 0)||((mb + ma2 - ma1)>=(*C).rows))
    printf(" %% error, writing to out of bound rows C(%d:%d,:) when C is only %d by %d in mplugin\n",mb,mb+ma2-ma1,C->rows,C->cols);
  if ((nb < 0)||((nb + na2 - na1)>=(*C).cols))
    printf(" %% error, writing to out of bound cols C(:,%d:%d) when C is only %d by %d in mplugin\n",nb,nb+na2-na1,C->rows,C->cols);
  if (C==A)
    {
      tmp = mmake();
      mplusd(C,0,tmp);
      for (index1=ma1;index1<=ma2;index1++)
	for (index2=na1;index2<=na2;index2++)
	  sentry(tmp,mb+index1-ma1,nb+index2-na1,gentry(A,index1,index2));
      mcopy(tmp,C);
      mtfree(tmp);tfree(tmp);
    }
  else
    {
      for (index1=ma1;index1<=ma2;index1++)
	for (index2=na1;index2<=na2;index2++)
	  sentry(C,mb+index1-ma1,nb+index2-na1,gentry(A,index1,index2));
    }
}

void mplugout(struct matrix *A,int ma1,int ma2,int na1,int na2,struct matrix *C)
{
  /* this extracts submatrix A(ma1:ma2,na1:na2) and stores it as C */
  int index1,index2;
  struct matrix *tmp;
  if ((ma2 > (*A).rows)||(na2 > (*A).cols))
    printf(" %% error, copying out of bound elements in mplugout\n");
  if (C==A)
    {
      tmp = mmake();
      minit0(tmp,ma2-ma1+1,na2-na1+1);
      for (index1=ma1;index1<=ma2;index1++)
	for (index2=na1;index2<=na2;index2++)
	  sentry(tmp,index1-ma1,index2-na1,gentry(A,index1,index2));
      mcopy(tmp,C);
      mtfree(tmp);tfree(tmp);
    }
  else
    {
      mtfree(C);
      minit0(C,ma2-ma1+1,na2-na1+1);
      for (index1=ma1;index1<=ma2;index1++)
	for (index2=na1;index2<=na2;index2++)
	  sentry(C,index1-ma1,index2-na1,gentry(A,index1,index2));
    }
}

double mmax(struct matrix *A)
{
  int rowindex,colindex;
  double temp=0,maxsofar=0;
  maxsofar = gentry(A,0,0);
  for (rowindex=0;rowindex<A->rows;rowindex++){
    for (colindex=0;colindex<A->cols;colindex++){
      temp = gentry(A,rowindex,colindex);
      if (temp > maxsofar){ maxsofar = temp;}}}
  return maxsofar;
}

double mmin(struct matrix *A)
{
  int rowindex,colindex;
  double temp=0,maxsofar=0;
  maxsofar = gentry(A,0,0);
  for (rowindex=0;rowindex<A->rows;rowindex++){
    for (colindex=0;colindex<A->cols;colindex++){
      temp = gentry(A,rowindex,colindex);
      if (temp < maxsofar){ maxsofar = temp;}}}
  return maxsofar;
}

void mtrans(struct matrix *A,struct matrix *C)
{
  int index1,index2;
  struct matrix *tmp;
  if (C==A)
    {
      tmp = mmake();
      (*tmp).rows = (*A).cols;
      (*tmp).cols = (*A).rows;
      (*tmp).mtrx = (double *) tmalloc(sizeof(double) * (*tmp).rows * (*tmp).cols);
      for (index1=0;index1<(*A).rows;index1++)
	for (index2=0;index2<(*A).cols;index2++)
	  sentry(tmp,index2,index1,gentry(A,index1,index2));
      mcopy(tmp,C);
      mtfree(tmp);tfree(tmp);
    }
  else
    {
      mtfree(C);
      (*C).rows = (*A).cols;
      (*C).cols = (*A).rows;
      (*C).mtrx = (double *) tmalloc(sizeof(double) * (*C).rows * (*C).cols);
      for (index1=0;index1<(*A).rows;index1++)
	for (index2=0;index2<(*A).cols;index2++)
	  sentry(C,index2,index1,gentry(A,index1,index2));
    }
}

void mtimesm(struct matrix *A,struct matrix *B,struct matrix *C)
{
  int index1,index2,index3;
  struct matrix *tmp;
  if ((C==A)||(C==B))
    {
      if ((*A).cols != (*B).rows)
	printf(" %% error: mtimesm called with matrices of incompatible dimension.\n");
      tmp = mmake();
      minit0(tmp,(*A).rows,(*B).cols);
      for (index1=0;index1<(*tmp).rows;index1++)
	for (index2=0;index2<(*tmp).cols;index2++)
	  for (index3=0;index3<(*A).cols;index3++)
	    sentry(tmp,index1,index2,gentry(tmp,index1,index2) + gentry(A,index1,index3)*gentry(B,index3,index2));
      mcopy(tmp,C);
      mtfree(tmp);tfree(tmp);
    }
  else
    {
      if ((*A).cols != (*B).rows)
	printf(" %% error: mtimesm called with matrices of incompatible dimension.\n");
      mtfree(C);
      minit0(C,(*A).rows,(*B).cols);
      for (index1=0;index1<(*C).rows;index1++)
	for (index2=0;index2<(*C).cols;index2++)
	  for (index3=0;index3<(*A).cols;index3++)
	    sentry(C,index1,index2,gentry(C,index1,index2) + gentry(A,index1,index3)*gentry(B,index3,index2));
    }
}

void mtimesm_block(struct matrix *A,struct matrix *B,struct matrix *C,int n,int K)
{
  /* This function multiplies matrices in a pseudo-block fashion
     for example, if
     A = [A1;A2;A3;A4] is an n*K by n matrix (K=4)
     and
     B = [B1;B2;B3;B4] is an n*K by n matrix (K=4)
     then this function computes
     C = [A1B1;A2B2;A3B3;A4B4] an n*K by n matrix (K=4)
     if on the other hand we have
     A = [A] an n by n matrix
     or
     B = [B] an n by n matrix
     then this function computes
     C = [AB1;AB2;AB3;AB4] or [A1B;A2B;A3B;A4B] an n*K by n matrix (K=4)
  */
  int m,index1;
  struct matrix *tmp,*tmp1,*tmp2,*tmp3;
  tmp = mmake();
  tmp1 = mmake();
  tmp2 = mmake();
  tmp3 = mmake();
  m = (*B).cols;
  if ((*A).cols != n)
    printf(" %% error: mtimesm_block called with left matrix of incompatible column dimension\n");
  if (((*B).rows != n)&&((*B).rows != n*K))
    printf(" %% error: mtimesm_block called with right matrix of incompatible row dimension %d \n",B->rows);
  if ((C==A)||(C==B))
    {
      minit0(tmp,n*K,m);
      if (((*A).rows == n*K)&&((*B).rows == n*K))
	for (index1=0;index1<K;index1++)
	  {
	    mplugout(A,index1*n,index1*n+n-1,0,n-1,tmp1);
	    mplugout(B,index1*n,index1*n+n-1,0,m-1,tmp2);
	    mtimesm(tmp1,tmp2,tmp3);
	    mplugin(tmp3,0,n-1,0,m-1,tmp,index1*n,0);
	  }
      else if (((*A).rows == n)&&((*B).rows == n*K))
	for (index1=0;index1<K;index1++)
	  {
	    mplugout(B,index1*n,index1*n+n-1,0,m-1,tmp2);
	    mtimesm(A,tmp2,tmp3);
	    mplugin(tmp3,0,n-1,0,m-1,tmp,index1*n,0);
	  }
      else if (((*A).rows == n*K)&&((*B).rows == n))
	for (index1=0;index1<K;index1++)
	  {
	    mplugout(A,index1*n,index1*n+n-1,0,n-1,tmp1);
	    mtimesm(tmp1,B,tmp3);
	    mplugin(tmp3,0,n-1,0,m-1,tmp,index1*n,0);
	  }
      else if (((*A).rows == n)&&((*B).rows == n))
	mtimesm(A,B,tmp);
      else
	printf(" %% error, something is wrong in mtimesm_block\n");
      mcopy(tmp,C);
    }
  else
    {
      mtfree(C);minit0(C,n*K,m);
      if (((*A).rows == n*K)&&((*B).rows == n*K))
	for (index1=0;index1<K;index1++)
	  {
	    mplugout(A,index1*n,index1*n+n-1,0,n-1,tmp1);
	    mplugout(B,index1*n,index1*n+n-1,0,m-1,tmp2);
	    mtimesm(tmp1,tmp2,tmp3);
	    mplugin(tmp3,0,n-1,0,m-1,C,index1*n,0);
	  }
      else if (((*A).rows == n)&&((*B).rows == n*K))
	for (index1=0;index1<K;index1++)
	  {
	    mplugout(B,index1*n,index1*n+n-1,0,m-1,tmp2);
	    mtimesm(A,tmp2,tmp3);
	    mplugin(tmp3,0,n-1,0,m-1,C,index1*n,0);
	  }
      else if (((*A).rows == n*K)&&((*B).rows == n))
	for (index1=0;index1<K;index1++)
	  {
	    mplugout(A,index1*n,index1*n+n-1,0,n-1,tmp1);
	    mtimesm(tmp1,B,tmp3);
	    mplugin(tmp3,0,n-1,0,m-1,C,index1*n,0);
	  }
      else if (((*A).rows == n)&&((*B).rows == n))
	mtimesm(A,B,C);
      else
	printf(" %% error, something is wrong in mtimesm_block\n");
    }
  mtfree(tmp);tfree(tmp);
  mtfree(tmp1);tfree(tmp1);
  mtfree(tmp2);tfree(tmp2);
  mtfree(tmp3);tfree(tmp3);
}

void minteg_block(struct matrix *T,struct matrix *A,struct matrix *C,int n,int K)
{
  /* This matrix integrates A via integration matrix T in a pseudo-block fashion
     for example, if we have
     A = [A1;A2;A3;A4] an n*K by n matrix (K=4)
     and
     T = [T] a 1 by K matrix
     then this function computes
     C = [C]_{kl} = T_{j}Aj_{kl} an n by n matrix
     however, if 
     T = [T] is a K by K matrix
     then this function computes
     C = [C1;C2;C3;C4] an n*K by n matrix
     where
     Cj = [Cj]_{lm} = T_{jk}Ak_{lm} an n by n matrix
  */
  int index1,index2;
  struct matrix *tmp1,*tmp2,*tmp3;
  tmp1 = mmake();
  tmp2 = mmake();
  tmp3 = mmake();
  if ((*A).rows != n*K)
    printf(" %% error, A has incompatible number of rows in minteg_block\n");
  if ((*T).cols != K)
    printf(" %% error, T has incompatible number of columns in minteg_block\n");
  if ((*T).rows == 1)
    {
      minit0(tmp2,n,(*A).cols);
      for (index1=0;index1<K;index1++)
	{
	  mtfree(tmp1);
	  mplugout(A,index1*n,index1*n+n-1,0,(*A).cols-1,tmp1);
	  mtimesd(tmp1,gentry(T,0,index1),tmp1);
	  mplusm(tmp1,tmp2,tmp2);
	}
      mcopy(tmp2,C);
    }
  else if ((*T).rows == K)
    {
      minit0(tmp3,n*K,(*A).cols);
      for (index1=0;index1<K;index1++)
	{
	  minit0(tmp2,n,(*A).cols);
	  for (index2=0;index2<K;index2++)
	    {
	      mtfree(tmp1);
	      mplugout(A,index2*n,index2*n+n-1,0,(*A).cols-1,tmp1);
	      mtimesd(tmp1,gentry(T,index1,index2),tmp1);
	      mplusm(tmp1,tmp2,tmp2);
	    }
	  mplugin(tmp2,0,n-1,0,(*tmp2).cols-1,tmp3,index1*n,0);
	}
      mcopy(tmp3,C);
    }
  else
    printf(" %% error, someting is wrong with rows of T in minteg_block\n");
  mtfree(tmp1);tfree(tmp1);
  mtfree(tmp2);tfree(tmp2);
  mtfree(tmp3);tfree(tmp3);
}


void mtimesd(struct matrix *A,double d,struct matrix *C)
{
  int index1,index2;
  if (C==A){
    for (index1=0;index1<A->rows;index1++){ for (index2=0;index2<A->cols;index2++){ sentry(A,index1,index2,gentry(A,index1,index2)*d);}}}
  else
    {
      mtfree(C);
      (*C).rows=(*A).rows;
      (*C).cols=(*A).cols;
      (*C).mtrx = (double *) tmalloc(sizeof(double) * (*C).rows * (*C).cols);
      for (index1=0;index1<(*C).rows;index1++)
	for (index2=0;index2<(*C).cols;index2++)
	  sentry(C,index1,index2,gentry(A,index1,index2)*d);
    }
}	      

void mplusm(struct matrix *A,struct matrix *B,struct matrix *C)
{
  int verbose=0;
  int index1,index2;
  double temp1,temp2;
  struct matrix *tmp;
  if (verbose && ((A->rows != B->rows)||(A->cols != B->cols))) 
    printf(" %% warning, A = %d x %d and B = %d x %d in mplusm\n",A->rows,A->cols,B->rows,B->cols);
  if ((C==A)||(C==B))
    {
      tmp = mmake();
      (*tmp).rows= maximum((*A).rows,(*B).rows);
      (*tmp).cols= maximum((*A).cols,(*B).cols);
      (*tmp).mtrx = (double *) tmalloc(sizeof(double) * (*tmp).rows * (*tmp).cols);
      for (index1=0;index1<(*tmp).rows;index1++)
	for (index2=0;index2<(*tmp).cols;index2++)
	  {
	    temp1 = (((index1<(*A).rows)&&(index2<(*A).cols)) ? gentry(A,index1,index2) : 0);
	    temp2 = (((index1<(*B).rows)&&(index2<(*B).cols)) ? gentry(B,index1,index2) : 0);
	    sentry(tmp,index1,index2,temp1+temp2);
	  }
      mcopy(tmp,C);
      mtfree(tmp);tfree(tmp);
    }
  else
    {
      mtfree(C);
      (*C).rows= maximum((*A).rows,(*B).rows);
      (*C).cols= maximum((*A).cols,(*B).cols);
      (*C).mtrx = (double *) tmalloc(sizeof(double) * (*C).rows * (*C).cols);
      for (index1=0;index1<(*C).rows;index1++)
	for (index2=0;index2<(*C).cols;index2++)
	  {
	    temp1 = (((index1<(*A).rows)&&(index2<(*A).cols)) ? gentry(A,index1,index2) : 0);
	    temp2 = (((index1<(*B).rows)&&(index2<(*B).cols)) ? gentry(B,index1,index2) : 0);
	    sentry(C,index1,index2,temp1+temp2);
	  }
    }
}

void mplusm_block(struct matrix *A,struct matrix *B,struct matrix *C,int n,int K)
{
  /* this function is expressly intended for pseudo-block matrix addition of the form
     A = n by ?
     B = n*K by ?
     or
     A = n*K by ?
     B = n by ?
  */
  int index1,m;
  struct matrix *tmp1,*tmp2;
  tmp1 = mmake();
  tmp2 = mmake();
  m = (*A).cols;
  if (m != (*B).cols)
    printf(" %% error, A and B of incompatible column dimension in mplusm_block\n");
  if (((*A).rows == n)&&((*B).rows == n*K))
    {
      minit0(tmp2,n*K,(*A).cols);
      for (index1=0;index1<K;index1++)
	{
	  mtfree(tmp1);
	  mplugout(B,index1*n,index1*n+n-1,0,m-1,tmp1);
	  mplusm(A,tmp1,tmp1);
	  mplugin(tmp1,0,n-1,0,m-1,tmp2,index1*n,0);
	}
      mcopy(tmp2,C);
    }
  else if (((*A).rows == n*K)&&((*B).rows == n))
    {
      minit0(tmp2,n*K,(*A).cols);
      for (index1=0;index1<K;index1++)
	{
	  mtfree(tmp1);
	  mplugout(A,index1*n,index1*n+n-1,0,m-1,tmp1);
	  mplusm(tmp1,B,tmp1);
	  mplugin(tmp1,0,n-1,0,m-1,tmp2,index1*n,0);
	}
      mcopy(tmp2,C);
    }
  else
    printf(" %% error, something is wrong in mplusm_block\n");
  mtfree(tmp1);tfree(tmp1);
  mtfree(tmp2);tfree(tmp2);  
}

void msubtm(struct matrix *A,struct matrix *B,struct matrix *C)
{
  int verbose=0;
  int index1,index2;
  double temp1,temp2;
  struct matrix *tmp;
  if (verbose && ((A->rows != B->rows)||(A->cols != B->cols)))
    printf(" %% warning, A = %d x %d and B = %d x %d in msubtm\n",A->rows,A->cols,B->rows,B->cols);
  if ((C==A)||(C==B))
    {
      tmp = mmake();
      (*tmp).rows= maximum((*A).rows,(*B).rows);
      (*tmp).cols= maximum((*A).cols,(*B).cols);
      (*tmp).mtrx = (double *) tmalloc(sizeof(double) * (*tmp).rows * (*tmp).cols);
      for (index1=0;index1<(*tmp).rows;index1++)
	for (index2=0;index2<(*tmp).cols;index2++)
	  {
	    temp1 = (((index1<(*A).rows)&&(index2<(*A).cols)) ? gentry(A,index1,index2) : 0);
	    temp2 = (((index1<(*B).rows)&&(index2<(*B).cols)) ? gentry(B,index1,index2) : 0);
	    sentry(tmp,index1,index2,temp1-temp2);
	  }
      mcopy(tmp,C);
      mtfree(tmp);tfree(tmp);
    }
  else
    {
      mtfree(C);
      (*C).rows= maximum((*A).rows,(*B).rows);
      (*C).cols= maximum((*A).cols,(*B).cols);
      (*C).mtrx = (double *) tmalloc(sizeof(double) * (*C).rows * (*C).cols);
      for (index1=0;index1<(*C).rows;index1++)
	for (index2=0;index2<(*C).cols;index2++)
	  {
	    temp1 = (((index1<(*A).rows)&&(index2<(*A).cols)) ? gentry(A,index1,index2) : 0);
	    temp2 = (((index1<(*B).rows)&&(index2<(*B).cols)) ? gentry(B,index1,index2) : 0);
	    sentry(C,index1,index2,temp1-temp2);
	  }
    }
}

void mplusd(struct matrix *A,double d,struct matrix *C)
{
  int index1,index2;
  struct matrix *tmp;
  if (d==0)
    {
      if (C==A)
	{
	  /* do nothing */
	}
      else
	{
	  mtfree(C);
	  (*C).rows = (*A).rows;
	  (*C).cols = (*A).cols;
	  (*C).mtrx = (double *) tmalloc(sizeof(double) * (*C).rows * (*C).cols);
	  for (index1=0;index1<(*C).rows;index1++)
	    for (index2=0;index2<(*C).cols;index2++)
	      sentry(C,index1,index2,gentry(A,index1,index2));
	}
    }
  else
    {
      if (C==A)
	{
	  tmp = mmake();
	  (*tmp).rows = (*A).rows;
	  (*tmp).cols = (*A).cols;
	  (*tmp).mtrx = (double *) tmalloc(sizeof(double) * (*tmp).rows * (*tmp).cols);
	  for (index1=0;index1<(*tmp).rows;index1++)
	    for (index2=0;index2<(*tmp).cols;index2++)
	      sentry(tmp,index1,index2,gentry(A,index1,index2)+d);
	  mcopy(tmp,C);
	  mtfree(tmp);tfree(tmp);
	}
      else
	{
	  mtfree(C);
	  (*C).rows = (*A).rows;
	  (*C).cols = (*A).cols;
	  (*C).mtrx = (double *) tmalloc(sizeof(double) * (*C).rows * (*C).cols);
	  for (index1=0;index1<(*C).rows;index1++)
	    for (index2=0;index2<(*C).cols;index2++)
	      sentry(C,index1,index2,gentry(A,index1,index2)+d);
	}
    }
}

void mfeval(struct matrix *A, double (*fun)(double),struct matrix *C)
{
  int index1,index2;
  struct matrix *tmp;
  if (C==A)
    {
      tmp = mmake();
      (*tmp).rows=(*A).rows;
      (*tmp).cols=(*A).cols;
      (*tmp).mtrx = (double *) tmalloc(sizeof(double) * (*tmp).rows * (*tmp).cols);
      for (index1=0;index1<(*tmp).rows;index1++)
	for (index2=0;index2<(*tmp).cols;index2++)
	  sentry(tmp,index1,index2,(*fun)(gentry(A,index1,index2)));
      mcopy(tmp,C);
      mtfree(tmp);tfree(tmp);
    }
  else
    {
      mtfree(C);
      (*C).rows=(*A).rows;
      (*C).cols=(*A).cols;
      (*C).mtrx = (double *) tmalloc(sizeof(double) * (*C).rows * (*C).cols);
      for (index1=0;index1<(*C).rows;index1++)
	for (index2=0;index2<(*C).cols;index2++)
	  sentry(C,index1,index2,(*fun)(gentry(A,index1,index2)));
    }
}

void mvnoc(struct matrix *A,struct matrix *B,struct matrix *C)
{
  /* This function convolves two vectors *A and *B
     If *A and *B are not one dimensional matrices, this function treats them as such
     This corresponds to polynomial multiplication where
     P_{A} = \sumj A_{j} x^{j}
     P_{B} = \sumj B_{j} x^{j}
     P_{C} = \sumj C_{j} x^{j} = \sumj (\sum_{i=0}^{j} A_{i} B_{j-i}) x^{j}
  */
  int lengthA,lengthB,index1,index2;
  struct matrix *tmp;
  lengthA = (*A).rows * (*A).cols;
  lengthB = (*B).rows * (*B).cols;
  if ((C==A)||(C==B))
    {
      tmp = mmake();
      minit0(tmp,lengthA+lengthB-1,1);
      for (index1=0;index1<(*tmp).rows;index1++)
	for (index2=0;index2<=index1;index2++)
	  (*tmp).mtrx[index1] = (*tmp).mtrx[index1] + ((index2<lengthA)&&((index1-index2)<lengthB) ? ((*A).mtrx[index2] * (*B).mtrx[index1-index2]) : 0);
      mcopy(tmp,C);
      mtfree(tmp);tfree(tmp);
    }
  else
    {
      mtfree(C);
      minit0(C,lengthA+lengthB-1,1);
      for (index1=0;index1<(*C).rows;index1++)
	for (index2=0;index2<=index1;index2++)
	  (*C).mtrx[index1] = (*C).mtrx[index1] + ((index2<lengthA)&&((index1-index2)<lengthB) ? ((*A).mtrx[index2] * (*B).mtrx[index1-index2]) : 0);
    }
}

void mpint(struct matrix *A,struct matrix *C)
{
  int index1;
  struct matrix *tmp;
  if (C==A)
    {
      tmp = mmake();
      (*tmp).rows = 1 + ((*A).rows*(*A).cols);
      (*tmp).cols = 1;
      (*tmp).mtrx = (double *) tmalloc(sizeof(double) * (*tmp).rows * (*tmp).cols);
      (*tmp).mtrx[0] = 0;
      for (index1=1;index1<(*tmp).rows;index1++)
	(*tmp).mtrx[index1] = ((*A).mtrx[index1-1])/index1;
      mcopy(tmp,C);
      mtfree(tmp);tfree(tmp);
    }
  else
    {
      mtfree(C);
      (*C).rows = 1 + ((*A).rows*(*A).cols);
      (*C).cols = 1;
      (*C).mtrx = (double *) tmalloc(sizeof(double) * (*C).rows * (*C).cols);
      (*C).mtrx[0] = 0;
      for (index1=1;index1<(*C).rows;index1++)
	(*C).mtrx[index1] = ((*A).mtrx[index1-1])/index1;
    }
}

void mpder(struct matrix *A,struct matrix *C)
{
  int index1;
  struct matrix *tmp;
  if (((*A).rows*(*A).cols)==1)
    {
      mtfree(C);
      minit0(C,1,1);
    }
  else
    {
      if (C==A)
	{
	  tmp = mmake();
	  (*tmp).rows = ((*A).rows * (*A).cols) - 1;
	  (*tmp).cols = 1;
	  (*tmp).mtrx = (double *) tmalloc(sizeof(double) * (*tmp).rows * (*tmp).cols);
	  for (index1=0;index1<(*tmp).rows;index1++)
	    (*tmp).mtrx[index1] = (*A).mtrx[index1+1] * (index1+1);
	  mcopy(tmp,C);
	  mtfree(tmp);tfree(tmp);
	}
      else
	{
	  mtfree(C);
	  (*C).rows = ((*A).rows * (*A).cols) - 1;
	  (*C).cols = 1;
	  (*C).mtrx = (double *) tmalloc(sizeof(double) * (*C).rows * (*C).cols);
	  for (index1=0;index1<(*C).rows;index1++)
	    (*C).mtrx[index1] = (*A).mtrx[index1+1] * (index1+1);
	}
    }
}

double mpeval(struct matrix *A,double d)
{
  int index1;
  double temp1;
  temp1=0;
  for (index1=((*A).cols*(*A).rows)-1;index1>0;index1--)
    temp1 = d*((*A).mtrx[index1] + temp1);
  temp1 = temp1 + (*A).mtrx[0];
  return temp1;
}

void mplu(struct matrix *A,struct matrix *Pfinal,struct matrix *LUfinal)
{
  int m,index1,index2,index3,temp1;
  double temp2,temp3;
  struct matrix *P,*L,*U,*I;
  /* Here we store *L and *U in the same matrix *LUfinal.
     Of course, since LU might very well be A, 
     we construct L and U and overwrite later.
  */
  if ((*A).rows != (*A).cols)
    printf(" %% error, nonsquare matrix passed into mplu\n");
  m = (*A).rows;
  P = mmake();
  L = mmake();
  U = mmake();
  I = mmake();
  mtfree(P);
  minit0(P,m,1);
  for (index1=0;index1<m;index1++)
    sentry(P,index1,0,index1);
  minit1(L,m,m);
  mplusd(A,0,U);
  minit1(I,m,m);
  for (index1=0;index1<m-1;index1++)
    {
      temp1 = index1;
      for (index2=index1+1;index2<m;index2++)
	if (fabs(gentry(U,index2,index1))>fabs(gentry(U,temp1,index1)))
	  temp1 = index2;
      for (index2=index1;index2<m;index2++)
	{
	  temp2 = gentry(U,index1,index2);
	  sentry(U,index1,index2,gentry(U,temp1,index2));
	  sentry(U,temp1,index2,temp2);
	}
      for (index2=0;index2<index1;index2++)
	{
	  temp2 = gentry(L,index1,index2);
	  sentry(L,index1,index2,gentry(L,temp1,index2));
	  sentry(L,temp1,index2,temp2);
	}
      temp2 = gentry(P,index1,0);
      sentry(P,index1,0,gentry(P,temp1,0));
      sentry(P,temp1,0,temp2);
      for (index2=(index1+1);index2<m;index2++)
	{
	  sentry(L,index2,index1,gentry(U,index2,index1)/gentry(U,index1,index1));
	  temp3 = gentry(L,index2,index1);
	  for (index3=index1;index3<m;index3++)
	    sentry(U,index2,index3,gentry(U,index2,index3) - temp3*gentry(U,index1,index3));
	}
    }
  msubtm(L,I,L);
  mplusm(L,U,LUfinal);
  mplusd(P,0,Pfinal);
  mtfree(P);tfree(P);
  mtfree(L);tfree(L);
  mtfree(U);tfree(U);
  mtfree(I);tfree(I);
}

void mplusolve(struct matrix *P,struct matrix *LU, struct matrix *B,struct matrix *X)
{
  int m,n,index1,index2,index3;
  double temp1=1;
  struct matrix *b,*x,*tmp1,*tmp2,*tmpX;
  //printf("starting mplusolve ");memprintf(0);
  b = mmake();
  x = mmake();
  tmp1 = mmake();
  tmp2 = mmake();
  tmpX = mmake();
  m = (*B).rows;
  n = (*B).cols;
  for (index1=0;index1<m;index1++)
    temp1 = temp1 * gentry(LU,index1,index1);
  if (temp1==0)
    printf(" %% error, matrix LU may be singular in mplusolve\n");
  /* we go through and solve each row of *B */
  mtfree(tmpX);
  minit0(tmpX,m,n);
  for (index3=0;index3<n;index3++)
    {
      mtfree(b);mtfree(x);mtfree(tmp1);mtfree(tmp2);
      mplugout(B,0,m-1,index3,index3,b);
      /* first we set tmp1 = P^{-1} * b */
      minit0(tmp1,m,1);
      for (index1=0;index1<m;index1++)
	sentry(tmp1,index1,0,(*b).mtrx[(int) gentry(P,index1,0)]);
      /* then we solve L * tmp2 = tmp1 */
      mtfree(tmp2);
      minit0(tmp2,m,1);
      for (index1=0;index1<m;index1++)
	{
	  temp1 = gentry(tmp1,index1,0);
	  for (index2=0;index2<index1;index2++)
	    temp1 = temp1 - gentry(LU,index1,index2)*gentry(tmp2,index2,0);
	  sentry(tmp2,index1,0,temp1);
	}
      /* then we solve U * x = tmp2 */
      mtfree(x);
      minit0(x,m,1);
      for (index1=m-1;index1>=0;index1--)
	{
	  temp1 = gentry(tmp2,index1,0);
	  for (index2=index1+1;index2<m;index2++)
	    temp1 = temp1 - gentry(LU,index1,index2)*gentry(x,index2,0);
	  sentry(x,index1,0,temp1/gentry(LU,index1,index1));
	}
      mplugin(x,0,m-1,0,0,tmpX,0,index3);
    }
  mtfree(b);tfree(b);
  mtfree(x);tfree(x);
  mtfree(tmp1);tfree(tmp1);
  mtfree(tmp2);tfree(tmp2);
  mcopy(tmpX,X);
  mtfree(tmpX);tfree(tmpX);
  //printf("ending mplusolve ");memprintf(0);
}

void mplusolve_block(struct matrix *P,struct matrix *LU,struct matrix *B,struct matrix *X,int n,int K)
{
  /* This matrix solves the system (P\inv) LU X = B
     in a pseudo-block fashion.
     for details, see the comment in mtimesm_block
  */
  int index1,nb=(*B).cols;
  struct matrix *tmp1;
  struct matrix *tmp2;
  struct matrix *tmp3;
  struct matrix *tmp4;
  tmp1 = mmake();
  tmp2 = mmake();
  tmp3 = mmake();
  tmp4 = mmake();
  if ((*P).cols != 1)
    printf(" %% error, matrix P has incompatible column dimension in mplusolve_block\n");
  if ((*LU).cols != n)
    printf(" %% error, matrix LU has incompatible column dimension in mplusolve_block\n");
  if (((*P).rows != n)&&((*P).rows != n*K))
    printf(" %% error, matrix P has incompatible row dimension in mplusolve_block\n");
  if (((*LU).rows != n)&&((*LU).rows != n*K))
    printf(" %% error, matrix LU has incompatible row dimension in mplusolve_block\n");
  if ((*P).rows != (*LU).rows)
    printf(" %% error, matrices P and U have incompatible row dimension in mplusolve_block\n");
  if (((*B).rows != n)&&((*B).rows != n*K))
    printf(" %% error, matrix B has incompatible row dimension in mplusolve_block\n");
  minit0(X,n*K,nb);
  if (((*P).rows == n*K)&&((*B).rows == n*K))
    for (index1=0;index1<K;index1++)
      {
	mtfree(tmp1);mtfree(tmp2);mtfree(tmp3);mtfree(tmp4);
	mplugout(P,index1*n,index1*n+n-1,0,0,tmp1);
	mplugout(LU,index1*n,index1*n+n-1,0,n-1,tmp2);
	mplugout(B,index1*n,index1*n+n-1,0,nb-1,tmp3);
	mplusolve(tmp1,tmp2,tmp3,tmp4);
	mplugin(tmp4,0,n-1,0,nb-1,X,index1*n,0);
      }
  else if (((*P).rows == n*K)&&((*B).rows == n))
    for (index1=0;index1<K;index1++)
      {
	mtfree(tmp1);mtfree(tmp2);mtfree(tmp3);mtfree(tmp4);
	mplugout(P,index1*n,index1*n+n-1,0,0,tmp1);
	mplugout(LU,index1*n,index1*n+n-1,0,n-1,tmp2);
	mplusolve(tmp1,tmp2,B,tmp4);
	mplugin(tmp4,0,n-1,0,nb-1,X,index1*n,0);
      }
  else if (((*P).rows == n)&&((*B).rows == n*K))
    for (index1=0;index1<K;index1++)
      {
	mtfree(tmp1);mtfree(tmp2);mtfree(tmp3);mtfree(tmp4);
	mplugout(B,index1*n,index1*n+n-1,0,nb-1,tmp3);
	mplusolve(P,LU,tmp3,tmp4);
	mplugin(tmp4,0,n-1,0,nb-1,X,index1*n,0);
      }
  else if (((*P).rows == n)&&((*B).rows == n))
    mplusolve(P,LU,B,X);
  else
    printf(" %% error, something is wrong in mplusolve_block\n");
  mtfree(tmp1);tfree(tmp1);
  mtfree(tmp2);tfree(tmp2);
  mtfree(tmp3);tfree(tmp3);
  mtfree(tmp4);tfree(tmp4);
}

void minv(struct matrix *A,struct matrix *B)
{
  /* this inverts A and sticks the answer into B */
  struct matrix *I=mmake(),*P=mmake(),*LU=mmake();
  minit1(I,A->rows,A->cols);
  mplu(A,P,LU); mplusolve(P,LU,I,B);
  mtfree(I);mtfree(P);mtfree(LU);
  tfree(I);tfree(P);tfree(LU);
}

double mfrobnorm(struct matrix *A)
{
  int index1;
  double temp1;
  temp1=0;
  for (index1=0;index1<((*A).rows*(*A).cols);index1++)
    temp1 = temp1 + (*A).mtrx[index1]*(*A).mtrx[index1];
  temp1 = sqrt(temp1);
  return temp1;
}

void mqrh(struct matrix *A,struct matrix *vR)
{
  int m,n,index1,index2,index3;
  struct matrix *tmp1,*tmp2,*tmp3,*tmp4,*v,*R;
  /* Here we are using the householder QR factorization.
     This produces reflection vectors v_{1},\ldots,v_{n}
     as well as the matrix R.
     We place all the reflection vectors, as well as R into the matrix vR
     the upper triangular part of vR is equal to R
     the subdiagonal entries of the kth column of vR is equal to v_{k}
  */
  m = (*A).rows;n = (*A).cols; if (n>m) printf(" %% error, n>m in call to mqrh\n");
  tmp1 = mmake();
  tmp2 = mmake();
  tmp3 = mmake();
  tmp4 = mmake();
  v = mmake();
  minit0(v,1 + m,n);
  R = mmake();
  mplusd(A,0,R);
  for (index1=0;index1<n;index1++)
    {
      mtfree(tmp1);
      minit0(tmp1,m-index1,1);
      for (index2=index1;index2<m;index2++)
	sentry(tmp1,index2-index1,0,gentry(R,index2,index1));
      sentry(tmp1,0,0,mfrobnorm(tmp1)*(gentry(tmp1,0,0)>0 ? 1.0 : -1.0) + gentry(tmp1,0,0));
      mtimesd(tmp1,1.0/mfrobnorm(tmp1),tmp1);
      mtfree(tmp2);
      minit0(tmp2,m-index1,n-index1);
      for (index2=0;index2<m-index1;index2++)
	for (index3=0;index3<n-index1;index3++)
	  sentry(tmp2,index2,index3,gentry(R,index2+index1,index3+index1));
      mtfree(tmp3);
      mplusd(tmp2,0,tmp3);
      mtfree(tmp4);
      mtrans(tmp1,tmp4);
      mtimesd(tmp4,-2.0,tmp4);
      mtimesm(tmp4,tmp2,tmp2);
      mtimesm(tmp1,tmp2,tmp2);
      mplusm(tmp3,tmp2,tmp2);
      for (index2=0;index2<m-index1;index2++)
	for (index3=0;index3<n-index1;index3++)
	  sentry(R,index2+index1,index3+index1,gentry(tmp2,index2,index3));
      for (index2=0;index2<m-index1;index2++)
	sentry(v,1 + index1 + index2,index1,gentry(tmp1,index2,0));
    }
  mtfree(vR);
  mplusm(v,R,vR);
  mtfree(tmp1);tfree(tmp1);
  mtfree(tmp2);tfree(tmp2);
  mtfree(tmp3);tfree(tmp3);
  mtfree(tmp4);tfree(tmp4);
}

double mkget(struct matrix *A)
{
  /* This just gets the damn condition number */
  int n=0;
  double temp;
  struct matrix *tmp1,*tmp2,*Ainv;
  tmp1=mmake();
  tmp2=mmake();
  Ainv=mmake();
  if ((*A).rows != (*A).cols)
    printf(" %% error, nonsquare matrix passed to mkget\n");
  else
    n=(*A).rows;
  minit1(Ainv,n,n);
  mplu(A,tmp1,tmp2);    
  mplusolve(tmp1,tmp2,Ainv,Ainv);
  temp = (mfrobnorm(A)*mfrobnorm(Ainv));
  mtfree(tmp1);tfree(tmp1);
  mtfree(tmp2);tfree(tmp2);
  mtfree(Ainv);tfree(Ainv);
  return temp;
}

double mkest(struct matrix *A)
{
  /* This is Hager's condition number estimator for matrix A */
  int n,maxindex,iteration,stopflag,index1;
  double normA=0,normAinv=0;
  struct matrix *x,*W,*z,*ztx,*Pa,*LUa,*Pat,*LUat;
  x = mmake();
  W = mmake();
  z = mmake();
  ztx = mmake();
  Pa = mmake();
  LUa = mmake();
  Pat = mmake();
  LUat = mmake();
  if ((*A).rows != (*A).cols)
    printf(" %% error, nonsquare matrix passed to mkest\n");
  n = (*A).rows;
  /* First we estimate \|A\|_{1} */
  mtfree(x);minit0(x,n,1);
  for (index1=0;index1<n;index1++)
    sentry(x,index1,0,1.0/((double) n));
  iteration=0;
  stopflag=0;
  while ((stopflag==0)&&(iteration<(2*n)))
    {
      mtfree(W);mtfree(z);mtfree(ztx);
      mtimesm(A,x,W);
      minit0(z,n,1);
      for (index1=0;index1<n;index1++)
	sentry(z,index1,0,(gentry(W,index1,0) >= 0 ? 1.0 : -1.0));
      mtrans(z,z);mtimesm(z,A,z);mtrans(z,z);
      maxindex=0;
      for (index1=1;index1<n;index1++)
	maxindex = (fabs(gentry(z,index1,0)) > fabs(gentry(z,maxindex,0)) ? index1 : maxindex);
      mtrans(z,ztx);
      mtimesm(ztx,x,ztx);
      if (gentry(z,maxindex,0) <= gentry(ztx,0,0))
	{
	  normA = 0;
	  for (index1=1;index1<n;index1++)
	    normA = normA + fabs(gentry(W,index1,0));
	  stopflag = 1;
	}
      else
	{
	  mtfree(x);minit0(x,n,1);sentry(x,maxindex,0,1.0);
	  iteration++;
	}
    }
  if (iteration>(2*n))
    printf(" %% error, iteration exceeds 2*n in estimating normA for mkest\n");
  /* Now we estimate \|A^{-1}\|_{1} 
     We are sloppy here, and take 2 PLU decompositions...
     one of A and one of A^{T}
  */
  mplu(A,Pa,LUa);
  mtrans(A,x);
  mplu(x,Pat,LUat);
  mtfree(x);minit0(x,n,1);
  for (index1=0;index1<n;index1++)
    sentry(x,index1,0,1.0/((double) n));
  iteration=0;
  stopflag=0;
  while ((stopflag==0)&&(iteration<(2*n)))
    {
      mtfree(W);mtfree(z);mtfree(ztx);
      mplusolve(Pa,LUa,x,W);
      minit0(z,n,1);
      for (index1=0;index1<n;index1++)
	sentry(z,index1,0,(gentry(W,index1,0) >= 0 ? 1.0 : -1.0));
      mplusolve(Pat,LUat,z,z);
      maxindex=0;
      for (index1=1;index1<n;index1++)
	maxindex = (fabs(gentry(z,index1,0)) > fabs(gentry(z,maxindex,0)) ? index1 : maxindex);
      mtrans(z,ztx);
      mtimesm(ztx,x,ztx);
      if (gentry(z,maxindex,0) <= gentry(ztx,0,0))
	{
	  normAinv = 0;
	  for (index1=1;index1<n;index1++)
	    normAinv = normAinv + fabs(gentry(W,index1,0));
	  stopflag = 1;
	}
      else
	{
	  mtfree(x);minit0(x,n,1);sentry(x,maxindex,0,1.0);
	  iteration++;
	}
    }
  if (iteration>(2*n))
    printf(" %% error, iteration exceeds 2*n in estimating normAinv for mkest\n");
  /* dispose of all temporary matrices */
  mtfree(x);tfree(x);
  mtfree(W);tfree(W);
  mtfree(z);tfree(z);
  mtfree(ztx);tfree(ztx);
  mtfree(Pa);tfree(Pa);
  mtfree(LUa);tfree(LUa);
  mtfree(Pat);tfree(Pat);
  mtfree(LUat);tfree(LUat);
  /* Finally we compute the actual condition number */
  return(normA*normAinv);    
}

void mtfree(struct matrix *A)
{
  if (A != NULL)
    {
      if ((*A).mtrx != NULL)
	{
	  tfree((*A).mtrx);
	  (*A).mtrx = NULL;
	}
      (*A).rows = 0;
      (*A).cols = 0;
    }
}

void mprintf(struct matrix *A)
{
  int index1,index2;
  printf(" %% matrix stored at %d\n",(A==NULL ? 0 : (int) A));
  printf(" %% matrix rows and columns stored at %d,%d\n",(A==NULL ? 0 : (int) &(A->rows)),(A==NULL ? 0 : (int) &(A->cols)));
  printf(" %% matrix array address stored at %d\n",(A->mtrx==NULL ? 0 : (int) &(A->mtrx)));
  printf(" %% matrix entries stored at %d-%d\n",(A->mtrx==NULL ? 0 : (int) &(A->mtrx[0])),(A->mtrx==NULL ? 0 : (int) &(A->mtrx[A->rows*A->cols-1])));
  if (A->mtrx != NULL)
    {
      printf(" %% rows = %d, cols = %d \n",(*A).rows,(*A).cols);
      for (index1=0;index1<(*A).rows;index1++)
	{
	  printf(" %% | ");
	  for (index2=0;index2<(*A).cols;index2++)
	    {
	      printf("%0.2e",gentry(A,index1,index2));
	      printf("\t | ");
	    }
	  printf("\n");
	}
      printf("\n");
    }
  else
    printf(" %% no rows or columns yet\n");
}

void mspy(struct matrix *A,double epsilon)
{
  int index1,index2;
  printf("rows = %d, cols = %d \n",(*A).rows,(*A).cols);
  for (index1=0;index1<(*A).rows;index1++)
    {
      printf("|");
      for (index2=0;index2<(*A).cols;index2++)
        {
          if (fabs(gentry(A,index1,index2))>epsilon)
	    printf("#");
	  else
	    printf(" ");
        }
      printf("|\n");
    }
  printf("\n");
}

void mplotf(struct matrix *A,int ci,int nbins)
{
  /* quick-plot the values of column ci of A into nbins bins*/
  int rowindex=0,binindex=0,plot_bin=0,z_bin=0;
  double Acmax=0,Acmin=0,temp=0;
  struct matrix *Ac=mmake();
  minit0(Ac,A->rows,1);
  for (rowindex=0;rowindex<A->rows;rowindex++){
    sentry(Ac,rowindex,0,gentry(A,rowindex,ci));}
  Acmax = mmax(Ac);
  Acmin = mmin(Ac);
  z_bin = floor(nbins*(0-Acmin)/(Acmax-Acmin));
  printf(" %% column %d of matrix %d in interval [%e,%e]\n",ci,(int)A,Acmin,Acmax);
  for (rowindex=0;rowindex<Ac->rows;rowindex++){
    temp = gentry(Ac,rowindex,0);
    printf(" %%    ");
    for (binindex=0;binindex<=nbins;binindex++){
      plot_bin = floor(nbins*(temp-Acmin)/(Acmax-Acmin));
      if (binindex==plot_bin){ printf("*");}
      else if (binindex==z_bin){ printf("+");}
      else if (binindex!=plot_bin && binindex!=z_bin){ printf(".");}}
    printf("\n");}
  mtfree(Ac);tfree(Ac);
}

void aprintf(struct matrix *A)
{
  int index1;
  for (index1=0;index1<((*A).rows)*((*A).cols);index1++)
    {
      printf("index # %d has entry:  ",index1);
      printf("%f",(*A).mtrx[index1]);
      printf("\n");
    }
  printf("\n");
}

/* file_end */

