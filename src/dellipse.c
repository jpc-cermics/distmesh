// Copyright (C) 2004-2012 Per-Olof Persson. See COPYRIGHT.TXT for details.

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <nsp/interf.h> 
#include <nsp/objects.h> 
#include <nsp/matrix.h> 

#define sqr(x) ((x)*(x))
#define dot2(a,b) (a[0]*b[0]+a[1]*b[1])
#define length(x,y) (sqrt(sqr(x)+sqr(y)))

static int roots(double *pol,double *rr,double *ri, int n);
static double dellipse(double x0,double y0,double a,double b, int *ok);

extern void C2F(dhseqr)(char*,char*, int*, int*, int*,double*, int*,double*,
			double*,void*, int*,double*,void*,int*);

int int_dellipse(Stack stack, int rhs, int opt, int lhs)
{
  NspMatrix *P,*Axes,*Res;
  int ip, ok;
  CheckStdRhs (2,2);
  CheckLhs (0,1);
  if ((P = GetRealMat (stack, 1)) == NULLMAT) return RET_BUG;
  if ((Axes = GetRealMat (stack, 2)) == NULLMAT) return RET_BUG;
  CheckLength(NspFname(stack),2,Axes,2);
  if ((Res = nsp_matrix_create(NVOID,'r',P->m,1)) == NULLMAT) return RET_BUG;
  for (ip = 0; ip < P->m; ip++)
    {
      Res->R[ip]= dellipse( P->R[ip], P->R[ip+P->m], Axes->R[0],Axes->R[1],&ok);
      if ( ok == FAIL) return RET_BUG;
    }
  MoveObj(stack,1,NSP_OBJECT(Res));
  return 1;
}

static double dellipse(double x0,double y0,double a,double b, int *ok)
{
  double t1,t2,t4,t6,t8,t9,t11,t15,t16;
  t1=a*a; t2=b*b; t4=y0*y0; t6=x0*x0; t8=t1*t1;  
  t9=t2*t1; t11=t2*t2; t15=t8*t2; t16=t11*t1;
  double pol[5]={ 1.0, 2.0*t1+2.0*t2, -t2*t4-t1*t6+t8+4.0*t9+t11,
                  -2.0*t9*t6-2.0*t9*t4+2.0*t15+2.0*t16,
                  -t16*t6+t8*t11-t15*t4 };
  double rr[4],ri[4];
  
  if ( roots(pol,rr,ri,4) == FAIL)
    {
      *ok=FAIL;
      return 0.0;
    }
  else
    {
      *ok= OK;
    }
  double t=-HUGE_VAL;
  for (int i=0; i<4; i++)
    t=Max(t,rr[i]);

  double x=sqr(a)*x0/(t+sqr(a));
  double y=sqr(b)*y0/(t+sqr(b));
  double d=t*sqrt(sqr(x)/sqr(sqr(a))+sqr(y)/sqr(sqr(b)));

  return d;
}

static int roots(double *pol,double *rr,double *ri, int n)
{
  double *H,*work;
  int o=1;
  int info;
  char chE='E',chN='N';
  H=(double*) calloc(n*n,sizeof(double));
  work=(double*)calloc(n,sizeof(double));

  memset(H,0,n*n*sizeof(double));
  for (int i=0; i<n-1; i++)
    H[1+(n+1)*i]=1.0;
  for (int i=0; i<n; i++)
    H[n*i]=-pol[i+1]/pol[0];
  
  C2F(dhseqr)(&chE,&chN,&n,&o,&n,H,&n,rr,ri,0,&n,work,&n,&info);
  
  free(work); free(H);
  if (info!=0)
    {
      Scierror("Roots not found.\n");
      return FAIL;
    }
  return OK;
}
