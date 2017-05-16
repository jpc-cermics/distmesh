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
static double dellipsoid(double x0,double y0,double z0,double a,double b,double c, int *ok);
extern void C2F(dhseqr)(char*,char*,int*,int*,int*,double*,int*,double*,
			double*,void*,int*,double*,void*,int*);

int int_dellipsoid(Stack stack, int rhs, int opt, int lhs)
{
  NspMatrix *P,*Axes,*Res;
  int ip, ok;
  CheckStdRhs (2,2);
  CheckLhs (0,1);
  if ((P = GetRealMat (stack, 1)) == NULLMAT) return RET_BUG;
  if ((Axes = GetRealMat (stack, 2)) == NULLMAT) return RET_BUG;
  CheckLength(NspFname(stack),2,Axes,3);
  if ((Res = nsp_matrix_create(NVOID,'r',P->m,1)) == NULLMAT) return RET_BUG;
  for (ip = 0; ip < P->m; ip++)
    {
      Res->R[ip]= dellipsoid( P->R[ip], P->R[ip+P->m],P->R[ip+2*P->m], Axes->R[0],Axes->R[1],Axes->R[2],&ok);
      if ( ok == FAIL) return RET_BUG;
    }
  MoveObj(stack,1,NSP_OBJECT(Res));
  return 1;
}

static double dellipsoid(double x0,double y0,double z0,double a,double b,double c, int *ok)
{
  // Begin Maple Generated
  double t1,t2,t3,t5,t6,t7,t8,t9,t10,t11,t14,t15,
         t16,t17,t18,t20,t22,t24,t41,t42,t60;

  t1=a*a; t2=b*b; t3=c*c; t5=x0*x0; t6=t1*t5; t7=t1*t1; t8=t1*t2; t9=4.0*t8;
  t10=t2*t2; t11=t1+t2; t14=t3*t3; t15=y0*y0; t16=t2*t15; t17=z0*z0;
  t18=t3*t17; t20=t7*t2; t22=t1*t10; t24=t7+t9+t10; t41=t7*t10; t42=t20+t22;
  t60=t41+4.0*t42*t3+t24*t14-t18*t7-4.0*t18*t8-t18*t10-t16*t7-
      4.0*t16*t1*t3-t16*t14-t6*t10-4.0*t6*t2*t3-t6*t14;
  double pol[7]= {
    1.0,
    2.0*t1+2.0*t2+2.0*t3,
    -t6+t7+t9+t10+4.0*t11*t3+t14-t16-t18,
    2.0*t20+2.0*t22+2.0*t24*t3+2.0*t11*t14-2.0*t16*t1-2.0*t16*t3-2.0*t6*t2-2.0*t6*t3-2.0*t18*t1-2.0*t18*t2,
    t60,
    2.0*t41*t3+2.0*t42*t14-2.0*t18*t20-2.0*t18*t22-2.0*t16*t7*t3-2.0*t16*t1*t14-2.0*t6*t10*t3-2.0*t6*t2*t14,
    t41*t14-t6*t10*t14-t18*t41-t16*t7*t14
  };
  // End Maple Generated

  double rr[6],ri[6];
  if ( roots(pol,rr,ri,6)==FAIL)
    {
      *ok=FAIL;return 0.0;
    }
  else
    {
      *ok=OK;
    }
  double t=-HUGE_VAL;
  for (int i=0; i<6; i++)
    {
      t=Max(t,rr[i]);
    }

  double x=sqr(a)*x0/(t+sqr(a));
  double y=sqr(b)*y0/(t+sqr(b));
  double z=sqr(c)*z0/(t+sqr(c));
  double d=t*sqrt(sqr(x)/sqr(sqr(a))+sqr(y)/sqr(sqr(b))+sqr(z)/sqr(sqr(c)));

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
