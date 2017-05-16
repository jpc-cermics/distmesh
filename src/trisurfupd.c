// Copyright (C) 2004-2012 Per-Olof Persson. See COPYRIGHT.TXT for details.

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <nsp/interf.h> 
#include <nsp/objects.h> 
#include <nsp/matrix.h> 

static const char mod3x1[3]={1,2,0};
static const char mod3x2[3]={2,0,1};

#define sqr(x) ((x)*(x))
#define dot2(a,b) (a[0]*b[0]+a[1]*b[1])
#define length(x) (sqrt(sqr(x[0])+sqr(x[1])+sqr(x[2])))
#define dot(x,y)  (x[0]*y[0]+x[1]*y[1]+x[2]*y[2])
#define cross(v1,v2,n) { n[0]=v1[1]*v2[2]-v1[2]*v2[1];			\
    n[1]=v1[2]*v2[0]-v1[0]*v2[2];					\
    n[2]=v1[0]*v2[1]-v1[1]*v2[0];}

static void tupdate(double *p,int *t,int *t2t, int *t2n, int np,int nt);

int int_trisurfupd(Stack stack, int rhs, int opt, int lhs)
{
  NspMatrix *P,*T,*T2t,*T2n;
  int nt,np;
  CheckRhs (4,4);
  CheckLhs (0,1);
  if ((T = GetRealMatCopyInt (stack, 1)) == NULLMAT) return RET_BUG;
  if ((T2t = GetRealMatCopyInt (stack, 2)) == NULLMAT) return RET_BUG;
  if ((T2n = GetRealMatCopy (stack, 3)) == NULLMAT) return RET_BUG;
  if ((P = GetRealMatCopy (stack, 4)) == NULLMAT) return RET_BUG;
  nt = T->n;
  np = T2n->n;
  tupdate(P->R,T->I,T2t->I, T2n->I,np,nt);
  MoveObj(stack,1,NSP_OBJECT(P));
  return 1;
}

#if 0
static double triarea(double *p1,double *p2,double *p3)
{
  double d12x=p2[0]-p1[0];
  double d12y=p2[1]-p1[1];
  double d13x=p3[0]-p1[0];
  double d13y=p3[1]-p1[1];
  return (d12x*d13y-d12y*d13x)/2.0;
}
#endif

static void trinormal3(double *p1,double *p2,double *p3,double *n)
{
  double d12[3]={ p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2] };
  double d13[3]={ p3[0]-p1[0], p3[1]-p1[1], p3[2]-p1[2] };
  cross(d12,d13,n);

  double nnorm=length(n);
  n[0]/=nnorm; n[1]/=nnorm; n[2]/=nnorm;
}

static double triqual3(double *p1,double *p2,double *p3)
{
  double n[3];
  double d12[3]={ p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2] };
  double d13[3]={ p3[0]-p1[0], p3[1]-p1[1], p3[2]-p1[2] };
  double d23[3]={ p3[0]-p2[0], p3[1]-p2[1], p3[2]-p2[2] };
  cross(d12,d13,n);
  double vol=length(n)/2.0;
  double den=dot(d12,d12)+dot(d13,d13)+dot(d23,d23);
  return 6.928203230275509*vol/den;
}

static void tupdate(double *p,int *t,int *t2t, int *t2n, int np,int nt)
{
  for (int t1=0; t1<nt; t1++)
    for (int n1=0; n1<3; n1++) {
      int t2=t2t[n1+3*t1];
      if (t2>=0) {
        double q1=triqual3(&p[3*t[0+3*t1]],&p[3*t[1+3*t1]],&p[3*t[2+3*t1]]);
        double q2=triqual3(&p[3*t[0+3*t2]],&p[3*t[1+3*t2]],&p[3*t[2+3*t2]]);
        double minqold=Min(q1,q2);
        if (minqold<0.9) {
          int n2=t2n[n1+3*t1];
          int tix11=mod3x1[n1];
          int tix12=mod3x2[n1];
          int tix21=mod3x1[n2];
          int tix22=mod3x2[n2];
          
          int newt[2][3]={{t[0+3*t1],t[1+3*t1],t[2+3*t1]},
                          {t[0+3*t2],t[1+3*t2],t[2+3*t2]}};
          
          // Swap edge
          newt[0][tix12]=newt[1][n2];
          newt[1][tix22]=newt[0][n1];
          
          double q3=triqual3(&p[3*newt[0][0]],&p[3*newt[0][1]],&p[3*newt[0][2]]);
          double q4=triqual3(&p[3*newt[1][0]],&p[3*newt[1][1]],&p[3*newt[1][2]]);
          double minqnew= Min(q3,q4);

          if (minqnew>minqold+.025) {
            int flip;
            double normal1[3],normal2[3],normal3[3],normal4[3];
            trinormal3(&p[3*t[0+3*t1]],&p[3*t[1+3*t1]],&p[3*t[2+3*t1]],normal1);
            trinormal3(&p[3*t[0+3*t2]],&p[3*t[1+3*t2]],&p[3*t[2+3*t2]],normal2);
            trinormal3(&p[3*newt[0][0]],&p[3*newt[0][1]],&p[3*newt[0][2]],normal3);
            trinormal3(&p[3*newt[1][0]],&p[3*newt[1][1]],&p[3*newt[1][2]],normal4);
            flip=dot(normal1,normal2)>0 && dot(normal3,normal4)>0;
            if (flip) {
              int nbt;
              char nbn;
              
              // Insert new triangles
              memcpy(t+3*t1,newt[0],3*sizeof(int));
              memcpy(t+3*t2,newt[1],3*sizeof(int));
              
              // Update t2t and t2n
              nbt=t2t[tix21+3*t2];
              nbn=t2n[tix21+3*t2];
              t2t[n1+3*t1]=nbt;
              t2n[n1+3*t1]=nbn;
              if (nbt>=0) {
                t2t[nbn+3*nbt]=t1;
                t2n[nbn+3*nbt]=n1;
              }
              
              nbt=t2t[tix11+3*t1];
              nbn=t2n[tix11+3*t1];
              t2t[n2+3*t2]=nbt;
              t2n[n2+3*t2]=nbn;
              if (nbt>=0) {
                t2t[nbn+3*nbt]=t2;
                t2n[nbn+3*nbt]=n2;
              }
              
              t2t[tix11+3*t1]=t2;
              t2n[tix11+3*t1]=tix21;
              t2t[tix21+3*t2]=t1;
              t2n[tix21+3*t2]=tix11;
            }
          }
        }
      }
    }
}
