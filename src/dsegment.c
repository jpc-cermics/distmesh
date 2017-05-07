// Copyright (C) 2004-2012 Per-Olof Persson. See COPYRIGHT.TXT for details.

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <nsp/interf.h> 
#include <nsp/objects.h> 
#include <nsp/matrix.h> 

#define _p(i,j) p[(i)+np*(j)]
#define _pv(i,j) pv[(i)+nvs*(j)]
#define _ds(i,j) ds[(i)+np*(j)]

#define sqr(x) ((x)*(x))
#define dot2(a,b) (a[0]*b[0]+a[1]*b[1])
#define length(x,y) (sqrt(sqr(x)+sqr(y)))

int int_dsegment(Stack stack, int rhs, int opt, int lhs)
{
  NspMatrix *P,*Vs,*Res;
  int np,nvs, iv, ip;
  double *p, *pv,*ds;
  
  CheckRhs (2,2);
  CheckLhs (0,1);
  if ((P = GetRealMat (stack, 1)) == NULLMAT) return RET_BUG;
  if ((Vs = GetRealMat (stack, 2)) == NULLMAT) return RET_BUG;
  np = P->m;
  nvs = Vs->m;
  p = P->R;
  pv = Vs->R;

  if ((Res = nsp_matrix_create(NVOID,'r',np,nvs-1)) == NULLMAT) return RET_BUG;
  ds = Res->R;
  
  for ( iv=0; iv<nvs-1; iv++)
    {
      for ( ip=0; ip<np; ip++)
	{
	  double v[2]={_pv(iv+1,0)-_pv(iv,0),
		       _pv(iv+1,1)-_pv(iv,1)};
	  double w[2]={_p(ip,0)-_pv(iv,0),
		       _p(ip,1)-_pv(iv,1)};
	  double c1=dot2(v,w);
	  double c2=dot2(v,v);
	  if (c1<=0)
	    _ds(ip,iv)=length(_p(ip,0)-_pv(iv,0),_p(ip,1)-_pv(iv,1));
	  else if (c1>=c2)
	    _ds(ip,iv)=length(_p(ip,0)-_pv(iv+1,0),_p(ip,1)-_pv(iv+1,1));
	  else
	    _ds(ip,iv)=length(_p(ip,0)-(_pv(iv,0)+c1/c2*v[0]),
			     _p(ip,1)-(_pv(iv,1)+c1/c2*v[1]));
	}
    }
  MoveObj(stack,1,NSP_OBJECT(Res));
  return 1;
}
