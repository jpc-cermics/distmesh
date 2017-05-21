//   Example: (Uniform Mesh on Unit Sphere)
function y=fd(p); y=dsphere(p,0,0,0,1);endfunction;
[p,t]=distmeshsurface(fd, huniform,0.2,1.1*[-1,-1,-1;1,1,1]);
xclick();

//   Example: (Graded Mesh on Unit Sphere)
function y=fd(p) ;y=dsphere(p,0,0,0,1);endfunction;
function y=fh(p) ;y=0.05+0.5*dsphere(p,0,0,1,0);endfunction;
[p,t]=distmeshsurface(fd,fh,0.15,1.1*[-1,-1,-1;1,1,1]);
xclick();

//   Example: (Uniform Mesh on Torus)
function y=fd(p) ; y=(sum(p.^2,2)+.8^2-.2^2).^2-4*.8^2*(p(:,1).^2+p(:,2).^2);endfunction;
[p,t]=distmeshsurface(fd,huniform,0.1,[-1.1,-1.1,-.25;1.1,1.1,.25]);
xclick();

//   Example: (Uniform Mesh on Ellipsoid)
function y=fd(p) ;y=p(:,1).^2/4+p(:,2).^2/1+p(:,3).^2/1.5^2-1;endfunction;
[p,t]=distmeshsurface(fd,huniform,0.2,[-2.1,-1.1,-1.6; 2.1,1.1,1.6]);
xclick()
// trunc graphics 
xclear();trisurf(t([1:700,1500:$],:),p);

//   Example: (Uniform Mesh on Truncated Ellipsoid)
function y=fd(p);
  z1= p(:,1) >=0.5;
  y=p(:,1).^2/4+p(:,2).^2/1+p(:,3).^2/1.5^2-1;
  y(find(~z1))= %inf;
endfunction;

[p,t]=distmeshsurface(fd,huniform,0.2,[-2.1,-1.1,-1.6; 2.1,1.1,1.6]);
xclick()
xclear();trisurf(t([1:700,1500:$],:),p);

//   Example: (Uniform Mesh on Hyperboloid)
function y=fd(p)
  z1= p(:,1) <= 4 && p(:,2) <= 4 && p(:,3) <= 4;
  z2= p(:,1) >= -4 && p(:,2) >= -4 && p(:,3) >= -4;
  y=p(:,1).^2 - p(:,2).^2+p(:,3).^2 -2 ;
  y(find( ~z1 || ~z2 ))= %inf;
endfunction;

[p,t]=distmeshsurface(fd,huniform,0.4,5*[-1,-1,-1;1,1,1]);
xclick()

