//   Example: (Uniform Mesh on Unit Sphere)

function y=fd(p); y=dsphere(p,0,0,0,1);endfunction
[p,t]=distmeshsurface(fd, huniform,0.2,1.1*[-1,-1,-1;1,1,1]);

//   Example: (Graded Mesh on Unit Sphere)

function y=fd(p) ;y=dsphere(p,0,0,0,1);endfunction
function y=fh(p) ;y=0.05+0.5*dsphere(p,0,0,1,0);endfunction
[p,t]=distmeshsurface(fd,fh,0.15,1.1*[-1,-1,-1;1,1,1]);

//   Example: (Uniform Mesh on Torus)
function y=fd(p) ; y=(sum(p.^2,2)+.8^2-.2^2).^2-4*.8^2*(p(:,1).^2+p(:,2).^2);endfunction
[p,t]=distmeshsurface(fd,huniform,0.1,[-1.1,-1.1,-.25;1.1,1.1,.25]);

//   Example: (Uniform Mesh on Ellipsoid)

function y=fd(p) ;y=p(:,1).^2/4+p(:,2).^2/1+p(:,3).^2/1.5^2-1;endfunction
[p,t]=distmeshsurface(fd,huniform,0.2,[-2.1,-1.1,-1.6; 2.1,1.1,1.6]);
  
