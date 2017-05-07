//MESHDEMO2d Distmesh2d examples.
//   Copyright (C) 2004-2012 Per-Olof Persson. See COPYRIGHT.TXT for details.
// rand('state', 1); // Always the same results

function fstats(p,t) 
  fprintf('%d nodes, %d elements, min quality %.2f\n',  ...
	  size(p, 1), size(t, 1), min(simpqual(p, t)));
endfunction 

fprintf('Uniform Mesh on Unit Circle\n');

function y=fd(p)
  y=sqrt(sum(p .^ 2, 2)) - 1;
endfunction

[p,t] = distmesh2d(fd, huniform, 0.2, [-1,-1;1,1], []);

fstats(p, t);
fprintf('(press any key)\n\n');pause

// ----------------

fprintf('Rectangle with circular hole, refined at circle boundary\n');
function y=fd(p)
  y=ddiff(drectangle(p, -1, 1, -1, 1), dcircle(p, 0, 0, 0.5));
endfunction

function y=fh(p) 
  y=0.05 + 0.3 * dcircle(p, 0, 0, 0.5);
endfunction
[p,t] = distmesh2d(fd, fh, 0.05, [-1,-1;1,1], [-1,-1;-1,1;1,-1;1,1]);
fstats(p, t);
fprintf('(press any key)\n\n');pause

// ----------------

fprintf('Polygon\n');
pv = [-0.4,-0.5;0.4,-0.2;0.4,-0.7;1.5,-0.4;0.9,0.1;
      1.6,0.8;0.5,0.5;0.2,1;0.1,0.4;-0.7,0.7;-0.4,-0.5];
[p,t] = distmesh2d(dpoly, huniform, 0.1, [-1,-1;2,1], pv, pv);
fstats(p, t);
fprintf('(press any key)\n\n');pause

fprintf('Ellipse\n');
function y=fd(p) 
  y=p(:, 1) .^ 2 / 2 ^ 2 + p(:, 2) .^ 2 / 1 ^ 2 - 1;
endfunction
[p,t] = distmesh2d(fd, huniform, 0.2, [-2,-1;2,1], []);
fstats(p, t);
fprintf('(press any key)\n\n');pause

fprintf('Square, with size function point and line sources\n');
function y=fd(p) y=drectangle(p, 0, 1, 0, 1);
endfunction
function y=fh(p) 
  y=min(min(0.01 + 0.3 * abs(dcircle(p, 0, 0, 0)), ...
	  0.025 + 0.3 * abs( dpoly( p, [0.3,  0.7;0.7, 0.5]))),  0.15);
endfunction

[p,t] = distmesh2d(fd, fh, 0.01, [0,0;1,1], [0,0;1,0;0,1;1,1]);
fstats(p, t);
fprintf('(press any key)\n\n');pause

fprintf('NACA0012 airfoil\n');

hlead = 0.01;htrail = 0.04;hmax = 2;circx = 2;circr = 4;
a = .12 / .2 * [0.2969,-0.1260,-0.3516,0.2843,-0.1036];

function y=fd(p) 
  y=ddiff(dcircle(p, circx, 0, circr), (abs(p(:, 2)) - polyval( ...
                                                                    [ ...
                                                                    a( ...
                                                                    (5 : -1 : 2)), ...
                                                                    0],  ...
                                                                    p(:, 1))) .^ 2 - a( ...
                                                                    1) ^ 2 * p(:,  ...
                                                                    1));
endfunction
function y=fh(p) 
  y=min(min(hlead + 0.3 * dcircle(p, 0, 0, 0), htrail + 0.3 * dcircle( ...
                                                                    p, 1, 0,  ...
                                                                    0)),  ...
             hmax);
endfunction

fixx = 1 - htrail * cumsum(1.3 .^ ((0 : 4))');
fixy = a(1) * sqrt(fixx) + polyval([a((5 : -1 : 2)),0], fixx);
fix = [[circx + [-1,1,0,0] * circr;0,0,circr * [-1,1]]';0,0;1,0;fixx, ...
       fixy;fixx,-fixy];
box = [circx - circr,-circr;circx + circr,circr];
h0 = min([hlead,htrail,hmax]);

[p,t] = distmesh2d(fd, fh, h0, box, fix);

fstats(p, t);
fprintf('(press any key)\n\n');pause
