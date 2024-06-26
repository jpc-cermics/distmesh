function [t2t, t2n]=mkt2t(t)
// MKT2T  Compute element connectivities from element indices.
// [T2T,T2N]=MKT2T(T)
//   Copyright (C) 2004-2012 Per-Olof Persson. See COPYRIGHT.TXT for details.

// we have nt lines each lines describe a polyhedron 
// we create the edges 
  
  nt = size(t, 1);
  dim = size(t, 2) - 1;
  
  select dim 
   case 1
    edges = [t(:, 2);
	     t(:, 1)];
   case 2
    edges = [t(:, [2,3]);
	     t(:, [3,1]);
	     t(:, [1,2])];
   case 3
    edges = [t(:, [2,3,4]);
	     t(:, [3,4,1]);
	     t(:, [4,1,2]);
	     t(:, [1,2,3])];
  end
  
  ts = [repmat((1 : nt), 1, dim + 1); ...
        kron((1 : (dim + 1)), ones(1, nt))]';
  
  edges = sort(edges,type="c",dir="i");
  
  // Using mtlb=%t in unique 

  [foo,_foo,jx] = unique(edges, which = 'rows', mtlb_mode =%t);
  
  [jx,ix] = sort(jx, dir="i");
  ts = ts(ix, :);
  
  ix = find(diff(jx) == 0);
  ts1 = ts(ix, :);
  ts2 = ts(ix + 1, :);
  
  t2t = zeros(nt, dim + 1);
  t2t(ts1(:, 1) + nt * (ts1(:, 2) - 1)) = ts2(:, 1);
  t2t(ts2(:, 1) + nt * (ts2(:, 2) - 1)) = ts1(:, 1);
  
  if nargout >= 2
    t2n = zeros(nt, dim + 1);
    t2n(ts1(:, 1) + nt * (ts1(:, 2) - 1)) = ts2(:, 2);
    t2n(ts2(:, 1) + nt * (ts2(:, 2) - 1)) = ts1(:, 2);
  end
endfunction
