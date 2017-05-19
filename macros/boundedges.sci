function e=boundedges(p,t)
//BOUNDEDGES Find boundary edges from triangular mesh
//   E=BOUNDEDGES(P,T)
  
//   Copyright (C) 2004-2012 Per-Olof Persson. See COPYRIGHT.TXT for details.
// Form all edges, non-duplicates are boundary edges

  edges = [t(:, [1,2]);
           t(:, [1,3]);
           t(:, [2,3])];
  node3 = [t(:, 3);t(:, 2);t(:, 1)];
  edges = sort(edges, type="c",dir="i");
  // unique + hist c replace by nsp unique with occ 
  // [foo,ix,jx] = unique(edges, which = 'rows', mtlb_mode=%t);
  // vec = histc(jx, (1 : max(jx)));
  [foo,ix,vec] = unique(edges, which = 'rows');
  qx = find(vec == 1);
  e = edges(ix(qx), :);
  node3 = node3(ix(qx));
  // Orientation
  v1 = p(e(:, 2), :) - p(e(:, 1), :);
  v2 = p(node3, :) - p(e(:, 1), :);
  ix = find(v1(:, 1)  .* v2(:, 2) - v1(:, 2)  .* v2(:, 1) > 0);
  e(ix, [1,2]) = e(ix, [2,1]);
endfunction
