function simpplot(p,t,expr,bcol,icol,nodes,tris)
//   Copyright (C) 2004-2012 Per-Olof Persson. See COPYRIGHT.TXT for details.

  function [I]=ismember(A,B)
  // quick implementation
    I=zeros(size(A));
    for i=1:size(B,'*');
      I(find(A==B(i)))=1;
    end
  endfunction
  
  dim = size(p, 2);
  select dim 
   case 2
    if nargin < 4 || isempty(bcol),
      bcol = [.8,.9,1];
    end
    if nargin < 5 || isempty(icol),
      icol = [0,0,0];
    end
    if nargin < 6,
      nodes = 0;
    end
    if nargin < 7,
      tris = 0;
    end
    printf("Revoir facecolor et edgecolor\n");
    triplot(t, p(:, 1), p(:, 2));// 0 * p(:, 1))// 'facecolor', bcol, 'edgecolor','k');
    if %f then 
      if nodes == 1
	line(p(:, 1), p(:, 2), 'linest', 'none', 'marker', '.', 'col', icol, 'markers', 24);
      elseif nodes == 2
	for ip = 1: size(p, 1)
	  txtpars = {'fontname','times','fontsize',12};
	  text(p(ip, 1), p(ip, 2), num2str(ip), txtpars{:});
	end
      end
      if tris == 2
	for it = 1: size(t, 1)
	  pmid = mean(p(t(it, :), :), 1);
	  txtpars = {'fontname','times','fontsize',12,'horizontala','center'};
	  text(pmid(1), pmid(2), num2str(it), txtpars{:});
	end
      end
    end
    // view(2)
    //axis('equal')
    // axis('off')
    // ax = axis;axis(ax * 1.001);
   case 3
    if nargin < 4 || isempty(bcol),
      bcol = [.8,.9,1];
    end
    if nargin < 5 || isempty(icol),
      icol = [.9,.8,1];
    end
    if size(t, 2) == 4
      tri1 = surftri(p, t);
      if nargin > 2 && ~isempty(expr)
	incl = find(evstr(expr)); // eval -> evstr
	t = t(any(ismember(t, incl), 2), :);
	tri1 = tri1(any(ismember(tri1, incl), 2), :);
	tri2 = surftri(p, t);
	tri2 = setdiff(tri2, tri1, which= 'rows');
	trisurf(tri2,[p,zeros(size(p,1),1)]);
	hold('on');
      end
    else 
      tri1 = t;
      if nargin > 2 && ~isempty(expr)
	incl = find(eval(expr));
	tri1 = tri1(any(ismember(tri1, incl), 2), :);
      end
    end
    trisurf(tri1, p);
    hold('off')
    //set(h, 'facecolor', bcol, 'edgecolor', 'k');
    //axis('equal')
    //cameramenu
  else 
    error('Unimplemented dimension.');
  end
endfunction

