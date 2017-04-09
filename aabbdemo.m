function aabbdemo(varargin)
%AABBDEMO run a series of aabb-tree demos.
%   AABBDEMO(II) runs the II-th demo, where +1 <= II <= +3. 
%   Demos illustrate the functionallity of the MAKETREE
%   routine -- performing spatial queries on collections of
%   bounding-boxes.
%
%   See also MAKETREE, DRAWTREE, QUERYSET, MAPVERT, MAPRECT

%   Darren Engwirda : 2014 --
%   Email           : engwirda@mit.edu
%   Last updated    : 08/04/2017

    if (nargin>=1) 
        id = varargin{1}; 
    else
        id =   +1;
    end

%------------------------------------------------- call demo
    switch (id)
        case 1, demo1;
        case 2, demo2;
        case 3, demo3;
        otherwise
    error('aabbdemo:invalidInput','Invalid demo selection.');
    end

end

function demo1
%-----------------------------------------------------------
    fprintf(1,[...
'   AABBTREE offers d-dimensional aabb-tree construction &\n'... 
'   search for collections of spatial objects. These trees\n'...
'   are useful when seeking to implement efficient spatial\n'...
'   queries -- determining intersections between collecti-\n'...
'   ons of spatial objects. \n\n'...
'   Given a collection of spatial objects, an aabb-tree p-\n'...
'   artitions the bounding-boxes of the elements in the c-\n'...
'   ollection (the aabb''s) into a "tree" (hierarchy) of  \n'...
'   rectangular "nodes". In contrast to other types of ge-\n'...
'   ometric trees (quadtrees, kd-trees, etc) the nodes in \n'...
'   an aabb-tree enclose aabb''s -- not points -- and may \n'...
'   overlap as a result. Objects in the collection are co-\n'...
'   ntained in a single node only. \n\n']);

    filename = mfilename('fullpath');
    filepath = fileparts( filename );

    addpath([filepath,'/mesh-file']);

   [geom] = readmsh([filepath,'/test-data/airfoil.msh']);
    
    pp = geom.point.coord(:,1:2);
    tt = geom.tria3.index(:,1:3);
   
    bi = pp(tt(:,1),:); bj = pp(tt(:,1),:);
    for ii = 2 : size(tt,2)    
        bi = min(bi,pp(tt(:,ii),:)) ;
        bj = max(bj,pp(tt(:,ii),:)) ;
    end
    
    tr = maketree([bi,bj]) ;
    
    fc = [.95,.95,.55];
    ec = [.25,.25,.25];
    
    figure; 
    subplot(1,2,1); hold on;
    patch('faces',tt,'vertices',pp,'facecolor',fc,...
        'edgecolor',ec,'facealpha',+.3);
    axis image off;
    set(gca,'units','normalized','position',[0.01,0.05,.48,.90]);
    subplot(1,2,2); hold on;
    drawtree(tr);
    axis image off;
    set(gca,'units','normalized','position',[0.51,0.05,.48,.90]);
    
end

function demo2
%-----------------------------------------------------------
    fprintf(1,[...
'   AABBTREE is a d-dimensional library, storing objects  \n'... 
'   and performing search operations in R^d. AABBTREE sim-\n'...
'   ply requires an description of the d-dimensional boun-\n'...
'   ding-boxes of a given collection. It is not limited to\n'...
'   simplexes (triangles, tetrahedrons, etc).\n\n']);

    filename = mfilename('fullpath');
    filepath = fileparts( filename );

    addpath([filepath,'/mesh-file']);

   [geom] = readmsh([filepath,'/test-data/veins.msh']);
    
    pp = geom.point.coord(:,1:3);
    tt = geom.tria3.index(:,1:3);
   
    bi = pp(tt(:,1),:); bj = pp(tt(:,1),:);
    for ii = 2 : size(tt,2)    
        bi = min(bi,pp(tt(:,ii),:)) ;
        bj = max(bj,pp(tt(:,ii),:)) ;
    end
    
    op.vtol = .67;
    
    tr = maketree([bi,bj],op) ;
    
    fc = [.95,.95,.55];
    ec = [.25,.25,.25];
    
    figure; 
    subplot(1,2,1); hold on;
    patch('faces',tt,'vertices',pp,'facecolor',fc,...
        'edgecolor',ec,'facealpha',+1.);
    axis image off;
    set(gca,'units','normalized','position',[0.01,0.05,.48,.90]);
    view(80,15);
    light; camlight;
    subplot(1,2,2); hold on;
    drawtree(tr);
    axis image off;
    view(80,15);
    light; camlight;
    set(gca,'units','normalized','position',[0.51,0.05,.48,.90]);
    
end

function demo3
%-----------------------------------------------------------
    fprintf(1,[...
'   AABBTREE facilitates efficient spatial queries through\n'...
'   "localisation" -- reducing a large O(M * N) comparison\n'...
'   to a sequence of small O(m * n) operations, where m<<M\n'...
'   and n<<N. By partitioning the data about the aabb-tree\n'...
'   itself, intersection tests can be carried out over sm-\n'...
'   local subsets, rather than between every pair of obje-\n'...
'   cts in the collection. \n\n'...
'   In the following example, the intersections between a \n'...
'   set of points and a set of circles is computed.       \n'...
'   The "slow" algorithm simply tests every point against \n'...
'   every circle (an O(N^2) operation). \n\n'...
'   The "fast" algorithm relies on an aabb-tree to partit-\n'...
'   ion the data, and then computes the intersections loc-\n'...
'   ally (an approx. O(N*log(N) operation)). The speed-up \n'...
'   is around a factor of 10 (on my machine).\n\n']);

    nc = +10000;
    np = +50000;

    pc = randcirc(nc,2,0.02);
    pi = rand(np,size(pc,2)-1);
    
    fprintf(1,...
'   "Slow" algorithm: \n');
    
    tic
   [ii_slow,ip_slow,cj_slow   ] = slowfindcirc(pc,pi);
    toc
    
    fprintf(1,...
'   "Fast" algorithm: \n');
    
    tic
   [ii_fast,ip_fast,cj_fast,tr] = fastfindcirc(pc,pi);
    toc
   
    fprintf(1,...
'   Equivalent results? \n');

    if (size(ii_slow) == size(ii_fast))
        
        same = true ;
        
        for ii = +1 : size(ip_slow,1)

            c1 = ...
        cj_slow(ip_slow(ii,1):ip_slow(ii,2));
            c2 = ...
        cj_fast(ip_fast(ii,1):ip_fast(ii,2));

            c1 = sort (c1) ;
            c2 = sort (c2) ;

            if (length(c1) == length(c2))
                if (any(c1 ~= c2))
                    same = false; break ; 
                end
            else
                same = false; break ;
            end

        end
        
    else
        same = false ;
    end
       
    if (same)
        fprintf(1,'   TRUE \n') ;
    else
        fprintf(1,'  FALSE \n') ;
    end
    
    vp = ceil(linspace(+1,min(2500,np),2500));
    vc = ceil(linspace(+1,min(5000,nc),5000));
    
    figure; 
    subplot(1,2,1); hold on;
    drawcirc(pc(vc,:));
    plot(pi(vp,1),pi(vp,2),'r.');
    axis image off;
    set(gca,'units','normalized','position',[0.01,0.05,.48,.90]);
    subplot(1,2,2); hold on;
    drawtree(tr);
    axis image off;
    set(gca,'units','normalized','position',[0.51,0.05,.48,.90]);

end

function [ii,ip,cj] = slowfindcirc(pc,pi)
%SLOWFINDCIRC find the points enclosed by a set of circles.
%   [II,IP,CJ] = SLOWFINDCIRC(PC,PI) computes the pairwise
%   point-circle intersections between the points PI and 
%   the circles PC. PC(:,1:2) are the circle centres and 
%   PC(:,3) are the circle radii. 
%   [II,IP,CJ] is the set of intersections in compressed
%   "sparse-style" indexing. Each point II(K) intersects
%   with the list of circles CJ(IP(K,1):IP(K,2)).
%
%   This is the "slow" brute-force variant.

    ip = 1:size(pi,1);
    ic = 1:size(pc,1);
    
   [pj,cj] = incircle(ip,ic,pi,[pc(:,1:2),pc(:,3).^2]);
   
    pj = pj( :);
    cj = cj( :);
   
%-- re-index to the sparse-style representation 
%-- of QUERYSET
   [pj,ix] = sort (pj) ; 
    cj = cj(ix);
    ix = find(diff(pj)>+0) ;
    
%-- the points in II intersect with >= 1 circle
    ni = length (pj) ;
    ii = pj([ix;ni]) ;
    
    nj = length (cj) ;
    ni = length (ii) ;
    
%-- the points in II intersect with the circles 
%-- CJ(IP(K,1):IP(K,2)) {for point II(K)}.
    ip = zeros(ni,2) ;
    ip(:,1) = [+1; ix+1] ;
    ip(:,2) = [ix; nj+0] ;
   
end

function [ii,ip,cj,tr] = fastfindcirc(pc,pi)
%FASTFINDCIRC find the points enclosed by a set of circles.
%   [II,IP,CJ] = FASTFINDCIRC(PC,PI) computes the pairwise
%   point-circle intersections between the points PI and 
%   the circles PC. PC(:,1:2) are the circle centres and 
%   PC(:,3) are the circle radii. 
%   [II,IP,CJ] is the set of intersections in compressed
%   "sparse-style" indexing. Each point II(K) intersects
%   with the list of circles CJ(IP(K,1):IP(K,2)).
%
%   This is the "fast" aabb-indexed variant.

    nd = size (pc,2)-1;
    
%-- compute the set of aabb's for circ.
    bb = zeros(size(pc,1),nd*+2);
    for id = +1 : nd
        bb(:,id+nd*0) = ...
            pc(:,id)-pc(:,nd+1) ;
        bb(:,id+nd*1) = ...
            pc(:,id)+pc(:,nd+1) ;
    end
%-- compute aabb-tree for set of aabb's    
    
    op.nobj = +64;

    tr = maketree(bb,op);

%-- compute tree-to-vert. indexing maps
    tm = mapvert (tr,pi);
   
%-- the points in II intersect with >= 1 circle
%-- the points in II intersect with the circles 
%-- CJ(IP(K,1):IP(K,2)) {for point II(K)}.
   [ii,ip,cj] = ...
    queryset(tr,tm,@incircle,pi,[pc(:,1:2),pc(:,3).^2]) ;
    
end

function [pj,cj] = incircle(ip,ic,pi,pc)
%INCIRCLE pairwise point-circle comparison kernel function.
%   [PJ,CJ] = INCIRCLE(IP,IC,PI,PC) compute the pairwise in-
%   tersections between the points PI(IP,:) and the circles
%   PC(IC,:). PC(:,3) are the squared circle radii.
%   [PJ,CJ] are pairs of intersections, such that the point 
%   PJ(K) intersects with the circle CJ(K).

    li = cell(length(ic),1);
    lj = cell(length(ic),1);

    pt = pi(ip,:) ;
    
    for ii = +1 : length(ic)
       
        di = (pt(:,1)-pc(ic(ii),1)).^2 + ...
             (pt(:,2)-pc(ic(ii),2)).^2 ;
           
        li{ii} = find(di<=pc(ic(ii),3));
        lj{ii} = ...
            ii * ones(length(li{ii}),1);
        
    end

    pj = ip(vertcat(li{:})) ;
    cj = ic(vertcat(lj{:})) ;
   
end

function [pc] = randcirc(nc,nd,rr)
%RANDCIRC make a set of NC randomised d-circles in R^ND with
%mean radius RR.

    pc = [rand(nc,nd), rr*(rand(nc,1))];

end

function drawcirc(pc)
%DRAWCIRC draw a set of NC d-circles.

    fc = [.95,.95,.55];
    ec = [.25,.25,.25];

    switch (size(pc,2))
        case 3
%--------------------------------------------------- circles
    
        tt = linspace(+.0,+2.*pi,24);
        xx = cos(tt)';
        yy = sin(tt)';
        nt = length(tt);
            
        ee = [(1:nt-1)',(2:nt-0)'];
        
        xc = cell(size(pc,1),1);
        yc = cell(size(pc,1),1);
        jc = cell(size(pc,1),1);
        
        for ic = +1 : size(pc,1)
            
        xc{ic} = pc(ic,3)*xx + pc(ic,1);
        yc{ic} = pc(ic,3)*yy + pc(ic,2);

        jc{ic} = ee + (ic-1)*nt;
        end
  
        pp =[vertcat(xc{:}), ...
             vertcat(yc{:})] ;
        ff = vertcat(jc{:});
        
        case 4
%--------------------------------------------------- spheres
    %%!!todo:

    end

    patch('faces',ff,'vertices',pp,'facecolor',fc,...
        'edgecolor',ec,'facealpha',+.3);
    
end

