function aabbdemo(id)
%AABBDEMO run a series of aabb-tree demos.
%   AABBDEMO(II) runs the II-th demo, where +1 <= II <= +3. 
%   Demos illustrate the functionallity of the MAKETREE
%   routine -- performing spatial queries on collections of
%   bounding-boxes.
%
%   See also MAKETREE, SCANTREE, FINDPTS, FINDRAY, FINDBOX

%   Darren Engwirda : 2014 --
%   Email           : engwirda@mit.edu
%   Last updated    : 20/12/2014

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

    ffid = fopen('test/airfoil.node_2d','r');
    data = fscanf(ffid,'%e,%e,%i \r\n');
    fclose(ffid);
    
    pp = [data(1:3:end), ...
          data(2:3:end)] ;
      
    ffid = fopen('test/airfoil.tria_2d','r');
    data = fscanf(ffid,'%u,%u,%u,%u \r\n');
    fclose(ffid);
    
    tt = [data(1:4:end), ...
          data(2:4:end), ...
          data(3:4:end)] ;
    tt = tt+1;
    
    
    bi = pp(tt(:,1),:); bj = pp(tt(:,1),:);
    for ii = 2 : size(tt,2)    
        bi = min(bi,pp(tt(:,ii),:)) ;
        bj = max(bj,pp(tt(:,ii),:)) ;
    end
    tr = maketree([bi,bj]);
    
    
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

    ffid = fopen('test/femur.node_3d','r');
    data = fscanf(ffid,'%e,%e,%e,%i,%e \r\n');
    fclose(ffid);
    
    pp = [data(1:5:end), ...
          data(2:5:end), ...
          data(3:5:end)] ;
      
    ffid = fopen('test/femur.tria_3d','r');
    data = fscanf(ffid,'%u,%u,%u,%i \r\n');
    fclose(ffid);
    
    tt = [data(1:4:end), ...
          data(2:4:end), ...
          data(3:4:end)] ;
    tt = tt+1;

    bi = pp(tt(:,1),:); bj = pp(tt(:,1),:);
    for ii = 2 : size(tt,2)    
        bi = min(bi,pp(tt(:,ii),:)) ;
        bj = max(bj,pp(tt(:,ii),:)) ;
    end
    tr = maketree([bi,bj]);
    
    
    fc = [.95,.95,.55];
    ec = [.25,.25,.25];
    
    figure; 
    subplot(1,2,1); hold on;
    patch('faces',tt,'vertices',pp,'facecolor',fc,...
        'edgecolor',ec,'facealpha',+1.);
    axis image off;
    set(gca,'units','normalized','position',[0.01,0.05,.48,.90]);
    view(-150,50);
    light; camlight;
    subplot(1,2,2); hold on;
    drawtree(tr);
    axis image off;
    view(-150,50);
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
'   every circle (an O(N^2) operation). \n'...
'   The "fast" algorithm relies on an aabb-tree to partit-\n'...
'   ion the data, and then computes the intersections loc-\n'...
'   ally (an approx. O(N*log(N) operation)). \n\n']);


    nc = +10000;
    np = +10000;

    pc = randcirc(nc,2,0.05);
    pi = randn(np,size(pc,2)-1);
    
    fprintf(1,...
'   "Slow" algorithm: \n');
    
    tic
    cj     = slowfindcirc(pc,pi);
    toc
    
    fprintf(1,...
'   "Fast" algorithm: \n');
    
    tic
   [ci,tr] = fastfindcirc(pc,pi);
    toc
   
    fprintf(1,...
'   Equivalent results? \n');
    
    same = true;
    for ip = +1 : nc
       
        li = sort(ci{ip}) ;
        lj = sort(cj{ip}) ;
        
        if (~isempty(li) && ...
            ~isempty(lj))
        if any(li-lj ~= +0)
            same = false;
        end
        end
    end
    if (same)
    fprintf(1,...
'   TRUE \n');
    end
    
    
    vp = ceil(linspace(+1,min(1000,np),1000));
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

function [cl] = slowfindcirc(pc,pi)
%SLOWFINDCIRC find the points enclosed by a set of circles.
%   [CL] = SLOWFINDCIRC(PC,PI) returns an NC-by-1 cell
%   array CL, representing the points in PI enclosed by each
%   circle in PC. Specifically, CL{II} is a list of points
%   enclosed by the II-th circle. 

    rs = pc(:,3).^2;

    cl = cell(size(pc,1),1);
    for ic = +1 : size(pc,1)
       
        di =(pi(:,1)-pc(ic,1)).^2 + ...
            (pi(:,2)-pc(ic,2)).^2 ;
           
        cl{ic} = find(di<=rs(ic)) ;
    end
    
end

function [cl,tr] = fastfindcirc(pc,pi)
%FASTFINDCIRC find the points enclosed by a set of circles.
%   [CL,TR] = FASTFINDCIRC(PC,PI) returns an NC-by-1 cell
%   array CL, representing the points in PI enclosed by each
%   circle in PC. Specifically, CL{II} is a list of points
%   enclosed by the II-th circle. Also returns the aabb-tree
%   TR used to form the query.

    nd = size (pc,2)-1;
    
%------------------------- compute set of aabb's for circles
    bb = zeros(size(pc,1),nd*+2);
    for id = +1 : nd
        bb(:,id+nd*0) = ...
            pc(:,id)-pc(:,nd+1) ;
        bb(:,id+nd*1) = ...
            pc(:,id)+pc(:,nd+1) ;
    end
%----------------------- compute aabb-tree for set of aabb's    
    tr = maketree(bb,[]);

%----------------------- compute tree-to-point mapping lists
    tm = findpts (tr,pi);
   
    rs = pc(:,3).^2 ;
   
%---------------------------- loop over non-empty tree nodes
    cl = cell(size(pc,1),1);
    for ip = +1 : size(tm.ii,1)
       
        ni = tm.ii(ip);
        
        cj = tr.ll{ni}; % list of circle per node
        pj = tm.ll{ip}; % list of points per node
        
        pk = pi(pj,:) ;
   
    %------------------------ find points within each circle
        for jc = +1 : length(cj)
            
            ic = cj(jc);
            
            di =(pk(:,1)-pc(ic,1)).^2 + ...
                (pk(:,2)-pc(ic,2)).^2 ;

            cl{ic} = pj(di<=rs(ic));
        end     
    end
    
end

function [pc] = randcirc(nc,nd,rr)
%RANDCIRC make a set of NC randomised d-circles in R^ND with
%mean radius RR.

    pc = [randn(nc,nd), rr*(rand(nc,1))];

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
            
        xc{ic} = pc(ic,3)*xx - pc(ic,1);
        yc{ic} = pc(ic,3)*yy - pc(ic,2);

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

