%% swan_structured_hsv.m
% Coons boundary -> centers, OPTIONAL Lloyd with fixed boundary samples,
% harmonic (Tutte) parameterization -> structured nx-by-ny grid,
% plot gridlines on black background with HSV color by (x,y).

clear; clc; close all

%% ---------------- knobs ----------------
nx = 40; ny = 40;               % structured resolution
nEdge = 800;                    % points to draw the swan boundary polygon
edgeTol = 1e-12; EPS = 1e-12;
doSmoothing = true;            % smooth the preimage stair boundaries
doLloyd     = true;            % run Lloyd on sites (with fixed boundary samples)
maxIterL    = 30;              % modest Lloyd iters
omega       = 0.85;            % Lloyd under-relaxation

%% ---------------- swan boundary polygon (x,y) ----------------
[bx, by] = boundaryPolygon(nEdge);
swan = polyshape(bx,by);
xmin=min(bx); xmax=max(bx); ymin=min(by); ymax=max(by);
pad  = 2*max(xmax-xmin, ymax-ymin);
far  = [xmin-pad ymin-pad; xmax+pad ymin-pad; xmax+pad ymax+pad; xmin-pad ymax+pad];

%% ---------------- original mapped quads -> keep mask ----------------
u = linspace(0,1,nx+1); v = linspace(0,1,ny+1);
cells = struct('ij',{},'P',{},'area',{},'bbox',{});
k = 0;
for i=1:nx
  for j=1:ny
    P = [ coons(u(i),v(j));
          coons(u(i+1),v(j));
          coons(u(i+1),v(j+1));
          coons(u(i),v(j+1)) ];
    A = polyAreaSigned(P);
    if A < -EPS, P = flipud(P); A = -A; end
    if A <= EPS, continue; end
    if selfIntersectQuad(P, EPS), continue; end
    C = mean(P,1);
    if ~isinterior(swan,C(1),C(2)), continue; end
    k=k+1; cells(k).ij=[i j]; cells(k).P=P; cells(k).area=A; cells(k).bbox=bbox(P);
  end
end

% prune overlaps among non-neighbors
validC = true(1,numel(cells));
for a=1:numel(cells)
  if ~validC(a), continue; end
  Pa=cells(a).P; Ba=cells(a).bbox; ia=cells(a).ij(1); ja=cells(a).ij(2);
  for b=a+1:numel(cells)
    if ~validC(b), continue; end
    ib=cells(b).ij(1); jb=cells(b).ij(2);
    if areNeighbors(ia,ja,ib,jb), continue; end
    Bb=cells(b).bbox;
    if ~bboxesOverlap(Ba,Bb,EPS), continue; end
    if polysIntersectStrict(Pa,cells(b).P,EPS)
      if cells(a).area >= cells(b).area, validC(b)=false; else, validC(a)=false; break; end
    end
  end
end
cells = cells(validC);

% keep mask in (u,v)
keep = false(ny,nx);
for t=1:numel(cells), ij=cells(t).ij; keep(ij(2),ij(1))=true; end

%% ---------------- preimage boundary (raw -> optional smooth) -----------
[northRaw, southRaw, eastRaw, westRaw] = extract_raw_boundary_lists(keep,u,v,edgeTol);
if doSmoothing
  north = smooth_turns(northRaw, 1e-12, edgeTol);
  south = smooth_turns(southRaw, 1e-12, edgeTol);
  east  = smooth_turns(eastRaw,  1e-12, edgeTol);
  west  = smooth_turns(westRaw,  1e-12, edgeTol);
else
  north = northRaw; south = southRaw; east = eastRaw; west = westRaw;
end

%% ---------------- structured Coons grid in (u,v) -> map to (x,y) -------
xi  = linspace(0,1,nx+1);
eta = linspace(0,1,ny+1);
[Xi,Eta] = meshgrid(xi,eta);
[U,V] = coons_uv_from_boundary(Xi,Eta, south, east, north, west); % (ny+1)x(nx+1)
[X,Y] = coons_map_array(U,V);

% centers: exactly nx*ny
centers = zeros(nx*ny,2); idx=0;
for j=1:ny
  for i=1:nx
    idx = idx+1;
    P = [ X(j,i)     Y(j,i)    ;
          X(j,i+1)   Y(j,i+1)  ;
          X(j+1,i+1) Y(j+1,i+1);
          X(j+1,i)   Y(j+1,i)  ];
    centers(idx,:) = mean(P,1);
  end
end
fprintf('Structured centers: %d (expected %d)\n', size(centers,1), nx*ny);

% quick visualization of the original Coons grid and centers
figure('Color','w'); tl=tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
nexttile; hold on
plot(bx,by,'k-','LineWidth',1.5)
plot_gridlines(X,Y);
axis equal tight; box on
title(sprintf('Structured Coons grid — %dx%d cells',nx,ny)); xlabel x; ylabel y;

nexttile; hold on
plot(bx,by,'k-','LineWidth',1.5)
scatter(centers(:,1),centers(:,2),6,'k','filled')
axis equal tight; box on
title('Cell centers (nx*ny)'); xlabel x; ylabel y;

%% ---------------- OPTIONAL: Lloyd on centers with fixed boundary samples
if doLloyd
  % boundary samples (exclude corners to avoid duplicates)
  sH = linspace(0,1,nx+2); sH = sH(2:end-1);
  sV = linspace(0,1,ny+2); sV = sV(2:end-1);

  northPts = mk_samples(@top,    sH);   % nx x 2
  southPts = mk_samples(@bottom, sH);   % nx x 2
  eastPts  = mk_samples(@right,  sV);   % ny x 2
  westPts  = mk_samples(@left,   sV);   % ny x 2

  boundarySites = [northPts; southPts; eastPts; westPts];

  sites   = [centers; boundarySites];
  isFixed = false(size(sites,1),1);
  isFixed(size(centers,1)+1:end) = true;

  % deduplicate if any centers lie exactly on boundary
  [sites, ia] = unique(round(sites,12), 'rows', 'stable');
  isFixed = isFixed(ia);

  for it = 1:maxIterL
    [Vv,Cv] = voronoin([sites; far]);  % padded to avoid unbounded cells
    newSites = sites; areas = nan(size(sites,1),1);
    for p = 1:size(sites,1)
      if isFixed(p), continue, end
      id = Cv{p};
      if isempty(id) || any(id==1) || any(isinf(id)), continue, end
      Pcell = polyshape(Vv(id,1), Vv(id,2),'Simplify',true);
      Q = intersect(Pcell, swan);
      if isempty(Q.Vertices), continue, end
      [cx,cy] = centroid(Q);
      if ~isfinite(cx) || ~isfinite(cy) || ~isinterior(Q,cx,cy)
        [cx,cy] = project_inside_polyshape(Q,[cx,cy]);
      end
      newSites(p,:) = (1-omega)*sites(p,:) + omega*[cx,cy];
      areas(p) = area(Q);
    end
    sites = newSites;
    cv = std(areas(~isFixed),'omitnan')/mean(areas(~isFixed),'omitnan');
    fprintf('Lloyd it %2d: cv=%.3f (interior only)\n', it, cv);
  end

  % show Lloyd Voronoi for sanity
  figure('Color','w'); hold on
  plot(bx,by,'k-','LineWidth',1.5)
  [Vv,Cv] = voronoin([sites; far]);
  for p=1:size(sites,1)
    id=Cv{p}; if isempty(id)||any(id==1)||any(isinf(id)), continue, end
    Pcell=polyshape(Vv(id,1),Vv(id,2),'Simplify',true);
    Q=intersect(Pcell,swan);
    if isempty(Q.Vertices), continue, end
    plot(Q,'FaceColor',[0.9 0.95 1],'FaceAlpha',0.35,'EdgeColor',[0.2 0.4 0.9])
  end
  scatter(sites(~isFixed,1), sites(~isFixed,2), 6, 'k', 'filled')
  scatter(sites(isFixed,1),  sites(isFixed,2), 12, [0.85 0 0], 'filled')
  axis equal tight; box on
  title('Voronoi with fixed boundary sites'); xlabel x; ylabel y;
else
  sites = centers;
  isFixed = false(size(sites,1),1);
end

%% ===== Approximate unstructured mesh with a structured nx-by-ny grid =====
% Harmonic (Tutte) parameterization (boundary -> rectangle), then sample & invert

% Boundary polyline (drop closing duplicate if present)
if (bx(1)==bx(end)) && (by(1)==by(end))
    B = [bx(1:end-1) by(1:end-1)];
else
    B = [bx by];
end

% Build CDT of sites + boundary, keep triangles inside the swan
m  = size(sites,1);
nb = size(B,1);
Pts = [sites; B];

C = [(m+1:m+nb-1).' (m+2:m+nb).'];   % boundary edge constraints
C = [C; m+nb, m+1];

DT = delaunayTriangulation(Pts, C);
T  = DT.ConnectivityList;

ctr = (Pts(T(:,1),:) + Pts(T(:,2),:) + Pts(T(:,3),:))/3;
inside = isinterior(swan, ctr(:,1), ctr(:,2));
T  = T(inside,:);

% ---- COMPRESS to vertices actually used by T -----------------------------
used = false(size(Pts,1),1);
used(unique(T(:))) = true;                % vertices referenced by T
map  = zeros(size(Pts,1),1);
map(used) = 1:nnz(used);

PtsU = Pts(used,:);                       % compressed coordinates
T     = map(T);                           % reindex connectivity

% Used boundary vertices and their rectangle UV targets
bndIdxAll = (m+1):(m+nb);                 % indices (in Pts) of boundary polyline
bndKeep   = used(bndIdxAll);              % which boundary verts are used by T
bndIdxU   = map(bndIdxAll(bndKeep));      % their indices in the compressed set

% Map used boundary to unit square (choose 4 corners at quarter-perimeter)
nbUsed = nnz(bndKeep);
idxCorners = [1, round(1+nbUsed/4), round(1+nbUsed/2), round(1+3*nbUsed/4)];
walk = @(a,b,n) 1 + mod((a-1) : (a + mod(b-a+n,n) - 1), n);

uB = nan(nbUsed,1); vB = nan(nbUsed,1);
% south: (0,0)->(1,0)
k = walk(idxCorners(1), idxCorners(2), nbUsed); t = linspace(0,1,numel(k))';
uB(k) = t; vB(k) = 0;
% east: (1,0)->(1,1)
k = walk(idxCorners(2), idxCorners(3), nbUsed); t = linspace(0,1,numel(k))';
uB(k) = 1; vB(k) = t;
% north: (1,1)->(0,1)
k = walk(idxCorners(3), idxCorners(4), nbUsed); t = linspace(0,1,numel(k))';
uB(k) = 1 - t; vB(k) = 1;
% west: (0,1)->(0,0)
k = walk(idxCorners(4), idxCorners(1), nbUsed); t = linspace(0,1,numel(k))';
uB(k) = 0; vB(k) = 1 - t;

% Assemble U,V with Dirichlet data on used boundary vertices
nU = size(PtsU,1);
U  = nan(nU,1);  V = nan(nU,1);
U(bndIdxU) = uB; V(bndIdxU) = vB;

% ---- Uniform graph Laplacian on the used vertex graph --------------------
E = [T(:,[1 2]); T(:,[2 3]); T(:,[3 1])];
E = sort(E,2); E = unique(E,'rows');
A = sparse(E(:,1), E(:,2), 1, nU, nU); A = A + A.';   % symmetric adjacency
deg = sum(A,2);
L   = spdiags(deg,0,nU,nU) - A;

free = isnan(U);  bdry = ~free;
if any(free)
    Aint = L(free,free);  bU = -L(free,bdry) * U(bdry);
    if rcond(full(Aint)) < 1e-12
        U(free) = pinv(full(Aint)) * bU;
    else
        U(free) = Aint \ bU;
    end
end
free = isnan(V);  bdry = ~free;
if any(free)
    Aint = L(free,free);  bV = -L(free,bdry) * V(bdry);
    if rcond(full(Aint)) < 1e-12
        V(free) = pinv(full(Aint)) * bV;
    else
        V(free) = Aint \ bV;
    end
end

% Triangulation in UV space (all finite now)
triUV = triangulation(T, [U V]);

% ---- Sample a structured (nx+1)×(ny+1) UV grid and invert to XY ---------
xi  = linspace(0,1,nx+1);
eta = linspace(0,1,ny+1);
[Xi,Eta] = meshgrid(xi,eta);
UVq = [Xi(:) Eta(:)];

ti = pointLocation(triUV, UVq);
Xs = nan(numel(UVq(:,1)),1); Ys = Xs;
validUV = ~isnan(ti);
bc  = cartesianToBarycentric(triUV, ti(validUV), UVq(validUV,:));
TT  = T(ti(validUV),:);
XY0 = PtsU(TT(:,1),:); XY1 = PtsU(TT(:,2),:); XY2 = PtsU(TT(:,3),:);
XY  = bc(:,1).*XY0 + bc(:,2).*XY1 + bc(:,3).*XY2;
Xs(validUV) = XY(:,1); Ys(validUV) = XY(:,2);
Xs = reshape(Xs, ny+1, nx+1);
Ys = reshape(Ys, ny+1, nx+1);

%% ---- Colored structured grid on black background (NaN-robust) ----------
xmin = min(Xs,[],'all','omitnan');  xmax = max(Xs,[],'all','omitnan');
ymin = min(Ys,[],'all','omitnan');  ymax = max(Ys,[],'all','omitnan');
nxp = size(Xs,2); nyp = size(Xs,1);

nxden = max(xmax - xmin, eps);
nyden = max(ymax - ymin, eps);
normx = @(x) (x - xmin)/nxden;         % map to [0,1]
normy = @(y) (y - ymin)/nyden;
clamp = @(t) max(0,min(1,t));
hwrap = @(h) h - floor(h);             % wrap hue to [0,1)

figure('Color','k'); ax = gca; hold(ax,'on');
set(ax,'Color','k','XColor',[.9 .9 .9],'YColor',[.9 .9 .9]);

% boundary in gray if present
if exist('bx','var') && exist('by','var') && ~isempty(bx) && ~isempty(by)
  plot(bx,by,'Color',[.6 .6 .6],'LineWidth',1.25)
end

% vertical lines (i fixed): hue from x, brightness from y
for i = 1:nxp
  xline = Xs(:,i);  yline = Ys(:,i);
  validLine = ~isnan(xline) & ~isnan(yline);
  if nnz(validLine) < 2, continue; end
  h = hwrap( normx( mean(xline(validLine),'omitnan') ) );
  v = clamp( 0.35 + 0.65*normy( mean(yline(validLine),'omitnan') ) );
  s = 0.95;
  rgb = hsv2rgb([h s v]); rgb = rgb(1,:);
  plot(xline, yline, '-', 'LineWidth', 1.15, 'Color', rgb);
end

% horizontal lines (j fixed): hue from y, brightness from x
for j = 1:nyp
  xline = Xs(j,:);  yline = Ys(j,:);
  validLine = ~isnan(xline) & ~isnan(yline);
  if nnz(validLine) < 2, continue; end
  h = hwrap( normy( mean(yline(validLine),'omitnan') ) );
  v = clamp( 0.35 + 0.65*normx( mean(xline(validLine),'omitnan') ) );
  s = 0.95;
  rgb = hsv2rgb([h s v]); rgb = rgb(1,:);
  plot(xline, yline, '-', 'LineWidth', 1.15, 'Color', rgb);
end

axis equal tight; box on
title('Structured grid with HSV by (x,y)','Color',[.95 .95 .95])
xlabel('x','Color',[.9 .9 .9]); ylabel('y','Color',[.9 .9 .9]);

%% ========================= local funcs =========================
function p = bottom(s), p=[s,0]; end
function p = right(s),  p=[1+2*s-2*s.^2, s]; end
function p = top(s),    p=[s, 1-3*s+3*s.^2]; end
function p = left(s),   p=[0,s]; end

function x = coons(u,v)
  Bu=bottom(u); Tu=top(u); Lv=left(v); Rv=right(v);
  P00=bottom(0); P10=bottom(1); P01=left(1); P11=right(1);
  bilinear=(1-u)*(1-v)*P00 + u*(1-v)*P10 + (1-u)*v*P01 + u*v*P11;
  x=(1-v)*Bu + v*Tu + (1-u)*Lv + u*Rv - bilinear;
end

function [bx,by]=boundaryPolygon(nPerEdge)
  s=linspace(0,1,nPerEdge+1); s(end)=[];
  Bcells=arrayfun(@(t) bottom(t), s, 'UniformOutput', false); B=vertcat(Bcells{:});
  Rcells=arrayfun(@(t) right(t),  s, 'UniformOutput', false); R=vertcat(Rcells{:});
  Tcells=arrayfun(@(t) top(t),    fliplr(s), 'UniformOutput', false); T=vertcat(Tcells{:});
  Lcells=arrayfun(@(t) left(t),   fliplr(s), 'UniformOutput', false); L=vertcat(Lcells{:});
  P=[B;R;T;L;bottom(0)];
  bx=P(:,1); by=P(:,2);
end

function A=polyAreaSigned(P)
  x=P(:,1); y=P(:,2);
  A=0.5*sum(x.*y([2:end 1]) - y.*x([2:end 1]));
end

function val=selfIntersectQuad(P,EPS)
  val = segmentsIntersectStrict(P(1,:),P(2,:),P(3,:),P(4,:),EPS) || ...
        segmentsIntersectStrict(P(2,:),P(3,:),P(4,:),P(1,:),EPS);
end

function B=bbox(P), B=[min(P(:,1)) min(P(:,2)) max(P(:,1)) max(P(:,2))]; end
function tf=bboxesOverlap(B1,B2,EPS)
  tf = ~(B1(3)<B2(1)-EPS || B2(3)<B1(1)-EPS || B1(4)<B2(2)-EPS || B2(4)<B1(2)-EPS);
end
function tf=areNeighbors(i1,j1,i2,j2)
  di=abs(i1-i2); dj=abs(j1-j2); tf=(di+dj==1) || (di==0 && dj==0);
end
function tf=polysIntersectStrict(P,Q,EPS)
  if ~bboxesOverlap(bbox(P),bbox(Q),EPS), tf=false; return; end
  for i=1:4, p1=P(i,:); p2=P(mod(i,4)+1,:);
    for j=1:4, q1=Q(j,:); q2=Q(mod(j,4)+1,:);
      if segmentsIntersectStrict(p1,p2,q1,q2,EPS), tf=true; return; end
    end
  end
  [inP,onP]=inpolygon(P(1,1),P(1,2),Q(:,1),Q(:,2));
  [inQ,onQ]=inpolygon(Q(1,1),Q(1,2),P(:,1),P(:,2));
  tf=(inP && ~onP) || (inQ && ~onQ);
end
function tf=segmentsIntersectStrict(a,b,c,d,EPS)
  o1=orient(a,b,c); o2=orient(a,b,d); o3=orient(c,d,a); o4=orient(c,d,b);
  tf=(o1*o2<-EPS) && (o3*o4<-EPS);
end
function s=orient(a,b,c)
  s=(b(1)-a(1))*(c(2)-a(2)) - (b(2)-a(2))*(c(1)-a(1));
end

% --- boundary extraction + colinear pruning
function [northZ,southZ,eastZ,westZ]=extract_raw_boundary_lists(mask,u,v,edgeTol)
  [ny,nx]=size(mask);
  yTop=zeros(1,nx);
  for i=1:nx
    col=mask(:,i); j=find(col,1,'last');
    if isempty(j), yTop(i)=v(1); else, yTop(i)=v(j+1); end
  end
  northZ = build_stair_poly_xfirst(u,yTop,edgeTol);
  yBot=zeros(1,nx);
  for i=1:nx
    col=mask(:,i); j=find(col,1,'first');
    if isempty(j), yBot(i)=v(1); else, yBot(i)=v(j); end
  end
  southZ = build_stair_poly_xfirst(u,yBot,edgeTol);
  xRight=zeros(1,ny);
  for j=1:ny
    row=mask(j,:); i=find(row,1,'last');
    if isempty(i), xRight(j)=u(1); else, xRight(j)=u(i+1); end
  end
  eastZ = build_stair_poly_yfirst(v,xRight,edgeTol);
  xLeft=zeros(1,ny);
  for j=1:ny
    row=mask(j,:); i=find(row,1,'first');
    if isempty(i), xLeft(j)=u(1); else, xLeft(j)=u(i); end
  end
  westZ = build_stair_poly_yfirst(v,xLeft,edgeTol);
end

function z=build_stair_poly_xfirst(u,y,edgeTol)
  N=numel(y); verts=zeros(0,2); ycur=y(1);
  verts(end+1,:)= [u(1), ycur];
  for i=1:N-1
    ynext=y(i+1); unext=u(i+1);
    if abs(ynext - ycur) < eps
      if abs(ycur-0)<=edgeTol || abs(ycur-1)<=edgeTol
        verts(end+1,:)=[unext, ycur];
      end
    else
      verts(end+1,:)=[unext, ycur];
      verts(end+1,:)=[unext, ynext];
      ycur=ynext;
    end
  end
  if ~(abs(ycur-0)<=edgeTol || abs(ycur-1)<=edgeTol)
    if ~(verts(end,1)==u(end) && verts(end,2)==ycur)
      verts(end+1,:)=[u(end), ycur];
    end
  end
  if ~isempty(verts), d=[true; any(diff(verts,1,1),2)]; verts=verts(d,:); end
  z=complex(verts(:,1), verts(:,2));
end

function z=build_stair_poly_yfirst(v,x,edgeTol)
  N=numel(x); verts=zeros(0,2); xcur=x(1);
  verts(end+1,:)=[xcur, v(1)];
  for j=1:N-1
    xnext=x(j+1); vnext=v(j+1);
    if abs(xnext - xcur) < eps
      if abs(xcur-0)<=edgeTol || abs(xcur-1)<=edgeTol
        verts(end+1,:)=[xcur, vnext];
      end
    else
      verts(end+1,:)=[xcur, vnext];
      verts(end+1,:)=[xnext, vnext];
      xcur=xnext;
    end
  end
  if ~(abs(xcur-0)<=edgeTol || abs(xcur-1)<=edgeTol)
    if ~(verts(end,1)==xcur && verts(end,2)==v(end))
      verts(end+1,:)=[xcur, v(end)];
    end
  end
  if ~isempty(verts), d=[true; any(diff(verts,1,1),2)]; verts=verts(d,:); end
  z=complex(verts(:,1), verts(:,2));
end

function zOut = smooth_turns(zIn, imagTol, edgeTol)
  z = zIn;
  for k = 2:numel(z)-1
    xk = real(z(k)); yk = imag(z(k));
    if abs(xk-0)<=edgeTol || abs(xk-1)<=edgeTol || abs(yk-0)<=edgeTol || abs(yk-1)<=edgeTol
      continue
    end
    a = z(k)-z(k-1); b = z(k+1)-z(k);
    if abs(a)==0 || abs(b)==0, continue; end
    r = b/a;
    if abs(real(r)) <= imagTol && abs(imag(r)) > 0
      z(k) = 0.5*(z(k-1)+z(k+1));
    end
  end
  zOut = z;
end

function [U,V] = coons_uv_from_boundary(Xi,Eta, southZ, eastZ, northZ, westZ)
  B = pl_eval(southZ, Xi); T = pl_eval(northZ, Xi);
  L = pl_eval(westZ, Eta); R = pl_eval(eastZ, Eta);
  P00 = pl_eval(southZ, 0); P10 = pl_eval(southZ, 1);
  P01 = pl_eval(westZ, 1); P11 = pl_eval(eastZ, 1);
  bilinear=(1-Xi).*(1-Eta).*P00 + Xi.*(1-Eta).*P10 + (1-Xi).*Eta.*P01 + Xi.*Eta.*P11;
  Z = (1-Eta).*B + Eta.*T + (1-Xi).*L + Xi.*R - bilinear;
  U=real(Z); V=imag(Z);
end

function Z = pl_eval(zlist, s)
  z=zlist(:); if numel(z)<2, Z=z*ones(size(s)); return; end
  seg=diff(z); L=[0;cumsum(abs(seg))]; Ltot=L(end);
  if Ltot==0, Z=z(end)*ones(size(s)); return; end
  t=s.*Ltot; Z=zeros(size(s));
  for k=1:numel(t)
    x=t(k);
    if x<=0, Z(k)=z(1); continue; end
    if x>=Ltot, Z(k)=z(end); continue; end
    i=find(L<=x,1,'last'); if i==numel(z), i=i-1; end
    a=(x-L(i))/(L(i+1)-L(i));
    Z(k)=z(i)+a*(z(i+1)-z(i));
  end
end

function [X,Y] = coons_map_array(U,V)
  X=zeros(size(U)); Y=zeros(size(V));
  for k=1:numel(U)
    p=coons(U(k),V(k)); X(k)=p(1); Y(k)=p(2);
  end
end

function plot_gridlines(X,Y)
  [nr,nc]=size(X);
  for j=1:nr, plot(X(j,:),Y(j,:),'-','LineWidth',0.7); end
  for i=1:nc, plot(X(:,i),Y(:,i),'-','LineWidth',0.7); end
end

function p_in = project_inside_polyshape(P, p)
  [vx,vy] = boundary(P); V=[vx vy]; V=V(~any(isnan(V),2),:);
  if isempty(V), [cx,cy]=centroid(P); p_in=[cx,cy]; return, end
  if any(V(1,:)~=V(end,:)), V=[V; V(1,:)]; end
  p = p(:).'; bestd2=inf; p_proj=V(1,:);
  for k=1:size(V,1)-1
    a=V(k,:); b=V(k+1,:); ab=b-a; den=dot(ab,ab);
    if den==0, q=a; else, t=max(0,min(1,dot(p-a,ab)/den)); q=a+t*ab; end
    d2=sum((p-q).^2); if d2<bestd2, bestd2=d2; p_proj=q; end
  end
  [cx,cy]=centroid(P); v=[cx,cy]-p_proj; nv=hypot(v(1),v(2));
  if nv<1e-12, p_in=[cx,cy]; return, end
  v=v/nv; step=max(1e-6,min(1e-2,sqrt(bestd2)+1e-3));
  p_try=p_proj+step*v; if isinterior(P,p_try(1),p_try(2)), p_in=p_try; return, end
  lo=0; hi=1; p_in=p_proj;
  while hi-lo>1e-6
    mid=0.5*(lo+hi); q=p_proj+mid*([cx,cy]-p_proj);
    if isinterior(P,q(1),q(2)), p_in=q; lo=mid; else, hi=mid; end
  end
  if ~isinterior(P,p_in(1),p_in(2)), p_in=[cx,cy]; end
end

% ---- helper: sample points along boundary funcs (works on older MATLAB)
function M = mk_samples(f, ss)
  C = arrayfun(@(t) f(t), ss, 'UniformOutput', false);
  M = vertcat(C{:});                     % N×2
end
