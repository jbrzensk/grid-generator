function swan_Splus_grid_areaheat
% SWAN_SPLUS_GRID_AREAHEAT (toolbox-free)
% 1) Build S_plus = {(xi,eta) in [0,1]^2 : inside Swan AND J>0}
% 2) Generate structured grid inside S_plus (J-safe vertical trimming)
% 3) Map the grid to the Swan
% 4) Plot: (a) grid in S_plus, (b) mapped grid, (c) heatmap on [0,1]^2 of
%    the AREA of the mapped quadrilaterals (Ny-1 by Nx-1), shown over the unit square.

%% ---- knobs ----
NxS = 161;                 % sampling for mask (xi)
NyS = 121;                 % sampling for mask (eta)
Nx  = 81;                  % structured grid columns
Ny  = 61;                  % structured grid rows
poly_pts = 1000;           % Swan polygon density (for inpolygon)
supersample = 2;           % marching-squares upsample (>=1)
simplifyTol = 5e-4;        % boundary simplification tolerance (param space)
J_tau  = 3e-4;             % treat J <= J_tau as unsafe for trimming
Mgeom  = 12*Ny;            % samples/column to find inside runs (geometry)
MJscan = 32*Ny;            % samples/column to trim by J

%% ---- Swan boundary polygon (physical) ----
s = linspace(0,1,poly_pts);
px = [s,  1+2*s-2*s.^2,  fliplr(s),       zeros(size(s))];
py = [0*s,    s,        fliplr(1-3*s+3*s.^2),  fliplr(s)];

f1=@(s)([s,0]);
f2=@(s)([1+2*s-2*s^2,s]);
f3=@(s)([s,1-3*(s)+3*(s)^2]);
f4=@(s)([0,s]);

%% ---- TFI on parameter grid & masks on cell centers ----
[Xi,Et] = meshgrid(linspace(0,1,NxS), linspace(0,1,NyS));
[Xp,Yp] = coons_eval_swan(Xi,Et);

% "inside Swan?" at cell centers (physical)
Xc = 0.25*(Xp(1:end-1,1:end-1)+Xp(1:end-1,2:end)+Xp(2:end,1:end-1)+Xp(2:end,2:end));
Yc = 0.25*(Yp(1:end-1,1:end-1)+Yp(1:end-1,2:end)+Yp(2:end,1:end-1)+Yp(2:end,2:end));
BW_in = inpolygon(Xc, Yc, px, py);

% numeric J on TFI nodes -> average to cell centers for robust masking
dxi  = 1/(NxS-1);  deta = 1/(NyS-1);
[Xxi,Xet] = grads(Xp,dxi,deta);
[Yxi,Yet] = grads(Yp,dxi,deta);
Jnode = Xxi.*Yet - Xet.*Yxi;
Jc = 0.25*(Jnode(1:end-1,1:end-1)+Jnode(1:end-1,2:end)+Jnode(2:end,1:end-1)+Jnode(2:end,2:end));
BW_J = (Jc > J_tau);

% intersect: S_plus mask on cell centers
BW = BW_in & BW_J;
BW = fillholes_nip(BW);
BW = majority3x3_nip(BW);
BW = areaopen_nip(BW,3);
BW = keep_largest_component_nip(BW);
if ~any(BW(:)), error('S ∩ {J>0} empty; reduce J_tau or refine NxS/NyS.'); end

%% ---- boundary of S_plus via padded marching-squares ----
BWpad = false(size(BW,1)+2, size(BW,2)+2);
BWpad(2:end-1, 2:end-1) = BW;
if supersample>1
    BWc = kron(BWpad, ones(supersample));
    C = contourc(double(BWc), [0.5 0.5]);
    scale = 1/supersample; BWsample=BWc;
else
    C = contourc(double(BWpad), [0.5 0.5]);
    scale = 1;             BWsample=BWpad;
end
off = 1.5;
loops = contour_to_loops(C);

keep=false(1,numel(loops)); xi_loops=cell(size(loops)); eta_loops=xi_loops;
for k=1:numel(loops)
    xc = mean(loops{k}(:,1)); yc = mean(loops{k}(:,2));
    r = max(1,min(size(BWsample,1),round(yc)));
    c = max(1,min(size(BWsample,2),round(xc)));
    if BWsample(r,c)~=1, continue; end       % keep only BW==1 loops
    xidx = loops{k}(:,1)*scale; yidx = loops{k}(:,2)*scale;
    xi  = (xidx - off)/(NxS-1); eta = (yidx - off)/(NyS-1);
    xi  = min(max(xi,0),1);    eta = min(max(eta,0),1);
    P = rdp_simplify([xi(:) eta(:)], simplifyTol);
    if signed_area(P(:,1),P(:,2))<0, P=flipud(P); end
    if any(P(1,:)~=P(end,:)), P(end+1,:)=P(1,:); end
    xi_loops{k}=P(:,1); eta_loops{k}=P(:,2); keep(k)=true;
end
kept=find(keep);
if isempty(kept), error('Could not extract boundary of S_plus.'); end
if numel(kept)==1
    xi_b=xi_loops{kept}; eta_b=eta_loops{kept};
else
    A = cellfun(@(x,y) abs(signed_area(x,y)), xi_loops(kept), eta_loops(kept));
    [~,ix] = max(A); xi_b=xi_loops{kept(ix)}; eta_b=eta_loops{kept(ix)};
end

%% ---- structured grid inside S_plus (J-safe per column) ----
[Xi_keep, Et_keep] = grid_in_polygon_vertical_Jsafe( ...
    xi_b, eta_b, Nx, Ny, Mgeom, MJscan, @coonsJ_swan, J_tau);

%% ---- map to Swan ----
[Xk, Yk] = coons_eval_swan(Xi_keep, Et_keep);
[Xlap,Ylap]=laplaceGridFromReference(Xk,Yk,f1,f2,f3,f4);


%% ---- AREA heatmap on full [0,1]^2 (no normalization) -------------------
% Build a uniform Ny-by-Nx node grid over [0,1]^2 for the heatmap.
xiU  = linspace(0,1,Nx);
etaU = linspace(0,1,Ny);
[XiU, EtU] = meshgrid(xiU, etaU);
[XU, YU]   = coons_eval_swan(XiU, EtU);

% Compute physical area for each param cell; store as Ny-1 by Nx-1
Area = zeros(Ny-1, Nx-1);
for j=1:Ny-1
  for i=1:Nx-1
    xq = [XU(j,i), XU(j,i+1), XU(j+1,i+1), XU(j+1,i)];
    yq = [YU(j,i), YU(j,i+1), YU(j+1,i+1), YU(j+1,i)];
    Area(j,i) = 0.5*abs( sum(xq.*circshift(yq,[0 -1])) - sum(yq.*circshift(xq,[0 -1])) );
  end
end
Area=Area*(Nx-1)*(Ny-1)
% Plot on cell centers for a true heatmap over [0,1]^2
xiC  = 0.5*(xiU(1:end-1)+xiU(2:end));
etaC = 0.5*(etaU(1:end-1)+etaU(2:end));
[XiC, EtC] = meshgrid(xiC, etaC);

%% ---- plots: exactly 3 subplots ----
figure('Color','w'); tiledlayout(1,3,'Padding','compact','TileSpacing','compact');

% (1) structured grid in S^+
nexttile; hold on; axis equal; box on;
title('Grid in S^+'); xlabel('\xi'); ylabel('\eta');
plot([0 1 1 0 0],[0 0 1 1 0],'k-');
plot(xi_b,eta_b,'r-','LineWidth',1.6);
for j=1:Ny, plot(Xi_keep(j,:),Et_keep(j,:),'b-'); end
for i=1:Nx, plot(Xi_keep(:,i),Et_keep(:,i),'b-'); end
legend({'[0,1]^2','\partial S^+','grid'},'Location','best');

% (2) mapped grid in Swan
nexttile; hold on; axis equal; box on;
title('Mapped grid (Swan)'); xlabel('x'); ylabel('y');
plot(px,py,'k-','LineWidth',1.2);
for j=1:Ny, plot(Xlap(j,:),Ylap(j,:),'k-'); end
for i=1:Nx, plot(Xlap(:,i),Ylap(:,i),'k-'); end
legend({'boundary','mapped grid'},'Location','best');

% (3) heatmap of mapped quad AREAS on [0,1]^2 (cell centers)
nexttile; hold on; axis equal tight; box on;
title('Mapped quad area on [0,1]^2'); xlabel('\xi'); ylabel('\eta');
surf(XiC, EtC, Area, 'EdgeColor','none'); view(2);
colormap(gca, parula);
amin = min(Area(:)); amax = max(Area(:));
if ~isfinite(amin) || ~isfinite(amax) || amax<=amin
    amin = 0; amax = 1;
end
caxis([amin amax]);
cb = colorbar; cb.Label.String = 'area(F(cell))';


% (4) laplace grid




end

%% =================== Swan Coons map & Jacobian ==========================
function [X,Y] = coons_eval_swan(XI,ETA)
% Coons patch for Swan edges:
%   x = xi*(1 + 2*eta - 2*eta^2),  y = eta*(1 - 3*xi + 3*xi^2)
X = XI .* (1 + 2*ETA - 2*ETA.^2);
Y = ETA .* (1 - 3*XI  + 3*XI.^2);
end

function J = coonsJ_swan(XI,ETA)
xxi = 1 + 2*ETA - 2*ETA.^2;
xet = XI .* (2 - 4*ETA);
yxi = ETA .* (6*XI - 3);
yet = 1 - 3*XI + 3*XI.^2;
J = xxi.*yet - xet.*yxi;
end

%% =================== Column builder with J-safe trim ====================
function [XI, ETA] = grid_in_polygon_vertical_Jsafe(xp,yp,Nx,Ny,Mgeom,MJ,Jfun,Jtau)
% Build vertical columns inside polygon; choose inside run (prefer eta=0),
% then trim to longest contiguous sub-run where J>Jtau; place Ny points.
if xp(1)~=xp(end) || yp(1)~=yp(end), xp(end+1)=xp(1); yp(end+1)=yp(1); end
xi_cols = linspace(0,1,Nx);
XI  = repmat(xi_cols,Ny,1);
ETA = nan(Ny,Nx);
eta_s = linspace(0,1,Mgeom).'; tol=1e-12; warned=false;

for i=1:Nx
    xi = xi_cols(i);
    in = inpolygon(xi*ones(Mgeom,1), eta_s, xp, yp);
    e  = diff([false; in; false]);
    st = find(e==1); en = find(e==-1)-1;

    if isempty(st)
        if i==1
            ok=false;
            for ii=2:Nx
                in2 = inpolygon(xi_cols(ii)*ones(Mgeom,1), eta_s, xp, yp);
                e2  = diff([false; in2; false]); rs=find(e2==1); re=find(e2==-1)-1;
                if ~isempty(rs)
                    [a,b] = choose_run(rs,re,eta_s,tol);
                    [a,b,okJ] = trim_by_J(xi_cols(ii), a, b, MJ, Jfun, Jtau);
                    if okJ, ETA(:,i)=linspace(a,b,Ny).'; ok=true; break, end
                end
            end
            if ~ok, error('No valid J>Jtau column found; relax J_tau or refine.'); end
        else
            if ~warned, warning('Some columns miss S; reusing previous column.'); warned=true; end
            ETA(:,i) = ETA(:,i-1);
        end
        continue
    end

    [a,b] = choose_run(st,en,eta_s,tol);
    [a,b,okJ] = trim_by_J(xi, a, b, MJ, Jfun, Jtau);
    if ~okJ
        [~,kmax] = max(eta_s(en)-eta_s(st));
        a2 = eta_s(st(kmax)); b2 = eta_s(en(kmax));
        [a2,b2,okJ] = trim_by_J(xi, a2, b2, MJ, Jfun, Jtau);
        if okJ
            ETA(:,i)=linspace(a2,b2,Ny).';
        else
            if i==1, ETA(:,i)=linspace(0,1,Ny).'; else, ETA(:,i)=ETA(:,i-1); end
        end
    else
        ETA(:,i)=linspace(a,b,Ny).';
    end
end
end

function [a,b] = choose_run(st,en,eta_s,tol)
k = find(eta_s(st)<=tol,1);         % prefer touching eta=0
if isempty(k), [~,k]=max(eta_s(en)-eta_s(st)); end
a = eta_s(st(k)); b = eta_s(en(k));
end

function [a_out,b_out,ok] = trim_by_J(xi, a, b, M, Jfun, Jtau)
eta = linspace(a,b,M).';
J   = Jfun(xi*ones(M,1), eta);
keep = (J > Jtau);
ok = any(keep);
if ~ok, a_out=a; b_out=b; return, end
edges = diff([false; keep; false]); st=find(edges==1); en=find(edges==-1)-1;
[~,k] = max(eta(en)-eta(st));
a_out = eta(st(k)); b_out = eta(en(k));
% small safety margin away from J≈0
margin = max(2/M, 3e-3);
a_out = a_out + margin;  b_out = b_out - margin;
if b_out <= a_out, ok=false; end
end

%% =================== Numeric grads on rect grid =========================
function [Xxi,Xet] = grads(X,dxi,deta)
Xxi = zeros(size(X)); Xet = Xxi;
Xxi(:,2:end-1) = (X(:,3:end)-X(:,1:end-2))/(2*dxi);
Xxi(:,1)       = (X(:,2)-X(:,1))/dxi;
Xxi(:,end)     = (X(:,end)-X(:,end-1))/dxi;
Xet(2:end-1,:) = (X(3:end,:)-X(1:end-2,:))/(2*deta);
Xet(1,:)       = (X(2,:)-X(1,:))/deta;
Xet(end,:)     = (X(end,:)-X(end-1,:))/deta;
end

%% =================== Morphology (no toolboxes) ==========================
function BWf = fillholes_nip(BW)
BW = logical(BW); Z = ~BW;
edge=false(size(BW)); edge([1 end],:)=Z([1 end],:); edge(:,[1 end])=edge(:,[1 end])|Z(:,[1 end]);
K=[0 1 0;1 1 1;0 1 0]; changed=true;
while changed
    grow = conv2(double(edge),K,'same')>0 & Z & ~edge;
    changed=any(grow(:)); edge=edge|grow;
end
holes = Z & ~edge; BWf = BW | holes;
end

function BW2 = majority3x3_nip(BW)
S = conv2(double(BW), ones(3), 'same'); BW2=BW; BW2(S>=5)=true; BW2(S<=4)=false;
end

function BW2 = areaopen_nip(BW, minsz)
BW=logical(BW); [nr,nc]=size(BW); vis=false(nr,nc); BW2=BW; nbr=[0 1;0 -1;1 0;-1 0];
for r=1:nr, for c=1:nc
  if BW2(r,c) && ~vis(r,c)
    q=zeros(numel(BW2),2); head=1; tail=1; q(1,:)=[r c]; vis(r,c)=true; comp=[r c];
    while head<=tail
      rr=q(head,1); cc=q(head,2); head=head+1;
      for k=1:4
        r2=rr+nbr(k,1); c2=cc+nbr(k,2);
        if r2>=1 && r2<=nr && c2>=1 && c2<=nc && BW2(r2,c2) && ~vis(r2,c2)
          tail=tail+1; q(tail,:)=[r2 c2]; vis(r2,c2)=true; comp(end+1,:)=[r2 c2]; %#ok<AGROW>
        end
      end
    end
    if size(comp,1)<minsz
      BW2(sub2ind([nr nc], comp(:,1), comp(:,2))) = false;
    end
  end
end, end
end

function BW1 = keep_largest_component_nip(BW)
BW=logical(BW); [nr,nc]=size(BW); vis=false(nr,nc); BW1=false(size(BW));
nbr=[0 1;0 -1;1 0;-1 0]; best=0; bestComp=[];
for r=1:nr, for c=1:nc
  if BW(r,c) && ~vis(r,c)
    q=zeros(numel(BW),2); head=1; tail=1; q(1,:)=[r c]; vis(r,c)=true; comp=[r c];
    while head<=tail
      rr=q(head,1); cc=q(head,2); head=head+1;
      for k=1:4
        r2=rr+nbr(k,1); c2=cc+nbr(k,2);
        if r2>=1 && r2<=nr && c2>=1 && c2<=nc && BW(r2,c2) && ~vis(r2,c2)
          tail=tail+1; q(tail,:)=[r2 c2]; vis(r2,c2)=true; comp(end+1,:)=[r2 c2]; %#ok<AGROW>
        end
      end
    end
    if size(comp,1)>best, best=size(comp,1); bestComp=comp; end
  end
end, end
if ~isempty(bestComp)
    BW1(sub2ind([nr nc],bestComp(:,1),bestComp(:,2))) = true;
end
end

function loops = contour_to_loops(C)
i=1; loops={};
while i < size(C,2)
    n = C(2,i); loops{end+1} = C(:, i+1:i+n).'; %#ok<AGROW>
    i = i + n + 1;
end
end

function A = signed_area(x,y)
A = 0.5*sum( x(:).*circshift(y(:),-1) - y(:).*circshift(x(:),-1) );
end

function P = rdp_simplify(P, eps)
if size(P,1)<=2, return; end
dmax=0; idx=0; A=P(1,:); B=P(end,:); AB=B-A; AB2=sum(AB.^2);
for i=2:size(P,1)-1
  AP=P(i,:)-A; t=max(0,min(1,(AP*AB')/AB2)); proj=A+t*AB;
  d=norm(P(i,:)-proj); if d>dmax, dmax=d; idx=i; end
end
if dmax>eps
  P1=rdp_simplify(P(1:idx,:),eps); P2=rdp_simplify(P(idx:end,:),eps);
  P=[P1(1:end-1,:); P2];
end
end

function [x_lap, y_lap] = laplaceGridFromReference(xref, yref, s1, s2, s3, s4)
%LAPLACEGRIDFROMREFERENCE  Weighted-length (Eq. 12) Laplace grid from a reference grid.
%   [X_LAP, Y_LAP] = LAPLACEGRIDFROMREFERENCE(XREF,YREF,S1,S2,S3,S4)
%   takes a reference grid (XREF,YREF) of size I-by-J and four boundary
%   curves S1..S4: [x,y] = S_k(t), t in [0,1], and returns the interior
%   grid (X_LAP,Y_LAP) that minimizes the reference-weighted length
%   functional (Eq. 12) subject to Dirichlet boundary data.
%
%   Reference weights use the reference edge lengths:
%     l_h(i,j) = ||(xref(i+1,j),yref(i+1,j)) - (xref(i,j),yref(i,j))||
%     l_v(i,j) = ||(xref(i,j+1),yref(i,j+1)) - (xref(i,j),yref(i,j))||
%   Each interior node solves a weighted 5-point Laplacian:
%     (wE+wW+wN+wS) u_ij
%       - wE u_{i+1,j} - wW u_{i-1,j} - wN u_{i,j+1} - wS u_{i,j-1} = rhs
%   with Dirichlet contributions folded into rhs. The same SPD matrix is
%   solved for X and for Y.
%
%   Boundary convention (counter-clockwise):
%     s1: bottom  (left->right), s2: right (bottom->top),
%     s3: top     (right->left), s4: left  (top->bottom).

  [I,J] = size(xref);
  assert(isequal(size(yref),[I,J]), 'xref and yref must be same size.');
  assert(I>=2 && J>=2, 'Grid must be at least 2-by-2.');

  % --- Build and fill boundary from curves (counter-clockwise) ---
  X = nan(I,J); Y = nan(I,J);

  % helper to safely evaluate [x,y]=curve(t) for vector t
  function [xs,ys] = eval_curve(curve, tt)
    xs = zeros(numel(tt),1); ys = zeros(numel(tt),1);
    for k = 1:numel(tt)
      p = curve(tt(k));
      xs(k) = p(1); ys(k) = p(2);
    end
  end

  tb = linspace(0,1,I);       % bottom param (i = 1..I)
  tr = linspace(0,1,J);       % right  param (j = 1..J)
  tt = linspace(0,1,I);       % top    param (i = 1..I)  (right->left)
  tl = linspace(0,1,J);       % left   param (j = 1..J)  (top->bottom)

  [xb,yb] = eval_curve(s1, tb);           % bottom
  [xr,yr] = eval_curve(s2, tr);           % right
  [xt,yt] = eval_curve(s3, tt);           % top (right->left along i=1..I)
  [xl,yl] = eval_curve(s4, tl);           % left (top->bottom along j=1..J)

  % bottom j=1, left->right
  X(:,1) = xb(:);         Y(:,1) = yb(:);
  % right  i=I, bottom->top
  X(I,:) = xr(:).';       Y(I,:) = yr(:).';
  % top    j=J, right->left -> place xt(1) at i=1, xt(end) at i=I
  X(:,J) = xt(:);         Y(:,J) = yt(:);
  % left   i=1, top->bottom
  X(1,:) = xl(:).';       Y(1,:) = yl(:).';

  % Optional: check/average corner consistency (simple average if mismatch)
  corners = [X(1,1)    Y(1,1);    X(I,1)   Y(I,1); ...
             X(I,J)    Y(I,J);    X(1,J)   Y(1,J)];
  if any(any(isnan(corners)))
      error('Boundary curves must define all four sides (including corners).');
  end
  % (If you prefer strict consistency, replace with assertions.)

  % --- Reference weights from (xref,yref) ---
  % Horizontal edge lengths l_h: size (I-1) x J
  lh = sqrt( (xref(2:I,:) - xref(1:I-1,:)).^2 + (yref(2:I,:) - yref(1:I-1,:)).^2 );
  % Vertical edge lengths l_v: size I x (J-1)
  lv = sqrt( (xref(:,2:J) - xref(:,1:J-1)).^2 + (yref(:,2:J) - yref(:,1:J-1)).^2 );

  epsw = 1e-14;                        % protect against zero-length refs
  wh = 1 ./ max(lh.^2, epsw);          % horizontal weights = 1 / (l_h)^2
  wv = 1 ./ max(lv.^2, epsw);          % vertical   weights = 1 / (l_v)^2

  % --- Assemble SPD matrix A for interior nodes ---
  ni = I-2; nj = J-2;
  if ni<=0 || nj<=0
    % No interior: just return boundary
    x_lap = X; y_lap = Y; return;
  end
  N = ni*nj;

  % Map (i,j) with i=2..I-1, j=2..J-1 to 1..N
  idx = @(ii,jj) (jj-2)*ni + (ii-1);

  rows = zeros(5*N,1); cols = rows; vals = rows; nnzCtr = 0;
  bx = zeros(N,1); by = zeros(N,1);

  for j = 2:J-1
    for i = 2:I-1
      k = idx(i,j);
      wE = wh(i  ,j  );       % edge (i,j) -> (i+1,j)
      wW = wh(i-1,j  );       % edge (i-1,j) -> (i,j)
      wN = wv(i  ,j  );       % edge (i,j) -> (i,j+1)
      wS = wv(i  ,j-1);       % edge (i,j-1) -> (i,j)
      d  = wE + wW + wN + wS; % diagonal

      % Diagonal
      nnzCtr = nnzCtr+1; rows(nnzCtr)=k; cols(nnzCtr)=k; vals(nnzCtr)=d;

      % East neighbor
      if i+1 <= I-1
        kk = idx(i+1,j); nnzCtr=nnzCtr+1; rows(nnzCtr)=k; cols(nnzCtr)=kk; vals(nnzCtr) = -wE;
      else
        bx(k) = bx(k) + wE*X(i+1,j); by(k) = by(k) + wE*Y(i+1,j);
      end
      % West neighbor
      if i-1 >= 2
        kk = idx(i-1,j); nnzCtr=nnzCtr+1; rows(nnzCtr)=k; cols(nnzCtr)=kk; vals(nnzCtr) = -wW;
      else
        bx(k) = bx(k) + wW*X(i-1,j); by(k) = by(k) + wW*Y(i-1,j);
      end
      % North neighbor
      if j+1 <= J-1
        kk = idx(i,j+1); nnzCtr=nnzCtr+1; rows(nnzCtr)=k; cols(nnzCtr)=kk; vals(nnzCtr) = -wN;
      else
        bx(k) = bx(k) + wN*X(i,j+1); by(k) = by(k) + wN*Y(i,j+1);
      end
      % South neighbor
      if j-1 >= 2
        kk = idx(i,j-1); nnzCtr=nnzCtr+1; rows(nnzCtr)=k; cols(nnzCtr)=kk; vals(nnzCtr) = -wS;
      else
        bx(k) = bx(k) + wS*X(i,j-1); by(k) = by(k) + wS*Y(i,j-1);
      end
    end
  end

  % Trim and build sparse A
  rows = rows(1:nnzCtr); cols = cols(1:nnzCtr); vals = vals(1:nnzCtr);
  A = sparse(rows, cols, vals, N, N);

  % --- Solve A * u = b for x and y interiors ---
  ux = A \ bx; uy = A \ by;

  % --- Write back interior and return ---
  X(2:I-1, 2:J-1) = reshape(ux, [ni, nj]);
  Y(2:I-1, 2:J-1) = reshape(uy, [ni, nj]);

  x_lap = X; y_lap = Y;
end

