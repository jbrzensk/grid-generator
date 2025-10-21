%% cvt_bfgs_warmstart_nostats.m 
% CVT via good seeding (DIY Halton), Lloyd warm-start, then BFGS.
% No toolboxes required.

clear; clc; close all
 
%% ----------------------- knobs -----------------------
Nint = 400;                 % number of generators
useSwanDomain = true;       % false -> unit square
nEdge = 200;                % Swan polygon sampling

warmIters = 3;             % Lloyd warm-start iterations
omega     = 1.0;            % 1 = pure Lloyd step

% BFGS params
maxIter   = 150;            % BFGS iterations
c1        = 1e-4;           % Armijo parameter
beta      = 0.5;            % backtracking factor
t0        = 1.0;            % initial trial step
gtol      = 1e-10;          % stop if ||g||/sqrt(2N) <= gtol
stol      = 1e-12;          % stop if relative step <= stol
resetIfBad= true;           % reset H if curvature test fails

rng(1);

%% ----------------------- domain -----------------------
if useSwanDomain
  [bx,by] = boundaryPolygon(nEdge);
  domainPoly = polyshape(bx,by);
  xmin=min(bx); xmax=max(bx); ymin=min(by); ymax=max(by);
else
  xmin=0; xmax=1; ymin=0; ymax=1;
  domainPoly = polyshape([xmin xmax xmax xmin],[ymin ymin ymax ymax]);
end
pad  = 2*max(xmax-xmin, ymax-ymin);
far  = [xmin-pad ymin-pad; xmax+pad ymin-pad; xmax+pad ymax+pad; xmin-pad ymax+pad]; % bounds Voronoi

%% --------- uniform (DIY Halton) seeding inside domain ----------
sites = halton_in_polyshape_nostats(Nint, domainPoly, [xmin xmax ymin ymax]);

%% --------- Lloyd warm-start (monotone, robust) -------------
for k=1:warmIters
  cent = clipped_centroids(sites, domainPoly, far);
  sites = (1-omega)*sites + omega*cent;
end

%% --------------------- BFGS finish ---------------------
N = size(sites,1); D = 2*N; H = eye(D);
E_hist = nan(maxIter+1,1); e_hist = nan(maxIter+1,1); step_hist = nan(maxIter,1);

[cent, areas] = clipped_centroids_areas(sites, domainPoly, far);
[g, E0]       = grad_and_energy(sites, cent, areas);
E_hist(1)=E0; e_hist(1)=rms_site_to_centroid(sites, cent);

for k = 1:maxIter
  xk = sites(:); gk = g(:);
  if norm(gk)/sqrt(D) <= gtol
    E_hist = E_hist(1:k); e_hist = e_hist(1:k); step_hist = step_hist(1:k-1);
    break
  end

  pdir = -H*gk; if gk.'*pdir >= 0, pdir = -gk; end

  Ek = E_hist(k); t=t0; accepted=false;
  while t>1e-16 && ~accepted
    xtrial = xk + t*pdir; sites_t = reshape(xtrial,[],2);
    [cent_t, areas_t] = clipped_centroids_areas(sites_t, domainPoly, far);
    [~, Et] = grad_and_energy(sites_t, cent_t, areas_t);
    if isfinite(Et) && (Et <= Ek + c1*t*(gk.'*pdir))
      accepted = true;
    else
      t = beta*t;
    end
  end
  if ~accepted
    warning('Line search failed; stopping.'); 
    E_hist = E_hist(1:k); e_hist = e_hist(1:k); step_hist = step_hist(1:k-1);
    break
  end

  xkp1 = xk + t*pdir; sites = reshape(xkp1,[],2);
  [cent, areas] = clipped_centroids_areas(sites, domainPoly, far);
  [gp1, Ep1]    = grad_and_energy(sites, cent, areas); gkp1 = gp1(:);

  E_hist(k+1)=Ep1; e_hist(k+1)=rms_site_to_centroid(sites, cent);
  step_hist(k)=norm(xkp1-xk)/max(1e-16,norm(xk));

  % BFGS with Powell damping
  s = xkp1 - xk; y = gkp1 - gk; sy = s.'*y;
  if sy <= 1e-12
    if resetIfBad
      gamma = (y.'*y)/max(sy,1e-12); H = (1/max(gamma,1e-12))*eye(D);
    end
  else
    Hy = H*y; theta = 1;
    if y.'*Hy < 0.2*sy, theta = (0.8*sy)/(sy - y.'*Hy); end
    ytil = theta*y + (1-theta)*Hy; rho = 1/max(s.'*ytil,1e-12); I = eye(D);
    H = (I - rho*(s*ytil.'))*H*(I - rho*(ytil*s.')) + rho*(s*s.');
  end

  if step_hist(k) <= stol
    E_hist = E_hist(1:k+1); e_hist = e_hist(1:k+1); step_hist = step_hist(1:k);
    break
  end
end
iters_done = numel(E_hist)-1;
fprintf('Warm Lloyd %d iters, then BFGS %d iters. Final e=%.3e\n', warmIters, iters_done, e_hist(end));

%% ----------------------- plots -----------------------
figure('Color','w'); 
semilogy(0:iters_done, e_hist,'-o','LineWidth',1.6,'MarkerSize',4); grid on
xlabel('BFGS iteration k'); ylabel('e_k = RMS ||p-c||');
title('BFGS after Lloyd warm-start');

figure('Color','w'); 
plot(domainPoly,'FaceColor',[.97 .97 1],'FaceAlpha',.25,'EdgeColor',[.5 .5 .7]); hold on
[Vv,Cv]=voronoin([sites; far]);
for p=1:size(sites,1)
  id=Cv{p}; if isempty(id)||any(id==1)||any(isinf(id)), continue, end
  Q=intersect(polyshape(Vv(id,1),Vv(id,2),'Simplify',true), domainPoly);
  if isempty(Q.Vertices), continue, end
  plot(Q,'FaceColor',[0.88 0.93 1],'FaceAlpha',0.45,'EdgeColor',[0.25 0.45 0.9])
end
scatter(sites(:,1),sites(:,2),12,'k','filled'); axis equal tight; box on
title('CVT after warm-start + BFGS'); xlabel x; ylabel y;

%% ======================= helpers =======================
function P = halton_in_polyshape_nostats(N, poly, bbox)
  % Low-discrepancy uniform points inside 'poly' via Halton(2,3) + rejection
  lb = [bbox(1) bbox(3)];
  ub = [bbox(2) bbox(4)];

  P = zeros(N,2); i = 1;
  skip = 1024;                 % burn-in
  shift = rand(1,2);           % Cranley–Patterson shift

  batch = 0;
  while i <= N
    batch = batch + 500;
    H = halton23(batch, skip); skip = skip + batch;
    H = mod(bsxfun(@plus, H, shift), 1);          %#ok<*NBRAK>
    XY = lb + H.*(ub - lb);
    in = isinterior(poly, XY(:,1), XY(:,2));
    keep = XY(in,:);
    n = min(size(keep,1), N - i + 1);
    if n>0
      P(i:i+n-1,:) = keep(1:n,:);
      i = i + n;
    end
  end

  % tiny jitter to avoid exact symmetries
  diam = max(ub - lb);
  P = P + 1e-4*diam*randn(size(P));
end

function H = halton23(n, skip)
  % Basic Halton sequence in bases 2 and 3 (no toolboxes)
  H = zeros(n,2);
  for t = 1:n
    k = skip + t;
    H(t,1) = radicalInverse(k, 2);
    H(t,2) = radicalInverse(k, 3);
  end
end

function v = radicalInverse(k, base)
  % van der Corput radical inverse in given base
  v = 0; f = 1/base;
  while k > 0
    v = v + mod(k, base) * f;
    k = floor(k / base);
    f = f / base;
  end
end

function [g, E] = grad_and_energy(sites, cent, areas)
  % gradient and surrogate energy E = sum m ||p-c||^2
  valid = all(isfinite(cent),2) & isfinite(areas);
  diff  = sites(valid,:) - cent(valid,:);
  m     = areas(valid);
  gfull = zeros(size(sites));
  gfull(valid,:) = 2*diff.*m;                % 2 m_i (p_i - c_i)
  g = gfull;
  E = sum(m .* sum(diff.^2,2));
end

function e = rms_site_to_centroid(sites, cent)
  valid = all(isfinite(cent),2);
  dif = cent(valid,:) - sites(valid,:);
  e = sqrt(mean(sum(dif.^2,2)));
end

function cent = clipped_centroids(sites, domainPoly, far)
  [Vv,Cv]=voronoin([sites; far]); N=size(sites,1); cent=nan(N,2);
  for p=1:N
    id=Cv{p}; if isempty(id)||any(id==1)||any(isinf(id)), continue, end
    Q=intersect(polyshape(Vv(id,1),Vv(id,2),'Simplify',true), domainPoly);
    if isempty(Q.Vertices), continue, end
    [cx,cy]=centroid(Q);
    if ~isfinite(cx) || ~isfinite(cy) || ~isinterior(Q,cx,cy)
      [cx,cy]=project_inside_polyshape(Q,[cx,cy]);
    end
    cent(p,:)=[cx,cy];
  end
end

function [cent, areas] = clipped_centroids_areas(sites, domainPoly, far)
  [Vv,Cv]=voronoin([sites; far]); N=size(sites,1); cent=nan(N,2); areas=nan(N,1);
  for p=1:N
    id=Cv{p}; if isempty(id)||any(id==1)||any(isinf(id)), continue, end
    Q=intersect(polyshape(Vv(id,1),Vv(id,2),'Simplify',true), domainPoly);
    if isempty(Q.Vertices), continue, end
    areas(p)=area(Q);                               % polyshape area()
    [cx,cy]=centroid(Q);
    if ~isfinite(cx) || ~isfinite(cy) || ~isinterior(Q,cx,cy)
      [cx,cy]=project_inside_polyshape(Q,[cx,cy]);
    end
    cent(p,:)=[cx,cy];
  end
end

% ---- Swan-ish Coons map for reference (kept from earlier code) ----
function [bx,by]=boundaryPolygon(nPerEdge)
  s=linspace(0,1,nPerEdge+1); s(end)=[];
  Bcells=arrayfun(@(t) bottom(t), s, 'UniformOutput', false); B=vertcat(Bcells{:});
  Rcells=arrayfun(@(t) right(t),  s, 'UniformOutput', false); R=vertcat(Rcells{:});
  Tcells=arrayfun(@(t) top(t),    fliplr(s), 'UniformOutput', false); T=vertcat(Tcells{:});
  Lcells=arrayfun(@(t) left(t),   fliplr(s), 'UniformOutput', false); L=vertcat(Lcells{:});
  P=[B;R;T;L;bottom(0)]; bx=P(:,1); by=P(:,2);
end
function p=bottom(s), p=[s,0]; end
function p=right(s),  p=[1+2*s-2*s.^2, s]; end
function p=top(s),    p=[s, 1-3*s+3*s.^2]; end
function p=left(s),   p=[0,s]; end
function [X,Y]=coons_map_array(U,V)
  X=zeros(size(U)); Y=zeros(size(V));
  for k=1:numel(U), q=coons(U(k),V(k)); X(k)=q(1); Y(k)=q(2); end
end
function x=coons(u,v)
  Bu=bottom(u); Tu=top(u); Lv=left(v); Rv=right(v);
  P00=bottom(0); P10=bottom(1); P01=left(1); P11=right(1);
  bilinear=(1-u)*(1-v)*P00 + u*(1-v)*P10 + (1-u)*v*P01 + u*v*P11;
  x=(1-v)*Bu + v*Tu + (1-u)*Lv + u*Rv - bilinear;
end

% ---- push a point inside a polyshape if numerical drift occurs ----
function p_in = project_inside_polyshape(P, p)
  [vx,vy]=boundary(P); V=[vx vy]; V=V(~any(isnan(V),2),:);
  if isempty(V), [cx,cy]=centroid(P); p_in=[cx,cy]; return, end
  if any(V(1,:)~=V(end,:)), V=[V; V(1,:)]; end
  p=p(:).'; bestd2=inf; p_proj=V(1,:);
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


