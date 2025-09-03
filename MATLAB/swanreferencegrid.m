function swanreferencegrid
% SWANREFERENCEGRID (toolbox-free)
% 1) Build preimage mask (cell-center test)
% 2) Extract boundary of S by padded marching squares (keeps BW==1 loops)
% 3) Build an Nx-by-Ny structured grid INSIDE S using robust inpolygon sampling
% 4) Map it to the Swan with the closed-form Coons patch and plot

%% --- knobs ---
NxS = 161;                 % sampling for mask (xi)
NyS = 121;                 % sampling for mask (eta)
Nx  = 81;                  % structured grid columns inside S
Ny  = 61;                  % structured grid rows inside S
poly_pts = 800;            % Swan polygon density
supersample = 2;           % >=1; upsample mask before contouring (nearest)
simplifyTol = 1e-3;        % RDP polyline simplification tolerance
Mscan = 8*Ny;              % vertical samples per column for inpolygon scan

%% --- Swan polygon (physical, CCW) ---
sp = linspace(0,1,poly_pts);
px = [sp,  1+2*sp-2*sp.^2,  fliplr(sp),       zeros(size(sp))];
py = [0*sp,    sp,         fliplr(1-3*sp+3*sp.^2),  fliplr(sp) ];

%% --- TFI on [0,1]^2 and inside/outside mask on cells ---
[Xi,Et] = meshgrid(linspace(0,1,NxS), linspace(0,1,NyS));
[Xp,Yp] = coons_eval_swan(Xi,Et);

Xc = 0.25*(Xp(1:end-1,1:end-1)+Xp(1:end-1,2:end)+Xp(2:end,1:end-1)+Xp(2:end,2:end));
Yc = 0.25*(Yp(1:end-1,1:end-1)+Yp(1:end-1,2:end)+Yp(2:end,1:end-1)+Yp(2:end,2:end));
BW = inpolygon(Xc, Yc, px, py);            % (NyS-1) x (NxS-1) logical

% light cleanup (no toolboxes)
BW = fillholes_nip(BW);
BW = majority3x3_nip(BW);
BW = areaopen_nip(BW,3);

%% --- Boundary of S via padded marching squares; keep only BW==1 loops ---
BWpad = false(size(BW,1)+2, size(BW,2)+2);
BWpad(2:end-1, 2:end-1) = BW;

if supersample > 1
    BWc = kron(BWpad, ones(supersample));   % nearest upsample
    C = contourc(double(BWc), [0.5 0.5]);
    scale = 1/supersample;  BWsample = BWc;
else
    C = contourc(double(BWpad), [0.5 0.5]);
    scale = 1;               BWsample = BWpad;
end
off = 1.5;                                  % 1-cell pad

loops = contour_to_loops(C);
keep = false(1,numel(loops)); xi_loops = cell(1,numel(loops)); eta_loops = xi_loops;

for k = 1:numel(loops)
    % classify by centroid pixel in the (padded/upsampled) mask
    xc = mean(loops{k}(:,1)); yc = mean(loops{k}(:,2));
    r = max(1, min(size(BWsample,1), round(yc)));
    c = max(1, min(size(BWsample,2), round(xc)));
    if BWsample(r,c) ~= 1, continue; end      % this loop bounds BW==0

    % map to parameter space
    xidx = loops{k}(:,1) * scale;  yidx = loops{k}(:,2) * scale;
    xi   = (xidx - off) / (NxS-1);
    eta  = (yidx - off) / (NyS-1);
    xi   = min(max(xi,0),1);  eta = min(max(eta,0),1);

    % simplify, CCW, close
    P = rdp_simplify([xi(:) eta(:)], simplifyTol);
    if signed_area(P(:,1),P(:,2)) < 0, P = flipud(P); end
    if any(P(1,:) ~= P(end,:)), P(end+1,:) = P(1,:); end
    xi_loops{k}  = P(:,1);  eta_loops{k} = P(:,2);  keep(k)=true;
end

kept = find(keep);
if isempty(kept), error('No BW==1 boundary loops found.'); end
if numel(kept)==1
    xi_b = xi_loops{kept}; eta_b = eta_loops{kept};
else
    A = cellfun(@(x,y) abs(signed_area(x,y)), xi_loops(kept), eta_loops(kept));
    [~,ix] = max(A); xi_b = xi_loops{kept(ix)}; eta_b = eta_loops{kept(ix)};
end

%% --- Structured grid INSIDE S using robust inpolygon sampling (per column)
[Xi_keep, Et_keep] = grid_in_polygon_vertical_sampled(xi_b, eta_b, Nx, Ny, Mscan);

%% --- Map to Swan via Coons ---
[Xk, Yk] = coons_eval_swan(Xi_keep, Et_keep);

%% --- Plots ---
figure('Color','w'); tl = tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

% (a) parameter space
nexttile; hold on; axis equal; box on;
title('Preimage S and structured (\xi,\eta) grid'); xlabel('\xi'); ylabel('\eta');
plot([0 1 1 0 0],[0 0 1 1 0],'k-');
plot(xi_b, eta_b, 'r-', 'LineWidth', 1.6);
for j=1:Ny, plot(Xi_keep(j,:), Et_keep(j,:),'b-'); end
for i=1:Nx, plot(Xi_keep(:,i), Et_keep(:,i),'b-'); end
legend({'[0,1]^2','boundary of S','structured grid'},'Location','best');

% (b) physical space
nexttile; hold on; axis equal; box on;
title('Mapped grid in Swan'); xlabel('x'); ylabel('y');
plot(px,py,'k-','LineWidth',1.2);
for j=1:Ny, plot(Xk(j,:), Yk(j,:),'k-'); end
for i=1:Nx, plot(Xk(:,i), Yk(:,i),'k-'); end
legend({'Swan boundary','mapped grid'},'Location','best');

fprintf('Boundary verts: %d   Grid: %dx%d   Mask coverage: %.1f%%\n', ...
        numel(xi_b), Nx, Ny, 100*mean(BW(:)));
end

%% ===== local helpers (no toolboxes) =====
function [XI, ETA] = grid_in_polygon_vertical_sampled(xp,yp,Nx,Ny,M)
% Structured grid inside polygon (xp,yp) via vertical columns.
% For each xi_i, sample M points eta∈[0,1], keep the LONGEST inside run.
% Prefer a run that touches eta=0 if present.
if xp(1)~=xp(end) || yp(1)~=yp(end), xp(end+1)=xp(1); yp(end+1)=yp(1); end
xi_cols = linspace(0,1,Nx); XI = repmat(xi_cols,Ny,1); ETA = nan(Ny,Nx);
eta_s   = linspace(0,1,M).'; tol = 1e-10; warned=false;

for i=1:Nx
    x0 = xi_cols(i) * ones(M,1);
    in = inpolygon(x0, eta_s, xp, yp);               % inside mask along the column
    % find contiguous runs where in==true
    edges = diff([false; in; false]);
    run_st = find(edges==1);  run_en = find(edges==-1)-1;
    if isempty(run_st)
        % no intersection; reuse previous column if exists
        if i==1
            % look ahead to find first good column
            found=false;
            for ii=2:Nx
                x2 = xi_cols(ii)*ones(M,1);
                in2 = inpolygon(x2, eta_s, xp, yp);
                e2 = diff([false; in2; false]);
                if any(e2==1)
                    rs = find(e2==1); re = find(e2==-1)-1;
                    % prefer a run touching bottom if any
                    spans = eta_s(re) - eta_s(rs);
                    touchB = find(eta_s(rs)<=tol, 1);
                    if ~isempty(touchB), k=touchB; else [~,k]=max(spans); end
                    ETA(:,i) = linspace(eta_s(rs(k)), eta_s(re(k)), Ny).'; found=true; break
                end
            end
            if ~found, error('No vertical inside-run found anywhere; S may be empty.'); end
        else
            if ~warned, warning('Some columns miss S; reusing previous column.'); warned=true; end
            ETA(:,i) = ETA(:,i-1);
        end
    else
        spans = eta_s(run_en) - eta_s(run_st);
        % prefer one that touches bottom, else the longest
        k = find(eta_s(run_st)<=tol, 1);
        if isempty(k), [~,k] = max(spans); end
        ETA(:,i) = linspace(eta_s(run_st(k)), eta_s(run_en(k)), Ny).';
    end
end
end

function [X,Y] = coons_eval_swan(XI,ETA)
% Closed-form Coons for these edges:
%   x(ξ,η) = ξ * (1 + 2η - 2η^2)
%   y(ξ,η) = η * (1 - 3ξ + 3ξ^2)
X = XI .* (1 + 2*ETA - 2*ETA.^2);
Y = ETA .* (1 - 3*XI  + 3*XI.^2);
end

function BWf = fillholes_nip(BW)
% Fill 0-holes by border flood-fill on zeros (4-neighborhood).
BW = logical(BW); Z = ~BW;
edge=false(size(BW)); edge([1 end],:)=Z([1 end],:); edge(:,[1 end])=edge(:,[1 end])|Z(:,[1 end]);
K = [0 1 0; 1 1 1; 0 1 0]; changed=true;
while changed
    grow = conv2(double(edge),K,'same')>0 & Z & ~edge;
    changed = any(grow(:)); edge = edge | grow;
end
holes = Z & ~edge; BWf = BW | holes;
end

function BW2 = majority3x3_nip(BW)
S = conv2(double(BW), ones(3), 'same'); BW2 = BW; BW2(S>=5)=true; BW2(S<=4)=false;
end

function BW2 = areaopen_nip(BW, minsz)
BW = logical(BW); [nr,nc]=size(BW); vis=false(nr,nc); BW2=BW; nbr=[0 1;0 -1;1 0;-1 0];
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
    if size(comp,1) < minsz
      BW2(sub2ind([nr nc], comp(:,1), comp(:,2))) = false;
    end
  end
end, end
end

function loops = contour_to_loops(C)
i=1; loops={};
while i < size(C,2)
    n = C(2,i); loops{end+1} = C(:,i+1:i+n).'; %#ok<AGROW>
    i = i + n + 1;
end
end

function A = signed_area(x,y)
A = 0.5*sum( x(:).*circshift(y(:),-1) - y(:).*circshift(x(:),-1) );
end

function P = rdp_simplify(P, eps)
% Ramer–Douglas–Peucker polyline simplification.
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
