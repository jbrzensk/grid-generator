%% test_softC_heat_blob_2nd_vs_4th
clc; clear; close all;

thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir);
addpath(fullfile(thisDir,'Grids'));

%% --- Original parameters ---
Ns = 81;
Nt = 81;
N_outer = 0010;

alpha  = 1.0;
Tfinal = 1;
CFL    = 0.10;

grid_fun = @makeGridHardC;

%% --- Build optimized soft C grid quietly ---
figsBefore = findall(0,'Type','figure');
evalc('[C_opt, E_opt, Xopt, Yopt, Jc, pUopt, pVopt, Uopt, Vopt, info] = gridGenerator(grid_fun, Ns, Nt, N_outer);');
figsAfter = findall(0,'Type','figure');
close(setdiff(figsAfter, figsBefore));

%% --- Cell geometry ---
[Xcell,Ycell,Acell,folded] = cellGeometryPhysical(Xopt,Yopt);
NcY = Nt-1;
NcX = Ns-1;

fprintf('\nGRID:\n');
fprintf('  cells        = %d x %d\n', NcX, NcY);
fprintf('  folded cells = %d\n', nnz(folded));

%% --- Operators ---
[G2,D2,K2] = buildMimeticOperators2nd(Xopt,Yopt);
[G4,D4,K4] = buildMimeticOperators4th(Xopt,Yopt);

%% --- CFL timestep ---
faceSum2 = max(-full(diag(K2)),0);
faceSum4 = max(-full(diag(K4)),0);

dt2 = min(Acell(:) ./ max(alpha*faceSum2,1e-14));
dt4 = min(Acell(:) ./ max(alpha*faceSum4,1e-14));

dt = CFL * min(dt2,dt4);
nSteps = ceil(Tfinal/dt);
dt = Tfinal/nSteps;

fprintf('\nTIME STEP:\n');
fprintf('  dt     = %.6e\n', dt);
fprintf('  nSteps = %d\n', nSteps);

%% --- Initial heat blob ---
x0 = mean(Xcell(:));
y0 = mean(Ycell(:));
sigma = 0.25;

u0 = exp(-((Xcell - x0).^2 + (Ycell - y0).^2) / sigma^2);

u2 = u0;
u4 = u0;

plotCellField(Xopt,Yopt,u0,'Initial temperature');

%% --- Time stepping ---
for n = 1:nSteps
    rhs2 = heatRHS(u2,K2,Acell,alpha);
    u2 = u2 + dt*rhs2;

    k1 = heatRHS(u4,K4,Acell,alpha);
    k2 = heatRHS(u4 + 0.5*dt*k1,K4,Acell,alpha);
    k3 = heatRHS(u4 + 0.5*dt*k2,K4,Acell,alpha);
    k4 = heatRHS(u4 + dt*k3,K4,Acell,alpha);

    u4 = u4 + dt*(k1 + 2*k2 + 2*k3 + k4)/6;
end

%% --- Compare ---
err = u2 - u4;

fprintf('\n2ND VS 4TH ORDER:\n');
fprintf('  L2 difference   = %.6e\n', sqrt(sum(Acell(:).*err(:).^2)/sum(Acell(:))));
fprintf('  Linf difference = %.6e\n', max(abs(err(:))));

plotCellField(Xopt,Yopt,u2,'2nd-order heat solution');
plotCellField(Xopt,Yopt,u4,'4th-order D*G RK4 reference solution');
plotCellField(Xopt,Yopt,err,'Difference: 2nd - 4th');

%% =======================================================================
function rhs = heatRHS(u,K,Acell,alpha)
rhs = alpha * (K*u(:)) ./ Acell(:);
rhs = reshape(rhs,size(u));
end

function [Xcell,Ycell,Acell,folded] = cellGeometryPhysical(X,Y)
[Nt,Ns] = size(X);
NcY = Nt-1; NcX = Ns-1;

Xcell = zeros(NcY,NcX);
Ycell = zeros(NcY,NcX);
Acell = zeros(NcY,NcX);
folded = false(NcY,NcX);

for j = 1:NcY
    for i = 1:NcX
        px = [X(j,i), X(j,i+1), X(j+1,i+1), X(j+1,i)];
        py = [Y(j,i), Y(j,i+1), Y(j+1,i+1), Y(j+1,i)];

        signedA = 0.5*sum(px.*py([2 3 4 1]) - py.*px([2 3 4 1]));

        Xcell(j,i) = mean(px);
        Ycell(j,i) = mean(py);
        Acell(j,i) = abs(signedA);
        folded(j,i) = signedA <= 0;
    end
end
end

function [G,D,K] = buildMimeticOperators2nd(X,Y)
[Nt,Ns] = size(X);
NcY = Nt-1; NcX = Ns-1;
Ncells = NcY*NcX;

cellIdx = @(j,i) j + (i-1)*NcY;
[Xc,Yc,~,~] = cellGeometryPhysical(X,Y);

nVerticalFaces = NcY*(NcX-1);
nHorizontalFaces = (NcY-1)*NcX;
Nfaces = nVerticalFaces + nHorizontalFaces;

gRows=[]; gCols=[]; gVals=[];
dRows=[]; dCols=[]; dVals=[];
face = 0;

for j = 1:NcY
    for i = 1:NcX-1
        face = face+1;
        Lcell = cellIdx(j,i);
        Rcell = cellIdx(j,i+1);

        Lface = hypot(X(j+1,i+1)-X(j,i+1), Y(j+1,i+1)-Y(j,i+1));
        dcent = hypot(Xc(j,i+1)-Xc(j,i), Yc(j,i+1)-Yc(j,i));

        gRows = [gRows; face; face];
        gCols = [gCols; Lcell; Rcell];
        gVals = [gVals; -1/max(dcent,1e-14); 1/max(dcent,1e-14)];

        dRows = [dRows; Lcell; Rcell];
        dCols = [dCols; face; face];
        dVals = [dVals; Lface; -Lface];
    end
end

for j = 1:NcY-1
    for i = 1:NcX
        face = face+1;
        Bcell = cellIdx(j,i);
        Tcell = cellIdx(j+1,i);

        Lface = hypot(X(j+1,i+1)-X(j+1,i), Y(j+1,i+1)-Y(j+1,i));
        dcent = hypot(Xc(j+1,i)-Xc(j,i), Yc(j+1,i)-Yc(j,i));

        gRows = [gRows; face; face];
        gCols = [gCols; Bcell; Tcell];
        gVals = [gVals; -1/max(dcent,1e-14); 1/max(dcent,1e-14)];

        dRows = [dRows; Bcell; Tcell];
        dCols = [dCols; face; face];
        dVals = [dVals; Lface; -Lface];
    end
end

G = sparse(gRows,gCols,gVals,Nfaces,Ncells);
D = sparse(dRows,dCols,dVals,Ncells,Nfaces);
K = D*G;
end

function [G,D,K] = buildMimeticOperators4th(X,Y)
[Nt,Ns] = size(X);
NcY = Nt-1; NcX = Ns-1;
Ncells = NcY*NcX;

cellIdx = @(j,i) j + (i-1)*NcY;
[Xc,Yc,~,~] = cellGeometryPhysical(X,Y);

nVerticalFaces = NcY*(NcX-1);
nHorizontalFaces = (NcY-1)*NcX;
Nfaces = nVerticalFaces + nHorizontalFaces;

gRows=[]; gCols=[]; gVals=[];
dRows=[]; dCols=[]; dVals=[];
face = 0;

for j = 1:NcY
    for i = 1:NcX-1
        face = face+1;

        Lcell = cellIdx(j,i);
        Rcell = cellIdx(j,i+1);

        Lface = hypot(X(j+1,i+1)-X(j,i+1), Y(j+1,i+1)-Y(j,i+1));

        if i >= 2 && i <= NcX-2
            ids = [cellIdx(j,i-1), cellIdx(j,i), cellIdx(j,i+1), cellIdx(j,i+2)];
            d = hypot(Xc(j,i+1)-Xc(j,i), Yc(j,i+1)-Yc(j,i));
            coeffs = [1 -8 8 -1] / (12*max(d,1e-14));
        else
            ids = [Lcell, Rcell];
            d = hypot(Xc(j,i+1)-Xc(j,i), Yc(j,i+1)-Yc(j,i));
            coeffs = [-1 1] / max(d,1e-14);
        end

        gRows = [gRows; face*ones(numel(ids),1)];
        gCols = [gCols; ids(:)];
        gVals = [gVals; coeffs(:)];

        dRows = [dRows; Lcell; Rcell];
        dCols = [dCols; face; face];
        dVals = [dVals; Lface; -Lface];
    end
end

for j = 1:NcY-1
    for i = 1:NcX
        face = face+1;

        Bcell = cellIdx(j,i);
        Tcell = cellIdx(j+1,i);

        Lface = hypot(X(j+1,i+1)-X(j+1,i), Y(j+1,i+1)-Y(j+1,i));

        if j >= 2 && j <= NcY-2
            ids = [cellIdx(j-1,i), cellIdx(j,i), cellIdx(j+1,i), cellIdx(j+2,i)];
            d = hypot(Xc(j+1,i)-Xc(j,i), Yc(j+1,i)-Yc(j,i));
            coeffs = [1 -8 8 -1] / (12*max(d,1e-14));
        else
            ids = [Bcell, Tcell];
            d = hypot(Xc(j+1,i)-Xc(j,i), Yc(j+1,i)-Yc(j,i));
            coeffs = [-1 1] / max(d,1e-14);
        end

        gRows = [gRows; face*ones(numel(ids),1)];
        gCols = [gCols; ids(:)];
        gVals = [gVals; coeffs(:)];

        dRows = [dRows; Bcell; Tcell];
        dCols = [dCols; face; face];
        dVals = [dVals; Lface; -Lface];
    end
end

G = sparse(gRows,gCols,gVals,Nfaces,Ncells);
D = sparse(dRows,dCols,dVals,Ncells,Nfaces);
K = D*G;
end

function plotCellField(X,Y,u,ttl)
[Nt,Ns] = size(X);
NcY = Nt-1; NcX = Ns-1;

figure('Color','w'); hold on; axis equal; box on;

for j = 1:NcY
    for i = 1:NcX
        px = [X(j,i), X(j,i+1), X(j+1,i+1), X(j+1,i)];
        py = [Y(j,i), Y(j,i+1), Y(j+1,i+1), Y(j+1,i)];

        patch(px,py,u(j,i),'EdgeColor',[0.65 0.65 0.65]);
    end
end

axis equal tight;
colorbar;
title(ttl);
xlabel('x'); ylabel('y');
end