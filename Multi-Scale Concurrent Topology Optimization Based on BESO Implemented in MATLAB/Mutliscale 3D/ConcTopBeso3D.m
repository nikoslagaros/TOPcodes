function [c,CH,x,y] = ConcTopBeso3D(lx,ly,lz,nelx,nely,nelz,nlx,nly,nlz, ...
    volfx,volfy,penal,rxmin,rymin,er,p)
% USER-DEFINED LOOP PARAMETERS
maxloop = 200;    % Maximum number of iterations
tolx = 0.01;      % Termination criterion
displayflag = 1;  % Display structure flag
% USER-DEFINED MATERIAL PROPERTIES
E0 = 2;           % Young's modulus of solid material
Emin = 1e-8;      % Young's modulus of void-like material
xmin = 1e-8;
nu = 0.3;         % Poisson's ratio
% USER-DEFINED LOAD DOFs
il = nelx/2; jl = nely; kl = nelz/2;                    % Coordinates
loadnid = kl*(nelx+1)*(nely+1)+il*(nely+1)+(nely+1-jl); % Node IDs
loaddof = 3*loadnid(:) - 1;                             % DOFs
% USER-DEFINED SUPPORT FIXED DOFs
iif = [0 0 nelx nelx]; jf = [0 0 0 0]; kf = [0 nelz 0 nelz];  % Coordinates
fixednid = kf*(nelx+1)*(nely+1)+iif*(nely+1)+(nely+1-jf);     % Node IDs
fixeddof = [3*fixednid(:); 3*fixednid(:)-1; 3*fixednid(:)-2]; % DOFs
% PREPARE FINITE ELEMENT ANALYSIS
nele = nelx*nely*nelz;
ndof = 3*(nelx+1)*(nely+1)*(nelz+1);
F = sparse(loaddof,1,-1,ndof,1);
U = zeros(ndof,size(F,2));
freedofs = setdiff(1:ndof,fixeddof);
h1 = hexahedron(lx/nelx/2,ly/nely/2,lz/nelz/2);
h2 = hexahedron(1e-6/nlx/2,1e-6/nly/2,1e-6/nlz/2);
nodegrd = reshape(1:(nely+1)*(nelx+1),nely+1,nelx+1);
nodeids = reshape(nodegrd(1:end-1,1:end-1),nely*nelx,1);
nodeidz = 0:(nely+1)*(nelx+1):(nelz-1)*(nely+1)*(nelx+1);
nodeids = repmat(nodeids,size(nodeidz))+repmat(nodeidz,size(nodeids));
edofVec = 3*nodeids(:)+1;
edofMat = repmat(edofVec,1,24)+ ...
    repmat([0 1 2 3*nely + [3 4 5 0 1 2] -3 -2 -1 ...
    3*(nely+1)*(nelx+1)+[0 1 2 3*nely + [3 4 5 0 1 2] -3 -2 -1]],nele,1);
iK = kron(edofMat,ones(24,1))';
jK = kron(edofMat,ones(1,24))';
% PREPARE FILTER MACRO
iH = ones(nele*(2*(ceil(rxmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for k1 = 1:nelz
    for i1 = 1:nelx
        for j1 = 1:nely
            e1 = (k1-1)*nelx*nely + (i1-1)*nely+j1;
            for k2 = max(k1-(ceil(rxmin)-1),1):min(k1+(ceil(rxmin)-1),nelz)
                for i2 = max(i1-(ceil(rxmin)-1),1):min(i1+(ceil(rxmin)-1),nelx)
                    for j2 = max(j1-(ceil(rxmin)-1),1):min(j1+(ceil(rxmin)-1),nely)
                        e2 = (k2-1)*nelx*nely + (i2-1)*nely+j2;
                        k = k+1;
                        iH(k) = e1;
                        jH(k) = e2;
                        sH(k) = max(0,rxmin-sqrt((i1-i2)^2+(j1-j2)^2+(k1-k2)^2));
                    end
                end
            end
        end
    end
end
Hx = sparse(iH,jH,sH);
Hsx = sum(Hx,2);
% PREPARE FILTER MICRO
iH = ones(nlx*nly*nlz*(2*(ceil(rymin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for k1 = 1:nlz
    for i1 = 1:nlx
        for j1 = 1:nly
            e1 = (k1-1)*nlx*nly + (i1-1)*nly+j1;
            for k2 = max(k1-(ceil(rymin)-1),1):min(k1+(ceil(rymin)-1),nlz)
                for i2 = max(i1-(ceil(rymin)-1),1):min(i1+(ceil(rymin)-1),nlx)
                    for j2 = max(j1-(ceil(rymin)-1),1):min(j1+(ceil(rymin)-1),nly)
                        e2 = (k2-1)*nlx*nly + (i2-1)*nly+j2;
                        k = k+1;
                        iH(k) = e1;
                        jH(k) = e2;
                        sH(k) = max(0,rymin-sqrt((i1-i2)^2+(j1-j2)^2+(k1-k2)^2));
                    end
                end
            end
        end
    end
end
Hy = sparse(iH,jH,sH);
Hsy = sum(Hy,2);
% INITIALIZE ITERATION
x = ones(nely,nelx,nelz);
y = repmat(ones,[nly,nlx,nlz]);
for i=1:nly
    for j=1:nlx
        for k=1:nlz
            if sqrt((i-nlx/2-0.5)^2+(j-nly/2-0.5)^2+(k-nlz/2-0.5)^2) < min(min(nlx,nly),nlz)/3
                y(i,j,k) = 0;
            end
        end
    end
end
loop = 0; 
change = 1;
volx = 1;
voly = 1;
% START ITERATION
while change > tolx && loop < maxloop
    loop = loop+1;
    % UPDATE VOLUME
    volx = max(volx*(1-er),volfx);
    voly = max(voly*(1-er),volfy);
    if loop > 1; olddc = dc; end
    if loop > 1; olddcy = dcy; end
    % INTERPOLATION
    [E,dE] = interpolate(y,E0,Emin,penal);
    % Homogenization
    [CH,DCH] = homo3DY(1e-6,1e-6,1e-6,E,dE,nu,1e-8,h2);
    % FE-ANALYSIS
    KE = h1.elemStiffness(CH);
    sK = reshape(KE(:)*(xmin+x(:)'.^penal*(1-xmin)),576*nelx*nely*nelz,1);
    K = sparse(iK(:),jK(:),sK(:)); K = (K+K')/2;
    U(freedofs) = p.solve(K(freedofs,freedofs),F(freedofs));
    % OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS MACRO
    ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),[nely,nelx,nelz]);
    c(loop) = sum(sum(sum((xmin+x.^penal*(1-xmin)).*ce)));
    dc = penal*(1-xmin)*x.^(penal-1).*ce;
    % SENSITIVITY ANALYSIS MICRO
    f = @(DCh) h1.elemStiffness(DCh);
    dKe = cellfun(f,DCH,'UniformOutput',false);
    f = @(dKE) reshape(sum((U(edofMat)*dKE).*U(edofMat),2),[nely,nelx,nelz]);
    cey = cellfun(f,dKe,'UniformOutput',false);
    dcy = cell2mat(cellfun(@sum,cellfun(@sum,cellfun(@sum, cey, ...
          'UniformOutput',false),'UniformOutput',false),'UniformOutput',false));
    % FILTERING AND MODIFICATION OF SENSITIVITIES MACRO
    dc(:) = Hx*(dc(:))./Hsx;
    if loop > 1; dc = (dc+olddc)/2.; end
    % FILTERING AND MODIFICATION OF SENSITIVITIES MICRO
    dcy(:) = Hy*(dcy(:)./Hsy);
    if loop > 1; dcy = (dcy+olddcy)/2.; end
    % BESO UPDATE MACRO
    l1 = min(dc(:)); l2 = max(dc(:));
    while (l2-l1)/(l1+l2) > 1e-3
        lmid = 0.5*(l2+l1);
        x = max(0.001,sign(dc-lmid));
        if sum(x(:)) > volx*nele, l1 = lmid; else, l2 = lmid; end
    end
    % BESO UPDATE MICRO
    l1 = min(dcy(:)); l2 = max(dcy(:));
    while (l2-l1)/(l1+l2) > 1e-4
        lmid = 0.5*(l2+l1);
        y = max(0.001,sign(dcy-lmid));
        if sum(y(:)) > voly*nlx*nly*nlz, l1 = lmid; else, l2 = lmid; end
    end
    if loop > 10
        change  = abs(sum(c(loop-9:loop-5))-sum(c(loop-4:loop)))/sum(c(loop-4:loop));
    end
    % PRINT RESULTS
    fprintf(' It.:%5i Obj.:%11.4f Volx.:%7.3f Voly.:%7.3f ch.:%7.3f\n', ...
        loop,c(loop),mean(x(:)),mean(y(:)),change);
    % PLOT DENSITIES
    if displayflag, figure(1); clf; display_3D(x); end
    if displayflag, figure(2); clf; display_3D(y); end
end
figure(1); clf; display_3D(x);
figure(2); clf; display_3D(y);
end