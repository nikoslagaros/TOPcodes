function [c,CH,x,y] = ConcTopOptBESOV2(lx,ly,nelx,nely,nlx,nly,volf, ...
    penal,rxmin,rymin,er)
%% MATERIAL PROPERTIES
E0 = 2;
Emin = 1.5;
xmin = 1e-8;
nu = 0.3;
%% PREPARE FINITE ELEMENT ANALYSIS
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
% DEFINE LOADS AND SUPPORTS
% F = sparse(2,1,-1,2*(nely+1)*(nelx+1),1);
% fixeddofs = union(1:2:2*(nely+1),2*(nelx+1)*(nely+1));
% DEFINE LOADS AND SUPPORTS (MBBeam)
F = sparse(2,1,-1,2*(nely+1)*(nelx+1),1);
fixeddofs = union(1:2:2*(nely+1),2*(nelx+1)*(nely+1));
% % DEFINE LOADS AND SUPPORTS (Long Cantilever)
% F = sparse(2*(nelx*(nely+1)+nely/2+1),1,-1,2*(nely+1)*(nelx+1),1);
% fixeddofs = 1:2*(nely+1);
U = zeros(2*(nely+1)*(nelx+1),size(F,2));
alldofs = 1:2*(nely+1)*(nelx+1);
freedofs = setdiff(alldofs,fixeddofs);
%% PREPARE FILTER GLOBAL
iH = ones(nelx*nely*(2*(ceil(rxmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for i1 = 1:nelx
  for j1 = 1:nely
    e1 = (i1-1)*nely+j1;
    for i2 = max(i1-(ceil(rxmin)-1),1):min(i1+(ceil(rxmin)-1),nelx)
      for j2 = max(j1-(ceil(rxmin)-1),1):min(j1+(ceil(rxmin)-1),nely)
        e2 = (i2-1)*nely+j2;
        k = k+1;
        iH(k) = e1;
        jH(k) = e2;
        sH(k) = max(0,rxmin-sqrt((i1-i2)^2+(j1-j2)^2));
      end
    end
  end
end
Hx = sparse(iH,jH,sH);
Hsx = sum(Hx,2);
%% PREPARE FILTER UNIT CELL
iH = ones(nlx*nly*(2*(ceil(rymin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for i1 = 1:nlx
  for j1 = 1:nly
    e1 = (i1-1)*nly+j1;
    for i2 = max(i1-(ceil(rymin)-1),1):min(i1+(ceil(rymin)-1),nlx)
      for j2 = max(j1-(ceil(rymin)-1),1):min(j1+(ceil(rymin)-1),nly)
        e2 = (i2-1)*nly+j2;
        k = k+1;
        iH(k) = e1;
        jH(k) = e2;
        sH(k) = max(0,rymin-sqrt((i1-i2)^2+(j1-j2)^2));
      end
    end
  end
end
Hy = sparse(iH,jH,sH);
Hsy = sum(Hy,2);
%% INITIALIZE ITERATION
x = ones(nely,nelx);
y = ones(nly,nlx);
for i = 1:nlx
    for j = 1:nly
        if sqrt((i-nlx/2-0.5)^2+(j-nly/2-0.5)^2) < min(nlx,nly)/3
            y(j,i) = 0;
        end
    end
end
loop = 0;
change = 1;
vol = 1;
%% START ITERATION
while change > 0.001 && loop < 200
  loop = loop + 1;
  %% UPDATE VOLUME
  vol = max(vol*(1-er),volf);
  if loop > 1; olddc = dc; end
  if loop > 1; olddcy = dcy; end
  %% INTERPOLATION
  [E, dE] = interpolate(y,E0,Emin,penal);
  %% HOMOGENIZATION
  [CH,DCH] = homogenize(1e-6,1e-6,E,nu,dE,90);
  %% FE-ANALYSIS
  KE = elementMatVec(lx/nelx/2,ly/nely/2,90,CH);
  sK = reshape(KE(:)*(xmin+x(:)'.^penal*(1-xmin)),64*nelx*nely,1);
  K = sparse(iK,jK,sK); K = (K+K')/2;
  U(freedofs,:) = K(freedofs,freedofs)\F(freedofs,:);
  %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS MACRO
  ce = reshape(0.5*sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
  c(loop) = sum(sum((xmin+x.^penal*(1-xmin)).*ce));
  dc = penal*(1-xmin)*x.^(penal-1).*ce;
  %% SENSITIVITY ANALYSIS MICRO
  f = @(DCh) elementMatVec(lx/nelx/2,ly/nely/2,90,DCh);
  dKe = cellfun(f,DCH,'UniformOutput',false);
  f= @(dKE) reshape(0.5*sum((U(edofMat)*dKE).*U(edofMat),2),nely,nelx);
  cey = cellfun(f,dKe,'UniformOutput',false);
  dcy = cell2mat(cellfun(@sum,cellfun(@sum,cellfun(@sum, cey, 'UniformOutput',false),'UniformOutput',false),'UniformOutput',false));
  %% FILTERING/MODIFICATION OF SENSITIVITIES MACRO
  dc(:) = Hx*(dc(:))./Hsx;
  if loop > 1; dc = (dc+olddc)/2.; end
  %% FILTERING/MODIFICATION OF SENSITIVITIES MICRO
  dcy(:) = Hy*(dcy(:)./Hsy);
  if loop > 1; dcy = (dcy+olddcy)/2.; end
  %% BESO UPDATE
  l1a = min(min(dc)); l1b = min(min(dcy)); l1 = min(l1a,l1b);
  l2a = max(max(dc)); l2b = max(max(dcy)); l2 = max(l2a,l2b);
  while (l2-l1)/(l1+l2) > 1e-3
    lmid = 0.5*(l2+l1);
    x = max(0.001,sign(dc-lmid));
    y = max(0.001,sign(dcy-lmid));
    if sum(x(:))/nelx/nely * sum(y(:)+(1 - y(:))*1/8)/nlx/nly > vol, l1 = lmid; else, l2 = lmid; end
  end
  if loop > 10
      change  = abs(sum(c(loop-9:loop-5))-sum(c(loop-4:loop)))/sum(c(loop-4:loop));
  end
  %% PRINT RESULTS
  fprintf(' It.:%5i Obj.:%11.4f Volx.:%7.3f Voly.:%7.3f ch.:%7.3f\n',loop,c(loop), ...
    mean(x(:)),mean(y(:)+(1 - y(:))*1/8),change);
  %% PLOT DENSITIES
  figure(1)
  colormap(gray); imagesc(1-x); caxis([0 1]); axis equal; axis off; drawnow;
  figure(2)
  map = [0 0 0
    0.7 0.7 0.7];
  colormap(map); imagesc(1-y); caxis([0 1]); axis equal; axis off; drawnow;
end
end