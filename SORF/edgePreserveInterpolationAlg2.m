function [agirlikMatrisi] = edgePreserveInterpolationAlg2(refImg, agirlikMatrisi,distPower)

agirlikMatrisiYedek = agirlikMatrisi;      
[wm,hm]=size(refImg);
if nargin < 3
    distPower = 1;
end
%% ilk 4 satir
% sol üst köþe
i=2; 
j=2;
[b1,b2] = meshgrid(i-1:i+1,j-1:j+1);
[c1,c2] = meshgrid(i:3:i+3,j:3:j+3);
b11 = b1(:);
b22 = b2(:);
c11 = c1(:);
c22 = c2(:);
for ii=1:1:length(b11)
    distVals(:,ii) = sqrt( (b11(ii) - c11).^2 + (b22(ii)-c22).^2 );
end

distVals = distVals./repmat(sum(distVals,1),[4,1]);

nPxls = refImg(i-1:i+1,j-1:j+1);
cPxls = refImg(i:3:i+3,j:3:j+3);
                
gDiff = abs( repmat(nPxls(:)',[4,1]) - repmat(cPxls(:),[1,9]) );
gDiffN = gDiff./repmat(sum(gDiff,1),[4,1]);
        
gSim = 1./(gDiffN + repmat(min(gDiffN+(gDiffN==0)),[4,1])*0.1);
gSimN = gSim./repmat(sum(gSim,1),[4,1]);
        
simMatris = gSimN.*distVals;
        
simMatrisN = simMatris./(repmat(sum(simMatris,1),[4,1]));
        
cHSPxls = agirlikMatrisi(i:3:i+3,j:3:j+3);
        
a_cHSPxls = cHSPxls(:);
        
newPxls = simMatrisN'*a_cHSPxls;
        
centerPxl = agirlikMatrisi(i,j);
agirlikMatrisi(i-1:i+1,j-1:j+1) = reshape(newPxls,[3,3]);
agirlikMatrisi(i,j) = centerPxl;

clear i j b1 b2 c1 c2 b11 b22 c11 c22 ii distVals nPxls cPxls gDiff ...
    gDiffN gSim gSimN simMatris simMatrisN cHSPxls a_cHSPxls newPxls ...
    centerPxl 
% ----------- ilk üç satýr --------- %
i=2; 
j=5;
[b1,b2] = meshgrid(i-1:i+1,j-1:j+1);
[c1,c2] = meshgrid(i:3:i+3,j-3:3:j+3);
b11 = b1(:);
b22 = b2(:);
c11 = c1(:);
c22 = c2(:);
distVals = [];
for ii=1:1:length(b11)
    distVals(:,ii) = sqrt( (b11(ii) - c11).^2 + (b22(ii)-c22).^2 );
end

distVals = distVals./repmat(sum(distVals,1),[6,1]);

for i=2
    for j=5:3:hm-4
        
        nPxls = refImg(i-1:i+1,j-1:j+1);
        cPxls = refImg(i:3:i+3,j-3:3:j+3);
                
        gDiff = abs( repmat(nPxls(:)',[6,1]) - repmat(cPxls(:),[1,9]) );
        gDiffN = gDiff./repmat(sum(gDiff,1),[6,1]);
        
        gSim = 1./(gDiffN + repmat(min(gDiffN+(gDiffN==0)),[6,1])*0.1);
        gSimN = gSim./repmat(sum(gSim,1),[6,1]);
        
        simMatris = gSimN.*distVals;
        
        simMatrisN = simMatris./(repmat(sum(simMatris,1),[6,1]));
        
        cHSPxls = agirlikMatrisi(i:3:i+3,j-3:3:j+3);
        
        a_cHSPxls = cHSPxls(:);
        
        newPxls = simMatrisN'*a_cHSPxls;
        
        centerPxl = agirlikMatrisi(i,j);
        agirlikMatrisi(i-1:i+1,j-1:j+1) = reshape(newPxls,[3,3]);
        agirlikMatrisi(i,j) = centerPxl;
        
    end
end

clear i j b1 b2 c1 c2 b11 b22 c11 c22 ii distVals nPxls cPxls gDiff ...
    gDiffN gSim gSimN simMatris simMatrisN cHSPxls a_cHSPxls newPxls ...
    centerPxl 

% sað üst köþe
i=2; 
j=hm-1;
[b1,b2] = meshgrid(i-1:i+1,j-1:j+1);
[c1,c2] = meshgrid(i:3:i+3,j-3:3:j);
b11 = b1(:);
b22 = b2(:);
c11 = c1(:);
c22 = c2(:);
for ii=1:1:length(b11)
    distVals(:,ii) = sqrt( (b11(ii) - c11).^2 + (b22(ii)-c22).^2 );
end

distVals = distVals./repmat(sum(distVals,1),[4,1]);

nPxls = refImg(i-1:i+1,j-1:j+1);
cPxls = refImg(i:3:i+3,j-3:3:j);
                
gDiff = abs( repmat(nPxls(:)',[4,1]) - repmat(cPxls(:),[1,9]) );
gDiffN = gDiff./repmat(sum(gDiff,1),[4,1]);
        
gSim = 1./(gDiffN + repmat(min(gDiffN+(gDiffN==0)),[4,1])*0.1);
gSimN = gSim./repmat(sum(gSim,1),[4,1]);
        
simMatris = gSimN.*distVals;
        
simMatrisN = simMatris./(repmat(sum(simMatris,1),[4,1]));
        
cHSPxls = agirlikMatrisi(i:3:i+3,j-3:3:j);
        
a_cHSPxls = cHSPxls(:);
        
newPxls = simMatrisN'*a_cHSPxls;
        
centerPxl = agirlikMatrisi(i,j);
agirlikMatrisi(i-1:i+1,j-1:j+1) = reshape(newPxls,[3,3]);
agirlikMatrisi(i,j) = centerPxl;

clear i j b1 b2 c1 c2 b11 b22 c11 c22 ii distVals nPxls cPxls gDiff ...
    gDiffN gSim gSimN simMatris simMatrisN cHSPxls a_cHSPxls newPxls ...
    centerPxl 

%% son 4 satýr
% sol alt köþe
i=wm-1; 
j=2;
[b1,b2] = meshgrid(i-1:i+1,j-1:j+1);
[c1,c2] = meshgrid(i-3:3:i,j:3:j+3);
b11 = b1(:);
b22 = b2(:);
c11 = c1(:);
c22 = c2(:);
for ii=1:1:length(b11)
    distVals(:,ii) = sqrt( (b11(ii) - c11).^2 + (b22(ii)-c22).^2 );
end

distVals = distVals./repmat(sum(distVals,1),[4,1]);

nPxls = refImg(i-1:i+1,j-1:j+1);
cPxls = refImg(i-3:3:i,j:3:j+3);
                
gDiff = abs( repmat(nPxls(:)',[4,1]) - repmat(cPxls(:),[1,9]) );
gDiffN = gDiff./repmat(sum(gDiff,1),[4,1]);
        
gSim = 1./(gDiffN + repmat(min(gDiffN+(gDiffN==0)),[4,1])*0.1);
gSimN = gSim./repmat(sum(gSim,1),[4,1]);
        
simMatris = gSimN.*distVals;
        
simMatrisN = simMatris./(repmat(sum(simMatris,1),[4,1]));
        
cHSPxls = agirlikMatrisi(i-3:3:i,j:3:j+3);
        
a_cHSPxls = cHSPxls(:);
        
newPxls = simMatrisN'*a_cHSPxls;
        
centerPxl = agirlikMatrisi(i,j);
agirlikMatrisi(i-1:i+1,j-1:j+1) = reshape(newPxls,[3,3]);
agirlikMatrisi(i,j) = centerPxl;

clear i j b1 b2 c1 c2 b11 b22 c11 c22 ii distVals nPxls cPxls gDiff ...
    gDiffN gSim gSimN simMatris simMatrisN cHSPxls a_cHSPxls newPxls ...
    centerPxl 
% ----------- son üç satýr --------- %
i=wm-1; 
j=5;
[b1,b2] = meshgrid(i-1:i+1,j-1:j+1);
[c1,c2] = meshgrid(i-3:3:i,j-3:3:j+3);
b11 = b1(:);
b22 = b2(:);
c11 = c1(:);
c22 = c2(:);
distVals = [];
for ii=1:1:length(b11)
    distVals(:,ii) = sqrt( (b11(ii) - c11).^2 + (b22(ii)-c22).^2 );
end

distVals = distVals./repmat(sum(distVals,1),[6,1]);

for i=wm-1
    for j=5:3:hm-4
        
        nPxls = refImg(i-1:i+1,j-1:j+1);
        cPxls = refImg(i-3:3:i,j-3:3:j+3);
                
        gDiff = abs( repmat(nPxls(:)',[6,1]) - repmat(cPxls(:),[1,9]) );
        gDiffN = gDiff./repmat(sum(gDiff,1),[6,1]);
        
        gSim = 1./(gDiffN + repmat(min(gDiffN+(gDiffN==0)),[6,1])*0.1);
        gSimN = gSim./repmat(sum(gSim,1),[6,1]);
        
        simMatris = gSimN.*distVals;
        
        simMatrisN = simMatris./(repmat(sum(simMatris,1),[6,1]));
        
        cHSPxls = agirlikMatrisi(i-3:3:i,j-3:3:j+3);
        
        a_cHSPxls = cHSPxls(:);
        
        newPxls = simMatrisN'*a_cHSPxls;
        
        centerPxl = agirlikMatrisi(i,j);
        agirlikMatrisi(i-1:i+1,j-1:j+1) = reshape(newPxls,[3,3]);
        agirlikMatrisi(i,j) = centerPxl;
        
    end
end

clear i j b1 b2 c1 c2 b11 b22 c11 c22 ii distVals nPxls cPxls gDiff ...
    gDiffN gSim gSimN simMatris simMatrisN cHSPxls a_cHSPxls newPxls ...
    centerPxl 

% sað alt köþe
i=wm-1; 
j=hm-1;
[b1,b2] = meshgrid(i-1:i+1,j-1:j+1);
[c1,c2] = meshgrid(i-3:3:i,j-3:3:j);
b11 = b1(:);
b22 = b2(:);
c11 = c1(:);
c22 = c2(:);
for ii=1:1:length(b11)
    distVals(:,ii) = sqrt( (b11(ii) - c11).^2 + (b22(ii)-c22).^2 );
end

distVals = distVals./repmat(sum(distVals,1),[4,1]);

nPxls = refImg(i-1:i+1,j-1:j+1);
cPxls = refImg(i-3:3:i,j-3:3:j);
                
gDiff = abs( repmat(nPxls(:)',[4,1]) - repmat(cPxls(:),[1,9]) );
gDiffN = gDiff./repmat(sum(gDiff,1),[4,1]);
        
gSim = 1./(gDiffN + repmat(min(gDiffN+(gDiffN==0)),[4,1])*0.1);
gSimN = gSim./repmat(sum(gSim,1),[4,1]);
        
simMatris = gSimN.*distVals;
        
simMatrisN = simMatris./(repmat(sum(simMatris,1),[4,1]));
        
cHSPxls = agirlikMatrisi(i-3:3:i,j-3:3:j);
        
a_cHSPxls = cHSPxls(:);
        
newPxls = simMatrisN'*a_cHSPxls;
        
centerPxl = agirlikMatrisi(i,j);
agirlikMatrisi(i-1:i+1,j-1:j+1) = reshape(newPxls,[3,3]);
agirlikMatrisi(i,j) = centerPxl;

clear i j b1 b2 c1 c2 b11 b22 c11 c22 ii distVals nPxls cPxls gDiff ...
    gDiffN gSim gSimN simMatris simMatrisN cHSPxls a_cHSPxls newPxls ...
    centerPxl 

%% ilk 3 sütun

i=5; 
j=2;
[b1,b2] = meshgrid(i-1:i+1,j-1:j+1);
[c1,c2] = meshgrid(i-3:3:i+3,j:3:j+3);
b11 = b1(:);
b22 = b2(:);
c11 = c1(:);
c22 = c2(:);
distVals = [];
for ii=1:1:length(b11)
    distVals(:,ii) = sqrt( (b11(ii) - c11).^2 + (b22(ii)-c22).^2 );
end

distVals = distVals./repmat(sum(distVals,1),[6,1]);

for j=2
    for i=5:3:wm-4
            
        nPxls = refImg(i-1:i+1,j-1:j+1);
        cPxls = refImg(i-3:3:i+3,j:3:j+3);
                
        gDiff = abs( repmat(nPxls(:)',[6,1]) - repmat(cPxls(:),[1,9]) );
        gDiffN = gDiff./repmat(sum(gDiff,1),[6,1]);
        
        gSim = 1./(gDiffN + repmat(min(gDiffN+(gDiffN==0)),[6,1])*0.1);
        gSimN = gSim./repmat(sum(gSim,1),[6,1]);
        
        simMatris = gSimN.*distVals;
        
        simMatrisN = simMatris./(repmat(sum(simMatris,1),[6,1]));
        
        cHSPxls = agirlikMatrisi(i-3:3:i+3,j:3:j+3);
        
        a_cHSPxls = cHSPxls(:);
        
        newPxls = simMatrisN'*a_cHSPxls;
        
        centerPxl = agirlikMatrisi(i,j);
        agirlikMatrisi(i-1:i+1,j-1:j+1) = reshape(newPxls,[3,3]);
        agirlikMatrisi(i,j) = centerPxl;
        
    end
end

clear i j b1 b2 c1 c2 b11 b22 c11 c22 ii distVals nPxls cPxls gDiff ...
    gDiffN gSim gSimN simMatris simMatrisN cHSPxls a_cHSPxls newPxls ...
    centerPxl 

%% son 3 sütun

i=5; 
j=hm-1;
[b1,b2] = meshgrid(i-1:i+1,j-1:j+1);
[c1,c2] = meshgrid(i-3:3:i+3,j-3:3:j);
b11 = b1(:);
b22 = b2(:);
c11 = c1(:);
c22 = c2(:);
distVals = [];
for ii=1:1:length(b11)
    distVals(:,ii) = sqrt( (b11(ii) - c11).^2 + (b22(ii)-c22).^2 );
end

distVals = distVals./repmat(sum(distVals,1),[6,1]);

for j=hm-1
    for i=5:3:wm-4
            
        nPxls = refImg(i-1:i+1,j-1:j+1);
        cPxls = refImg(i-3:3:i+3,j-3:3:j);
                
        gDiff = abs( repmat(nPxls(:)',[6,1]) - repmat(cPxls(:),[1,9]) );
        gDiffN = gDiff./repmat(sum(gDiff,1),[6,1]);
        
        gSim = 1./(gDiffN + repmat(min(gDiffN+(gDiffN==0)),[6,1])*0.1);
        gSimN = gSim./repmat(sum(gSim,1),[6,1]);
        
        simMatris = gSimN.*distVals;
        
        simMatrisN = simMatris./(repmat(sum(simMatris,1),[6,1]));
        
        cHSPxls = agirlikMatrisi(i-3:3:i+3,j-3:3:j);
        
        a_cHSPxls = cHSPxls(:);
        
        newPxls = simMatrisN'*a_cHSPxls;
        
        centerPxl = agirlikMatrisi(i,j);
        agirlikMatrisi(i-1:i+1,j-1:j+1) = reshape(newPxls,[3,3]);
        agirlikMatrisi(i,j) = centerPxl;
    end
end

clear i j b1 b2 c1 c2 b11 b22 c11 c22 ii distVals nPxls cPxls gDiff ...
    gDiffN gSim gSimN simMatris simMatrisN cHSPxls a_cHSPxls newPxls ...
    centerPxl 



%% pikseller

i=5; j=5;
[b1,b2] = meshgrid(i-1:i+1,j-1:j+1);
[c1,c2] = meshgrid(i-3:3:i+3,j-3:3:j+3);
b11 = b1(:);
b22 = b2(:);
c11 = c1(:);
c22 = c2(:);
for i=1:1:length(b11)
    distVals(:,i) = sqrt( (b11(i) - c11).^2 + (b22(i)-c22).^2 );
end

distVals = distVals./repmat(sum(distVals,1),[9,1]);

for i=5:3:wm-4
    for j=5:3:hm-4
        
        nPxls = refImg(i-1:i+1,j-1:j+1);
        cPxls = refImg(i-3:3:i+3,j-3:3:j+3);
                
        gDiff = abs( repmat(nPxls(:)',[9,1]) - repmat(cPxls(:),[1,9]) );
        gDiffN = gDiff./repmat(sum(gDiff,1),[9,1]);
        
        gSim = 1./(gDiffN + repmat(min(gDiffN+(gDiffN==0)),[9,1])*0.1);
        gSimN = gSim./repmat(sum(gSim,1),[9,1]);
        
        simMatris = gSimN.*distVals;
        
        simMatrisN = simMatris./(repmat(sum(simMatris,1),[9,1]));
        
        cHSPxls = agirlikMatrisi(i-3:3:i+3,j-3:3:j+3);
        
        a_cHSPxls = cHSPxls(:);
        
        newPxls = simMatrisN'*a_cHSPxls;
        
        centerPxl = agirlikMatrisi(i,j);
        agirlikMatrisi(i-1:i+1,j-1:j+1) = reshape(newPxls,[3,3]);
        agirlikMatrisi(i,j) = centerPxl;
        
    end
end

