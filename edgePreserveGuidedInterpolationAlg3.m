function [agirlikMatrisi] = edgePreserveGuidedInterpolationAlg3(refImg, agirlikMatrisi,distPower)

agirlikMatrisiYedek = agirlikMatrisi;      
[wm,hm]=size(refImg);
if nargin < 3
    distPower = 1;
end

%% satirlar
i=2; 
j=2;
[b1,b2] = meshgrid(i,j+1:j+2);
[c1,c2] = meshgrid(i,j:3:j+3);
b11 = b1(:);
b22 = b2(:);
c11 = c1(:);
c22 = c2(:);
distVals = [];
for ii=1:1:length(b11)
    distVals(:,ii) = sqrt( (b11(ii) - c11).^2 + (b22(ii)-c22).^2 );
end

distVals = distVals./repmat(sum(distVals,1),[2,1]);

for i=2:3:wm-2
    for j=2:3:hm-2
        
        nPxls = refImg(i,j+1:j+2);
        cPxls = refImg(i,j:3:j+3);
                
        gDiff = abs( repmat(nPxls(:)',[2,1]) - repmat(cPxls(:),[1,2]) );
        gDiffN = gDiff./repmat(sum(gDiff,1)+eps,[2,1]);
        
        gSim = 1./(gDiffN + repmat(min(gDiffN+(gDiffN==0)),[2,1])*0.1);
        gSimN = gSim./repmat(sum(gSim,1),[2,1]);
        
        simMatris = gSimN.*(distVals.^distPower);
        
        simMatrisN = simMatris./(repmat(sum(simMatris,1),[2,1]));
        
        cHSPxls = agirlikMatrisi(i,j:3:j+3);
        
        a_cHSPxls = cHSPxls(:);
        
        newPxls = simMatrisN'*a_cHSPxls;
        
        if sum(isnan(newPxls))>=1
            'dur';
        end
        
        agirlikMatrisi(i,j+1:j+2) = newPxls';
        
    end
end

clear i j b1 b2 c1 c2 b11 b22 c11 c22 ii distVals nPxls cPxls gDiff ...
    gDiffN gSim gSimN simMatris simMatrisN cHSPxls a_cHSPxls newPxls ...
    centerPxl 

%% sütunlar
i=2; 
j=2;
[b1,b2] = meshgrid(i+1:i+2,j);
[c1,c2] = meshgrid(i:3:i+3,j);
b11 = b1(:);
b22 = b2(:);
c11 = c1(:);
c22 = c2(:);
distVals = [];
for ii=1:1:length(b11)
    distVals(:,ii) = sqrt( (b11(ii) - c11).^2 + (b22(ii)-c22).^2 );
end

distVals = distVals./repmat(sum(distVals,1),[2,1]);

for i=2:3:wm-2
    for j=2:3:hm-2
        
        nPxls = refImg(i+1:i+2,j);
        cPxls = refImg(i:3:i+3,j);
                
        gDiff = abs( repmat(nPxls(:)',[2,1]) - repmat(cPxls(:),[1,2]) );
        gDiffN = gDiff./repmat(sum(gDiff,1)+eps,[2,1]);
        
        gSim = 1./(gDiffN + repmat(min(gDiffN+(gDiffN==0)),[2,1])*0.1);
        gSimN = gSim./repmat(sum(gSim,1),[2,1]);
        
        simMatris = gSimN.*(distVals.^distPower);
        
        simMatrisN = simMatris./(repmat(sum(simMatris,1),[2,1]));
        
        cHSPxls = agirlikMatrisi(i:3:i+3,j);
        
        a_cHSPxls = cHSPxls(:);
        
        newPxls = simMatrisN'*a_cHSPxls;
        
        if sum(isnan(newPxls))>=1
            'dur';
        end
        
        agirlikMatrisi(i+1:i+2,j) = newPxls';
        
    end
end

clear i j b1 b2 c1 c2 b11 b22 c11 c22 ii distVals nPxls cPxls gDiff ...
    gDiffN gSim gSimN simMatris simMatrisN cHSPxls a_cHSPxls newPxls ...
    centerPxl 


%% diagonal pikseller

i=3; j=3;
[b1,b2] = meshgrid(i:i+1,j:j+1);
[c1,c2] = meshgrid(i-1:3:i+2,j-1:3:j+2);
b11 = b1(:);
b22 = b2(:);
c11 = c1(:);
c22 = c2(:);
distVals = [];
for ii=1:1:length(b11)
    distVals(:,ii) = sqrt( (b11(ii) - c11).^2 + (b22(ii)-c22).^2 );
end

distVals = distVals./repmat(sum(distVals,1),[4,1]);

for i=3:3:wm-3
    for j=3:3:hm-3
        
        nPxls = refImg(i:i+1,j:j+1);
        cPxls = refImg(i-1:3:i+2,j-1:3:j+2);
                
        gDiff = abs( repmat(nPxls(:)',[4,1]) - repmat(cPxls(:),[1,4]) );
        gDiffN = gDiff./repmat(sum(gDiff,1)+eps,[4,1]);
        
        gSim = 1./(gDiffN + repmat(min(gDiffN+(gDiffN==0)),[4,1])*0.1);
        gSimN = gSim./repmat(sum(gSim,1),[4,1]);
        
        simMatris = gSimN.*(distVals.^distPower);
        
        simMatrisN = simMatris./(repmat(sum(simMatris,1),[4,1]));
        
        cHSPxls = agirlikMatrisi(i-1:3:i+2,j-1:3:j+2);
        
        a_cHSPxls = cHSPxls(:);
        
        newPxls = simMatrisN'*a_cHSPxls;
        agirlikMatrisi(i:i+1,j:j+1) = reshape(newPxls,[2,2]);
        
    end
end

