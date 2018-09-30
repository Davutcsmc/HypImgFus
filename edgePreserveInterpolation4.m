function [refImageInterpolated] = edgePreserveInterpolation4(guideImage,refImage,distPower)

[wm,hm]=size(guideImage);
if nargin < 3
    distPower = 1;
end

refImageInterpolated = refImage;
if (size(refImage,1) ~= wm || size(refImage,2) ~= hm)
    refImageInterpolated = imresize(refImage,[wm hm], 'bicubic','Antialiasing',true);
    refImageInterpolated(2:end-1,2:end-1)= 0;
    refImageInterpolated(2:3:end-1,2:3:end-1) = refImage;
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

distValsEdgeRegion = distVals./repmat(sum(distVals,1),[2,1]);

i=5;
j=5;
[b1,b2] = meshgrid(i,j+1:j+2);
[c1,c2] = meshgrid(i-3:3:i+3,j:3:j+3);
b11 = b1(:);
b22 = b2(:);
c11 = c1(:);
c22 = c2(:);
distVals = [];
for ii=1:1:length(b11)
    distVals(:,ii) = sqrt( (b11(ii) - c11).^2 + (b22(ii)-c22).^2 );
end
distValsMiddleRegion = distVals./repmat(sum(distVals,1),[6,1]);

for i=2:3:wm-1
    for j=2:3:hm-2
        
        %%% kenarbölgesi
        if ((i-3)<=0) || ((i+3)>wm)
            
            simMatris = distValsEdgeRegion.^distPower;
            
            simMatrisN = simMatris./(repmat(sum(simMatris,1),[2,1]));
            
            cHSPxls = refImageInterpolated(i,j:3:j+3);
            
            a_cHSPxls = cHSPxls(:);
            
            newPxls = simMatrisN'*a_cHSPxls;
            
            if sum(isnan(newPxls))>=1
                'nan value detected';
            end
            refImageInterpolated(i,j+1:j+2) = newPxls';
        else %%% orta alan bölgesi
            
            simMatris =distValsMiddleRegion.^distPower;
            
            simMatrisN = simMatris./(repmat(sum(simMatris,1),[6,1]));
            
            cHSPxls = refImageInterpolated(i-3:3:i+3,j:3:j+3);
            
            a_cHSPxls = cHSPxls(:);
            
            newPxls = simMatrisN'*a_cHSPxls;
            
            if sum(isnan(newPxls))>=1
                'nan value detected';
            end
            
            refImageInterpolated(i,j+1:j+2) = newPxls';
        end
    end
end

clear i j b1 b2 c1 c2 b11 b22 c11 c22 ii distVals nPxls cPxls gDiff ...
    gDiffN gSim gSimN simMatris simMatrisN cHSPxls a_cHSPxls newPxls ...
    centerPxl distValsEdgeRegion distValsMiddleRegion

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

distValsEdgeRegion = distVals./repmat(sum(distVals,1),[2,1]);

i=5;
j=5;
[b1,b2] = meshgrid(i+1:i+2,j);
[c1,c2] = meshgrid(i:3:i+3,j-3:3:j+3);
b11 = b1(:);
b22 = b2(:);
c11 = c1(:);
c22 = c2(:);
distVals = [];
for ii=1:1:length(b11)
    distVals(:,ii) = sqrt( (b11(ii) - c11).^2 + (b22(ii)-c22).^2 );
end
distValsMiddleRegion = distVals./repmat(sum(distVals,1),[6,1]);


for i=2:3:wm-2
    for j=2:3:hm-1
        
        %%% kenarbölgesi
        if ((j-3)<=0) || ((j+3)>hm)
            
            simMatris = distValsEdgeRegion.^distPower;
            
            simMatrisN = simMatris./(repmat(sum(simMatris,1),[2,1]));
            
            cHSPxls = refImageInterpolated(i:3:i+3,j);
            
            a_cHSPxls = cHSPxls(:);
            
            newPxls = simMatrisN'*a_cHSPxls;
            
            if sum(isnan(newPxls))>=1
                'nan value detected';
            end
            refImageInterpolated(i+1:i+2,j) = newPxls';
            
        else %%%% orta alan bölgesi
            
            simMatris = distValsMiddleRegion.^distPower;
            
            simMatrisN = simMatris./(repmat(sum(simMatris,1),[6,1]));
            
            cHSPxls = refImageInterpolated(i:3:i+3,j-3:3:j+3);
            
            a_cHSPxls = cHSPxls(:);
            
            newPxls = simMatrisN'*a_cHSPxls;
            
            if sum(isnan(newPxls))>=1
                'nan value detected';
            end
            
            refImageInterpolated(i+1:i+2,j) = newPxls';
        end
    end
end

clear i j b1 b2 c1 c2 b11 b22 c11 c22 ii distVals nPxls cPxls gDiff ...
    gDiffN gSim gSimN simMatris simMatrisN cHSPxls a_cHSPxls newPxls ...
    centerPxl distValsEdgeRegion distValsMiddleRegion


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
        
        simMatris = distVals.^distPower;
        
        simMatrisN = simMatris./(repmat(sum(simMatris,1),[4,1]));
        
        cHSPxls = refImageInterpolated(i-1:3:i+2,j-1:3:j+2);
        
        a_cHSPxls = cHSPxls(:);
        
        newPxls = simMatrisN'*a_cHSPxls;
        
        if sum(isnan(newPxls))>=1
            'nan value detected';
        end
        
        refImageInterpolated(i:i+1,j:j+1) = reshape(newPxls,[2,2]);
        
    end
end

