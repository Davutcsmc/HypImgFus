function [agirlikMatrisi] = edgePreserveGuidedInterpolationAlg1(refImg, agirlikMatrisi,distPower)

agirlikMatrisiYedek = agirlikMatrisi;

[wm,hm]=size(refImg);
if nargin < 3
    distPower = 1;
end

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

for i=5:3:wm-4
    for j=5:3:hm-4
        
        nPxls = refImg(i-1:i+1,j-1:j+1);
        cPxls = refImg(i-3:3:i+3,j-3:3:j+3);
                
        gDiff = abs( repmat(nPxls(:)',[9,1]) - repmat(cPxls(:),[1,9]) );
        
        diffMatris = gDiff.*(distVals.^4);
        
        sumX = sum(diffMatris,1);
        
        diffMatrisN = diffMatris./repmat(sumX,[9,1]);
        
        simMatrisN = (1 - diffMatrisN)./repmat(sum(1-diffMatrisN),[9,1]);
        
%         nHSPxls = agirlikMatrisi(i-1:i+1,j-1:j+1);
        cHSPxls = agirlikMatrisi(i-3:3:i+3,j-3:3:j+3);
        
%         a_nHSPxls = nHSPxls(:);
        a_cHSPxls = cHSPxls(:);
        
        newPxls = simMatrisN'*a_cHSPxls;
        
        centerPxl = agirlikMatrisi(i,j);
        agirlikMatrisi(i-1:i+1,j-1:j+1) = reshape(newPxls,[3,3]);
        agirlikMatrisi(i,j) = centerPxl;
        
    end
end

