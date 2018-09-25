function [hsiR]=funcResampleImg(hsi,resR)

[coeffM,paydaNorm]=funcCoeffMat(resR);

band = hsi(:,:,1);
hsiR = zeros(size(band,1)*resR,size(band,2)*resR,size(hsi,3));
for itrDim = 1:1:size(hsi,3)
    band = hsi(:,:,itrDim);
    bandr = zeros(size(band,1)*resR,size(band,2)*resR);
    %%%%% sub sampling - manual
    for iii=1:1:size(band,1)*resR
        iix = (iii-1)*(1/resR)+1; iiy = iii*(1/resR); 
        for jjj=1:1:size(band,2)*resR
            jjx = (jjj-1)*(1/resR)+1; jjy = jjj*(1/resR); 
            blockA = band( iix:iiy,jjx:jjy );
            bandr(iii,jjj) = sum(sum(( blockA.*coeffM ) ))./paydaNorm;
        end
    end, clear iii jjj blockA
    hsiR(:,:,itrDim) = bandr;
end, clear itrDim iii jjj iix jjx iiy jjy

% X= hsi(:,:,20);
% Y=imfilter(X, coeffM, 'circular');
% Y=imfilter(X, coeffM, 'symmetric');
% Y1=Y(1:1/resR:end, 1:1/resR:end,:);
% Y2=Y(2:1/resR:end, 2:1/resR:end,:);
% Y3 = hsiR(:,:,20);

end