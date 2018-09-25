function [imgSNR]=funcAddSNR(img, SNRRatio)

if (SNRRatio ~= 0)
    imgSNR = zeros(size(img));
    for itrDim = 1:1:size(img,3)
        band = img(:,:,itrDim);
        bandSNR = awgn(band,SNRRatio);
        imgSNR(:,:,itrDim) = bandSNR;
    end  
else
    imgSNR = img;
end

end