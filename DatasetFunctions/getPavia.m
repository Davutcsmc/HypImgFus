function [hsiResampled,hsi,rH,cH,dimH, msi,hsi2msiCoeff,wavelength]=getPavia(strDataType,resR)

load('..\Datasets\Pavia\PaviaU');

% a = paviaU(163:242,101:180,50);
% imshow(a,[]);
hsi = paviaU(162:245,100:183,:);
hsi = ( hsi - min(hsi(:)) )/( max(hsi(:)) - min(hsi(:)) );

[rH,cH,dimH] = size(hsi);

%% resampled hsi
% resR = 0.5;
% resR = 3;
if (1/resR == 2 || 1/resR == 4 || 1/resR == 8 || 1/resR==3)
    [hsiResampled] = funcResampleImg(hsi,resR);
% elseif (resR == 1/3)
%     %     [hsiResampled] = funcResampleImg(hsi,resR);    
%     rHr = rH/3;
%     cHr = cH/3;
%     hsiResampled = zeros(rHr,cHr,dimH);
%     for itrDim = 1:1:dimH
%         band = hsi(:,:,itrDim);
%         %%%%% sub sampling - manual
%         bandr = imresize(band,[rHr,cHr],'Method','bicubic','Antialiasing',true);
%         hsiResampled(:,:,itrDim) = bandr;
%     end, clear itrDim band bandr
end

%% create msi
[msi, hsi2msiCoeff,wavelength] = funcHSI2MSI(hsi,strDataType,'Pavia');

end
