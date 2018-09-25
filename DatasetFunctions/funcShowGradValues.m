function[]=funcShowGradValues(gradMsiBands,gradHsiBands,mR,mC,hR,hC,figId)

if nargin < 7
    figId = 21;
end

gradMsiBandsM.B = reshape(gradMsiBands.B,[mR-1,mC-1]);
gradMsiBandsM.G = reshape(gradMsiBands.G,[mR-1,mC-1]);
gradMsiBandsM.R = reshape(gradMsiBands.R,[mR-1,mC-1]);


for i=1:1:size(gradHsiBands.B,2)
    gradHsiBandsM.B(:,:,i) = reshape(gradHsiBands.B(:,i),[hR-1,hC-1]);
    gradHsiBandsM.G(:,:,i) = reshape(gradHsiBands.G(:,i),[hR-1,hC-1]);
    gradHsiBandsM.R(:,:,i) = reshape(gradHsiBands.R(:,i),[hR-1,hC-1]);
end

artm = 1;
figure(figId), 
for i=1:1:size(gradHsiBands.B,2)+1
    subplot(4,4,artm); hold on,
    if i==1
        imshow(cat(3,gradMsiBandsM.B,gradMsiBandsM.G,gradMsiBandsM.R));
        title(strcat('MSI'));
    else
        imshow(cat(3,gradHsiBandsM.B(:,:,i-1),gradHsiBandsM.G(:,:,i-1),gradHsiBandsM.R(:,:,i-1)));
        title(strcat('HSI ( BGR - ',num2str(i-1),' )'));
    end
    artm=artm+1;
end

end