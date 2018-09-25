

[wh,hh,dh] = size(I_HS);
[wm,hm,dm] = size(I_PAN);


I_PANds = zeros(wh,hh,dm);
for i=1:1:dm    
    I_PANds(:,:,i) = imresize(I_PAN(:,:,i),[wh,hh],'bicubic','Antialiasing',true);
end

ccVals = zeros(dh,dm);
for i=1:1:dh
    bandH = I_HS(:,:,i);
    for j=1:1:dm
        bandP = I_PANds(:,:,j);
        ccVals(i,j) = corr2(bandH,bandP);
    end
end

figure(11), hold on,
plot(ccVals(:,1),'b','LineWidth',2);
plot(ccVals(:,2),'g','LineWidth',2);
plot(ccVals(:,3),'r','LineWidth',2);
plot(ccVals(:,4),'k','LineWidth',2);

[~,inds] = max(ccVals,[],2);
