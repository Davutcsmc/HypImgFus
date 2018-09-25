

i = 70;
bandL = I_HS(:,:,i);
bandO = I_REF(:,:,i);
bandS = I_GFPCA(:,:,i);
bandR = I_HS_mGF_Res(:,:,i);

subplot(221), imshow(imresize(bandL,2,'nearest'),[]); title(strcat('LR   : ',num2str(i)));
subplot(222), imshow(bandO,[]);title(strcat('Orj   : ',num2str(i)));
subplot(223), imshow(bandS,[]);title(strcat('GFPCA : ',num2str(i)));
subplot(224), imshow(bandR,[]);title(strcat('SORF  : ',num2str(i)));






h1 = figure(1);
for i=1:1:204
    
    clf;
    bandS = I_BayesSparse(:,:,i);
    bandR = I_HS_mGF(:,:,i);
    
    subplot(121), imshow(bandS,[]);
    subplot(122), imshow(bandR,[]);
    
    
    title(strcat('band : ',num2str(i)));
    
    waitforbuttonpress;
end


