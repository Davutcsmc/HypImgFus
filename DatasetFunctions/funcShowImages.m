function []= funcShowImages(hsiData,strDataName)

%% set band indices
switch(strDataName)
   
    case {'Salinas','salinas','SALINAS'}
        bR= 30;   bG= 20; bB =10;  
    case {'Pavia','pavia','PAVIA'}
        bR= 45;   bG= 35; bB =30;
    case {'Sentetic','sentetic','SENTETIC'}
        bR= 30;   bG= 20; bB =10;
    otherwise
        warning('please define avalid data name');
        return;
end

%% rgb - hsi
hsiRGB = cat(3,hsiData.hsi(:,:,bR),hsiData.hsi(:,:,bG),hsiData.hsi(:,:,bB));
[hsiRGBN] = funcNormalize(hsiRGB,'all'); % "all", "individual"

%% rgb - hsiOrj 
hsiOrjRGB = cat(3,hsiData.hsiOrj(:,:,bR),hsiData.hsiOrj(:,:,bG),hsiData.hsiOrj(:,:,bB));
[hsiOrjRGBN] = funcNormalize(hsiOrjRGB,'all'); % "all", "individual"

%% rgb - msi
[msiN] = funcNormalize(hsiData.msi,'all'); % "all", "individual"
msiNrgb = cat(3,msiN(:,:,3),msiN(:,:,2),msiN(:,:,1));

%% show images
figure(12), 
subplot(131), imshow(msiNrgb); title('MSI-RGB')
subplot(132), imshow(hsiOrjRGBN); title('HSI-Orj')
subplot(133), imshow(hsiRGBN); title('HSI')

end





