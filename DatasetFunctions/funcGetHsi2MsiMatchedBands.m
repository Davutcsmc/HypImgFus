function [hsi2msiBands]=funcGetHsi2MsiMatchedBands(strDataName)

hsi2msiBands = 0;
switch(strDataName)
   
    case {'Salinas','salinas','SALINAS'}
        [hsi2msiBands]=getSalinasIkonosBands();
    case {'Pavia','pavia','PAVIA'}
        [hsi2msiBands]=getPaviaIkonosBands();
    case {'Sentetic','sentetic','SENTETIC'}
        [hsi2msiBands]=getSenteticIkonosBands();
    otherwise
        warning('please define avalid data name');
        return;
end
end

function [hsi2msiBands]=getSalinasIkonosBands()

load('C:\Users\davut.cesmeci\Google Drive\2016_fusion\FusionWorks\codes\data\Salinas\Salinas_corrected');

hsi = salinas_corrected(20:103,65:148,:);
[rH,cH,dimH] = size(hsi);

%% resampled hsi
resR = 0.25;
[hsiResampled] = funcResampleImg(hsi,resR);

%% create msi
[msi] = funcHSI2MSI(hsi,'salinas');
[hsi2msiBands]=funcHSI2MSIMatchingBands(hsi,'Salinas');

end

function [hsi2msiBands]=getPaviaIkonosBands()
load('C:\Users\davut.cesmeci\Google Drive\2016_fusion\FusionWorks\codes\data\Pavia\PaviaU');

% a = paviaU(163:242,101:180,50);
% imshow(a,[]);
hsi = paviaU(162:245,100:183,:);
[rH,cH,dimH] = size(hsi);

%% resampled hsi
resR = 0.25;
[hsiResampled] = funcResampleImg(hsi,resR);

%% create hsi2msi bands
[hsi2msiBands]=funcHSI2MSIMatchingBands(hsi,'Pavia');
end


function [hsi2msiBands]=getSenteticIkonosBands()
if (nargin== 0)
    SNRRatio = 40;
end

dataDir = 'C:\Users\davut.cesmeci\Google Drive\2016_fusion\FusionWorks\codes\data\senteticData\';
load(strcat(dataDir,'signsUSGS'));
load(strcat(dataDir,'spectResp'));

% % ----------------------------------------------------------
% % 1. durum :(hepsi farklý)--  2. durum (3-4 MSI da benzer)--
% 1-alunite                 -- 1- alunite                   --
% 2-muscovite               -- 2- muscovite                 --
% 3-kaolinite               -- 3- kaolinite                 --
% 4-dumortierite            -- 4- halloysite                --
% 5-brick                   -- 5- brick                     --
% % -----------------------------------------------------------
durum = 1; %%%% durum = 2;
if isequal(durum,1) %%%% hepsi farklý
    inds = [ 1 8 6 3 2 ];
    signNames   = signatureLabels(inds);
    signVectors = signsUSGS(inds,:);
elseif isequal(durum,2) %%%% 3-4 MSI da benzer HSI da farklý
    inds = [ 1 8 6 4 2 ];
    signNames   = signatureLabels(inds);
    signVectors = signsUSGS(inds,:);
else
    error('Davut: define a valid ''durum'' value...');
end

Signs.SignNames = signNames;
Signs.SignVec = signVectors;
Signs.SpectResp = spectResp;

dimH = size(signVectors,2);
clear signsUSGS signatureLabels inds
maprgb = double( imread(strcat(dataDir,'groundTruthMapLarge.tif')));
maprgb = imresize(maprgb,0.5,'nearest');
gt = maprgb;
% maprgb = maprgb(:,2:401,:);
[row,column,dimMap] = size(maprgb);
maprgbV = reshape(maprgb,row*column,dimMap);

%% repair ground-Truth

imBrick = zeros(24,24,3);
for ii=1:1:8
    for jj=9:1:24
        imBrick(ii,jj,:) = [0 255 255]';
    end
end
for ii=9:1:24
    for jj=1:1:24
        imBrick(ii,jj,:) = [0 255 255]';
    end
end
for ii=1:1:24
    for jj=1:1:24
        if isequal( squeeze( imBrick(ii,jj,:)),[0 255 255]' )
            maprgb(27+ii,135+jj,:) = imBrick(ii,jj,:);
        end
    end
end,clear row column dimMap ii jj imBrick
% dataSize = '8kat';
dataSize = '4kat';
if strcmp(dataSize,'4kat')
    maprgb = maprgb(21:end-20,21:end-20,:);
    [row,column,dimMap] = size(maprgb);
else
    error('Davut: define sparial res ratio as 8 ...')
end
imagesc(uint8(maprgb));
maprgbV = reshape(maprgb,row*column,dimMap);

%% create hsi
[uniqueLabels] = unique(maprgbV,'rows');
if strcmp(dataSize,'kucuk')
    imzaSirasi = [1 4 2 3];
elseif strcmp(dataSize,'cokKucuk')
    imzaSirasi = [1 2 ];
else
    imzaSirasi = [1 4 5 2 3];
end
hsiV = ones(size(maprgbV,1),dimH)*NaN;
labelMapV = ones(size(maprgbV,1),1)*NaN;
labelFracV= zeros(size(maprgbV,1),size(uniqueLabels,1));
for i=1:1:size(uniqueLabels,1)
    
    vecSearch= double( uniqueLabels(i,:) );
    [vecPos] = strmatch(vecSearch,maprgbV);
    
    hsiV(vecPos,:) = repmat(signVectors(imzaSirasi(i),:),length(vecPos),1);
    labelMapV(vecPos,:) = repmat(imzaSirasi(i),length(vecPos),1);
    labelFracV(vecPos,imzaSirasi(i)) = repmat(1,length(vecPos),1);
end, clear i vecSearch vecPos
hsi = reshape(hsiV,row,column,dimH);
labelMap = reshape(labelMapV,row,column);
labelFrac = reshape(labelFracV,row,column,size(uniqueLabels,1));

clear hsiV i vecSearch vecPos imzaSirasi uniqueLabels labelMapV maprgbV labelFracV
hsi = hsi / 1e4;
[rH, cH, dimH] = size(hsi);
close all;

%% resampled hsi
resR  = 0.25;
[hsiResample]=funcResampleImg(hsi,resR);
[hsiResampleSNR] = funcAddSNR(hsiResample,SNRRatio);

%% create msi
[hsi2msiBands]=funcHSI2MSIMatchingBands(hsiResampleSNR,'Sentetic');

end

function [hsi2msiBands]=funcHSI2MSIMatchingBands(hsi,strDataName)

switch(strDataName)
   
    case {'Salinas','salinas','SALINAS'}
        spectResp = 380:10:2410;
    case {'Pavia','pavia','PAVIA'}
        spectResp = 434:4:844;
    case {'Sentetic','sentetic','SENTETIC'}
        dataDir = 'C:\Users\davut.cesmeci\Google Drive\2016_fusion\FusionWorks\codes\data\senteticData\';
        load(strcat(dataDir,'spectResp'));
    otherwise
        warning('please define avalid data name');
        return;
end

%%%% IKONOS wavelengths
cB = 500; cG = 550; cR = 670; cIR = 770; stdVal = 20;
scaler = 1/(sqrt(2*pi)*stdVal);
gaussB = scaler * exp( -1/2 * ((spectResp - cB) /(stdVal)).^2 );
gaussG = scaler * exp( -1/2 * ((spectResp - cG) /(stdVal)).^2 );
gaussR = scaler * exp( -1/2 * ((spectResp - cR) /(stdVal)).^2 );
gaussIR= scaler * exp( -1/2 * ((spectResp - cIR)/(stdVal)).^2 );

gaussB  = gaussB  / sum(gaussB);
gaussG  = gaussG  / sum(gaussG);
gaussR  = gaussR  / sum(gaussR);
gaussIR = gaussIR / sum(gaussIR);

% plot(spectResp,gaussB), hold on,
% plot(spectResp,gaussG,'g')
% plot(spectResp,gaussR,'r')
% plot(spectResp,gaussIR,'m')
[row,col,dim] = size(hsi);
msiB = [];
msiG = [];
msiR = [];
msiIR = [];
for i=1:1:length(spectResp)
    band  = hsi(:,:,i);
    if gaussB(i) > 0.001
        msiB  = cat(3,msiB,band);
    end
    if gaussR(i) > 0.001
        msiG  = cat(3,msiG,band);
    end
    if gaussG(i) > 0.001
        msiR  = cat(3,msiR,band);
    end
    if gaussIR(i) > 0.001
        msiIR  = cat(3,msiIR,band);
    end
end
hsi2msiBands.B = msiB;
hsi2msiBands.G = msiG;
hsi2msiBands.R = msiR;
hsi2msiBands.IR = msiIR;
clear i band msiB msiG msiR msiIR cB cG cR cIR stdVal scaler

end