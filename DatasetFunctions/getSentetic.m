function [hsiSNRResample,hsiSNR,hsi, rH,cH,dimH, msiSNR, gt, labelMap, labelFrac, Signs, SNRRatio, hsi2msiCoeff,wavelength]=getSentetic(SNRRatio,strDataType,resR)

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
    maprgb = maprgb(11:end-10,11:end-10,:);
    [row,column,dimMap] = size(maprgb);
    gt = gt(11:end-10,11:end-10,:);
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
[hsiSNR] = funcAddSNR(hsi,SNRRatio);
% resR  = 0.5;
[hsiSNRResample]=funcResampleImg(hsiSNR,resR);

%% create msi
[msi, hsi2msiCoeff, wavelength] = funcHSI2MSI(hsi,strDataType,'sentetic');
[msiSNR] = funcAddSNR(msi,SNRRatio);
% msirgb = cat(3,msi(:,:,3),msi(:,:,2),msi(:,:,1));
% imshow(msirgb,[]);


