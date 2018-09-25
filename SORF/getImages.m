function [dataset] = getImages(strDataName,strDataType,resR)

if nargin < 2
    strDataType = 'MS';
end

I_REF = NaN;

p = fileparts(mfilename('fullpath'));

PRECISION    = 'double';
OFFSET       = 0 ;
INTERLEAVE   = 'bsq';
BYTEORDER    = 'ieee-le';

if (strcmp(strDataName,'data1')) %% data1
    FILENAME_REF = [p '\REF']; % where you put the data
    SIZE_REF     = [395,185,176];
    I_REF        = multibandread(FILENAME_REF, SIZE_REF, PRECISION, OFFSET, INTERLEAVE, BYTEORDER);
    
elseif (strcmp(strDataName,'Sentetic') || ...
        strcmp(strDataName,'Salinas') || ...
        strcmp(strDataName,'Pavia'))
    
    [hsiData, I_REF] = getMyData(strDataName,strDataType,resR);
    
else
    error('Davut : Define a valid "DataName"');
end

if (isnan(I_REF))
    error('LocalError : I_REF is NaN');
end

ratio = 1/resR;
% ratio = 3;
% ratio = 5;
% overlap = 1:41; % commun bands (or spectral domain) between I_PAN and I_HS
overlap = hsiData.hsi2msiCoeff;
size_kernel=[ratio ratio];
sig = (1/(2*(2.7725887)/ratio^2))^0.5;
% start_pos(1)=1; % The starting point of downsampling
% start_pos(2)=1; % The starting point of downsampling
if mod(ratio,2)
    start_pos(1)= floor(ratio/2)+1; % The starting point of downsampling
    start_pos(2)= floor(ratio/2)+1; % The starting point of downsampling
else
    start_pos(1)= 1; % The starting point of downsampling
    start_pos(2)= 1;
end

I_MS = hsiData.msi;
I_PAN = mean(hsiData.msi,3);

[coeffM,paydaNorm]=funcCoeffMat(resR);
KerBlu = coeffM./repmat(paydaNorm,[1/resR 1/resR]);
% [I_HS,KerBlu]=conv_downsample(I_REF,ratio,size_kernel,sig,start_pos);
I_HS = hsiData.hsi;

dataset.ratio = ratio;
dataset.overlap = overlap;
dataset.size_kernel = size_kernel;
dataset.sig = sig;
dataset.start_pos = start_pos;
dataset.start_pos = start_pos;
dataset.I_REF = I_REF;
dataset.I_MS = I_MS;
dataset.I_PAN = I_PAN;
dataset.I_HS = I_HS;
dataset.KerBlu = KerBlu;
dataset.INTERLEAVE = INTERLEAVE;
dataset.wavelength = hsiData.wavelength;

end

function [hsiData, I_REF] = getMyData(strDataName,strDataType,ratio)
addpath('..\dataFunctions');

hsiData = getSourceData(strDataName,strDataType,ratio);

if (~isnan(hsiData.hsiSNR))
    I_REF = hsiData.hsiSNR;
else
    I_REF = hsiData.hsiOrj;
end
end

