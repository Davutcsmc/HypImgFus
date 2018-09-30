
clc, clearvars, close all

%% Get Dataset

selectData = 'Pavia'; % 'Pavia' or 'Salinas'

if strcmp(selectData,'Pavia')
    strDataType = 'MS'; % PAN, MS
    strDataName = 'Pavia'; % Salinas, Pavia, Sentetic, data1
elseif strcmp(selectData,'Salinas')
    strDataType = 'MS'; % PAN, MS
    strDataName = 'Salinas'; % Salinas, Pavia, Sentetic, data1
elseif strcmp(selectData,'Sentetic')
    strDataType = 'MS'; % PAN, MS
    strDataName = 'Sentetic'; % Salinas, Pavia, Sentetic, data1
end
rVal = 1/3;
[dataset] = getImages(strDataName,strDataType,rVal);

ratio       = dataset.ratio;
overlap     = dataset.overlap;
size_kernel =dataset.size_kernel;
sig         = dataset.sig;
start_pos   = dataset.start_pos;
KerBlu      = dataset.KerBlu;
INTERLEAVE  = dataset.INTERLEAVE;

I_REF       = dataset.I_REF;
I_PAN       = dataset.I_PAN;
if (strcmp(strDataType ,'MS'))
    I_PAN   = dataset.I_MS;
end
I_HS        = dataset.I_HS;

%% Hyperspectral pansharpening

L = nextpow2(max(I_HS(:)));
addpath(genpath('../QualityIndices'));
addpath(genpath('../MethodsOnLiterature'));

%%%% B G R IR weighted ratio for GF
tic
[ I_HS_mGF_Res ] = GF_BGRIR_Residual( I_HS,I_PAN, dataset );
disp(strcat('Comp. time (I_HS_mGF_Res): ',num2str(toc)));
QI_GF_mGF_Res = QualityIndices(I_HS_mGF_Res(5:end-4,5:end-4,:),I_REF(5:end-4,5:end-4,:),ratio);

%%%% B G R IR weighted ratio for GF
distPower = 1;
tic
[ I_HS_mGF_Res21 ] = GF_BGRIR_Residual2( I_HS,I_PAN, dataset, distPower );
disp(strcat('Comp. time (I_HS_mGF_Res21): ',num2str(toc)));
QI_GF_mGF_Res21 = QualityIndices(I_HS_mGF_Res21(5:end-4,5:end-4,:),I_REF(5:end-4,5:end-4,:),ratio);

%%%% B G R IR weighted ratio for GF
distPower = 1;
tic
[ I_HS_mGF_Res31 ] = GF_BGRIR_Residual3( I_HS,I_PAN, dataset, distPower );
disp(strcat('Comp. time (I_HS_mGF_Res31): ',num2str(toc)));
QI_GF_mGF_Res31 = QualityIndices(I_HS_mGF_Res31(5:end-4,5:end-4,:),I_REF(5:end-4,5:end-4,:),ratio);

% GFPCA
tic
I_GFPCA = GFPCA(I_HS,I_PAN,4, 8, 0.001^2);
disp(strcat('Comp. time (GFPCA): ',num2str(toc)));
QI_GFPCA = QualityIndices(I_GFPCA(5:end-4,5:end-4,:),I_REF(5:end-4,5:end-4,:),ratio);

%% end of pansharpining
%% plot qualities

minCC = min( [ min( [ QI_GF_mGF_Res.ccMap QI_GF_mGF_Res21.ccMap ])] );
maxCC = max( [ max( [ QI_GF_mGF_Res.ccMap QI_GF_mGF_Res21.ccMap ])] );
figure(1111),hold on,
plot(dataset.wavelength,QI_GFPCA.ccMap,'g','LineWidth',2), hold on,
plot(dataset.wavelength,QI_GF_mGF_Res.ccMap,'k','LineWidth',2),
plot(dataset.wavelength,QI_GF_mGF_Res21.ccMap,'r','LineWidth',2),
plot(dataset.wavelength,QI_GF_mGF_Res31.ccMap,'b','LineWidth',2),
legend('GFPCA','oransal GFR','oransal GFR2','oransal GFR3'), title('CC')
axis([min(dataset.wavelength) max(dataset.wavelength) 0.9 maxCC+0.01]);

minRMSE = min( [ min( [ QI_GF_mGF_Res.rmseBands; QI_GF_mGF_Res21.rmseBands ] ) ] );
maxRMSE = max( [ max( [ QI_GF_mGF_Res.rmseBands; QI_GF_mGF_Res21.rmseBands ] ) ] );

figure(11112), hold on,
plot(dataset.wavelength,QI_GFPCA.rmseBands,'g','LineWidth',2), hold on,
plot(dataset.wavelength,QI_GF_mGF_Res.rmseBands,'k','LineWidth',2),
plot(dataset.wavelength,QI_GF_mGF_Res21.rmseBands,'r','LineWidth',2),
plot(dataset.wavelength,QI_GF_mGF_Res31.rmseBands,'b','LineWidth',2),
legend('GFPCA','oransal GFR','oransal GFR2','oransal GFR3'), title('RMSE')
axis([min(dataset.wavelength) max(dataset.wavelength) minRMSE maxRMSE+0.1]);

%% end of plotting qualities
%%


