function [msi,hsi2msiCoeff,spectResp]=funcHSI2MSI(hsi,strDataType,strDataName)

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

if (strcmp(strDataType,'MS'))
    
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
    
    hsi2msiCoeff = [gaussB' gaussG' gaussR' gaussIR'];
    % plot(spectResp,gaussB), hold on,
    % plot(spectResp,gaussG,'g')
    % plot(spectResp,gaussR,'r')
    % plot(spectResp,gaussIR,'m')
    
    msiB = 0; msiG = 0; msiR = 0; msiIR = 0;
    for i=1:1:length(spectResp)
        band  = hsi(:,:,i);
        msiB  = msiB  + gaussB(i)  * band;
        msiG  = msiG  + gaussG(i)  * band;
        msiR  = msiR  + gaussR(i)  * band;
        msiIR = msiIR + gaussIR(i) * band;
    end
    msi = cat(3,msiB,msiG,msiR,msiIR);
    clear i band msiB msiG msiR msiIR cB cG cR cIR stdVal scaler
elseif (strcmp(strDataType,'PAN'))
    coeff = 1/40;    
    hsi2msiCoeff = zeros(length(spectResp),1);
    hsi2msiCoeff(1:40) = coeff;
    
    msi = 0;
    %%%% PAN image
    for i=1:1:40
        band  = hsi(:,:,i);
        msi  = msi +  hsi2msiCoeff(i) * band;
    end  
end
end