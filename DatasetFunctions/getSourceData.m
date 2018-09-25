function [dataSet]=getSourceData(strDataName,strDataType,ratio)

if (nargin==0)
    disp('please set a data name as input: Salinas, Pavia or Sentetic ');
end
    
dataSet(1).hsi      = 0;
dataSet(1).hsiOrj   = 0;
dataSet(1).hsiSNR   = NaN;
dataSet(1).W        = 0;
dataSet(1).H        = 0;
dataSet(1).Dim      = 0;
dataSet(1).GT       = 0;
dataSet(1).msi      = 0;
dataSet(1).LabelMap = 0;
dataSet(1).LabelFrac= 0;
dataSet(1).Signs    = 0;
dataSet(1).SNRRatio = 0;
dataSet(1).hsi2msiCoeff = 0;
dataSet(1).wavelength = 0;

switch(lower(strDataName))
   
    case {'salinas'}
        [dataSet(1).hsi,dataSet(1).hsiOrj,dataSet(1).W,dataSet(1).H,dataSet(1).DimH, dataSet(1).msi, dataSet(1).GT, dataSet(1).hsi2msiCoeff,dataSet(1).wavelength]=...
            getSalinas(strDataType,ratio);        
    case {'pavia'}
        [dataSet(1).hsi,dataSet(1).hsiOrj,dataSet(1).W,dataSet(1).H,dataSet(1).DimH, dataSet(1).msi, dataSet(1).hsi2msiCoeff,dataSet(1).wavelength]=...
            getPavia(strDataType,ratio);
    case {'sentetic'}
        [dataSet(1).hsi,dataSet(1).hsiSNR,dataSet(1).hsiOrj,dataSet(1).W,dataSet(1).H,dataSet(1).DimH, dataSet(1).msi,dataSet(1).GT, ...
            dataSet(1).LabelMap, dataSet(1).LabelFrac, dataSet(1).Signs, dataSet(1).SNRRatio, dataSet(1).hsi2msiCoeff, dataSet(1).wavelength ]=...
            getSentetic(40,strDataType,ratio);            
    otherwise
        warning('please define avalid data name');
        return;
end

end






