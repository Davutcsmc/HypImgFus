function [msiN]=funcNormalize(msi,normType)

if (nargin==0)
    disp('please set a data name as input: Salinas, Pavia or Sentetic ');
elseif (nargin==1)
    normType = 'all';
end

msiN = zeros(size(msi));
switch (normType)    
    case {'all','All','ALL'}
        msiN = (( msi - min(msi(:)))./( max(msi(:)) - min(msi(:)) ));
    case{'individual','Individual','INDIVIDUAL'}        
        for i=1:1:size(msi,3)
            msiN(:,:,i) = (( msi(:,:,i) - min(min(min(msi(:,:,i))))))./( max(max(max(msi(:,:,i)))) - min(min(min(msi(:,:,i)))) );
        end        
    otherwise 
        warning('please define a valid normType');
        return;
end
