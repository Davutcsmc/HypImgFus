function [out,outBands] = RMSE(ref,tar)
%--------------------------------------------------------------------------
% Root mean squared error (RMSE)
%
% USAGE
%   out = RMSE(ref,tar)
%
% INPUT
%   ref : reference HS data (rows,cols,bands)
%   tar : target HS data (rows,cols,bands)
%
% OUTPUT
%   out : RMSE (scalar)
%
%--------------------------------------------------------------------------
[rows,cols,bands] = size(ref);
out = (sum(sum(sum((tar-ref).^2)))/(rows*cols*bands)).^0.5;
outBands= ones(bands,1)*NaN;
for i=1:1:bands
    outBands(i,1) = (sum(sum(sum((tar(:,:,i)-ref(:,:,i)).^2)))/(rows*cols)).^0.5;
end

