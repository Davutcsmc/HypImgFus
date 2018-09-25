function [gradMap]=funcVecGradient(MSI,gradientType,SElement,normDegree)

if nargin < 4
    normDegree = 2;
end

if nargin < 3
    SElement = [1 1 1; 1 1 1; 1 1 1];
    vecSElement = SElement(:);
end

if nargin < 2
    gradientType = 'CMG';
end

[row,col,dim] = size(MSI);
[rowE,colE] = size(SElement);
rEh = floor(rowE/2);
cEh = floor(colE/2);

gradMap = zeros(row,col);

MSIL = zeros(size(MSI,1)+2,size(MSI,2)+2,size(MSI,3));
MSIL(2:end-1,2:end-1,:) = MSI;
MSIL(1,2:end-1,:) = MSI(1,:,:);
MSIL(end,2:end-1,:) = MSI(end,:,:);
MSIL(:,1,:) = MSIL(:,2,:);
MSIL(:,end,:) = MSIL(:,end-1,:);

for i=rEh+1:1:row+rEh
    for j=cEh+1:1:col+cEh
        
        E = MSIL(i-rEh:i+rEh,j-cEh:j+cEh,:);
        vecE = reshape(E,[(2*rEh+1)*(2*cEh+1),dim]);
        
        vecES = vecE(vecSElement==1,:);
        
        if (strcmp(gradientType,'RCMG'))
            [maxValInds]=RCMG(vecES,normDegree);
        else 
            [maxValInds]=CMG(vecES,normDegree);
        end
        gradMap(i-1,j-1) = maxValInds(3);
    end
end


end

%%%%%%%  ------- CMG 
function [maxValInds]=CMG(vecES,normDegree)

pxlNbr = size(vecES,1);
maxDistInds = zeros(pxlNbr,3);
itr = 1;
for i=1:1:pxlNbr-1
    for j=i+1:1:pxlNbr
        maxDistInds(itr,1:2) = [i,j];
        maxDistInds(itr,3) = sqrt( sum(abs(( vecES(i,:) - vecES(j,:) ).^normDegree )) );
        itr=itr+1;
    end
end

[~,maxInds] = max(maxDistInds(:,3));
[maxValInds] = maxDistInds(maxInds,:);
end

%%%%%%%  ------- RCMG 
function [maxValInds2]=RCMG(vecES,normDegree)

pxlNbr = size(vecES,1);
maxDistInds = zeros(pxlNbr,3);
itr = 1;
for i=1:1:pxlNbr-1
    for j=i+1:1:pxlNbr
        maxDistInds(itr,1:2) = [i,j];
        maxDistInds(itr,3) = sqrt( sum(abs(( vecES(i,:) - vecES(j,:) ).^normDegree )) );
        itr=itr+1;
    end
end

[~,maxInds] = max(maxDistInds(:,3));
[maxValInds1] = maxDistInds(maxInds,:);
maxInds1 = maxValInds1(1:2);

[ind1] = find(maxDistInds(:,1)==maxInds1(1));
[ind2] = find(maxDistInds(:,1)==maxInds1(2));

[ind3] = find(maxDistInds(:,2)==maxInds1(1));
[ind4] = find(maxDistInds(:,2)==maxInds1(2));

maxDistInds([ind1;ind2;ind3;ind4],3) = 0;

[~,maxInds] = max(maxDistInds(:,3));
[maxValInds2] = maxDistInds(maxInds,:);
% maxInds2 = maxValInds2(1:2);

end

