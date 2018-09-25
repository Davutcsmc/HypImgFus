function [gradXYNVec] = funcMultibandGradyen(hsi2msi)

if (isstruct(hsi2msi))
B = hsi2msi.B;
G = hsi2msi.G;
R = hsi2msi.R;
IR = hsi2msi.IR;
else
   bandNames(1) = {'B'};
   bandNames(2) = {'G'};
   bandNames(3) = {'R'};
   bandNames(4) = {'IR'};
   for i=1:1:min(size(hsi2msi,3),4)
       eval(strcat(bandNames{i},'=hsi2msi(:,:,i);'));
   end
end

[gradBxyNVec] = calculateGradyen(B);
[gradGxyNVec] = calculateGradyen(G);
[gradRxyNVec] = calculateGradyen(R);
[gradIRxyNVec] = calculateGradyen(IR);

gradXYNVec.B = gradBxyNVec;
gradXYNVec.G = gradGxyNVec;
gradXYNVec.R = gradRxyNVec;
gradXYNVec.IR = gradIRxyNVec;

end

function [gradXYNVec] = calculateGradyen(B)

gradBx = B(1:end-1,:,:) - B(2:end,:,:);
gradBy = B(:,1:end-1,:) - B(:,2:end,:);
gradBxy= sqrt(gradBx(:,1:end-1,:).^2 + gradBy(1:end-1,:,:).^2);

gradBxyN = funcNormalize(gradBxy);

gradXYNVec = reshape(gradBxyN,size(gradBxyN,1)*size(gradBxyN,2),size(gradBxyN,3));
end