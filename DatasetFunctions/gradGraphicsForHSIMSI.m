function [] = gradGraphicsForHSIMSI(hsi,msi)

figure(21),

subplot(4,4,1); imagesc(squeeze(msi(:,:,1))); title('MSI');
for i=1:1:14
    subplot(4,4,i+1); imagesc(squeeze(hsi(:,:,i+4))); title(strcat('HSI - ', num2str(4+i)));
end

[gradMVec] = gradVec(msi);

[gradHVecB] = gradVec(hsi(:,:,5:18));
[gradHVecG] = gradVec(hsi(:,:,10:23));
[gradHVecR] = gradVec(hsi(:,:,23:36));
[gradHVecIR] = gradVec(hsi(:,:,34:46));

figure(22), 

j= 1;
for i=1:2:13
    subplot(2,4,j); hold on, grid on;
    plot3(gradMVec(:,1),gradMVec(:,2),gradMVec(:,3),'*','LineWidth',2,'Color','b');
    plot3(gradHVecB(:,i),gradHVecG(:,i),gradHVecR(:,i),'*','LineWidth',2,'Color','r');
    plot3(gradHVecB(:,i+1),gradHVecG(:,i+1),gradHVecR(:,i+1),'*','LineWidth',2,'Color','g');
    title(strcat('MSI HSI B:',num2str(i+4),' - ',num2str(i+5), ', G:',num2str(i+9),' - ',num2str(i+10), ', R:',num2str(i+22),' - ',num2str(i+23)));
    xlabel('Blue Band'); ylabel('Green Band'); zlabel('Red Band');
    j=j+1;
end

end

function [gradMVec] = gradVec(B)

gradMx = B(1:end-1,:,:) - B(2:end,:,:);
gradMy = B(:,1:end-1,:) - B(:,2:end,:);
gradMxy= sqrt(gradMx(:,1:end-1,:).^2 + gradMy(1:end-1,:,:).^2);

gradMxyN = funcNormalize(gradMxy);

gradMVec = reshape(gradMxyN,size(gradMxyN,1)*size(gradMxyN,2),size(gradMxyN,3));

end