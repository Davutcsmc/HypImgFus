function [ I_HS_seq ] = GF_GRatioSpatial( I_HS,I_PAN, dataset )
%GFSEQUENTÝAL Summary of this function goes here
%   Detailed explanation goes here

    [~,maxInds] = max(dataset.overlap);
    
    %%%% HS
    I_HS( I_HS<0 ) = 0;
    minV = min(I_HS(:)); maxV = max(I_HS(:));
    [wh,hh,dh] = size(I_HS);
    sRatio = I_HS ./( repmat(I_HS(:,:,maxInds(2)),[1,1,dh]) +1e-2);
%     vsRatio = reshape(sRatio,[wh*hh,dh-1])';
    %     figure(11), plot(vsRatio);
%     vsMedRatio = medfilt1(vsRatio,3);
    
%     sMedRatio = reshape(vsMedRatio',[wh,hh,dh-1]);
    
    %%%% HS bicubic    
    [wm,hm,dm] = size(I_PAN);
    I_HSbicubic = imresize(I_HS,[wm hm], 'bicubic','Antialiasing',true);
%     sCubicRatio = I_HSbicubic(:,:,1:end-1) ./( I_HSbicubic(:,:,2:end) +1e-2);
%     vsCubicRatio = reshape(sCubicRatio,[wm*hm,dh-1])';
%     %     figure(11), plot(vsRatiobicubic);
%     vsMedCubicRatio = medfilt1(vsCubicRatio,3);
    
    rows1 = 2:3:wm; cols1 = 2:3:hm; Mbicubic = zeros(wm,hm);
    Mbicubic(rows1, cols1) = 1;
    clear rows1 cols1;
    
    %% plot ratio between bands
    %     for i=1:1:w*h
    %         figure(1), clf, hold on,
    %         plot(vsRatio(:,i),'LineWidth',2,'Color','b');
    %         plot(vsMedRatio(:,i),'LineWidth',2,'Color','r');
    %         legend('orj','median')
    %         drawnow;
    %         k=waitforbuttonpress;
    %     end
    %% end of plotting
    %%     
    r=8;
    B_GF = guidedfilter(I_PAN(:,:,1),I_HSbicubic(:,:,maxInds(1)),r,eps); % Fusion of HSI and pan image
    G_GF = guidedfilter(I_PAN(:,:,2),I_HSbicubic(:,:,maxInds(2)),r,eps); % Fusion of HSI and pan image
    R_GF = guidedfilter(I_PAN(:,:,3),I_HSbicubic(:,:,maxInds(3)),r,eps); % Fusion of HSI and pan image
    IR_GF= guidedfilter(I_PAN(:,:,4),I_HSbicubic(:,:,maxInds(4)),r,eps); % Fusion of HSI and pan image
    
    %% imshow images
%     h1 = figure(1);
%     subplot(131), imshow(I_HSbicubic(:,:,maxInds(1)),[]);
%     subplot(132), imshow(I_PAN(:,:,1),[]);
%     subplot(133), imshow(B_GF,[]);
%     
%     h2 = figure(2);
%     subplot(131), imshow(I_HSbicubic(:,:,maxInds(2)),[]);
%     subplot(132), imshow(I_PAN(:,:,2),[]);
%     subplot(133), imshow(G_GF,[]);
%     
%     h3 = figure(3);
%     subplot(131), imshow(I_HSbicubic(:,:,maxInds(3)),[]);
%     subplot(132), imshow(I_PAN(:,:,3),[]);
%     subplot(133), imshow(R_GF,[]);
%     
%     h4 = figure(4);
%     subplot(131), imshow(I_HSbicubic(:,:,maxInds(4)),[]);
%     subplot(132), imshow(I_PAN(:,:,4),[]);
%     subplot(133), imshow(IR_GF,[]);
    %%  end of showing images
    %% 
    I_HS_MaskB = I_HSbicubic(:,:,maxInds(1)).*Mbicubic;
    I_HS_MaskG = I_HSbicubic(:,:,maxInds(2)).*Mbicubic;
    I_HS_MaskR = I_HSbicubic(:,:,maxInds(3)).*Mbicubic;
    I_HS_MaskIR= I_HSbicubic(:,:,maxInds(4)).*Mbicubic;
    %%%% residuals    
    I_HS_MaskBR = I_HS_MaskB - B_GF.*Mbicubic;
    I_HS_MaskGR = I_HS_MaskG - G_GF.*Mbicubic;
    I_HS_MaskRR = I_HS_MaskR - R_GF.*Mbicubic;
    I_HS_MaskIRR= I_HS_MaskIR- IR_GF.*Mbicubic;
       
    %%%% residual aradeðerleme
    rows1 = 2:3:wm; cols1 = 2:3:hm;
    I_HS_MaskBRorj = I_HS_MaskBR(rows1,cols1); 
    I_HS_MaskGRorj = I_HS_MaskGR(rows1,cols1); 
    I_HS_MaskRRorj = I_HS_MaskRR(rows1,cols1); 
    I_HS_MaskIRRorj = I_HS_MaskIRR(rows1,cols1); 
    
    I_HS_BRcubic = imresize(I_HS_MaskBRorj,[wm hm], 'bicubic','Antialiasing',true);
    I_HS_GRcubic = imresize(I_HS_MaskGRorj,[wm hm], 'bicubic','Antialiasing',true);
    I_HS_RRcubic = imresize(I_HS_MaskRRorj,[wm hm], 'bicubic','Antialiasing',true);
    I_HS_IRRcubic = imresize(I_HS_MaskIRRorj,[wm hm], 'bicubic','Antialiasing',true);
    clear rows1 cols1 I_HS_MaskBRorj I_HS_MaskGRorj I_HS_MaskRRorj I_HS_MaskIRRorj 
    
    %%%% estimations
%     I_HS_MaskBE = B_GF.*Mbicubic + I_HS_MaskBR;
%     I_HS_MaskGE = G_GF.*Mbicubic + I_HS_MaskGR;
%     I_HS_MaskRE = R_GF.*Mbicubic + I_HS_MaskRR;
%     I_HS_MaskIRE= IR_GF.*Mbicubic + I_HS_MaskIRR;
    
%     I_HS_BE = B_GF + I_HS_BRcubic;
%     I_HS_GE = G_GF + I_HS_GRcubic;
%     I_HS_RE = R_GF + I_HS_RRcubic;
%     I_HS_IRE= IR_GF + I_HS_IRRcubic;

    I_HS_BE = B_GF;
    I_HS_GE = G_GF;
    I_HS_RE = R_GF;
    I_HS_IRE= IR_GF;
    
    %% imshow images
%     figure(1);
%     subplot(141), imshow( I_REF(:,:,maxInds(1)),[] );
%     subplot(142), imshow( I_HSbicubic(:,:,maxInds(1)),[] );
%     subplot(143), imshow( B_GF,[] );
%     subplot(144), imshow( I_HS_BE,[] );
%     
%     figure(2);
%     subplot(141), imshow( I_REF(:,:,maxInds(2)),[] );
%     subplot(142), imshow( I_HSbicubic(:,:,maxInds(2)),[] );
%     subplot(143), imshow( G_GF,[] );
%     subplot(144), imshow( I_HS_GE,[] );
%     
%     figure(3);
%     subplot(141), imshow( I_REF(:,:,maxInds(3)),[] );
%     subplot(142), imshow( I_HSbicubic(:,:,maxInds(3)),[] );
%     subplot(143), imshow( R_GF,[] );
%     subplot(144), imshow( I_HS_RE,[] );
%     
%     figure(4);
%     subplot(141), imshow( I_REF(:,:,maxInds(4)),[] );
%     subplot(142), imshow( I_HSbicubic(:,:,maxInds(4)),[] );
%     subplot(143), imshow( IR_GF,[] );
%     subplot(144), imshow( I_HS_IRE,[] );
    %% end of show images
    %%    
    %%%% ardýþýk oranlar üzerinden B seçilerek sonuç çýkarýlýyor    
    sRatio_cubic = zeros(wm,hm,size(sRatio,3));
    for i=1:1:size(sRatio,3)
        band = sRatio(:,:,i);
        sRatio_cubic(:,:,i) = imresize(band,[wm hm], 'bicubic','Antialiasing',true);
        
    end
    
    I_HS_seq =zeros(wm,hm,dh);
    I_HS_seq(:,:,maxInds(1)) = I_HS_GE;
    
    bantInit = I_HS_GE;
    for i=1:1:size(sRatio,3)
        I_HS_seq(:,:,i) = bantInit.*sRatio_cubic(:,:,i);
    end
    
    %% imshow bands and plot signatures
%     
%     figure(1),
%     for i=1:1:dh
%         subplot(121), imshow(I_REF(:,:,i),[]); title(num2str(i));
%         subplot(122), imshow(I_HS_seq(:,:,i),[]); truesize([400,400]); title(num2str(i));        
%         drawnow;
%         pause(0.5);
%     end
%     
%      h1 = figure(1); imshow(I_HS_seq(:,:,13),[]); truesize([400,400]);
%      h2 = figure(2);
%      figure(h1);
%      for i=1:1:6
%          pause(0.1);
%          figure(h1),
%          [ypos,xpos] = ginput(1);
%          ypos= round(ypos); xpos= round(xpos);
%          
%          pause(0.1);
%          figure(h2), clf, hold on,
%          plot(squeeze(dataset.I_REF(xpos,ypos,:)),'b','LineWidth',2);
%          plot(squeeze(I_HS_seq(xpos,ypos,:)),'r','LineWidth',2);
%          
%      end
    %% end of showing and plotting
    %%
    

end

