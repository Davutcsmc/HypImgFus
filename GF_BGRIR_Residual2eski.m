function [ I_HS_BGRIR ] = GF_BGRIR_Residual2eski( I_HS,I_PAN, dataset )
%GFSEQUENTÝAL Summary of this function goes here
%   Detailed explanation goes here

    [~,maxInds] = max(dataset.overlap);
    
    %%%% HS
    I_HS( I_HS<0 ) = 0;
    minV = min(I_HS(:)); maxV = max(I_HS(:));
    [wh,hh,dh] = size(I_HS);
    sRatioB = I_HS ./( repmat(I_HS(:,:,maxInds(1)),[1,1,dh]) +1e-2);
    sRatioG = I_HS ./( repmat(I_HS(:,:,maxInds(2)),[1,1,dh]) +1e-2);
    sRatioR = I_HS ./( repmat(I_HS(:,:,maxInds(3)),[1,1,dh]) +1e-2);
    sRatioIR= I_HS ./( repmat(I_HS(:,:,maxInds(4)),[1,1,dh]) +1e-2);
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
    
    %%%% local clustering
    BGR_GF = cat(3,B_GF,G_GF,R_GF,IR_GF);
    labels = 1:wh*hh;
    labels = reshape(labels,[wh,hh]);
    labelmap = zeros(wm,hm);
    rowsl = 2:3:wm; colsl = 2:3:hm;
    labelmap(rowsl,colsl) = labels;
    
    labelmap2 = zeros(wm,hm);
    for i=1:1:wh
        if i==1
            x1 = 1;
            x2 = rowsl(i);
        elseif i==wh
            x1 = rowsl(i-1);
            x2 = wm;
        else
            x1 = rowsl(i-1);
            x2 = rowsl(i);
        end
        for j=1:1:hh
            if j==1
                y1 = 1;
                y2 = colsl(j);
            elseif j==hh
                y1 = colsl(j-1);
                y2 = hm;
            else
                y1 = colsl(j-1);
                y2 = colsl(j);
            end
            sampleArea = BGR_GF(x1:x2,y1:y2,:);
            [wL,hL,dL]= size(sampleArea);            
            localLabelArea  = labelmap(x1:x2,y1:y2,:);
            localLabels = localLabelArea(localLabelArea~=0);
            
            vSampleArea = reshape(sampleArea,[wL*hL,dL])';
            distMap = zeros(length(localLabels),wL*hL);
            for k=1:1:length(localLabels)
                [str,stn] = find(localLabelArea==localLabels(k));
                center = squeeze(sampleArea(str,stn,:));
                distMap(k,:) = sum((vSampleArea - repmat(center,[1,size(vSampleArea,2)])).^2);
            end
            [val,inds] = min(distMap,[],1);
            inds2 = zeros(1,length(inds));
            for k=1:1:length(localLabels)
                inds2(inds==k)=localLabels(k);
            end
            m_inds = reshape(inds2,[wL,hL]);
            labelmap2(x1:x2,y1:y2) = m_inds;
        end
    end
    
    %%%%show clustering result
    figure(11111), 
    for i=1:2:wh
        if i==1
            x1 = 1;
            x2 = rowsl(i);
        elseif i==wh
            x1 = rowsl(i-1);
            x2 = wm;
        else
            x1 = rowsl(i-1);
            x2 = rowsl(i+1);
        end
        for j=1:2:hh
            if j==1
                y1 = 1;
                y2 = colsl(j);
            elseif j==hh
                y1 = colsl(j-1);
                y2 = hm;
            else
                y1 = colsl(j-1);
                y2 = colsl(j+1);
            end                   
            localLabelArea  = labelmap2(x1:x2,y1:y2,:);
            uniqLabels = unique(localLabelArea);
            
            colorB = cat(3,B_GF,B_GF,B_GF);
            maxLocalV=max(colorB(:));
            colors = [ 1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1;0 1 1; 0.7 0.5 0.1; 0.1 0.5 0.7; 0.5 0.7 0.1]*maxLocalV;
            colormap = zeros([size(localLabelArea),3]);
            for k=1:1:length(uniqLabels)
                [str,stn] = find(localLabelArea==uniqLabels(k));
                for ii=1:1:length(str)
                    if (k > size(colors,1))
                        colormap(str(ii),stn(ii),:) = [0 0 0];
                    else
                        colormap(str(ii),stn(ii),:) = colors(k,:);
                    end
                end                
            end
            
            colorB(x1:x2,y1:y2,:) = colormap;
            
            imshow(uint8(colorB./maxLocalV*255),[]); 
            truesize([500 500]);
            drawnow, 
            pause(0.1);
            waitforbuttonpress;
        end
    end
    
    
    
    
    
    
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
     
    
    
%     I_HS_BRcubic = imresize(I_HS_MaskBRorj,[wm hm], 'bicubic','Antialiasing',true);
%     I_HS_GRcubic = imresize(I_HS_MaskGRorj,[wm hm], 'bicubic','Antialiasing',true);
%     I_HS_RRcubic = imresize(I_HS_MaskRRorj,[wm hm], 'bicubic','Antialiasing',true);
%     I_HS_IRRcubic = imresize(I_HS_MaskIRRorj,[wm hm], 'bicubic','Antialiasing',true);
%     clear rows1 cols1 I_HS_MaskBRorj I_HS_MaskGRorj I_HS_MaskRRorj I_HS_MaskIRRorj 
    
%     burada merkez piksellerinden bir clustering yapýlýp baþarým iyileþtirilmesi yapýlabilir.
    
    %%%% estimations
%     I_HS_MaskBE = B_GF.*Mbicubic + I_HS_MaskBR;
%     I_HS_MaskGE = G_GF.*Mbicubic + I_HS_MaskGR;
%     I_HS_MaskRE = R_GF.*Mbicubic + I_HS_MaskRR;
%     I_HS_MaskIRE= IR_GF.*Mbicubic + I_HS_MaskIRR;
    
    I_HS_BE = B_GF + I_HS_BRcubic;
    I_HS_GE = G_GF + I_HS_GRcubic;
    I_HS_RE = R_GF + I_HS_RRcubic;
    I_HS_IRE= IR_GF + I_HS_IRRcubic;
    
    %% imshow images
%     figure(1);
%     subplot(141), imshow( dataset.I_REF(:,:,maxInds(1)),[] );
%     subplot(142), imshow( I_HSbicubic(:,:,maxInds(1)),[] );
%     subplot(143), imshow( B_GF,[] );
%     subplot(144), imshow( I_HS_BE,[] );
%     
%     figure(2);
%     subplot(141), imshow( dataset.I_REF(:,:,maxInds(2)),[] );
%     subplot(142), imshow( I_HSbicubic(:,:,maxInds(2)),[] );
%     subplot(143), imshow( G_GF,[] );
%     subplot(144), imshow( I_HS_GE,[] );
%     
%     figure(3);
%     subplot(141), imshow( dataset.I_REF(:,:,maxInds(3)),[] );
%     subplot(142), imshow( I_HSbicubic(:,:,maxInds(3)),[] );
%     subplot(143), imshow( R_GF,[] );
%     subplot(144), imshow( I_HS_RE,[] );
%     
%     figure(4);
%     subplot(141), imshow( dataset.I_REF(:,:,maxInds(4)),[] );
%     subplot(142), imshow( I_HSbicubic(:,:,maxInds(4)),[] );
%     subplot(143), imshow( IR_GF,[] );
%     subplot(144), imshow( I_HS_IRE,[] );
    %% end of show images
    %%    
    %%%% ardýþýk oranlar üzerinden B seçilerek sonuç çýkarýlýyor    
    sRatioB_cubic = zeros(wm,hm,size(sRatioB,3));
    sRatioG_cubic = zeros(wm,hm,size(sRatioG,3));
    sRatioR_cubic = zeros(wm,hm,size(sRatioR,3));
    sRatioIR_cubic= zeros(wm,hm,size(sRatioIR,3));
    for i=1:1:size(sRatioB,3)
        bandB = sRatioB(:,:,i);
        bandG = sRatioG(:,:,i);
        bandR = sRatioR(:,:,i);
        bandIR= sRatioIR(:,:,i);
        sRatioB_cubic(:,:,i) = imresize(bandB,[wm hm], 'bicubic','Antialiasing',true);
        sRatioG_cubic(:,:,i) = imresize(bandG,[wm hm], 'bicubic','Antialiasing',true);
        sRatioR_cubic(:,:,i) = imresize(bandR,[wm hm], 'bicubic','Antialiasing',true);
        sRatioIR_cubic(:,:,i)= imresize(bandIR,[wm hm], 'bicubic','Antialiasing',true);
    end
    
    I_HS_BGRIR =zeros(wm,hm,dh);
%     I_HS_BGRIR(:,:,maxInds(1)) = I_HS_BE;
    
%     %%%% weights are calculating for B, G, R, IR    
%     bandWeights = max(dataset.overlap,1e-5);
%     closenessBandindiceMask = sum(dataset.overlap,2) <= 1e-4; 
%     closenessIndices = zeros(length(closenessBandindiceMask),1);
%     weights = zeros(length(closenessBandindiceMask),size(dataset.overlap,2));
%     for i=1:1:length(closenessBandindiceMask)
%         if (closenessBandindiceMask(i)==1)
%             [~,inds] = min(abs(i-maxInds));
%             closenessIndices(i) = inds;
%             weights(i,inds) = 1;
%         else
%             weights(i,:) = bandWeights(i,:) ./ repmat(sum(bandWeights(i,:),2),[1,size(dataset.overlap,2)]);
%         end
%     end
%         
%     for i=1:1:size(sRatioB,3)
%         I_HS_BGRIR(:,:,i) = I_HS_BE.*sRatioB_cubic(:,:,i)* weights(i,1) + ...
%             I_HS_GE.*sRatioG_cubic(:,:,i)* weights(i,2) + ...
%             I_HS_RE.*sRatioR_cubic(:,:,i)* weights(i,3) + ...
%             I_HS_IRE.*sRatioIR_cubic(:,:,i)* weights(i,4);
%     end
    
I_PANds = zeros(wh,hh,dm);
for i=1:1:dm
    I_PANds(:,:,i) = imresize(I_PAN(:,:,i),[wh,hh],'bicubic','Antialiasing',true);
end

ccVals = zeros(dh,dm);
for i=1:1:dh
    bandH = I_HS(:,:,i);
    for j=1:1:dm
        bandP = I_PANds(:,:,j);
        ccVals(i,j) = corr2(bandH,bandP);
    end
end

% figure(11), hold on,
% plot(ccVals(:,1),'b','LineWidth',2);
% plot(ccVals(:,2),'g','LineWidth',2);
% plot(ccVals(:,3),'r','LineWidth',2);
% plot(ccVals(:,4),'k','LineWidth',2);

[~,ccInds] = max(ccVals,[],2);

for i=1:1:dh
    if ccInds(i) == 1
        I_HS_BGRIR(:,:,i) = I_HS_BE.*sRatioB_cubic(:,:,i);
    elseif ccInds(i) == 2
        I_HS_BGRIR(:,:,i) = I_HS_GE.*sRatioG_cubic(:,:,i);
    elseif ccInds(i) == 3
        I_HS_BGRIR(:,:,i) = I_HS_RE.*sRatioR_cubic(:,:,i);
    else
        I_HS_BGRIR(:,:,i) = I_HS_IRE.*sRatioIR_cubic(:,:,i);
    end
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

