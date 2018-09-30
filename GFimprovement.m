
clc, clearvars, close all

%% Get Dataset

for iii=2:2
    
    if (iii==1)
        strDataType = 'MS'; % PAN, MS
        strDataName = 'Pavia'; % Salinas, Pavia, Sentetic, data1    
    elseif (iii==2)
        strDataType = 'MS'; % PAN, MS
        strDataName = 'Salinas'; % Salinas, Pavia, Sentetic, data1
    elseif (iii==3)
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
    
    distPower = 1;
    tic
    [ I_HS_mGF_Res41 ] = GF_BGRIR_Residual4( I_HS,I_PAN, dataset, distPower );
    disp(strcat('Comp. time (I_HS_mGF_Res41): ',num2str(toc)));
    QI_GF_mGF_Res41 = QualityIndices(I_HS_mGF_Res41(5:end-4,5:end-4,:),I_REF(5:end-4,5:end-4,:),ratio);
    
    distPower = 1;
    tic
    [ I_HS_mGF_Res61 ] = GF_BGRIR_Residual6( I_HS,I_PAN, dataset, distPower );
    disp(strcat('Comp. time (I_HS_mGF_Res61): ',num2str(toc)));
    QI_GF_mGF_Res61 = QualityIndices(I_HS_mGF_Res61(5:end-4,5:end-4,:),I_REF(5:end-4,5:end-4,:),ratio);
    
    distPower = 1;
    tic
    [ I_HS_mGF_Res71 ] = GF_BGRIR_Residual7( I_HS,I_PAN, dataset, distPower );
    disp(strcat('Comp. time (I_HS_mGF_Res71): ',num2str(toc)));
    QI_GF_mGF_Res71 = QualityIndices(I_HS_mGF_Res71(5:end-4,5:end-4,:),I_REF(5:end-4,5:end-4,:),ratio);
    
    
%     %%%% B G R IR weighted ratio for spatial GF
%     tic
%     [ I_HS_mGF ] = GF_BGRIR( I_HS,I_PAN, dataset );
%     disp(strcat('Comp. time (I_HS_mGF): ',num2str(toc)));
%     cd ..\uptodateReview\toolbox\Quality_Indices
%     QI_GF_mGF = QualityIndices(I_HS_mGF(5:end-4,5:end-4,:),I_REF(5:end-4,5:end-4,:),ratio);
%     cd ..\..\..\GF
    
    % GFPCA
    tic
    I_GFPCA = GFPCA(I_HS,I_PAN,4, 8, 0.001^2);
    disp(strcat('Comp. time (GFPCA): ',num2str(toc)));
    QI_GFPCA = QualityIndices(I_GFPCA(5:end-4,5:end-4,:),I_REF(5:end-4,5:end-4,:),ratio);

  
    % Bayesian Sparse
%     cd '..\uptodateReview\toolbox\ilerleme9\methods\BayesFusion'
%     setup;
%     tic;
%     overlapPan = 1:42;
%     [I_BayesSparse]= BayesianFusion(I_HS,mean(I_PAN,3),overlapPan,KerBlu,ratio,'Sparse',start_pos);
%     disp(strcat('Comp. time (BayesSparse): ',num2str(toc)));
%     close all;
%     cd ..\..\..\Quality_Indices
%     QI_BayesSparse = QualityIndices(I_BayesSparse(5:end-4,5:end-4,:),I_REF(5:end-4,5:end-4,:),ratio);
%     cd ..\..\..\GF
    %% end of pansharpining
    %% plot qualities
    
    minCC = min( [ min( [ QI_GF_mGF_Res.ccMap QI_GF_mGF_Res21.ccMap ])] );
    maxCC = max( [ max( [ QI_GF_mGF_Res.ccMap QI_GF_mGF_Res21.ccMap ])] );
    figure(1111),hold on,
%     plot(dataset.wavelength,QI_BayesSparse.ccMap,'b','LineWidth',2), hold on,
    plot(dataset.wavelength,QI_GFPCA.ccMap,'b','LineWidth',2), hold on,
%     plot(dataset.wavelength,QI_GF_mGF.ccMap,'r','LineWidth',2),
    plot(dataset.wavelength,QI_GF_mGF_Res.ccMap,'r','LineWidth',2),
    plot(dataset.wavelength,QI_GF_mGF_Res21.ccMap,'g','LineWidth',2),    
    plot(dataset.wavelength,QI_GF_mGF_Res31.ccMap,'k','LineWidth',2),
    plot(dataset.wavelength,QI_GF_mGF_Res41.ccMap,'m','LineWidth',2),
    plot(dataset.wavelength,QI_GF_mGF_Res61.ccMap,'c','LineWidth',2),  
    plot(dataset.wavelength,QI_GF_mGF_Res71.ccMap,'y','LineWidth',2),  
    legend('GFPCA','oransal GFR','oransal GFR2','oransal GFR3','oransal GFR4','oransal GFR6','oransal GFR7'), title('CC')
    axis([min(dataset.wavelength) max(dataset.wavelength) 0.9 maxCC+0.01]);
    
    minRMSE = min( [ min( [ QI_GF_mGF_Res.rmseBands; QI_GF_mGF_Res21.rmseBands ] ) ] );
    maxRMSE = max( [ max( [ QI_GF_mGF_Res.rmseBands; QI_GF_mGF_Res21.rmseBands ] ) ] );
    
    figure(11112), hold on,
%     plot(dataset.wavelength,QI_BayesSparse.rmseBands,'b','LineWidth',2),hold on,
    plot(dataset.wavelength,QI_GFPCA.rmseBands,'b','LineWidth',2), hold on,
%     plot(dataset.wavelength,QI_GF_mGF.rmseBands,'r','LineWidth',2),
    plot(dataset.wavelength,QI_GF_mGF_Res.rmseBands,'r','LineWidth',2),
    plot(dataset.wavelength,QI_GF_mGF_Res21.rmseBands,'g','LineWidth',2),
    plot(dataset.wavelength,QI_GF_mGF_Res31.rmseBands,'k','LineWidth',2),
    plot(dataset.wavelength,QI_GF_mGF_Res41.rmseBands,'m','LineWidth',2),
    plot(dataset.wavelength,QI_GF_mGF_Res61.rmseBands,'c','LineWidth',2),    
    plot(dataset.wavelength,QI_GF_mGF_Res71.rmseBands,'y','LineWidth',2),
    legend('GFPCA','oransal GFR','oransal GFR2','oransal GFR3','oransal GFR4','oransal GFR6','oransal GFR7'), title('RMSE')
    axis([min(dataset.wavelength) max(dataset.wavelength) minRMSE maxRMSE+0.1]);
     
%     figure(2222),
%     minSAMVal = min( [min(QI_GFPCA.SAMmap(:)),min(QI_GF_mGF.SAMmap(:)),min(QI_GF_mGF_Res.SAMmap(:))] );
%     maxSAMVal = max( [max(QI_GFPCA.SAMmap(:)),max(QI_GF_mGF.SAMmap(:)),max(QI_GF_mGF_Res.SAMmap(:))] );
%     subplot(221), imshow(QI_GFPCA.SAMmap,[minSAMVal,maxSAMVal]); title('GFPCA'),
%     subplot(222), imshow(QI_GF_mGF.SAMmap,[minSAMVal,maxSAMVal]);title('oransal GF');
%     subplot(223), imshow(QI_GF_mGF_Res.SAMmap,[minSAMVal,maxSAMVal]);title('oransal GFR');
%     %     subplot(223), imshow(QI_BayesSparse.SAMmap,[minSAMVal,maxSAMVal]);title('Bayes Sparse');
%     truesize([300 300]);
    
%         h3333 = figure(3333);
%         for i=1:1:size(I_HS,3)
%     
%             figure(h3333),
%             clf;
%             band1 = I_REF(:,:,i);
%             band2 = I_GFPCA(:,:,i);
%             band3 = I_HS_BRatio(:,:,i);
%             band4 = I_HS_GRatioSpatial(:,:,i);
%             band5 = I_HS_BGRIRRatio(:,:,i);
%             band6 = I_HS_BGRIRSpatial(:,:,i);
%     
%     
%             subplot(2,3,1), imshow(band1,[]), title(strcat('Orj ',num2str(i)));
%             subplot(2,3,2), imshow(band2,[]), title(strcat('GFPCA ',num2str(i)));
%             subplot(2,3,4), imshow(band3,[]), title(strcat('BRatio ',num2str(i)));
%             subplot(2,3,4), imshow(band4,[]), title(strcat('GRatioSpatial ',num2str(i)));
%             subplot(2,3,3), imshow(band5,[]), title(strcat('RGRIR ',num2str(i)));
%             subplot(2,3,5), imshow(band6,[]), title(strcat('RGRIRSpatial ',num2str(i)));
%             truesize([200 200]);
%             drawnow;
%             waitforbuttonpress;
%         end
    %
    %     h4444 = figure(1); imshow(I_REF(:,:,13),[]); truesize([400,400]);
    %     h5555 = figure(2);
    %     figure(h4444);
    %     for i=1:1:6
    %         pause(0.1);
    %         figure(h4444),
    %         [ypos,xpos] = ginput(1);
    %         ypos= round(ypos); xpos= round(xpos);
    %
    %         pause(0.1);
    %         figure(h5555), clf, hold on,
    %         plot(squeeze(dataset.I_REF(xpos,ypos,:)),'b','LineWidth',2);
    %         plot(squeeze(I_GFPCA(xpos,ypos,:)),'r','LineWidth',2);
    %         plot(squeeze(I_HS_BRatio(xpos,ypos,:)),'g','LineWidth',2);
    %         plot(squeeze(I_HS_GRatioSpatial(xpos,ypos,:)),'m','LineWidth',2);
    %         plot(squeeze(I_HS_BGRIRRatio(xpos,ypos,:)),'c','LineWidth',2);
    %         plot(squeeze(I_HS_BGRIRSpatial(xpos,ypos,:)),'k','LineWidth',2);
    %         legend('Orj','GFPCA','BRatio','GRatioSpatial','BGRIR','BGRIRSpatial')
    %     end
    
    %% end of plotting qualities
    %%
    
    
    %     % Bayesian Sparse
    %     cd 'BayesFusion'
    %     setup;
    %     tic;
    %     [I_BayesSparse]= BayesianFusion(I_HS,I_PAN,overlap,KerBlu,ratio,'Sparse',start_pos);
    %     t = toc;
    %     disp(strcat('Comp. time (BayesSparse): ',num2str(t)));
    %     cd ..
    %     QI_BayesSparse = QualityIndices(I_BayesSparse(5:end-4,5:end-4,:),I_REF(5:end-4,5:end-4,:),ratio);
    %     cd ..\..
    %     multibandwrite(I_BayesSparse,'Outputs/I_BayesSparse',INTERLEAVE);
    %     cd ilerleme9\Methods
    
    %     save(strcat(strDataName,'_',strDataType,'_GFPCA'));
end
% 
% for i=1:1:size(I_GFPCA,3)
% %     
% %     rmse1 = (I_REF(:,:,i) - I_REF(:,:,i)).^2;
% %     rmse2 = (I_HS_mGF_Res41(:,:,i) - I_REF(:,:,i)).^2;
% %     
% %     minV = min([rmse1(:); rmse2(:)]);
% %     maxV = max([rmse1(:); rmse2(:)]);
% %     
%     figure(11),
%     subplot(1,3,1), imshow(I_REF(:,:,i),[]); title(strcat('REF band : ',num2str(i) ))
%     subplot(1,3,2), imshow(I_HS(:,:,i),[]); truesize([300 300]); title(strcat('ORJ' )); 
%     subplot(1,3,3), imshow(I_HS_mGF_Res41(:,:,i),[]); truesize([300 300]); title(strcat('MY' )); 
%     drawnow(); 
%     waitforbuttonpress();
% end


