function[]=funcPlotGradValues(strDataName,gradMsiBands,gradHsiBands,figId)

if nargin < 4
    figId = 20;
end

switch(strDataName)
   
    case {'Salinas','salinas','SALINAS'}
        plotGradValuesSalinas(gradMsiBands,gradHsiBands,figId);
    case {'Pavia','pavia','PAVIA'}
        plotGradValuesPavia(gradMsiBands,gradHsiBands,figId);
    case {'Sentetic','sentetic','SENTETIC'}
        plotGradValuesSentetic(gradMsiBands,gradHsiBands,figId);
    otherwise
        warning('please define avalid data name');
        return;
end
end

%% SALINAS
function []=plotGradValuesSalinas(gradMsiBands,gradHsiBands,figId)
%%%% tek-çift kontrolu
nBand = size(gradHsiBands.B,2);
artm = 1;
figure(figId), 
for i=1:2:size(gradHsiBands.B,2)
    subplot(2,4,artm); hold on,
    plot3(gradMsiBands.B,gradMsiBands.G,gradMsiBands.R,'.','LineWidth',2,'Color','b'); grid on;
    plot3(gradHsiBands.B(:,i),gradHsiBands.G(:,i),gradHsiBands.R(:,i),'.','LineWidth',2,'Color','r');
    if (~isequal(i,nBand))
        plot3(gradHsiBands.B(:,i+1),gradHsiBands.G(:,i+1),gradHsiBands.R(:,i+1),'.','LineWidth',2,'Color','g');
        title(strcat('MSI - HSI( ',num2str(i),' , ',num2str(i+1),' )'));
    else
        title(strcat('MSI - HSI( ',num2str(i),' )'));
    end    
    xlabel('Blue'), ylabel('Green'), zlabel('Red');
    artm=artm+1;
end

end

%% SENTETIC
function []=plotGradValuesSentetic(gradMsiBands,gradHsiBands,figId)
artm = 1;
figure(figId), 
for i=1:2:size(gradHsiBands.B,2)
    subplot(2,4,artm); hold on,
    plot3(gradMsiBands.B,gradMsiBands.G,gradMsiBands.R,'.','LineWidth',2,'Color','b'); grid on;
    plot3(gradHsiBands.B(:,i),gradHsiBands.G(:,i),gradHsiBands.R(:,i),'.','LineWidth',2,'Color','r');
    if (~isequal(i,nBand))
        plot3(gradHsiBands.B(:,i+1),gradHsiBands.G(:,i+1),gradHsiBands.R(:,i+1),'.','LineWidth',2,'Color','g');
        title(strcat('MSI - HSI( ',num2str(i),' , ',num2str(i+1),' )'));
    else
        title(strcat('MSI - HSI( ',num2str(i),' )'));
    end 
    artm=artm+1;
end

end

%% PAVIA
function []=plotGradValuesPavia(gradMsiBands,gradHsiBands,figId)
artm = 1;
figure(figId),
for i=1:2:size(gradHsiBands.B,2)
    subplot(2,4,artm); hold on,
    plot3(gradMsiBands.B,gradMsiBands.G,gradMsiBands.R,'.','LineWidth',2,'Color','b'); grid on;
    plot3(gradHsiBands.B(:,i),gradHsiBands.G(:,i),gradHsiBands.R(:,i),'.','LineWidth',2,'Color','r');
    if (~isequal(i,nBand))
        plot3(gradHsiBands.B(:,i+1),gradHsiBands.G(:,i+1),gradHsiBands.R(:,i+1),'.','LineWidth',2,'Color','g');
        title(strcat('MSI - HSI( ',num2str(i),' , ',num2str(i+1),' )'));
    else
        title(strcat('MSI - HSI( ',num2str(i),' )'));
    end 
    artm=artm+1;
end

end