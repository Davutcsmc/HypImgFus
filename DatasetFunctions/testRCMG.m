clear all; close all; clc;

dataType = 'dcmall'; % dcmall %pine 

switch(dataType)    
    case {'Pine','pine','PINE'}
        load('C:\Users\davut.cesmeci\Google Drive\2016_fusion\FusionWorks\codes\data\Pine\pineData');
    case {'DCMall','dcmall','DCMALL'}
        load('C:\Users\davut.cesmeci\Google Drive\2016_fusion\FusionWorks\codes\data\DCmall\dcmall_1');
    otherwise
        warning('please define avalid data name');
        return;
end

[gradMap] = funcVecGradient(hsi,'RCMG');
imshow(gradMap,[]);

