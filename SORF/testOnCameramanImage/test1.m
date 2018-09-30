clear all; close all; clc;

img = double(imread('cameraman.tif'));
img = img(1:end-1,1:end-1);
[w,h] = size(img);

img2 = zeros(size(img));
rows = 2:3:w; cols = 2:3:h;

img2(rows,cols) = img(rows,cols);

agirlikMatrisi = edgePreserveInterpolationAlg1(img,img2);
agirlikMatrisi2= edgePreserveInterpolationAlg2(img,img2);
agirlikMatrisi3= edgePreserveInterpolationAlg3(img,img2);
agirlikMatrisi4= edgePreserveInterpolationAlg4(img,img2);

img3 = img2(rows,cols);
imgBicubic = imresize(img3,[255,255],'bicubic');

figure(1), 
subplot(231), imshow(img,[]); title('Orj');
subplot(232), imshow(imgBicubic,[]); title('Bicubic');
subplot(233), imshow(agirlikMatrisi2,[]); title('IntAlg-2');
subplot(234), imshow(agirlikMatrisi3,[]); title('IntAlg-3');
subplot(235), imshow(agirlikMatrisi4,[]); title('IntAlg-4');

