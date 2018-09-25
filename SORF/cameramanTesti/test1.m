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

img3 = img2(rows,cols);
imgBicubic = imresize(img3,[255,255],'bicubic');


