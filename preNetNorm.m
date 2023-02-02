function [imout] = preNetNorm(imin)

imin = im2double(imin);

u = [0.485, 0.456, 0.406];
st =[0.229, 0.224, 0.225];

for i = 1:3
U(:,:,i) = ones(size(imin,1),size(imin,2)).*u(i);
ST(:,:,i) = ones(size(imin,1),size(imin,2)).*st(i);
end


imnorm = (imin - U) ./ ST;

imsize = imresize(imnorm,[250,250]);

targetSize = [224 224];

win1 = centerCropWindow2d(size(imsize),targetSize);

imout = imcrop(imsize,win1);

