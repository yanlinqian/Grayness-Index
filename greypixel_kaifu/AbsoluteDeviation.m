function mag = AbsoluteDeviation(img,width)
% function mag= DerivGauss(img,sigma)

%=========================================================================%
dgau2D=-1*ones([width,width]);
dgau2D(round(width/2),round(width/2))=width*width-1;

mag = imfilter(img, dgau2D, 'conv','replicate');
mag=abs(mag);
%=========================================================================%
