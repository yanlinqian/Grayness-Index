function Lstd = LocalStd(img,wsize)
% function mag= DerivGauss(img,sigma)
% Design the filters - a derivative of 2D gaussian
% inputs:
%         img ---- one channel of the color-biased image.
%         wsize -- Window size
% outputs:
%         Lstd --- Local contrast (STD)
%
% Kaifu Yang <yang_kf@163.com>
% March 2015
%=========================================================================%

Mi = wsize(1); Mj = wsize(2);
Emap = padarray(img,[floor(Mi/2) floor(Mj/2)],'replicate','both');
[Ex Ey] = size(Emap);
contrast = colfilt(double(Emap), [Mi Mj], 'sliding', @std ); % get contrast = local std
Lstd = contrast(floor(Mi/2)+1:Ex-floor(Mi/2),floor(Mj/2)+1:Ey-floor(Mj/2));

%=========================================================================%