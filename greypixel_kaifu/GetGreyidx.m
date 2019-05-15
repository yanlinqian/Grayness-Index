function Greyidx = GetGreyidx(img,method,scale)
% function Greyidx = GetGreyidx(img,method,scale)
% Obtain the grey index of each pixel
% inputs:
%         img ------ Input color-biased image.
%         method --- The method to compute IIM, method = {'GPedge','GPstd'}
%         scale ---- Filter scale
% outputs:
%        Greyidx --- Grey index map
%
% Kaifu Yang <yang_kf@163.com>
% March 2015
%=========================================================================%

[rr cc dd] = size(img);
if dd~=3
    error('The input image must be RGB-format color image!');
end

R=img(:,:,1);
G=img(:,:,2);
B=img(:,:,3);

R(R==0)=eps;
G(G==0)=eps;
B(B==0)=eps;

if strcmp(method,'GPedge')
    sigma = scale;
    Mr = DerivGauss(log(R),sigma);
    Mg = DerivGauss(log(G),sigma);
    Mb = DerivGauss(log(B),sigma);
 
elseif strcmp(method,'GPstd')
    wsize = [scale scale];
    Mr = LocalStd(log(R),wsize);
    Mg = LocalStd(log(G),wsize);
    Mb = LocalStd(log(B),wsize);
end

data=[Mr(:),Mg(:),Mb(:)];
Ds= std(data,[],2);
Ds = Ds./(mean(data,2)+eps);

data1 = [R(:),G(:),B(:)];
Ps = Ds./(mean(data1,2)+eps);
Greyidx = reshape(Ps, [rr cc]);

Greyidx = Greyidx./(max(Greyidx(:))+eps);
Greyidx(Mr<eps & Mg<eps & Mb<eps)= max(Greyidx(:));

hh = fspecial('average',[7 7]);
Greyidx = imfilter(Greyidx,hh,'circular');

%=========================================================================%

