function [Greyidx,Greyidx_angular] = GetGreyidx_angular_darkbright(img,method,scale,prc)
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
% Ds= std(data,[],2);
% Ds = Ds./(mean(data,2)+eps);
% 
% lumi = [R(:),G(:),B(:)];
% Ps = Ds./(mean(lumi,2)+eps);
% Greyidx = reshape(Ps, [rr cc]);

%added by yanlin
data(data(:,1)==0,1)=eps;
data(data(:,2)==0,2)=eps;
data(data(:,3)==0,3)=eps;

data_normed=normr(data);
gt=normr(ones(size(data)));
angular_error=real(acos(dot(data_normed,gt,2)));

lumi = [R(:),G(:),B(:)]+1.0;


data = reshape(img, [], 3);
n = size(data,1);
% gray-world estimation
l = mean(data); l = l./norm(l); 
% projection on the gray-world estimation
data_p = sum(bsxfun(@times, data, l), 2); 

[data_p, idx] = sort(data_p);
dark_idx=idx(1:ceil(n*prc/100));
bright_idx=idx(floor(n*(100-prc)/100):end);
% data_selected = [data(idx(1:ceil(n*prc/100)),:); ...
% 		data(idx(floor(n*(100-prc)/100):end),:);];
weight_idx=zeros(size(angular_error));
weight_idx(dark_idx)=1;
weight_idx(bright_idx)=1;

% angular_error=angular_error./(mean(lumi,2)+eps);
%angular_error=angular_error./(abs(mean(lumi,2)-256/2)+eps);
angular_error=angular_error.*weight_idx;
Greyidx_angular = reshape(angular_error,[rr,cc]);

% mean_error=mean(real(angular_error))*180/pi
% median_error=median(angular_error)*180/pi
% min_error=min(angular_error)*180/pi
% max_error=max(angular_error)*180/pi

Greyidx = Greyidx_angular./(max(Greyidx_angular(:))+eps);
Greyidx(Mr<eps & Mg<eps & Mb<eps)= max(Greyidx(:));

Greyidx_angular(Mr<eps & Mg<eps & Mb<eps)= max(Greyidx_angular(:));

hh = fspecial('average',[7 7]);
Greyidx = imfilter(Greyidx,hh,'circular');

hh = fspecial('average',[7 7]);
Greyidx_angular = imfilter(Greyidx_angular,hh,'circular');

%=========================================================================%

