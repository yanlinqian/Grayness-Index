function [Greyidx,Greyidx_angular] = GetGreyidx_angular(img,method,scale)
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
%
% revisited by Yanlin Qian, yanlin.qian@tut.fi,2017
%
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
    
%     Fr = DerivGauss((R),sigma);
%     Fg = DerivGauss((G),sigma);
%     Fb = DerivGauss((B),sigma);
%     
%     delta_Fr = DerivGauss((Fr),sigma);
%     delta_Fg = DerivGauss((Fg),sigma);
%     delta_Fb = DerivGauss((Fb),sigma);
    
%     
%     F=Fr.*Fg.*Fb;
 
elseif strcmp(method,'GPstd')
    wsize = [scale scale];
    Mr = LocalStd(log(R),wsize);
    Mg = LocalStd(log(G),wsize);
    Mb = LocalStd(log(B),wsize);
    
elseif strcmp(method,'GPabso')
    Mr = AbsoluteDeviation(log(R),scale);
    Mg = AbsoluteDeviation(log(G),scale);
    Mb = AbsoluteDeviation(log(B),scale);
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


lumi = [R(:),G(:),B(:)];

mink_norm=1;
lumi_factor=power(sum(power(lumi,mink_norm),2),1/mink_norm)+eps;
[sort_lumi,sortid_lumi]=sort(lumi_factor,'descend');
threshold_lumi=sort_lumi(round(0.05*numel(lumi_factor)));


Greyidx_angular = reshape(angular_error,[rr,cc]);
angular_error=angular_error./lumi_factor;

Greyidx = Greyidx_angular./(max(Greyidx_angular(:))+eps);
Greyidx(Mr<eps & Mg<eps & Mb<eps)= max(Greyidx(:));

Greyidx_angular(Mr<eps & Mg<eps & Mb<eps)= max(Greyidx_angular(:));

hh = fspecial('average',[7 7]);
Greyidx = imfilter(Greyidx,hh,'circular');
hh = fspecial('average',[7 7]);
Greyidx_angular = imfilter(Greyidx_angular,hh,'circular');

%=========================================================================%

