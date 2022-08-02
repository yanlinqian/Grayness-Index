function [CorrImg MultiLum] = MultiLumConstancy_dgp(img,numGPs,Inum,mask,varargin)
% function [CorrImg MultiLum] = MultiLumConstancy(img,numGPs,Inum)
% inputs:
%         img   ------ Input color-biased image.
%         numGPs ----- The number of grey pixels.
%         Inum  ------ The number of illuminant.
% outputs:
%         CorrImg -----Corrected image.
%         MultiLum  --- pixel-wise illuminant.
% Main function for performing color constancy system using Grey Pixels in paper:
% Kaifu Yang, Shaobing Gao, and Yongjie Li*.
% Efficient Illuminant Estimation for Color Constancy Using Grey Pixels.CVPR, 2015.
%
% Contact:
% Visual Cognition and Computation Lab(VCCL),
% Key Laboratory for NeuroInformation of Ministry of Education,
% School of Life Science and Technology(SLST),
% University of Electrical Science and Technology of China(UESTC).
% Address: No.4, Section 2, North Jianshe Road,Chengdu,Sichuan,P.R.China, 610054
% Website: http://www.neuro.uestc.edu.cn/vccl/home.html

% Kaifu Yang <yang_kf@163.com>;
% March 2015
%========================================================================%



[ww hh dd] = size(img);
EvaLum =zeros(size(img));

R=img(:,:,1);
G=img(:,:,2);
B=img(:,:,3);

R(R==0)=eps;
G(G==0)=eps;
B(B==0)=eps;

% % Algorithm 1 -- using edge as IIM
%GreyEdge = GetGreyidx(img,'GPedge',.5);
%Greyidx = GreyEdge;

% Algorithm 2 -- using Local Contrast as IIM
%GreyStd = GetGreyidx(img,'GPstd',3);
%Greyidx = GreyStd;

    % Algorithm 3 -- using grayness function of DGP
    for i=1:2:nargin-4
    eval(sprintf('%s = varargin{%d+1};',varargin{i},i));
    end
    
    input_im=img;
    img_column=reshape(input_im,[],3);
    r=input_im(:,:,1); g=input_im(:,:,2); b=input_im(:,:,3);
    
    %add this if we are running on natural image, as they are with a lot noise.
    %averagiong
    %average_filter = fspecial('average',[7 7]);
    %r = imfilter(r,average_filter,'circular');g = imfilter(g,average_filter,'circular');b = imfilter(b,average_filter,'circular');
    
%     %mask 0 elements
%     mask=(r==0) | (g==0) | (b==0);
    r(r==0)=eps;  g(g==0)=eps;  b(b==0)=eps;norm1=r+g+b;

    %mask low contrast pixels     
    delta_r=DerivGauss(r,.5); delta_g=DerivGauss(g,.5); delta_b=DerivGauss(b,.5);
    mask_deltargb=(delta_r<=delta_threshold & delta_g<=delta_threshold & delta_b<=delta_threshold);

   
    log_r=log(r)-log(norm1); log_b=log(b)-log(norm1);
    
    delta_log_r=DerivGauss(log_r,.5); 
    delta_log_b=DerivGauss(log_b,.5);
    mask_deltargb=mask_deltargb | (delta_log_r==Inf) | (delta_log_b==Inf);
    

    data=[(delta_log_r(:)),(delta_log_b(:))];
    mink_norm=2;
    norm2_data=power(sum(power(data,mink_norm),2),1/mink_norm);
    norm2_data(mask_deltargb==1)=max(norm2_data(:));
    norm2_data(mask==0)=max(norm2_data(:));

    map_uniquelight=reshape(norm2_data,size(delta_log_r));

    %average
    average_filter = fspecial('average',[3 3]);
    map_uniquelight = imfilter(map_uniquelight,average_filter,'circular');
    
    Greyidx=map_uniquelight;
    
    
    %save GI
    if exist('name_img')
        
        Map       = jet(255);
        sequence_6=(map_uniquelight)/(max(map_uniquelight(:)));
        sequence_6=-log(sequence_6);
        sequence_6=1-sequence_6/max(sequence_6(:));
        sequence_6=uint8(sequence_6*255.0);
        sequence_6= ind2rgb(sequence_6, Map);
        imwrite(sequence_6,name_img)
    end

tt=sort(Greyidx(:));

Gidx = zeros(size(Greyidx));
Gidx(Greyidx<=tt(numGPs)) = 1;

[row col]=find(Gidx==1);
[gIdx,TempCent]=k_means([row col],Inum);
while size(TempCent,1)<Inum
    [gIdx,TempCent]=k_means([row col],Inum);
end


II=zeros(Inum,3);
Cent = zeros(Inum,2);
[xx,yy]=meshgrid(1:hh,1:ww);
Dist = zeros(ww*hh,Inum);
for ki=1:Inum
    Gidxk=zeros(size(Greyidx));
    Cr = row(gIdx==ki);
    Cc = col(gIdx==ki);
    for i = 1:length(Cr)
        Gidxk(Cr(i),Cc(i))=1;
    end
    RR = Gidxk.*R;
    GG = Gidxk.*G;
    BB = Gidxk.*B;
    II(ki,:) = [sum(RR(:)) sum(GG(:)) sum(BB(:))];
    Cent(ki,:) = floor(TempCent(ki,:));
    Dist(:,ki) = sqrt((Cent(ki,1)-yy(:)).^2 + (Cent(ki,2)-xx(:)).^2)...
        ./sqrt(ww^2 + hh^2);
end

tr = zeros(ww*hh,1);
tg = zeros(ww*hh,1);
tb = zeros(ww*hh,1);
Wi = exp(-(Dist./(2*sig^2)));
Swi=sum(Wi,2);
for jj=1:Inum
    Wi(:,jj)= Wi(:,jj)./Swi;
    tr = tr + Wi(:,jj)*II(jj,1);
    tg = tg + Wi(:,jj)*II(jj,2);
    tb = tb + Wi(:,jj)*II(jj,3);
end
tempr = reshape(tr,[ww hh]);
tempg = reshape(tg,[ww hh]);
tempb = reshape(tb,[ww hh]);

EvaLum(:,:,1) = tempr;
EvaLum(:,:,2) = tempg;
EvaLum(:,:,3) = tempb;


coff = sum(EvaLum,3);
EvaLum(:,:,1) = EvaLum(:,:,1)./coff; % normalization
EvaLum(:,:,2) = EvaLum(:,:,2)./coff; % normalization
EvaLum(:,:,3) = EvaLum(:,:,3)./coff; % normalization

EvaLum(:,:,1) = EvaLum(:,:,1).* mask; % normalization
EvaLum(:,:,2) = EvaLum(:,:,2).* mask; % normalization
EvaLum(:,:,3) = EvaLum(:,:,3).* mask; % normalization

CorrImg(:,:,1) = img(:,:,1)./EvaLum(:,:,1) .* mask;
CorrImg(:,:,2) = img(:,:,2)./EvaLum(:,:,2) .* mask;
CorrImg(:,:,3) = img(:,:,3)./EvaLum(:,:,3) .* mask;

EvaLum=reshape(normr(reshape(EvaLum,[],3)), size(EvaLum));
EvaLum(:,:,1) = EvaLum(:,:,1).* mask; % normalization
EvaLum(:,:,2) = EvaLum(:,:,2).* mask; % normalization
EvaLum(:,:,3) = EvaLum(:,:,3).* mask; % normalization
MultiLum = EvaLum;

%=========================================================================%

