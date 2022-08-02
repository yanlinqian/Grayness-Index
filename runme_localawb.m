% demonstrate usage of grayness index
% Copyright (c) by Yanlin Qian, yanlin.qian@tut.fi, 2018-3-22.
% Credits to Kaifu Yang for his first version gray pixel code.


%%path and load data
clear all, close all
addpath(genpath('./greypixel_kaifu'));
real_world_flag=1;   % real_world_flag=1 for real-world dataset / real_world_flag=0 laboratory dataset

%%reset your mimo dataset
mimo='D:\RESEARCH\mimo\';

if(real_world_flag)
    load(fullfile(mimo,'image_names_real.mat'));
    pathImages=fullfile(mimo,'/realworld/img/');
    pathGT=fullfile(mimo,'/realworld/groundtruth/');
    pathMasks=fullfile(mimo,'/realworld/masks/'); 
else
    load(fullfile(mimo,'image_names_lab.mat'));
    pathImages=fullfile(mimo,'/lab/img/');
    pathGT=fullfile(mimo,'/lab/groundtruth/');
    pathMasks=fullfile(mimo,'/lab/masks/');
end

nIm=size(image_names,2);  % number of images
errors=zeros(nIm,1);



%%
list_mean=zeros(1,1);
list_median=zeros(1,1);
for j=1:1
for ii=1:nIm  % loop over images;
    input_im = double(imread([pathImages,image_names{ii},'.png']));
    GT_im = double(imread([pathGT,image_names{ii},'.png']));
    mask = double(imread([pathMasks,image_names{ii},'.png']));
    %Npre=10;Inum=6;sig=0.2;
    Npre=10;Inum=6;sig=0.2;
    Npixels = size(input_im,1)*size(input_im,2);
    numGPs=floor(Npre*Npixels/100); 
    input_im=input_im/max(input_im(:));
    [CorrImg,EstIl] = MultiLumConstancy_dgp(input_im,numGPs,Inum,mask,'delta_threshold',10^(-4),'sig',sig);

    adist=angDistPixelwise(GT_im.*repmat(mask,[1,1,3]),EstIl);
    errors(ii)=mean(adist)/pi*180;          % error in degrees
    fprintf('the %dth image, angular error %.2f\n',ii,errors(ii))
    
end
fprintf(1,'The mean error=%f\n',mean(errors(errors==errors)));
fprintf(1,'The median error=%f\n',median(errors(errors==errors)));

list_mean(j)=mean(errors(errors==errors));
list_median(j)=median(errors(errors==errors));
end
mean(list_mean)
mean(list_median)

