% demo: example for multi-illuminant estimation
% Efficient Illuminant Estimation Using Grey Pixels.
% paper in CVPR 2015:
% Kaifu Yang, Shaobing Gao, and Yongjie Li*.
% Efficient Illuminant Estimation for Color Constancy Using Grey Pixels. CVPR, 2015.
% Website: http://www.neuro.uestc.edu.cn/vccl/home.html

% Kaifu Yang <yang_kf@163.com>
% March 2015
%=========================================================================%

clc;  clear;

% parameter setting
Inum = 2;   % number of luminants
Npre= 10;  % n% = 10%

img = double(imread('multi-Lums-original.png'));
GT_im = double(imread('multi-Lums-GT.png'));
mask = double(imread('multi-Lums-mask.png'));

figure;imshow(img./max(img(:)),[]);
figure;imshow(GT_im./max(GT_im(:)),[]);

Npixels = size(img,1)*size(img,2);
numGPs=floor(Npre*Npixels/100);

% Multi-Illuminant estimation
[CorrImg MultiLum] = MultiLumConstancy(img,numGPs,Inum);

figure;imshow(CorrImg./max(CorrImg(:)),[]);
figure;imshow(MultiLum./max(MultiLum(:)),[]);

adist=angDistPixelwise(GT_im.*repmat(mask,[1,1,3]),MultiLum);
Err = mean(adist(:))/pi*180;

fprintf('Average Angular Error = %f\n',Err);

fprintf(2,'======== THE END ========\n');
%=========================================================================%

