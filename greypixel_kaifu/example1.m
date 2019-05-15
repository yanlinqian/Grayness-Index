% demo: example for single-illuminant estimation
% Efficient Illuminant Estimation Using Grey Pixels.
% paper in CVPR 2015:
% Kaifu Yang, Shaobing Gao, and Yongjie Li*.
% Efficient Illuminant Estimation for Color Constancy Using Grey Pixels. CVPR, 2015.
% Website: http://www.neuro.uestc.edu.cn/vccl/home.html

% Kaifu Yang <yang_kf@163.com>
% March 2015
%=========================================================================%
clear; clc;

% parameter setting
Npre = 0.01; % n% = 0.01%
mask = []; % no mask needed for SFU321 lab dataset

load('sfulab_name'); % load the groundtruth_illuminants of SFU321 dataset
img = double(imread('demo_image284.tif')); % demo_image284
figure; imshow(img./max(img(:)),[]);

Npixels = size(img,1)*size(img,2);
numGPs=floor(Npre*Npixels/100); 
[outimg estimated_illuminants] = GPconstancy(img,numGPs,mask);
figure; imshow(outimg./max(outimg(:)),[]);

Err = angerr(estimated_illuminants,groundtruth_illuminants(284,:));
fprintf('Angular Error = %f\n',Err);

fprintf(2,'======== THE END ========\n');
%=========================================================================%
