% demonstrate usage of grayness index
% Copyright (c) by Yanlin Qian, yanlin.qian@tut.fi, 2018-3-22.
% Credits to Kaifu Yang for his first version gray pixel code.


%%path and load data
clear all, close all
addpath(genpath('./greypixel_kaifu'));
load('exampleimg.mat');


%% grayness index
%produce some intermediate figures for paper.
param.visualization.histogram=0
param.visualization.bigsequenceimage=0;
param.visualization.sequence_dir='./sequence/';
param.real_rgb=gt;
param.visualization.greypixel_comparison=0;
param.visualization.comp_dir='./comp/';

%whether use illumination prior
param.prior.use=0;

%whether use binclip to gray pixel counting
param.binclip.use=0;
    
%two parameters, tuned already
list_Npre=10.^[-1];%[-2,-1,0];
list_threshold=10.^[-4];%[-5,-4,-3];
Npre = list_Npre(1); % n% = 0.01%
delta_threshold=list_threshold(1);

%choose Npixels.
Npixels = size(input_im,1)*size(input_im,2);
numGPs=floor(Npre*Npixels/100); 

%mask saturated pixels and mask very dark pixels
mask=(max(input_im,[],3)>=0.95) | (sum(input_im,3)<=0.0315);
 

param.delta_threshold=delta_threshold;
param.mask=mask;
param.numGPs=numGPs;
param.runtime.i=i;
param.runtime.name_img='example.jpg';
EvaLum=GPconstancy_GI(input_im,param);
sprintf('pred error: %0.2f degree',acos(normr(EvaLum)*gt')*180/pi)

