function ei = illuminantEstimator(img, mask, prc)
%illuminantEstimator estimate the illuminant for a image
%
%   ei = illuminantEstimator(img, mask, prc)
%	
%   img : 3-channel normalized [0,1] image to estimate the illuminant 
%   mask: binary image mask indicates pixels not used in the estimation, 
%         e.g., color checker in the image and pixels possibly saturated
%   prc : percentage of pixels selected from both dark and bright side 
%
%
% Copyright (c) 2013 Dongliang Cheng
% School of Computing
% National University of Singapore
% http://www.nus.edu.sg/
%
% Permission is hereby granted, free of charge, to any person obtaining 
% a copy of this software and associated documentation files (the 
% "Software"), to deal in the Software with restriction for its use for 
% research purpose only, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included
% in all copies or substantial portions of the Software.
%
% The Software is provided "as is", without warranty of any kind.
%
% Please cite the following work if this program is used:
% [1] Cheng D.L., Prasad D., Brown M. S. "Illuminant Estimation for Color
% Constancy: why spatial domain methods work and the role of the color 
% distribution",  Journal of Optical Society of America - A (JOSA) 2014
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data = reshape(img, [], 3);
if mask
    data(mask(:),:)=[];
end

n = size(data,1);
% gray-world estimation
l = mean(data); l = l./norm(l); 
% projection on the gray-world estimation
data_p = sum(bsxfun(@times, data, l), 2); 

[data_p, idx] = sort(data_p);
data_selected = [data(idx(1:ceil(n*prc/100)),:); ...
		data(idx(floor(n*(100-prc)/100):end),:);];

sigma = data_selected'*data_selected;
[v, d] = eig(sigma);
[~, idx] = sort(diag(d));
v = v(:,idx);

ei = abs(v(:,3));