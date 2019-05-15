function [JDark,JDark_Ind] = darkChannel_our(im2)

[height, width, ~] = size(im2);

patchSize = 15; %the patch size is set to be 15 x 15
padSize = 7; % half the patch size to pad the image with for the array to 
%work (be centered at 1,1 as well as at height,1 and 1,width and height,width etc)

JDark = zeros(height, width); % the dark channel
JDark_Ind=zeros(height,width);
imJ = padarray(im2, [padSize padSize], Inf); % the new image
% imagesc(imJ); colormap gray; axis off image

for j = 1:height
    for i = 1:width
        % the patch has top left corner at (jj, ii)
        patch = imJ(j:(j+patchSize-1), i:(i+patchSize-1),:);
        % the dark channel for a patch is the minimum value for all
        % channels for that patch
        [JDark(j,i) JDark_Ind(j,i)] = min(patch(:));
     end
end
