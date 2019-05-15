% angerr: Computes the angular error between two illuminants
%
% err = angerr(l1,l2)
%
%   l1,l2: Two illuminant color vectors
% RESULT
%   err:   Angular error in degrees
function r = angerr(l1,l2)
  l1 = l1 / sqrt(sum(l1.^2));
  l2 = l2 / sqrt(sum(l2.^2));
  r = acosd(sum(l1(:).*l2(:)));
  
