function image = rgb2srgb(image, is_int)


if ndims(image) == 3
    [n1,n2,n3] = size(image);
    image = reshape(image,n1*n2,n3);
    do_reshape = 1;
elseif ndims(image) == 2
    do_reshape = 0;
end

ind = image<=0.0031308;
image(ind) = image(ind) .* 12.92;
image(~ind) = 1.055 * (image(~ind).^(1.0/2.4)) - 0.055;

image = image./max(image(:));
if is_int
    image = uint8(round(255 * image));
end


if do_reshape
    image = reshape(image,n1,n2,n3);
end
