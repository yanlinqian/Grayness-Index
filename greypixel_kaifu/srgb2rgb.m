function I = srgb2rgb(image)

if ~isfloat(image), 
    I = double(image);
else
    I = image;
end

if max(I(:)) > 1
    I  = double(I)./255;
end

ind = I <= 0.04045;
I(ind) = I(ind)./12.92;

I(~ind) = (((I(~ind)+0.055)/1.055).^2.4);
