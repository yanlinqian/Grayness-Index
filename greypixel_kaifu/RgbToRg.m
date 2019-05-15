function rg = RgbToRg(rgb)
%Turns an RGB matrix (3xn) specifying the color of illuminants
assert(size(rgb,1) == 3);
norm2=sum(rgb.^2,1);
rg = [(rgb(1,:).^2 ./ norm2(1,:)); (rgb(2,:).^2 ./ norm2(1,:))];
end

