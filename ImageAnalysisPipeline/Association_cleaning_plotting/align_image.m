%% Align source image to target by trying out all shifts within a given boundary.

function [shift1, shift2] = align_image(source, target, max_shift)

image_size = size(source);
sum_product_matrix = NaN*ones(2*max_shift + 1, 2*max_shift + 1);
for x1 = (-max_shift):max_shift
    for x2 = (-max_shift):max_shift
         sum_product_matrix(x1 + max_shift + 1, x2 + max_shift + 1) =...
             sum(sum(source((1+max_shift - x1):(image_size(1) - max_shift - x1),...
             (1+max_shift - x2):(image_size(2) - max_shift - x2)).*...
             target((1+max_shift):(image_size(1) - max_shift), (1+max_shift):(image_size(2) - max_shift))));
         
    end
end

[shift1, shift2] =  find(sum_product_matrix == max(max(sum_product_matrix)));
shift1 = shift1 - max_shift;
shift2 = shift2 - max_shift;
%% show score space.
figure,
[xx, yy] = meshgrid((-max_shift):max_shift, (-max_shift):max_shift);
surface(xx, yy,sum_product_matrix)
title('Shift score map.');


return