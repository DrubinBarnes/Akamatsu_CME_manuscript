function image = fit_contrast_median(image)
image_median = median(reshape(image, [1, numel(image)]));
image = (image - image_median)/(max(max(image)) - image_median);
image(image<0) = 0;
return