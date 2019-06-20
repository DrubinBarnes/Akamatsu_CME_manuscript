function image = fit_contrast_ignore_zero(image)
image_above0 = image(image>0);
image_min = min(min(image_above0));
image_max = max(max(image_above0));
image = (image - image_min)/(image_max - image_min);
image(image<0) = 0;
return