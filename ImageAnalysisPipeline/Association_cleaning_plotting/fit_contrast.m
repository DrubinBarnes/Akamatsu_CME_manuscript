function image = fit_contrast(image)
image = (image - min(min(image)))/(max(max(image)) - min(min(image)));
return