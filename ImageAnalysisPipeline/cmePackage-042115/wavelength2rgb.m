function rgbValue = wavelength2rgb(lambda)
% Converts wavelength to RGB approximation of the visible spectrum.
% INPUT :   lambda : wavelength in [m]
%
% Code based on http://www.midnightkite.com/color.html

% Francois Aguet, October 2010

lambda = lambda*1e9;
gammaVal = 0.8;

if lambda >= 380 && lambda < 440
    rgbValue = [-(lambda-440)/60 0 1];
end
if lambda >= 440 && lambda < 490
    rgbValue = [0 (lambda-440)/50 1];
end
if lambda >= 490 && lambda < 510
    rgbValue = [0 1 -(lambda-510)/20];
end
if lambda >= 510 && lambda < 580
    rgbValue = [(lambda-510)/70 1 0];
end
if lambda >= 580 && lambda < 645
    rgbValue = [1 -(lambda-645)/65 0];
end
if lambda >= 645 && lambda <= 780
    rgbValue = [1 0 0];
end

% Attenuate intensity near limits
if lambda > 700
    window = 0.3 + 0.7*(780-lambda)/80;
elseif lambda < 420
    window = 0.3 + 0.7*(lambda-380)/40;
else
    window = 1;
end

% Gamma
rgbValue = rgbValue * window * gammaVal; 