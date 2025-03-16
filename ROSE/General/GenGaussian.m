function img = GenGaussian(imgsize, parameters)
%Gaussian PSF Generation
%parameter
%[x0 y0 std int ]
% 1  2   3   4  
if length(imgsize) ==1
    imgsize = [imgsize imgsize];
end

x0 = parameters(1) + (imgsize(2)+1)/2;
y0 = parameters(2) + (imgsize(1)+1)/2;
psf_std = parameters(3);
psf_int = parameters(4);
bkg = 0;

img = ones(imgsize).*bkg;
img = AddGaussian2D(img, x0, y0, psf_int, psf_std);