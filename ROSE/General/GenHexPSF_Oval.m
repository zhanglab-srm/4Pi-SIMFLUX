function result = GenHexPSF_Oval(imgsize, x0, y0, rx, ry, std_psf, intlist, phaselist, PhaseWidth)
% pixelsize = 107;
% std_psf = 1.5;
% r = 5;


% phase = [30, 90, 150, 210, 270, 330];
% PhaseWidth = -15:1:15;

% intlist = rand(1,6);
% imgsize = size(xx,1);
% phase = phase / 180 * pi;
% PhaseWidth = PhaseWidth / 180 * pi;

if nargin < 9   %PahseWidth
    PhaseWidth = 0;
end
PhaseWidth = PhaseWidth / 180 * pi;

if nargin < 8   %phaselist
    phaselist = [30, 90, 150, 210, 270, 330];
    phaselist = phaselist / 180 * pi;
end

imgsize_half = (imgsize+1)/2;
x0 = imgsize_half + x0;
y0 = imgsize_half + y0;
img = zeros(imgsize, imgsize);

for m=1:length(phaselist)
    tint = intlist(m) / length(PhaseWidth);
    for n=1:length(PhaseWidth)
        tphase = phaselist(m) + PhaseWidth(n);
        img = AddGaussian2D(img, x0 + rx * cos(tphase), y0 + ry * sin(tphase),...
            tint, std_psf);
    end
    
end

result = img;
