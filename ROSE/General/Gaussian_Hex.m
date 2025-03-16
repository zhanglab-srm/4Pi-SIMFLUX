%Hex Gaussian PSF Generation
%parameter
%[x0 y0 std int1 int2 int3 int4 int5 int6 bkg r phase_offset]
function result=Gaussian_Hex(xp,xdata)
a = size(xdata, 1);
% xx=xdata(:,1:b/2);
% yy=xdata(:,b/2+1:b);
phase_offset = xp(12);
imgsize = a;
x0  = xp(1);
y0  = xp(2);
% r   = 5;
std_psf = xp(3);
intlist = xp(4:9);
phaselist = ([30, 90, 150, 210, 270, 330] + phase_offset)/180*pi;
bkg = xp(10);

r = xp(11);

result = GenHexPSF(imgsize, x0, y0, r, std_psf, intlist, phaselist);
result =result + bkg;
