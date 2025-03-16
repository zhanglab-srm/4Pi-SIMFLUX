%Hex Gaussian PSF Generation
%parameter
%[x0 y0 std int1 int2 int3 int4 int5 int6 bkg rx ry phase_offset dx1 dx2 dx3 dx4 dx5 dx6 dy1 dy2 dy3 dy4 dy5 dy6]
% 1  2   3   4    5     6    7    8   9    10  11 12 13           14  15  16  17  18  19  20  21  22  23  24  25                                  
function result=Gaussian_Hex_Oval_dpos(xp,xdata)
a = size(xdata, 1);
% xx=xdata(:,1:b/2);
% yy=xdata(:,b/2+1:b);
phase_offset = xp(13);
imgsize = a;
x0  = xp(1);
y0  = xp(2);
% r   = 5;
std_psf = xp(3);
intlist = xp(4:9);
phaselist = ([30, 90, 150, 210, 270, 330] + phase_offset)/180*pi;
bkg = xp(10);

rx = xp(11);
ry = xp(12);

result = GenHexPSF_Oval_dpos(imgsize, x0, y0, rx, ry, std_psf, intlist, phaselist, xp(14:19), xp(20:25), (-10:10)./180.*pi);
result =result + bkg;
