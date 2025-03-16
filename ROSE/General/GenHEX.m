function img = GenHEX(imgsize, parameters, modulationdepth, PhaseWidth, PhaseOffset)
%Hex Gaussian PSF Generation
%img = GenHEX(imgsize, parameters, modulationdepth, PhaseWidth, PhaseOffset)
%parameter
%[x0 y0 std int1 int2 int3 int4 int5 int6 bkg rx ry phase_offset]
% 1  2   3   4    5     6   7    8    9    10 11 12   13
%[x0 y0 std int phasex phasey rx (ry)], phase_offset = 15
% 1  2   3   4    5    6     7     8  
%
% PhaseOffset in degree

if nargin <3
    modulationdepth = 1;
end

if nargin <5
    PhaseOffset = 120;
end

if length(imgsize) ==1
    imgsize = [imgsize imgsize];
end

img = zeros(imgsize);
if length(parameters) == 13
    lp = parameters;
elseif length(parameters) == 8
    x0 = parameters(1);
    y0 = parameters(2);
    psf_std = parameters(3);
    psf_int = parameters(4);
    bkg=0;
    phasex = parameters(5); 
    phasey = parameters(6); 
    rx = parameters(7);
    ry = parameters(8);
    
    phaselist = [phasex-PhaseOffset, phasex, phasex+PhaseOffset, phasey-PhaseOffset, phasey, phasey+PhaseOffset]/180*pi;
    intlist = sin(phaselist).*modulationdepth + 1;
    intlist = intlist / sum(intlist)*psf_int;
%   [x0 y0 std int1 int2 int3 int4 int5 int6 bkg rx ry phase_offset]
%    1  2   3   4    5     6   7    8    9    10 11 12   13
    lp = [x0, y0, psf_std, intlist, bkg, rx, ry, 15];
    
elseif length(parameters) == 7
    x0 = parameters(1);
    y0 = parameters(2);
    psf_std = parameters(3);
    psf_int = parameters(4);
    bkg = 0;
    phasex = parameters(5); 
    phasey = parameters(6); 
    rx = parameters(7);
    ry = rx;
    
    phaselist = [phasex-PhaseOffset, phasex, phasex+PhaseOffset, phasey-PhaseOffset, phasey, phasey+PhaseOffset]/180*pi;
    intlist = sin(phaselist).*modulationdepth + 1;
    intlist = intlist / sum(intlist)*psf_int;
%   [x0 y0 std int1 int2 int3 int4 int5 int6 bkg rx ry phase_offset]
%    1  2   3   4    5     6   7    8    9    10 11 12   13
    lp = [x0, y0, psf_std, intlist, bkg, rx, ry, 15];
end

if nargin <4 || isempty(PhaseWidth)
    img = Gaussian_Hex_Oval(lp,img);
else
    img = Gaussian_Hex_Oval(lp,img, PhaseWidth);
end