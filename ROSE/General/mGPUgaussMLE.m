%mGPUgaussMLE  GPU based MLE of single molecule position, emission rate and background
%
%   [X Y N BG S CRLBx CRLBy CRLBn CRLBb CRLBs LogL]=GPUgaussMLE(data,PSFSigma,iterations,fittype)
%
%   This code performs a maximum likelihood estimate of particle position,
%   emission rate (photons/frame), and background rate (photons/pixels/frame).
%   The found values are used to calculate the Cramer-Rao Lower Bound (CRLB)
%   for each parameter and the CRLBs are returned along with the estimated parameters.
%   The input is one or more identically sized images, each of which is
%   assumed to contain a single fluorophore.  Since the maximum likelihood
%   estimate assumes Poisson statistics, images must be converted from arbitrary
%   ADC units to photon counts.
%
%       INPUTS:
%   data:       A 3D stack of images
%   PSFSigma:   microscope PSF sigma (sigma_0 for z fit, starting sigma for sigma fit)
%   iterations: number of iterations (default=10)
%   fittype:    1: position,bg,N only; 2: also PSF sigma; 3: z position;
%               4: PSF sigma_x,sigma_y (CRLBx not returned). Default=1
%
%       OUTPUTS:
%   X:      found X positions
%   Y:      found Y positions
%   N:      found Photons
%   BG:     found Background rate
%   S:      found Sigma
%   CRLBx:  x-position uncertainty
%   CRLBy:  y-position uncertainty
%   CRLBn:  Photons uncertainty
%   CRLBb:  background rate uncertainty
%   CRLBs:  sigma or z-position uncertainty
%   LogL:   LogLikelihood/Pixels

%   REFERENCE:
%   "Fast, single-molecule localization that
%   achieves theoretical optimal accuracy." Carlas S. Smith, Nikolai Joseph,
%   Bernd Rieger and Keith A. Lidke

function [X Y N BG S CRLBx CRLBy CRLBn CRLBb CRLBs LogL]=mGPUgaussMLE(data,PSFSigma,iterations,fittype)

if nargin<3
    iterations=10;
end
if nargin<4
    fittype=1;
end
if nargin<2
   error('Minial usage: mGPUgaussMLE(data,PSFSigma)');
end

%convert dipimage data to single
if isa(data,'dip_image')
    data=permute(single(data),[2 1 3]);
end

%convert other matlab datatyes to single
if ~ (isa(data,'double'))
    data=single(data);
end

if ~(isa(data,'single'))
    error('Input data must be type dipimage, single, or double')
end


[X Y N BG S CRLBx CRLBy CRLBn CRLBb CRLBs LogL]=GPUgaussMLE(data,PSFSigma,iterations,fittype);

%convert variances to standard deviations
CRLBx=sqrt(CRLBx);
CRLBy=sqrt(CRLBy);
CRLBn=sqrt(CRLBn);
CRLBb=sqrt(CRLBb);
CRLBs=sqrt(CRLBs);


end











