function [ DFT, interpSignal ] = generalSignalDFT( xVect, yVect, type )
% generalSignalDFT generalSignalDft evaluates the dft of a generic signal,
% even in non equally spaced form, performing an equally spaced 
% interpolation of signal data, and then using fft matlab built in
% function. Interpolation is performed using the interpolation type fed as input.
% Hanning windowing is further added to minimize leakage problematics 
% WARNING: generalSignalDFT filters out the zero frequency component
%
% PROTOTYPE:
%	[DFT, interpSignal] = generalSignalDFT(xVect, yVect, type)
%
% INPUT:
%   xVect[N]                Vector of time instants [s]
%	yVect[N]        		Signal values at time instants of xVect [-]
%	type[1]                 Interpolation setting, fed to interp1.m[string]
%
% OUTPUT :
%   DFT[struct]           	
%       wholeDom[N*900]         Full frequency domain of 1Sided transform [Hz]
%       whole1SideTran[N*900]   1 Sided transform of the signal [-]
%       redDom[-]               Reduced domain of the most valuable signals [Hz]
%       red1SideTran[-]         Reduced 1 sided transform [-]
%   interpSignal[N*900]         Rebuilt signal via interpolation of yVect, 
%                               following 'type' interpolation type [-]
%
% CONTRIBUTORS:
%   Fabio Spada
%
% VERSIONS
%   2021-02-11
%

xVectInterp = linspace(xVect(1), xVect(end), 900*length(xVect));
yVectInterp = interp1(xVect, yVect, xVectInterp, type);

L = length(xVectInterp);

interpSignal = [xVectInterp',yVectInterp'];

% zero frequency component elimination
yVectInterpFFT = yVectInterp - mean(yVectInterp);
hanningW = hanning(length(yVectInterpFFT));
hanningW = hanningW';
yVectInterpFFT = hanningW.*yVectInterpFFT;
yFFT = fft(yVectInterpFFT);

% inizialization of multiple dft represantations
% two sides dft
yDFT2 = abs(yFFT/L);
% one sided dft
yDFT1 = yDFT2(1: L/2+1);
yDFT1(2:end-1) = 2*yDFT1(2:end-1);

% inizialization of frequency discretized domain
% frequency resolution
fRes = 1/(xVectInterp(2) - xVectInterp(1));
% frequency vector
fVect = fRes*(0:(L/2))/L;


stem(fVect, yDFT1);
grid minor

% find the high-contribution-associated indexes
[~, ind] = find( yDFT1 > max(yDFT1)/20 );

% fixing the proper scale to represent the frequency content
xlim([0, fVect(ceil(length(fVect)/256))]);

fSig = fVect(ind);
ySig = yDFT1(ind);

DFT.wholeDom        = fVect;
DFT.whole1SideTran  = yDFT1;
DFT.redDom          = fSig; 
DFT.red1SideTran    = ySig;

end