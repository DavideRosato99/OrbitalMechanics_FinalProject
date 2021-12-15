function Y = filter_mov(X,T, kep,mu,n) 
% 
% 
% PROTOTYPE:
%   Y = filter_mov(X,T, kep,mu,n) 
% 
% INPUT
%   X[kk] generic vector of dimension kk defined with respect to a T vector
%   T[kk]
%   kep[6x1] keplerian elements [km,adimensional,rad,rad,rad,rad]
%   mu[1] gravitational constant, km^3/s^2
%   n[1] number of orbital periods on which the filtering is performed
% 
% OUTPUT
% Y[kk] filtered vector
% 
% CONTRIBUTORS
% Ciro Salvi
% 
% VERSIONS
% 2020-02-11

Y = movmean(X,n*2*pi*sqrt(kep(1)^3/mu),'SamplePoints',T,'endpoints', 'shrink');
end