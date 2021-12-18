function [rECI, vECI, tm] = teme2eci(rTEME, vTEME, timemjd2000)
% teme2eci - function to compute the rotation matrix to perform a
% transformation from TEME (True Equator Mean Equinox) reference frame to
% ECI Geocentric Equatorial Reference Frame (IAU-76/FK5 J2000, mean
% equator, mean equinox frame).
%
% PROTOTYPE
%   [rECI,vECI]=teme2eci(rTEME, vTEME, timemjd2000)
%
% INPUT:
%   rTEME        double  [Nx3]   positions vectors (TEME)            [km]
%   vTEME        double  [Nx3]   velocities vectors (TEME)         [km/s]
%   timemjd2000  double  [1x1]   modified Julian date (J2000)         [d]
%
% OUTPUT:
%   rECI         double  [Nx3]   positions vectors (ECI)             [km]
%   vECI         double  [Nx3]   velocities vectors (ECI)          [km/s]
%   tm           double  [3x3]   rotation matrix (TEME -> ECI)        [-]
%
% CALLED FUNCTIONS: precess, nutation
%
% REFERENCES:
%   David Vallado       2013, 231-233
%
% CONTRIBUTORS:
%   Rosato Davide               10618468
%   Saba Mohammadi Yengeje      10789462
%   Spinelli jason              10618465
%   Tagliati Alessia            10635119
%
% VERSIONS
%   2021-10-21: Release
%
% -------------------------------------------------------------------------

% Get size of the positions (and velocities) vector
N = size(rTEME, 1);

%% CALCULATE ROTATION MATRIX
ddpsi = -0.054522 * pi / (180*3600);
ddeps = -0.006209 * pi / (180*3600);
centjd = (timemjd2000)/(36525);

prec = precess(centjd);

[deltapsi, meaneps, nut] = nutation(centjd, ddpsi, ddeps);

eqeg = deltapsi* cos(meaneps);

eqeg = rem (eqeg, 2.0*pi);

eqe(1,1) =  cos(eqeg);
eqe(1,2) =  sin(eqeg);
eqe(1,3) =  0.0;
eqe(2,1) = -sin(eqeg);
eqe(2,2) =  cos(eqeg);
eqe(2,3) =  0.0;
eqe(3,1) =  0.0;
eqe(3,2) =  0.0;
eqe(3,3) =  1.0;

tm = prec * nut * eqe';

%% ROTATE VECTORS
rECI = zeros(N, 3);
vECI = zeros(N, 3);
for i = 1:N
    rECI(i, 1:3) = tm * rTEME(i, 1:3)';
    vECI(i, 1:3) = tm * vTEME(i, 1:3)';
end
