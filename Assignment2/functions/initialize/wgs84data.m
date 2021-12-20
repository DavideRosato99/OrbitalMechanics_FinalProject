function [tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2] = wgs84data
% wgs84data - the function gets constants for the propagator. note that mu 
%             is identified to facilitiate comparisons with newer models.
%
% INPUT: [-]
%    It uses by default wgs-84 constants
%
% OUTPUT:
%    tumin       - minutes in one time unit
%    mu          - earth gravitational parameter
%    radiusearthkm - radius of the earth in km
%    xke         - reciprocal of tumin
%    j2, j3, j4  - un-normalized zonal harmonic values
%    j3oj2       - j3 divided by j2
%
% AUTHOR:
%    David Vallado 719-573-2600   21 jul 2006
%
% REFERENCES:
%    norad spacetrack report #3
%    vallado, crawford, hujsak, kelso  2006
%    [tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2] = getgravc(whichconst);
%
% CHANGELOG:
%     2021-12-17, 2021-2022 Assignments changes for practical uses
%
% CONTRIBUTORS:
%   Rosato Davide               10618468
%   Saba Mohammadi Yengeje      10789462
%   Spinelli Jason              10618465
%   Tagliati Alessia            10635119
%
%  ------------------------------------------------------------------------


mu            = astroConstants(13);
radiusearthkm = astroConstants(23);
xke           = 60.0 / sqrt((radiusearthkm^3)/mu);
tumin         = 1.0 / xke;
j2            =   0.00108262998905;
j3            =  -0.00000253215306;
j4            =  -0.00000161098761;
j3oj2         =  j3 / j2;

end