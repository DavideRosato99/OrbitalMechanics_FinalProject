function [method, ao, con41, con42, cosio, cosio2, eccsq, omeosq,...
    posq, rp, rteosq, sinio, gsto, no_unkozai] = initl(xke, j2, ecco,...
    epoch, inclo, no_kozai)
% initl - the function initializes the spg4 propagator. all the initialization is
%     consolidated here instead of having multiple loops inside other routines.
%
% INPUTS:
%     ecco        - eccentricity                           [0.0 - 1.0]
%     epoch       - epoch time in days from jan 0, 1950. 0 hr
%     inclo       - inclination of satellite
%     no          - mean motion of satellite
%     satn        - satellite number
%
% OUTPUTS:
%     ainv        - 1.0 / a
%     ao          - semi major axis
%     con41       -
%     con42       - 1.0 - 5.0 cos(i)
%     cosio       - cosine of inclination
%     cosio2      - cosio squared
%     einv        - 1.0 / e
%     eccsq       - eccentricity squared
%     method      - flag for deep space                    'd', 'n'
%     omeosq      - 1.0 - ecco * ecco
%     posq        - semi-parameter squared
%     rp          - radius of perigee
%     rteosq      - square root of (1.0 - ecco*ecco)
%     sinio       - sine of inclination
%     gsto        - gst at time of observation               [rad]
%     no          - mean motion of satellite
%
% LOCALS:
%     ak          -
%     d1          -
%     del         -
%     adel        -
%     po          -
%
% COUPLING:
%     gstime      - find greenwich sidereal time from the julian date
%
% AUTHOR: 
%   Jeff Beck 
%   beckja@alumni.lehigh.edu
%   1.0 (aug 7, 2006) - update for paper dav
%   1.1 (nov 16, 2007)- update for better compliance
%
% ORIGINAL COMMENTS FROM VALLADO C++ VERSION:
%   David Vallado 719-573-2600   28 jun 2005
%
% REFERENCES:
%     hoots, roehrich, norad spacetrack report #3 1980
%     hoots, norad spacetrack report #6 1986
%     hoots, schumacher and glover 2004
%     vallado, crawford, hujsak, kelso  2006
%
% CHANGELOG:
%     2021-12-17, 2021-2022 Assignments changes for practical uses
%
% CONTRIBUTORS:
%   Rosato Davide               10618468
%   Saba Mohammadi Yengeje      10789462
%   Spinelli Jason              10618465
%   Tagliati Alessia            10635119
% ------------------------------------------------------------------------
% WGS-84 earth constants
x2o3 = 2/3;

% Calculate auxillary epoch quantities
eccsq  = ecco * ecco;
omeosq = 1.0 - eccsq;
rteosq = sqrt(omeosq);
cosio  = cos(inclo);
cosio2 = cosio^2;

% Un-kozai the mean motion
ak           = (xke / no_kozai)^x2o3;
d1          = 0.75 * j2 * (3.0 * cosio2 - 1.0) / (rteosq * omeosq);
del         = d1 / (ak * ak);
adel        = ak * (1.0 - del * del - del * (1.0 / 3.0 + 134.0 * del * del / 81.0));
del         = d1/(adel * adel);
no_unkozai  = no_kozai / (1.0 + del);

ao     = (xke / no_unkozai)^x2o3;
sinio  = sin(inclo);
po     = ao * omeosq;
con42  = 1.0 - 5.0 * cosio2;
con41  = -con42-cosio2-cosio2;
posq   = po * po;
rp     = ao * (1.0 - ecco);
method = 'n';

% SGP4 modern approach to finding sidereal time
gsto = mjd20002gmst(epoch); 
end