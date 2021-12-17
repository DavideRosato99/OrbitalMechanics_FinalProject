function GMST = mjd20002gmst(MJD2000)
% mjd20002gmst -  Convert a specified Julian Date Vector to Greenwhich Mean
%                 Sidereal Time.
%
% INPUT:
%   mjd2000     Date in MJD 2000. MJD2000 is defined as the number of days
%               since 01-01-2000, 12:00 noon.
%
% OUTPUT:
%   GMST        Greenwich Mean Sidereal Time in degrees from 0-360
%
% CALLED FUNCTIONS: -
%
% CONTRIBUTORS:
%   Rosato Davide               10618468
%   Saba Mohammadi Yengeje      10789462
%   Spinelli Jason              10618465
%   Tagliati Alessia            10635119
%
%
% REFERENCES:
%   Approximate Sidereal Time: http://www.usno.navy.mil/USNO/astronomica...
%   l-applications/astronomical-information-center/approx-sider-time
%
%   Universal Sidereal Times, The Astronomical Almanac For The Year 2004
%
%   Darin Koblick
%
% VERSIONS:
%    2010-11-07
%--------------------------------------------------------------------------

JD = mjd20002jd(MJD2000);
%Find the Julian Date of the previous midnight, JD0
JD0 = NaN(size(JD));
JDmin = floor(JD)-.5;
JDmax = floor(JD)+.5;
JD0(JD > JDmin) = JDmin(JD > JDmin);
JD0(JD > JDmax) = JDmax(JD > JDmax);
H = (JD-JD0).*24;       %Time in hours past previous midnight
D = JD - 2451545.0;     %Compute the number of days since J2000
D0 = JD0 - 2451545.0;   %Compute the number of days since J2000
T = D./36525;           %Compute the number of centuries since J2000
%Calculate GMST in hours (0h to 24h) ... then convert to degrees
GMST = mod(6.697374558 + 0.06570982441908.*D0  + 1.00273790935.*H + ...
    0.000026.*(T.^2),24).*15;

GMST = GMST * pi/180;