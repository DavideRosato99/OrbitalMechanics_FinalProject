function [rrmat, vvmat,kep] = kep2car(varargin)
% kep2car transformation from keplerian coordinates to cartesian
%
% PROTOTYPE:
%   [rrmat, vvmat,kep] = kep2car(varargin)
%
% INPUT:
%   a[1]            Semimajor axis [km]
%   e[1]            Eccentricity [-]
%   i[1]            Inclination [rad]
%   OM[1]          	Right ascension of ascending node [rad]
%   om[1]        	Argument of perigee [rad]
%   th[1]           True anomaly [rad]
%   mu[1]           Planetary constant [km^3/s^2]
%   type[1]         Unit of measure of input angles [ 'deg' || 'rad' ]
%or
%	kep[7x1]        Vector of keplerian elements [ a, e, i, OM, om, thvect, mu]
%	type[1]         Unit of measure of input angles [ 'deg' || 'rad' ]
%
%*if the unit of measure of input angles is rad it's not necessary type
%as an input
%
% OUTPUT :
%   rrmat[3x1]      Position components for each theta evaluated [km]
%   vvmat[3x1]      Velocity components for each theta evaluated [km/s]
%   kep [3x1]       Vector of keplerian elements [ a, e, i, OM, om, thvect, mu]
%
% CONTRIBUTORS:
%   Alessandro Staffolani
%   Fabio Spada
%   Ciro Salvi
%
% VERSIONS
%   2021-02-11
%

if( nargin <= 2 )
    kep = cell2mat(varargin(1));
    a       = kep(1);
    e       = kep(2);
    i       = kep(3);
    OM      = kep(4);
    om      = kep(5);
    thvect  = kep(6);
    mu      = kep(7);
    if( nargin == 2 )
        type  = varargin(2);
        if ( strcmp(type, 'deg') )
            i = i * pi/180;
            OM = OM * pi/180;
            om = om * pi/180;
            thvect = thvect * pi/180;
        end
    end
elseif( nargin>2 && nargin<=8 )
    a       = cell2mat(varargin(1));
    e       = cell2mat(varargin(2));
    i       = cell2mat(varargin(3));
    OM      = cell2mat(varargin(4));
    om      = cell2mat(varargin(5));
    thvect  = cell2mat(varargin(6));
    mu      = cell2mat(varargin(7));
    if( nargin == 8 )
        type  = varargin(8);
        if ( strcmp(type, 'deg') )
            i = i * pi/180;
            OM = OM * pi/180;
            om = om * pi/180;
            thvect = thvect * pi/180;
        end
    end
end

p = a * (1 - e^2);

ROM =  [cos(OM), -sin(OM),       0;
        sin(OM),  cos(OM),       0;
              0,        0,       1];

Ri =   [      1,        0,       0;
              0,   cos(i), -sin(i);
              0,   sin(i),  cos(i)];

Rom =  [cos(om), -sin(om),       0;
        sin(om),  cos(om),       0;
              0,        0,       1];

R313 = ROM * Ri * Rom ;

rrmat = [];
vvmat = [];

L = length(thvect);

for ii = 1 : L
    
    th = thvect(ii);
    r  = (a*(1 - e^2))/(1 + e*cos(th));
    rvect = r*[cos(th); sin(th);0];
    vvect = sqrt(mu/p) * [-sin(th); e + cos(th); 0];
    
    rr = R313 * rvect;
    vv = R313 * vvect;
    
    rrmat = [rrmat, rr];
    vvmat = [vvmat, vv];
end

end
