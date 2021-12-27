function PlotOrbit (a,e,i,Om,om,mu,th0,thf,dth,type, plotoptions)
% 
% PlotOrbit plots a steady orbit, given its keplerian parameters and the range of trume
% anomalies.
% 
% PROTOTYPE:
%    PlotOrbit (a,e,i,Om,om,mu,th0,thf,dth,type, plotoptions)
%
% INPUT:
%   a[1]            Semimajor axis [km]
%   e[1]            Eccentricity [-]
%   i[1]            Inclination [°/rad]
%   Om[1]          	Right ascension of ascending node [°/rad]
%   om[1]        	Argument of perigee [°/rad]
%   mu[1]           Planetary constant [km^3/s^2]
%   th0[1]          Left bound true anomaly [°/rad]
%   thf[1]          Right bound true anomaly [°/rad]
%   dth[1]          True anomaly plotting step [°/rad]
%   type[1]         Unit of measure of input angles [ 'deg' || 'rad' ]
%   plotoptions[struct] Struct with plotting options
%       style[1]    LineStyle
%       lw [1]      LineWidth'Linewidth'plotoptions.lw, 
%       color [-]   Color options; possible to feed an hexadecimal representation either an RGB vector
%
% CONTRIBUTORS:
%   Fabio Spada
%   Alessandro Staffolani
%   Ciro Salvi
%
% VERSION:
%   2021-02-11
%


if nargin > 6
    th = th0 : dth : thf;
elseif nargin == 6
    th = linspace (0,2*pi,100);
end


[rr,vv] = kep2car(a,e,i,Om,om,th,mu, type);

if nargin > 10
    plot3(rr(1,:),rr(2,:),rr(3,:), plotoptions.style, 'Linewidth',plotoptions.lw, 'Color', plotoptions.color );
else
    plot3(rr(1,:),rr(2,:),rr(3,:), 'Linewidth',1.5 );
end

hold on

end
