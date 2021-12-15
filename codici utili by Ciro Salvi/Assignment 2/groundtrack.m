%fg0 = theta Greenwich at t0
%we  = angular velocity Earth
%tspan = vector of integration times

function [alpha, delta, lat, long,a,Y,S,tspan,T] = groundtrack(a,e,i,OM,om,f0,mu,we,tspan,fg0,time_evaluation,k,m)

% 
% 
% PROTOTYPE:
%   [alpha, delta, lat, long,a,Y,S,tspan,T] = groundtrack(a,e,i,OM,om,f0,mu,we,tspan,fg0,time_evaluation,k,m)

% INPUT
%   a,e,i,OM,om,f0,fg0: initial keplerian elements [km,adimensional,rad,rad,rad,rad,rad]
%   mu: gravitational constant km^3/s^2
%   we: earth rotation deg/h
%   time_evaluation: struct set in the main
%   k,m: parameters of the repeating GT
% 
% OUTPUT
% 	[alpha, delta] latitude and longitude of the GT
%   [lat, long] latitude and longitude of the GT
%   a0: semimajor axis for repeating, km
%   Y: vector of the propagation, position (km), velocity(km/s)
%   tspan: time vector of the propagation
%   T: period of the repeating orbit (s)

% 
% CONTRIBUTORS
%    Ciro Salvi
%    Alessandro Staffolani 
%
% VERSIONS
% 2020-02-11


%we is in degree/h
we = deg2rad(we/3600);

if (nargin > 12) %for a repeating ground track with a given ratio k/m, k = number of satellite's periods, m = number of Earth's revolutions
    n = we*(k/m);
    a =(mu/(n^2))^(1/3);
end
    T = 2*pi*sqrt(a^3/mu);
if( nargin >= 11 && strcmp(time_evaluation, 'one-revolution') )
    tspan = linspace(0,T,10000);
elseif( nargin >= 11 && strcmp(time_evaluation, '12-revolution'))
    tspan = linspace(0,12*T,10000);
end
    
[rr, vv] = par2car(a,e,i,OM,om,f0,mu);
y0 = [rr; vv];

% Set options
options = odeset ( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

[~, Y] = ode113(@(t,y) ode_orbit(t,y,mu), tspan, y0, options);

R = [Y(:,1)' ;Y(:,2)' ;Y(:,3)'];
V = [Y(:,4)' ;Y(:,5)' ;Y(:,6)'];

for ii = 1:size(Y,1)
    r_mag(ii) = norm(R(:,ii));
    delta(ii) = asin(R(3,ii)/r_mag(ii));
    if(R(2,ii)/r_mag(ii) > 0)
        alpha(ii) = acos((R(1,ii)/r_mag(ii))/(cos(delta(ii))));
    else
        alpha(ii) = 2*pi - acos((R(1,ii)/r_mag(ii))/(cos(delta(ii))));
    end
    %longitude
    fg(ii) = fg0 + we*tspan(ii); %angle of Earth's rotation refere to Greenwich
    long(ii) = wrapToPi(alpha(ii)-fg(ii)); % vector of longitude
    if(ii~=1 && (isnan(long(ii-1)) == 0) && (sign(long(ii)) ~= sign(long(ii-1))) && (round(long(ii)) ~= 0))
        long(ii) = NaN;
    end    
    %latitude
    lat(ii) = delta(ii); % vector of latitude
    %x-y-z projected on the Earth's surface (spherical coordinates)
    R_e     = 6371.0087714; % equatorial radius (meters)
    S(ii,1) = R_e*cos(lat(ii))*cos(long(ii));
    S(ii,2) = R_e*cos(lat(ii))*sin(long(ii));
    S(ii,3) = R_e*sin(lat(ii));
end
tspan1 = tspan;
end











