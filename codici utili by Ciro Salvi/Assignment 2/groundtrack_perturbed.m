
function [alpha, delta, lat, long,a0,Y,S,tspan,T] = groundtrack_perturbed(a,e,i,OM,om,f0,mu,we,tspan,fg0,time_evaluation,J_2,R_e,k,m)
% 
% 
% PROTOTYPE:
%   [alpha, delta, lat, long,a0,Y,S,tspan,T] = groundtrack_perturbed(a,e,i,OM,om,f0,mu,we,tspan,fg0,time_evaluation,J_2,R_e,k,m)

% INPUT
%   a,e,i,OM,om,f0,fg0: initial keplerian elements [km,adimensional,rad,rad,rad,rad,rad]
%   mu: gravitational constant, km^3/s/2
%   we: earth rotation, deh/h
%   time_evaluation: struct set in the main
%   J_2: J2 constant
%   R_e: radius of earth, km
%   k,m: parameters of the repeating GT, adimensional
% 
% OUTPUT
% 	[alpha, delta] latitude and longitude of the GT
%   [lat, long] latitude and longitude of the GT
%   a0: semimajor axis for repeating
%   Y: vector of the propagation
%   tspan: time vector of the propagation
%   T: period of the repeating orbit

% 
% CONTRIBUTORS
%    Ciro Salvi
%    Alessandro Staffolani 
%
% VERSIONS
% 2020-02-11

we = deg2rad(we/3600);
a0 = a;
    OM_punto = @(a) -((3/2)*((sqrt(mu)*J_2*R_e^2)/(((1-e^2)^2)*(a^(7/2)))))*cos(i);
    om_punto = @(a) -((3/2)*((sqrt(mu)*J_2*R_e^2)/(((1-e^2)^2)*(a^(7/2)))))*(-2 + (5/2)*(sin(i))^2);
    M_punto  = @(a) -((3/2)*((sqrt(mu)*J_2*R_e^2)/(((1-e^2)^2)*(a^(7/2)))))*(1- (3/2)*(sin(i))^2);
    n        = @(a) sqrt(mu/(a^3));
    
if (nargin > 13) %for a repeating ground track with a given ratio k/m, k = number of satellite's periods, m = number of Earth's revolutions
    
    f = @(a) ((m/k) - (we - OM_punto(a))/(n(a) + M_punto(a) + om_punto(a)));
    a = fzero(f,a0);        % first guess of a
    a0 = a;
    T_mean = 2*pi/(n(a)+M_punto(a)+om_punto(a));
    % approaches to find right initial conditions
    %DELTA a converges rapidly: in 3 iterations, delta = 1.06e-4 km 
    for kkkk = 1:10
        [a,delta_a,T_mean] = iterative3 (a,e,i,OM,om,f0,mu,a0,we,tspan,fg0,time_evaluation,J_2,R_e,k,m,T_mean);
    end
    T = T_mean;
% %     second approach find theta is slower: 9 iterations to get to a_mean=a0
%     for kkkk = 1:50
%         [a_mean,e_mean,i_mean,thetastart,delta_a] = iterative (a,e,i,OM,om,f0,thetastart,mu,we,tspan,fg0,time_evaluation,J_2,R_e,k,m);
%         VETMEAN = [VETMEAN;a_mean,e,i];
%         f0 = thetastart
%     end
    
else 
    om_p = om_punto(a); M_p = M_punto(a); OM_p = OM_punto(a); nn = n(a);
    T = 2*pi/(nn+M_p+om_p);
end
   
if( nargin >= 13 && strcmp(time_evaluation, 'one-revolution') )
    tspan = linspace(0,T,10000);
elseif( nargin >= 13 && strcmp(time_evaluation, '12-revolution'))
    tspan = linspace(0,12*T,10000);
end

[rr, vv] = par2car(a,e,i,OM,om,f0,mu);
y0 = [rr; vv];


options = odeset ( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
[time, Y] = ode45(@(t,y) ode_orbit_perturbated(t,y,mu,R_e,J_2), tspan, y0, options);

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
end




