function [a,e,i,Om,om,theta] = car2kep(rr,vv,mu)
% PROTOTYPE:
%   [a,e,i,Om,om,theta] = car2kep(rr,vv,mu)
% 
% DESCRIPTION:
%   Returns Keplerian orbital elements given the s/c position and velocity
%   vectors in cartesian coordinates and the planetary constant mu
% 
% INPUT:
% rr     position vector in cartesian coordinates
% vv     velocity vector in cartesian coordinates
% mu     planetary constant
% 
% OUTPUT:
% a     semi-major axis
% e     eccentricity
% i     inclination
% Om     RAAN (Right Ascension of the Ascending Node)
% om     periapsis true anomaly
% theta     true anomaly
% 
% CALLED FUNCTIONS:
%   (none)

r = norm(rr);
v = norm(vv);
k = [0 0 1]';

a = -mu / (v^2 - 2*mu/r);

hh = cross(rr,vv);
h = norm(hh);

ee = cross(vv,hh)/mu - rr/r;
e = norm(ee);

i = acos(hh(3)/h);

N = cross(k,hh) / norm(cross(hh,k));

Om = acos(N(1));
if N(2) < 0
    Om = 2*pi - Om;
end

om = acos(dot(N,ee)/e);
if ee(3) < 0
    om = 2*pi - om;
end

theta = acos(dot(rr,ee)/(r*e));
if dot(rr,vv) < 0
    theta = 2*pi - theta;
end




end

