function [a,e,i,OM,om,th] = car2kep(rr,vv,mu)
% car2kep starting from the vector of carthesian data provides keplerian
% paramethers
%
% PROTOTYPE:
%   [a,e,i,OM,om,th] = car2kep(rr,vv,mu)
%
% INPUT:
%   rr[3x1]         Position vector [km]
%   vv[3x1]         Velocity vector [km/s]
%   mu[1]           Gravitational constant [km^3/s^2]
%
% OUTPUT :
%   a[1]            Semimajor axis [km]
%   e[1]            Eccentricity [-]
%   i[1]            Inclination [rad]
%   OM[1]          	Right ascension of ascending node [rad]
%   om[1]        	Argument of perigee [rad]
%   th[1]           True anomaly [rad]
%
% CONTRIBUTORS:
%   Fabio Spada
%   Alessandro Staffolani
%   Ciro Salvi
%
% VERSIONS:
%   2021-02-11
%

r = norm(rr);
v = norm(vv);
K = [0 0 1];

a = mu*r/(2*mu-r*v^2); %valore semiasse maggiore

H = cross(rr,vv);
h = norm(H);
ecc = ((cross((vv),H))/mu)-((rr)/r); %vettore eccentricità
e = norm(ecc);

i = acos(H(3)/h);
di = rad2deg(i);

%OM
N = (cross(K,H));
N = N/norm(N);
if N == [0 0 0];
    OM = 0;
else
    if (N(2) >= 0 )
        OM = acos(N(1));
    else
        OM = 2*pi - acos(N(1));
    end
end
dOM = rad2deg(OM);

% om
if N == [0 0 0]
    om = atan2(ecc(2)/e, ecc(1)/e);
else
    
    if (ecc(3) >= 0 )
        om = acos((N*ecc)/e);
    else
        om = 2*pi - acos(N*ecc/e);
    end
end
dom = rad2deg(om);

%theta
if ((vv)'*(rr) >= 0 )
    th = acos(((rr)'*ecc)/(r*e));
else
    th = 2*pi - acos(((rr)'*ecc)/(r*e));
end
dth = rad2deg(th);

end





