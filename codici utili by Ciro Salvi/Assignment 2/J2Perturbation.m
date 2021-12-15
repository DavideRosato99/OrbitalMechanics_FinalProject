function accelPert = J2Perturbation(Kep,RS, Att,rr,vv)

% 
% PROTOTYPE:
%   accelPert = J2Perturbation(Kep,RS, Att,rr,vv)
% 
% INPUT
%   Kep [6x1] keplerian elements [km,adimensional,rad,rad,rad,rad]
%   RS {struct} defined the reference frame 'tnh' or 'rsw'
%   rr[3x1] position vector, km 
%   vv[3x1] velocity vector, km/s
%   mu[1] gravitational constant, km^3/s^2
% 
% OUTPUT
% 	accelpert [3x1] acceleration given by the perturbation, km/s^2
% 
% 
% CONTRIBUTORS
%    Ciro Salvi
%    Fabio Spada
% 
% VERSIONS
% 2020-02-11

mu = astroConstants(13);
r = norm(rr);
v = norm(vv);
accelPert = [];
const = -1.5*astroConstants(9)*astroConstants(23)^2*mu/r^4;
h = sqrt(mu*Kep(1)*(1-Kep(2)^2));

    switch RS
        case 'rsw'
            accelPert = const * [1 - 3*sin(Kep(3))^2*sin(Kep(6) + Kep(5))^2;
                                sin(Kep(3))^2*sin(2*(Kep(6) + Kep(5)));
                                sin(2*Kep(3))*sin(Kep(6) + Kep(5))];
        case 'tnh'
            mat = [Kep(2)*sin(Kep(6)), -1-Kep(2)*cos(Kep(6));1+Kep(2)*cos(Kep(6)), Kep(2)*sin(Kep(6))]*h/v/Kep(1)/(1-Kep(2)^2);
            accelPer = const * [1 - 3*sin(Kep(3))^2*sin(Kep(6) + Kep(5))^2;
                                sin(Kep(3))^2*sin(2*(Kep(6) + Kep(5)))];
            accelPert = [mat\accelPer ; const * sin(2*Kep(3))*sin(Kep(6) + Kep(5))];

    end

end