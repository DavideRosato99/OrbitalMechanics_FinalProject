%rad
function [rrmat, vvmat] = par2car(a , e , i, OM, om, thvect, mu)

% 
% PROTOTYPE:
%   [rrmat, vvmat] = par2car(a , e , i, OM, om, thvect, mu)
% 
% INPUT
%   [a,e,i,OM,om,th] keplerian elements [km,adimensional,rad,rad,rad,rad]
%   mu[1] gravitational constant, km^3/s^2
% 
% OUTPUT
%   rrmat[3x1] positions vector, km
%   vvmat[3x1] velocities vector, km/s
% 	

% 
% CONTRIBUTORS
%   Alessandro Staffolani
%   Fabio Spada
%   Ciro Salvi
% 
% VERSIONS
% 2020-02-11

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
        r = (a*(1 - e^2))/(1 + e*cos(th));
        rvect = r*[cos(th); sin(th);0];
        vvect = sqrt(mu/p) * [-sin(th); e + cos(th); 0];
        
        rr = R313 * rvect;
        vv = R313 * vvect;
        
        rrmat = [rrmat, rr];
        vvmat = [vvmat, vv];
    end
    
end
        