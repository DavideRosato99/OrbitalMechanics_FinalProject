function dKepdt = OdeKepR2BP(t, Kep, type, RS,am_ratio,cd,R_E,W_E)

% PROTOTYPE:
%   dKepdt = OdeKepR2BP(t, Kep, type, RS,am_ratio,cd,R_E,W_E)
% 
% INPUT
%   type [struct] defines the kind of perturbations considered,'j2','drag',
%   'global', 'off'
%   t comes from the integration with ode113 or ode45, s
%   Kep [6x1] keplerian elements [km,adimensional,rad,rad,rad,rad]
%   RS {struct} defined the reference frame
%   am_ratio [1] a:m ratio of satellite, m^2/kg
%   cd[1] dragcoefficient of satellite
%   R_E[1] radius of earth,km
%   W_E[1] angular velocity of earth, deg/h
% 
% OUTPUT
% 	dKepdt [6x1] first derivative of the state

% 
% CONTRIBUTORS
%    Ciro Salvi
%    Fabio Spada
% 
% VERSIONS
% 2020-02-11

% RS:   - 'rsw'     RadialTangentOutOfPlane    

    a = Kep(1);
    e = Kep(2);
    i = Kep(3);
    Om = Kep(4);
    om = Kep(5);
    th = Kep(6);
    mu = astroConstants(13);
    
    [rr,vv] = kep2car(Kep(1), Kep(2), Kep(3), Kep(4), Kep(5), Kep(6), mu, 'rad');
    r = norm(rr);
    v = norm(vv);
    p = a*(1-e^2);
    h = sqrt(p*mu);
    
    Accel = Kep2AccelPerturb( type, t, Kep,RS,am_ratio,cd,R_E,W_E,rr,vv); 
    
    switch RS
        
        case 'rsw'
            
            accR = Accel(1);
            accS = Accel(2);
            accW = Accel(3);
            
            dKepdt = [  2*a^2/h*(e*sin(th)*accR + p/r*accS);
                        1/h *(p*sin(th)*accR + ((p+r)*cos(th) + r*e)*accS);
                        r*cos(th+om)/h*accW;
                        r*sin(th+om)/(h*sin(i))*accW;
                        1/(h*e)*(-p*cos(th)*accR + (p+r)*sin(th)*accS) - r*sin(th+om)*cos(i)/(h*sin(i))*accW;
                        h/r^2 + 1/(e*h) * (p*cos(th)*accR - (p+r)*sin(th)*accS)];
                    
        case 'tnh'
            
            at = Accel(1);
            an = Accel(2);
            ah = Accel(3);
            
            dKepdt = [  2*a^2*v/mu * at,
                        (2*(e+cos(th))*at-r*sin(th)*an/a)/v,
                        r*ah/h*cos(th+om),
                        r*ah/h*sin(th+om)/sin(i),
                        (2*at*sin(th)+(2*e+r/a*cos(th))*an)/e/v - r*ah/h*sin(th+om)*cos(i)/sin(i),
                        h/r^2 - (2*at*sin(th) + (2*e+r/a*cos(th))*an)/e/v];
    
    end

end