function [dx, parout] = ode_2bp(t,x,mu,typeSim,varargin)
% ode_2bp - Ode function to compute the motion of the satellite.
%
% PROTOTYPE
%   dx=ode_2bp(t,x,mu,'cart') will perform a computation of the orbit
%   without perturbations using cartesian derivatives equations.
%
%   dx=ode_2bp(t,x,mu,'gauss') will perform a computation of the orbit
%   without perturbations using Gauss derivatives equations.
%
%   dx=ode_2bp(t,x,mu,typeSim,date0) will perform a computation of the
%   orbit with J2 and Moon perturbations using the specified derivatives
%   equations in typeSim input.
%
% INPUT:
%   t        double  [1x1]   time                                 [s]
%   x        double  [6x1]   state vector               [km and km/s]
%   mu       double  [1x1]   gravitational parameter       [km^3/s^2]
%   typeSim  char    [1x1]   type of simulation                   [-]
%
% OUTPUT:
%   dx       double [6x1]   state vector derivative               [-]
%
% CALLED FUNCTIONS: -
%
% NOTE: time can be omitted
%
% CONTRIBUTORS:
%   Rosato Davide               10618468
%   Saba Mohammadi Yengeje      10789462
%   Spinelli jason              10618465
%   Tagliati Alessia            10635119
%
% VERSIONS
%   2021-10-21: Release
%
% -------------------------------------------------------------------------

J2 = astroConstants(9);
Re = astroConstants(23);
            
switch typeSim
    case 'cart' % here for cartesian equations
        r = norm(x(1:3));

        if not(isempty(varargin)) % here for perturbation allowed
            % Modified Julian date 2000
            date0 = varargin{1};
            [Y, Mo, D] = ymd(date0);
            [H, M, S] = hms(date0);
            S = S + t;
            date = datetime([Y, Mo, D, H, M, S]);
            [Y, Mo, D] = ymd(date);
            [H, M, S] = hms(date);
            date0mjd2000 = date2mjd2000([Y, Mo, D, H, M, S]);

            % Position of the Moon
            [r_Moon, ~] = ephMoon(date0mjd2000);
            % Moon perturbation
            aMOON = MoonPer(r_Moon', x(1:3));
            % J2 perturbation
            aJ2 = J2Pert(x(1:3), J2, Re, mu);

            %%% COMPUTE STATE VECTOR DERIVATIVE
            dx = [x(4:6); -mu/(r^3)*x(1:3) + aJ2 + aMOON];
            
            %%% SAVE IMPORTANT DATA
            parout.aMOON = aMOON;
            parout.aJ2 = aJ2;

        else

            %%% COMPUTE STATE VECTOR DERIVATIVE
            dx = [x(4:6); -mu/(r^3)*x(1:3)];
            
            %%% SAVE IMPORTANT DATA
            parout.aMOON = 0;
            parout.aJ2 = 0;

        end
        
    case 'gauss' % here for Gauss equations
        a = x(1);
        e = x(2);
        i = x(3);
        % OM = x(4);
        om = x(5);
        th = x(6);

        b = a*sqrt(1-e^2);
        p = b^2/a;
        n = sqrt(mu/a^3);
        h = n*a*b;
        r = p/(1+e*cos(th));
        v = sqrt((2*mu)/r-mu/a);
        th_star = th+om;

        [rr, vv] = par2car(x, mu);
        
        if not(isempty(varargin)) % here for perturbation allowed
            % Modified Julian date 2000
            date0 = varargin{1};
            [Y, Mo, D] = ymd(date0);
            [H, M, S] = hms(date0);
            S = S + t;
            date = datetime([Y, Mo, D, H, M, S]);
            [Y, Mo, D] = ymd(date);
            [H, M, S] = hms(date);
            date0mjd2000 = date2mjd2000([Y, Mo, D, H, M, S]);
            
            % Position of the Moon
            [r_Moon, ~] = ephMoon(date0mjd2000);
            % Moon perturbation
            aMOON = MoonPer(r_Moon', rr);
            % J2 perturbation
            aJ2 = J2Pert(rr, J2, Re, mu);
            ap_car = aJ2 + aMOON;
            
            %%% SAVE IMPORTANT DATA
            parout.aMOON = aMOON;
            parout.aJ2 = aJ2;
            
        else
            
            ap_car = 0;
            
            %%% SAVE IMPORTANT DATA
            parout.aMOON = 0;
            parout.aJ2 = 0;

        end

        % "Cartesian" to TNH reference frame rotation matrix
        t_vers = vv/norm(vv);
        h_vers = cross(rr,vv)/norm(cross(rr,vv));
        n_vers = cross(h_vers,t_vers);
        A = [t_vers, n_vers, h_vers]';
        
        ap_tnh = A*ap_car;

        da = (2*a^2*v)/mu * ap_tnh(1);
        de = 1/v * (2*(e+cos(th))*ap_tnh(1) - r/a*sin(th)*ap_tnh(2));
        di = r*cos(th_star)/h*ap_tnh(3);
        dOM = r*sin(th_star)/(h*sin(i))*ap_tnh(3);
        dom = 1/(e*v)*(2*sin(th)*ap_tnh(1)+(2*e+r/a*cos(th))*ap_tnh(2))- ...
            (r*sin(th_star)*cos(i))/(h*sin(i))*ap_tnh(3);
        dth = h/r^2-1/(e*v)*(2*sin(th)*ap_tnh(1)+(2*e+r/a*cos(th))*ap_tnh(2));

        dx = [da, de, di, dOM, dom, dth]';
end


end