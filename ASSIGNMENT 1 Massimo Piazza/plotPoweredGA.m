function [t_SOI, r_p, h_GA, DeltaV, e_minus, e_plus] = plotPoweredGA(fb_ID, mjd_fb, VV_minus, VV_plus, dt, zoomFactor)
% PROTOTYPE:
%  [t_SOI, r_p, h_GA, DeltaV, e_minus, e_plus] = plotPoweredGA(fb_ID, mjd_fb, VV_minus, VV_plus, dt, zoomFactor)
% 
% DESCRIPTION:
%   Returns relevant parameters of the Powered Gravity Assist and the plot
%   of the maneuver
% 
% INPUT:
%     fb_ID: Integer number identifying the celestial body (< 10)
%                   1:   Mercury
%                   2:   Venus
%                   3:   Earth
%                   4:   Mars
%                   5:   Jupiter
%                   6:   Saturn
%                   7:   Uranus
%                   8:   Neptune
%                   9:   Pluto
%     mjd_fb: time of the flyby [Modified Julian Days 2000]
%     VV_minus:  s/c velocity vector before the flyby(heliocentric frame)
%     VV_plus: s/c velocity vector before the flyby(heliocentric frame)
%     dt:  time  step-size 
%     zoomFactor 

% 
% OUTPUT:
%     t_SOI: Total permanence time inside the SOI
%     DeltaV: Delta velocity of the Power Gravity Assist
%     h_GA: Altitude of the closest approach with the planet
%     r_p:  Periapsis radius of the hyperbolas 
%     e_minus: eccentricity of the incoming asymptote
%     e_plus: eccentricity of the outcoming asymptote
%
% CALLED FUNCTIONS:
%     astroConstants.m
%     poweredGA.m
%     plotPlanet.m
%     uplanet.m
%     kep2car.m


mu_s = astroConstants(4); % [km^3/s^2]
mu_p = astroConstants(10 + fb_ID);
R_p = astroConstants(20 + fb_ID);
G = astroConstants(1);
m_p = mu_p/G;
m_s = mu_s/G;

kep_fb = uplanet(mjd_fb, fb_ID);
[RR2, VV_p] = kep2car(kep_fb(1), kep_fb(2), kep_fb(3), kep_fb(4), kep_fb(5), kep_fb(6), mu_s);
[DeltaV, h_GA, r_p, v_p_minus, v_p_plus] = poweredGA(fb_ID, VV_p, VV_minus, VV_plus);

vv_inf_minus = VV_minus - VV_p;
v_inf_minus = norm(vv_inf_minus);

vv_inf_plus = VV_plus - VV_p;
v_inf_plus = norm(vv_inf_plus);
delta_minus  = @(r_p) 2*asin(1./(1+ r_p * v_inf_minus^2/mu_p));
delta_plus   = @(r_p) 2*asin(1./(1+ r_p * v_inf_plus^2/mu_p));


e_minus = 1/sin(delta_minus(r_p)/2);
e_plus = 1/sin(delta_plus(r_p)/2);
a_minus_abs = -r_p/(1-e_minus); % note that for an hyperbola: a<0
a_plus_abs = -r_p/(1-e_plus);



r_s2fb = norm(RR2);
R_SOI = r_s2fb*(m_p/m_s)^(2/5);

% Incoming Hyperbola
a = a_minus_abs;
e = e_minus;
v_inf = v_inf_minus;

Delta = abs(a)*sqrt(e^2-1); % Impact parameter
h = Delta * v_inf;
theta = acos( 1/e * (h^2/(mu_p*R_SOI) - 1) );
F = 2*atanh(sqrt((e-1)/(e+1))*tan(theta/2));

deltaT_1 = sqrt(abs(a)^3/mu_p)*(e*sinh(F)-F);

% Outocming Hyperbola
a = a_plus_abs;
e = e_plus;
v_inf = v_inf_plus;

Delta = abs(a)*sqrt(e^2-1); % Impact parameter
h = Delta * v_inf;
theta = acos( 1/e * (h^2/(mu_p*R_SOI) - 1) );
F = 2*atanh(sqrt((e-1)/(e+1))*tan(theta/2));

deltaT_2 = sqrt(abs(a)^3/mu_p)*(e*sinh(F)-F);


% Total permanence time inside the SOI
t_SOI = deltaT_1 + deltaT_2;





% ***************** Hyperbolic arcs ******************

% Set options
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );


% Perform the integration
tspan = [0:dt:deltaT_1];
rr_p = [r_p 0 0]';
vv_p_minus = [0 v_p_minus 0]';
y0 = [rr_p; vv_p_minus];
[time_vec, StateMat_dep] = ode113( @(t,y) ode_2body(t,y, mu_p), -tspan, y0, options );
X_minus = StateMat_dep(:,1); Y_minus = StateMat_dep(:,2); Z_minus = StateMat_dep(:,3);
X_min_minus = min(X_minus);


tspan = [0:dt:deltaT_2];
rr_p = [r_p 0 0]';
vv_p_plus = [0 v_p_plus 0]';
y0 = [rr_p; vv_p_plus];
[time_vec, StateMat_dep] = ode113( @(t,y) ode_2body(t,y, mu_p), tspan, y0, options );
X_plus = StateMat_dep(:,1); Y_plus = StateMat_dep(:,2); Z_plus = StateMat_dep(:,3);
X_min_plus = min(X_plus);

% ***************** Asymptotes ******************

x_at_y0_minus = r_p+a_minus_abs;
x_at_y0_plus = r_p+a_plus_abs;
Y_asy_minus = @(x)   tan(pi/2 - delta_minus(r_p)/2) * (x-x_at_y0_minus);
Y_asy_plus =  @(x)  -tan(pi/2 - delta_plus(r_p)/2)  * (x-x_at_y0_plus);



% PLOT
h = figure;
view([0,90])
Planet = plotPlanet(fb_ID, [0,0,0], h);
axis equal
grid on
hold on

for i = 1:2
    
    if i == 2
        axes('position',[.675 .175 .3 .3])
        box on % put box around new pair of axes
        axis equal
        grid off
        hold on
        Planet = plotPlanet(fb_ID, [0,0,0], h);
        dt = dt/zoomFactor;
        % Perform the zoomed integration
        deltaT_zoom = max(deltaT_1, deltaT_2);
        tspan = [0:dt:deltaT_zoom/zoomFactor];
        rr_p = [r_p 0 0]';
        vv_p_minus = [0 v_p_minus 0]';
        y0 = [rr_p; vv_p_minus];
        [time_vec, StateMat_dep] = ode113( @(t,y) ode_2body(t,y, mu_p), -tspan, y0, options );
        X_minus = StateMat_dep(:,1); Y_minus = StateMat_dep(:,2); Z_minus = StateMat_dep(:,3);
        X_min_minus = min(X_minus);


        tspan = [0:dt:deltaT_zoom/zoomFactor];
        rr_p = [r_p 0 0]';
        vv_p_plus = [0 v_p_plus 0]';
        y0 = [rr_p; vv_p_plus];
        [time_vec, StateMat_dep] = ode113( @(t,y) ode_2body(t,y, mu_p), tspan, y0, options );
        X_plus = StateMat_dep(:,1); Y_plus = StateMat_dep(:,2); Z_plus = StateMat_dep(:,3);
        X_min_plus = min(X_plus);
        
        if x_at_y0_minus > x_at_y0_plus % Better scaled visualization by cutting the 'longest' asymptote
            x_at_y0_minus = x_at_y0_plus;
        else
            x_at_y0_plus = x_at_y0_minus;
        end
    end
    
    plot3(X_minus,Y_minus,Z_minus, 'LineWidth', 2)
    plot3(X_plus,Y_plus,Z_plus, 'LineWidth', 2)
        

        X = [X_min_minus:x_at_y0_minus];
                if i == 2
                    X = X(1:round(0.5*length(X)));
                end
        plot3(X, Y_asy_minus(X), zeros(length(X),1), '--k', 'LineWidth', .5) 

        X = [X_min_plus:x_at_y0_plus];
                if i == 2
                    X = X(1:round(0.5*length(X)));
                end
        plot3(X, Y_asy_plus(X), zeros(length(X),1), '--k', 'LineWidth', .5, 'HandleVisibility','off')
        % Just label one of the two asymptotes -> HandleVisibility:off for
        % the second one

    if i == 1
        axis equal
        legend({'Incoming Hyperbola', 'Outcoming Hyperbola', 'Asymptotes'}, 'Interpreter', 'latex')
    elseif i == 2
        axis tight
    end


end


end

