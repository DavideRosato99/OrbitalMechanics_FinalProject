%% EX.1.a

close all
clear all
addpath('textures')
% Data: 
% (vectors in the ecliptic frame; assume circular Earth orbit around the
% Sun)

mu_p = astroConstants(13); % [km^3/s^2]
mu_s = astroConstants(4); % [km^3/s^2]
AU = astroConstants(2); % [km]

vv_inf_minus = [15.1 0 0]'; %km/s
v_inf = norm(vv_inf_minus);
Delta = 9200; % [km]
rr_s2p = [AU; 0; 0]
V_p = sqrt(mu_s/AU); % we assume a circular Earth orbit
VV_p = V_p * [0 1 0]'


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Equations for the hyperbola
a = - mu_p/v_inf^2;
delta= 2*atan(-a/Delta);
e = 1/sin(delta/2);
r_p = a * (1-e);
Delta_V = 2*v_inf*sin(delta/2)




VV_minus = VV_p + vv_inf_minus


% Either assign assist_type as a:
% - word
% - angle in [rad] that identifies the location of the incoming asymptote
assist_types = {'leading', 'trailing', 'under', 'over'};
assist_type = assist_types{3};
assist_type = pi/2



VV_minus
norm2helioPlane1 = cross(rr_s2p, VV_minus);
norm2helioPlane1 = norm2helioPlane1/norm(norm2helioPlane1)

    switch assist_type
        case 'leading'
            rot_dir = -norm2helioPlane1;
        case 'trailing'
            rot_dir = +norm2helioPlane1;
        case 'under'
            rot_dir = - VV_p;
        case 'over'
            rot_dir = + VV_p;      
        otherwise
            theta = assist_type; % counterclockwise angle that defines the location
                                 % of incoming asymptote, considering a rotation around
                                 % the position vector and s.t. theta=0 -> leading-side
            rot_dir = rotVecAroundVecByAngle(norm2helioPlane1, rr_s2p, -theta);  
    end
rot_dir = rot_dir/norm(rot_dir)

vv_inf_plus = rotVecAroundVecByAngle(vv_inf_minus, rot_dir, delta)
VV_plus = VV_p + vv_inf_plus



norm2helioPlane2 = cross(rr_s2p, VV_plus);
norm2helioPlane2 = norm2helioPlane2/norm(norm2helioPlane2)

h = figure;
hold on
grid on
axis equal
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
dt = 10000;
tspan = [-200*86400:dt:0];
y0 = [rr_s2p; -VV_minus];
[time_vec, StateMat] = ode113( @(t,y) ode_2body(t,y, mu_s), tspan, y0, options );
X_arr = StateMat(:,1); Y_arr = StateMat(:,2); Z_arr = StateMat(:,3);
%VX = StateMat_arr(:,4); VY = StateMat_arr(:,5); VZ = StateMat_arr(:,6);
plot3(X_arr,Y_arr,Z_arr, 'LineWidth', 1.5)

tspan = [0:dt:200*86400];
y0 = [rr_s2p; VV_plus];
[time_vec, StateMat] = ode113( @(t,y) ode_2body(t,y, mu_s), tspan, y0, options );
X_arr = StateMat(:,1); Y_arr = StateMat(:,2); Z_arr = StateMat(:,3);
%VX = StateMat_arr(:,4); VY = StateMat_arr(:,5); VZ = StateMat_arr(:,6);
plot3(X_arr,Y_arr,Z_arr, 'LineWidth', 1.5)


view([54, 32])
sunGlobe = plotPlanet(10, [0,0,0], h, 30);








%% EX.1.b
vv_inf_minus = [15.1 0 0]'; % km/s
V_p = sqrt(mu_s/AU); % we assume a circualr Earth orbit
VV_p = V_p * [0 1 0]';
Delta = 9200;
inc_asymptote_angularPos = deg2rad(50)
h = figure;
view([-140, 20])
for inc_asymptote_angularPos = deg2rad([0:10:360])
    % for Delta = [8000:1000:15000]
        [VV_minus, VV_plus] = flybyPlotHelio(vv_inf_minus, VV_p, Delta, rr_s2p , mu_s, mu_p, inc_asymptote_angularPos)
        drawnow
    % end
end




%%  EX.2 - POWERED Flyby
clear all
mu_p = astroConstants(13); % [km^3/s^2]
mu_s = astroConstants(4); % [km^3/s^2]
AU = astroConstants(2); % [km]
R_planet = 6378; % max radius (i.e. equatorial)
h_atm = 200;

VV_minus = [31.5, 4.69, 0]';
VV_plus = [38.58, 0, 0]';

% rr_s2p = [0, -1, 0]';
V_p = sqrt(mu_s/AU); % we assume a circualr Earth orbit
VV_p = V_p * [1 0 0]';


vv_inf_minus = VV_minus - VV_p
v_inf_minus = norm(vv_inf_minus)

vv_inf_plus = VV_plus - VV_p
v_inf_plus = norm(vv_inf_plus)

delta = acos(dot(vv_inf_minus, vv_inf_plus)/(v_inf_minus*v_inf_plus))



delta_minus  = @(r_p) 2*asin(1./(1+ r_p * v_inf_minus^2/mu_p));
delta_plus   = @(r_p) 2*asin(1./(1+ r_p * v_inf_plus^2/mu_p));
r_p_SolveFun = @(r_p) (delta_minus(r_p) + delta_plus(r_p))/2 - delta;
r_p_min = R_planet + h_atm;
options = optimoptions('fsolve','TolFun',1e-14,'Display','off');

r_p = fsolve(r_p_SolveFun, r_p_min, options)  % Periapsis radiu of the hyperbolas [km]



v_p_minus = sqrt(v_inf_minus^2 + 2*mu_p/r_p)
v_p_plus = sqrt(v_inf_plus^2 + 2*mu_p/r_p)

delta_V_poweredFB = abs(v_p_plus - v_p_minus)

h_ga = r_p - R_planet


% ****************************************************
% ***************** Hyperbolic arcs ******************
% ****************************************************

% Set options
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

dt = 10; % Step size [s]

% Set time span
tspan = [0:dt:5000];

rr_p = [r_p 0 0]';
vv_p_minus = [0 v_p_minus 0]';
y0 = [rr_p; vv_p_minus];

% Perform the integration
[time_vec, StateMat_dep] = ode113( @(t,y) ode_2body(t,y, mu_p), -tspan, y0, options );
X = StateMat_dep(:,1); Y = StateMat_dep(:,2); Z = StateMat_dep(:,3);

h = figure;
view([0,90])
earthGlobe = plotPlanet(3, [0,0,0], h, (astroConstants(23)/astroConstants(3)));
axis equal
grid on
hold on
plot3(X,Y,Z, 'LineWidth', 2)
X_min_minus = min(X);


rr_p = [r_p 0 0]';
vv_p_plus = [0 v_p_plus 0]';
y0 = [rr_p; vv_p_plus];
[time_vec, StateMat_dep] = ode113( @(t,y) ode_2body(t,y, mu_p), tspan, y0, options );
X = StateMat_dep(:,1); Y = StateMat_dep(:,2); Z = StateMat_dep(:,3);
plot3(X,Y,Z, 'LineWidth', 2)
X_min_plus = min(X);





% ****************************************************
% ***************** Asymptotes ******************
% ****************************************************

e_minus = 1/sin(delta_minus(r_p)/2)
e_plus = 1/sin(delta_plus(r_p)/2)
a_minus_abs = -r_p/(1-e_minus); % note that for an hyperbola: a<0
a_plus_abs = -r_p/(1-e_plus);


x_at_y0_minus = r_p+a_minus_abs
x_at_y0_plus = r_p+a_plus_abs
Y_minus = @(x)   tan(pi/2 - delta_minus(r_p)/2) * (x-x_at_y0_minus);
Y_plus =  @(x)  -tan(pi/2 - delta_plus(r_p)/2)  * (x-x_at_y0_plus);

X = [X_min_minus:x_at_y0_minus];
plot3(X, Y_minus(X), zeros(length(X),1), '--k', 'LineWidth', .5)

X = [X_min_plus:x_at_y0_plus];
plot3(X, Y_plus(X), zeros(length(X),1), '--k', 'LineWidth', .7, 'HandleVisibility','off')

legend('Incoming Hyperbola', 'Outcoming Hyperbola', 'Asymptotes')





