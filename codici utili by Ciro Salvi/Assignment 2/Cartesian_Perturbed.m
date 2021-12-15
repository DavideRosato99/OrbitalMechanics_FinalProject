% 
% PROTOTYPE:
%   % Cartesian_perturbed
% 
% OUTPUT
% this script propagates given initial conditions with Cartesian propagation
% method, with the possibility to define plotting and evaluation preferences
% in the PlanetaryMission_group_14 script,
%  
% 
% CONTRIBUTORS
% Ciro Salvi
% Suhailah Alkhawashke
% 
% VERSIONS
% 2020-02-11% Cartesian_Perturbed

if moln_comparison
kep = kepmol_cut(1,:);
[r0,v0] = kep2car(kep(1) , kep(2) , kep(3), kep(4), kep(5), kep(6), mu, 'rad');
else
[r0,v0] = kep2car(a,e,i,OM,om,f0,mu,'rad');
end
x0 = [r0 ; v0];  
% r0_mag = norm(r0);

%% Initial evaluations and integration
if strcmp(J2_pert , 'on')
    J2 = J_2;
else J2 = 0;
end

omg_E = w_E*pi/180/3600;
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 , 'events', @check_earth_radius_cartesian);
tspan = [0,T_sat*n_revolutions];

[ T_Cartesian, x_Cartesian ] = ode113( @(t,x) Perturbed2BP_fnc(t,x,mu,J_2,R_e,omg_E,cd,am_ratio,J2_pert,drag_pert), tspan, x0, options);
r_Cartesian = x_Cartesian(:,1:3);
v_Cartesian = x_Cartesian(:,4:6);


%% Keplerian elements evaluation
if check_kep
kep_Cartesian = [];
for kk = 1 : size(x_Cartesian,1)
    [aprov,eprov,iprov,OMprov,omprov,thprov] = car2kep(r_Cartesian(kk,:)',v_Cartesian(kk,:)',mu);
    Kepprov = [aprov,eprov,iprov,OMprov,omprov,thprov];
    kep_Cartesian = [kep_Cartesian;Kepprov];
end
kep_Cartesian(:,6) = unwrap(kep_Cartesian(:,6));
kep_Cartesian(:,5) = unwrap(kep_Cartesian(:,5));
kep_Cartesian(:,4) = unwrap(kep_Cartesian(:,4));
omegadot_Cartesian_j2 = 3*J2*R_e^2*sqrt(mu)/2./(1-kep_Cartesian(:,2).^2).^2./kep_Cartesian(:,1).^(7/2).*(2-2.5*(sin(kep_Cartesian(:,3))).^2);   % rad
OMEGAdot_Cartesian_j2 = 3*J2*R_e^2*sqrt(mu)/2./(1-kep_Cartesian(:,2).^2).^2./kep_Cartesian(:,1).^(7/2).*(cos(kep_Cartesian(:,3)));              % rad
end

%% Kep elements graphs
if check_graphs_kep
n=2; % number of periods for movmean
figure()
plot( T_Cartesian/86400 ,kep_Cartesian(:,1),T_Cartesian/86400,filter_mov(kep_Cartesian(:,1),T_Cartesian, kep_Cartesian,mu,n))
grid on
xlabel ('time')
ylabel ('a [km]')
title ('semimajor axis evolution (and filtered)')

figure()
plot( T_Cartesian/86400,kep_Cartesian(:,2),T_Cartesian/86400,filter_mov(kep_Cartesian(:,2),T_Cartesian, kep_Cartesian,mu,n) )
grid on
xlabel ('time')
ylabel ('e')
title ('eccentricity evolution (and filtered)')

figure()
plot( T_Cartesian/86400,kep_Cartesian(:,3)*180/pi,T_Cartesian/86400,filter_mov(kep_Cartesian(:,3),T_Cartesian, kep_Cartesian,mu,n)*180/pi )
grid on
xlabel ('time')
ylabel ('i [deg]')
title ('inclination evolution (and filtered)')

figure()
plot( T_Cartesian/86400,kep_Cartesian(:,4)*180/pi,T_Cartesian/86400,filter_mov(kep_Cartesian(:,4),T_Cartesian, kep_Cartesian,mu,n)*180/pi )
grid on
xlabel ('time')
ylabel ('OM [deg]')
title ('RAAN evolution (and filtered)')

figure()
plot( T_Cartesian/86400,kep_Cartesian(:,5)*180/pi,T_Cartesian/86400,filter_mov(kep_Cartesian(:,5),T_Cartesian, kep_Cartesian,mu,n)*180/pi )
grid on
xlabel ('time')
ylabel ('om [deg]')
title ('argument of periapsis evolution (and filtered)')

figure()
plot( T_Cartesian/86400,kep_Cartesian(:,6)*180/pi,T_Cartesian/86400,filter_mov(kep_Cartesian(:,6),T_Cartesian, kep_Cartesian,mu,n)*180/pi )
grid on
xlabel ('time')
ylabel ('theta [deg]')
title ('true anomaly evolution (and filtered)')

figure()
plot( T_Cartesian,omegadot_Cartesian_j2*180/pi,T_Cartesian,OMEGAdot_Cartesian_j2*180/pi )
grid on
xlabel ('time')
ylabel ('theta [deg]')
title ('RAAN and argument of perigee derivative')
legend ('Perigee Precession','Nodal Regression')
end


%% draw the orbit
if check_draw
Earth_plot
axis equal
hold on
% plot3(r_Cartesian(:,1),r_Cartesian(:,2),r_Cartesian(:,3),'-' )
period_nume = floor(T_Cartesian./T_sat)+1;
no_period = period_nume(end);
scatter3(r_Cartesian(:,1),r_Cartesian(:,2),r_Cartesian(:,3),2,period_nume*256/no_period )
grid on 
xlabel('r_x');
ylabel('r_y');
zlabel('r_z');
hold off
end



%% energy / eccentricity / specific angular momentum
if check_energy
    
    r_mag_Cartesian = vecnorm(r_Cartesian,2,2);

% eps = zeros(length(r_mag),1);
for n=(1:length(r_mag_Cartesian))
RZH = J2*R_e^2/2/(norm(r_Cartesian(n,:)))^2*(3*(r_Cartesian(n,3)/(norm(r_Cartesian(n,:))))^2-1);
eps_Cartesian(n) = ((norm(v_Cartesian(n,:))^2)/2)-(mu/norm(r_Cartesian(n,:)))*(1-RZH);
end

h_Cartesian = cross(r_Cartesian,v_Cartesian);
h_mag_Cartesian =vecnorm(h_Cartesian,2,2);

evect_Cartesian = zeros(length(r_mag_Cartesian),3);
% e_mag = zeros(10000,1);
for n=(1:length(r_mag_Cartesian))
evect_Cartesian(n,:) = (cross(v_Cartesian(n,:),h_Cartesian(n,:))./mu)-(r_Cartesian(n,:)./norm(r_Cartesian(n,:)));
e_mag__Cartesian(n,:) = norm(evect_Cartesian(n,:));
end

% e_dot_h = zeros(10000,1);
for n=(1:length(r_mag_Cartesian))
e_dot_h_Cartesian(n,:) = dot(evect_Cartesian(n,:),h_Cartesian(n,:));
end
MAX_e_dot_h_Cartesian = max(abs(e_dot_h_Cartesian));
end    

%% plot energy / eccentricity / specific angular momentum
if check_graphs
figure()
plot(T_Cartesian/86400,h_Cartesian);
hold on
plot(T_Cartesian/86400,h_mag_Cartesian);
hold off
xlabel('time [days]');
ylabel('h_x, h_y, h_z, ||h||');
title('Specific Angular Momentum vs. Time');
legend ('h_x','h_y','h_z','||h||');

figure()
plot(T_Cartesian/86400,evect_Cartesian);
hold on 
plot(T_Cartesian/86400,e_mag__Cartesian);
hold off
xlabel('time [days]');
ylabel('Eccentricity (e_x, e_y, e_z, ||e||)');
title('Eccentricity vs. Time');
legend ('e_x','e_y','e_z','||e||');

figure();
plot(T_Cartesian/86400,eps_Cartesian);
xlabel('time [days]');
ylabel('Specific Energy ({\epsilon} [km^2/s^2])');
title('Specific Energy vs. Time');
legend ('{\epsilon}');

figure();
plot(T_Cartesian/86400,v_Cartesian(:,1));
hold on
plot(T_Cartesian/86400,v_Cartesian(:,2));
xlabel('time [days]');
ylabel('v_x, v_y [km/s]');
title('Radial and Transversal Velocity vs. Time');
legend ('v_x','v_y');

figure();
plot(T_Cartesian/86400,abs(e_dot_h_Cartesian));
xlabel('time [days]');
ylabel('\bf |e.h| \rm [km^2/s^2])');
title('e.h dot product vs. Time');
legend ('\bf e.h');
end

clear alpha cdata check_draw check_graphs drag_pert erad erot GMST0 J2_pert n_revolutions v0 x0 tspan  npanels omg_E prad