% 
% PROTOTYPE:
%   % Gauss_perturbed
% 
% OUTPUT
% this script propagates given initial conditions with gauss propagation
% method, with the possibility to define plotting and evaluation preferences
% in the PlanetaryMission_group_14 script,
%  
% 
% CONTRIBUTORS
% Ciro Salvi
% 
% VERSIONS
% 2020-02-11

if strcmp(J2_pert , 'on')
    J2 = J_2;
    if strcmp(drag_pert , 'on')
        Perturbation = 'global';
    else Perturbation = 'j2';
    end
else J2 = 0;
    if strcmp(drag_pert , 'on')
        Perturbation = 'drag';
    else Perturbation = 'off';
    end
end
%% initial evaluations and integration
kep = [a,e,i,OM,om,f0];

if moln_comparison
kep = kepmol_cut(1,:);
kep(end) = kep(end)/180*pi;
end

options = odeset ( 'RelTol', 1e-13, 'AbsTol', 1e-14 , 'events', @check_earth_radius_gauss );
tspan = [0,number_revolutions*2*pi*sqrt(kep(1)^3/mu)];
[ T_Gauss, kep_Gauss ] = ode113(@(t,y) OdeKepR2BP(t, y, Perturbation , 'rsw',am_ratio,cd,R_e,w_E), tspan, kep, options);
omegadot_Gauss_j2 = 3*J2*R_e^2*sqrt(mu)/2./(1-kep_Gauss(:,2).^2).^2./kep_Gauss(:,1).^(7/2).*(2-2.5*(sin(kep_Gauss(:,3))).^2);   % rad
OMEGAdot_Gauss_j2 = 3*J2*R_e^2*sqrt(mu)/2./(1-kep_Gauss(:,2).^2).^2./kep_Gauss(:,1).^(7/2).*(cos(kep_Gauss(:,3)));              % rad
%% Kep elements graphs
if (check_graphs_kep)
n=2; % number of periods for movmean
figure()
plot( T_Gauss/86400,kep_Gauss(:,1),T_Gauss/86400,filter_mov(kep_Gauss(:,1),T_Gauss, kep_Gauss,mu,n))
grid on
xlabel ('time')
ylabel ('a [km]')
title ('semimajor axis evolution (and filtered)')

figure()
plot( T_Gauss/86400,kep_Gauss(:,2),T_Gauss/86400,filter_mov(kep_Gauss(:,2),T_Gauss, kep_Gauss,mu,n) )
grid on
xlabel ('time')
ylabel ('e')
title ('eccentricity evolution (and filtered)')

figure()
plot( T_Gauss/86400,kep_Gauss(:,3)*180/pi,T_Gauss/86400,filter_mov(kep_Gauss(:,3),T_Gauss, kep_Gauss,mu,n)*180/pi )
grid on
xlabel ('time')
ylabel ('i [deg]')
title ('inclination evolution (and filtered)')

figure()
plot( T_Gauss/86400,kep_Gauss(:,4)*180/pi,T_Gauss/86400,filter_mov(kep_Gauss(:,4),T_Gauss, kep_Gauss,mu,1.5)*180/pi )
grid on
xlabel ('time')
ylabel ('OM [deg]')
title ('RAAN evolution (and filtered)')

figure()
plot( T_Gauss/86400,kep_Gauss(:,5)*180/pi,T_Gauss/86400,filter_mov(kep_Gauss(:,5),T_Gauss, kep_Gauss,mu,1.5) *180/pi)
grid on
xlabel ('time')
ylabel ('om [deg]')
title ('argument of periapsis evolution (and filtered)')

figure()
plot( T_Gauss/86400,kep_Gauss(:,6)*180/pi,T_Gauss/86400,filter_mov(kep_Gauss(:,6),T_Gauss, kep_Gauss,mu,1.5)*180/pi )
grid on
xlabel ('time')
ylabel ('theta [deg]')
title ('true anomaly evolution (and filtered)')

figure()
plot( T_Gauss,omegadot_Gauss_j2*180/pi,T_Gauss,OMEGAdot_Gauss_j2*180/pi )
grid on
xlabel ('time')
ylabel ('theta [deg]')
title ('RAAN and argument of perigee derivative')
legend ('Perigee Precession','Nodal Regression')
end

%% Cartesian coordinates evaluation
if check_draw || check_graphs || check_energy
    v_Gauss = [];
    r_Gauss = [];
    for ii = 1:(length(kep_Gauss(:,6)))
        [position, velocity] = kep2car(kep_Gauss(ii,1) , kep_Gauss(ii,2) , kep_Gauss(ii,3), kep_Gauss(ii,4), kep_Gauss(ii,5), kep_Gauss(ii,6), mu,'rad');
        r_Gauss = [r_Gauss;position'];
        v_Gauss = [v_Gauss;velocity'];
    end
end

%% draw the orbit
if check_draw
Earth_plot
period_nume = floor(T_Gauss./T_sat)+1;
no_period = period_nume(end);
scatter3(r_Gauss(:,1),r_Gauss(:,2),r_Gauss(:,3),2,period_nume*256/no_period )
% plot3 (r_Gauss(1:1000,1),r_Gauss(1:1000,2),r_Gauss(1:1000,3));
xlabel('r_x');
ylabel('r_y');
zlabel('r_z');
grid on
end

%% energy / eccentricity / specific angular momentum
if check_energy
    r_mag_Gauss = vecnorm(r_Gauss,2,2);

if drag_pert == 'on'
    Drag_dissipation_Gauss = [];
for n=(1:length(r_mag_Gauss)-1)
v_rel = [v_Gauss(n,1)+(w_E*r_Gauss(n,2)); v_Gauss(n,2)-(w_E*r_Gauss(n,1)); v_Gauss(n,3)]; 
Drag_dissipation_Gauss = [Drag_dissipation_Gauss -0.5*am_ratio*cd*rho(norm(r_Gauss(n,:))-R_e)*norm(v_rel)^2*1000*norm(r_Gauss(n+1,:)-r_Gauss(n,:))];
end
end 

for n=(1:length(r_mag_Gauss))
% eps_Gauss(n) = ((norm(v_Gauss(n,:))^2)/2)-(mu/norm(r_Gauss(n,:)));
RZH = J2*R_e^2/2/(norm(r_Gauss(n,:)))^2*(3*(r_Gauss(n,3)/(norm(r_Gauss(n,:))))^2-1);
eps_Gauss(n) = ((norm(v_Gauss(n,:))^2)/2)-(mu/norm(r_Gauss(n,:)))*(1-RZH);
end

h_Gauss = cross(r_Gauss,v_Gauss);
h_mag_Gauss =vecnorm(h_Gauss,2,2);

evect_Gauss = zeros(length(r_mag_Gauss),3);
% e_mag = zeros(10000,1);
for n=(1:length(r_mag_Gauss))
evect_Gauss(n,:) = (cross(v_Gauss(n,:),h_Gauss(n,:))./mu)-(r_Gauss(n,:)./norm(r_Gauss(n,:)));
end

% e_dot_h = zeros(10000,1);
for n=(1:length(r_mag_Gauss))
e_dot_h_Gauss(n,:) = dot(evect_Gauss(n,:),h_Gauss(n,:));
end
MAX_e_dot_h_Gauss = max(abs(e_dot_h_Gauss));
end    
   
%% plot energy / eccentricity / specific angular momentum
if check_graphs

figure()
plot(T_Gauss/86400,h_Gauss);
hold on
plot(T_Gauss/86400,h_mag_Gauss);
hold off
xlabel('time [days]');
ylabel('h_x, h_y, h_z, ||h||');
title('Specific Angular Momentum vs. Time');
legend ('h_x','h_y','h_z','||h||');

figure()
plot(T_Gauss/86400,evect_Gauss);
hold on 
plot(T_Gauss/86400,kep_Gauss(:,2));
hold off
xlabel('time [days]');
ylabel('Eccentricity (e_x, e_y, e_z, ||e||)');
title('Eccentricity vs. Time');
legend ('e_x','e_y','e_z','||e||');

figure();
plot(T_Gauss/86400,eps_Gauss)%,T_Gauss(1:end-1),Drag_dissipation_Gauss);
xlabel('time [days]');
ylabel('Specific Energy ({\epsilon} [km^2/s^2])');
title('Specific Energy vs. Time');
legend ('{\epsilon}');

figure();
plot(T_Gauss/86400,v_Gauss(:,1));
hold on
plot(T_Gauss/86400,v_Gauss(:,2));
xlabel('time [days]');
ylabel('v_x, v_y [km/s]');
title('Radial and Transversal Velocity vs. Time');
legend ('v_x','v_y');

figure();
plot(T_Gauss/86400,abs(e_dot_h_Gauss));
xlabel('time [days]');
ylabel('\bf |e.h| \rm [km^2/s^2])');
title('e.h dot product vs. Time');
legend ('\bf e.h');

% v_rel4check = [v_Gauss(:,1)+(omg_E*r_Gauss(:,2)) , v_Gauss(:,2)-(omg_E*r_Gauss(:,1)) , v_Gauss(:,3)]; 
% v_rel4check = [v_Cartesian(:,1)+(omg_E*r_Cartesian(:,2)) , v_Cartesian(:,2)-(omg_E*r_Cartesian(:,1)) , v_Cartesian(:,3)]; 
% 
% for nn = 1 : size (r_Gauss,1)-1
%     eps(nn) = ((norm(v_Gauss(nn,:))^2)/2)-(mu/norm(r_Gauss(nn,:)));
%     eps_drag(nn) = -0.5*cd*am_ratio*rho(norm(r_Gauss(nn,:))-R_e)*norm(v_rel4check(nn,:))^2*norm((r_Gauss(nn+1)-r_Gauss(nn)))*1000;
% end
% plot (T_Gauss(1:end-1),eps_drag,T_Gauss, -mu/2/a*ones(length(T_Gauss),1))


end

clear alpha cdata erad erot GMST0 hp ii n npanels position prad tspan check_draw check_graphs