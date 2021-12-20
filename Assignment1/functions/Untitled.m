clear all
close all
clc

a = 3.2776e4;                 % [km] Orbit semi-major axis
e = 0.1;                   % [-] Orbit eccentricity
i = 134.2783;                 % [deg] Orbit inclination

orb = [a e deg2rad(i) 0 0 0];
muE = astroConstants(13);
[rr0, vv0] = par2car(orb, muE);

Torbit = 2*pi * sqrt(a^3/muE);

Y0 = [rr0; vv0];


tspan = linspace(0, Torbit, 500);

options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14);
[T, Y] = ode113(@ode_2bp, tspan, Y0, options, muE, 'cart');


%%
figure
% Trajectory
plot3(Y(:,1), Y(:,2), Y(:,3)); hold on

% Satellite
sat = plot3(Y(1,1), Y(1,2), Y(1,3), 'ro', 'MarkerSize', 2);
zoom(1.1); 
omega = rad2deg(acos(Y(1,2) / norm(Y(1,:))));
view(omega, 10)

% Earth
[Xe, Ye, Ze] = sphere(1000);
Xe = 3671*Xe; Ye = 3671*Ye; Ze = 3671*Ze;
earth = imread('earth.png');  
surf(Xe, Ye, Ze, 'CData', flipud(earth), 'FaceColor', 'texture', 'edgecolor', 'none');

% Galaxy
XLIM = xlim;
YLIM = ylim;
ZLIM = zlim;

zc2 = ((YLIM(2)^2 * ZLIM(2)^2) - (YLIM(1)^2 * ZLIM(2)^2))/(YLIM(2)^2 - YLIM(1)^2);
yc2 = (YLIM(1)^2 * zc2)/(zc2 - ZLIM(2)^2);

[Xm, Ym, Zm] = ellipsoid(Y(1,1), Y(1,2), Y(1,3), sqrt(zc2), sqrt(yc2/10000), sqrt(zc2));
xMAX = min(abs(xlim));
yMAX = min(abs(ylim));
zMAX = min(abs(zlim));
maxR = min([xMAX yMAX zMAX]);
milkyWay = imread('milkyWay.jpg');  
% M = surf(Xm, Ym, Zm, 'CData', milkyWay, 'FaceColor', 'texture', 'edgecolor', 'none');
M = surf(Xm, Ym, Zm);

axis equal
error('ciao')
%%
for i = 1:10
    delete(sat)
    delete(M);
    sat = plot3(Y(i,1), Y(i,2), Y(i,3), 'ro', 'MarkerSize', 2); hold on
    
    zoom(1.1); 
    omega = rad2deg(acos(Y(i,2) / norm(Y(i,:))));
    view(omega, 10)
    
    xMAX = min(abs(xlim));
    yMAX = min(abs(ylim));
    zMAX = min(abs(zlim));
    
    maxR = min([xMAX, yMAX, zMAX]);
    [Xm, Ym, Zm] = sphere(1000);
    Xm = maxR*Xm; Ym = maxR*Ym; Zm = maxR*Zm;
    
    milkyWay = imread('milkyWay.jpg');  
    M = surf(Xm, Ym, Zm, 'CData', milkyWay, 'FaceColor', 'texture', 'edgecolor', 'none');
    grid on;
   
    
    
    drawnow limitrate
end