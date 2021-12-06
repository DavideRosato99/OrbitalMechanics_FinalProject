clc;
clear;

mu = astroConstants(13);      % Gravitational parameter [km^3/s^2];


a_1 = 12500; e_1 = 0; i_1 = 0; OM_1 = 0; om_1 = 0;  th_1 = deg2rad(120);
a_2 = 9500; e_2 = 0.3; i_2 = 0; OM_2 = 0; om_2 = 0;  th_2 = deg2rad(250);
ToF = 3300; % Time of flight [s];

[RR1, VV1] = kep2car(a_1, e_1, i_1, OM_1, om_1, th_1, mu);
[RR2, VV2] = kep2car(a_2, e_2, i_2, OM_2, om_2, th_2, mu);
orbitType=0;
Nrev=0;
Ncase=0;
optionsLMR=2;


[a,p,e,ERROR,VVT1,VVT2,TPAR,theta] = lambertMR( RR1, RR2 , ToF, mu, orbitType,Nrev,Ncase, optionsLMR ); 
RR1 = RR1(:); RR2 = RR2(:); VVT1 = VVT1(:); VVT2 = VVT2(:);

y0 = [RR1; VVT1];

% Set options
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

% Set time span
tspan = [0:1:ToF];

% Perform the integration
[T, StateMat] = ode113( @(t,y) ode_2body(t,y, mu), tspan, y0, options );

X = StateMat(:,1); Y = StateMat(:,2); Z = StateMat(:,3);
VX = StateMat(:,4); VY = StateMat(:,5); VZ = StateMat(:,6);

figure
plot3(X,Y,Z)
grid on
hold on
earth_sphere
axis equal
DeltaV_T1 = norm(VVT1-VV1);
DeltaV_T2 = norm(VV2-VVT2);
DeltaV_tot = DeltaV_T1 + DeltaV_T2











