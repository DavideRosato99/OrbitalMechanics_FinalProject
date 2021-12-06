r1=[-21800; 37900; 0];
r2=[27300; 27700; 0];
TOF=15*3600+6*60+40;
muE=astroConstants(13);
orbitType=0;
Nrev=0;
Ncase=0;
optionsLMR=2;
[A,P,E,ERROR,v1,v2,TPAR,THETA] = lambertMR(r1,r2,TOF,muE,orbitType,Nrev,Ncase,optionsLMR)
 %[r1,v1] = plot the orbit with your code of module 1

r1 = r1(:); r2 = r2(:); v1 = v1(:); v2 = v2(:);



y0 = [r1; v1];

% Set options
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

% Set time span
tspan = [0:100:TOF];

% Perform the integration
[T, StateMat] = ode113( @(t,y) ode_2body(t,y, muE), tspan, y0, options )

X = StateMat(:,1); Y = StateMat(:,2); Z = StateMat(:,3);
VX = StateMat(:,4); VY = StateMat(:,5); VZ = StateMat(:,6);



figure
plot3(X,Y,Z)
grid on
hold on
earth_sphere
axis equal
