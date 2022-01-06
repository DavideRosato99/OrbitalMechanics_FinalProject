function [] = plotMission(data)
% plotMission - Function to plot the trajectory followed by the satellite
%
% PROTOTYPE
%   [] = plotMission(data)
%
% INPUT:
%   data     
%
% OUTPUT: [-]
%
% CALLED FUNCTIONS: 
%   astroConstants
%   date2mjd2000
%   ephMoon
%   J2Pert
%   MoonPer
%   par2car
%
% CONTRIBUTORS:
%   Rosato Davide               10618468
%   Saba Mohammadi Yengeje      10789462
%   Spinelli Jason              10618465
%   Tagliati Alessia            10635119
%
% VERSIONS
%   2021-10-21: Release
%
% -------------------------------------------------------------------------

% depDate      = data.starting.depDate;
depDate      = [2041, 8, 1, 3, 22, 17];
depDate      = date2mjd2000(depDate);
FlyByDate    = [2042, 6, 9, 10, 40, 51];% data.starting.depDate;
FlyByDate    = date2mjd2000(FlyByDate);
% arrDate      = data.starting.arrDate;
arrDate      = [2042, 12, 6, 6, 59, 50];
arrDate      = date2mjd2000(arrDate);

minHfl = data.optimization.minHfl;

depPlanID    = data.starting.depPlanID;   
depColor     = [1 85/255 0];
flyByPlanID  = data.starting.flyByPlanID; 
fbColor     = [85/255 85/255 1];
arrPlanID    = data.starting.arrPlanID;   
arrColor     = [47/255 85/255 186/255];

muE = data.constants.muE;
muS = data.constants.muS;
AU  = data.constants.AU;
Re  = data.constants.Re;  

%% determine position of the planet
[depkep,~] = uplanet(depDate, depPlanID);
[flybykep,~] = uplanet(FlyByDate, flyByPlanID);
[arrkep,~] = uplanet(arrDate, arrPlanID);

Td = 2*pi*sqrt((depkep(1)^3)/muS);
Tfb = 2*pi*sqrt((flybykep(1)^3)/muS);
Ta = 2*pi*sqrt((arrkep(1)^3)/muS);
tspanD = linspace(0,Td,300);
tspanFB = linspace(0,Tfb,300);
tspanA = linspace(0,Ta,300);

for i = 1:300
   [depkep,~] = uplanet(depDate + tspanD(i)/(24*3600), depPlanID);
   [rrDin(:,i), vvDin(:,i)] = par2car(depkep, muS);
   [flybykep,~] = uplanet(FlyByDate + tspanFB(i)/(24*3600), flyByPlanID);
   [rrDfb(:,i), vvDfb(:,i)] = par2car(flybykep, muS);
   [arrkep,~] = uplanet(arrDate + tspanA(i)/(24*3600), arrPlanID);
   [rrDarr(:,i), vvDarr(:,i)] = par2car(arrkep, muS);
end

%% orbit planet
% % figure
% aminor = 13000000;                 % Torus minor radius
% % Rmajor = mean(vecnorm(rrDin));   % Torus major radius
% theta  = linspace(-pi, pi, 64)   ; % Poloidal angle
% phi    = linspace(0., 2.*pi, 64) ; % Toroidal angle
% [t, p] = meshgrid(phi, theta);
% x = (Rmajor + aminor.*cos(p)) .* cos(t);
% y = (Rmajor + aminor.*cos(p)) .* sin(t);
% z = aminor.*sin(p);
% H1 = surf(x, y, z);
% axis equal
% % rotate(H1,[1 0 0],0,[0 0 0]);
% set(H1,'edgecolor','none','facecolor',depColor,'FaceAlpha',0.2)

%% lambert's arc

[depkep,~] = uplanet(depDate, depPlanID);
[rr1, vv1] = par2car(depkep, muS);
[flybykep,~] = uplanet(FlyByDate, flyByPlanID);
[rr2, vv2] = par2car(flybykep, muS);
[arrkep,~] = uplanet(arrDate, arrPlanID);
[rr3, vv3] = par2car(arrkep, muS);



TOF1 = (FlyByDate -depDate)*24*3600;
TOF2 = (arrDate -FlyByDate)*24*3600;


[a1, p1, e1, error1, vi1, vf1, Tpar1, theta1, a2, p2, e2, ...
    error2, vi2, vf2, Tpar2, theta2, DV1, DV2, errorFB, r_p, h_ga, delta, ...
    delta_V_powFB, e_minus, e_plus, a_minus, a_plus] = deltaVtot(...
    rr1, rr2, rr3, vv1, vv2, vv3, TOF1, TOF2, minHfl, flyByPlanID);

%% first leg
firstlegin = car2par(rr1,vi1',muS);

firstlegend= car2par(rr2,vf1',muS);
if firstlegin(6) < firstlegend(6)
    thspan = linspace(firstlegin(6),firstlegend(6),300);
else
    thspan = linspace(firstlegin(6),firstlegend(6)+2*pi,300);
end


orb = firstlegin;
for i = 1:length(thspan)
   orb(6) = thspan(i);
   [rr_firstleg(:,i),~] = par2car(orb,muS); 
     
end

plot3(rr_firstleg(1,:),rr_firstleg(2,:),rr_firstleg(3,:)); hold on



%% second leg

secondlegin = car2par(rr2,vi2',muS);

secondlegend= car2par(rr3,vf2',muS);
if secondlegin(6) <secondlegend(6)
    thspan = linspace(secondlegin(6),secondlegend(6),300);
else
    thspan = linspace(secondlegin(6),secondlegend(6)+2*pi,300);
end


orb = secondlegin;
for i = 1:length(thspan)
   orb(6) = thspan(i);
   [rr_secondleg(:,i),~] = par2car(orb,muS); 
     
end

plot3(rr_secondleg(1,:),rr_secondleg(2,:),rr_secondleg(3,:));

%%
plot3(rrDin(1,:),rrDin(2,:),rrDin(3,:),'--k');
axis equal; grid on; hold on;
plot3(rrDfb(1,:),rrDfb(2,:),rrDfb(3,:),'--k');
plot3(rrDarr(1,:),rrDarr(2,:),rrDarr(3,:),'--k');
legend('First Leg', 'Second Leg')

%%
% load('planet_texture');
cdataS = imread('Sun_texture.jpg');
cdataM = imread('Mars_texture.jpg');
cdataE = imread('Earth_texture.jpg');
cdataV = imread('Venus_texture.jpg');

[Xs,Ys,Zs] = sphere(100);
Xs = astroConstants(3)*Xs*15;
Ys = astroConstants(3)*Ys*15;
Zs = astroConstants(3)*Zs*15;
S = surf(Xs,Ys,Zs,'FaceColor', 'none', 'EdgeColor', 'none');
set(S, 'FaceColor', 'texturemap', 'CData', cdataS);

r1 = astroConstants(20+depPlanID);
r2 = astroConstants(20+flyByPlanID);
r3 = astroConstants(20+arrPlanID);


[X1,Y1,Z1] = sphere(100);
X1 = r1*X1*1000 ;
Y1 = r1*Y1*1000 ;
Z1 = r1*Z1*1000 ;
depp = surf(X1 + rr1(1),Y1 + rr1(2),Z1 + rr1(3),'FaceColor', 'none', 'EdgeColor', 'none');
set(depp, 'FaceColor', 'texturemap', 'CData', cdataM);

[X2,Y2,Z2] = sphere(100);
X2 = r2*X2*1100 ;
Y2 = r2*Y2*1100 ;
Z2 = r2*Z2*1100 ;
fbp = surf(X2 + rr2(1),Y2 + rr2(2),Z2 + rr2(3),'FaceColor', 'none', 'EdgeColor', 'none');
set(fbp, 'FaceColor', 'texturemap', 'CData', cdataE);

[X3,Y3,Z3] = sphere(100);
X3 = r3*X3*900 ;
Y3 = r3*Y3*900 ;
Z3 = r3*Z3*900 ;
arrp = surf(X3 + rr3(1),Y3 + rr3(2),Z3 + rr3(3),'FaceColor', 'none', 'EdgeColor', 'none');
set(arrp, 'FaceColor', 'texturemap', 'CData', cdataV);

%% flyby 
vInfMin = vf1;
vInfPlus = vi2;
vMin = norm(vInfMin); 
vPlus = norm(vInfPlus);

delta = acos(dot(vInfMin,vInfPlus)/(vMin*vPlus));

% extract attractor mass parameter and radius
mu = astroConstants(10+flyByPlanID);
rP = astroConstants(20+flyByPlanID);

% valuate epericentre height
rp = rP + minHfl;

% evaluate both branches eccentricities, semimajor axes
eMin = 1 + rp*vMin^2/mu;
ePlus = 1 + rp*vPlus^2/mu; 
aMin = -mu/vMin^2; 
aPlus = -mu/vPlus^2;

% evaluate hyperbolic plane directions
planeDir = cross(vInfMin, vInfPlus)/(vMin*vPlus*sin(delta));
% norm(planeDir)
% dot(planeDir, vInfMin)
% dot(planeDir, vInfPlus)

% evaluate plane parameters of the hyperbolic trajectory
i = acos(planeDir(3));
N = cross([0 0 1], planeDir)/norm(cross([0 0 1], planeDir));
if ~isequal(N, [0 0 0])
    if N(2) >= 0
        Om = acos(N(1));
    else
        Om = 2*pi - acos(N(1));
    end
else
    Om = 0;
end

% evaluate impact parameters
DELTAMin = rp*sqrt(1 + 2*mu /(rp*vMin^2));
DELTAPlus = rp*sqrt(1 + 2*mu /(rp*vPlus^2));

% evaluate two branches specific angular momentums
hMin = DELTAMin*vMin; 
hPlus = DELTAPlus*vPlus; 
hMin = hMin*planeDir;
hPlus = hPlus*planeDir; 

% as vInfMin is the velocity at infinity direction, it coincides with radial direction at infinite
rHatInfMin = -vInfMin/norm(vMin);
rHatInfPlus = vInfPlus/norm(vPlus);

% eccentricity vectors evaluations
eMinVect = cross(vInfMin, hMin)/mu - rHatInfMin; 
ePlusVect = cross(vInfPlus, hPlus)/mu - rHatInfPlus; 

% pericentre anomaly evaluation
if eMinVect(3) >= 0
    om = acos(dot(N, eMinVect)/norm(eMin)); 
else
    om = 2*pi - acos(dot(N, eMinVect)/norm(eMin));
end

% infinite true anomalies evaluation
thetaInfMin = acos(-1/eMin); 
thetaInfPlus = acos(-1/ePlus); 

thvec1 = linspace(-thetaInfMin+0.4,0,10000);
thvec2 = linspace(0,thetaInfPlus-0.4,10000);

for ii = 1:10000
    
    orb1 = [aMin, eMin, i, Om, om, thvec1(ii)];
    [rr1(:,ii),~] = par2car(orb1,muE);     
    orb2 = [aPlus, ePlus, i, Om, om, thvec2(ii)];
    [rr2(:,ii),~] = par2car(orb2,muE);
    
end

Rsoi = astroConstants(2) * (muE/muS)^(2/5)

figure
plot3(rr1(1,:),rr1(2,:),rr1(3,:), 'LineWidth', 1.5); hold on;
plot3(rr2(1,:),rr2(2,:),rr2(3,:), 'LineWidth', 1.5);

[flybykep,~] = uplanet(FlyByDate, flyByPlanID);
[rr6, vv6] = par2car(flybykep, muS);
quiver3(0, 0, 0, -rr6(1)/10000, -rr6(2)/10000, -rr6(3)/10000+6500);

[X2,Y2,Z2] = sphere(100);
X2 = r2*X2 ;
Y2 = r2*Y2 ;
Z2 = r2*Z2 ;
fbp = surf(X2 ,Y2 ,Z2 ,'FaceColor', 'none', 'EdgeColor', 'none');
set(fbp, 'FaceColor', 'texturemap', 'CData', cdataE);
hold on; axis equal; grid on;

% [X2,Y2,Z2] = sphere(100);
% X2 = Rsoi*X2 ;
% Y2 = Rsoi*Y2 ;
% Z2 = Rsoi*Z2 ;
% surf(X2 ,Y2 ,Z2 ,'FaceColor', [.1, .1, .1], 'EdgeColor', 'none', 'FaceAlpha', 0.2)


F = @(th, e) 2*atanh(sqrt((e-1)/(e+1)) * tan(th/2));
dt = @(th1, th2, e, a) sqrt(muE/(-a^3)) * ( (e*sinh(F(th2, e)) - F(th2, e)) - (e*sinh(F(th1, e)) - F(th1, e)) );

[a, index1] = min(abs(vecnorm(rr1) - Rsoi))
[a, index2] = min(abs(vecnorm(rr2) - Rsoi))

% plot3(rr1(1, index1),rr1(2, index1),rr1(3, index1), '*')

thvec1(index1)

dt(0, 2*pi - thvec1(index1), eMin, aMin)
dt(0, thvec2(index2), ePlus, aPlus)



csad = 1;





































































































































