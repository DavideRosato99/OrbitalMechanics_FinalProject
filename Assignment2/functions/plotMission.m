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
% NOTE: time can be omitted
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
depDate      = [2041 08 01 01 08 00];
depDate      = date2mjd2000(depDate);
FlyByDate    = [2042 06 09 12 41 21];% data.starting.depDate;
FlyByDate    = date2mjd2000(FlyByDate);
% arrDate      = data.starting.arrDate;
arrDate      = [2042 12 06 05 19 24];
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
plot3(rrDin(1,:),rrDin(2,:),rrDin(3,:),'--k');
axis equal; grid on; hold on;
plot3(rrDfb(1,:),rrDfb(2,:),rrDfb(3,:),'--k');
plot3(rrDarr(1,:),rrDarr(2,:),rrDarr(3,:),'--k');

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
firstlegin = car2par(rr1,vi1',muS)

firstlegend= car2par(rr2,vf1',muS)
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

plot3(rr_firstleg(1,:),rr_firstleg(2,:),rr_firstleg(3,:));



%% second leg

secondlegin = car2par(rr2,vi2',muS)

secondlegend= car2par(rr3,vf2',muS)
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
figure; 
r = astroConstants(20+flyByPlanID);
[X,Y,Z] = sphere(100);
X = r*X ;
Y = r*Y ;
Z = r*Z ;
fbp = surf(X ,Y ,Z ,'FaceColor', 'none', 'EdgeColor', 'none');
set(fbp, 'FaceColor', 'texturemap', 'CData', cdataE);
hold on
axis equal
rsoi = astroConstants(2)*(muE/muS)^(2/5);

% [X,Y,Z] = sphere(100);
% X = rsoi*X ;
% Y = rsoi*Y ;
% Z = rsoi*Z ;
% fbp = surf(X ,Y ,Z,'FaceColor', [0.1 0.1 0.1], 'EdgeColor', 'none','FaceAlpha',0.2);

if firstlegin(6) < firstlegend(6)
    thspan = linspace(firstlegin(6),firstlegend(6),1000);
else
    thspan = linspace(firstlegin(6),firstlegend(6)+2*pi,1000);
end

tspan = linspace(0,TOF1,1000);
orb = firstlegin;
for i = 1:length(thspan)
    date = depDate + tspan(i)/(24*3600);
    [flybykep,~] = uplanet(date, flyByPlanID);
    [rrE, ~] = par2car(flybykep, muS);
    
    orb(6) = thspan(i);
    [rr_firstleg1(:,i),~] = par2car(orb,muS); 
    rr_firstleg(:,i) = rr_firstleg1(:,i) - rrE; 
 
end
[~,index] = min(abs(vecnorm(rr_firstleg) - rsoi))

% plot3(rr_firstleg(1,index-10:index),rr_firstleg(2,index-10:index),rr_firstleg(3,index-10:index));


[kepminus] = car2par(rr_firstleg(:,index),vf1',muE);
OM = kepminus(4);
om = kepminus(5);
ii  = kepminus(3);
R_OM = [cos(OM) sin(OM) 0; -sin(OM) cos(OM) 0; 0 0 1];
R_i =  [1 0 0; 0 cos(ii) sin(ii); 0 -sin(ii) cos(ii)];
R_om = [cos(om) sin(om) 0; -sin(om) cos(om) 0; 0 0 1];

T_ECI_PF = R_om*R_i*R_OM;
T_PF_ECI = T_ECI_PF';

rr=zeros(3,300);
thspan1 = [linspace(deg2rad(270),2*pi,300)];%,linspace(deg2rad(0),deg2rad(90),300)]
for i = 1:length(thspan1)
    r(i) = (a_minus*(1-e_minus^2))/(1+e_minus*cos(thspan1(i)));
    r_PF = r(i)*[cos(thspan1(i)); sin(thspan1(i)); 0];
    rr(:,i) = T_PF_ECI*r_PF;
end

plot3(rr(1,:),rr(3,:),rr(3,:),'LineWidth',1.5);
figure
plot(rad2deg(thspan1),r)

manna = 7;















































































































































