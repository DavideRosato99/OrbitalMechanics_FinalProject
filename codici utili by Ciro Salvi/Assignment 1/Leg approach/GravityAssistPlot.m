function GravityAssistPlot(vInfMin, vInfPlus, hp, idAtt)
% GravityAssistPlot plots the gravity assist manoeuvre with the auxiliary fly by
% trajectories used in powFlyby.m to evaluate hp, given vInfMin and vInfPlus
%
% PROTOTYPE:
%    GravityAssistPlot(vInfMin, vInfPlus, hp, idAtt)
%
% INPUT:
%   vInfMin[3x1]      Infinite velocity before flyby [km/s]
%   vInfPlus[3x1]     Infinite velocity after flyby [km/s]
%	hp[1]             Height of flyby perigee [km]
%   idAtt[1]          indices of the planets to be plotted
%                        1:   Mercury
%                        2:   Venus
%                        3:   Earth
%                        4:   Mars
%                        5:   Jupiter
%                        6:   Saturn
%                        7:   Uranus
%                        8:   Neptune
%                        9:   Pluto
%
% CONTRIBUTORS:
%   Fabio Spada
%
% VERSIONS
%   2020-02-11
%
%% Real Trajectory

% extract velocities norms
vMin = norm(vInfMin); 
vPlus = norm(vInfPlus);

delta = acos(dot(vInfMin,vInfPlus)/(vMin*vPlus));

% extract attractor mass parameter and radius
mu = astroConstants(10+idAtt);
rP = astroConstants(20+idAtt);

% valuate epericentre height
rp = rP + hp;

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

% true trajectory
f = figure('color', 'none');
plotoptions.style = '-';
plotoptions.lw = 2; 
plotoptions.color = [0.07,0.62,1.00];
PlotOrbit(aMin, eMin, i, Om, om, mu, -thetaInfMin + 0.5, 0 , 0.001 ,'rad',plotoptions);
plotoptions.color = [1.00,0.41,0.16];
PlotOrbit(aPlus, ePlus, i, Om, om, mu, 0, thetaInfPlus - 0.5, 0.001, 'rad',plotoptions);

grid minor 
axis equal

view(planeDir)

% [~, vInfCheckMin] = kep2car(aMin, eMin, i, Om, om, -thetaInfMin, mu,'rad');
% [~, vInfCheckPlus] = kep2car(aPlus, ePlus, i, Om, om, thetaInfPlus,  mu,'rad');

Planet_plot(2, [0 0 0], 1);


% quiver3(0,0,0, planeDir(1)*5e5, planeDir(2)*5e5, planeDir(3)*5e5);
% quiver3(0,0,0, vInfMin(1)*2.5e5, vInfMin(2)*2.5e5, vInfMin(3)*2.5e5);
% quiver3(0,0,0, vInfPlus(1)*2.5e5, vInfPlus(2)*2.5e5, vInfPlus(3)*2.5e5);
% quiver3(0,0,0, vInfCheckMin(1)*2.5e5, vInfCheckMin(2)*2.5e5, vInfCheckMin(3)*2.5e5);
% quiver3(0,0,0, vInfCheckPlus(1)*2.5e5, vInfCheckPlus(2)*2.5e5, vInfCheckPlus(3)*2.5e5);

%% powFlyby.m trajectories - 1st trajectory

eGuess = 1/sin(delta/2);

% entry condition respectful trajectory
rpGuess1 = (eGuess - 1)*mu/vMin^2;
eMin1 = 1 + rpGuess1*vMin^2/mu;
ePlus1 = 1 + rpGuess1*vPlus^2/mu;


% evaluate impact parameters
DELTAMin1 = rpGuess1*sqrt(1 + 2*mu /(rpGuess1*vMin^2));
DELTAPlus1 = rpGuess1*sqrt(1 + 2*mu /(rpGuess1*vPlus^2));

% evaluate two branches specific angular momentums
hMin1 = DELTAMin1*vMin; 
hPlus1 = DELTAPlus1*vPlus; 
hMin1 = hMin1*planeDir;
hPlus1 = hPlus1*planeDir; 

% eccentricity vectors evaluations
eMinVect1 = cross(vInfMin, hMin1)/mu - rHatInfMin; 

% pericentre anomaly evaluation
if eMinVect1(3) >= 0
    om1 = acos(dot(N, eMinVect1)/norm(eMin1)); 
else
    om1 = 2*pi - acos(dot(N, eMinVect1)/norm(eMin1));
end

% first trajectory drawing
plotoptions.style = ':';
plotoptions.lw = 1.5; 
plotoptions.color = [0.07,0.62,1.00];
PlotOrbit(aMin, eMin1, i, Om, om1, mu, -thetaInfMin + 0.5, 0 , 0.001 ,'rad', plotoptions);
plotoptions.color = [1.00,0.41,0.16];
PlotOrbit(aPlus, ePlus1, i, Om, om1, mu, 0, thetaInfPlus - 0.5, 0.001, 'rad', plotoptions);

%% powFlyby.m trajectories - 2nd trajectory
% exit condition respectful trajectory

rpGuess2 = (eGuess - 1)*mu/vPlus^2;
eMin2 = 1 + rpGuess2*vMin^2/mu;
ePlus2 = 1 + rpGuess2*vPlus^2/mu;

% evaluate impact parameters
DELTAMin2 = rpGuess2*sqrt(1 + 2*mu /(rpGuess2*vMin^2));
DELTAPlus2 = rpGuess2*sqrt(1 + 2*mu /(rpGuess2*vPlus^2));

% evaluate two branches specific angular momentums
hMin2 = DELTAMin2*vMin; 
hPlus2 = DELTAPlus2*vPlus; 
hMin2 = hMin2*planeDir;
hPlus2 = hPlus2*planeDir; 

% eccentricity vectors evaluations
ePlusVect2 = cross(vInfPlus, hPlus2)/mu - rHatInfPlus; 

% pericentre anomaly evaluation
if ePlusVect2(3) >= 0
    om2 = acos(dot(N, ePlusVect2)/norm(ePlus2)); 
else
    om2 = 2*pi - acos(dot(N, ePlusVect2)/norm(ePlus2));
end

% second trajectory drawing
plotoptions.style = '--';
plotoptions.lw = 1.5; 
plotoptions.color = [0.07,0.62,1.00];
PlotOrbit(aMin, eMin2, i, Om, om2, mu, -thetaInfMin + 0.1, 0 , 0.001 ,'rad',plotoptions);
plotoptions.color = [1.00,0.41,0.16];
PlotOrbit(aPlus, ePlus2, i, Om, om2, mu, 0, thetaInfPlus - 0.1, 0.001, 'rad',plotoptions);

%% graphics modifications

Ax = gca;
%Ax.Color = 'none';
Ax.MinorGridColor = [0.15 0.15 0.15];
Ax.XColor = [0.15 0.15 0.15];
Ax.YColor = [0.15 0.15 0.15];
Ax.ZColor = [0.15 0.15 0.15];
Ax.TickLabelInterpreter = 'latex';
Ax.MinorGridAlpha = 0.2;
Ax.OuterPosition = [-0.171192314745643,-0.025427831894398,1.086094728812638,1.023369949880217];
Ax.InnerPosition = [-0.03,0.087142862592426,0.841723414829795,0.834046509152377];
Ax.Position = [-0.03,0.087142862592426,0.841723414829795,0.834046509152377];

X = xlabel('$X_{VCI} $ [km]');
X.Interpreter = 'latex';

Y = ylabel('$Y_{VCI} $ [km]');
Y.Interpreter = 'latex';

Z = zlabel('$Z_{VCI} $ [km]');
Z.Interpreter = 'latex';

Leg = legend;
Leg.Color = 'none';
Leg.Interpreter = 'latex';
Leg.TextColor = [0.15 0.15 0.15];
Leg.String = {'Entry leg - Real'; 'Exit leg - Real'; 'Entry leg - Guess 1'; 'Exit leg - Guess 1'; 'Entry leg - Guess 2'; 'Exit leg - Guess 2'; 'Venus heliocentric velocity'};
Leg.FontSize = 9;
Leg.Position = [0.659530714280776,0.686030651940252,0.313145081501401,0.266904751550583];
Leg.Box = 'off';