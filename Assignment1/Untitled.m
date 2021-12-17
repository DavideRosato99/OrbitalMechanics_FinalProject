clear all
close all
clc


%% SET PATH
% global path
filePath = fileparts(mfilename('fullpath'));
currentPath = pwd;
if not(strcmp(filePath, currentPath))
    cd (filePath);
    currentPath = filePath;
end

addpath(genpath(currentPath));

% functions path
addpath(genpath('functions'))

%%

dateIn = [2027 12 1 0 0 0];
dateFin = [2067 12 1 0 0 0];

deltaTime = datetime(dateFin) - datetime(dateIn)

H = 350640;

% Mars
    [kepMars, ksun] = uplanet(date2mjd2000(dateIn), 4);
    Tmars = round(2*pi*sqrt(kepMars(1)^3/ksun));
    % Earth
    [kepEarth, ksun] = uplanet(date2mjd2000(dateIn), 3);
    Tearth = round(2*pi*sqrt(kepEarth(1)^3/ksun));
    % Venus
    [kepVenus, ksun] = uplanet(date2mjd2000(dateIn), 2);
    Tvenus = round(2*pi*sqrt(kepVenus(1)^3/ksun));
    
    Tme = (Tmars*Tearth) / abs(Tmars-Tearth);
    Tev = (Tvenus*Tearth) / abs(Tvenus-Tearth);
    
    mcm2 = lcm(round(Tev), round(Tme))

date = dateIn;

ctr = 1;

for h = 0:24:H
    datei = [date(1) date(2) date(3) date(4)+h date(5) date(6)];
    dateBello = datetime(datei);
    [Yl, Mol, Dl] = ymd(dateBello);
    [Hl, Ml, Sl] = hms(dateBello);
    
    datemjd2000 = date2mjd2000([Yl Mol Dl Hl Ml Sl]);
    % Mars
    [kepMars, ksun] = uplanet(datemjd2000, 4);
    [rrM(:, ctr), vvM(:, ctr)] = par2car(kepMars, ksun);
    % Earth
    [kepEarth, ksun] = uplanet(datemjd2000, 3);
    [rrE(:, ctr), vvE(:, ctr)] = par2car(kepEarth, ksun);
    % Venus
    [kepVenus, ksun] = uplanet(datemjd2000, 2);
    [rrV(:, ctr), vvV(:, ctr)] = par2car(kepVenus, ksun);

    ctr = ctr+1;
end

%%
figure
plot3(rrE(1, :), rrE(2, :), rrE(3, :)); hold on; axis equal; grid on;
plot3(rrM(1, :), rrM(2, :), rrM(3, :));
plot3(rrV(1, :), rrV(2, :), rrV(3, :));


E = plot3(rrE(1, 1), rrE(2, 1), rrE(3, 1), 'bo', 'MarkerSize', 3);
M = plot3(rrM(1, 1), rrM(2, 1), rrM(3, 1), 'ro', 'MarkerSize', 3);
V = plot3(rrV(1, 1), rrV(2, 1), rrV(3, 1), 'ko', 'MarkerSize', 3);

for i = 2:size(rrM, 2)
    delete(E)
    E = plot3(rrE(1, i), rrE(2, i), rrE(3, i), 'bo', 'MarkerSize', 3);
    delete(M)
    M = plot3(rrM(1, i), rrM(2, i), rrM(3, i), 'ro', 'MarkerSize', 3);
    delete(V)
    V = plot3(rrV(1, i), rrV(2, i), rrV(3, i), 'ko', 'MarkerSize', 3);
    
    drawnow limitrate
end

