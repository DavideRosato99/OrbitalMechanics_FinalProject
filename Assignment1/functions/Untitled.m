clear all
close all
clc

J2 = astroConstants(9);
Re = astroConstants(23);
mu = astroConstants(13);

k = 2;
m = 1;
e = 0.74;
i = 63.4*pi/180;
omE = deg2rad(15.04)/3600;
n = omE * (k/m);


Omdot = @(a) -((3/2) * (sqrt(mu)*J2*Re^2)/((1-e^2)^2 * a^(7/2)))*cos(i);
omdot = @(a) -((3/2) * (sqrt(mu)*J2*Re^2)/((1-e^2)^2 * a^(7/2)))*((5/2)*sin(i)^2 - 2);
M0dot = @(a) ((3/2) * (sqrt(mu)*J2*Re^2)/((1-e^2)^(3/2) * a^(7/2)))*(1 - (3/2)*sin(i)^2);
n = @(a) sqrt(mu/a^3);

funPert = @(a) ((m/k) - (omE - Omdot(a))./(n(a) + omdot(a) + M0dot(a)));


[X,FVAL] = fsolve(funPert, 26600, optimoptions('fsolve','FunctionTolerance',1e-16,...
'OptimalityTolerance',1e-16))
