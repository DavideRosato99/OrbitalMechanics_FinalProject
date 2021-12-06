%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                       ORBITAL MECHANICS                                 %
%                    Academic year 2020/2021                              %
%                                                                         %
%                  Module 3: Orbital manoeuvres                           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% DESCRIPTION: Sample script illustating the use of lambertMR
%
% Camilla Colombo, 11/11/2016
% Juan Luis Gonzalo, 17/11/2020: Additional comments
%

clc;
clearvars;

% Read Help of lambertMR for details. You can open file lambertMR.m in the
% editor, or use command 'help lambertMR'
%
% lambertMR(RI,RF,TOF,MU,orbitType,Nrev,Ncase,optionsLMR)
%
% As inputs you need only:
%  - RI: Vector containing the initial position in Cartesian coordinates [L]
%  - RF: Vector containing the final position vector in Cartesian coordinates [L]
%  - TOF: Transfer time, time of flight [T]
%  - MU: Planetary constant of the primary [L^3/T^2]
% For the other parameters set:
%  - orbitType = 0: Direct orbit (0: direct; 1: retrograde)
%  - Nrev = 0: Zero-revolution case
%  - Ncase = 0: Not used for the zero-revolution case
%  - optionsLMR = 0: No display


muSun = 0.39860e6;      % Sun's gravitational parameter [km^3/s^2];
ToF = 50*86400;         % Time in [s];
r1 = [10000,0,0];       % Initial position vector [km]
r2 = [20000,100,0];     % Final position vector [km]
[A,P,E,ERROR,VI,VF,TPAR,THETA] = lambertMR( r1, r2, ToF, muSun, 0, 0, 0 );