clear all; close all; clc;
%
% InterplanetaryMission_group14.m
% 
% The following main allows for analyzing both approaches for the inteplanetary problem;
% The switch is the following:
strategy = 'LegApproach' % 'LegApproach' | 'GlobalApproach'
% 
% Some of the scripts recall already evaluated results; if a new analysis is
% preferred to be performed, triggering this switch is required:
switchNew = 0; 
% 
% The Leg approach consists of five phases, each can be triggered by its related flag.
% N.B.: the five sections can be executed independently only if switchNew is 0. 
switchCoarseLeg = 0; 
switchPCRefineLeg = 0; 
switchPCWindowLeg = 0; 
switchGAResultsLeg = 0; 
switchResPlotLeg = 0; 
% 
% The Global Approach consists of four phases, each can be triggered by its related flag.
% N.B.: the four sections can be executed independently only if switchNew is 0.
switchCoarseGlobal = 0;
switchResMergingGlobal = 0; 
switchRefinementGlobal = 0; 
switchTransferPlottingGlobal = 0; 
%
%
% CONTRIBUTORS:
% Fabio Spada
% Alessandro Staffolani
% 

addpath(genpath('.\Functions'))

%% LEG APPROACH %%% ----- COARSE RESULTS ----- %%% 

if strcmp(strategy, 'LegApproach') && (switchCoarseLeg)
Interplanetary_mission.m
end

%% LEG APPROACH %%% ----- PORKCHOP REFINE + DFT ANALYSIS ----- %%%  

if strcmp(strategy, 'LegApproach') && (switchPCRefineLeg)
porkchopRefined2leg1.m
porkchopRefined2leg2.m 
end

%% LEG APPROACH %%% ----- PORKCHOP WINDOW  ----- %%% 

if strcmp(strategy, 'LegApproach') && (switchPCWindowLeg)
overlap_main.m
end


%% LEG APPROACH %%% ----- GA RESULTS  ----- %%% 
if strcmp(strategy, 'LegApproach') && (switchPCWindowLeg)
ga_mainMod4.m
end

%% LEG APPROACH %%% ----- RESULTS PLOTTING  ----- %%% 
if strcmp(strategy, 'LegApproach') && (switchPCWindowLeg)
gifminima.m
end


%% GLOBAL APPROACH %%% ----- COARSE ANALYSIS ------
if strcmp(strategy, 'GlobalApproach') && (switchCoarseGlobal)
ga_mainMod3.m
end

%% GLOBAL APPROACH %%% ----- RESULTS MERGING AND PLOTTING -----
if strcmp(strategy, 'GlobalApproach') && (switchResMergingGlobal)
final_results.m
end
%% GLOBAL APPROACH %%% ----- REFINEMENT -----
if strcmp(strategy, 'GlobalApproach') && (switchRefinementGlobal)
refine_minima.m
end
%% GLOBAL APPROACH %%% ----- TRANSFER PLOTTING -----
if strcmp(strategy, 'GlobalApproach') && (switchTransferPlottingGlobal)
gif_mission.m
end

