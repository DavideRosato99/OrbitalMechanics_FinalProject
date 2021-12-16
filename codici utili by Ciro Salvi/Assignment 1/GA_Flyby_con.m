function [DV_MIN, Dv_min_TOF_1, Dv_min_TOF_2, r1_arc, r2_arc, r3_arc, v_venus, t_venus, vInfMin, vInfPlus, dv_ga, Data, hp] = GA_Flyby_con( ga_it, t_start_j_1, t_end_j_2 )
% GA_Flyby_con starting from a window of dates this function minimize an objective function (dv_optMod2) with ga
%
% PROTOTYPE:
%	[DV_MIN, Dv_min_TOF_1, Dv_min_TOF_2, r1_arc, r2_arc, r3_arc, v_venus, t_venus,vInfMin, vInfPlus, dv_ga, Data, hp] = GA_Flyby_con( ga_it, t_start_j_1, t_end_j_2 )
%
% INPUT:
%	ga_it[1]                Number of ga iteration imposed
%	t_start_j_1[1]          Lower boundary value [mjd2000]
%   t_end_j_2[1]            Higher boundary value [mjd2000]
%
% OUTPUT :
%   DV_MIN[1]           	Minimum DV found by the optimizer [km/s]
%   Dv_min_TOF_1[1]     	TOF corresponding to the first transfer arc [days]
%   Dv_min_TOF_2[1]     	TOF corresponding to the second transfer arc [days]
%   r1_arc[3x1]        	 	Position of Mercury at the departure date from Mercury [km]
%   r2_arc[3x1]        		Position of Venus at the departure and arrival date on Venus [km]
%   r3_arc[3x1]        		Position of Jupiter at the arrival date on Jupiter [km]
%   v_venus[3x1]       	 	Venus velocity at flyby [km/s]
%   t_venus[1]          	Venus day at flyby [mjd2000]
%   vInfMin[3x1]            Infinite velocity before flyby [km/s]
%   vInfPlus[3x1]           Infinite velocity after flyby [km/s]
%   dv_ga[1]                DV to be given for the flyby [km/s]
%   Data[(ga_it*50)x4]      Data evaluated by ga: dates population [3xn][mjd2000] & Dv [km/s]
%	hp[1]                   Height of flyby perigee [km]
%
% CONTRIBUTORS:
%   Fabio Spada
%   Alessandro Staffolani
%   Suhailah Alkhawashke
%
% VERSIONS
%   2021-02-11
%

%Fitness function
ObjectiveFunction = @dv_optMod2;

% Number of variables : x(1) tof_1, x(2) tof_2, x(3) tof_3
nvars = 3;

% Costrained system of equation for TOFs
% Ax <= b :
% - Spacecraft can't arrive before  departing
%               x(1)-x(2)+0 <= 1e-5
%               0+x(2)-x(3) <= 1e-5
% - Spacecraft can't fly for more than a reasonable ToF
% Mercury-Venus              -x(1) + x(2) + 0 <= 30*30  % thirty months as ToF
% Venus-Jupiter                 0  - x(2) + x(3) <= 90*30 % ninety months as ToF

A = [1 -1 0; 0 1 -1 ; -1 1 0; 0 -1 1];
% b = [0; 0; 30*30; 90*30];
b = [1e-5; 1e-5; 30*30; 90*30];

% Upper bound and lower bound are first day and latest day of the given
% window for the mission
% min_TOF_1 = 40;
% min_TOF_2 = 450;
LB = [       t_start_j_1   (t_start_j_1+40)  (t_start_j_1+40+450) ];
UB = [(t_end_j_2-40-450)    (t_end_j_2-450)             t_end_j_2 ];

% Preallocation
vect_t(ga_it,3) = 0;
vect_dv(ga_it) = 0;
options = optimoptions('ga','Display','off','ConstraintTolerance',1e-9);

% We create a set of ga_it elements (obtained by running
% ga_it times the genetic algorithm)
Data = [];
flagErr = 0;
for i = 1:ga_it
    [t, dv_val, ~, ~,set_x,dv_evaluated] = ga(ObjectiveFunction,nvars,A,b,[],[],LB,UB,[],options);
    vect_t(i,:) = t;
    vect_dv(i) = dv_val;
    Data = [ Data; set_x, dv_evaluated ];
    
end

if ~isequal( isnan(Data(:, 4)), ones(size(Data, 1), 1))
    % We evaluate the smallest element of the population
    DV_MIN = min(vect_dv);
    T_DV_MIN = vect_t(vect_dv == DV_MIN,:);
    
    % Outputs DV_min and the position of the planets evaluated at the days
    % founded by the fmincon optimizator
    [DV_MIN, r1_arc, r2_arc, r3_arc, v_venus, vInfMin, vInfPlus,dv_ga, hp] = dv_optMod2(T_DV_MIN);
    
    % Evaluate the transfer times for Mercury -> Venus and Venus -> Jupiter
    Dv_min_TOF_1 = (T_DV_MIN(2)-T_DV_MIN(1))*86400;
    Dv_min_TOF_2 = (T_DV_MIN(3)-T_DV_MIN(2))*86400;
    
    % Day for the best flyby (relative to Venus)
    t_venus = T_DV_MIN(2);
    
else    
    fprintf('In the total population no option satisfied fly by constraints. \nAll output were initialized to NaN \n');
    DV_MIN = NaN;
    Dv_min_TOF_1 = NaN;
    Dv_min_TOF_2 = NaN;
    r1_arc = NaN;
    r2_arc = NaN;
    r3_arc = NaN;
    v_venus = NaN;
    t_venus = NaN;
    dv_ga = NaN;
    hp = NaN;
    vInfMin = NaN;
    vInfPlus = NaN;
end
end