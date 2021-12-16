function [DV_MIN, Dv_min_TOF_1, Dv_min_TOF_2, r1_arc, r2_arc, r3_arc, v_venus, t_venus, vInfMin, vInfPlus, dv_ga, Data, hp] = GA_ref( ga_it_ref,ref_days,options_GAref, Data_min_unc)
% GA function used for the refinement of the results obtained from the
% Global GA analysis of the unconstrained values 

% PROTOTYPE:
%   [DV_MIN, Dv_min_TOF_1, Dv_min_TOF_2, r1_arc, r2_arc, r3_arc, v_venus, t_venus, vInfMin, vInfPlus, dv_ga, Data, hp] = GA_ref( ga_it_ref,ref_days,options_GAref, Data_min_unc)
% 
% INPUT:
%   ga_it_ref [1]       number of iterations for the ga function execution 
%   ref_days [1]        semi-amplitude value taken for the refinement [days]
%   options_GAref []    Optimization options used by the GA
%   Data_min_unc [1x4]  vector of Data of a given minimum  [1x3][mjd2000] and [km/s] 
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
%   Suhailah Alkhawashke
%
% VERSIONS:
%   2021-02-01: First version
%

% Extract the departure, flyby, and arrival dates
t_dep_opt = Data_min_unc(1);
t_Flyby_opt = Data_min_unc(2);
t_arr_opt = Data_min_unc(3);

%Fitness function
ObjectiveFunction = @dv_optMod2;

% Number of variables : x(1) tof_1, x(2) tof_2, x(3) tof_3
nvars_ref = 3;

% Constrained system of equation for TOFs (arrival after departure)
% Ax <= b :
%     x(1)-x(2)+0=0
%     0+x(2)-x(3)=0

A = [1, -1, 0; 0 , 1, -1];
b = [1e-5; 1e-5];

% Upper bound and lower bound are first day and last day of the given
% window for the mission, for the refinement they are taken by considering the ref_days 

LB = [t_dep_opt-ref_days   t_Flyby_opt-ref_days  t_arr_opt-ref_days];
UB = [t_dep_opt+ref_days   t_Flyby_opt+ref_days  t_arr_opt+ref_days];

% Preallocation
vect_t(ga_it_ref,3) = 0;
vect_dv(ga_it_ref) = 0;

% create a set of ga_it elements (obtained by running
% ga_it times the genetic algorithm)
Data = [];
for i = 1:ga_it_ref
    [t, dv_val, ~, ~,set_x,dv_evaluated] = ga(ObjectiveFunction,nvars_ref,A,b,[],[],LB,UB,[],options_GAref);
    vect_t(i,:) = t;
    vect_dv(i) = dv_val;
    Data = [ Data; set_x, dv_evaluated ];
end


if ~isequal( isnan(Data(:, 4)), ones(size(Data, 1), 1)) 

% evaluate the smallest element of the population
DV_MIN = min(vect_dv);
T_DV_MIN = vect_t(vect_dv == DV_MIN,:);

% Outputs DV_min and the position of the planets
[DV_MIN, r1_arc, r2_arc, r3_arc, v_venus,vInfMin, vInfPlus, dv_ga, hp] = dv_optMod2(T_DV_MIN);

% finding transfer times from Mercury to Venus and from Venus to Jupiter
Dv_min_TOF_1 = (T_DV_MIN(2)-T_DV_MIN(1))*86400;
Dv_min_TOF_2 = (T_DV_MIN(3)-T_DV_MIN(2))*86400;

% Day for the best flyby
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
