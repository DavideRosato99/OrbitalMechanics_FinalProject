function [Data_min, Dv_min_TOF_1, Dv_min_TOF_2, r1_arc, r2_arc, r3_arc, v_venus, t_venus, vInfMin, vInfPlus, dv_ga, Data, hp] = GA_FlybyMod3( ga_it, t_merc, t_ven_1, t_ven_2, t_jup )
% GA_FlybyMod3 starting from a set of dates this function minimize an objective
% function (dv_optMod2) with ga shifting only the dates on Venus
%
% PROTOTYPE:
%	[Data_min, Dv_min_TOF_1, Dv_min_TOF_2, r1_arc, r2_arc, r3_arc, v_venus, t_venus, vInfMin, vInfPlus, dv_ga, Data, hp] = GA_FlybyMod3( ga_it, t_merc, t_ven_1, t_ven_2, t_jup )
%
% INPUT:
%	ga_it[1]                   Number of ga iteration imposed
%	t_merc[1]        		Mercury departure date [mjd2000]
%	t_ven_1[1]             Lower boundary of Venus flyby date [mjd2000]
%	t_ven_2[1]             Lower boundary of Venus flyby date [mjd2000]
%	t_merc[1]       	    Jupiter arrival date [mjd2000]

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
%   Data[4x(ga_it*50)]      Data evaluated by ga: dates population [3xn][mjd2000] & Dv [km/s]
%	hp[1]                   Height of flyby perigee [km]
%
% CONTRIBUTORS:
%   Alessandro Staffolani
%
% VERSIONS
%   2021-02-11
%

%Fitness function
ObjectiveFunction = @dv_optMod2;

% Number of variables : x(1) tof_1, x(2) tof_2, x(3) tof_3
nvars = 3;
% Inequality condition: t_ven_1 <= x(2) <= t_ven_2
%             0-x(2)+0 = -t_ven_1
%             0+x(2)+0 =  t_ven_2
A = [ 0, -1, 0; 0, 1, 0 ];
b = [ -t_ven_1; t_ven_2 ];
% Equality condition:
%             x(1)+0+0 = t_merc
%             0+0+x(3) = t_jup
Aeq = [ 1, 0, 0; 0, 0, 1 ];
beq = [ t_merc; t_jup ];

% Upper bound and lower bound are first day and latest day of the given
% window for the mission
% min_TOF_1 = 40;
% min_TOF_2 = 450;
% LB = [ -Inf t_ven_1 -Inf];
% UB = [  Inf t_ven_2  Inf];

% Preallocation
vect_t(ga_it,3) = 0;
vect_dv(ga_it) = 0;
options = optimoptions('ga','Display','off'); %,'ConstraintTolerance',1e-9

% We create a set of ga_it elements (obtained by running
% ga_it times the genetic algorithm)
Data = [];
for i = 1:ga_it
    [t, dv_val, ~, ~,set_x,dv_evaluated] = ga(ObjectiveFunction,nvars,A,b,Aeq,beq,[],[],[],options);
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
[DV_MIN, r1_arc, r2_arc, r3_arc, v_venus, vInfMin, vInfPlus, dv_ga, hp, DV] = dv_optMod2(T_DV_MIN);

% Evaluate the transfer times for Mars -> Saturn and Saturn -> Neptune
Dv_min_TOF_1 = (T_DV_MIN(2)-T_DV_MIN(1))*86400;
Dv_min_TOF_2 = (T_DV_MIN(3)-T_DV_MIN(2))*86400;

% Day for the best flyby (relative to saturn)
t_venus = T_DV_MIN(2);

Data_min = [ T_DV_MIN, DV_MIN, DV ];

else    
    fprintf('In the total population no option satisfied fly-by constraints. \nAll output were initialized to NaN \n');
    % é scritto col culo ma poi lo cambiamo perché non si capisce IN MEMORIAM
    Data_min = [NaN,NaN,NaN,NaN,NaN,NaN,NaN];
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