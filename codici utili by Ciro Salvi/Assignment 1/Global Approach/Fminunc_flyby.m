function [ Data_min, Dv_min_TOF_1, Dv_min_TOF_2, r_merc, r_ven, r_jup, v_venus, t_venus, dv_flyby, fminunc_struct] = Fminunc_flyby( x0 )
% Fminunc_Flyby starting from an initial guess this function minimize an objective function (dv_optMod2) with fminunc
%
% PROTOTYPE:
%	[ Data_min, Dv_min_TOF_1, Dv_min_TOF_2, r_merc, r_ven, r_jup, v_venus, t_venus, dv_flyby,fminunc_struct] = Fminunc_flyby( x0 )
%
% INPUT:
%   x[3x1]              	Vector of dates (initial guess):
%	  x(1)        			Mercury departure date [mjd2000]
%     x(2)        			Venus flyby date [mjd2000]
%     x(3)        	    	Jupiter arrival date [mjd2000]
%
% OUTPUT :
%   DV_MIN[1]           	Minimum DV found by the optimizer [km/s]
%   Dv_min_TOF_1[1]     	TOF corresponding to the first transfer arc [days]
%   Dv_min_TOF_2[1]     	TOF corresponding to the second transfer arc [days]
%   r_merc[3x1]        	 	Position of Mercury at the departure date from Mercury [km]
%   r_ven[3x1]        		Position of Venus at the departure and arrival date on Venus [km]
%   r_jup[3x1]        		Position of Jupiter at the arrival date on Jupiter [km]
%   v_venus[3x1]       	 	Venus velocity at flyby [km/s]
%   t_venus[1]          	Venus day at flyby [mjd2000]
%   dv_flyby[1]         	DV to be given for the flyby [km/s]
%	fminunc_struct[struct]  Struct of Optimizer info
%
% CONTRIBUTORS:
%   Alessandro Staffolani
%
% VERSIONS:
%   2021-02-11
%

% Options to remove output
options = optimoptions('fminunc');
options.Display = 'iter'; % 'none' || 'off' || 'iter' || 'iter-detailed' || 'notify' || 'notify-detailed' || 'final' || 'final-detailed'
options.MaxFunctionEvaluations = 5.000000e+03;
options.OptimalityTolerance = 1e-9;
options.StepTolerance = 1e-16;
options.DiffMinChange = 1e-1;
options.Algorithm = 'quasi-newton';
% options.Algorithm = 'trust-region';
% options.SpecifyObjectiveGradient = true;
% options.PrecondBandWidth = Inf; % uses a direct factorization (Cholesky) rather than the conjugate gradients (CG)

% Call the optimizator to find best days
[x,fval,exitflag,output] = fminunc(@dv_optMod2,x0,options);
fminunc_struct.fval = fval;
fminunc_struct.exitflag = exitflag;
fminunc_struct.output = output;

% Outputs DV_min and the position of the planets evaluated at the days
% founded by the fmincon optimizator
[DV_MIN, r_merc, r_ven, r_jup, v_venus, dv_flyby] = dv_optMod2(x);

% Evaluate the transfer times for Mercury -> Venus and Venus -> Jupiter
Dv_min_TOF_1 = (x(2)-x(1))*86400;
Dv_min_TOF_2 = (x(3)-x(2))*86400;

% Day for the best flyby (relative to saturn)
t_venus = x(2);
Data_min = [ x, DV_MIN ]; % [ Julian day, Julian day, Julian day, km/s ] Data of the minimum evaluated
end