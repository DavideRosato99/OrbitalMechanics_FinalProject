function [data] = optimizationGA(data, settings)

%%
minDate = data.timeWindows.depDate;
maxDate = data.timeWindows.arrDate;
depPlanID   = data.starting.depPlanID;
flyByPlanID = data.starting.flyByPlanID;
arrPlanID   = data.starting.arrPlanID;
minHfl = data.optimization.minHfl;
TOF1max = data.timeWindows.maxTOF1days;
TOF2max = data.timeWindows.maxTOF2days;

%%
minDatemjd2000 = date2mjd2000(minDate);

%% GA
A = [1 1 1];
B = [(datenum(maxDate) - datenum(minDate))*24*60];

LB(1:3) = 0;
UB(1) = (datenum(maxDate) - datenum(minDate))*24*60;
UB(2) = TOF1max*24*60;
UB(3) = TOF2max*24*60;

fitnessFcn = @(x) fitness(x, A, B, minHfl, minDatemjd2000, depPlanID, flyByPlanID, arrPlanID);

options = optimoptions('ga', 'MaxStallGenerations', 50, 'FunctionTolerance', ...
    1e-6, 'MaxGenerations', 10, 'NonlinearConstraintAlgorithm', 'penalty',...
    'PopulationSize', 5000, 'PlotFcn', {@(options,state,flag) gaplotsafe(options,state,flag),'gaplotbestindiv', 'gaplotbestf'}, ...
    'Display', 'iter', 'EliteCount', floor(0.01*5000), 'UseParallel', false);
[xOpt, DVopt, exitflag] = ga(fitnessFcn, 3, [], [], [], [], LB, UB, [], 1:3, options)
% delete('gaprogress.txt');

end


%% FUNCTIONS
%%% FITNESS FUNCTION
function DV = fitness(x, A, B, minHfl, minDatemjd2000, depPlanID, flyByPlanID, arrPlanID)

contr = sum(A.*x);

if contr <= B
    [kepDep, ksun] = uplanet(minDatemjd2000 + x(1)/(24*60), depPlanID);
    [rrDep, vvDep] = par2car(kepDep, ksun);

    [kepFlyBy, ksun] = uplanet(minDatemjd2000 + x(2)/(24*60) + x(1)/(24*60), flyByPlanID);
    [rrFlyBy, vvFlyBy] = par2car(kepFlyBy, ksun);

    [kepArr, ksun] = uplanet(minDatemjd2000 + x(3)/(24*60) + x(2)/(24*60) + x(1)/(24*60), arrPlanID);
    [rrArr, vvArr] = par2car(kepArr, ksun);

    [~, ~, ~, error1, ~, ~, ~, ~, ~, ~, ~, ...
    error2, ~, ~, ~, ~, DV1, DV2, errorFB, ~, ~, ~, ...
    delta_V_powFB, ~, ~, ~, ~] = deltaVtot(...
    rrDep, rrFlyBy, rrArr, vvDep, vvFlyBy, vvArr, x(2)*60,...
    x(3)*60, minHfl, flyByPlanID);

    if error1 == 0 && error2 == 0 && errorFB == 0
        DV = DV1 + DV2 + delta_V_powFB;
    else
        DV = NaN;
    end
else
    DV = NaN;
end

end







