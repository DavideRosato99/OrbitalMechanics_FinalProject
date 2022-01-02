function [data] = optimizationGA(data, settings)

%%
minDate = data.timeWindows.depDate;
maxDate = data.timeWindows.arrDate;
depPlanID   = data.starting.depPlanID;
flyByPlanID = data.starting.flyByPlanID;
arrPlanID   = data.starting.arrPlanID;
minHfl = data.optimization.minHfl;

nDep  = 30;
nTOF1 = 30;
nTOF2 = 30;

%%
minDatemjd2000 = date2mjd2000(minDate);

%% GA
A = [1 1 1];
B = [(datenum(maxDate) - datenum(minDate))*24*60];

LB(1:3) = 0;
UB(1:3) = (datenum(maxDate) - datenum(minDate))*24*60;

fitnessFcn = @(x) fitness(x, A, B, minHfl, minDatemjd2000, depPlanID, flyByPlanID, arrPlanID);
% nonLinConFcn = @(x) nonlinfcn(x, A, B, minHfl, minDatemjd2000, depPlanID, flyByPlanID, arrPlanID);

options = optimoptions('ga', 'MaxStallGenerations', 50, 'FunctionTolerance', ...
    1e-6, 'MaxGenerations', 200, 'NonlinearConstraintAlgorithm', 'penalty',...
    'PopulationSize', 500000, 'PlotFcn', {'gaplotbestindiv', 'gaplotbestf'}, ...
    'Display', 'iter', 'EliteCount', floor(0.01*500000), 'UseParallel', true);
[xOpt, DVopt, exitflag] = ga(fitnessFcn, 3, [], [], [], [], LB, UB, [], 1:3, options)


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

%%% NON LINEAR CONSTRAINT FUNCTION
function [c, ceq] = nonlinfcn(x, minHfl, minDatemjd2000, depPlanID, flyByPlanID, arrPlanID)

ceq = [];
[kepDep, ksun] = uplanet(minDatemjd2000 + x(1)/24, depPlanID);
[rrDep, vvDep] = par2car(kepDep, ksun);

[kepFlyBy, ksun] = uplanet(minDatemjd2000 + x(2)/24 + x(1)/24, flyByPlanID);
[rrFlyBy, vvFlyBy] = par2car(kepFlyBy, ksun);

[kepArr, ksun] = uplanet(minDatemjd2000 + x(3)/24 + x(2)/24 + x(1)/24, arrPlanID);
[rrArr, vvArr] = par2car(kepArr, ksun);
                            
[~, ~, ~, error1, ~, ~, ~, ~, ~, ~, ~, ...
error2, ~, ~, ~, ~, ~, ~, errorFB, ~, ~, ~, ...
~, ~, ~, ~, ~] = deltaVtot(rrDep, rrFlyBy, rrArr, vvDep, vvFlyBy, vvArr,...
x(2)*3600, x(3)*3600, minHfl, flyByPlanID);

if error1 == 0 && error2 == 0 && errorFB == 0
    c = - 1;
else
    c = error1 + error2 + errorFB;
end

end

%%% SAVE DATA OUT OF FMINCON ITERATIONS FUNCTION
function stop = fminconOutData(x, optimValues, state)
stop = false;

x

optimValues
state



end



