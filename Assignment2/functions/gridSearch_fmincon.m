function [data] = gridSearch_fmincon(data, settings)

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
maxDatemjd2000 = date2mjd2000(maxDate);

%%
bestGlobal = 10e9;
error = 10;
flag  = false;

departureStart = 0;
TOF1Start      = 0;
TOF2Start      = 0;

departureEnd = (datenum(maxDate) - datenum(minDate));
TOF1End      = (datenum(maxDate) - datenum(minDate));
TOF2End      = (datenum(maxDate) - datenum(minDate));

ctr = 1;

while error >= 0.5 && not(flag) && ctr < 2
    
    t = tic;
    best = 10e9;
    bestI = NaN;
    bestJ = NaN;
    bestK = NaN;
    
    DVTOT = nan(nDep, nTOF1, nTOF2);
    
    ctrI = 1;
    for i = linspace(departureStart, departureEnd, nDep)
        if minDatemjd2000 + i <= maxDatemjd2000
            ctrJ = 1;
            for j = linspace(TOF1Start, TOF1End, nTOF1)
                if minDatemjd2000 + i + j <= maxDatemjd2000
                    ctrK = 1;
                    for k = linspace(TOF2Start, TOF2End, nTOF2)
                        if minDatemjd2000 + i + j + k <= maxDatemjd2000
                            [kepDep, ksun] = uplanet(minDatemjd2000 + i, depPlanID);
                            [rrDep, vvDep] = par2car(kepDep, ksun);

                            [kepFlyBy, ksun] = uplanet(minDatemjd2000 + j + i, flyByPlanID);
                            [rrFlyBy, vvFlyBy] = par2car(kepFlyBy, ksun);

                            [kepArr, ksun] = uplanet(minDatemjd2000 + k + j + i, arrPlanID);
                            [rrArr, vvArr] = par2car(kepArr, ksun);
                            
                            [~, ~, ~, error1, ~, ~, ~, ~, ~, ~, ~, ...
                            error2, ~, ~, ~, ~, DV1, DV2, errorFB, ~, ~, ~, ...
                            delta_V_powFB, ~, ~, ~, ~] = deltaVtot(...
                            rrDep, rrFlyBy, rrArr, vvDep, vvFlyBy, vvArr, j*24*3600,...
                            k*24*3600, minHfl, flyByPlanID);
                            
                            if error1 == 0 && error2 == 0 && errorFB == 0
                                DVTOT(ctrI, ctrJ, ctrK) = DV1 + DV2 + delta_V_powFB;

                                if DVTOT(ctrI, ctrJ, ctrK) <= best
                                    best = DVTOT(ctrI, ctrJ, ctrK);
                                    bestI = i;
                                    bestJ = j;
                                    bestK = k;
                                end
                            end
                        
                        end
                        ctrK = ctrK + 1;
                    end
                end
                ctrJ = ctrJ + 1;
            end
        end
        ctrI = ctrI + 1;
    end
    
    DVTOTglobal{ctr} = DVTOT;
    departureStartGlobal{ctr} = departureStart;
    departureEndGlobal{ctr}   = departureEnd;
    TOF1StartGlobal{ctr} = TOF1Start;
    TOF1EndGlobal{ctr}   = TOF1End;
    TOF2StartGlobal{ctr} = TOF2Start;
    TOF2EndGlobal{ctr}   = TOF2End;
    compTime(ctr) = toc(t);
    
    ctr = ctr + 1;
   
    error = abs(bestGlobal - best)
    
    if best < bestGlobal
        
        bestGlobal = best
        bestGlobalI = bestI;
        bestGlobalJ = bestJ;
        bestGlobalK = bestK;
        
        deltaDep  = departureEnd - departureStart;
        deltaTOF1 = TOF1End - TOF1Start;
        deltaTOF2 = TOF2End - TOF2Start;
        
        if not( bestGlobalI - deltaDep/4 < departureStart )
            departureStart = bestGlobalI - deltaDep/4;
        end
        if not ( bestGlobalI + deltaDep/4 > departureEnd )
            departureEnd   = bestGlobalI + deltaDep/4;
        end
        
        if not( bestGlobalJ - deltaTOF1/4 < TOF1Start )
            TOF1Start = bestGlobalJ - deltaTOF1/4;
        end
        if not ( bestGlobalJ + deltaTOF1/4 > TOF1End )
            TOF1End   = bestGlobalJ + deltaTOF1/4;
        end
        
        if not( bestGlobalK - deltaTOF2/4 < TOF2Start  )
            TOF2Start = bestGlobalK - deltaTOF2/4;
        end
        if not ( bestGlobalK + deltaTOF2/4 > TOF2End )
            TOF2End   = bestGlobalK + deltaTOF2/4;
        end
        
    else
        
        flag = true;
        ctr = ctr-1;
        
    end
    
end

%% FMINCON

x0 = [bestGlobalI, bestGlobalJ, bestGlobalK];

A = [1 1 1];
B = [(datenum(maxDate) - datenum(minDate))];

LB(1:3) = 0;
UB(1:3) = (datenum(maxDate) - datenum(minDate));


fitnessFcn = @(x) fitness(x, minHfl, minDatemjd2000, depPlanID, flyByPlanID, arrPlanID);
nonLinConFcn = @(x) nonlinfcn(x, minHfl, minDatemjd2000, depPlanID, flyByPlanID, arrPlanID);

options = optimoptions('fmincon', 'Algorithm', 'interior-point',...
    'MaxFunctionEvaluations', 5e6, 'OptimalityTolerance', 1e-9,...
    'StepTolerance', 1e-30, 'ConstraintTolerance', 1e-9, ...
    'DiffMinChange', 1e-5, 'FunctionTolerance', 1e-10, 'OutputFcn', @fminconOutData);
[xOpt, DVopt] = fmincon(fitnessFcn, x0, A, B, [], [], LB, UB, nonLinConFcn, options)


end


%% FUNCTIONS
%%% FITNESS FUNCTION
function DV = fitness(x, minHfl, minDatemjd2000, depPlanID, flyByPlanID, arrPlanID)
    
[kepDep, ksun] = uplanet(minDatemjd2000 + x(1), depPlanID);
[rrDep, vvDep] = par2car(kepDep, ksun);

[kepFlyBy, ksun] = uplanet(minDatemjd2000 + x(2) + x(1), flyByPlanID);
[rrFlyBy, vvFlyBy] = par2car(kepFlyBy, ksun);

[kepArr, ksun] = uplanet(minDatemjd2000 + x(3) + x(2) + x(1), arrPlanID);
[rrArr, vvArr] = par2car(kepArr, ksun);
                            
[~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ...
~, ~, ~, ~, ~, DV1, DV2, ~, ~, ~, ~, ...
delta_V_powFB, ~, ~, ~, ~] = deltaVtot(...
rrDep, rrFlyBy, rrArr, vvDep, vvFlyBy, vvArr, x(2)*24*3600,...
x(3)*24*3600, minHfl, flyByPlanID);

DV = DV1 + DV2 + delta_V_powFB;

end

%%% NON LINEAR CONSTRAINT FUNCTION
function [c, ceq] = nonlinfcn(x, minHfl, minDatemjd2000, depPlanID, flyByPlanID, arrPlanID)

ceq = [];
[kepDep, ksun] = uplanet(minDatemjd2000 + x(1), depPlanID);
[rrDep, vvDep] = par2car(kepDep, ksun);

[kepFlyBy, ksun] = uplanet(minDatemjd2000 + x(2) + x(1), flyByPlanID);
[rrFlyBy, vvFlyBy] = par2car(kepFlyBy, ksun);

[kepArr, ksun] = uplanet(minDatemjd2000 + x(3) + x(2) + x(1), arrPlanID);
[rrArr, vvArr] = par2car(kepArr, ksun);
                            
[~, ~, ~, error1, ~, ~, ~, ~, ~, ~, ~, ...
error2, ~, ~, ~, ~, ~, ~, errorFB, ~, ~, ~, ...
~, ~, ~, ~, ~] = deltaVtot(rrDep, rrFlyBy, rrArr, vvDep, vvFlyBy, vvArr,...
x(2)*24*3600, x(3)*24*3600, minHfl, flyByPlanID);

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



