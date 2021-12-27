function data = prova(data, settings)
%
%
%
%


%% SET UP VARIABLES
%%% STARTING
depDate     = data.starting.depDate;
arrDate     = data.starting.arrDate;
depPlanID   = data.starting.depPlanID;
flyByPlanID = data.starting.flyByPlanID;
arrPlanID   = data.starting.arrPlanID;
%%% TIME WINDOWS
nDep  = data.timeWindows.nDep;
nTOF1 = data.timeWindows.nTOF1;
nTOF2 = data.timeWindows.nTOF2;

%% FIRST INITIALIZATION
TOTdays = (datenum(arrDate) - datenum(depDate));
secFromDep = linspace(0, TOTdays*24*3600, nDep);

%% PLANETS CHARACTERIZATION
datemjd2000 = date2mjd2000(depDate);
[kepDep, ksun] = uplanet(datemjd2000, depPlanID);
[kepFlyBy, ~]  = uplanet(datemjd2000, flyByPlanID);
[kepArr, ~]    = uplanet(datemjd2000, arrPlanID);

%% SYNODIC PERIOD
% Planets period
periodDepPlanet   = 2*pi * sqrt(kepDep(1)^3 / ksun);
periodFlyByPlanet = 2*pi * sqrt(kepFlyBy(1)^3 / ksun);
periodArrPlanet   = 2*pi * sqrt(kepArr(1)^3 / ksun);

% Synodic period
synPeriod1 = (periodFlyByPlanet * periodDepPlanet)/(abs(periodFlyByPlanet - periodDepPlanet))
synPeriod2 = (periodArrPlanet * periodFlyByPlanet)/(abs(periodArrPlanet - periodFlyByPlanet));

%% TOF VECTORS
arrFlyByVec = linspace(0, TOTdays*24*3600, nTOF1);
arrArrVec   = linspace(0, TOTdays*24*3600, nTOF2);

%% CALCULATION
DVfirstLeg = nan(nDep, nTOF1);
DVsecondLeg = nan(nDep, nTOF1, nTOF2);

nDep
for i = 1:nDep
    disp(i)
    depDateLoop = depDate;
    depDateLoop(6) = depDateLoop(6) + secFromDep(i);
    
    [Yl, Mol, Dl] = ymd(datetime(depDateLoop));
    [Hl, Ml, Sl]  = hms(datetime(depDateLoop));
    
    datemjd2000l = date2mjd2000([Yl Mol Dl Hl Ml Sl]);
    
    % Departure planet
    [kepDep, ksun] = uplanet(datemjd2000l, depPlanID);
    [rrDep, vvDep] = par2car(kepDep, ksun);
    
    for j = 1:nTOF1
        
        FlyByDateLoop = depDate;
        FlyByDateLoop(6) = FlyByDateLoop(6) + arrFlyByVec(j);
        
        controlCons = datenum(FlyByDateLoop) - datenum(depDateLoop);
        
        if controlCons >= 0
            
            controlDate = datenum(arrDate) - datenum(FlyByDateLoop);
        
            if controlDate >= 0
                
                [YL, MoL, DL] = ymd(datetime(FlyByDateLoop));
                [HL, ML, SL]  = hms(datetime(FlyByDateLoop));
                datemjd2000L = date2mjd2000([YL MoL DL HL ML SL]);
                
                [kepFlyBy, ksun] = uplanet(datemjd2000L, flyByPlanID);
                [rrFlyBy, vvFlyBy] = par2car(kepFlyBy, ksun);

                TOF1 = (datenum(FlyByDateLoop) - datenum(depDateLoop))*24*3600;
                [~, ~, ~, err, vI, vF, ~, ~] = lambertMR(rrDep, rrFlyBy, TOF1, ksun, 0, 0, 0, 1);


                if err == 0

                    DV1 = norm(vI' - vvDep);
                    DV2 = norm(vF' - vvFlyBy);
                    DVfirstLeg(i, j) = DV1 + DV2;

                    for k = 1:nTOF2

                        arrDateLoop = depDate;
                        arrDateLoop(6) = arrDateLoop(6) + arrArrVec(k);

                        [YLL, MoLL, DLL] = ymd(datetime(arrDateLoop));
                        [HLL, MLL, SLL]  = hms(datetime(arrDateLoop));

                        datemjd2000Ll = date2mjd2000([YLL MoLL DLL HLL MLL SLL]);

                        controlDateL = datenum(arrDate) - datenum(arrDateLoop);

                        if controlDateL >= 0

                            [kepArr, ksun] = uplanet(datemjd2000Ll, arrPlanID);
                            [rrArr, vvArr] = par2car(kepArr, ksun);

                            TOF2 = (datenum(arrDateLoop) - datenum(depDate))*24*3600;
                            [~, ~, ~, err, vI, vF, ~, ~] = lambertMR(rrFlyBy, rrArr, TOF2, ksun, 0, 0, 0, 1);

                            if err == 0

                                DV1 = norm(vI' - vvFlyBy);
                                DV2 = norm(vF' - vvArr);
                                DVsecondLeg(i, j, k) = DV1 + DV2;

                            end

                        end


                    end

                end


            end
        end
            
        
    end
    
end


%%
data.timeWindows.depDays = secFromDep/(24*3600);
data.timeWindows.DVfirstLeg = DVfirstLeg;
data.timeWindows.FlyByDays = arrFlyByVec/(24*3600);
data.timeWindows.DVsecondLeg = DVsecondLeg;
data.timeWindows.arrDays = arrArrVec/(24*3600);

%%
vec1 = '';
vec2 = linspace(0, TOTdays, 15);

for i = 1:length(vec2)
    date = depDate;
    date(3) = date(3) + floor(vec2(i));
    [YL, MoL, DL] = ymd(datetime(date));
    [HL, ML, SL]  = hms(datetime(date));
    vec1{end + 1} = datestr([YL MoL DL HL ML SL]);
end

vec3 = '';
vec4 = linspace(0, TOTdays, 10);

for i = 1:length(vec4)
    date = depDate;
    date(3) = date(3) + floor(vec4(i));
    [YL, MoL, DL] = ymd(datetime(date));
    [HL, ML, SL]  = hms(datetime(date));
    vec3{end + 1} = datestr([YL MoL DL HL ML SL]);
end


%%
daysSyn1 = synPeriod1/(24*3600);

if settings.timeWindows.plot
    
    [A, B] = meshgrid(data.timeWindows.depDays, data.timeWindows.FlyByDays);
    contour(A, B, data.timeWindows.DVfirstLeg', 100)
    xticks(vec2); xticklabels(vec1); xtickangle(-45)
    yticks(vec4); yticklabels(vec3)
    hold on
    plot([0 daysSyn1], [daysSyn1 0], 'r');
    plot([0 2*daysSyn1], [2*daysSyn1 0], 'r');
    plot([0 3*daysSyn1], [3*daysSyn1 0], 'r');
    plot([0 4*daysSyn1], [4*daysSyn1 0], 'r');
    plot([0 5*daysSyn1], [5*daysSyn1 0], 'r');
    plot([0 6*daysSyn1], [6*daysSyn1 0], 'r');
    plot([0 7*daysSyn1], [7*daysSyn1 0], 'r');
    plot([0 8*daysSyn1], [8*daysSyn1 0], 'r');
    plot([0 9*daysSyn1], [9*daysSyn1 0], 'r');
    plot([0 10*daysSyn1], [10*daysSyn1 0], 'r');
    
    figure
    plot(data.timeWindows.FlyByDays, data.timeWindows.DVfirstLeg(1, :))
    hold on
    plot([2*daysSyn1 2*daysSyn1], [0 max(data.timeWindows.DVfirstLeg(1, :))], 'r');
    plot([3*daysSyn1 3*daysSyn1], [0 max(data.timeWindows.DVfirstLeg(1, :))], 'r');
    plot([4*daysSyn1 4*daysSyn1], [0 max(data.timeWindows.DVfirstLeg(1, :))], 'r');
    plot([5*daysSyn1 5*daysSyn1], [0 max(data.timeWindows.DVfirstLeg(1, :))], 'r');
    
    
    
end


end