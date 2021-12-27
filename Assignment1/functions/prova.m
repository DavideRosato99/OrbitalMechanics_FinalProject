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
deltaDdep = data.timeWindows.deltaDdep;
nTOF1     = data.timeWindows.nTOF1;
nTOF2     = data.timeWindows.nTOF2;

%% FIRST INITIALIZATION
TOTdays = (datenum(arrDate) - datenum(depDate));
daysFromDep = 0 : deltaDdep : TOTdays;

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
synPeriod1 = (periodFlyByPlanet * periodDepPlanet)/(abs(periodFlyByPlanet - periodDepPlanet));
synPeriod2 = (periodArrPlanet * periodFlyByPlanet)/(abs(periodArrPlanet - periodFlyByPlanet));

%% TOF VECTORS
TOF1vec = linspace(0, synPeriod1, nTOF1);
TOF2vec = linspace(0, synPeriod2, nTOF2);

%% CALCULATION
DVfirstLeg = nan(length(daysFromDep), nTOF1);
DTfirstLeg = nan(length(daysFromDep), nTOF1);
DVsecondLeg = nan(length(daysFromDep), nTOF1, nTOF2);

length(daysFromDep)
parfor i = 1:length(daysFromDep)
    disp(i)
    depDateLoop = depDate;
    depDateLoop(3) = depDateLoop(3) + daysFromDep(i);
    
    [Yl, Mol, Dl] = ymd(datetime(depDateLoop));
    [Hl, Ml, Sl]  = hms(datetime(depDateLoop));
    
    datemjd2000l = date2mjd2000([Yl Mol Dl Hl Ml Sl]);
    
    % Departure planet
    [kepDep, ksun] = uplanet(datemjd2000l, depPlanID);
    [rrDep, vvDep] = par2car(kepDep, ksun);
    
    for j = 1:nTOF1
        
        FlyByDateLoop = depDateLoop;
        FlyByDateLoop(6) = FlyByDateLoop(6) + TOF1vec(j);
        
        [YL, MoL, DL] = ymd(datetime(FlyByDateLoop));
        [HL, ML, SL]  = hms(datetime(FlyByDateLoop));
        
        datemjd2000L = date2mjd2000([YL MoL DL HL ML SL]);
        
        controlDate = datenum(arrDate) - datenum(FlyByDateLoop);
        
        if controlDate >= 0
            
            [kepFlyBy, ksun] = uplanet(datemjd2000L, flyByPlanID);
            [rrFlyBy, vvFlyBy] = par2car(kepFlyBy, ksun);
            
            [~, ~, ~, err, vI, vF, ~, ~] = lambertMR(rrDep, rrFlyBy, TOF1vec(j), ksun, 0, 0, 0, 1);
            
            
            if err == 0
                
                DV1 = norm(vI' - vvDep);
                DV2 = norm(vF' - vvFlyBy);
                DVfirstLeg(i, j) = DV1 + DV2;
                DTfirstLeg(i, j) = daysFromDep(i)*24*3600 + TOF1vec(j);
                
                for k = 1:nTOF2
                    
                    arrDateLoop = FlyByDateLoop;
                    arrDateLoop(6) = arrDateLoop(6) + TOF2vec(k);
                    
                    [YLL, MoLL, DLL] = ymd(datetime(arrDateLoop));
                    [HLL, MLL, SLL]  = hms(datetime(arrDateLoop));

                    datemjd2000Ll = date2mjd2000([YLL MoLL DLL HLL MLL SLL]);

                    controlDateL = datenum(arrDate) - datenum(arrDateLoop);
                    
                    if controlDateL >= 0
                        
                        [kepArr, ksun] = uplanet(datemjd2000Ll, arrPlanID);
                        [rrArr, vvArr] = par2car(kepArr, ksun);

                        [~, ~, ~, err, vI, vF, ~, ~] = lambertMR(rrFlyBy, rrArr, TOF2vec(j), ksun, 0, 0, 0, 1);
                        
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

%% MATRICES PARSING
DVLeg1 = [];
controlArrFlyBy = [];

for i = 1:length(daysFromDep)
    for j = 1:nTOF1
        arrFlyBy = daysFromDep(i)*24*3600 + TOF1vec(j);
        
        if any(controlArrFlyBy == arrFlyBy)
        else
            controlArrFlyBy(end+1) = arrFlyBy;
            [row, col] = find(DTfirstLeg == arrFlyBy);
            
            vec = [];
            
            for k = 1:length(row)
                vec(end + 1) = DVfirstLeg(row(k), col(k));
            end
            
            DVLeg1(i, end+1) = min(vec);
            
        end
        
    end
end

[~, indexes] = sort(controlArrFlyBy);
DVLeg1 = DVLeg1(:, indexes);

%%
data.timeWindows.DVLeg1 = DVLeg1;
data.timeWindows.DVsecondLeg = DVsecondLeg;


end