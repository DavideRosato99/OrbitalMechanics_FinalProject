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
synPeriod1 = (periodFlyByPlanet * periodDepPlanet)/(abs(periodFlyByPlanet - periodDepPlanet));
synPeriod2 = (periodArrPlanet * periodFlyByPlanet)/(abs(periodArrPlanet - periodFlyByPlanet));

%% TOF VECTORS
arrFlyByVec = linspace(0, TOTdays*24*3600, nTOF1);
arrArrVec   = linspace(0, TOTdays*24*3600, nTOF2);

%% LARGE GRID-SEARCH CALCULATIONS
DVfirstLeg = nan(nDep, nTOF1);
DVsecondLeg = nan(nDep, nTOF1, nTOF2);
DVsecondLeg2 = nan(nTOF1, nTOF2);

nDep
flag = false;
parfor i = 1:nDep
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
                [~, ~, ~, err, vI_1, vF_1, ~, ~] = lambertMR(rrDep, rrFlyBy, TOF1, ksun, 0, 0, 0, 1);

                if err == 0

                    for k = 1:nTOF2

                        arrDateLoop = depDate;
                        arrDateLoop(6) = arrDateLoop(6) + arrArrVec(k);

                        [YLL, MoLL, DLL] = ymd(datetime(arrDateLoop));
                        [HLL, MLL, SLL]  = hms(datetime(arrDateLoop));

                        datemjd2000Ll = date2mjd2000([YLL MoLL DLL HLL MLL SLL]);
                        
                        controlConsL = datenum(arrDateLoop) - datenum(FlyByDateLoop);
                        
                        if controlConsL >= 0
                        
                            controlDateL = datenum(arrDate) - datenum(arrDateLoop);

                            if controlDateL >= 0

                                [kepArr, ksun] = uplanet(datemjd2000Ll, arrPlanID);
                                [rrArr, vvArr] = par2car(kepArr, ksun);

                                TOF2 = (datenum(arrDateLoop) - datenum(FlyByDateLoop))*24*3600;
                                [~, ~, ~, err, vI_2, vF_2, ~, ~] = lambertMR(rrFlyBy, rrArr, TOF2, ksun, 0, 0, 0, 1);

                                if err == 0
                                    
                                    DV1_1 = norm(vI_1' - vvDep);
                                    DV2_1 = norm(vF_1' - vvFlyBy);
                                    DVfirstLeg(i, j) = DV1_1;

                                    DV1_2 = norm(vI_2' - vvFlyBy);
                                    DV2_2 = norm(vF_2' - vvArr);
                                    DVsecondLeg(i, j, k) = DV2_2;

                                end

                            end
                        
                        end


                    end

                end


            end
        end
            
        
    end
    
end


%% OBTAIN MINIMA FOR SECOND LEG
for j = 1:nTOF1
    for k = 1:nTOF2
        
        DV = DVsecondLeg(:, j, k);
        minDV = min(DV);
        DVsecondLeg2(j, k) = minDV;
        
    end
end


%%
minDV     = zeros(floor(TOTdays*24*3600/synPeriod2)+1, 1);
minDVdays = zeros(floor(TOTdays*24*3600/synPeriod2)+1, 1);
minima    = min(DVsecondLeg2, [], 2);

for i = 1:floor(TOTdays*24*3600/synPeriod2)+1
    indexes = find(arrFlyByVec >= (i-1)*synPeriod2 + 1 & arrFlyByVec < i*synPeriod2);
    [minDV(i), ind] = min(minima(indexes));
    ciao = arrFlyByVec(indexes);
    minDVdays(i) = ciao(ind)/(24*3600);
end

maxMinDV = max(minDV);
minMinDV = min(minDV);

f = @(DV) (maxMinDV - DV).^4/(maxMinDV - minMinDV)^4;

figure
plot(minDVdays, minDV); grid on


dateOptFlyBy = (sum(minDVdays(not(isnan(minDV))).*f(minDV(not(isnan(minDV))))))/(sum(f(minDV(not(isnan(minDV))))));


figure
plot(minDVdays, f(minDV)); hold on; grid on
xline(dateOptFlyBy)

%%%
DVfirstLeg2 = nan(nDep, 1);
secFromDep2 = linspace(0, dateOptFlyBy*24*3600, nDep);
for i = 1:nDep
    depDateLoop2 = depDate;
    depDateLoop2(6) = depDateLoop2(6) + secFromDep2(i);
    
    [Yl, Mol, Dl] = ymd(datetime(depDateLoop2));
    [Hl, Ml, Sl]  = hms(datetime(depDateLoop2));
    
    datemjd2000l = date2mjd2000([Yl Mol Dl Hl Ml Sl]);
    
    % Departure planet
    [kepDep, ksun] = uplanet(datemjd2000l, depPlanID);
    [rrDep, vvDep] = par2car(kepDep, ksun);
    
    % FlyBy planet
    FlyByDateLoop2 = depDate;
    FlyByDateLoop2(6) = FlyByDateLoop2(6) + dateOptFlyBy*24*3600;
    [YL, MoL, DL] = ymd(datetime(FlyByDateLoop2));
    [HL, ML, SL]  = hms(datetime(FlyByDateLoop2));
    datemjd2000L = date2mjd2000([YL MoL DL HL ML SL]);
                
    [kepFlyBy, ksun] = uplanet(datemjd2000L, flyByPlanID);
    [rrFlyBy, vvFlyBy] = par2car(kepFlyBy, ksun);
    
    TOF1 = (dateOptFlyBy*24*3600) - secFromDep2(i);
    [~, ~, ~, err, vI_1, vF_1, ~, ~] = lambertMR(rrDep, rrFlyBy, TOF1, ksun, 0, 0, 0, 1);
    if err == 0
        DV1_1 = norm(vI_1' - vvDep);
        DV2_1 = norm(vF_1' - vvFlyBy);
        DVfirstLeg2(i) = DV1_1;
    end
    
end

for i = 1:floor(dateOptFlyBy*24*3600/synPeriod1)+1
    indexes = find(secFromDep2 >= (i-1)*synPeriod1 + 1 & secFromDep2 < i*synPeriod1);
    [minDV2(i), ind] = min(DVfirstLeg2(indexes));
    ciao = secFromDep2(indexes);
    minDVdays2(i) = ciao(ind)/(24*3600);
end

maxMinDV = max(minDV2);
minMinDV = min(minDV2);

f = @(DV) (maxMinDV - DV).^4/(maxMinDV - minMinDV)^4;

figure
plot(minDVdays2, minDV2); grid on


dateOptDep = (sum(minDVdays2(not(isnan(minDV2))).*f(minDV2(not(isnan(minDV2))))))/(sum(f(minDV2(not(isnan(minDV2))))));


figure
plot(minDVdays2, f(minDV2)); hold on; grid on
xline(dateOptDep)


dateOptDep = dateOptDep;
date = depDate;
date(6) = date(6) + dateOptDep*24*3600;
[YLLL, MoLLL, DLLL] = ymd(datetime(date));
[HLLL, MLLL, SLLL]  = hms(datetime(date));
dateOptDepString = datestr([YLLL MoLLL DLLL HLLL MLLL SLLL])

%%%
arrArrVec2 = linspace(dateOptFlyBy*24*3600, TOTdays*24*3600, nTOF2);
DVsecondLeg3 = nan(nDep, 1);
for i = 1:nTOF2
    depDateLoop3 = depDate;
    depDateLoop3(6) = depDateLoop3(6) + arrArrVec2(1);
    
    [Yl, Mol, Dl] = ymd(datetime(depDateLoop3));
    [Hl, Ml, Sl]  = hms(datetime(depDateLoop3));
    
    datemjd2000l = date2mjd2000([Yl Mol Dl Hl Ml Sl]);
    
    % FlyBy planet
    [kepFlyBy, ksun] = uplanet(datemjd2000l, flyByPlanID);
    [rrFlyBy, vvFlyBy] = par2car(kepFlyBy, ksun);
    
    % Arr planet
    ArrDateLoop3 = depDate;
    ArrDateLoop3(6) = ArrDateLoop3(6) + arrArrVec2(i);
    [YL, MoL, DL] = ymd(datetime(ArrDateLoop3));
    [HL, ML, SL]  = hms(datetime(ArrDateLoop3));
    datemjd2000L = date2mjd2000([YL MoL DL HL ML SL]);
                
    [kepArr, ksun] = uplanet(datemjd2000L, arrPlanID);
    [rrArr, vvArr] = par2car(kepArr, ksun);
    
    TOF2 = arrArrVec2(i) - dateOptFlyBy*24*3600;
    [~, ~, ~, err, vI_2, vF_2, ~, ~] = lambertMR(rrFlyBy, rrArr, TOF2, ksun, 0, 0, 0, 1);
    if err == 0
        DV1_2 = norm(vI_2' - vvFlyBy);
        DV2_2 = norm(vF_2' - vvArr);
        DVsecondLeg3(i) = DV2_2;
    end
    
end

for i = 1:floor((TOTdays - dateOptFlyBy)*24*3600/synPeriod2)+1
    indexes = find(arrArrVec2 >= (i-1)*synPeriod2 + dateOptFlyBy*24*3600 + 1 & arrArrVec2 < i*synPeriod2 + dateOptFlyBy*24*3600);
    [minDV3(i), ind] = min(DVsecondLeg3(indexes));
    ciao = arrArrVec2(indexes);
    minDVdays3(i) = ciao(ind)/(24*3600);
end

maxMinDV = max(minDV3);
minMinDV = min(minDV3);

f = @(DV) (maxMinDV - DV).^2/(maxMinDV - minMinDV)^2;

figure
plot(minDVdays3, minDV3); grid on


dateOptArr = (sum(minDVdays3(not(isnan(minDV3))).*f(minDV3(not(isnan(minDV3))))))/(sum(f(minDV3(not(isnan(minDV3))))));


figure
plot(minDVdays3, f(minDV3)); hold on; grid on
xline(dateOptArr)


dateOptArr = dateOptArr;
date = depDate;
date(6) = date(6) + dateOptArr*24*3600;
[YLL, MoLL, DLL] = ymd(datetime(date));
[HLL, MLL, SLL]  = hms(datetime(date));
dateOptArrString = datestr([YLL MoLL DLL HLL MLL SLL])


%%
data.timeWindows.depDays = secFromDep/(24*3600);
data.timeWindows.DVfirstLeg = DVfirstLeg;
data.timeWindows.FlyByDays = arrFlyByVec/(24*3600);
data.timeWindows.DVsecondLeg = DVsecondLeg;
data.timeWindows.DVsecondLeg2 = DVsecondLeg2;
data.timeWindows.arrDays = arrArrVec/(24*3600);
data.timeWindows.depDate = [YLLL MoLLL DLLL HLLL MLLL SLLL];
data.timeWindows.arrDate = [YLL MoLL DLL HLL MLL SLL];
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
daysSyn2 = synPeriod2/(24*3600);

if settings.timeWindows.plot
    
    figure
    [A, B] = meshgrid(data.timeWindows.depDays, data.timeWindows.FlyByDays);
    contour(A, B, data.timeWindows.DVfirstLeg', 100)
    xticks(vec2); xticklabels(vec1); xtickangle(-45)
    yticks(vec4); yticklabels(vec3)
    xlabel('Mars Departure Date'); ylabel('Earth Arrival Date')
    title('$\Delta$V First Leg Porkchop')
    
    figure
    [A, B] = meshgrid(data.timeWindows.FlyByDays, data.timeWindows.arrDays);
    contour(A, B, data.timeWindows.DVsecondLeg2', 100)
    xticks(vec2); xticklabels(vec1); xtickangle(-45)
    yticks(vec4); yticklabels(vec3)
    xlabel('Earth Departure Date'); ylabel('Venus Arrival Date')
    title('$\Delta$V Second Leg Porkchop')
    
    figure
    minima1 = min(data.timeWindows.DVfirstLeg(1:end-5,:), [], 2);
    plot(data.timeWindows.depDays(1:end-5), minima1, '-o', 'MarkerSize', 5, 'LineWidth', 1)
    hold on; grid on;
    for i = 1:floor(TOTdays/daysSyn1)
        xline(i*daysSyn1, '--k');
    end
    xticks(vec2); xticklabels(vec1); xtickangle(-45)
    xlabel('Mars Departure Date'); ylabel('min $\Delta$V [$\frac{km}{s}$]')
    title('min $\Delta$V First Leg')
    
    figure
    minima2 = min(data.timeWindows.DVsecondLeg2(1:end-5,:), [], 2);
    plot(data.timeWindows.FlyByDays(1:end-5), minima2, '-o', 'MarkerSize', 5, 'LineWidth', 1)
    hold on; grid on;
    for i = 1:floor(TOTdays/daysSyn2)
        xline(i*daysSyn2, '--k');
    end
    xticks(vec4); xticklabels(vec3); xtickangle(-45)
    xlabel('Earth Departure Date'); ylabel('min $\Delta$V [$\frac{km}{s}$]')
    title('min $\Delta$V Second Leg')
    
    figure
    minima3 = min(data.timeWindows.DVsecondLeg(1:end-5,:), [], [2 3]);
    plot(data.timeWindows.depDays(1:end-5), minima3, '-o', 'MarkerSize', 5, 'LineWidth', 1)
    hold on; grid on;
    for i = 1:floor(TOTdays/daysSyn2)
        xline(i*daysSyn2, '--k');
    end
    xticks(vec2); xticklabels(vec1); xtickangle(-45)
    xlabel('Mars Departure Date'); ylabel('min $\Delta$V [$\frac{km}{s}$]')
    title('min $\Delta$V Second Leg')
   
    
    
    
end


end