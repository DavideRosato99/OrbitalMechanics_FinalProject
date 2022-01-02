function [data] = gridSearch(data, settings)

%%
minDate = data.timeWindows.depDate;
maxDate = data.timeWindows.arrDate;
depPlanID   = data.starting.depPlanID;
flyByPlanID = data.starting.flyByPlanID;
arrPlanID   = data.starting.arrPlanID;
minHfl = data.optimization.minHfl;
TOF1max = data.timeWindows.maxTOF1days;
TOF2max = data.timeWindows.maxTOF2days;
nDep = data.gridSearch.nDep;
nTOF1 = data.gridSearch.nTOF1;
nTOF2 = data.gridSearch.nTOF2;

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
TOF1End      = TOF1max;
TOF2End      = TOF2max;

ctr = 1;

t = tic;
while error >= 0.001 && not(flag)
    
    DVTOT = nan(nDep, nTOF1, nTOF2);
    
    vec1 = linspace(departureStart, departureEnd, nDep);
    vec2 = linspace(TOF1Start, TOF1End, nTOF1);
    vec3 = linspace(TOF2Start, TOF2End, nTOF2);
    
    fid = fopen('parforbest.txt', 'w');
    fprintf(fid, '%5.15f %d %d %d\n', 10000000000, 0, 0, 0);
    fclose(fid);
    
    fprintf(strcat("Calculate ", num2str(ctr), "th generation:            "));
    parforProgress(nDep);
    parfor i = 1:nDep
        if minDatemjd2000 + vec1(i) <= maxDatemjd2000
            for j = 1:nTOF1
                if minDatemjd2000 + vec1(i) + vec2(j) <= maxDatemjd2000
                    for k = 1:nTOF2
                        if minDatemjd2000 + vec1(i) + vec2(j) + vec3(k) <= maxDatemjd2000
                            [kepDep, ksun] = uplanet(minDatemjd2000 + vec1(i), depPlanID);
                            [rrDep, vvDep] = par2car(kepDep, ksun);
                            [kepFlyBy, ksun] = uplanet(minDatemjd2000 + vec2(j) + vec1(i), flyByPlanID);
                            [rrFlyBy, vvFlyBy] = par2car(kepFlyBy, ksun);

                            [kepArr, ksun] = uplanet(minDatemjd2000 + vec3(k) + vec2(j) + vec1(i), arrPlanID);
                            [rrArr, vvArr] = par2car(kepArr, ksun);
                            
                            [~, ~, ~, error1, ~, ~, ~, ~, ~, ~, ~, ...
                            error2, ~, ~, ~, ~, DV1, DV2, errorFB, ~, ~, ~, ...
                            delta_V_powFB, ~, ~, ~, ~] = deltaVtot(...
                            rrDep, rrFlyBy, rrArr, vvDep, vvFlyBy, vvArr, vec2(j)*24*3600,...
                            vec3(k)*24*3600, minHfl, flyByPlanID);
                           
                            if error1 == 0 && error2 == 0 && errorFB == 0
                                DVTOT(i, j, k) = DV1 + DV2 + delta_V_powFB;
                                
                                fid = fopen('parforbest.txt', 'r');
                                best = fscanf(fid, '%f %d %d %d');
                                while isempty(best)
                                    best = fscanf(fid, '%f %d %d %d');
                                end
                                fclose(fid);
                                
                                if DVTOT(i, j, k) <= best(1)
                                    fid = fopen('parforbest.txt', 'w');
                                    fprintf(fid, '%5.15f %d %d %d\n', DVTOT(i, j, k), i, j, k);
                                    fclose(fid);
                                end
                            end
                        end
                    end
                end
            end
        end
        fprintf(repmat('\b', 1, 5));
        parforProgress
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
    
    fid = fopen('parforbest.txt', 'r');
    best = fscanf(fid, '%f %d %d %d');
    fclose(fid);
    delete('parforbest.txt');
   
    error = abs(bestGlobal - best(1))
    best(1)
    
    if best(1) < bestGlobal
        
        bestGlobal = best(1);
        bestGlobalI = best(2);
        bestGlobalJ = best(3);
        bestGlobalK = best(4);
        
        deltaDep  = departureEnd - departureStart;
        deltaTOF1 = TOF1End - TOF1Start;
        deltaTOF2 = TOF2End - TOF2Start;
        
        if not( vec1(bestGlobalI) - deltaDep/4 < departureStart )
            departureStart = vec1(bestGlobalI) - deltaDep/4
        end
        if not ( vec1(bestGlobalI) + deltaDep/4 > departureEnd )
            departureEnd   = vec1(bestGlobalI) + deltaDep/4
        end
        
        if not( vec2(bestGlobalJ) - deltaTOF1/4 < TOF1Start )
            TOF1Start = vec2(bestGlobalJ) - deltaTOF1/4
        end
        if not ( vec2(bestGlobalJ) + deltaTOF1/4 > TOF1End )
            TOF1End   = vec2(bestGlobalJ) + deltaTOF1/4
        end
        
        if not( vec3(bestGlobalK) - deltaTOF2/4 < TOF2Start  )
            TOF2Start = vec3(bestGlobalK) - deltaTOF2/4
        end
        if not ( vec3(bestGlobalK) + deltaTOF2/4 > TOF2End )
            TOF2End   = vec3(bestGlobalK) + deltaTOF2/4
        end
        
    else
        
        flag = true;
        ctr = ctr-1;
        
    end
    
end

comptTime = toc(t);

date = minDate;
date(6) = date(6) + vec1(bestGlobalI)*24*3600;
[Y, Mo, D] = ymd(datetime(date));
[H, M, S] = hms(datetime(date));
fprintf(strcat("Optimal Departure Date:  ", datestr([Y Mo D H M S]), "\n"))

date = minDate;
date(6) = date(6) + vec1(bestGlobalI)*24*3600 + vec2(bestGlobalJ)*24*3600;
[Y, Mo, D] = ymd(datetime(date));
[H, M, S] = hms(datetime(date));
fprintf(strcat("Optimal Fly-By Date:  ", datestr([Y Mo D H M S]), "\n"))

date = minDate;
date(6) = date(6) + vec1(bestGlobalI)*24*3600 + vec2(bestGlobalJ)*24*3600 + vec3(bestGlobalK)*24*3600;
[Y, Mo, D] = ymd(datetime(date));
[H, M, S] = hms(datetime(date));
fprintf(strcat("Optimal Arrival Date:  ", datestr([Y Mo D H M S]), "\n"))

%% PLOT
if settings.gridSearch.plot
    %%%
    figure
    maxVal = 0;
    minVal = 10e9;
    
    for i = 1:ctr-1
        maxI = max(max(max(DVTOTglobal{i})));
        minI = min(min(min(DVTOTglobal{i})));
        if maxI >= maxVal
            maxVal = maxI;
        end
        if minI <=  minVal
            minVal = minI;
        end
    end
    
    for i = 1:ctr-1
        x = linspace(departureStartGlobal{i}, departureEndGlobal{i}, nDep);
        y = linspace(TOF1StartGlobal{i}, TOF1EndGlobal{i}, nTOF1);
        z = linspace(TOF2StartGlobal{i}, TOF2EndGlobal{i}, nTOF2);
        colormap = parula(300);
        for ii = 1:length(x)
            for jj = 1:length(y)
                for kk = 1:length(z)
                    if not(isnan(DVTOTglobal{i}(ii, jj, kk)))
                        plot3(x(ii), y(jj), z(kk), 'o', 'color',...
                            colormap(floor((DVTOTglobal{i}(ii, jj, kk) - minVal)/(maxVal - minVal)*299) + 1, :),...
                            'MarkerFaceColor', colormap(floor((DVTOTglobal{i}(ii, jj, kk) - minVal)/(maxVal - minVal)*299) + 1, :),...
                            'MarkerSize', 0.1);
                        hold on;
                    end
                end
            end
        end
        
        patch([x(1) x(1) x(1) x(1)] , [y(1) y(end) y(end) y(1)], [z(1) z(1) z(end) z(end)], [z(1) z(1) z(end) z(end)], 'FaceColor', [1, 1, 1], 'FaceAlpha', 0, 'LineWidth', 0.1, 'LineStyle', '-')
        patch([x(end) x(end) x(end) x(end)] , [y(1) y(end) y(end) y(1)], [z(1) z(1) z(end) z(end)], [z(1) z(1) z(end) z(end)], 'FaceColor', [1, 1, 1], 'FaceAlpha', 0, 'LineWidth', 0.1, 'LineStyle', '-')
        patch([x(1) x(end) x(end) x(1)] , [y(1) y(1) y(1) y(1)], [z(1) z(1) z(end) z(end)], [z(1) z(1) z(end) z(end)], 'FaceColor', [1, 1, 1], 'FaceAlpha', 0, 'LineWidth', 0.1, 'LineStyle', '-')
        patch([x(1) x(end) x(end) x(1)] , [y(end) y(end) y(end) y(end)], [z(1) z(1) z(end) z(end)], [z(1) z(1) z(end) z(end)], 'FaceColor', [1, 1, 1], 'FaceAlpha', 0, 'LineWidth', 0.1, 'LineStyle', '-')
        patch([x(1) x(end) x(end) x(1)] , [y(1) y(1) y(end) y(end)], [z(1) z(1) z(1) z(1)], [z(1) z(1) z(1) z(1)], 'FaceColor', [1, 1, 1], 'FaceAlpha', 0, 'LineWidth', 0.1, 'LineStyle', '-')
        patch([x(1) x(end) x(end) x(1)] , [y(1) y(1) y(end) y(end)], [z(end) z(end) z(end) z(end)], [z(end) z(end) z(end) z(end)], 'FaceColor', [1, 1, 1], 'FaceAlpha', 0, 'LineWidth', 0.1, 'LineStyle', '-')
    end
    xlabel('Time from first date [days]');
    ylabel('TOF 1 [days]');
    zlabel('TOF 2 [days]');
    title('Evaluations - Grid Search')
    axis equal; grid on;
        
    for i = 1:ctr-1
        minMin = min(min(min(DVTOTglobal{i})));
        x = linspace(departureStartGlobal{i}, departureEndGlobal{i}, nDep);
        y = linspace(TOF1StartGlobal{i}, TOF1EndGlobal{i}, nTOF1);
        z = linspace(TOF2StartGlobal{i}, TOF2EndGlobal{i}, nTOF2);
        for ii = 1:length(x)
            for jj = 1:length(y)
                for kk = 1:length(z)
                    if not(isnan(DVTOTglobal{i}(ii, jj, kk)))
                        if DVTOTglobal{i}(ii, jj, kk) == minMin
                            plot3(x(ii), y(jj), z(kk), 'o', 'color', 'r',...
                                'MarkerFaceColor', 'r', ...
                                'MarkerSize', 6);
                        end
                        hold on;
                    end
                end
            end
        end
    end
    
    %%%
    figure
    for i = 1:ctr-1
        minMin = min(min(min(DVTOTglobal{i})));
        x = linspace(departureStartGlobal{i}, departureEndGlobal{i}, nDep);
        y = linspace(TOF1StartGlobal{i}, TOF1EndGlobal{i}, nTOF1);
        z = linspace(TOF2StartGlobal{i}, TOF2EndGlobal{i}, nTOF2);
        for ii = 1:length(x)
            for jj = 1:length(y)
                for kk = 1:length(z)
                    if not(isnan(DVTOTglobal{i}(ii, jj, kk)))
                        if DVTOTglobal{i}(ii, jj, kk) == minMin
                            minMinima{i}= [x(ii), y(jj), z(kk)];
                        end
                    end
                end
            end
        end
    end
    colormap = parula(ctr-1);
    plot3(minMinima{1}(1), minMinima{1}(2), minMinima{1}(3), 'o', 'color', colormap(1, :),...
        'MarkerFaceColor', colormap(1, :), ...
        'MarkerSize', 6); hold on; grid on;
    for i = 2:ctr-1
        plot3(minMinima{i}(1), minMinima{i}(2), minMinima{i}(3), 'o', 'color', colormap(i, :),...
            'MarkerFaceColor', colormap(i, :), ...
            'MarkerSize', 6);
        U = minMinima{i}(1) - minMinima{i-1}(1);
        V = minMinima{i}(2) - minMinima{i-1}(2);
        W = minMinima{i}(3) - minMinima{i-1}(3);
        quiver3(minMinima{i-1}(1), minMinima{i-1}(2), minMinima{i-1}(3), ...
            U, V, W, 'k');
    end
    xlabel('Time from first date [days]');
    ylabel('TOF 1 [days]');
    zlabel('TOF 2 [days]');
    title('Best individuals - Grid Search')
    axis equal; grid on;
    
    
    
    %%%
    figure
    meanVal = zeros(ctr-1, 1);
    minVal  = zeros(ctr-1, 1);
    stdDevUp  = zeros(ctr-1, 1);
    stdDevDown  = zeros(ctr-1, 1);
    for i = 1:ctr-1
        DVTOTparsed = DVTOTglobal{i}(not(isnan(DVTOTglobal{i}(:, :, :))));
        
        meanVal(i) = mean(DVTOTparsed);
        minVal(i)  = min(DVTOTparsed);
        
        diff = DVTOTparsed - meanVal(i);
        stdDevUp(i) = sqrt(sum(diff(diff >= 0).^2) / length(diff(diff >= 0)));
        stdDevDown(i) = sqrt(sum(diff(diff < 0).^2) / length(diff(diff < 0)));
    end
    
    subplot(2,1,1)
    plot([1 : ctr-1], minVal, '-o', 'LineWidth', 1); grid on;
    ylabel('$\Delta$V [$\frac{km}{s}$]'); title('Output - Grid Search')
    subplot(2,1,2)
    errorbar([1 : ctr-1], meanVal, stdDevDown, stdDevUp, 'LineWidth', 1); grid on;
    xlim([0.5 ctr-0.5]);
    legend('Mean and std Dev.', 'Best individual');
    xlabel('Generation [-]'); ylabel('$\Delta$V [$\frac{km}{s}$]')
    
    
end






