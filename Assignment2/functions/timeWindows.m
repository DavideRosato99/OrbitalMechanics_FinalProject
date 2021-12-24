function [DVfirstLeg, DVsecondLeg, DVTOT, daysSpan, TOF1span, TOF2span,...
    DVTOTreal, dtTOTreal] = timeWindows(depID, flybyID, arrID, depDate, arrDate, deltaD)

TOTdaysWindow = datenum(arrDate) - datenum(depDate);

daysSpan = 0 : deltaD : TOTdaysWindow;

N = length(daysSpan);

% First leg time span
TOF1span = 0 : 3*365 : TOTdaysWindow;
Ntof1 = length(TOF1span);
% Secon leg time span
TOF2span = 0 : 3*365 : TOTdaysWindow;
Ntof2 = length(TOF2span);

% Initialize
DVfirstLeg  = NaN(N, Ntof1);
DVsecondLeg = NaN(N, Ntof1, Ntof2);
DVTOT       = NaN(N, Ntof1, Ntof2);
dtTOT       = NaN(N, Ntof1, Ntof2);
dtTOTreal   = [];
date = depDate;
for i = 1:N
    
    d = daysSpan(i);
    date(3) = depDate(3) + d;
    [Yl, Mol, Dl] = ymd(datetime(date));
    [Hl, Ml, Sl] = hms(datetime(date));
    
    datemjd2000 = date2mjd2000([Yl Mol Dl Hl Ml Sl]);
    
    % Departure planet
    [kepDep, ksun] = uplanet(datemjd2000, depID);
    [rrDep, vvDep] = par2car(kepDep, ksun);
    
    dateL = date;
    %%% TOFs loops
    for j = 1:Ntof1
        dateL(3) = date(3) + TOF1span(j);
        [YL, MoL, DL] = ymd(datetime(dateL));
        [HL, ML, SL] = hms(datetime(dateL));
        datemjd2000L = date2mjd2000([YL MoL DL HL ML SL]);
        
        deltaTime = datenum(arrDate) - datenum(dateL);
        if deltaTime >= 0
            % Flyby planet
            [kepFlyby, ~] = uplanet(datemjd2000L, flybyID);
            [rrFlyby, vvFlyby] = par2car(kepFlyby, ksun);

            [~, ~, ~, err, vI, vF, ~, ~] = lambertMR(rrDep, rrFlyby, TOF1span(j)*24*60*60, ksun, 0, 0, 0, 1);

            if err ~= 0
                DVfirstLeg(i, j) = NaN;
                DVTOT(i, j, :) = NaN;
                dtTOT(i, j, :) = NaN;
            else
                DV1 = norm(vI' - vvDep);
                DV2 = norm(vF' - vvFlyby);
                DVfirstLeg(i, j) = DV1 + DV2;

                for k = 1:Ntof2
                    dateLL(3) = dateL(3) + TOF2span(k);
                    [YLL, MoLL, DLL] = ymd(datetime(dateLL));
                    [HLL, MLL, SLL] = hms(datetime(dateLL));
                    datemjd2000LL = date2mjd2000([YLL MoLL DLL HLL MLL SLL]);
                    
                    deltaTime = datenum(arrDate) - datenum(dateLL);
                    if deltaTime >= 0
                        % Arrival planet
                        [kepArr, ~] = uplanet(datemjd2000LL, arrID);
                        [rrArr, vvArr] = par2car(kepArr, ksun);
                        [~, ~, ~, err, vI, vF, ~, ~] = lambertMR(rrFlyby, rrArr, TOF2span(k)*24*60*60, ksun, 0, 0, 0, 1);

                        if err ~= 0
                            DVsecondLeg(i, j, k) = NaN;
                        else
                            DV1 = norm(vI' - vvFlyby);
                            DV2 = norm(vF' - vvArr);
                            DVsecondLeg(i, j, k) = DV1 + DV2;
                            DVTOT(i, j, k) = DVfirstLeg(i, j) + DVsecondLeg(i, j, k);
                            dtTOT(i, j, k) = daysSpan(i) + TOF1span(j) + TOF2span(k);
                            if not(any(dtTOTreal == dtTOT(i, j, k)))
                                dtTOTreal(end+1) = dtTOT(i, j, k);
                            end
                        end
                    else
                        DVsecondLeg(i, j, k) = NaN;
                        DVTOT(i, j, k) = NaN;
                        dtTOT(i, j, k) = NaN;
                    end
                end
            end
        else
            DVfirstLeg(i, j) = NaN;
        end
    end
    
    
    
end

%% PARSING DATA
%%% PARSE dtTOT
dtTOTreal = sort(dtTOTreal);
DVTOTreal = nan(N, length(dtTOTreal));
for i = 1:N
    for j = 1:Ntof1
        for k = 1:Ntof2
            if not(isnan(dtTOT(i, j, k)))
                index = find(dtTOTreal == dtTOT(i, j, k));
                [a, b] = find(squeeze(dtTOT(i, :, :)) == dtTOTreal(index));
                DVvec = [];
                
                for ii = 1:length(a)
                    DVvec = [DVvec, DVTOT(i, a(ii), b(ii))];
                end

                DVTOTreal(i, index) = min(DVvec);
                
            end
        end
    end
end

a = 1;








