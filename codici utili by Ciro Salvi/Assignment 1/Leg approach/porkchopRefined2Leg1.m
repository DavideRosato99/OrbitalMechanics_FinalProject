% porkchopRefined2Leg1 allows a sampling of the minima along the axis of departure dates of the
% first leg. The synodic period between Mercury and Venus is expoited to shift the rectangular window
% in which the minimum is found.
% The frequency analysis of the signal of minima is then performed.
% 
% PROTOTYPE: 
%   porkchopRefined2Leg1
%
% CONTRIBUTORS:
%   Fabio Spada
%

clc; clear all; close all;

checkPorkchop = 1;
checkPorkchopRefine = 1;

%% Mercury-Venus _ First impulse

idDep = 1;
idArr = 2;

% suggested step: 4 months per departure date, 8 months per arrival date

if (checkPorkchopRefine)                    % precisely plotting the porkchop may be computationally expensive for a large window:
    t_start_j_1 = [2027 6 1 0 0 0];             % if the aim is obtaining good looking porkchop plots, the window is set as a small one,
    t_end_j_1 = [2029 6 1 0 0 0];               % otherwise a large window can be chosen for gathering data on a wider window
    t_start_j_2 = [2027 6 1 0 0 0];
    t_end_j_2 = [2029 6 1 0 0 0];
else
    t_start_j_1 = [2027 6 1 0 0 0];             % if the aim is obtaining good looking porkchop plots, the window is set as a small one,
    t_end_j_1 = [2107 6 1 0 0 0];               % otherwise a large window can be chosen for gathering data on a wider window
    t_start_j_2 = [2027 6 1 0 0 0];
    t_end_j_2 = [2107 6 1 0 0 0];
end

% % Mercury Venus Synodic period = 87.97*224.7/(224.7 - 87.97) = 4.819 30days-lasting months
step1 = 4.819 * 30;                 % number of months * ndays per month
step2 = 10 * 30;                    % number of months * ndays per month
if (checkPorkchopRefine)            % better refinement for a better looking plot
    s_size_1 = ceil( 4*step1 );     % number of days * 4 ---> 6 h step
    s_size_2 = ceil( 4*step2 );     % number of days * 4 ---> 6 h step
else
    s_size_1 = ceil( 2*step1 );     % number of days * 2 ---> 12 h step
    s_size_2 = ceil( 2*step2 );     % number of days * 2 ---> 12 h step
end

t_start_j_1 = date2mjd2000(t_start_j_1);
t_end_j_1 = date2mjd2000(t_end_j_1);
t_start_j_2 = date2mjd2000(t_start_j_2);
t_end_j_2 = date2mjd2000(t_end_j_2);


% The aim is building up a vector of optimal solutions in terms of first deltav of the
% transfer Mercury-Venus. This can be done analyzing the blocks of the porkchop plot so to
% get the minimum per each nucleus

tDepmjd = t_start_j_1;
tArrmjd = t_start_j_2;
tDepMinCost = [];
tArrMinCost = [];
ToFMinCost = [];
MinCost = [];
ii = 0;

tDepEndmjd = tDepmjd + step1;
tArrEndmjd = tArrmjd + step2;

while tArrEndmjd < t_end_j_2
    
    ii = ii+1;
    
    [ ~, T_1, T_2, DV1, ~, ~, ~, ~, ~, ~, ~ ] = deltaV1( idDep, idArr, tDepmjd , tDepEndmjd, tArrmjd, tArrEndmjd , s_size_1, s_size_2);
    
    tPlot.Dep = T_1;
    tPlot.Arr = T_2;
    
    DV1(DV1==0) = NaN;
    
    minDV1 = min(min(DV1));
    [indDep, indArr] = find(DV1 == minDV1);
    
    tDepMinCost(ii) = T_1(indDep);
    tArrMinCost(ii) = T_2(indArr);
    MinCost(ii) = minDV1;
    ToFMinCost(ii) = [tArrMinCost(ii)-tDepMinCost(ii)];
    
    tPlot.tArrMinCost = T_2(indArr);
    tPlot.tDepMinCost = T_1(indDep);
    
    tDepmjd = tDepEndmjd;
    tArrmjd = tDepEndmjd;
    
    tDepEndmjd = tDepmjd + step1;
    tArrEndmjd = tArrmjd + step2;
    
    if checkPorkchop
        PorkchopPlot1(tPlot, DV1', min(min(DV1)), 20, 'Dates');
    end
    
end

%% Cost Plot with arrival date

for bb = 1 : length( tArrMinCost )
    tArrMinCost_g(bb) = datenum(mjd20002date(tArrMinCost(bb)));
end

figure; hold on; grid minor;
plot( tArrMinCost_g, MinCost, 'o--', 'LineWidth', 1, 'MarkerSize', 5 );
title('Venus date vs Minimum Costs');
xlabel('Venus date');
ylabel('[km/s]');
title( 'minimum \Deltav of the leg 1 ');
datetick('x', 'yyyy mmm dd','keeplimits');
xtickangle(-45)

dcm = datacursormode;
dcm.Enable = 'on';
dcm.SnapToDataVertex = 'on';
dcm.DisplayStyle = 'datatip';
dcm.UpdateFcn = @displayCoordinates_j;

%% Cost Plot with departure date

for bb = 1 : length( tDepMinCost )
    tDepMinCost_g(bb) = datenum(mjd20002date(tDepMinCost(bb)));
end

figure; hold on; grid minor;
plot( tDepMinCost_g, MinCost, 'o--', 'LineWidth', 1, 'MarkerSize', 5 );
title('Mercury date vs Minimum Costs');
xlabel('Mercury date');
ylabel('[km/s]');
title( 'minimum \Deltav of the leg 1 ');
datetick('x', 'yyyy mmm dd','keeplimits');
xtickangle(-45)

dcm = datacursormode;
dcm.Enable = 'on';
dcm.SnapToDataVertex = 'on';
dcm.DisplayStyle = 'datatip';
dcm.UpdateFcn = @displayCoordinates_j;

%% ToF plot with departure date

figure; hold on; grid minor;
plot( tDepMinCost_g, ToFMinCost, 'o--', 'LineWidth', 1, 'MarkerSize', 5 );
title('mercury date vs ToF');
xlabel('Mercury departure date');
ylabel('[Time]');
title( 'Minimum-cost-transfer related ToF ');
datetick('x', 'yyyy mmm dd','keeplimits');
xtickangle(-45)

dcm = datacursormode;
dcm.Enable = 'on';
dcm.SnapToDataVertex = 'on';
dcm.DisplayStyle = 'datatip';
dcm.UpdateFcn = @displayCoordinates_j;

%% Frequency Analysis / Arrival data / pchip

tArrMinCostSec = tArrMinCost*24*3600; % swap to seconds

% OPTIMAL BEHAVIOUR
[DFTArrPchip, ArrSig] = generalSignalDFT(tArrMinCostSec, MinCost, 'pchip' );
title('Minimum cost vs Arrival date approx DFT - pchip');

% Verify interpolation correctness
figure
plot(tArrMinCostSec, MinCost,'o', ArrSig(:, 1), ArrSig(:, 2));
legend('Real signal','Interpolated signal');
grid on

%% Frequency Analysis / Departure data

tDepMinCostSec = tDepMinCost*24*3600; % swap to seconds

% OPTIMAL BEHAVIOUR
[DFTDepPchip, DepSig] = generalSignalDFT(tDepMinCostSec, MinCost, 'pchip' );
title('Minimim cost vs departure date approx DFT - pchip');

% Verify interpolation correctness
figure
plot(tDepMinCostSec, MinCost,'o', DepSig(:, 1), DepSig(:, 2));
legend('Real signal','Interpolated signal');
grid minor

%% Overlap variables definition

% Minima related vaiables definition
ILeg.tDepMinCost = tDepMinCost;
ILeg.tArrMinCost = tArrMinCost;
ILeg.MinCost = MinCost;
ILeg.ToFMinCost = ToFMinCost;


%% Callback function
function txt = displayCoordinates_j(~,info)
x = info.Position(1);
y = info.Position(2);
txt = {['Date = ', datestr(datetime(mjd20002date(x)), ' dd mmm yyyy HH:MM:SS')], ['Delta V = ' num2str(y) ' [km/s]']};
end
