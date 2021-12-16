% porkchopRefined2Leg2 allows a sampling of the minima along the axis of departure dates of the
% second leg. The synodic period between Venus and Jupiter is expoited to shift the rectangular window
% in which the minimum is found.
% The frequency analysis of the signal of minima is then performed.
%
% PROTOTYPE: 
%   porkchopRefined2Leg2
%
% CONTRIBUTORS:
%   Fabio Spada
%

clc; clear all; close all;
checkPorkchop = 0;

%% Venus-Jupiter _ Second impulse

idDep = 2;
idArr = 5;

% 1) the departure window was cut at the beginning, where an early departure from venus could
% determine too high costs for the first transfer;
% 2) the arrival window was postponed of 1.5 years wrt the departure window, as 576 days were 
% required to get reasonable costs (from 10km/s below)
t_start_j_1 = [2027 10 1 0 0 0];
t_end_j_1 = [2065 12 1 0 0 0];

% to-be-improven control using dates
t_start_j_2 = [2029 4 1 0 0 0];
t_end_j_2 = [2067 6 1 0 0 0];

% % Venus Jupiter Synodic period = 4332.59*224.7/(4332.59-224.7) = 236.9910 days = 7.8997
% 30days-lasting months
step1 = 7.8997 * 30;              % number of months * ndays per month
step2 = 30 * 30;                % number of months * ndays per month
s_size_1 = ceil(2 * step1);     % number of days * 2 ---> 12 h step
s_size_2 = ceil(2 * step2);     % number of days * 2 ---> 12 h step 

t_start_j_1 = date2mjd2000(t_start_j_1);
t_end_j_1 = date2mjd2000(t_end_j_1);
t_start_j_2 = date2mjd2000(t_start_j_2);
t_end_j_2 = date2mjd2000(t_end_j_2);


% The aim is building up a vector of optimal solutions in terms of second deltav of the
% transfer Venus-Jupiter. This can be done analyzing the blocks of the porkchop plot so to
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

    [ ~, T_1, T_2, ~, DV2, ~, ~, ~, ~, ~, ~ ] = deltaV1( idDep, idArr, tDepmjd , tDepEndmjd, tArrmjd, tArrEndmjd , s_size_1, s_size_2);
    
    tPlot.Dep = T_1;
    tPlot.Arr = T_2;
    
    DV2(DV2==0) = NaN;
    
    minDV2 = min(min(DV2));
    [indDep, indArr] = find(DV2 == minDV2);
    
    tDepMinCost(ii) = T_1(indDep);
    tArrMinCost(ii) = T_2(indArr);
    MinCost(ii) = minDV2;
    ToFMinCost(ii) = [tArrMinCost(ii) - tDepMinCost(ii)];
    
    tPlot.tArrMinCost = T_2(indArr);
    tPlot.tDepMinCost = T_1(indDep);
    
    tDepmjd = tDepEndmjd; 
    tArrmjd = tDepEndmjd + (t_start_j_2 - t_start_j_1);
    
    tDepEndmjd = tDepmjd + step1;
    tArrEndmjd = tArrmjd + step2;
    
    if checkPorkchop
        PorkchopPlot1(tPlot, DV2', min(min(DV2)), 5);
    end
end

%% Cost Plot with arrival date

for bb = 1 : length( tArrMinCost )
    tArrMinCost_g(bb) = datenum(mjd20002date(tArrMinCost(bb)));
end

    figure; hold on; grid minor;
    plot( tArrMinCost_g, MinCost, 'o--', 'LineWidth', 1, 'MarkerSize', 5, 'Color',[0.8500, 0.3250, 0.0980] ); 
    xlabel('Jupiter date');
    ylabel('[km/s]');
    title( 'minimum \Deltav of the leg 2 ');
    datetick('x', 'yyyy mmm dd','keeplimits');
    xtickangle(-45)
    
    ax = gca;
    ax.TickLabelInterpreter = 'Latex';
    
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
    plot( tDepMinCost_g, MinCost, 'o--', 'LineWidth', 1, 'MarkerSize', 5 , 'Color',[0.8500, 0.3250, 0.0980]); 
    xlabel('Venus date');
    ylabel('[km/s]');
    title( 'minimum \Deltav of the leg 2 ');
    datetick('x', 'yyyy mmm dd','keeplimits');
    xtickangle(-45)
    
    ax = gca;
    ax.TickLabelInterpreter = 'Latex';
    
    dcm = datacursormode;
    dcm.Enable = 'on';
    dcm.SnapToDataVertex = 'on';
    dcm.DisplayStyle = 'datatip';
    dcm.UpdateFcn = @displayCoordinates_j;

%% ToF plot with departure date

    figure; hold on; grid minor;
    plot( tDepMinCost_g, ToFMinCost, 'o--', 'LineWidth', 1, 'MarkerSize', 5, 'Color',[0.8500, 0.3250, 0.0980]); 
    xlabel('Venus departure date');
    ylabel('[km/s]');
    title( 'Minimum-cost-transfer related ToF ');
    datetick('x', 'yyyy mmm dd','keeplimits');
    xtickangle(-45)
    
    ax = gca;
    ax.TickLabelInterpreter = 'Latex';
    
    dcm = datacursormode;
    dcm.Enable = 'on';
    dcm.SnapToDataVertex = 'on';
    dcm.DisplayStyle = 'datatip';
    dcm.UpdateFcn = @displayCoordinates_j;
    

%% Frequency Analysis / Arrival data / pchip

tArrMinCostSec = tArrMinCost*24*3600; % swap to seconds

% build up interpolation vector 
tArrMinCostInterp = linspace(tArrMinCostSec(1), tArrMinCostSec(end), 100*length(tArrMinCostSec));

% build up equally spaced minimum vector (Necessary to be equally spaced for FFT)
MinCostInterpArr = interp1(tArrMinCostSec, MinCost, tArrMinCostInterp, 'pchip');

% Verify interpolation correctness
figure
plot(tArrMinCostSec, MinCost,'o', tArrMinCostInterp, MinCostInterpArr);
legend('Real signal','Interpolated signal');
xlabel('Arrival time');
title('Minimum Costs interpolation - Arrival Time');
grid on

%% Frequency Analysis / Departure data / pchip

tDepMinCostSec = tDepMinCost*24*3600; % swap to seconds

% OPTIMAL BEHAVIOUR
[DFTDepPchip, DepSig] = generalSignalDFT(tDepMinCostSec, MinCost, 'pchip' );
title('Minimim cost approx DFT - pchip');

% Verify interpolation correctness
figure
plot(tDepMinCostSec, MinCost,'o', DepSig(:, 1), DepSig(:, 2));
legend('Real signal','Interpolated signal');
xlabel('Departure time');
title('Minimum Costs interpolation - Departure Time');
grid on


%% Overlap variables definition

% Minima related vaiables definition
IILeg.tDepMinCost = tDepMinCost;
IILeg.tArrMinCost = tArrMinCost;
IILeg.MinCost = MinCost; 
IILeg.ToFMinCost = ToFMinCost;

%% Callback function
function txt = displayCoordinates_j(~,info)
x = info.Position(1);
y = info.Position(2);
txt = {['Date = ', datestr(datetime(mjd20002date(x)), ' dd mmm yyyy HH:MM:SS')], ['Delta V = ' num2str(y) ' [km/s]']};
end




