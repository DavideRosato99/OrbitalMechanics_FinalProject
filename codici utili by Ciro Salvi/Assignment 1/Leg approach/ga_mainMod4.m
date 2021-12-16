% ga_mainMod4 starting from a set of dates this script analyze the possible
% flyby option only shifting the dates on Venus
%
% PROTOTYPE: 
%   ga_mainMod4
%
% CONTRIBUTORS:
%   Alessandro Staffolani
%

if (isequal(exist('switchNew'),0) || switchNew == 0)
close all; clear all; clc
% load('box_matched_1_1.mat'); % all Nan results
% load('box_matched_5_5.mat');
% load('box_matched_8_2.mat');
load('box_matched_15_5.mat');
end

profile on

DATA = [];
ga_it = 15;
% good_ii = [ 1 7 21 29 44 46 52 63 ]; % 5_5
% good_ii = [ 4 35 ]; % 8_2
% good_ii = [ 1 3 9 23 26 28 34 50 52 58 71 74 76 ]; % 15_5
for ii = 1 : size(int_box,1)  % good_ii || 1 : size(int_box,1)
    ii
    t_merc_n  = int_box(ii,1); % [datenum]
    t_ven_1_n = int_box(ii,2); % [datenum]
    t_ven_2_n = int_box(ii,3); % [datenum]
    t_jup_n   = int_box(ii,4); % [datenum]
    t_merc_j  = jd2mjd2000(juliandate(datetime(t_merc_n,'ConvertFrom','datenum')));   % [Julian day] arrival/departure date on Venus
    t_ven_1_j = jd2mjd2000(juliandate(datetime(t_ven_1_n,'ConvertFrom','datenum')));  % [Julian day] arrival/departure date on Venus
    t_ven_2_j = jd2mjd2000(juliandate(datetime(t_ven_2_n,'ConvertFrom','datenum')));  % [Julian day] arrival/departure date on Venus
    t_jup_j   = jd2mjd2000(juliandate(datetime(t_jup_n,'ConvertFrom','datenum')));    % [Julian day] arrival/departure date on Venus

    % perform genetic algorithm based analysis
    [Data_min, Dv_min_TOF_1, Dv_min_TOF_2, r1_arc, r2_arc, r3_arc, v_ven, t_ven,vInfMin, vInfPlus, dv_ga, Data, hp] = GA_FlybyMod3( ga_it, t_merc_j, t_ven_1_j, t_ven_2_j, t_jup_j );
    DATA = [ DATA; Data_min ];
end

ga_info = profile('info');
[~,index] = sortrows([ga_info.FunctionTable.TotalTime].'); ga_info.FunctionTable = ga_info.FunctionTable(index(end:-1:1)); clear index
fprintf('Elapsed time: %s\n', duration(seconds(ga_info.FunctionTable(1).TotalTime),'Format','hh:mm:ss'));

%% PLOT
switch_plot      = 0; % plot the possible options with linear ticked axes of dates
switch_plot_zoom = 0; % plot the possible options with selected axes of dates
switch_plot_all  = 0; % plot all the points evalutated [ NaN, over the threshold, good options ]
switchPlotFlyby  = 0; % plot both unconstrained and constrained flyby solutions
switch_minima    = 1;
%%
D = Data;
% D = DataVectCon;
id = find( isnan(D(:,4)) );
Dv_nan = D(id,:);  % erase the rows corresponding to a nan value of Delta_v
D(id,:)=[];
dv_threshold = 50; % [km/s] higher value of Delta_v acceptable
id_ = find( D(:,4) > dv_threshold );
Dv_trs = D(id_,:);
D(id_,:)=[];       %  erase the rows corresponding to a too high value of Delta_v
[B,I] = sort(D);
Bx = sortrows(D,1);
for ii = 1:size(D,1)
    for jj = 1:3
        D_idx(I(ii,jj),jj) = ii;
        D_idx(I(ii,jj),jj) = ii;
        D_idx(I(ii,jj),jj) = ii;
        dates(ii,jj) =  {['',datestr(datetime(mjd20002date(B(ii,jj))), ' dd mmm yyyy HH:MM:SS')]};
    end
end
D_idx(:,4) = D(:,4); % add Delta_v values in the last column

DDate = [];
for ii = 1:size(D,1)
    DDate(ii, 1) = datenum(mjd20002date(D(ii,1)));
    DDate(ii, 2) = datenum(mjd20002date(D(ii,2)));
    DDate(ii, 3) = datenum(mjd20002date(D(ii,3)));
end
for ii = 1:size(Dv_nan,1)
    Dv_nan(ii, 1) = datenum(mjd20002date(Dv_nan(ii,1)));
    Dv_nan(ii, 2) = datenum(mjd20002date(Dv_nan(ii,2)));
    Dv_nan(ii, 3) = datenum(mjd20002date(Dv_nan(ii,3)));
end
for ii = 1:size(Dv_trs,1)
    Dv_trs(ii, 1) = datenum(mjd20002date(Dv_trs(ii,1)));
    Dv_trs(ii, 2) = datenum(mjd20002date(Dv_trs(ii,2)));
    Dv_trs(ii, 3) = datenum(mjd20002date(Dv_trs(ii,3)));
end
[min_dv,min_dv_idx] = min(D(:,4)); % find the Delta_v minimum and its index
%% plot with datenum
if( switch_plot==1 && size(D,1)~=0 )
    set_color = '#808080';
    figure('Color',set_color);
    grid on; hold on;
    scatter3( DDate(:,1), DDate(:,2), DDate(:,3), 15, D(:,4), '*')    % draw the scatter plot
    scatter3( DDate(min_dv_idx,1), DDate(min_dv_idx,2), DDate(min_dv_idx,3), 20, D(min_dv_idx,4), 'r*')    % draw the minimum Delta_V
    set(gca,'color','none') % no background for the graph
    
    datetick('x', 'yyyy mmm dd','keeplimits');
    xtickangle(-45)
    datetick('y', 'yyyy mmm dd','keeplimits');
    ytickangle(45)
    datetick('z', 'yyyy mmm dd','keeplimits');
    
    hcb = colorbar;
    hcb.Title.String = '$\Delta v$ [km/s]';
    hcb.Title.Interpreter = 'latex';
    hcb.Title.FontSize = 15;
    view(3)
    title(sprintf('Travel options evaluated with ga\n from   [%s]   to   [%s]', datetime(t_merc), datetime(t_jup) ));
    ylabel('Venus''s date');
    xlabel('Mercury''s date');
    zlabel('Jupiter''s date');
    
    dcm = datacursormode;
    dcm.Enable = 'on';
    dcm.SnapToDataVertex = 'on';
    dcm.DisplayStyle = 'datatip';
    dcm.UpdateFcn = @displayCoordinates_3d;
end
%% plot only date evaluated

if( switch_plot_zoom==1 && size(D,1)~=0 )
    set_color = '#808080';
    
    figure('Color',set_color);
    grid on; hold on;
    scatter3( D_idx(:,1), D_idx(:,2), D_idx(:,3), 15, D_idx(:,4), '*')    % draw the scatter plot
    set(gca,'color','none') % no background for the graph
    scatter3( D_idx(min_dv_idx,1), D_idx(min_dv_idx,2), D_idx(min_dv_idx,3), 20, D(min_dv_idx,4), 'r*')    % draw the minimum Delta_V
    
    t_s = floor(size(D,1)/5); %tick span
    tick_vect = 0:t_s:(size(D,1)-1);
    xticks(tick_vect);
    xticklabels(dates(tick_vect+1,1))
    xtickangle(-45)
    yticks(tick_vect);
    yticklabels(dates(tick_vect+1,2))
    ytickangle(45)
    zticks(tick_vect);
    zticklabels(dates(tick_vect+1,3))
    ztickangle(0)
    
    hcb = colorbar;
    hcb.Title.String = '$\Delta v$ [km/s]';
    hcb.Title.Interpreter = 'latex';
    hcb.Title.FontSize = 15;
    view(3)
    title(sprintf('Travel options evaluated with ga\n from   [%s]   to   [%s]', datetime(t_merc), datetime(t_jup) ));
    ylabel('Venus''s date');
    xlabel('Mercury''s date');
    zlabel('Jupiter''s date');
end
%% plot all the population of points evaluated

if( switch_plot_all )
    set_color = '#808080';
    figure('Color',set_color);
    grid on; hold on;
    scatter3( Dv_nan(:,1), Dv_nan(:,2), Dv_nan(:,3), 5, 'w.','LineWidth',3);         % scatter plot of nan points evaluated
    scatter3( Dv_trs(:,1), Dv_trs(:,2), Dv_trs(:,3), 10, 'ko','LineWidth',3);        % scatter plot of points with too high Delta_v value
    if( size(D,1) )
        scatter3( DDate(:,1), DDate(:,2), DDate(:,3), 25, D(:,4), 'o','LineWidth',2);end % scatter plot of possible options
    set(gca,'color','none') % no background for the graph
    xlabel('Mercury''s date');
    ylabel('Venus''s date');
    zlabel('Jupiter''s date');
    title(sprintf('Overall points evaluated with ga\n from   [%s]   to   [%s]', datetime(t_merc), datetime(t_jup) ));
    view(3)
    datetick('x', 'yyyy mmm dd','keeplimits');
    xtickangle(-45)
    datetick('y', 'yyyy mmm dd','keeplimits');
    ytickangle(45)
    datetick('z', 'yyyy mmm dd','keeplimits');
    legend({'nan','over the threshold','acceptable solution'},'Location','north','NumColumns',3)
    hcb = colorbar;
    hcb.Title.String = '$\Delta v$ [km/s]';
    hcb.Title.Interpreter = 'latex';
    hcb.Title.FontSize = 15;
end
%% plot needed
if( switch_plot_all )
set_color = '#808080';
figure('Color',set_color);
grid on; hold on;
scatter3( ones(size(D,1),1)*datenum(mjd20002date(t_merc_j)), DDate(:,2), ones(size(D,1),1)*datenum(mjd20002date(t_jup_j)), 15, D(:,4), '*') % draw the scatter plot
% scatter3( datenum(mjd20002date(t_merc_j)), datenum(mjd20002date(t_ven_j)), datenum(mjd20002date(t_jup_j)), 20, D(min_dv_idx,4), 'wo')       % draw the initial guess Delta_V
scatter3( datenum(mjd20002date(t_merc_j)), DDate(min_dv_idx,2), datenum(mjd20002date(t_jup_j)), 20, D(min_dv_idx,4), 'r*')                  % draw the minimum Delta_V
% scatter3( DDate(min_dv_idx,1), DDate(min_dv_idx,2), DDate(min_dv_idx,3), 20, D(min_dv_idx,4), 'r*')    % draw the minimum Delta_V
set(gca,'color','none') % no background for the graph

datetick('x', 'yyyy mmm dd','keeplimits');
xtickangle(-45)
datetick('y', 'yyyy mmm dd','keeplimits');
ytickangle(45)
datetick('z', 'yyyy mmm dd','keeplimits');

hcb = colorbar;
hcb.Title.String = '$\Delta v$ [km/s]';
hcb.Title.Interpreter = 'latex';
hcb.Title.FontSize = 15;
view(3)
title(sprintf('Travel options shifting Venus date evaluated with ga\n Mercury''s date   [%s]   Jupiter''s date   [%s]', datetime(t_merc), datetime(t_jup) ));
ylabel('Venus''s date');
xlabel('Mercury''s date');
zlabel('Jupiter''s date');
dcm = datacursormode;
dcm.Enable = 'on';
dcm.SnapToDataVertex = 'on';
dcm.DisplayStyle = 'datatip';
dcm.UpdateFcn = @displayCoordinates_3d;
rotate3d on
end
%% plot minima
if( switch_minima )
id = find( isnan(DATA(:,1)) );
DATA(id,:)=[];

    for ii = 1:size(DATA,1)
    DATA(ii, 5) = datenum(mjd20002date(DATA(ii,1)));
    DATA(ii, 6) = datenum(mjd20002date(DATA(ii,2)));
    DATA(ii, 7) = datenum(mjd20002date(DATA(ii,3)));
end

set_color = '#808080';
figure('Color',set_color);
grid on; hold on;
scatter3( DATA(:,5), DATA(:,6), DATA(:,7), 15, DATA(:,4), '*') % draw the scatter plot
set(gca,'color','none') % no background for the graph

datetick('x', 'yyyy mmm dd','keeplimits');
xtickangle(-45)
datetick('y', 'yyyy mmm dd','keeplimits');
ytickangle(45)
datetick('z', 'yyyy mmm dd','keeplimits');

hcb = colorbar;
hcb.Title.String = '$\Delta v$ [km/s]';
hcb.Title.Interpreter = 'latex';
hcb.Title.FontSize = 15;
view(3)
ylabel('Venus''s date');
xlabel('Mercury''s date');
zlabel('Jupiter''s date');
dcm = datacursormode;
dcm.Enable = 'on';
dcm.SnapToDataVertex = 'on';
dcm.DisplayStyle = 'datatip';
dcm.UpdateFcn = @displayCoordinates_3d;
rotate3d on
end
%% Auxiliary functions
function txt = displayCoordinates_3d(~,info)
x = info.Position(1);
y = info.Position(2);
z = info.Position(3);
txt = {['Departure date {Mercury} = ', datestr(x, ' dd mmm yyyy HH:MM:SS')], ['Venus''s date = ', datestr(y, ' dd mmm yyyy HH:MM:SS')], ['Arrival date {Jupiter} = ', datestr(z, ' dd mmm yyyy HH:MM:SS')]};
end