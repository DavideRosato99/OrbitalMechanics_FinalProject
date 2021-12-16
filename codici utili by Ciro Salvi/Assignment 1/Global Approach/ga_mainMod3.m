% ga_mainMod3 a given transfer window is analyzed using shifting
% sub-window, assigned the shifting step and the window width
%
% PROTOTYPE:
% ga_mainMod3
%
% CONTRIBUTORS:
%   Fabio Spada
%   Suhailah Alkhawashke
%


close all; clear all; clc
profile on
%%
% First windows to catch the first date guess
tstart_dep = [ 2027, 6, 1, 0, 0, 0];      % [Gregorian date] first value of the departure window
tend_arr   = [ 2054, 6, 1, 0, 0, 0];      % [Gregorian date] final value of the arrival window

% Transform the date in Modified Julian day 2000 number from Gregorian date
t_start_j_1 = date2mjd2000(tstart_dep);   % [Julian day] first value of the departure window
t_end_j_2   = date2mjd2000(tend_arr);     % [Julian day] final value of the arrival window

% initialize vectors of useful variables
DvminVect = [];
DataVect = [];
hpVect = [];
vInfMinVect = [];
vInfPlusVect = [];

DvminVectCon = [];
DataVectCon = [];
hpVectCon = [];
vInfMinConVect = [];
vInfPlusConVect = [];

timeStep = 30; % day step between each starting point
nYears = 6;
timeSpan = 365.25*nYears;
nGroup = ceil((t_end_j_2-t_start_j_1-timeSpan)/timeStep);                        
ga_it = 10;         

n_it = ga_it*nGroup*2;

tWind = [t_start_j_1 t_start_j_1+timeSpan];

for ii = 1 : nGroup
    
    fprintf('Group N°%d \n', ii) 
    
    [DV_MINCon, Dv_min_TOF_1Con, Dv_min_TOF_2Con, r1_arcCon, r2_arcCon, r3_arcCon, v_venCon, t_venCon,vInfMinCon, vInfPlusCon, dv_gaCon, DataCon, hpCon] = GA_Flyby_con( ga_it, tWind(1), tWind(2) );
    [DV_MIN, Dv_min_TOF_1, Dv_min_TOF_2, r1_arc, r2_arc, r3_arc, v_ven, t_ven,vInfMin, vInfPlus, dv_ga, Data, hp] = GA_FlybyMod2( ga_it, tWind(1), tWind(2) );
    
    fprintf('%.2f %% of process completed \n', ii/nGroup*100);
    
    DvminVect = [DvminVect; DV_MIN];
    DataVect = [DataVect; Data];
    hpVect = [hpVect; hp];
    vInfMinVect = [vInfMinVect; vInfMin];
    vInfPlusVect = [vInfPlusVect; vInfPlus];

    DvminVectCon = [DvminVectCon; DV_MINCon];
    DataVectCon = [DataVectCon; DataCon];
    hpVectCon = [hpVectCon; hpCon];
    vInfMinConVect = [vInfMinConVect; vInfMinCon];
    vInfPlusConVect = [vInfPlusConVect; vInfPlusCon];

    tWind = tWind + timeStep;
    
end

if( ~isempty(DataVect))
    fprintf('%d evaluated points with no contraints have not a NaN value\n', sum(~isnan(DataVect(:,4))) );end
if( ~isempty(DataVectCon))
    fprintf('%d evaluated points with double contraints haven''t a NaN value\n', sum(~isnan(DataVectCon(:,4))) );end
ga_info = profile('info');
[~,index] = sortrows([ga_info.FunctionTable.TotalTime].'); ga_info.FunctionTable = ga_info.FunctionTable(index(end:-1:1)); clear index
fprintf('Elapsed time: %s\n', duration(seconds(ga_info.FunctionTable(1).TotalTime),'Format','hh:mm:ss'));

%% compare min
if( ~isempty(DataVect) && ~isempty(DataVectCon) )
    [~,id_min] = min( DataVect(:,4) );
    [~, id_min_con] = min( DataVectCon(:,4) );
    Compare = table( [ datetime(mjd20002date(DataVect(id_min,1))); datetime(mjd20002date(DataVectCon(id_min_con,1))) ],...
                     [ datetime(mjd20002date(DataVect(id_min,2))); datetime(mjd20002date(DataVectCon(id_min_con,2))) ],...
                     [ datetime(mjd20002date(DataVect(id_min,3))); datetime(mjd20002date(DataVectCon(id_min_con,3))) ],...
                     [                         DataVect(id_min,4);                         DataVectCon(id_min_con,4) ],...
                     'VariableNames',{'Mercury','Venus','Jupiter','DV'}, 'RowNames',{'Single constraint', 'Double constraints'} );
    disp(Compare);
end
%% PLOT
switch_plot      = 1; % plot the possible options with linear ticked axes of dates
switch_plot_zoom = 1; % plot the possible options with selected axes of dates
switch_plot_all  = 1; % plot all the points evalutated [ NaN, over the threshold, good options ]
switchPlotFlyby  = 0;  
%%
clear D_idx    % <<<<<< If the plotting are based on already evaluated data
D = DataVect;
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

nGroupMin = ceil(min_dv_idx/(50*ga_it));
hpMin = hpVect(nGroupMin);
vInfMinMin = vInfMinVect(nGroupMin);
vInfPlusMin = vInfPlusVect(nGroupMin);

%% plot flyby solution

if (switchPlotFlyby)
    GravityAssistPlot(vInfMinMin, vInfPlusPlus, hpMin, 2);
    deltavInfNorm = norm(vInfPlus - vInfMin)
end

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
title(sprintf('Travel options evaluated with ga\n from   [%s]   to   [%s]', datetime(tstart_dep), datetime(tend_arr) ));
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
title(sprintf('Travel options evaluated with ga\n from   [%s]   to   [%s]', datetime(tstart_dep), datetime(tend_arr) ));
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
title(sprintf('Overall points evaluated with ga\n from   [%s]   to   [%s]', datetime(tstart_dep), datetime(tend_arr) ));
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

ga_info = profile('info');
[~,index] = sortrows([ga_info.FunctionTable.TotalTime].'); ga_info.FunctionTable = ga_info.FunctionTable(index(end:-1:1)); clear index
fprintf('Time elapsed: %s\n', duration(seconds(ga_info.FunctionTable(1).TotalTime),'Format','hh:mm:ss'));

%%

D = DataVectCon;
id = find( isnan(D(:,4)) );
Dv_nan = D(id,:);  % erase the rows corresponding to a nan value of Delta_v
D(id,:)=[];
dv_threshold = 70; % [km/s] higher value of Delta_v acceptable
id_ = find( D(:,4) > dv_threshold );
Dv_trs = D(id_,:);
D(id_,:)=[];       %  erase the rows corresponding to a too high value of Delta_v
[B,I] = sort(D);
Bx = sortrows(D,1);
D_idx = [];
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

nGroupConMin = ceil(min_dv_idx/(50*ga_it));
hpConMin = hpVect(nGroupConMin);
vInfMinConMin = vInfMinConVect(nGroupConMin);
vInfPlusConMin = vInfPlusConVect(nGroupConMin);

%% CONSTRAINED plot with datenum
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
title(sprintf('Travel options evaluated with ga\n from   [%s]   to   [%s]', datetime(tstart_dep), datetime(tend_arr) ));
ylabel('Venus''s date');
xlabel('Mercury''s date');
zlabel('Jupiter''s date');

dcm = datacursormode;
dcm.Enable = 'on';
dcm.SnapToDataVertex = 'on';
dcm.DisplayStyle = 'datatip';
dcm.UpdateFcn = @displayCoordinates_3d;
end
%% CONSTRAINED plot only date evaluated

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
title(sprintf('Travel options evaluated with ga\n from   [%s]   to   [%s]', datetime(tstart_dep), datetime(tend_arr) ));
ylabel('Venus''s date');
xlabel('Mercury''s date');
zlabel('Jupiter''s date');
end
%% CONSTRAINED plot all the population of points evaluated 

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
title(sprintf('Overall points evaluated with ga\n from   [%s]   to   [%s]', datetime(tstart_dep), datetime(tend_arr) ));
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

function txt = displayCoordinates_3d(~,info)
x = info.Position(1);
y = info.Position(2);
z = info.Position(3);
txt = {['Departure date {Mercury} = ', datestr(x, ' dd mmm yyyy HH:MM:SS')], ['Venus''s date = ', datestr(y, ' dd mmm yyyy HH:MM:SS')], ['Arrival date {Jupiter} = ', datestr(z, ' dd mmm yyyy HH:MM:SS')]};
end