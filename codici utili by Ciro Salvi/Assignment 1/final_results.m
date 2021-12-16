% final_results represent the data evaluated with the global ga 
%
% PROTOTYPE:
% final_results
%
% CONTRIBUTORS:
% Alessandro Staffolani
%
if (isequal(exist('switchNew'),0) || switchNew == 0)
close all; clear all; clc

% DataVect loading
load( 'workspace-27-54.mat', 'DataVect' );
DataVect_27_54 = DataVect;
load( 'workspace-46-67.mat', 'DataVect' );
DataVect_46_67 = DataVect;
DataUnc = [ DataVect_27_54; DataVect_46_67 ];

% DataVectCon loading
load( 'workspace-27-54.mat', 'DataVectCon' );
DataVectCon_27_54 = DataVectCon;
load( 'workspace-46-67.mat', 'DataVectCon' );
DataVectCon_46_67 = DataVectCon;
DataCon = [ DataVectCon_27_54; DataVectCon_46_67 ];

else

DataUnc = DataVect;
DataCon = DataVectCon;

end

%% IterationMinimaIsolation

% N.B. this problem accounts for 'PopulationSize' for ga optimizator set to 50, the default value
PopulationSize = 50; % CTT, must coincide with 'PopulationSize' used for ga

DataUncMin = IterationMinimaIsolate(DataUnc, PopulationSize);
DataConMin = IterationMinimaIsolate(DataCon, PopulationSize);

%% PLOT
% First windows to catch the first date guess
tstart_dep = [ 2027, 6, 1, 0, 0, 0];      % [Gregorian date] first value of the departure window
tend_arr   = [ 2067, 6, 1, 0, 0, 0];      % [Gregorian date] final value of the arrival window

switch_plot      = 1; % plot the possible options with linear ticked axes of dates
switch_plot_zoom = 1; % plot the possible options with selected axes of dates
switch_plot_all  = 1; % plot all the points evalutated [ NaN, over the threshold, good options ]
%%
D =[];
D = DataUnc;
id = find( isnan(D(:,4)) );
Dv_nan = D(id,:);  % erase the rows corresponding to a nan value of Delta_v
D(id,:)=[];
dv_threshold = 50; % [km/s] higher value of Delta_v acceptable
id_ = find( D(:,4) > dv_threshold );
Dv_trs = D(id_,:);
D(id_,:)=[];       %  erase the rows corresponding to a too high value of Delta_v
[B,I] = sort(D);
B1x = sortrows(D,1);
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
min_dates = D(min_dv_idx, [1:3]);

%% plot with datenum
if( switch_plot==1 && size(D,1)~=0 )
set_color = '#808080';
figure('Color',set_color); 
grid on; hold on;
scatter3( DDate(:,1), DDate(:,2), DDate(:,3), 15, D(:,4), '*','LineWidth',2)    % draw the scatter plot
scatter3( DDate(min_dv_idx,1), DDate(min_dv_idx,2), DDate(min_dv_idx,3), 20, D(min_dv_idx,4), 'r*','LineWidth',2)    % draw the minimum Delta_V
% scatter3( DDate(:,1), DDate(:,2)-DDate(:,1), DDate(:,3)-DDate(:,2), 15, D(:,4), '*','LineWidth',2)    % draw the scatter plot
% scatter3( DDate(min_dv_idx,1), DDate(min_dv_idx,2)-DDate(min_dv_idx,1), DDate(min_dv_idx,3)-DDate(min_dv_idx,2), 20, D(min_dv_idx,4), 'r*','LineWidth',2)    % draw the minimum Delta_V
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
% ylabel('ToF to Venus');
% zlabel('ToF to Jupiter');

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
scatter3( D_idx(:,1), D_idx(:,2), D_idx(:,3), 15, D_idx(:,4), '*','LineWidth',2)    % draw the scatter plot
% scatter3( D_idx(:,1), D_idx(:,2)-D_idx(:,1), D_idx(:,3)-D_idx(:,2), 15, D_idx(:,4), '*','LineWidth',2)    % draw the scatter plot
set(gca,'color','none') % no background for the graph
scatter3( D_idx(min_dv_idx,1), D_idx(min_dv_idx,2), D_idx(min_dv_idx,3), 20, D(min_dv_idx,4), 'r*','LineWidth',2)    % draw the minimum Delta_V
% scatter3( D_idx(min_dv_idx,1), D_idx(min_dv_idx,2)-D_idx(min_dv_idx,1), D_idx(min_dv_idx,3)-D_idx(min_dv_idx,2), 20, D(min_dv_idx,4), 'r*','LineWidth',2)    % draw the minimum Delta_V

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
% ylabel('ToF to Venus');
% zlabel('ToF to Jupiter');
end
%% plot all the population of points evaluated

if( switch_plot_all )
set_color = '#808080';
figure('Color',set_color); 
grid on; hold on; 
scatter3( Dv_nan(:,1), Dv_nan(:,2), Dv_nan(:,3), 5, 'w.','LineWidth',3);         % scatter plot of nan points evaluated
% scatter3( Dv_nan(:,1), Dv_nan(:,2)-Dv_nan(:,1), Dv_nan(:,3)-Dv_nan(:,2), 10, 'wo','LineWidth',2);         % scatter plot of nan points evaluated
scatter3( Dv_trs(:,1), Dv_trs(:,2), Dv_trs(:,3), 10, 'ko','LineWidth',3);        % scatter plot of points with too high Delta_v value
% scatter3( Dv_trs(:,1), Dv_trs(:,2)-Dv_trs(:,1), Dv_trs(:,3)-Dv_trs(:,2), 10, 'k+','LineWidth',2);        % scatter plot of points with too high Delta_v value
if( size(D,1) )
scatter3( DDate(:,1), DDate(:,2), DDate(:,3), 25, D(:,4), 'o','LineWidth',2);
% scatter3( DDate(:,1), DDate(:,2)-DDate(:,1), DDate(:,3)-DDate(:,2), 25, D(:,4), '+','LineWidth',2);
end % scatter plot of possible options
set(gca,'color','none') % no background for the graph
xlabel('Mercury''s date');
ylabel('Venus''s date');
zlabel('Jupiter''s date');
% ylabel('ToF to Venus');
% zlabel('ToF to Jupiter');
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

% ga_info = profile('info');
% [~,index] = sortrows([ga_info.FunctionTable.TotalTime].'); ga_info.FunctionTable = ga_info.FunctionTable(index(end:-1:1)); clear index
% fprintf('Time enlapsed: %s\n', duration(seconds(ga_info.FunctionTable(1).TotalTime),'Format','hh:mm:ss'));

%%
D =[]; D_idx=[];
D = DataCon;
id = find( isnan(D(:,4)) );
Dv_nan = D(id,:);  % erase the rows corresponding to a nan value of Delta_v
D(id,:)=[];
dv_threshold = 50; % [km/s] higher value of Delta_v acceptable
id_ = find( D(:,4) > dv_threshold );
Dv_trs = D(id_,:);
D(id_,:)=[];       %  erase the rows corresponding to a too high value of Delta_v
[B,I] = sort(D);
B2x = sortrows(D,1);
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
min_dates = D(min_dv_idx, [1:3]);

%% CONSTRAINED plot with datenum
if( switch_plot==1 && size(D,1)~=0 )
set_color = '#808080';
figure('Color',set_color); 
grid on; hold on;
scatter3( DDate(:,1), DDate(:,2)-DDate(:,1), DDate(:,3)-DDate(:,2), 15, D(:,4), '+','LineWidth',2)    % draw the scatter plot
scatter3( DDate(min_dv_idx,1), DDate(min_dv_idx,2)-DDate(min_dv_idx,1), DDate(min_dv_idx,3)-DDate(min_dv_idx,2), 20, D(min_dv_idx,4), 'r*','LineWidth',2)    % draw the minimum Delta_V
set(gca,'color','none') % no background for the graph

datetick('x', 'yyyy mmm dd','keeplimits');
xtickangle(-45)
% datetick('y', 'yyyy mmm dd','keeplimits');
% ytickangle(45)
% datetick('z', 'yyyy mmm dd','keeplimits');

hcb = colorbar;
hcb.Title.String = '$\Delta v$ [km/s]';
hcb.Title.Interpreter = 'latex';
hcb.Title.FontSize = 15;
view(3)
title(sprintf('Travel options evaluated with ga\n from   [%s]   to   [%s]', datetime(tstart_dep), datetime(tend_arr) ));
xlabel('Mercury''s date');
% ylabel('Venus''s date');
% zlabel('Jupiter''s date');
ylabel('ToF to Venus');
zlabel('ToF to Jupiter');

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
% scatter3( D_idx(:,1), D_idx(:,2), D_idx(:,3), 15, D_idx(:,4), '*','LineWidth',2)    % draw the scatter plot
scatter3( D_idx(:,1), D_idx(:,2)-D_idx(:,1), D_idx(:,3)-D_idx(:,2), 15, D_idx(:,4), '*','LineWidth',2)    % draw the scatter plot
set(gca,'color','none') % no background for the graph
% scatter3( D_idx(min_dv_idx,1), D_idx(min_dv_idx,2), D_idx(min_dv_idx,3), 20, D(min_dv_idx,4), 'r*','LineWidth',2)    % draw the minimum Delta_V
scatter3( D_idx(min_dv_idx,1), D_idx(min_dv_idx,2)-D_idx(min_dv_idx,1), D_idx(min_dv_idx,3)-D_idx(min_dv_idx,2), 20, D(min_dv_idx,4), 'r*','LineWidth',2)    % draw the minimum Delta_V

t_s = floor(size(D,1)/5); %tick span
tick_vect = 0:t_s:(size(D,1)-1);
xticks(tick_vect);
xticklabels(dates(tick_vect+1,1))
xtickangle(-45)
% yticks(tick_vect);
% yticklabels(dates(tick_vect+1,2))
% ytickangle(45)
% zticks(tick_vect);
% zticklabels(dates(tick_vect+1,3))
% ztickangle(0)

hcb = colorbar;
hcb.Title.String = '$\Delta v$ [km/s]';
hcb.Title.Interpreter = 'latex';
hcb.Title.FontSize = 15;
view(3)
title(sprintf('Travel options evaluated with ga\n from   [%s]   to   [%s]', datetime(tstart_dep), datetime(tend_arr) ));
ylabel('Venus''s date');
% xlabel('Mercury''s date');
% zlabel('Jupiter''s date');
ylabel('ToF to Venus');
zlabel('ToF to Jupiter');
end
%% CONSTRAINED plot all the population of points evaluated 

if( switch_plot_all )
set_color = '#808080';
figure('Color',set_color); 
grid on; hold on; 
% scatter3( Dv_nan(:,1), Dv_nan(:,2), Dv_nan(:,3), 5, 'w.','LineWidth',3);         % scatter plot of nan points evaluated
scatter3( Dv_nan(:,1), Dv_nan(:,2)-Dv_nan(:,1), Dv_nan(:,3)-Dv_nan(:,2), 10, 'wo','LineWidth',2);         % scatter plot of nan points evaluated
% scatter3( Dv_trs(:,1), Dv_trs(:,2), Dv_trs(:,3), 10, 'ko','LineWidth',3);        % scatter plot of points with too high Delta_v value
scatter3( Dv_trs(:,1), Dv_trs(:,2)-Dv_trs(:,1), Dv_trs(:,3)-Dv_trs(:,2), 10, 'k+','LineWidth',2);        % scatter plot of points with too high Delta_v value
if( size(D,1) )
% scatter3( DDate(:,1), DDate(:,2), DDate(:,3), 25, D(:,4), 'o','LineWidth',2);
scatter3( DDate(:,1), DDate(:,2)-DDate(:,1), DDate(:,3)-DDate(:,2), 25, D(:,4), '+','LineWidth',2);
end % scatter plot of possible options
set(gca,'color','none') % no background for the graph
xlabel('Mercury''s date');
% ylabel('Venus''s date');
% zlabel('Jupiter''s date');
ylabel('ToF to Venus');
zlabel('ToF to Jupiter');
title(sprintf('Overall points evaluated with ga\n from   [%s]   to   [%s]', datetime(tstart_dep), datetime(tend_arr) ));
view(3)
datetick('x', 'yyyy mmm dd','keeplimits');
xtickangle(-45)
% datetick('y', 'yyyy mmm dd','keeplimits');
% ytickangle(45)
% datetick('z', 'yyyy mmm dd','keeplimits');
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