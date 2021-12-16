% refine_minima this script call the refinement methods
%
% PROTOTYPE: 
%   refine_minima
%
% CONTRIBUTORS:
%   Alessandro Staffolani
%   Fabio Spada
%   Suhailah Alkhawashke
%
% VERSIONS:
%   2021-02-11
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

threshold = 23.5;                               % depend the analysis of this threshold choice
[id_unc] = find( DataUncMin(:,4)<threshold );
Data_min_unc = DataUncMin(id_unc,:);
[id_con] = find( DataConMin(:,4)<threshold );
Data_min_con = DataConMin(id_con,:);
windAmpl = 10; 

%% REFINEMENT SECTION - GridSearch
grid_search_switch = 0;
% behaviour: 
% for windAmpl = 10; timeStep = 1 the results appear to be improving for higher
% cost options but to be worsening for the least expensive options

if grid_search_switch
timeStep = 1;

DataMinUncRef = zeros(length(Data_min_unc),8);
DataMinConRef = zeros(length(Data_min_con),8);
cTimeVectUnc = [];
cTimeVectCon = [];

% unconstrained data
for ii = 1 : length(id_unc)
    TVect = Data_min_unc(ii, 1:3);
    [MinRefined, ~,~,~,~,~,cTime] = GridSearchRefine(TVect, windAmpl, timeStep);
    DataMinUncRef(ii, :) = MinRefined(1:8);
    cTimeVectUnc = [cTimeVectUnc; cTime];
end

% constrained data
for ii = 1 : length(id_con)
    TVect = Data_min_con(ii, 1:3);
    [MinRefined, ~,~,~,~,~,cTime] = GridSearchRefine(TVect, windAmpl, timeStep);
    DataMinConRef(ii, :) = MinRefined(1:8);
    cTimeVectCon = [cTimeVectCon; cTime];
end
out_gs_unc.DataMinUncRef = DataMinUncRef;
out_gs_con.DataMinConRef = DataMinConRef;
out_gs_unc.cTimeVectUnc = cTimeVectUnc;
out_gs_con.cTimeVectCon = cTimeVectCon;
end

%% REFINEMENT SECTION - GA analysis
GA_switch = 0;

if (GA_switch)
    
    % input
    ref_days = 10;
    ga_it_ref = 10;
    
    for E = 1:2
        
        options_GAref = optimoptions('ga','Display','off','Elitecount',E);
        
        fprintf('Elitecount %d \n PopulationSize %d \n',E);
        
% <<<<< Unconstained >>>>>

        fprintf('\n***Unconstrained***\n');
        tic
        GAr_abs_min = [];
        
        for mm = 1:size(Data_min_unc,1)
            % initialize vectors of useful variables
            GAr_DvminVect = [];
            GAr_DataVect = [];
            GAr_hpVect = [];
            GAr_vInfMinVect = [];
            GAr_vInfPlusVect = [];
            
            for ii = 1 : ga_it_ref
                
                fprintf('GA Refinement iteration N°%d \n', ii)
                
                [GAr_DV_MIN, GAr_Dv_min_TOF_1, GAr_Dv_min_TOF_2, GAr_r1_arc, GAr_r2_arc, GAr_r3_arc, GAr_v_venus, GAr_t_venus, GAr_vInfMin, GAr_vInfPlus, GAr_dv_ga, GAr_Data, GAr_hp] = GA_ref( ga_it_ref,ref_days,options_GAref,Data_min_unc(mm,:));
                
                fprintf('%.2f %% of process completed \n', ii/ga_it_ref*100);
                
                options_vect = repmat(E,size(GAr_Data,1),1);
                
                GAr_DvminVect = [GAr_DvminVect; GAr_DV_MIN];
                GAr_DataVect = [GAr_DataVect; GAr_Data,options_vect];
                GAr_hpVect = [GAr_hpVect; GAr_hp];
                GAr_vInfMinVect = [GAr_vInfMinVect, GAr_vInfMin]; % There was a typo here
                GAr_vInfPlusVect = [GAr_vInfPlusVect, GAr_vInfPlus]; % There was a typo here
                
            end
            
            [~, ind] = min(GAr_DataVect(:,4));
            GAr_abs_min = [GAr_abs_min;GAr_DataVect(ind,:)]; clear ind
            
            if( ~isempty(GAr_DataVect))
                fprintf('%d evaluated points with no contraints have not a NaN value\n', sum(~isnan(GAr_DataVect(:,4))) );end
        end
        c_time_GAunc = toc;
        
% <<<<< Constrained >>>>>
        fprintf('\n***Constrained***\n');
        tic
        GAr_abs_minCon = [];
        
        for mm = 1:size(Data_min_con,1)
            % initialize vectors of useful variables
            GAr_DvminVectCon = [];
            GAr_DataVectCon = [];
            GAr_hpVectCon = [];
            GAr_vInfMinConVect = [];
            GAr_vInfPlusConVect = [];
            
            for ii = 1 : ga_it_ref
                
                fprintf('GA Refinement iteration N°%d \n', ii)
                
                [GAr_DV_MINCon, GAr_Dv_min_TOF_1Con, GAr_Dv_min_TOF_2Con, GAr_r1_arcCon, GAr_r2_arcCon, GAr_r3_arcCon, GAr_v_venCon, GAr_t_venCon,GAr_vInfMinCon, GAr_vInfPlusCon, GAr_dv_gaCon, GAr_DataCon, GAr_hpCon] = GA_ref_con( ga_it_ref,ref_days,options_GAref,Data_min_con(mm,:));
                
                fprintf('%.2f %% of process completed \n', ii/ga_it_ref*100);
                
                options_vect_con = repmat(E,size(GAr_DataCon,1),1);
                
                GAr_DvminVectCon = [GAr_DvminVectCon; GAr_DV_MINCon];
                GAr_DataVectCon = [GAr_DataVectCon; GAr_DataCon,options_vect_con];
                GAr_hpVectCon = [GAr_hpVectCon; GAr_hpCon];
                GAr_vInfMinConVect = [GAr_vInfMinConVect; GAr_vInfMinCon];
                GAr_vInfPlusConVect = [GAr_vInfPlusConVect; GAr_vInfPlusCon];
                
            end
            
            [~, indC] = min(GAr_DataVectCon(:,4));
            GAr_abs_minCon = [GAr_abs_minCon;GAr_DataVectCon(indC,:)]; clear indC
            
            if( ~isempty(GAr_DataVectCon))
                fprintf('%d evaluated points with double contraints haven''t a NaN value\n', sum(~isnan(GAr_DataVectCon(:,4))) );end
            
        end
        c_time_GAcon = toc;
    end
    
    out_GAref_unc.GAr_abs_min = GAr_abs_min;
    out_GAref_con.GAr_abs_minCon = GAr_abs_minCon;
    
end

%% REFINEMENT SECTION - Gradient based
gradient_based_switch = 0;
if gradient_based_switch
    attempts = 200;
    unconstrained data
    for ii = 1 : length(id_unc)
        [ D_min_fcon, D_min_func, c_time_fmincon, c_time_fminunc ] = fmin_con_unc_refinement( Data_min_unc(ii,1:3), windAmpl, attempts );
        Data_min_unc_fcon(ii,:) = D_min_fcon(:);
        Data_min_unc_func(ii,:) = D_min_func(:);
        comp_time_unc_fmincon(ii) = c_time_fmincon;
        comp_time_unc_fminunc(ii) = c_time_fminunc;
    end
    out_gb_unc.Data_min_unc_fcon = Data_min_unc_fcon;
    out_gb_unc.Data_min_unc_func = Data_min_unc_func;
    out_gb_unc.comp_time_unc_fmincon = comp_time_unc_fmincon;
    out_gb_unc.comp_time_unc_fminunc = comp_time_unc_fminunc;
    
    % constrained data
    for ii = 1 : length(id_con)
        [ D_min_fcon, D_min_func, c_time_fmincon, c_time_fminunc ] = fmin_con_unc_refinement( Data_min_con(ii,1:3), windAmpl, attempts );
        Data_min_con_fcon(ii,:) = D_min_fcon(:);
        Data_min_con_func(ii,:) = D_min_func(:);
        comp_time_con_fmincon(ii) = c_time_fmincon;
        comp_time_con_fminunc(ii) = c_time_fminunc;
    end    
    out_gb_con.Data_min_con_fcon = Data_min_con_fcon;
    out_gb_con.Data_min_con_func = Data_min_con_func;
    out_gb_con.comp_time_con_fmincon = comp_time_con_fmincon;
    out_gb_con.comp_time_con_fminunc = comp_time_con_fminunc;
end
%% COMPARISON within methods
% deltavshift % computational time % deltav
switch_compare = 1;
if switch_compare
%loading grid_search struct
load( 'out_gs_unc.mat');
DataMinUncRef = out_gs_unc.DataMinUncRef;
cTimeVectUnc = out_gs_unc.cTimeVectUnc;
load( 'out_gs_con.mat');
DataMinConRef = out_gs_con.DataMinConRef;
cTimeVectCon = out_gs_con.cTimeVectCon;

%loading ga struct
load( 'out_GAref_unc.mat' );
GAr_abs_min = out_GAref_unc.GAr_abs_min;
c_time_GAunc = out_GAref_unc.c_time_GAunc;
load( 'out_GAref_con.mat');
GAr_abs_minCon = out_GAref_con.GAr_abs_minCon;
c_time_GAcon = out_GAref_con.c_time_GAcon;

%loading gradient_based struct
load( 'out_gb_unc.mat');
Data_min_unc_fcon = out_gb_unc.Data_min_unc_fcon;
Data_min_unc_func = out_gb_unc.Data_min_unc_func;
comp_time_unc_fmincon = out_gb_unc.comp_time_unc_fmincon;
comp_time_unc_fminunc = out_gb_unc.comp_time_unc_fminunc;
load( 'out_gb_con.mat');
Data_min_con_fcon = out_gb_con.Data_min_con_fcon;
Data_min_con_func = out_gb_con.Data_min_con_func;
comp_time_con_fmincon = out_gb_con.comp_time_con_fmincon;
comp_time_con_fminunc = out_gb_con.comp_time_con_fminunc;

% constrained_data compare
figure; grid on; hold on;
plot(Data_min_con(:,4),'-.o', 'DisplayName', '$$Global\:GA$$');
plot(DataMinConRef(:,1),'-.o', 'DisplayName', '$$Grid\:search$$');
plot(Data_min_con_fcon(:,4),'-.o', 'DisplayName', '$$Fmincon$$');
plot(Data_min_con_func(:,4),'-.o', 'DisplayName', '$$Fminunc$$');
plot(GAr_abs_minCon(:,4),'-.o', 'DisplayName', '$$GA$$');
legend('Interpreter','latex');
set(gca,'TickLabelInterpreter','latex')
xlabel('$$ID\:Minima\:selected$$','Interpreter','latex');
ylabel('$$\Delta V\;\left [ \frac{km}{s} \right ]$$','Interpreter','latex');
title('\textbf{$$Refinement\: of\: constrained\: data\: evaluated\: with\: global\: ga $$}','Interpreter','latex')
figure; grid on; hold on;
plot(cTimeVectCon(:,1),'-.o', 'DisplayName', '$$Grid\:search$$');
plot(comp_time_con_fmincon,'-.o', 'DisplayName', '$$Fmincon$$');
plot(comp_time_con_fminunc,'-.o', 'DisplayName', '$$Fminunc$$');
legend('Interpreter','latex');
set(gca,'TickLabelInterpreter','latex')
xlabel('$$ID\:Minima\:selected$$','Interpreter','latex');
ylabel('$$time\;\left [ sec \right ]$$','Interpreter','latex');
title('\textbf{$$computational\; time\; taken\; by\; each\; method $$}','Interpreter','latex')

% unconstrained_data compare
figure; grid on; hold on;
plot(Data_min_unc(:,4),'-.o', 'DisplayName', '$$Global\:GA$$');
plot(DataMinUncRef(:,1),'-.o', 'DisplayName', '$$Grid\:search$$');
plot(Data_min_unc_fcon(:,4),'-.o', 'DisplayName', '$$Fmincon$$');
plot(Data_min_unc_func(:,4),'-.o', 'DisplayName', '$$Fminunc$$');
plot(GAr_abs_min(:,4),'-.o', 'DisplayName', '$$GA$$');
legend('Interpreter','latex');
set(gca,'TickLabelInterpreter','latex')
xlabel('$$ID\:Minima\:selected$$','Interpreter','latex');
ylabel('$$\Delta V\;\left [ \frac{km}{s} \right ]$$','Interpreter','latex');
title('\textbf{$$Refinement\: of\: unconstrained\: data\: evaluated\: with\: global\: ga $$}','Interpreter','latex')
figure; grid on; hold on;
plot(cTimeVectUnc(:,1),'-.o', 'DisplayName', '$$Grid\:search$$');
plot(comp_time_unc_fmincon,'-.o', 'DisplayName', '$$Fmincon$$');
plot(comp_time_unc_fminunc,'-.o', 'DisplayName', '$$Fminunc$$');
legend('Interpreter','latex');
set(gca,'TickLabelInterpreter','latex')
xlabel('$$ID\:Minima\:selected$$','Interpreter','latex');
ylabel('$$time\;\left [ sec \right ]$$','Interpreter','latex');
title('\textbf{$$computational\; time\; taken\; by\; each\; method $$}','Interpreter','latex')

Compare = table( [          (sum(cTimeVectUnc));          (sum(cTimeVectCon)) ],...
                 [ (sum(comp_time_unc_fmincon)); (sum(comp_time_con_fmincon)) ],...
                 [ (sum(comp_time_unc_fminunc)); (sum(comp_time_con_fminunc)) ],...
                 [               (c_time_GAunc);               (c_time_GAcon) ],...
                     'VariableNames',{'Grid_search','Fmincon','Fminunc','GA'}, 'RowNames',{'Single constraint', 'Double constraints'} );
disp(Compare);

end
%%
switch_plot = 0;
if switch_plot
figure; grid minor; hold on;

plot( Data_min_unc(:,1), Data_min_unc(:,4), '-xr' )
plot( Data_min_con(:,1), Data_min_con(:,4), '-xb' )
datetick('x', 'yyyy mmm dd','keeplimits');
xtickangle(45)
legend('unconstrained','constrained')

dcm = datacursormode;
dcm.Enable = 'on';
dcm.SnapToDataVertex = 'on';
dcm.DisplayStyle = 'datatip';
dcm.UpdateFcn = @displayCoordinates_2d;
%%
set_color = '#808080';

figure('Color',set_color);
grid on; hold on;
% scatter3( Data_min_unc(:,1), Data_min_unc(:,2), Data_min_unc(:,3), 15, Data_min_unc(:,4), '*','LineWidth',2)    % draw the scatter plot
scatter3( Data_min_unc(:,1), Data_min_unc(:,2)-Data_min_unc(:,1), Data_min_unc(:,3)-Data_min_unc(:,2), 15, Data_min_unc(:,4), '*','LineWidth',2)    % draw the scatter plot
% scatter3( Data_min_con(:,1), Data_min_con(:,2), Data_min_con(:,3), 15, Data_min_con(:,4), '*','LineWidth',2)    % draw the scatter plot
scatter3( Data_min_con(:,1), Data_min_con(:,2)-Data_min_con(:,1), Data_min_con(:,3)-Data_min_con(:,2), 15, Data_min_con(:,4), '*','LineWidth',2)    % draw the scatter plot
set(gca,'color','none') % no background for the graph

% t_s = floor(size(D,1)/5); %tick span
% tick_vect = 0:t_s:(size(D,1)-1);
% xticks(tick_vect);
% xticklabels(dates(tick_vect+1,1))
% xtickangle(-45)
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
xlabel('Mercury''s date');
% ylabel('Venus''s date');
% zlabel('Jupiter''s date');
ylabel('ToF to Venus');
zlabel('ToF to Jupiter');
end
function txt = displayCoordinates_2d(~,info)
x = info.Position(1);
dv = info.Position(2);
txt = {['Departure date {Mercury} = ', datestr(x, ' dd mmm yyyy HH:MM:SS')],['\DeltaV = ', num2str(dv),' [km/s]']};
end
