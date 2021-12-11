%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             ASSIGNMENT 1                                %
%                       Interplanetary Explorer Mission                   %
%                                                                         %
%                               GROUP 5                                   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% INITIAL SETUP

close all
clear all
set(0,'defaultTextInterpreter','latex')
set(0,'DefaultAxesFontSize', 14);

addpath('time'); % Subfolder containing time-conversion functions
addpath('textures'); % Subfolder containing planet-texture images


mu_s = astroConstants(4); % Gravitational parameter of the Sun [km^3/s^2]


dep_ID = 8;  % Neptune ID [-]
fb_ID = 5;   % Jupiter ID [-]
arr_ID = 4;  % Mars ID    [-]
sun_ID = 10;

mu_d  = astroConstants(10 + dep_ID); % Gravitational parameter of departure planet [km^3/s^2]
mu_fb = astroConstants(10 + fb_ID);  % Gravitational parameter of flyby planet     [km^3/s^2]
mu_a  = astroConstants(10 + arr_ID); % Gravitational parameter of arrival planet   [km^3/s^2]


kepElems = uplanet(0, dep_ID);
T_1_JD = 1/(86400) * 2*pi*sqrt(kepElems(1)^3/mu_s);

kepElems  = uplanet(0, fb_ID);
T_2_JD = 1/(86400)* 2*pi*sqrt(kepElems(1)^3/mu_s);

kepElems = uplanet(0, arr_ID);
T_3_JD = 1/(86400) * 2*pi*sqrt(kepElems(1)^3/mu_s);




T_syn_DwrtFB_JD = T_1_JD*T_2_JD / abs(T_1_JD-T_2_JD); % Synodic period of Departure planet w.r.t Flyby planet
T_syn_FBwrtA_JD = T_2_JD*T_3_JD / abs(T_2_JD-T_3_JD); % Synodic period of Flyby planet w.r.t Arrival planet

% Departure Window

dep_min_date = [2033, 05, 01, 00, 00, 00]'   % [year, month, day, hours, minutes, seconds]
dep_min = date2mjd2000(dep_min_date);
dep_max_date = [2033, 07, 01, 00, 00, 00]'
dep_max = date2mjd2000(dep_max_date);

% TOF1 Window (Departure to Flyby)
TOF1_min = 0.1*T_syn_DwrtFB_JD
TOF1_max = 1.3*T_syn_DwrtFB_JD


% TOF2 Window (Flyby to Arrival)
TOF2_min = 0.1*T_syn_FBwrtA_JD
TOF2_max = 3*T_syn_FBwrtA_JD

% Time step-size
dt = date2mjd2000([2000, 01, 06, 12, 00, 00]); % [MJD2000]




% ******************* SELECTED TRAJECTORY *******************
date1 = [2033, 05, 29, 23, 23, 04]';   % [year, month, day, hours, minutes, seconds]
t1_sel = date2mjd2000(date1);

date2 = [2050, 01, 11, 15, 21, 04]';
t2_sel = date2mjd2000(date2);
TOF1_sel = t2_sel-t1_sel;

date3 = [2054, 05, 15, 14, 35, 30]';
t3_sel = date2mjd2000(date3);
TOF2_sel = t3_sel-t2_sel;
% ******************* SELECTED TRAJECTORY *******************



departure_W  = [dep_min:dt:dep_max]';
TOF1_W   = [TOF1_min:dt:TOF1_max]';
TOF2_W   = [TOF2_min:dt:TOF2_max]';



fprintf('3D array size: %i x %i x %i\nElements No.: %i\n', length(departure_W), length(TOF1_W), length(TOF2_W),...
                                                         prod([length(departure_W), length(TOF1_W), length(TOF2_W)]))
                                                     
            
prompt = sprintf('\n\n\n1. Run Genetic Algorithm\n2. Run Triple-Loop\n3. Flyby Detail\n4. Porkchop plots\n5. Evaluate selected trajectory\n6. Render movie\n0. Exit\nSelect an option: ');                                                     
prompt_GA = sprintf('\n\n\n1. Single dep. window domain\n2. Split dep. window into sub-domains\nSelect an option: ');                                                     


myInp_1 = 1;

while myInp_1
    
    myInp_1 = input(prompt);

    switch myInp_1
        
        case 1

            
        myInp_2 = input(prompt_GA);
        switch myInp_2
            
            case 1

                                                     
%% GA: Single Dep. Window Domain

lowerBound = [dep_min TOF1_min TOF2_min];
upperBound = [dep_max TOF1_max TOF2_max];
options = optimoptions('ga','PlotFcn', @gaplotbestf, 'FunctionTolerance', 1e-4, 'PopulationSize', 2000);
[timeVec_genAlg, DeltaV_genAlg, exitflag, output, population, scores] = ga(@fitnessFcnDeltaV,3,[],[],[],[],lowerBound,upperBound, [],options);

fprintf('Minimum DeltaV arc w/ cost: %f km/s\n', DeltaV_genAlg)
fprintf('Departure: %s\n', date2string(mjd20002date(timeVec_genAlg(1))))
fprintf('Flyby: %s\n', date2string(mjd20002date(sum(timeVec_genAlg(1:2)))))
fprintf('Arrival: %s\n', date2string(mjd20002date(sum(timeVec_genAlg))))

            case 2

%% GA: Subdomains
for i = 1:4
N_blocks = 80;
departure_W_blockBounds = linspace(min(departure_W), max(departure_W), N_blocks)';

    DATA = {'Departure','Flyby','Arrival','Delta_V'};
for b = 1:length(departure_W_blockBounds)-1
    
    dep_min = departure_W_blockBounds(b);
    dep_max = departure_W_blockBounds(b+1);
    
    lowerBound = [dep_min TOF1_min TOF2_min];
    upperBound = [dep_max TOF1_max TOF2_max];
    options = optimoptions('ga','PlotFcn', @gaplotbestf, 'FunctionTolerance', 1e-4, 'PopulationSize', 2000);
    [timeVec_genAlg, DeltaV_genAlg, exitflag, output, population, scores] = ga(@fitnessFcnDeltaV,3,[],[],[],[],lowerBound,upperBound, [],options);
    
    DATA{b+1,1} = sprintf(date2string(mjd20002date(timeVec_genAlg(1))));
    DATA{b+1,2} = sprintf(date2string(mjd20002date(sum(timeVec_genAlg(1:2)))));
    DATA{b+1,3} = sprintf(date2string(mjd20002date(sum(timeVec_genAlg))));
    DATA{b+1,4} = DeltaV_genAlg;
    
    DATA{b+1,:}
    
end


% Convert cell to a table and use first row as variable names
myTable = cell2table(DATA(2:end,:),'VariableNames',DATA(1,:))
 
% Write the table to a CSV file
writetable(myTable, sprintf('wideRangeLaunchWindow_GASubDomains (%s).csv', string(datetime('now'))));


end


%% Merge Multiple Genetic-Alg. Runs
directory_GA = '4T_syn';

Files=dir( sprintf('%s/*.csv', directory_GA));
for k=1:length(Files)
   fileName = Files(k).name;
   fileDir = sprintf('%s/%s', directory_GA, fileName);
   DeltaV(:,k) = csvread(fileDir, 1,3);
   
   fid = fopen(fileDir, 'rt');
   H = textscan(fid, '%s ', 'delimiter', ',', 'MultipleDelimsAsOne',1);
   H = H{1};
   H = reshape(H, 4, numel(H)/4)';
   myCell{k} = H;
end

[DV, min_idx] = min(DeltaV, [], 2)

subCellSize = size(myCell{1});

    DATA = {'Departure','Flyby','Arrival','Delta_V'};
for b = 1:subCellSize(1)-1
    H = myCell{min_idx(b)};
    DATA(b+1,:) = H(b+1,:);
    
end


% Convert cell to a table and use first row as variable names
myTable = cell2table(DATA(2:end,:),'VariableNames',DATA(1,:))
 
% Write the table to a CSV file
writetable(myTable, sprintf('GA_combinedResults_%s.csv', directory_GA));
            
        end
        
        case 2


%% Brute-force triple loop


% Matrix of DeltaV values to achieve the insertion of the S/C into
% the 1st Lambert's arc (dep2flyby)
DeltaV1_mat = zeros(length(departure_W), length(TOF1_W));

% Matrix of DeltaV values to perform the powered flyby, given initial and
% final conditions, which are respectively specified by the incoming and
% the outcoming heliocentric trajectories (i.e. the two Lambert's arcs)
DeltaV2_mat = zeros(length(departure_W), length(TOF1_W), length(TOF2_W));

% Matrix of DeltaV values to achieve the insertion of the S/C into
% the final orbit (flyby2arrival)
DeltaV3_mat = zeros(length(departure_W), length(TOF1_W));
DeltaV_tot_tens = zeros(length(departure_W), length(TOF1_W), length(TOF2_W));


lastMsgNumel = 0; % Initialize counter of last message-to-delete elements (right before cycling)
for i = 1:length(departure_W)
    t1 = departure_W(i); % Departure time
    
    kep_dep = uplanet(t1, dep_ID);
    [RR1, VV1] = kep2car(kep_dep(1), kep_dep(2), kep_dep(3), kep_dep(4), kep_dep(5), kep_dep(6), mu_s);
    
    for j = 1:length(TOF1_W)
        
        TOF1 = TOF1_W(j);
        t2 = t1 + TOF1; % Gravity Assist time
        ToF_seconds = TOF1 * (24*3600); % [s]
        kep_arr = uplanet(t2, fb_ID);
        [RR2, VV2] = kep2car(kep_arr(1), kep_arr(2), kep_arr(3), kep_arr(4), kep_arr(5), kep_arr(6), mu_s);
        
        % Lambert-solver: 1st call
        [a,p,e,ERROR,VVT1a,VVT2a,TPAR,theta] = lambertMR( RR1, RR2 , ToF_seconds, mu_s, 0, 0, 2 );
        VVT1a = VVT1a(:); VVT2a = VVT2a(:);

        DeltaV_1 = norm(VVT1a-VV1);
        DeltaV1_mat(i,j) = DeltaV_1;
        
            fprintf(repmat('\b',1,lastMsgNumel));
            % Note that k is 'imprecisely' set to zero
            progressString = sprintf('Progress: %.1f/100', floor( 1000 * ( (i-1)*length(TOF1_W)*length(TOF2_W) + (j-1)*length(TOF2_W) + 0) / (length(departure_W)*length(TOF1_W)*length(TOF2_W)))/10);
            fprintf(progressString); % print new progress-message
            lastMsgNumel = numel(progressString); % update elements counter of last message-to-delete
        
        
        parfor k = 1:length(TOF2_W)
            
            TOF2 = TOF2_W(k);
            t3 = t2 + TOF2; % Arrival time
            ToF_seconds = TOF2 * (24*3600); % [s]
            kep_arr = uplanet(t3, arr_ID);
            [RR3, VV3] = kep2car(kep_arr(1), kep_arr(2), kep_arr(3), kep_arr(4), kep_arr(5), kep_arr(6), mu_s);

            % Lambert-solver: 2nd call
            [a,p,e,ERROR,VVT1b,VVT2b,TPAR,theta] = lambertMR( RR2, RR3 , ToF_seconds, mu_s, 0, 0, 2 );
            VVT1b = VVT1b(:); VVT2b = VVT2b(:);
            
            
            [DeltaV_2, ~] = poweredGA(fb_ID, VV2, VVT2a, VVT1b);


            DeltaV_3 = norm(VVT2b-VV3);
            DeltaV3_mat(i,j,k) = DeltaV_3;
            
            
            DeltaV_tot = DeltaV_1 + DeltaV_2 + DeltaV_3;
            DeltaV_tot_tens(i,j,k) = DeltaV_tot;
            
     
%             fprintf(repmat('\b',1,lastMsgNumel)); % delete latest progress-message       
%             % The progress formula can be easily derived by taking into
%             % account that the current element number throughout the
%             % scrolling of the third order tensor DeltaV_tot_mat is given
%             % by: elemNo = [(i-1)*l(j)*l(k) + (j-1)*l(k) + k]
%             progressString = sprintf('Progress: %.1f/100', floor( 1000 * ( (i-1)*length(TOF1_W)*length(TOF2_W) + (j-1)*length(TOF2_W) + k) / (length(departure_W)*length(TOF1_W)*length(TOF2_W)))/10);
%             fprintf(progressString); % print new progress-message
%             lastMsgNumel = numel(progressString); % update elements counter of last message-to-delete
            
        end
                       
           
    end
        
    
end

    %clearvars -except DeltaV_tot_tens departure_W TOF1_W TOF2_W DeltaV1_mat DeltaV2_mat DeltaV3_mat
    filename = sprintf('tripleLoop (%s)', string(datetime('now')));
    save(filename, 'DeltaV_tot_tens', 'departure_W', 'TOF1_W', 'TOF2_W', 'DeltaV1_mat', 'DeltaV2_mat', 'DeltaV3_mat');


    
        case 3
    
    
%% Flyby Detail

mu_s = astroConstants(4); % Gravitational parameter of the Sun [km^3/s^2]


dep_ID = 8;  % Neptune ID [-]
fb_ID = 5;   % Jupiter ID [-]
arr_ID = 4;  % Mars ID    [-]
sun_ID = 10;

mu_d  = astroConstants(10 + dep_ID); % Gravitational parameter of departure planet [km^3/s^2]
mu_fb = astroConstants(10 + fb_ID);  % Gravitational parameter of flyby planet     [km^3/s^2]
mu_a  = astroConstants(10 + arr_ID); % Gravitational parameter of arrival planet   [km^3/s^2]



VV_minus = [ -11.7668403164062;
             -13.2841477532301;
             0.411952768117762];


VV_plus = [ -9.21171471407473;
            -0.772235992445798;
             0.643264211137162];

         
date2 = [2050, 01, 11, 15, 21, 04]';
t2 = date2mjd2000(date2);

dt = 10000;
zoomFactor = 400;
[t_SOI, r_p, h_GA, DeltaV, e_minus, e_plus] = plotPoweredGA(fb_ID, t2, VV_minus, VV_plus, dt, zoomFactor)

h = gcf;
set(h, 'Units', 'Normalized', 'OuterPosition', [.05 .11 .9 .78]);
fileType = 'png';
bitmapRes = 5000;
ppi = round(bitmapRes/15); % pixel-per-inch
printFigure(h, 'flybyDetail', fileType, bitmapRes, ppi)


        case 4


%% Porkchop plots for each separate arc (assuming uncoupled trajectories and neglecting PGA burn)

dep_min_date = [2019, 01, 01, 00, 00, 00]';   % [year, month, day, hours, minutes, seconds]
tDep_min = date2mjd2000(dep_min_date);
dep_max_date = [2050, 01, 01, 00, 00, 00]';
tDep_max = date2mjd2000(dep_max_date);
% dep_min = -T_1_JD;
%dep_max = dep_min + T_1_JD;

% TOF1 Window (Departure to Flyby)
TOF1_min = 1*T_syn_DwrtFB_JD
TOF1_max = 1.3*T_syn_DwrtFB_JD

% TOF2 Window (Flyby to Arrival)
TOF2_min = 1*T_syn_FBwrtA_JD;
TOF2_max = 3*T_syn_FBwrtA_JD;

tFB_min = tDep_min + TOF1_min;
tFB_max = tDep_max + TOF1_max;

tArr_min = tFB_min + TOF2_min;
tArr_max = tFB_max + TOF2_max;




shift_MJD2000 = datenum(2000,01,01,12,00,00);
dt_MJD = date2mjd2000([2000, 01, 11, 12, 00, 00]); % 5 days [MJD2000]

% PNG printing setup
fileType = 'png';
bitmapRes = 5000;
ppi = round(bitmapRes/10); % pixel-per-inch


% 1st Arc
[h1, h2] = porkChopPlot2DoF(dep_ID, fb_ID, tDep_min, tDep_max, tFB_min, tFB_max, dt_MJD, 1);
figure(h1)
hold on
plot3(t1_sel+shift_MJD2000, t2_sel+shift_MJD2000, 100, '*', 'MarkerSize', 16, 'LineWidth', 2, 'MarkerEdgeColor', 'r' )

figure(h2)
hold on
plot(t1_sel+shift_MJD2000, t2_sel+shift_MJD2000, '*', 'MarkerSize', 16, 'LineWidth', 2, 'MarkerEdgeColor', 'r' )
printFigure(h1, 'porkchop_dep2FB_1', fileType, bitmapRes, ppi)
printFigure(h2, 'porkchop_dep2FB_2', fileType, bitmapRes, ppi)


% 2nd Arc
[h1, h2] = porkChopPlot2DoF(fb_ID, arr_ID, tFB_min, tFB_max, ... 
                      tArr_min, tArr_max, dt_MJD, 2);       
figure(h1)
hold on
plot3(t2_sel+shift_MJD2000, t3_sel+shift_MJD2000, 100, '*', 'MarkerSize', 16, 'LineWidth', 2, 'MarkerEdgeColor', 'r' )

figure(h2)
hold on
plot(t2_sel+shift_MJD2000, t3_sel+shift_MJD2000, '*', 'MarkerSize', 16, 'LineWidth', 2,'MarkerEdgeColor', 'r' )
printFigure(h1, 'porkchop_FB2arr_1', fileType, bitmapRes, ppi)
printFigure(h2, 'porkchop_FB2arr_2', fileType, bitmapRes, ppi)


        case 5


%% Optimal Trajectory Evaluation

date1 = [2033, 05, 29, 23, 23, 04]';   % [year, month, day, hours, minutes, seconds]
t1 = date2mjd2000(date1);

date2 = [2050, 01, 11, 15, 21, 04]';
t2 = date2mjd2000(date2);
TOF1 = t2-t1;

date3 = [2054, 05, 15, 14, 35, 30]';
t3 = date2mjd2000(date3);
TOF2 = t3-t2;



[DeltaV_1, DeltaV_2, DeltaV_3, DeltaV_tot]  = interplTrajEvaluate(dep_ID, fb_ID, arr_ID, t1, TOF1, TOF2)


        case 6

%% Movie
isAnimatedPlot = 1;
isMovingCamera = 1;
isMovieRecorded = 1;
dt = 100000;

[myMovie] = plotLambertArcPlanet(dep_ID, fb_ID, arr_ID, t1_sel, t2_sel, t3_sel, dt, isAnimatedPlot, isMovingCamera, isMovieRecorded);
if isMovieRecorded
    writeMP4Movie(myMovie, 'interplanetaryTrajectory')
end


    end

end
