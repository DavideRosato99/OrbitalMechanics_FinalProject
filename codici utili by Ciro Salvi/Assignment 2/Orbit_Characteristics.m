% 
% PROTOTYPE:
%   Orbit_Characteristics
% 
% OUTPUT
% this script evaluates all the GT for the assignment, and with the inputs
% defined in the PlanetaryMission_group_14 script, it plots the desired one
% 
% CONTRIBUTORS
% Alessandro Staffolani
% Ciro Salvi
% 
% VERSIONS
% 2020-02-11


%% --- Ground tracks ---
[alpha, delta, lat, long] = groundtrack(a,e,i,OM,om,f0,mu,w_E,tspan,fg0,time_evaluation);                                                                  %Original
[alpharep, deltarep, latrep, longrep, a_rep, ~,~, tspan_rep, T_rep ] = groundtrack(a,e,i,OM,om,f0,mu,w_E,tspan,fg0,time_evaluation,k,m);                   %Repeating GT for unperturbed
[alphaper, deltaper, latper, longper, a_per, ~,~, tspan_per, T_per ] = groundtrack_perturbed(a,e,i,OM,om,f0,mu,w_E,tspan,fg0,time_evaluation,J_2,R_e); %original GT for secular J_2
[alphaperrep, deltaperrep, latperrep, longperrep, a_perrep, Y_perrep,~, tspan_perrep, T_perrep ] =  groundtrack_perturbed(a,e,i,OM,om,f0,mu,w_E,tspan,fg0,time_evaluation,J_2,R_e,k,m);  %Repeating GT for secular J_2

lat  = rad2deg(lat);  latrep  = rad2deg(latrep);   latper  = rad2deg(latper); latperrep  = rad2deg(latperrep);
long = rad2deg(long); longrep = rad2deg(longrep);  longper = rad2deg(longper); longperrep = rad2deg(longperrep);
time.T_rep = duration(seconds(T_rep),'Format','dd:hh:mm:ss');
time.T_per = duration(seconds(T_per),'Format','dd:hh:mm:ss');

%% --- Results ---
Orbit_data = array2table( [a, a_rep, a_perrep; T_sat, T_rep, T_perrep], 'rowNames',{'Semimajor axis [km]','Period [s]'}, 'VariableNames',{'Original','Repeating_unperturbed','Repeating_secularJ2'} )
% Period     = table( time.T_sat, time.T_rep, time.T_per, 'rowNames',{'Period [hh:mm:ss]'}, 'VariableNames',{'Original','Repeating GT for unperturbed','Repeating GT for secular J_2'})
%% Plot
% switch_plot = 'perturbed'; %% 'all' || 'original' || 'repeating' || 'perturbed' || 'original_repeating' || 'repeating_perturbed' || 'original_perturbed'
figure('WindowState','maximized'); grid on; hold on;

title(['Ground track of: "',name,'"'])
A = imread('EarthTexture.jpg');
Z = image([-180 180], [90 -90],A);
uistack(Z,'bottom');
xticks([-180:30:180]);
yticks([-90:30:90]);
xticklabels({'180W', '150W', '120W', '90W', '60W', '30W', '0', '30E', '60E', '90E', '120E', '150E', '180E'})
yticklabels({'90S', '60S', '30S', '0', '30N', '60N', '90N'})
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')

switch switch_plot
    case 'all'
        plot(long,lat,'cy','LineWidth',0.5);
        plot(long(1),lat(1),'og','MarkerSize',10,'LineWidth',2);
        plot(long(end),lat(end),'xcy','MarkerSize',10,'LineWidth',2);
        
        plot(longrep,latrep,'r','LineWidth',0.5);
        plot(longrep(1),latrep(1),'og','MarkerSize',10,'LineWidth',2);
        plot(longrep(end),latrep(end),'xr','MarkerSize',10,'LineWidth',2);
        
        plot(longper,latper,'y','LineWidth',0.5);
        plot(longper(1),latper(1),'og','MarkerSize',10,'LineWidth',2);
        plot(longper(end),latper(end),'xy','MarkerSize',10,'LineWidth',2);
        
        plot(longperrep,latperrep,'b','LineWidth',0.5);
        plot(longperrep(1),latperrep(1),'og','MarkerSize',10,'LineWidth',2);
        plot(longperrep(end),latperrep(end),'xb','MarkerSize',10,'LineWidth',2);
        
        legend({ 'Original','Start', 'End',...
            'Repeating GT for unperturbed','Start repeated', 'End repeated',...
            'Original GT for secular J_2','Start perturbed', 'End perturbed',...
            'Repeating GT for secular J_2','Start perturbed repeating', 'End perturbed repeating'},...
            'Location','south','NumColumns',4);
    case 'original'
        plot(long,lat,'cy','LineWidth',0.5);
        plot(long(1),lat(1),'og','MarkerSize',10,'LineWidth',2);
        plot(long(end),lat(end),'xcy','MarkerSize',10,'LineWidth',2);
        legend({'Original','Start', 'End'},'Location','northeastoutside');
    case 'repeating'
        plot(longrep,latrep,'r','LineWidth',0.5);
        plot(longrep(1),latrep(1),'og','MarkerSize',10,'LineWidth',2);
        plot(longrep(end),latrep(end),'xr','MarkerSize',10,'LineWidth',2);
        legend({'Repeating GT for unperturbed','Start', 'End'},'Location','northeastoutside');
    case 'perturbed'
        plot(longper,latper,'y','LineWidth',0.5);
        plot(longper(1),latper(1),'og','MarkerSize',10,'LineWidth',2);
        plot(longper(end),latper(end),'xy','MarkerSize',10,'LineWidth',2);
        legend({'Original GT for secular J_2','Start', 'End'},'Location','northeastoutside');
    case 'repeating_j2'
        plot(longperrep(1:end),latperrep(1:end),'y','LineWidth',0.5);
        plot(longperrep(1),latperrep(1),'og','MarkerSize',10,'LineWidth',2);
        plot(longperrep(end),latperrep(end),'xy','MarkerSize',10,'LineWidth',2);
        legend({'Repeating GT for secular J_2','Start perturbed', 'End perturbed'},...
            'Location','northeast');
    case 'original_repeating'
        plot(long,lat,'cy','LineWidth',0.5);
        plot(long(1),lat(1),'og','MarkerSize',10,'LineWidth',2);
        plot(long(end),lat(end),'xcy','MarkerSize',10,'LineWidth',2);
        plot(longrep,latrep,'r','LineWidth',0.5);
        plot(longrep(1),latrep(1),'og','MarkerSize',10,'LineWidth',2);
        plot(longrep(end),latrep(end),'xr','MarkerSize',10,'LineWidth',2);
        legend({ 'Original','Start', 'End',...
            'Repeating GT for unperturbed','Start repeated', 'End repeated'},...
            'Location','south','NumColumns',2);        
    case 'repeating_perturbed'
        plot(longrep,latrep,'r','LineWidth',0.5);
        plot(longrep(1),latrep(1),'og','MarkerSize',10,'LineWidth',2);
        plot(longrep(end),latrep(end),'xr','MarkerSize',10,'LineWidth',2);
        plot(longper,latper,'y','LineWidth',0.5);
        plot(longper(1),latper(1),'og','MarkerSize',10,'LineWidth',2);
        plot(longper(end),latper(end),'xy','MarkerSize',10,'LineWidth',2);
        legend({'Repeating GT for unperturbed','Start repeated', 'End repeated',...
            'Repeating GT for secular J_2','Start perturbed', 'End perturbed'},...
            'Location','south','NumColumns',2);
        
    case 'original_perturbed'
        plot(long,lat,'cy','LineWidth',0.5);
        plot(long(1),lat(1),'og','MarkerSize',10,'LineWidth',2);
        plot(long(end),lat(end),'xcy','MarkerSize',10,'LineWidth',2);
        plot(longper,latper,'y','LineWidth',0.5);
        plot(longper(1),latper(1),'og','MarkerSize',10,'LineWidth',2);
        plot(longper(end),latper(end),'xy','MarkerSize',10,'LineWidth',2);
        legend({ 'Original','Start', 'End',...
            'Repeating GT for secular J_2','Start perturbed', 'End perturbed'},...
            'Location','south','NumColumns',2);
end
%% Plot Orbit
if( check_orbit )

Earth_plot
grid on; hold on;
PlotOrbit (a,e,i,OM,om,mu,0,2*pi,1e-3);     %Original
PlotOrbit (a_rep,e,i,OM,om,mu,0,2*pi,1e-3); %Repeating GT for unperturbed
PlotOrbit (a_perrep,e,i,OM,om,mu,0,2*pi,1e-3); %Repeating GT for secular J_2
hold on
axis equal
view(-35,15)
% legend({'', 'Original','Repeating GT for unperturbed','Repeating GT for secular J_2'},...
%          'Location','south','NumColumns',3);
title(['Possible orbits of: "',name,'"'])
end
