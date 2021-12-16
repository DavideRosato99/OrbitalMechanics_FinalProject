% gif_mission this script represent the best minimum achieved from the ga global approach
% 
% PROTOTYPE: 
%   gif_mission
%
% CONTRIBUTORS:
%   Alessandro Staffolani
%
% VERSIONS:
%   2021-02-11
%

close all; clear all; clc
DATA = [17986.5015486790  ,18732.3699104108  ,19714.9074062153,22.5596023954640]; % best minimum from ga global
for ii = 1:3
    DATA_n(ii) = datenum(mjd20002date(DATA(ii)));
end
%% Evaluation of data
T_t_1 = DATA(1); % [Julian day] T_merc
T_t_2 = DATA(2); % [Julian day] T_ven
T_t_3 = DATA(3); % [Julian day] T_jup
%planet_ID
p1 = 1;
p2 = 2;
p3 = 5;
muS = astroConstants(4);    % [km^3/s^2] Sun's gravitational parameter
R_1 = 0.4;                  % [AU] Mercury's distance from the Sun
R_2 = 0.7;                  % [AU] Venus's distance from the Sun
R_3 = 5.2;                  % [AU] Jupiter's distance from the Sun
% Physical parameter
%[s_sys_3] = solar_system_3( p1, p2, p3, R_1, R_2, R_3, T_t_1, T_t_2, T_t_3, T_t_4, r_t_1, r_t_2, r_t_3, r_t_4, v_p1, v_p21, v_p22, v_p3, Sun_fct, p1_fct, p2_fct, p3_fct, SOI_fct, tr_dotted )
 
AU = 1.495978707e8;                % [km] Astronomical Unit
% Sun
mu_S = astroConstants(4);          % [km^3/s^2] Sun's planetary constants
mean_radius_S = astroConstants(3); % [km] mean radius Sun

[kep_1] = uplanet( T_t_1, p1 ); % planet_1
[ r_t_1, v_p1 ] = kep2car( kep_1(1), kep_1(2), kep_1(3), kep_1(4), kep_1(5), kep_1(6), muS );
[kep_2] = uplanet( T_t_2, p2 );  % planet_2
[ r_t_2, v_p2 ] = kep2car( kep_2(1), kep_2(2), kep_2(3), kep_2(4), kep_2(5), kep_2(6), muS );
[kep_3] = uplanet( T_t_3, p3 );  % planet_3
[ r_t_3, v_p3 ] = kep2car( kep_3(1), kep_3(2), kep_3(3), kep_3(4), kep_3(5), kep_3(6), muS );

% Planet_1
mu_p1 = astroConstants(p1+10);     % [km^3/s^2] Planet_1's planetary constants
mr_p1 = astroConstants(p1+20);     % [km] mean radius Planet_1
%[kep_1_arr(:),~] = uplanet( T1, p1 );
a_p1 = kep_1(1);               % [km] semimajor axis Planet_1
% Planet_2
mu_p2 = astroConstants(p2+10);     % [km^3/s^2] Planet_2's planetary constants
mr_p2 = astroConstants(p2+20);     % [km] mean radius Planet_2
%[kep_2_arr(:),~] = uplanet( T_1, p2 );
a_p2 = kep_2(1);               % [km] semimajor axis Planet_2
% Planet_3
mu_p3 = astroConstants(p3+10);     % [km^3/s^2] Planet_3's planetary constants
mr_p3 = astroConstants(p3+20);     % [km] mean radius Planet_3
%[kep_3_arr(:),~] = uplanet( T_1, p3 );
a_p3 = kep_3(1);               % [km] semimajor axis Planet_3

Ratio_s1 = mean_radius_S / mr_p1;  % ratio between Planet_1's radius and Planet_2's radius
Ratio_12 = mr_p1 / mr_p2;          % ratio between Planet_1's radius and Planet_2's radius
Ratio_13 = mr_p1 / mr_p3;          % ratio between Planet_1's radius and Planet_3's radius

r_SOI_1 = norm(R_1*AU)*(mu_p1/mu_S)^(2/5); % [km] radius of the Sphere of Influence
r_SOI_2 = norm(R_2*AU)*(mu_p2/mu_S)^(2/5); % [km] radius of the Sphere of Influence
r_SOI_3 = norm(R_3*AU)*(mu_p3/mu_S)^(2/5); % [km] radius of the Sphere of Influence
%% Orbit
% Set options
options = odeset ( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

% Set transfer trajectory 1 data
Dt1 = T_t_2 - T_t_1;                  % [days] duration of the transfer orbit 1
tspan_t1 = [ 0 Dt1*3600*24 ];         % arc of ellipse of transfer orbit 1
[a_t1,~,e_t1,~,v_t_1,v_t_2,~,~] = lambertMR( r_t_1, r_t_2, Dt1*24*3600 , mu_S, 0, 0, 0, 1);
y0_t1     = [ r_t_1; v_t_1' ];        % Initial conditions of the transfer orbit 1
y0_t_dot1 = [ r_t_2; v_t_2' ];        % Final conditions of the transfer orbit 1

% Set transfer trajectory 2 data
Dt2 = T_t_3 - T_t_2;                  % [days] duration of the transfer orbit 2
tspan_t2 = [ 0 Dt2*3600*24 ];         % arc of ellipse of transfer orbit 2
[a_t2,~,e_t2,~,v_t_2,v_t_3,~,~] = lambertMR( r_t_2, r_t_3, Dt2*24*3600 , mu_S, 0, 0, 0, 1);
y0_t2     = [ r_t_2; v_t_2' ];        % Initial conditions of the transfer orbit 2
y0_t_dot2 = [ r_t_3; v_t_3' ];        % Final conditions of the transfer orbit 2

% Set planets data
T_p1 = (2*pi*sqrt(a_p1^3/mu_S));  % [sec] period of Planet_1
T_p2 = (2*pi*sqrt(a_p2^3/mu_S));  % [sec] period of Planet_2
T_p3 = (2*pi*sqrt(a_p3^3/mu_S));  % [sec] period of Planet_3
tspan_p1 = [ 0 T_p1 ];
tspan_p2 = [ 0 T_p2 ];
tspan_p3 = [ 0 T_p3 ];
y0_p1 = [ r_t_1; v_p1 ];          % Initial conditions of Planet_1
y0_p2 = [ r_t_2; v_p2 ];          % Initial conditions of Planet_2
y0_p3 = [ r_t_3; v_p3 ];          % Initial conditions of Planet_3

%evaluate where is the Planet 1 at arrival time 1
[kep_1_arr(:),~] = uplanet( T_t_2, p1 );
[ r_1_arr1(:),~ ] = kep2car( kep_1_arr(1), kep_1_arr(2), kep_1_arr(3), kep_1_arr(4), kep_1_arr(5), kep_1_arr(6), mu_S );
%evaluate where is the Planet 2 at departure time 1
[kep_2_dep(:),~] = uplanet( T_t_1, p2 );
[ r_2_dep1(:),v_2_dep1(:)] = kep2car( kep_2_dep(1), kep_2_dep(2), kep_2_dep(3), kep_2_dep(4), kep_2_dep(5), kep_2_dep(6), mu_S );
%evaluate where is the Planet 3 at departure time 1
[kep_3_dep(:),~] = uplanet( T_t_1, p3 );
[ r_3_dep1(:),v_3_dep1(:) ] = kep2car( kep_3_dep(1), kep_3_dep(2), kep_3_dep(3), kep_3_dep(4), kep_3_dep(5), kep_3_dep(6), mu_S );
%evaluate where is the Planet 3 at arrival time 1
[kep_3_arr(:),~] = uplanet( T_t_2, p3 );
[ r_3_arr1(:),~] = kep2car( kep_3_arr(1), kep_3_arr(2), kep_3_arr(3), kep_3_arr(4), kep_3_arr(5), kep_3_arr(6), mu_S );

%evaluate where is the Planet 1 at departure time 2
[kep_1_dep(:),~] = uplanet( T_t_2, p1 );
[ r_1_dep2(:),~] = kep2car( kep_1_dep(1), kep_1_dep(2), kep_1_dep(3), kep_1_dep(4), kep_1_dep(5), kep_1_dep(6), mu_S );
%evaluate where is the Planet 1 at arrival time 2
[kep_1_arr(:),~] = uplanet( T_t_3, p1 );
[ r_1_arr2(:),~] = kep2car( kep_1_arr(1), kep_1_arr(2), kep_1_arr(3), kep_1_arr(4), kep_1_arr(5), kep_1_arr(6), mu_S );
%evaluate where is the Planet 2 at arrival time 2
[kep_2_arr(:),~] = uplanet( T_t_3, p2 );
[ r_2_arr2(:),~] = kep2car( kep_2_arr(1), kep_2_arr(2), kep_2_arr(3), kep_2_arr(4), kep_2_arr(5), kep_2_arr(6), mu_S );
%evaluate where is the Planet 3 at departure time 2
[kep_3_dep(:),~] = uplanet( T_t_2, p3 );
[ r_3_dep2(:),~] = kep2car( kep_3_dep(1), kep_3_dep(2), kep_3_dep(3), kep_3_dep(4), kep_3_dep(5), kep_3_dep(6), mu_S );

% Propagation of the orbits
% Orbits of planets
[time1, Y_1] = ode45(@(t,y) ode_orbit(t,y,mu_S), tspan_p1, y0_p1, options);
[time2, Y_2] = ode45(@(t,y) ode_orbit(t,y,mu_S), tspan_p2, y0_p2, options);
[time3, Y_3] = ode45(@(t,y) ode_orbit(t,y,mu_S), tspan_p3, y0_p3, options);
% transfer orbits
[time1t, Y_t1] = ode45(@(t,y) ode_orbit(t,y,mu_S), tspan_t1, y0_t1, options);
[time2t, Y_t2] = ode45(@(t,y) ode_orbit(t,y,mu_S), tspan_t2, y0_t2, options);
% Time propagation
% Tspan_global = linspace(0, (Dt1+Dt2)*3600*24,42000);
Tspan_global = linspace(0, (T_p3/2),42000);
% Orbits of planets
[time1g, Y_1g] = ode45(@(t,y) ode_orbit(t,y,mu_S), Tspan_global, y0_p1, options);
[time2g, Y_2g] = ode45(@(t,y) ode_orbit(t,y,mu_S), Tspan_global, [ r_2_dep1'; v_2_dep1' ], options);
[time3g, Y_3g] = ode45(@(t,y) ode_orbit(t,y,mu_S), Tspan_global, [ r_3_dep1'; v_3_dep1' ], options);
% transfer orbits
[time1tg, Y_t1g] = ode45(@(t,y) ode_orbit(t,y,mu_S), Tspan_global, y0_t1, options);
nt1 = sum(time1tg<=(Dt1*3600*24));
[time2tg, Y_t2g] = ode45(@(t,y) ode_orbit(t,y,mu_S), time1tg((nt1+1):end), y0_t2, options);
nt2 = sum(time2tg<=((Dt1+Dt2)*3600*24));

%% Plot Orbit
static_plot = 0;
if static_plot
    Sun_fct = 40;
    p1_fct  = 3e4;
    p2_fct  = 1/5;
    p3_fct  = 1/(4e2);
    
    % figure;
    % figure('WindowState','maximized','Color','black'); %'WindowState','maximized',
    % [0.07 0.62 1] [1.00 0.41 0.16] [0.39 0.83 0.07] [1.00,0.82,0.41]
    hold on; grid on
    plot3 ( Y_1(:,1),Y_1(:,2),Y_1(:,3), 'Color',[0.07 0.62    1],'LineWidth',1);    % Planet_1 orbit '#0072BD',
    plot3 ( Y_2(:,1),Y_2(:,2),Y_2(:,3), 'Color',[1.00 0.41 0.16],'LineWidth',1);    % Planet_2 orbit '#D95319',
    plot3 ( Y_3(:,1),Y_3(:,2),Y_3(:,3), 'Color',[0.39 0.83 0.07],'LineWidth',1);    % Planet_3 orbit '#4DBEEE'
    plot3 ( Y_t1(:,1),Y_t1(:,2),Y_t1(:,3), 'Color','#4DBEEE','LineWidth',1); % Transfer orbit 1
    plot3 ( Y_t2(:,1),Y_t2(:,2),Y_t2(:,3), 'Color',[1.00,0.82,0.41],'LineWidth',1); % Transfer orbit 2 '#EDB120'
    
    factor_p1  = p1_fct * Sun_fct/Ratio_s1;                    % Size factor for planet_1
    Planet_plot ( 10,[ 0 0 0 ], Sun_fct );                     % Sun plot
    
    Planet_plot (  p1, r_t_1,    factor_p1 );                  % Planet_1 plot departure 1
    Planet_plot (  p1, r_1_arr1, factor_p1 );                  % Planet_1 plot arrival 1
    Planet_plot (  p1, r_1_dep2, factor_p1 );                  % Planet_1 plot departure 2
    Planet_plot (  p1, r_1_arr2, factor_p1 );                  % Planet_1 plot arrival 2
    
    Planet_plot (  p2, r_2_dep1, p2_fct*factor_p1/Ratio_12 );  % Planet_2 plot departure 1
    Planet_plot (  p2, r_t_2,    p2_fct*factor_p1/Ratio_12 );  % Planet_2 plot arrival 1
    Planet_plot (  p2, r_t_3,    p2_fct*factor_p1/Ratio_12 );  % Planet_2 plot departure 2
    Planet_plot (  p2, r_2_arr2, p2_fct*factor_p1/Ratio_12 );  % Planet_2 plot arrival 2
    
    Planet_plot (  p3, r_3_dep1, p3_fct*factor_p1/Ratio_13 );  % Planet_3 plot departure 1
    Planet_plot (  p3, r_3_arr1, p3_fct*factor_p1/Ratio_13 );  % Planet_3 plot arrival 1
    Planet_plot (  p3, r_3_dep2, p3_fct*factor_p1/Ratio_13 );  % Planet_3 plot departure 2
    Planet_plot (  p3, r_t_3,    p3_fct*factor_p1/Ratio_13 );  % Planet_3 plot arrival 2
    set(gca,'color','none','TickLabelInterpreter','latex') % no background for the graph
    
    view(35,20);
    axis equal
    % if( nargin >= 23 ) % SOI plot
    % SOI_plot ( r_SOI_2*SOI_fct*Sun_fct, 0.5, r_t_2);end
    
    
    % s_sys_3.r_SOI_1 = r_SOI_1;
    % s_sys_3.r_SOI_2 = r_SOI_2;
    % s_sys_3.r_SOI_3 = r_SOI_3;
    % s_sys_3.e_t1 = e_t1;
    % s_sys_3.e_t2 = e_t2;
end
%% animation

check_animation = 1;
check_save_gif  = 0;
frame_factor = 0.0000001;
if(check_animation)
    filename = 'test_4.gif';
    %     colordef black;
    %     h = figure('WindowState','maximized'); %'WindowState','maximized',
    h = figure('WindowState','maximized','Color','black'); %'WindowState','maximized',
    % gifinfo.DelayTime = frame_factor*(max(modv)-modv);
    gifinfo.DelayTime = frame_factor;
    delay = gifinfo.DelayTime;
    count = 0; %counter
    for n =  1 :1000: size(Y_1g,1)
        clf
        hold on
        grid on
        %
        Sun_fct = 30;
        p1_fct  = 3e4;
        p2_fct  = 1/5;
        p3_fct  = 1/(4e2);
        if( n<=nt1 )
            plot3 ( Y_t1g(1:n,1),Y_t1g(1:n,2),Y_t1g(1:n,3), 'Color','#4DBEEE','LineWidth',1); % Transfer orbit 1
        else
            plot3 ( Y_t1g(1:nt1,1),Y_t1g(1:nt1,2),Y_t1g(1:nt1,3), 'Color','#4DBEEE','LineWidth',1); % Transfer orbit 1
        end
        if( n>nt1 && n<=(nt1+nt2) )
            plot3 ( Y_t2g(1:(n-nt1),1),Y_t2g(1:(n-nt1),2),Y_t2g(1:(n-nt1),3), 'Color',[1.00,0.82,0.41],'LineWidth',1); % Transfer orbit 2 '#EDB120'
        elseif(n>(nt1+nt2))
            plot3 ( Y_t2g(1:(nt2),1),Y_t2g(1:(nt2),2),Y_t2g(1:(nt2),3), 'Color',[1.00,0.82,0.41],'LineWidth',1); % Transfer orbit 2 '#EDB120'
        end
        factor_p1  = p1_fct * Sun_fct/Ratio_s1;                    % Size factor for planet_1
        Planet_plot ( 10,[ 0 0 0 ], Sun_fct );                     % Sun plot
        
        Planet_plot (  p1, Y_1g(n,1:3),    factor_p1 );            % Planet_1
        %         Planet_plot (  p1, r_t_1,    factor_p1 );                  % Planet_1 plot departure 1
        plot3 ( r_t_1(1), r_t_1(2), r_t_1(3), 'o', 'Color',[0.07 0.62    1], 'LineWidth',3,'MarkerSize',5); % Planet_1 plot departure 1
        if( n>=nt1  )
            %             Planet_plot (  p1, r_1_arr1, factor_p1 );           % Planet_1 plot arrival 1
            plot3 ( r_1_arr1(1), r_1_arr1(2), r_1_arr1(3), 'o', 'Color',[0.07 0.62    1], 'LineWidth',3,'MarkerSize',5); % Planet_1 plot arrival 1
        end
        if( n>=(nt1+nt2)  )
            %             Planet_plot (  p1, r_1_arr2, factor_p1 );           % Planet_1 plot arrival 2
            plot3 ( r_1_arr2(1), r_1_arr2(2), r_1_arr2(3), 'o', 'Color',[0.07 0.62    1], 'LineWidth',3,'MarkerSize',5); % Planet_1 plot arrival 2
        end
        
        Planet_plot (  p2, Y_2g(n,1:3), p2_fct*factor_p1/Ratio_12 );     % Planet_2
        %         Planet_plot (  p2, r_2_dep1,    p2_fct*factor_p1/Ratio_12 );     % Planet_2 plot departure 1
        plot3 ( r_2_dep1(1), r_2_dep1(2), r_2_dep1(3), 'o', 'Color',[1.00 0.41 0.16], 'LineWidth',3,'MarkerSize',5); % Planet_2 plot departure 1
        if( n>=nt1  )
            %             Planet_plot (  p2, r_t_2,    p2_fct*factor_p1/Ratio_12 ); % Planet_2 plot departure 2
            plot3 ( r_t_2(1), r_t_2(2), r_t_2(3), 'o', 'Color',[1.00 0.41 0.16], 'LineWidth',3,'MarkerSize',5); % Planet_2 plot departure 2
        end
        if( n>=(nt1+nt2)  )
            %             Planet_plot (  p2, r_2_arr2, p2_fct*factor_p1/Ratio_12 ); % Planet_2 plot arrival 2
            plot3 ( r_2_arr2(1), r_2_arr2(2), r_2_arr2(3), 'o', 'Color',[1.00 0.41 0.16], 'LineWidth',3,'MarkerSize',5); % Planet_2 plot arrival 2
        end
        
        Planet_plot (  p3, Y_3g(n,1:3), p3_fct*factor_p1/Ratio_13 );     % Planet_3
        %         Planet_plot (  p3, r_3_dep1, p3_fct*factor_p1/Ratio_13 );        % Planet_3 plot departure 1
        plot3 ( r_3_dep1(1), r_3_dep1(2), r_3_dep1(3), 'o', 'Color',[0.39 0.83 0.07], 'LineWidth',3,'MarkerSize',5); % Planet_3 plot departure 1
        if( n>=nt1  )
            %             Planet_plot (  p3, r_3_dep2, p3_fct*factor_p1/Ratio_13 ); % Planet_3 plot departure 2
            plot3 ( r_3_dep2(1), r_3_dep2(2), r_3_dep2(3), 'o', 'Color',[0.39 0.83 0.07], 'LineWidth',3,'MarkerSize',5); % Planet_3 plot departure 2
        end
        if( n>=(nt1+nt2)  )
            %             Planet_plot (  p3, r_t_3,    p3_fct*factor_p1/Ratio_13 );  % Planet_3 plot arrival 2
            plot3 ( r_t_3(1), r_t_3(2), r_t_3(3), 'o', 'Color',[0.39 0.83 0.07], 'LineWidth',3,'MarkerSize',5); % Planet_3 plot arrival 2
        end
        
        plot3 ( Y_1g(1:n,1),Y_1g(1:n,2),Y_1g(1:n,3), 'Color',[0.07 0.62    1],'LineWidth',1);    % Planet_1 orbit '#0072BD',
        plot3 ( Y_2g(1:n,1),Y_2g(1:n,2),Y_2g(1:n,3), 'Color',[1.00 0.41 0.16],'LineWidth',1);    % Planet_2 orbit '#D95319',
        plot3 ( Y_3g(1:n,1),Y_3g(1:n,2),Y_3g(1:n,3), 'Color',[0.39 0.83 0.07],'LineWidth',1);    % Planet_3 orbit '#4DBEEE'

        %
        set(gca,{'xcolor'},{'w'})
        set(gca,{'ycolor'},{'w'})
        set(gca,{'zcolor'},{'w'})
        set(gca,'color','none','TickLabelInterpreter','latex') % no background for the graph
        %
        current_date = datenum(mjd20002date(DATA(1)+(time1g(n)/(24*3600))));
        if( n<=nt1  )
            sc_speed = norm(Y_t1g(n,4:6));
        elseif( n>nt1 && n<=(nt1+nt2) )
            sc_speed = norm(Y_t2g((n-nt1),4:6));
        else
            sc_speed = norm(Y_3g(n,4:6));
        end
        text_plot1 = {'$$Date:$$ ';...
            ['$$\;Departure\;\;$$ ', datestr(DATA_n(1), ' dd mmm yyyy HH:MM:SS')];...
            ['$$\;Flyby\;\;\;\;\;\;\;\;\;\;$$ ', datestr(DATA_n(2), ' dd mmm yyyy HH:MM:SS')];...
            ['$$\;Arrival\;\;\;\;\;\;\;$$ ', datestr(DATA_n(3), ' dd mmm yyyy HH:MM:SS')];...
%             ['$$\;Current\:date\;$$ ', datestr(current_date, ' dd mmm yyyy HH:MM:SS')];...
%             ['$$\;v_{S/C}\;\;\;\;\;$$ ', num2str(sc_speed), '  [km/s]'];...
            };
        %             ['$$\;v_{\infty+}\;\;\;\;\;$$ ', num2str(DATA(n,9)), '  [km/s]'];...
        %             ['$$\;h_{pericentre}\;\;$$ ', num2str(DATA(n,10)), '  [km]'];...
        %             ['$$\;r_{SOI}\;\;\;\;$$ ', num2str(s_sys_3.r_SOI_2), '  [km]'];...
        %             };
        annotation('textbox',[0,0.57,0.38,0.41],'String',text_plot1,'FontSize',15,'Color','w','EdgeColor','none','FitBoxToText','on','Interpreter','latex');
        text_plot2 = {['$$\;Current\:date\;\;\;$$ ', datestr(current_date, ' dd mmm yyyy HH:MM:SS'),'$$\;\;\;\;\;v_{S/C}\;\;\;\;\;$$ ', num2str(sc_speed), ' [km/s]'];...
            };
         annotation('textbox',[0,0.05,0.4,0.0],'String',text_plot2,'FontSize',15,'Color','w','EdgeColor','none','FitBoxToText','on','Interpreter','latex');
        
        if(n<(size(Y_1g,1)/2))
            view(200+(n*180/size(Y_1g,1)), 20+(n*90/size(Y_1g,1)));
        else
            view(200+(n*180/size(Y_1g,1)), 20+(45)-(n*50/size(Y_1g,1)));
        end
        axis equal
        
        drawnow
        if(check_save_gif)
            % Capture the plot as an image
            frame = getframe(h);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            % Write to the GIF File
            if n == 1
                imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
            else
                imwrite(imind,cm,filename,'gif','WriteMode','append');
            end
        end
        count = count+1;
    end
end