% gifminima this script represent the set of minima achieved from the leg approach
% 
% PROTOTYPE: 
%   gifminima
%
% CONTRIBUTORS:
%   Alessandro Staffolani
%

if (isequal(exist('switchNew'),0) || switchNew == 0)
close all; clear all; clc
% load('workspace_minima_5_5.mat', 'DATA');
% load('workspace_minima_8_2.mat', 'DATA');
load('workspace_minima_15_5.mat', 'DATA');
end

id = find( isnan(DATA(:,1)) );
DATA(id,:)=[];
for ii = 1 : size(DATA,1)
    [Dv, ~, ~, ~, ~, vInfMin, vInfPlus, deltavp, h_p, DV] = dv_optMod2(DATA(ii,1:3));
    DATA(ii,5:7) = DV(1:3);
    DATA(ii,8) = norm(vInfMin);
    DATA(ii,9) = norm(vInfPlus);
    DATA(ii,10) = h_p;
end
%% animation
switch_plot = 'gif+plot_dv'; % 'gif' || 'gif+plot_dv'
dv_component = 1;
data_plot = 1;

check_animation = 1;
check_save_gif  = 0;
frame_factor = 0.01;
if(check_animation)
    filename = 'minima_1.gif';
    %     colordef black;
    %     h = figure('WindowState','maximized'); %'WindowState','maximized',
    h = figure('WindowState','maximized','Color','black'); %'WindowState','maximized',
    % gifinfo.DelayTime = frame_factor*(max(modv)-modv);
    gifinfo.DelayTime = frame_factor;
    delay = gifinfo.DelayTime;
    count = 0; %counter
    for ii = 1 : size(DATA,1)
        DATA_n(ii,1) = datenum(mjd20002date(DATA(ii,1)));
        DATA_n(ii,2) = datenum(mjd20002date(DATA(ii,2)));
        DATA_n(ii,3) = datenum(mjd20002date(DATA(ii,3)));
    end
    for n =  1 : size(DATA,1)  % best 10
        clf
        if strcmp(switch_plot,'gif+plot_dv')
            subplot(1,2,1)
        end
        hold on; grid on;
        [s_sys_3] = planet_representation( DATA(n,1), DATA(n,2), DATA(n,3) );
        set(gca,{'ycolor'},{'w'})
        set(gca,{'xcolor'},{'w'})
        set(gca,'color','none','TickLabelInterpreter','latex') % no background for the graph
        view(2)
        if strcmp(switch_plot,'gif+plot_dv')
            s2 = subplot(1,2,2);
            pos1 = get(s2, 'Position'); % gives the position of current sub-plot
            new_pos1 = pos1 + [-0.05,0.43,0.1,-0.4]; %[0,0.43,0.05,-0.4];
            set(s2, 'Position',new_pos1 ); % set new position of current sub-plot
            hold on; grid on;
            plot(DATA_n(:,1),DATA(:,4),'-w*','LineWidth',1,'DisplayName', '$$\Delta V$$');
            if dv_component
                plot(DATA_n(:,1),DATA(:,5),'-g*','LineWidth',1,'DisplayName', ' $$\Delta V_{1} $$');
                plot(DATA_n(:,1),DATA(:,6),'-y*','LineWidth',1,'DisplayName', ' $$\Delta V_{flyby} $$');
                plot(DATA_n(:,1),DATA(:,7),'-r*','LineWidth',1,'DisplayName', ' $$\Delta V_{2} $$');
%                 text(DATA_n(n,1),DATA(n,4)*0.95,['',num2str(DATA(n,4))],'Color','w','FontSize',10,'VerticalAlignment','top','Interpreter','latex')
%                 text(DATA_n(n,1),DATA(n,5)*0.95,['',num2str(DATA(n,5))],'Color','w','FontSize',10,'VerticalAlignment','top','Interpreter','latex')
%                 text(DATA_n(n,1),DATA(n,6)*1.05,['',num2str(DATA(n,6))],'Color','w','FontSize',10,'VerticalAlignment','bottom','Interpreter','latex')
%                 text(DATA_n(n,1),DATA(n,7)*0.95,['',num2str(DATA(n,7))],'Color','w','FontSize',10,'VerticalAlignment','top','Interpreter','latex')
            end
            plot(DATA_n(n,1),DATA(n,4),'go','LineWidth',2,'MarkerSize',10,'HandleVisibility','off');
            legend('TextColor','w','FontSize',14,'Location','northoutside','NumColumns',4,'Interpreter','latex');
            legend('boxoff')
            datetick('x', 'yyyy mmm dd','keeplimits');
            xtickangle(-20)
            set(gca,{'ycolor'},{'w'})
            set(gca,{'xcolor'},{'w'})
            set(gca,'color','none','TickLabelInterpreter','latex') % no background for the graph
            if( data_plot )
%             s2 = subplot(1,2,2);
%             pos1 = get(s2, 'Position'); % gives the position of current sub-plot
%             new_pos1 = pos1 + [-0.05,0.43,0.1,-0.4]; %[0,0.43,0.05,-0.4];
                text_plot1 = {'$$Date:$$ ';...
                    ['$$\;Departure\;\;$$ ', datestr(DATA_n(n,1), ' dd mmm yyyy HH:MM:SS')];...
                    ['$$\;Flyby\;\;\;\;\;\;\;\;\;\;$$ ', datestr(DATA_n(n,2), ' dd mmm yyyy HH:MM:SS')];...
                    ['$$\;Arrival\;\;\;\;\;\;\;$$ ', datestr(DATA_n(n,3), ' dd mmm yyyy HH:MM:SS')];...
                    '';...
                    '$$Transfer\:orbit:$$ ';...
                    ['$$\;e_{1}\;\;\;\;\;\;\;$$ ', num2str(s_sys_3.e_t1)];...
                    ['$$\;e_{2}\;\;\;\;\;\;\;$$ ', num2str(s_sys_3.e_t2)];...
                    ['$$\;ToF_{1}\;\;$$ ', num2str(DATA_n(n,2)-DATA_n(n,1)), '  [days]'];...
                    ['$$\;ToF_{2}\;\;$$ ', num2str(DATA_n(n,3)-DATA_n(n,2)), '  [days]'];...
%                     '$$Flyby:$$ ';...
%                     ['$$\;v_{\infty-}\;\;\;$$ ', num2str(DATA(n,8)), '  [km/s]'];...
%                     ['$$\;v_{\infty+}\;\;\;$$ ', num2str(DATA(n,9)), '  [km/s]'];...
%                     ['$$\;h_{pericentre}\;\;$$ ', num2str(DATA(n,10)), '  [km]'];...
%                     ['$$\;r_{SOI}\;\;$$ ', num2str(s_sys_3.r_SOI_2), '  [km]'];...
                    };
                annotation('textbox',[0.55,0.01,0.38,0.41],'String',text_plot1,'FontSize',15,'Color','w','EdgeColor','none','FitBoxToText','on','Interpreter','latex');
                text_plot2 = {'$$Flyby:$$ ';...
                    ['$$\;v_{\infty-}\;\;\;\;\;$$ ', num2str(DATA(n,8)), '  [km/s]'];...
                    ['$$\;v_{\infty+}\;\;\;\;\;$$ ', num2str(DATA(n,9)), '  [km/s]'];...
                    ['$$\;h_{pericentre}\;\;$$ ', num2str(DATA(n,10)), '  [km]'];...
                    ['$$\;r_{SOI}\;\;\;\;$$ ', num2str(s_sys_3.r_SOI_2), '  [km]'];...
                    };
                annotation('textbox',[0.8,0.01,0.38,0.41],'String',text_plot2,'FontSize',15,'Color','w','EdgeColor','none','FitBoxToText','on','Interpreter','latex');
            else
                pos1 = get(s2, 'Position'); % gives the position of current sub-plot
                %         new_pos1 = pos1 + [-0.1 0.15 0.18 -0.25];
                new_pos1 = pos1 + [0 0.2 0.05 -0.4];
                set(s2, 'Position',new_pos1 ); % set new position of current sub-plot
            end
        end
        drawnow
        if(check_save_gif)
            % Capture the plot as an image
            frame = getframe(h);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            % Write to the GIF File
            if n==1
                imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
            elseif n==10
                imwrite(imind,cm,filename,'gif','WriteMode','append', 'Loopcount',inf);
            else
                imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.00001);
            end
        end
        count = count+1;
    end
end