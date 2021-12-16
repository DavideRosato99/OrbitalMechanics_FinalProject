% Starting from a set of minima for each leg overlap_main
% evaluates the variation of cost for the first impulse and second impulse respectively
% of first and second leg shifting Venus encountering date. Selecting a
% threshold on the cost shift, some matching windows between the two leg are
% found
% 
% PROTOTYPE: 
%   overlap_main
%
% CONTRIBUTORS:
%   Fabio Spada
%   Alessandro Staffolani
%

if (isequal(exist('switchNew'),0) || switchNew == 0)
close all; clear all; clc
load('ILeg.mat');
load('IILeg.mat');
end

T1_1 = ILeg.tDepMinCost;
T1_2 = ILeg.tArrMinCost;
DV1_1 = ILeg.MinCost;

T2_1 = IILeg.tDepMinCost;
T2_2 = IILeg.tArrMinCost;
DV2_2 = IILeg.MinCost;
%%
for bb = 1 : length(T1_2) %convert arrival date
    T1_plot2(bb) = datenum(mjd20002date(T1_2(bb)));
end
for bb = 1 : length(T2_1) %convert departure date
    T2_plot1(bb) = datenum(mjd20002date(T2_1(bb)));
end

%% Range condition

range_dv_arr = 15;        % [km/s] shifting value around the minimum Delta_v
range_dv_dep = 5;         % [km/s] shifting value around the minimum Delta_v
% [8 2]  
% [15 5]  

check_arrival_range = 1;
checkPlot1 = 0;
check_departure_range = 1;
checkPlot2 = 0;

check_overlap = 1;
check_int_box = 1;
%% first leg
if( check_arrival_range )
    bound_1 = [];
    for ii = 1:length(DV1_1)
        vINorm = DV1_1(ii);
        vIVect = vINorm;
        
        tMerc = T1_1(ii);
        tVen = T1_2(ii);
        tVenLow = tVen;
        tVenUp = tVen;
        tVenVect = tVen;
        
        [kepI, ksun] = uplanet(tMerc, 1);
        [RI, VMerc]= kep2car([kepI, ksun]);
        
        while vINorm < DV1_1(ii) + range_dv_arr
            tVenLow = tVenLow - 1;                  % anticipate arrival date of one day
            TOF = (tVenLow - tMerc)*86400;
            kepF = uplanet(tVenLow, 2);
            RF = kep2car([kepF, ksun]);
            [~,~,~,~,VIH,~,~,~] = lambertMR(RI,RF,TOF,ksun,0,0,0,0);
            VIM = VIH' - VMerc;
            vINorm = norm(VIM);
            vIVect = [vINorm; vIVect];
            tVenVect = [tVenLow; tVenVect];
        end
        
        vINorm = DV1_1(ii);
        while vINorm < DV1_1(ii) + range_dv_arr
            tVenUp = tVenUp + 1;
            TOF = (tVenUp - tMerc)*86400;
            kepF = [uplanet(tVenUp, 2), ksun];
            RF = kep2car(kepF);
            [~,~,~,~,VIH,~,~,~] = lambertMR(RI,RF,TOF,ksun,0,0,0,0);
            VIM = VIH' - VMerc;
            vINorm = norm(VIM);
            vIVect = [ vIVect; vINorm];
            tVenVect = [tVenVect; tVenUp];
        end
        
        bound_1 = [bound_1; datenum(mjd20002date(tVenLow)), datenum(mjd20002date(tVenUp))];
        
        if (checkPlot1)
            for bb = 1 : length(tVenVect) %convert departure date
                tVenVect_g(bb) = datenum(mjd20002date(tVenVect(bb)));
            end
            tVen_g = datenum(mjd20002date(tVen));
            
            figure; hold on; grid minor;
            plot(tVenVect_g,vIVect, '--o', 'Color','#0072BD','LineWidth', 1, 'MarkerSize', 5 );
            plot(tVen_g, DV1_1(ii), 'r--*', 'LineWidth', 1, 'MarkerSize', 5 );
            datetick('x', 'yyyy mmm dd','keeplimits');
            xtickangle(-45);
        end
        clear tVenVect tVenVect_g vIVect
    end
end
%% second leg
if( check_departure_range )
    bound_2 = [];
    for ii = 1:length(DV2_2)
        vFNorm = DV2_2(ii);
        vFVect = vFNorm;
        
        tJup = T2_2(ii);
        tVen = T2_1(ii);
        tVenLow = tVen;
        tVenUp = tVen;
        tVenVect = tVen;
        
        [kepF, ksun] = uplanet(tJup, 5);
        [RF, VJup]= kep2car([kepF, ksun]);
        
        while vFNorm < DV2_2(ii) + range_dv_dep
            tVenLow = tVenLow - 1;
            TOF = (tJup - tVenLow)*86400;
            kepI = uplanet(tVenLow, 2);
            RI = kep2car([kepI, ksun]);
            [~,~,~,~,~,VFH,~,~] = lambertMR(RI,RF,TOF,ksun,0,0,0,0);
            VFJ = VFH' - VJup;
            vFNorm = norm(VFJ);
            vFVect = [vFNorm; vFVect];
            tVenVect = [tVenLow; tVenVect];
        end
        
        vFNorm = DV2_2(ii);
        while vFNorm < DV2_2(ii) + range_dv_dep
            tVenUp = tVenUp + 1;
            TOF = (tJup - tVenUp)*86400;
            kepI = [uplanet(tVenUp, 2), ksun];
            RI = kep2car(kepI);
            [~,~,~,~,~,VFH,~,~] = lambertMR(RI,RF,TOF,ksun,0,0,0,0);
            VFJ = VFH' - VJup;
            vFNorm = norm(VFJ);
            vFVect = [ vFVect; vFNorm];
            tVenVect = [tVenVect; tVenUp];
        end
        
        bound_2 = [bound_2; datenum(mjd20002date(tVenLow)), datenum(mjd20002date(tVenUp))];
        
        if (checkPlot2)
            for bb = 1 : length(tVenVect) %convert departure date
                tVenVect_g(bb) = datenum(mjd20002date(tVenVect(bb)));
            end
            tVen_g = datenum(mjd20002date(tVen));
            figure; hold on; grid minor;
            plot(tVenVect_g,vFVect, '--o', 'Color','#0072BD','LineWidth', 1, 'MarkerSize', 5 );
            plot(tVen_g, DV2_2(ii), 'r--*', 'LineWidth', 1, 'MarkerSize', 5 );
            datetick('x', 'yyyy mmm dd','keeplimits');
            xtickangle(-45);
        end
        clear tVenVect tVenVect_g vFVect
    end
end
%% Overlap graph

if( check_overlap )
    figure; hold on; grid minor;
    yyaxis right
    plot( T2_plot1, DV2_2, '--o', 'LineWidth', 1, 'MarkerSize', 5 ); % Leg_2
    
    yyaxis left
    plot( T1_plot2, DV1_1, '--o', 'LineWidth', 1, 'MarkerSize', 5 ); % Leg_1
    xlabel('$$Venus\: date$$','Interpreter','latex');
    ylabel('$$[km/s]$$','Interpreter','latex');
    title( '$$ Compare\: minimum \Delta V\: of\:  the\: leg $$','Interpreter','latex');
    datetick('x', 'yyyy mmm dd','keeplimits');
    xtickangle(-45)
    legend('$$ \Delta v_1 \:of\: first\: leg $$','$$ \Delta v_2\: of \:second\: leg $$','Interpreter','latex');
    yyaxis right
    if(check_departure_range)
        yl = ylim;
        for ii = 1 : size(bound_2,1)
            x_1 = bound_2(ii,1);
            x_2 = bound_2(ii,2);
            xBox = [x_1, x_1, x_2, x_2];
            yBox = [yl(1), yl(2), yl(2), yl(1)];
            p = patch(xBox, yBox, 'black','EdgeColor','none', 'FaceColor', 'r', 'FaceAlpha', 0.1);
            p.Annotation.LegendInformation.IconDisplayStyle = 'off';
        end
    end
    
    yyaxis left
    if(check_arrival_range)
        yl = ylim;
        for ii = 1 : size(bound_1,1)
            x_1 = bound_1(ii,1);
            x_2 = bound_1(ii,2);
            xBox = [x_1, x_1, x_2, x_2];
            yBox = [yl(1), yl(2), yl(2), yl(1)];
            p = patch(xBox, yBox, 'black','EdgeColor','none', 'FaceColor', 'b', 'FaceAlpha', 0.2);
            p.Annotation.LegendInformation.IconDisplayStyle = 'off';
        end
    end
    set(gca,'TickLabelInterpreter','latex')
    dcm = datacursormode;
    dcm.Enable = 'on';
    dcm.SnapToDataVertex = 'on';
    dcm.DisplayStyle = 'datatip';
    dcm.UpdateFcn = @displayCoordinates_j;
end

%% overlapping boxes extraction
% % Save the overlapping boxes in which the trasfer is possible
% % int_box = [ t_merc t_ven_1 t_ven_end t_jup ]
if( check_int_box )
    int_box = [];
    primo=0;secondo=0;terzo=0;quarto=0;
    for ii = 1 : size(bound_1,1)
        
        x_1 = bound_1(ii,1);
        x_2 = bound_1(ii,2);
        
        for jj = 1 : size(bound_2,1)
            if( x_1>=bound_2(jj,1) && x_2<=bound_2(jj,2) ) % the box of bound_1 is completely incluse in box bound_2
                primo = primo+1;
                int_box = [ int_box; datenum(mjd20002date(T1_1(ii))), x_1, x_2, datenum(mjd20002date(T2_2(jj))) ];
                %                 int_box(ii,2) = x_1;
                %                 int_box(ii,3) = x_2;
                %                 int_box(ii,1) = datenum(mjd20002date(T1_1(ii)));
                %                 int_box(ii,4) = datenum(mjd20002date(T2_2(jj)));
            elseif( x_1<=bound_2(jj,1) && x_2>=bound_2(jj,2) ) % the box of bound_2 is completely incluse in box bound_1
                secondo = secondo+1;
                int_box = [ int_box; datenum(mjd20002date(T1_1(ii))), bound_2(jj,1), bound_2(jj,2), datenum(mjd20002date(T2_2(jj))) ];
                %                 int_box(ii,2) = bound_2(jj,1);
                %                 int_box(ii,3) = bound_2(jj,2);
                %                 int_box(ii,1) = datenum(mjd20002date(T1_1(ii)));
                %                 int_box(ii,4) = datenum(mjd20002date(T2_2(jj)));
            elseif( x_1>=bound_2(jj,1) && x_2>=bound_2(jj,2) && x_1<=bound_2(jj,2) ) % the box of bound_1 is across the right side of box bound_2
                terzo = terzo+1;
                int_box = [ int_box; datenum(mjd20002date(T1_1(ii))), x_1, bound_2(jj,2), datenum(mjd20002date(T2_2(jj))) ];
                %                 int_box(ii,2) = x_1;
                %                 int_box(ii,3) = bound_2(jj,2);
                %                 int_box(ii,1) = datenum(mjd20002date(T1_1(ii)));
                %                 int_box(ii,4) = datenum(mjd20002date(T2_2(jj)));
            elseif( x_1<=bound_2(jj,1) && x_2<=bound_2(jj,2) && x_2>=bound_2(jj,1) ) % the box of bound_1 is across the left side of box bound_2
                quarto = quarto+1;
                int_box = [ int_box; datenum(mjd20002date(T1_1(ii))), bound_2(jj,1), x_2, datenum(mjd20002date(T2_2(jj))) ];
                %                 int_box(ii,2) = bound_2(jj,1);
                %                 int_box(ii,3) = x_2;
                %                 int_box(ii,1) = datenum(mjd20002date(T1_1(ii)));
                %                 int_box(ii,4) = datenum(mjd20002date(T2_2(jj)));
            end
        end
    end
    
    debug_vect=[primo,secondo,terzo,quarto];
    id = find( (int_box(:,1))==0 );
    int_box(id,:) = [];            % erase the rows no matching
end
function txt = displayCoordinates_j(~,info)
x = info.Position(1);
y = info.Position(2);
txt = {['Date = ', datestr(x, ' dd mmm yyyy HH:MM:SS')], ['Delta V = ' num2str(y) ' [km/s]']};
end