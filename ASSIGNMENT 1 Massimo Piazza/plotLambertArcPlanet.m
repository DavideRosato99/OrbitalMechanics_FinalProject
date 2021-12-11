function [myMovie] = plotLambertArcPlanet(dep_ID, fb_ID, arr_ID, t1, t2, t3, dt, isAnimatedPlot, isMovingCamera, isMovieRecorded)
% PROTOTYPE:
%   plotLambertArcPlanet(dep_ID, fb_ID, arr_ID, t1, t2, t3, dt,
%   isAnimatedPlot, isMovingCamera, isMovieRecorded)
% 
% DESCRIPTION:
%   Returns Plot and Movie of the interplanetary mission
% 
% INPUT:
%	dep_ID, fb_ID, arr_ID     Integer number identifying the celestial body
%	of departure, flyby and arrival (< 10)
%   t1, t2, t3     time of departure, flyby and arrival [Modified Julian
%   Days 2000] 
%   dt     time  step-size 
%   isAnimatedPlot, isMovingCamera, isMovieRecorded    Boolean parameters
% 
% OUTPUT:
%   Plot, Movie
% 
% CALLED FUNCTIONS:
%   astroConstants.m
%   lambertMR.m
%   kep2car.m
%   ode_2body.m
%   uplanet.m
%   plotPlanet.m
%   date2string.m

mu_sun = astroConstants(4);

ToF1_MJD2000 = t2 - t1; % [Julian Days]
ToF1 = ToF1_MJD2000 * (86400);

ToF2_MJD2000 = t3 - t2; % [Julian Days]
ToF2 = ToF2_MJD2000 * (86400);



h = figure;
view([30,20])
set(h, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1])
axis vis3d equal
grid on
hold on
addpath('textures');

myMovie = struct('cdata',[],'colormap',[]);
myMovie = [];


if ~isAnimatedPlot || nargin < 7
    isMovieRecorded = 0;
end

% Set integration options
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );


ToF_tot = 86400*(t3-t1);
tspan_tot = [0:dt:ToF_tot];


arc1Color = [6, 113, 201]/255;
arc2Color = [193, 54, 11]/255;

% Globe scales
sunScale_wrtSun = 100;  depScale_wrtSun = 75;  fbScale_wrtSun = 75;  arrScale_wrtSun = 75;
intermScale_wrtSun = 75;

markerSize = 6;


  kep_dep = uplanet(t1, dep_ID);
        [RR1, VV1] = kep2car(kep_dep(1), kep_dep(2), kep_dep(3), kep_dep(4), kep_dep(5), kep_dep(6), mu_sun);
  kep_fb = uplanet(t2, fb_ID);
        [RR2, VV2] = kep2car(kep_fb(1), kep_fb(2), kep_fb(3), kep_fb(4), kep_fb(5), kep_fb(6), mu_sun);
  kep_arr = uplanet(t3, arr_ID);
        [RR3, VV3] = kep2car(kep_arr(1), kep_arr(2), kep_arr(3), kep_arr(4), kep_arr(5), kep_arr(6), mu_sun);

% Compute Lambert's arcs
	% Lambert's arc 1
    [a,p,e,ERROR,VVT1,VVT2,TPAR,theta] = lambertMR( RR1, RR2 , ToF1, mu_sun, 0, 0, 2 );
    RR1 = RR1(:); RR2 = RR2(:); VVT1 = VVT1(:); VVT2 = VVT2(:);

    y0 = [RR1; VVT1];
    tspan = [0:dt:ToF1];
    [time_vec, StateMat_dep] = ode113( @(t,y) ode_2body(t,y, mu_sun), tspan, y0, options );
    X1 = StateMat_dep(:,1); Y1 = StateMat_dep(:,2); Z1 = StateMat_dep(:,3);
    
    
    % Lambert's arc 2
    [a,p,e,ERROR,VVT1fb,VVT2fb,TPAR,theta] = lambertMR( RR2, RR3 , ToF2, mu_sun, 0, 0, 2 );
    RR2 = RR2(:); RR3 = RR3(:); VVT1fb = VVT1fb(:); VVT2fb = VVT2fb(:);

    y0 = [RR2; VVT1fb];
    tspan = [0:dt:ToF2];
    [time_vec_fb, StateMat_fb] = ode113( @(t,y) ode_2body(t,y, mu_sun), tspan, y0, options );
    X2 = StateMat_fb(:,1); Y2 = StateMat_fb(:,2); Z2 = StateMat_fb(:,3);



% Initial orbit
kep0 = uplanet(t1, dep_ID);
[rr0_dep, vv0_dep] = kep2car(kep0(1), kep0(2), kep0(3), kep0(4), kep0(5), kep0(6), mu_sun);
y0 = [rr0_dep; vv0_dep];
T = 2*pi*sqrt(kep0(1)^3/mu_sun);
if isAnimatedPlot & (T < ToF_tot)
    tspan = tspan_tot;
else
    tspan = [0:dt:T];
end
[time_vec_dep, StateMat_dep] = ode113( @(t,y) ode_2body(t,y, mu_sun), tspan, y0, options );
X_dep = StateMat_dep(:,1); Y_dep = StateMat_dep(:,2); Z_dep = StateMat_dep(:,3);
plot3(X_dep,Y_dep,Z_dep, 'LineStyle', '-.', 'LineWidth', 1,  'Color', [20, 79, 221]/255)

% Intermediate orbits
if dep_ID > arr_ID
    interm_IDs = [dep_ID-1:-1:arr_ID+1];
else
    interm_IDs = [dep_ID+1:+1:arr_ID-1];
end

interm_IDs(interm_IDs==fb_ID) = []; % Remove the flyby planet in between

for k = interm_IDs
    
    j = k - min(interm_IDs) + 1;
    % Plot intermediate orbits
    kep_interm = uplanet(t1, k);
    [rr0_interm, vv0_interm] = kep2car(kep_interm(1), kep_interm(2), kep_interm(3), kep_interm(4), kep_interm(5), kep_interm(6), mu_sun);
    y0 = [rr0_interm; vv0_interm];
    T = 2*pi*sqrt(kep0(1)^3/mu_sun);
    if isAnimatedPlot & (T < ToF_tot)        
        tspan = tspan_tot;
    else
        tspan = [0:dt:T];
    end    
    [time_vec_interm, StateTens_interm(:,:,j)] = ode113( @(t,y) ode_2body(t,y, mu_sun), tspan, y0, options );
    X_interm = StateTens_interm(:,1,j); Y_interm = StateTens_interm(:,2,j); Z_interm = StateTens_interm(:,3,j);
    hVisibilityTrue = 'off';
    if k == interm_IDs(1);
        hVisibilityTrue = 'on';
    end
    plot3(X_interm,Y_interm,Z_interm,  '--k', 'LineWidth', .5, 'Color', [59 59 59]/255, 'HandleVisibility', hVisibilityTrue)
    
end


% Flyby-planet orbit
kep_fb = uplanet(t2, fb_ID);
kep0 = uplanet(t1, fb_ID);
[rr0_fb, vv0_fb] = kep2car(kep0(1), kep0(2), kep0(3), kep0(4), kep0(5), kep0(6), mu_sun);
y0 = [rr0_fb; vv0_fb];
T = 2*pi*sqrt(kep0(1)^3/mu_sun);
if isAnimatedPlot & (T < ToF_tot)
    tspan = tspan_tot;
else
    tspan = [0:dt:T];
end
[time_vec_fb, StateMat_fb] = ode113( @(t,y) ode_2body(t,y, mu_sun), tspan, y0, options );
X_fb = StateMat_fb(:,1); Y_fb = StateMat_fb(:,2); Z_fb = StateMat_fb(:,3);
plot3(X_fb,Y_fb,Z_fb, 'LineStyle', '-.', 'LineWidth', 1, 'Color', [60, 178, 70]/255)

% Final orbit
kep_arr = uplanet(t3, arr_ID);
kep0 = uplanet(t1, arr_ID);
[rr0_arr, vv0_arr] = kep2car(kep0(1), kep0(2), kep0(3), kep0(4), kep0(5), kep0(6), mu_sun);
y0 = [rr0_arr; vv0_arr];
T = 2*pi*sqrt(kep0(1)^3/mu_sun);
if isAnimatedPlot & (T < ToF_tot)
    tspan = tspan_tot;
else
    tspan = [0:dt:T];
end
[time_vec_arr, StateMat_arr] = ode113( @(t,y) ode_2body(t,y, mu_sun), tspan, y0, options );
X_arr = StateMat_arr(:,1); Y_arr = StateMat_arr(:,2); Z_arr = StateMat_arr(:,3);
plot3(X_arr,Y_arr,Z_arr, 'LineStyle', '-.', 'LineWidth', 1, 'Color', [206, 61, 21]/255)


light('Position',[0 0 0],'Style','local');
sunGlobe = plotPlanet(10, [0,0,0], h, sunScale_wrtSun);

legendPosition = [0.7 .8 0.05 0.05];
legendFontSize = 16;
if ~isAnimatedPlot || nargin < 7
    depGlobe = plotPlanet(dep_ID, RR1, h, depScale_wrtSun);
    fbGlobe  = plotPlanet(fb_ID,  RR2, h, fbScale_wrtSun);
    arrGlobe = plotPlanet(arr_ID, RR3, h, arrScale_wrtSun);
    

    plot3(X1,Y1,Z1, 'LineWidth', 3, 'Color', arc1Color)
    plot3(X2,Y2,Z2, 'LineWidth', 3, 'Color', arc2Color)
    
    h1_L = legend({'Initial Orbit', 'Intermediate-planets Orbits', 'Flyby-planet Orbit', 'Target Orbit', '1st Lambert''s Arc', '2nd Lambert''s Arc'}, 'Interpreter', 'latex', 'FontSize', legendFontSize);
    set(h1_L, 'Position', legendPosition, 'Units', 'normalized');
    
    
else
    h1_L = legend({'Initial Orbit', 'Intermediate-planets Orbits', 'Flyby-planet Orbit', 'Target Orbit'}, 'Interpreter', 'latex', 'FontSize', legendFontSize);
    set(h1_L, 'Position', legendPosition, 'Units', 'normalized');
end


camzoom(1.5)


%% Dynamic plot

X = X1; Y = Y1;   Z = Z1;
trajColor = arc1Color;
fbIndex = 0;
dynStep = 1000;
if isMovieRecorded
    dynStep = 30;
end

bitmapRes = 5000;
ppi = bitmapRes/20;

if isAnimatedPlot
    
    % Set up camera live-zooming
    R1 = norm([X_dep(1) Y_dep(1) Z_dep(1)]);
    R2 = norm([X_fb(end) Y_fb(end) Z_fb(end)]);
    zoomFact_tot = 0.45*R1/R2;
    nFrames = (length(tspan_tot)/dynStep);
    zoomFact_step = zoomFact_tot^(1/nFrames);
    
    
    % Set up camera live-rotation
    deltaAzimuth = -120; %[deg]
    deltaElevation = 10;   %[deg]
    rotationAxis = [0 0 1];

    if ~isMovingCamera % i.e. fixed camera view
        zoomFact_step = 1;
        deltaAzimuth = 0;
        deltaElevation = 0;
    end
	dAz = deltaAzimuth/nFrames;
    dEl = deltaElevation/nFrames;



    % Time-box
    timeBoxLocation = [.2 .7 .1 .1];
    timeBox = annotation('textbox',timeBoxLocation,'String', date2string(mjd20002date(t1), 'date'),  'Interpreter', 'latex');
    timeBox.FontSize = 24;
    timeBox.BackgroundColor = 'w';
   
    for k = 1:dynStep:length(tspan_tot)
    
        if k > 1
           delete(depGlobe);  delete(fbGlobe);  delete(arrGlobe);   delete(SC);   delete(traj_SC);
           delete( intermGlobe );
        end
        

        depGlobe = plotPlanet(dep_ID, StateMat_dep(k,1:3), h, depScale_wrtSun);
        fbGlobe  = plotPlanet(fb_ID,  StateMat_fb(k,1:3),  h, fbScale_wrtSun );
        arrGlobe = plotPlanet(arr_ID, StateMat_arr(k,1:3), h, arrScale_wrtSun);
        for r = interm_IDs
            s = r - min(interm_IDs) + 1;
            intermGlobe(s) = plotPlanet(r, StateTens_interm(k,1:3,s), h, intermScale_wrtSun);
        end

        if tspan_tot(k) >= ToF1
            if  ~fbIndex
                fbIndex = k;
                delete(fbGlobe);
                fbGlobe  = plotPlanet(fb_ID,  RR2, h, fbScale_wrtSun);
                delete(traj_SC);
                plot3(X,Y,Z, 'LineWidth', 3, 'Color', arc1Color);
                legend('Initial Orbit', 'Intermediate-planets Orbits', 'Flyby-planet Orbit', 'Target Orbit', '1st Lambert''s Arc')
                X = X2; Y = Y2; Z = Z2;
                trajColor = arc2Color;
                
                printFigure(h, 'heliocConfig_FB', 'png', bitmapRes, ppi) % Print frame @FB
            end
        end
        
        
        j = k - fbIndex + 1;

        date_MJD2000 = tspan_tot(k)/86400 + t1;
        
        
        if (length(tspan_tot)-k) < dynStep % last frame
            %delete(depGlobe); delete(arrGlobe); delete(SC); 
            delete(traj_SC);
            plot3(X,Y,Z, 'LineWidth', 3, 'Color', arc2Color);
            j = 1; X = RR3(1); Y = RR3(2); Z = RR3(3);
            date_MJD2000 = tspan_tot(end)/86400 + t1;
            legend('Initial Orbit', 'Intermediate-planets Orbits', 'Flyby-planet Orbit', 'Target Orbit', '1st Lambert''s Arc', '2nd Lambert''s Arc')
            printFigure(h, 'heliocConfig_arr', 'png', bitmapRes, ppi) % Print frame @arr
        end
        SC = plot3(X(j), Y(j), Z(j), 'o', 'HandleVisibility', 'off',...    % Plot the current position of the S/C
                            'MarkerSize', markerSize,'MarkerEdgeColor','r','MarkerFaceColor',[0.8,0.2,0.2]);
                        
        if (length(tspan_tot)-k) >= dynStep
            traj_SC = plot3(X(1:j), Y(1:j), Z(1:j), 'HandleVisibility', 'off', 'LineWidth', 3, 'Color', trajColor);
        end
        
        


        timeBox.String = sprintf('%s', date2string(mjd20002date(date_MJD2000), 'date'));   
        camzoom(zoomFact_step); % zoom camera
        camorbit(dAz, dEl,'data', rotationAxis) % rotate camera
        drawnow
        if isMovieRecorded
            myMovie = [myMovie getframe(h)];
        end
        
        if k == 1 % Print frame @dep (i.e. 1st frame after having properly plotted)
            printFigure(h, 'heliocConfig_dep', 'png', bitmapRes, ppi)
        end
        

    end
     
    
end




end
