function data = propagate(data, settings)

%%
mu = data.constants.muE;
TmaxComp = data.propagation.TmaxComp;
if not(isnan(data.starting.OM) && isnan(data.starting.om))
    orbIn  = data.starting.orbIn;         % [-] Given orbit initial parameters
end
%%% Optimal
if isnan(data.starting.OM) && isnan(data.starting.om)
    orbIn  = data.optimal.orbIn;          % [-] Given orbit initial parameters
end
date0 = data.starting.date;
deltaSpan1 = data.propagation.deltaSpan1;
deltaSpan2 = data.propagation.deltaSpan2;
Tmax = data.propagation.Tmax;
Tperiod = data.starting.Torbit;

%% SPAN
deltaSpan = [deltaSpan1, deltaSpan2];
N = length(deltaSpan);

%% COMPUTATIONAL TIME COMPARISON
time    = zeros(N, 2);
Tcart   = cell(N, 1);
Ycart   = cell(N, 1);
Tgauss  = cell(N, 1);
Ygauss  = cell(N, 1);
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14);
[rr0, vv0] = par2car(orbIn, mu);
x0 = [rr0; vv0];
fprintf('** ORBIT PROPAGATION **\n')
fprintf("Computational costs:                        ");
str = fprintf('0%%\n');

ctr = 1;
for i = 1:N
    timespan = linspace(0, TmaxComp, deltaSpan(i));
    t1 = cputime;
    [Tcart{i}, Ycart{i}] = ode113(@ode_2bp, timespan, x0, options, mu, 'cart', datetime(date0));
    t2 = cputime;
    time(i, 1) = t2-t1;
    fprintf(repmat('\b', 1, str));
    str = fprintf('%g%%\n', ctr/(2*N)*100);
    ctr = ctr + 1;
    t1 = cputime;
    [Tgauss{i}, Ygauss{i}] = ode113(@ode_2bp, timespan, orbIn, options, mu, 'gauss', datetime(date0));
    t2 = cputime;
    time(i, 2) = t2-t1;
    fprintf(repmat('\b', 1, str));
    str = fprintf('%g%%\n', ctr/(2*N)*100);
    ctr = ctr + 1;
end


%% ORBIT PROPAGATION
tspan = linspace(0, Tmax, 100 * floor(Tmax/Tperiod));
size(tspan)
%%% GAUSS
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14);
[TGauss, YGauss] = ode113(@ode_2bp, tspan, orbIn, options, mu, 'gauss', datetime(date0));
fprintf('ciao')
perturbationsGauss = recallOdeFcn(@ode_2bp, TGauss, YGauss, mu, 'gauss', datetime(date0));
fprintf('ciao')
%%% CARTESIAN
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14);
[TCart, YCart] = ode113(@ode_2bp, tspan, x0, options, mu, 'cart', datetime(date0));
perturbationsCart = recallOdeFcn(@ode_2bp, TCart, YCart, mu, 'cart', datetime(date0));

%% PLOT
if settings.propagation.plot
    %%% COMPUTATIONAL COST
    figure('Name', 'Computational cost', 'NumberTitle', 'off');
    cart = loglog(deltaSpan(6:end),time(6:end, 1), 'o-', 'LineWidth', 1); hold on; grid on;
    gauss = loglog(deltaSpan(6:end),time(6:end, 2), 'o-', 'LineWidth', 1);
    xlabel('N$_{step}$ [-]'); ylabel('Time [s]'); title('Computational Cost')
    pC = polyfit(log(deltaSpan(6:end)), log(time(6:end, 1)'), 1);
    polyC = loglog(deltaSpan(6:end), exp(polyval(pC, log(deltaSpan(6:end)))), '--');
    pG = polyfit(log(deltaSpan(6:end)), log(time(6:end, 2)'), 1);
    polyG = loglog(deltaSpan(6:end), exp(polyval(pG, log(deltaSpan(6:end)))), '--');

    legend([cart, gauss, polyC, polyG], {'cartesian', 'gauss', ...
        strcat('n$^{', num2str(pC(2)), '}$'), strcat('n$^{', num2str(pG(2)), '}$')},...
        'interpreter', 'latex')
    
    %%% PROPAGATED ORBIT
    figure('Name', 'Propagated orbit', 'NumberTitle', 'off');
    for i = 1:size(YGauss, 1)
        [rr(i, 1:3), vv(1:3)] = par2car(YGauss(i, :), mu);
    end
    totYears = floor(TGauss(end)/(365*24*60*60));
    totDays = floor(TGauss(end)/(24*60*60));
    color = parula(totDays);
    for i = 1:totDays
        [index] = find(TGauss <= i*24*60*60 & TGauss >= (i-1)*24*60*60);
        if i >= totDays-1
            plot3(rr([index(1):index(end)], 1), rr([index(1):index(end)], 2), rr([index(1):index(end)], 3), 'color', color(i, :)); hold on
        else
            plot3(rr([index(1):index(end)+1], 1), rr([index(1):index(end)+1], 2), rr([index(1):index(end)+1], 3), 'color', color(i, :)); hold on
        end
    end
    
    colormap(color)
    hcb = colorbar; 
    caxis([0 totYears])
    set(get(hcb, 'Title'), 'String', '$\Delta$ t [years]', 'interpreter', 'latex')
    axis equal
    Re = astroConstants(23);
    [Xe, Ye, Ze] = sphere(100);
    Xe = Re*Xe; Ye = Re*Ye; Ze = Re*Ze;
    earthPNG = imread('earth.png');
    surf(Xe, Ye, Ze, 'CData', flipud(earthPNG), 'FaceColor', 'texture', 'edgecolor', 'none');
    grid on
    
    %% RETRIEVE PARAMETGERS FROM CARTESIAN 
    orbCart = zeros(size(YCart, 1), 6);
    for i = 1:size(YCart, 1)
        orbCart(i, :) = car2par(YCart(i, 1:3), YCart(i, 4:6), mu);
        [YCartGauss(i, 1:3), ~] = par2car(YGauss(i, :), mu);
    end
    
    %%% SEMI-MAJOR AXIS ERROR
    figure('Name', 'Semi-major axis error', 'NumberTitle', 'off');
    plot(TCart/(365*24*3600), abs(orbCart(:, 1) - YGauss(:, 1))); grid on;
    xlabel('Time [years]'); ylabel('$err_a$ [km]'); title('Semi-major axis error')
    
    %%% ECCENTRICITY ERROR
    figure('Name', 'Eccentricity error', 'NumberTitle', 'off');
    plot(TCart/(365*24*3600), abs(orbCart(:, 2) - YGauss(:, 2))); grid on;
    xlabel('Time [years]'); ylabel('$err_e$ [-]'); title('Eccentricity error')
    
    %%% INCLINATION ERROR
    figure('Name', 'Inclination error', 'NumberTitle', 'off');
    plot(TCart/(365*24*3600), abs(orbCart(:, 3) - YGauss(:, 3))*180/pi); grid on;
    xlabel('Time [years]'); ylabel('$err_i$ [deg]'); title('Inclination axis error')
    
    %%% RAAN ERROR
    figure('Name', '$\Omega$ error', 'NumberTitle', 'off');
    plot(TCart/(365*24*3600), abs(orbCart(:, 4) - YGauss(:, 4))*180/pi); grid on;
    xlabel('Time [years]'); ylabel('$err_{\Omega}$ [deg]'); title('$\Omega$ error')
    
    %%% ARGUMENT OF PERICENTER ERROR
    figure('Name', '$\omega$ error', 'NumberTitle', 'off');
    plot(TCart/(365*24*3600), abs(orbCart(:, 5) - YGauss(:, 5))*180/pi); grid on;
    xlabel('Time [years]'); ylabel('$err_{\omega}$ [km]'); title('$\omega$ error')
    
    %%% PERTURBATIONS ERROR
    figure('Name', 'Perturbations error', 'NumberTitle', 'off');
    plot(TCart/(365*24*3600), abs((vecnorm(perturbationsCart.aJ2+ perturbationsCart.aMOON))...
        - (vecnorm(perturbationsGauss.aJ2+ perturbationsGauss.aMOON)))); grid on;
    xlabel('Time [years]'); ylabel('$err_{acc}$ [$\frac{km}{s^2}$]'); title('$acc_{pert}$ error')
     
%    %%%THETA ERROR
%     figure('Name', 'True Anomaly error', 'NumberTitle', 'off');
%     plot(TCart/(365*24*3600), abs(orbCart(:, 6) - YGauss(:, 6))*180/pi); grid on;
%     xlabel('Time [years]'); ylabel('$err_{\theta}$ [km]'); title('$\tehta$ axis error')
%     
    %%% SEMI-MAJOR AXIS
    figure('Name', 'Semi-major axis evolution', 'NumberTitle', 'off');
    plot(TGauss/(365*24*3600), YGauss(:, 1)); hold on; grid on;
    plot(TGauss/(365*24*3600), movmean(YGauss(:, 1), 5000), 'LineWidth',2);
    xlabel('Time [years]'); ylabel('a [km]'); title('Semi-major axis evolution')
    legend('Gauss', 'Filtered', 'interpreter', 'latex')
    
    %%% ECCENTRICITY
    figure('Name', 'Eccentricity evolution', 'NumberTitle', 'off');
    plot(TGauss/(365*24*3600), YGauss(:, 2)); hold on; grid on;
    plot(TGauss/(365*24*3600), movmean(YGauss(:, 2), 5000), 'LineWidth',2);
    xlabel('Time [years]'); ylabel('e [-]'); title('Eccentricity evolution')
    legend('Gauss', 'Filtered', 'interpreter', 'latex')
    
    %%% INCLINATION
    figure('Name', 'Inclination evolution', 'NumberTitle', 'off');
    plot(TGauss/(365*24*3600), YGauss(:, 3)*180/pi); hold on; grid on;
    plot(TGauss/(365*24*3600), movmean(YGauss(:, 3)*180/pi, 5000), 'LineWidth',2);
    xlabel('Time [years]'); ylabel('i [deg]'); title('Inclination evolution')
    legend('Gauss', 'Filtered', 'interpreter', 'latex')
    
    %%% RAAN
    figure('Name', 'RAAN evolution', 'NumberTitle', 'off');
    plot(TGauss/(365*24*3600), YGauss(:, 4)*180/pi); hold on; grid on;
    plot(TGauss/(365*24*3600), movmean(YGauss(:, 4)*180/pi, 5000), 'LineWidth',2);
    xlabel('Time [years]'); ylabel('$\Omega$ [deg]'); title('RAAN evolution')
    legend('Gauss', 'Filtered', 'interpreter', 'latex')
    
    %%% ARGUMENT OF PERICENTER
    figure('Name', 'Argument of pericenter evolution', 'NumberTitle', 'off');
    plot(TGauss/(365*24*3600), YGauss(:, 5)*180/pi); hold on; grid on;
    plot(TGauss/(365*24*3600), movmean(YGauss(:, 5)*180/pi, 5000), 'LineWidth',2);
    xlabel('Time [years]'); ylabel('$\omega$ [deg]'); title('Argument of pericenter evolution')
    legend('Gauss', 'Filtered', 'interpreter', 'latex')
    
    %%% DISTANCE FROM EARTH
    figure('Name', 'Distance from Earth', 'NumberTitle', 'off');
    plot(TCart/(365*24*3600), vecnorm(YCartGauss')); grid on; hold on;
    plot(TCart/(365*24*3600), movmean(vecnorm(YCartGauss'), 5000), 'LineWidth', 2);
    xlabel('Time [years]'); ylabel('$|r|$ [km]'); title('Distance from Earth')
    legend('Gauss', 'Filtered', 'interpreter', 'latex')
    
    %%% J2 PERTURBATION
    figure('Name', 'J2 perturbation', 'NumberTitle', 'off');
    plot(TCart/(365*24*3600), vecnorm(perturbationsGauss.aJ2)); grid on; hold on;
    plot(TGauss/(365*24*3600), movmean(vecnorm(perturbationsGauss.aJ2), 5000), 'LineWidth', 2);
    xlabel('Time [years]'); ylabel('a$_{J2}$ [$\frac{km}{s^2}$]'); title('J2 acceleration')
    legend('Gauss', 'Filtered', 'interpreter', 'latex')
    
    %%% MOON PERTURBATION
    figure('Name', 'Moon Perturbation', 'NumberTitle', 'off');
    plot(TCart/(365*24*3600), vecnorm(perturbationsGauss.aMOON)); grid on; hold on;
    plot(TGauss/(365*24*3600), movmean(vecnorm(perturbationsGauss.aMOON), 5000), 'LineWidth', 2);
    xlabel('Time [years]'); ylabel('a$_{Moon}$ [$\frac{km}{s^2}$]'); title('Moon acceleration')
    legend('Gauss', 'Filtered', 'interpreter', 'latex')
%     
     %%% SEMI-MAJOR AXIS
    figure('Name', 'Semi-major axis evolution', 'NumberTitle', 'off');
    plot(TCart/(365*24*3600), orbCart(:, 1)); hold on; grid on;
    plot(TGauss/(365*24*3600), YGauss(:, 1));
    xlabel('Time [years]'); ylabel('a [km]'); title('Semi-major axis evolution')
    legend('Cartesian', 'Gauss', 'interpreter', 'latex')
    
    %%% ECCENTRICITY
    figure('Name', 'Eccentricity evolution', 'NumberTitle', 'off');
    plot(TCart/(365*24*3600), orbCart(:, 2)); hold on; grid on;
    plot(TGauss/(365*24*3600), YGauss(:, 2));
    xlabel('Time [years]'); ylabel('e [-]'); title('Eccentricity evolution')
    legend('Cartesian', 'Gauss', 'interpreter', 'latex')
    
    %%% INCLINATION
    figure('Name', 'Inclination evolution', 'NumberTitle', 'off');
    plot(TCart/(365*24*3600), orbCart(:, 3)*180/pi); hold on; grid on;
    plot(TGauss/(365*24*3600), YGauss(:, 3)*180/pi);
    xlabel('Time [years]'); ylabel('i [deg]'); title('Inclination evolution')
    legend('Cartesian', 'Gauss', 'interpreter', 'latex')
    
    %%% RAAN
    figure('Name', 'RAAN evolution', 'NumberTitle', 'off');
    plot(TCart/(365*24*3600), orbCart(:, 4)*180/pi); hold on; grid on;
    plot(TGauss/(365*24*3600), YGauss(:, 4)*180/pi);
    xlabel('Time [years]'); ylabel('$\Omega$ [deg]'); title('RAAN evolution')
    legend('Cartesian', 'Gauss', 'interpreter', 'latex')
    
    %%% ARGUMENT OF PERICENTER
    figure('Name', 'Argument of pericenter evolution', 'NumberTitle', 'off');
    plot(TCart/(365*24*3600), orbCart(:, 5)*180/pi); hold on; grid on;
    plot(TGauss/(365*24*3600), YGauss(:, 5)*180/pi);
    xlabel('Time [years]'); ylabel('$\omega$ [deg]'); title('Argument of pericenter evolution')
    legend('Cartesian', 'Gauss', 'interpreter', 'latex')
    
end
