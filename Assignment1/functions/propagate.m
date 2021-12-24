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
%%% GAUSS
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14);
[TGauss, YGauss] = ode113(@ode_2bp, [0 Tmax], orbIn, options, mu, 'gauss', datetime(date0));
fprintf('ciao')
perturbationsGauss = recallOdeFcn(@ode_2bp, TGauss, YGauss, mu, 'gauss', datetime(date0));
fprintf('ciao')
%%% CARTESIAN
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14);
[TCart, YCart] = ode113(@ode_2bp, [0 Tmax], x0, options, mu, 'cart', datetime(date0));
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
    end
    
    %%% SEMI-MAJOR AXIS
    figure('Name', 'Semi-major axis evolution', 'NumberTitle', 'off');
    plot(TCart/(365*24*3600), orbCart(:, 1)); grid on; hold on;
    plot(TGauss/(365*24*3600), YGauss(:, 1));
    xlabel('Time [years]'); ylabel('a [km]'); title('Semi-major axis evolution')
    legend('Cartesian', 'Gauss', 'interpreter', 'latex')
    
    %%% ECCENTRICITY
    figure('Name', 'Eccentricity evolution', 'NumberTitle', 'off');
    plot(TCart/(365*24*3600), orbCart(:, 2)); grid on; hold on;
    plot(TGauss/(365*24*3600), YGauss(:, 2));
    xlabel('Time [years]'); ylabel('e [-]'); title('Eccentricity evolution')
    legend('Cartesian', 'Gauss', 'interpreter', 'latex')
    
    %%% INCLINATION
    figure('Name', 'Inclination evolution', 'NumberTitle', 'off');
    plot(TCart/(365*24*3600), rad2deg(orbCart(:, 3))); grid on; hold on;
    plot(TGauss/(365*24*3600), rad2deg(YGauss(:, 3)));
    xlabel('Time [years]'); ylabel('i [deg]'); title('Inclination evolution')
    legend('Cartesian', 'Gauss', 'interpreter', 'latex')
    
    %%% RAAN
    figure('Name', 'RAAN evolution', 'NumberTitle', 'off');
    plot(TCart/(365*24*3600), rad2deg(orbCart(:, 4))); grid on; hold on;
    plot(TGauss/(365*24*3600), rad2deg(YGauss(:, 4)));
    xlabel('Time [years]'); ylabel('$\Omega$ [deg]'); title('RAAN evolution')
    legend('Cartesian', 'Gauss', 'interpreter', 'latex')
    
    %%% ARGUMENT OF PERICENTER
    figure('Name', 'Argument of pericenter evolution', 'NumberTitle', 'off');
    plot(TCart/(365*24*3600), rad2deg(orbCart(:, 5))); grid on; hold on;
    plot(TGauss/(365*24*3600), rad2deg(YGauss(:, 5)));
    xlabel('Time [years]'); ylabel('$\omega$ [deg]'); title('Argument of pericenter evolution')
    legend('Cartesian', 'Gauss', 'interpreter', 'latex')
    
    %%% DISTANCE FROM EARTH
    figure('Name', 'Distance from Earth', 'NumberTitle', 'off');
    plot(TCart/(365*24*3600), vecnorm(YCart')); grid on; hold on;
    plot(TCart/(365*24*3600), movmean(vecnorm(YCart'), 10000));
    xlabel('Time [years]'); ylabel('|r| [km]'); title('Distance from Earth')
    legend('Cartesian', 'Gauss', 'interpreter', 'latex')
    
    %%% J2 PERTURBATION
    figure('Name', 'J2 perturbation', 'NumberTitle', 'off');
    plot(TCart/(365*24*3600), vecnorm(perturbationsCart.aJ2)); grid on; hold on;
    plot(TGauss/(365*24*3600), vecnorm(perturbationsGauss.aJ2));
    plot(TCart/(365*24*3600), movmean(vecnorm(perturbationsCart.aJ2), 10000));
    xlabel('Time [years]'); ylabel('a$_{J2}$ [$frac{km}{s^2}$]'); title('J2 acceleration')
    legend('Cartesian', 'Gauss', 'interpreter', 'latex')
    
    %%% MOON PERTURBATION
    figure('Name', 'Moon Perturbation', 'NumberTitle', 'off');
    plot(TCart/(365*24*3600), vecnorm(perturbationsCart.aMOON)); grid on; hold on;
    plot(TGauss/(365*24*3600), vecnorm(perturbationsGauss.aMOON));
    xlabel('Time [years]'); ylabel('a$_{Moon}$ [$frac{km}{s^2}$]'); title('Moon acceleration')
    legend('Cartesian', 'Gauss', 'interpreter', 'latex')
    
end
