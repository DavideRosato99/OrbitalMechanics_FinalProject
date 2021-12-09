function [alpha, delta, lon, lat] = groundTrack(rr0,vv0,   greenwich0,   n_orb,  mu, perturbedMotion)

mu   = astroConstants(13);
om_E = (2*pi + 2*pi/365.26) / (24*3600); % Earth's rotation velocity in rad (=15.04 deg/h)/ 1h in seconds 

[a,e,i,Om,om,theta] = car2kep(rr0,vv0,  mu);


% Initial conditions
y0 = [rr0(:); vv0(:)];

% Integration timespan
time_step = 10;
T = 2*pi*sqrt(a^3/mu);
t_span = [0:time_step:n_orb*T];

% Set integration options
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );


% Perform the integration
if perturbedMotion
    [T, Y] = ode113( @(t,y) ode_2bodyPerturb(t,y, mu), t_span, y0, options );
else
    [T, Y] = ode113( @(t,y) ode_2body(t,y, mu), t_span, y0, options );
end


rr_mat = Y(:,1:3);
[alpha, delta, l, m, n] = car2RA_Dec(rr_mat);

thetaG_vect = wrapTo2Pi( greenwich0 - om_E*t_span );

lon = rad2deg( wrapTo2Pi( atan2(m,l)+ thetaG_vect' )-pi) ; % THIS IS WRONG BY 180 degrees

lat = rad2deg( delta ); % vector of latitudes



%%  PLOT

    h = figure('Name','Ground Track');
    set(h, 'Units', 'Normalized', 'OuterPosition', [.15 .25 .7 .7]);
    hold on;
    axis equal;
    set(gca,'XTick',[-180:15:180],'XTickMode','manual');
    set(gca,'YTick',[-90:10:90],'YTickMode','manual');
    xlim([-180,180]); ylim([-90,90]);
        

    image_file = 'earth.jpg';
    cdata      = flip(imread(image_file));
    imagesc([-180,180],[-90, 90],cdata);

    
    plot(lon, lat, '.g');
    xlabel('Longitude \lambda   [deg]')
    ylabel('Latitude \phi   [deg]')



end

