function [a_rep,ra, dec, lon, lat] = Repeating_groundTrack_corretto(kep,   greenwich0,   n_orb,  mu, m, k)
%% 
% all angles must be given in radiant 
% missing the greenwich0 part
%%
p = greenwich0;

e  =  kep(2);
Om =  kep(3); 
om =  kep(4); 
i  =  kep(5);
theta_0 = kep(6);

om_E = deg2rad(15.04)/3600;  % Earth's rotation velocity in rad (=15.04 deg/h)/ 1h in seconds 
T = 2*pi/(om_E) * (m/k);
a_rep = (((T/(2*pi))^2)*mu)^(1/3);

time_step = 10;
t_vect = [0:time_step:n_orb*T];

[r_vect, v_vect] = kep2car(a_rep,e,i,Om,om,theta_0,mu);

theta= (2*pi+2*pi/365.26)/86400*t_vect;

X0 = [r_vect,v_vect];
options = odeset('Reltol',1e-13,'AbsTol',1e-14);

[T, X] = ode113( @(t,y) ode_2body(t,y, mu), t_vect, X0, options );



ra = zeros(size(t_vect));
dec = zeros(size(t_vect));
lat = zeros(size(t_vect));
lon = zeros(size(t_vect));

for j = 1:length(t_vect)
    R = [cos(theta(j)), sin(theta(j)), 0;
        -sin(theta(j)), cos(theta(j)), 0;
               0              0        1];
    r_local = R*X(j,1:3)';
    
    r = norm(X(j,1:3));
    
    dec(j) = asin(X(j,3)/r);
    
    if X(j,2) > 0
        ra(j) = acos(X(j,1)/r/cos(dec(j)));
    else
        ra(j) = 2*pi-acos(X(j,1)/r/cos(dec(j)));
    end
    
    
    r = norm(r_local);
    
    lat(j) = asin(r_local(3)/r);
    
    if r_local(2) > 0
        lon(j) = acos(r_local(1)/r/cos(lat(j)));
    else
        lon(j) = 2*pi-acos(r_local(1)/r/cos(lat(j)));
    end
    
    if lon(j) > pi
        lon(j) = lon(j) -2*pi;
    end
    
end

figure
hold on;
axis equal;
set(gca,'XTick',[-180:30:180],'XTickMode','manual');
set(gca,'YTick',[-90:30:90],'YTickMode','manual');
xlim([-180,180]); ylim([-90,90]);
        
image_file = 'earth.jpg';
cdata      = flip(imread(image_file));
imagesc([-180,180],[-90, 90],cdata);

plot(lon/pi*180,lat/pi*180,'.r');
end