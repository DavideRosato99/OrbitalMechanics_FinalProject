% ADD earth image
a=8350; 
e=0.1976;
i=deg2rad(60); 
Om=deg2rad(270);
om=deg2rad(45); 
theta_0=deg2rad(230);
greenwich0 = deg2rad(0);
n_orb=3.25;
mu=astroConstants(13);
[r_vect, v_vect] = kep2car(a,e,i,Om,om,theta_0,mu);

time_step = 10;
T = 2*pi*sqrt(a^3/mu);
t_vect = [0:time_step:n_orb*T];


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
hold on
plot(lon/pi*180,lat/pi*180,'.b')
axis equal
axis([-180 180 -90 +90])
xticks(-180:30:180)
yticks(-90:30:90)
xlabel('Longitude [°]')
ylabel('Latitude [°]')

