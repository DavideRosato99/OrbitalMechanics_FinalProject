%% Computation of groudtrack exercises

%% Exercise 1.1
a=8350; 
e=0.1976;
i=deg2rad(60); 
Om=deg2rad(270);
om=deg2rad(45); 
theta_0=deg2rad(230);
greenwich0 = deg2rad(0);
n_orb=3.25;
mu=astroConstants(13);
[rr_0, vv_0] = kep2car(a,e,i,Om,om,theta_0,mu);
[alpha, delta, lon, lat] = groundTrack(rr_0,vv_0,   greenwich0,   n_orb,   mu,1);


%% Exercise 1.2
a=26600; 
e=0.74; 
i=deg2rad(63.4);
Om=deg2rad(50); 
om=deg2rad(280); 
theta_0=deg2rad(0);
greenwich0 = deg2rad(0);
n_orb=30;
mu=astroConstants(13);
[rr_0, vv_0] = kep2car(a,e,i,Om,om,theta_0,mu);
[alpha, delta, lon, lat] = groundTrack(rr_0,vv_0,   greenwich0,   n_orb,   mu,1);


%% Exercise 1.3
a = 800 + 6371;
e = 0.0004953;
Om = deg2rad(0); 
om = deg2rad(40);  
theta_0 = 0;
greenwich0 = deg2rad(0);
mu=astroConstants(13);
n_orb=5;
i_vec = deg2rad([0, 30, 98])';
for k = 1:length(i_vec)
    i = i_vec(k);
    [rr_0, vv_0] = kep2car(a,e,i,Om,om,theta_0,mu);
    [alpha, delta, lon, lat] = groundTrack(rr_0,vv_0,   greenwich0,   n_orb,   mu, 1);
end
