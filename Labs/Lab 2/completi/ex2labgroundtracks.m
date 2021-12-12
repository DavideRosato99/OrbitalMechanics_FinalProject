%% Computation of groundtrack exercises

%% Exercise 1.1
a=8350;
e=0.1976;
i=deg2rad(60); 
Om=deg2rad(270);
om=deg2rad(45); 
theta_0=deg2rad(230);
greenwich0 = deg2rad(0);
n_orb=15;
k=12;
m=1;
mu=astroConstants(13);
kep=[ a e Om om i theta_0];
[ra, dec, lon, lat] = Single_groundTrack_corretto(kep,   greenwich0,   n_orb,  mu)
[a_rep,ra, dec, lon, lat] = R_groundTrack_corretto(kep,   greenwich0,   n_orb,  mu, m, k)

%% Exercise 1.2
a=26600; 
e=0.74; 
i=deg2rad(63.4);
Om=deg2rad(50); 
om=deg2rad(280); 
theta_0=deg2rad(0);
greenwich0 = deg2rad(0);
n_orb=30;
k=2;
m=1;
mu=astroConstants(13);
kep=[ a e Om om i theta_0]
[ra, dec, lon, lat] = Single_groundTrack_corretto(kep,   greenwich0,   n_orb,  mu)
[a_rep,ra, dec, lon, lat] = R_groundTrack_corretto(kep,   greenwich0,   n_orb,  mu, m, k)


%% Exercise 1.3a
a = 7171.01;
e =0;
Om = deg2rad(0); 
om = deg2rad(40);  
theta_0 = 0;
greenwich0 = deg2rad(0);
mu=astroConstants(13);
n_orb=30;
i=deg2rad(0)
k=20;
m=2;
kep=[ a e Om om i theta_0];
[ra, dec, lon, lat] = Single_groundTrack_corretto(kep,   greenwich0,   n_orb,  mu)
[a_rep,ra, dec, lon, lat] = R_groundTrack_corretto(kep,   greenwich0,   n_orb,  mu, m, k)

%% Exercise 1.3b
a = 7171.01;
e =0;
Om = deg2rad(0); 
om = deg2rad(40);  
theta_0 = 0;
greenwich0 = deg2rad(0);
mu=astroConstants(13);
n_orb=30;
i=deg2rad(30);
k=29;
m=2;
kep=[ a e Om om i theta_0];
[ra, dec, lon, lat] = Single_groundTrack_corretto(kep,   greenwich0,   n_orb,  mu)
[a_rep,ra, dec, lon, lat] = R_groundTrack_corretto(kep,   greenwich0,   n_orb,  mu, m, k)

%% Exercise 1.3c
a = 7171.01;
e =0;
Om = deg2rad(0); 
om = deg2rad(40);  
theta_0 = 0;
greenwich0 = deg2rad(0);
mu=astroConstants(13);
n_orb=30;
i=deg2rad(98);
k=15;
m=1;
kep=[ a e Om om i theta_0];
[ra, dec, lon, lat] = Single_groundTrack_corretto(kep,   greenwich0,   n_orb,  mu)
[a_rep,ra, dec, lon, lat] = R_groundTrack_corretto(kep,   greenwich0,   n_orb,  mu, m, k)
