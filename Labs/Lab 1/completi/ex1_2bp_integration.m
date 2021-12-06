clear
close all
clc
%%Two body problem integration
%% To evaluate perturbed solutions, uncomment lines

pert=0;%    0->no perturbation 1->perturbation

%Given data
mu=astroConstants(13);%km^3/s^2
r0=[ 26578.137, 0, 0]';%km
v0=[0; 2.221; 3.173];%km/s
s0=[r0; v0];
n_orbits=200;


%derivate values
a = 1/( 2/norm(r0) - dot(v0,v0)/mu);%semimajor axis - km
T=2*pi*sqrt((a^3)/(mu));%Orbital period

%Quality of the results settings

n_steps=10000;
tspan=linspace(0,n_orbits*T,n_steps);%Observation period sampling [s]
s_per_orbit=n_steps/n_orbits;%Samples per orbit

%ODE RESOLUTION

options=odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 ) ; % decreased tolerances 
%for accurate results in orbital mechanics characteristic dimension frame

if (pert==1)
    [T,S] = ode113( @(t,s) ode_2bp_pert2(t,s,mu) , tspan , s0 , options);
else
    [T,S] = ode113( @(t,s) ode_2bp(t,s,mu) , tspan , s0 , options);
end

%Variable pre-allocation

h=zeros(length(tspan),3);%angular momentum vector
eps=zeros(length(tspan),1);%total energy scalar
e=zeros(length(tspan),3);%eccentricity vector
edoth=zeros(length(tspan),1); %dot(e,h) scalar
norm_h=zeros(length(tspan),1); %norm(h) scalar


%%Study of integrals of motion
for i=1:length(tspan) %Saving values of parameters
    
    rv=[S(i,1) S(i,2) S(i,3)]';
    vv=[S(i,4) S(i,5) S(i,6)]';
    
    r=norm(rv);
    
    eps(i)= dot(vv,vv)/2-mu/r;
    e(i,:)=cross(vv,h(i,:))'/mu-rv/r;    
    h(i,:)=cross(rv,vv)';
    norm_h(i)=sqrt(dot(h(i,:),h(i,:)));    
    edoth(i)=dot(e(i,:),h(i,:));
end

%Overly complex visualisation

pa1=figure('Name','Panel 1')
tiledlayout(1,2)

nexttile
plot3(S(:,1),S(:,2),S(:,3)); %Unperturbed orbit
hold on
earth_sphere
axis equal


nexttile
first_orbit=[1:1:s_per_orbit+1];
last_orbit=[n_steps-s_per_orbit:1:n_steps];
plot3(S(first_orbit,1),S(first_orbit,2),S(first_orbit,3),'k');
hold on
plot3(S(last_orbit,1),S(last_orbit,2),S(last_orbit,3),'r');%Comparison between
%first and last orbit - useful for visualising perturbation effects
earth_sphere


pa2=figure('Name','Panel 2')
tiledlayout(2,2)

nexttile
plot(tspan(1:s_per_orbit*4),eps(1:s_per_orbit*4));
title('Specific Energy over 4*T')


nexttile
plot(tspan,edoth)
title('Magnitude of dot(e,h)')

if (pert==1)
 accel=zeros(length(tspan),3);
 j=1;
 while j<=length(tspan)
     accel(j,:)=szhJ2(S(j,:));
     
     j=j+1;
 end
 nexttile
 plot(tspan(1:s_per_orbit*4),accel(1:s_per_orbit*4,:))
 title('Perturbation accelerations over 4*T')

end
 
 
nexttile
plot(tspan(1:s_per_orbit*4),norm_h(1:s_per_orbit*4))
title('Magnitude and components of Specific Angular Momentum over 4xT')
hold on

for k=[1:1:3]
    plot(tspan(1:s_per_orbit*4),h(1:s_per_orbit*4,k))
end
legend ('norm(h)','hx','hy','hz')



OH=figure('Name','Orbit history')
i=1;
s_per_frame=2;%Discretisation steps per frame
while i<=length(tspan)

plot3(S(1:s_per_frame:i,1),S(1:s_per_frame:i,2),S(1:s_per_frame:i,3));
hold on
earth_sphere

axis([min(S(:,1)) max(S(:,1)) min(S(:,2)) max(S(:,2)) min(S(:,3)) max(S(:,3))])
figure(OH)

i=i+s_per_frame;
end

