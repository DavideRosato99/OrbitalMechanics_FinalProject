function [ DeltaV_powered,R_p,tof, ERROR ] = poweredGA_reloaded( V_IN,V_OUT,Pid, state_p,stamp )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
V_P=state_p(4:6)';
kep_p=car2kep(state_p,astroConstants(4));

mu_sun=astroConstants(4);
mu_planet=astroConstants(Pid+10);
rp_min=astroConstants(Pid+20);

ERROR=0;



V_INF_minus=V_IN-V_P; % relative infinite velocty od s/c wrt to planet at entering
V_INF_plus=V_OUT-V_P; % relative infinite velocty od s/c wrt to planet at exit

% norm of velocty
V_norm_INF_minus=norm(V_INF_minus);
V_norm_INF_plus=norm(V_INF_plus);

% powered gravity assist: V_inf_min is different from V_inf_plus

% conservation of angular momentum fro the single hyperbola
% vp_min*rp=delta_min*v_inf_min
% vp_plus*rp=delta_plus*v_inf_plus

% geometry of hyperbola
% delta is the angle between v_inf_minus and v_inf_plus
ab=V_norm_INF_minus*V_norm_INF_plus;
cosdelta=dot(V_INF_minus,V_INF_plus)/ab;
sindelta=sqrt(1-cosdelta^2);
% sindelta=(norm(cross(V_INF_minus,V_INF_plus)))/ab;
delta = atan2(sindelta,cosdelta);



f_rp=@(rp) delta - asin(1./(1+(rp.*V_norm_INF_minus^2)./mu_planet)) - asin(1./(1+(rp.*V_norm_INF_plus^2)./mu_planet));

 options = optimoptions('fsolve','TolFun',1e-14,'Display','off');
R_p = fsolve(f_rp,rp_min,options);
%   R_p=lsqnonlin(f_rp,1e4,rp_min,[],options);


% 
if R_p<rp_min
    ERROR=1;
    R_p=NaN;
    DeltaV_powered=NaN;
    tof=NaN;
%     return;
else

%% hyperbola minus
e_minus=1+(R_p*V_norm_INF_minus^2)/mu_planet;
a_h_minus=-mu_planet/((V_norm_INF_minus)^2); % from energy equation -mu/2a=v^2/2
% Delta_minus=-a_h_minus*sqrt(e_minus^2 - 1);
Delta_minus=R_p*sqrt((e_minus+1)/(e_minus-1));
vp_minus=Delta_minus*V_norm_INF_minus/R_p;
h_minus=R_p*vp_minus;
P_minus=(h_minus)^2/mu_planet;
% theta_minus=acos((R_p/P_minus-1)/e_minus);
r_soi = norm(state_p(1:3))*(mu_planet/mu_sun)^(2/5);
theta_minus=acos(((P_minus/r_soi)-1)/e_minus);
 t1=timecalculator(a_h_minus,e_minus,0,theta_minus,mu_planet);

%% hyperbola plus
e_plus=1+(R_p*V_norm_INF_plus^2)/mu_planet;
a_h_plus=-mu_planet/((V_norm_INF_plus)^2); % from energy equation -mu/2a=v^2/2
% Delta_plus=-a_h_plus*sqrt(e_plus^2 - 1);
Delta_plus=R_p*sqrt((e_plus+1)/(e_plus-1));
vp_plus=Delta_plus*V_norm_INF_plus/R_p;
h_plus=R_p*vp_plus;
P_plus=(h_plus)^2/mu_planet;
% theta_plus=acos((R_p/P_plus-1)/e_plus);
theta_plus=acos((P_plus/r_soi-1)/e_plus);
 t2 = timecalculator(a_h_plus,e_plus,0, theta_plus,mu_planet);

%% deltav provided at pericenter
DeltaV_powered=abs(vp_plus-vp_minus);


 tof=abs(t1)+abs(t2);

end

%% orbit representation
if stamp
Vers1=V_INF_minus/norm(V_INF_minus);
Vers2=V_INF_plus/norm(V_INF_plus);

AR1 = R_p*sqrt(1+2*mu_planet/(R_p*norm(V_INF_minus)^2));
AR2 = R_p*sqrt(1+2*mu_planet/(R_p*norm(V_INF_plus)^2));


N = cross(Vers1,Vers2);
HH = cross(state_p(1:3),state_p(4:6));
H = HH/norm(HH);
RP = state_p(1:3)/norm(state_p(1:3));

R_Om = [cos(kep_p(4)) sin(kep_p(4)) 0; -sin(kep_p(4)) cos(kep_p(4)) 0; 0 0 1];
R_inc = [1 0 0; 0 cos(kep_p(3)) sin(kep_p(3)); 0 -sin(kep_p(3)) cos(kep_p(3))];
R_om = [-cos(kep_p(5)+kep_p(6)) -sin(kep_p(5)+kep_p(6)) 0; sin(kep_p(5)+kep_p(6)) -cos(kep_p(5)+kep_p(6)) 0; 0 0 1];

R = R_om*R_inc*R_Om; % rotation matrix (Sun centered frame -> Jupiter centered frame)

Vp = R*state_p(4:6)';
vp = Vp/norm(Vp);
n = R*N;
vers1 = R*Vers1;
vers2 = R*Vers2;


ARv1 = cross(vers1,n);
ARv2 = cross(vers2,n);
Peri = ARv1+ARv2; %regola del parallelogramma
peri = Peri/norm(Peri);

vperi = cross(n,peri);
vp1v = vp_minus*vperi;
vpp2v = vp_plus*vperi;
rpv = R_p*peri;
[kep1] = car2kep([rpv,vp1v],mu_planet);
[kep2] = car2kep([rpv,vpp2v],mu_planet);

figure(5)
drawPlanet('Mars',[0,0,0],5,1);
 plotOrbit(kep1,mu_planet,[-theta_minus/2+pi 0+pi]);
hold on
plotOrbit(kep2,mu_planet,[0+pi theta_plus/2+pi]);
% p3 = quiver3(0,0,0,peri(1),peri(2),peri(3),R_p,'r');
% p4 = quiver3(0,0,0,n(1),n(2),n(3),-kep1(1),'k');
% p5 = quiver3(0,0,0,vp(1),vp(2),vp(3),-kep1(1),'g');
% p6 = quiver3(0,0,0,1,0,0,-kep1(1),'k--');
% p7 = quiver3(0,0,0,ARv1(1),ARv1(2),ARv1(3),AR1,'m');
% p8 = quiver3(0,0,0,ARv2(1),ARv2(2),ARv2(3),AR2,'color',[.7 .2 .2]);
 
axis equal
% legend('First branch of hyperbola (approaching)','Second branch of hyperbola (leaving)',...
%     'Pericenter','Hyperbola plane normal unit vector','Jupiter velocity','Sun direction',...
%     'Aiming radius of first branch','Aiming radius of second branch')

end
