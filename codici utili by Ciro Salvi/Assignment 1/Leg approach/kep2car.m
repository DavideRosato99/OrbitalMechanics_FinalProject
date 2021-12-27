function [r_ECI, v_ECI] = kep2car(a,e,i,OM,om,th, mu)
% Remember to put inputs in radiants with deg2rad   

p = a*(1-e^2);
h = sqrt(p*mu);
r = p / (1+e*cos(th));

r_PF = r*[cos(th), sin(th), 0]';
v_PF = (mu/h) * [-sin(th), (e+cos(th)), 0]';

% Rotation matrices: Earth-Centered Inertial --> Perifocal   (ECI->PF)
R_om = [cos(om)  sin(om)    0   ;
        -sin(om) cos(om)    0   ;
           0        0       1   ];
       
R_i  = [   1        0       0   ;
           0     cos(i)   sin(i);
           0    -sin(i)   cos(i)];
       
R_OM = [cos(OM)  sin(OM)    0   ;
        -sin(OM) cos(OM)    0   ;
           0        0       1   ];
    
R313 = R_om * R_i * R_OM; % ECI --> PF


% PF --> ECI
r_ECI = R313'* r_PF;
v_ECI = R313'* v_PF;


end

