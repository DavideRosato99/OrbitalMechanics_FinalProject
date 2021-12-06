function dy = ode_2bodyPerturb( t, y, mu )


rr = y(1:3);
vv = y(4:6);

r = norm(rr);
v = norm(vv);
x = rr(1);
y = rr(2);
z = rr(3);

R_e = 6378.137;
J2 = 0.00108263;
kJ2 = 1.5*J2*mu*R_e^2/r^4;

a_J2_x = kJ2 * x/r*(5*z^2/r^2-1);
a_J2_y = kJ2 * y/r*(5*z^2/r^2-1);
a_J2_z = kJ2 * z/r*(5*z^2/r^2-3);

% Set the derivatives of the state
dy = [  vv(1)                   ;
        vv(2)                   ;
        vv(3)                   ;
        -mu/r^3 * x  +  a_J2_x  ;
        -mu/r^3 * y  +  a_J2_y  ;
        -mu/r^3 * z  +  a_J2_z  ];
    
  
% Validity check:
%  epsilon = v^2/2 - mu/r;
%  hh = cross(rr, vv);                h = norm(hh);
%  ee = 1/mu * cross(vv, hh) - rr/r;  e = norm(ee);
% 
%  [epsilon h e]'
 
end

