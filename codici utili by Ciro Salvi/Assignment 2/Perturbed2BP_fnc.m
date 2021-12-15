
function [dy] = Perturbed2BP_fnc(~,x,mu,J2,R_e,W_E,cd,am_ratio,J2_pert,drag_pert)

%SUMMARY
%
% HIGHLIGHT:
% perturbed2BP can take into consideration the perturbation due to the 
% oblaitness of Earth and airdrag, thanks to the switches J2_pert, drag_pert:
% Switch them 'on' to consider the desired perturbation
%
% INPUT:
% t[1] Time (we omit is here by indicating "~") [T]  
% x[6x1] State of the orbit (position and velocity) [ L, L/T ]
% mu[1] gravitational parameter of primary body [1/T] km^3/s^2
% J2 [1] Second Zonal Harmonic J2
% R_e [1] Earth's equatorial radius [L] km
% W_E [1] rad/s
% cd Drag [1] Coefficient
% am_ratio [1] Area/mass ration 
%
% OUTPUT:
% dy[6x1] Derivative of the state [ L/T^2, L/T^3 ]

r = norm(x(1:3));
alt = r - R_e;
density = rho(alt);
% density = rho_table(alt);


if J2_pert == 'on'
    aJ2_x = ((3*J2*mu*R_e^2)/(2*r^4))*((x(1)/r)*((5*x(3)^2/r^2)-1));
    aJ2_y = ((3*J2*mu*R_e^2)/(2*r^4))*((x(2)/r)*((5*x(3)^2/r^2)-1));
    aJ2_z = ((3*J2*mu*R_e^2)/(2*r^4))*((x(3)/r)*((5*x(3)^2/r^2)-3));
else
    aJ2_x = 0;
    aJ2_y = 0;
    aJ2_z = 0;
end

if drag_pert == 'on'
    v_rel = [x(4)+(W_E*x(2)); x(5)-(W_E*x(1)); x(6)];        % 1) get v(rel): 
    a_Drg = -0.5*cd*am_ratio*density*norm(v_rel)^2*(v_rel/norm(v_rel))*1000; % 3) get the  a_drg
else a_Drg = zeros(3,1);
end

dy = [x(4);x(5);x(6);(-(mu/(r^3))*x(1))+ aJ2_x + a_Drg(1);(-(mu/(r^3))*x(2))+ aJ2_y + a_Drg(2);(-(mu/(r^3))*x(3))+ aJ2_z + a_Drg(3)];

end

