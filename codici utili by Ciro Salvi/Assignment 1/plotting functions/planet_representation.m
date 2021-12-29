function [s_sys_3] = planet_representation(T1,T2,T3)
% planet_representation This function call solar_system_3.m to represent
% three planet and the relative orbits
%
% PROTOTYPE:
%	[s_sys_3] = planet_representation(T1,T2,T3)
%
% INPUT:
%	  T1[1]        Mercury departure date [mjd2000]
%     T2[1]        Venus flyby date [mjd2000]
%     T3[1]        Jupiter arrival date [mjd2000]
%OUTPUT:
%   s_sys_3[struct]:
%    r_SOI_1[1]        radius of SOI of the planet_1 [km]
%    r_SOI_2[1]        radius of SOI of the planet_2 [km] 
%    r_SOI_3[1]        radius of SOI of the planet_3 [km] 
%    e_t1[1]           eccentricity of the first transfer orbit [-] 
%    e_t2[1]           eccentricity of the second transfer orbit [-] 
%
% CONTRIBUTORS:
%   Alessandro Staffolani
%
% VERSIONS:
%   2021-02-11
%

muS = astroConstants(4);    % [km^3/s^2] Sun's gravitational parameter
r_1 = 0.4;                  % [AU] Mercury's distance from the Sun
r_2 = 0.7;                  % [AU] Venus's distance from the Sun
r_3 = 5.2;                  % [AU] Jupiter's distance from the Sun
p1 = 1;
p2 = 2;
p3 = 5;
[kep_merc] = uplanet( T1, p1 ); % planet_1
[ r_merc, v_merc ] = kep2car( kep_merc(1), kep_merc(2), kep_merc(3), kep_merc(4), kep_merc(5), kep_merc(6), muS );
[kep_ven] = uplanet( T2, p2 );  % planet_2
[ r_ven, v_ven ] = kep2car( kep_ven(1), kep_ven(2), kep_ven(3), kep_ven(4), kep_ven(5), kep_ven(6), muS );
[kep_jup] = uplanet( T3, p3 );  % planet_3
[ r_jup, v_jup ] = kep2car( kep_jup(1), kep_jup(2), kep_jup(3), kep_jup(4), kep_jup(5), kep_jup(6), muS );

[s_sys_3] = solar_system_3( p1, p2, p3, r_1, r_2, r_3, T1, T2, T2, T3, r_merc, r_ven, r_ven, r_jup, v_merc, v_ven, v_ven, v_jup, 40, 3e4, 1/5, 1/(4e2),0,1);
% [s_sys_2] = solar_system_2( p1, p2, r_1, r_2, T1, T2, r_merc, r_ven, v_merc, v_ven, 20, 2e4, 1/5, 0, 1 ); 
end