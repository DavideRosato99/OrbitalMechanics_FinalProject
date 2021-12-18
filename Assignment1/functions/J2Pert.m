function [apJ2_vect] = J2Pert(r_vect,J2,R,mu)
% j2peracc - the function gives the perturbed acceleration due to J2 
%            perturbation
%    
% INPUT:
%   r_vec        double [1x3]   Position vector                  [km]
%   J2           double [1x1]   J2 parameter                     [-]
%   R            double [1x1]   Planet radius                    [km]
%   mu           double [1x1]   Planet gravitational constant    [km^3/s^2]
% 
% OUTPUT:
%   apJ2_vect    double [3x1]   Perturbed acceleration due J2    [km/s^2]
%
% CALLED FUNCTIONS: [-] 
%
% CONTRIBUTORS:
%   Rosato Davide               10618468
%   Saba Mohammadi Yengeje      10789462
%   Spinelli Jason              10618465
%   Tagliati Alessia            10635119
%
% VERSIONS
%   2021-10-21: Release
%
% -------------------------------------------------------------------------

r = norm(r_vect);

x = r_vect(1);
y = r_vect(2);
z = r_vect(3);

c = 3/2*(J2*mu*R^2)/r^4;

ap_x = c*(x/r*(5*((z^2)/(r^2))- 1));
ap_y = c*(y/r*(5*((z^2)/(r^2))- 1));
ap_z = c*(z/r*(5*((z^2)/(r^2))- 3));

apJ2_vect = [ap_x ap_y ap_z]';

end