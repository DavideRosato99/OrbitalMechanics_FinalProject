function dx = ode_2bp(t,x,mu)
% ode_2bp - Ode function to compute the motion of the satellite in an
% unperturbed orbit and keep track of other objects orbiting the Earth
% through TLEs.
%
% PROTOTYPE
%   dx=ode_2bp(t,x,mu,TLEs)
%
% INPUT:
%   t        double  [1x1]   time                                 [s]
%   x        double  [6x1]   state vector               [km and km/s]
%   mu       double  [1x1]   gravitational parameter       [km^3/s^2]
%
% OUTPUT:
%   dx       double [6x1]   state vector derivative              [-]
%
% CALLED FUNCTIONS: -
%
% NOTE: time can be omitted
%
% CONTRIBUTORS:
%   Rosato Davide               10618468
%   Saba Mohammadi Yengeje      10789462
%   Spinelli jason              10618465
%   Tagliati Alessia            10635119
%
% VERSIONS
%   2021-10-21: Release
%
% -------------------------------------------------------------------------

%% COMPUTE STATE VECTOR DERIVATIVE
r = norm(x(1:3));
dx = [x(4:6); -mu/(r^3)*x(1:3)];


end