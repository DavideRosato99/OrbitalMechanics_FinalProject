function [a1, p1, e1, error1, vi1, vf1, Tpar1, theta1, a2, p2, e2, ...
    error2, vi2, vf2, Tpar2, theta2, DV1, DV2, errorFB, r_p, h_ga, delta, ...
    delta_V_powFB, e_minus, e_plus, a_minus, a_plus] = deltaVtot(...
    rr1, rr2, rr3, vv1, vv2, vv3, TOF1, TOF2, minHfl, ID)
% deltaVtot - function to compute the charcateristics of the entire flight
% from a planet 1 to a planet 3, performing a powered fly-by at planet 2.
% All input quantities are in the ecliptical reference frame.
%
% PROTOTYPE
%   [a1,p1,e1,error1,vi1,vf1,Tpar1,theta1,a2,p2,e2,error2,vi2,vf2,Tpar2,...
%   theta2,DV1,DVfb,DV2,r_p,h_ga,Delta,delta_V_powFB,e_minus,e_plus,...
%   a_minus,a_plus]=deltaVtot(rr1,rr2,rr3,vv1,vv2,vv3,depDate,tof1,tof2,...
%   minHfl,ID1,ID2,ID3)
%
% INPUT:
%   rr1        double [3x1]   Departure planet position    [km]
%   rr2        double [3x1]   Fly-By planet position       [km]
%   rr3        double [3x1]   Arrival planet position      [km]
%   vv1        double [3x1]   Departure planet velocity    [km/s]
%   vv2        double [3x1]   Fly-By planet velocity       [km/s]
%   vv3        double [3x1]   Arrival planet velocity      [km/s]
%
% OUTPUT:
%   errorFB          double  [1x1]   error in computing the radius of
%                                    pericenter (1: error)           [-]
%   r_p              double  [1x1]   perigree radius                 [km]   
%   h_ga             double  [1x1]   altitude of perigree            [km]
%   Delta            double  [1x1]   turning angle                   [deg]
%   delta_V_powFB    double  [1x1]   cost of the Fly-By manouvre     [km/s]
%   e_minus          double  [1x1]   incoming eccentricity           [-]
%   e_plus           double  [1x1]   outgoing eccentricity           [-]
%   a_minus          double  [1x1]   incoming semi-major axis        [km]
%   a_plus           double  [1x1]   outgoing semi-major axis        [km]
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

ksun = astroConstants(4);

%% DEPARTURE
[a1, p1, e1, error1, vi1, vf1, Tpar1, theta1] = lambertMR(rr1, rr2, TOF1, ksun, 0, 0, 0, 1);
DV1 = norm(vi1' - vv1);

%% ARRIVAL
[a2, p2, e2, error2, vi2, vf2, Tpar2, theta2] = lambertMR(rr2, rr3, TOF2, ksun, 0, 0, 0, 1);
DV2 = norm(vf2' - vv3);

%% FLY-BY
[errorFB, r_p, h_ga, delta, delta_V_powFB, e_minus, e_plus, a_minus,...
    a_plus] = poweredFlyby(vi1, vf2, vv2, minHfl, ID);



