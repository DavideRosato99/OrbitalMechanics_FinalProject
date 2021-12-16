function [Dv, r_merc, r_ven, r_jup, v_ven, vInfMin, vInfPlus, deltavp, h_p, DV] = dv_optMod2(x)
% dv_optMod2 starting from the vector of three dates this function evaluate
% the entire interplanetary mission in terms of costs
%
% PROTOTYPE:
%   [Dv, r_merc, r_ven, r_jup, v_ven, vInfMin, vInfPlus, deltavp, h_p, DV] = dv_optMod2(x)
%
% INPUT:
%   x[3x1]              Vector of dates (initial guess):
%	  x(1)        		Mercury departure date [mjd2000]
%     x(2)        	    Venus flyby date [mjd2000]
%     x(3)        	    Jupiter arrival date [mjd2000]
%
% OUTPUT :
%   Dv[1]               Dv for the overall transfer [km/s]
%   r_merc[3x1]         Position of Mercury at the departure date from Mercury [km]
%   r_ven[3x1]          Position of Venus at the departure and arrival date on Venus [km]
%   r_jup[3x1]          Position of Jupiter at the arrival date on Jupiter [km]
%   v_venus[3x1]        Venus velocity at flyby [km/s]
%   vInfMin[3x1]        Flyby hyperbola entry velocity [km/s]
%   vInfPlus[3x1]       Flyby hyperbola exit velocity [km/s]
%   deltavp[1]          DV to be given for the flyby [km/s]
%   h_p[1]              height of perigee for the flyby [km]
%   DV[3x1]             Contibutions of Dv for manoeuvre [km/s]
%
% CONTRIBUTORS:
%   Fabio Spada
%   Alessandro Staffolani
%   Suhailah Alkhawashke
%
% VERSIONS
%   2021-02-11
%

% Integer number identifying Mercury, Venus and Jupiter for uplanet
ibody_merc = 1;
ibody_venus = 2;
ibody_jupiter = 5;

% From ephemeris compute position and velocity for the selected day

try
    
    [kep_dep_vect_merc,muS] = uplanet(x(1),ibody_merc);
    
    kep_1 = [ kep_dep_vect_merc, muS ];
    [r_merc,v_merc] = kep2car(kep_1);
    
    [kep_dep_vect_venus,~] = uplanet(x(2),ibody_venus);
    kep_2 = [ kep_dep_vect_venus, muS ];
    [r_ven,v_ven] = kep2car(kep_2);
    
    [kep_dep_vect_jupiter,~] = uplanet(x(3),ibody_jupiter);
    kep_3 = [ kep_dep_vect_jupiter, muS ];
    [r_jup,v_jup] = kep2car(kep_3);
    
    % TOF to go from Mercury to Venus, in days
    tof_1 = (x(2)-x(1))*86400;
    
    % Compute the Lambert arc for tof_1
    [~,~,~,~,VI_merc,VF_ven,~,~] = lambertMR(r_merc,r_ven,tof_1,muS);
    dv1_merc = norm(VI_merc - v_merc');
    dv2_ven = VF_ven - v_ven';
    
    vInfMin = dv2_ven;
    
    % TOF to go from Venus to Jupiter, in days
    tof_2 = (x(3)-x(2))*86400;
    
    % Compute the Lambert arc for tof_2
    [~,~,~,~,VI_ven,VF_jup,~,~] = lambertMR(r_ven,r_jup,tof_2,muS);
    dv1_ven = VI_ven - v_ven';
    dv2_jup = norm(v_jup' - VF_jup);
    
    vInfPlus = dv1_ven;
    % "Powered" DV
    [~, h_p, deltavp, ~, ~] = powFlyby( dv2_ven, dv1_ven, 2);
    
    % Compute the total DV
    Dv = dv1_merc + deltavp + dv2_jup;
    DV = [ dv1_merc, deltavp, dv2_jup ];

catch err
    
    Dv = nan;
    r_merc = nan;
    r_ven = nan;
    r_jup = nan;
    v_ven = nan;
    deltavp = nan;
    h_p = nan;
    
    if( isempty(x) )
        fprintf('x vector is empty, ga gave no input to dv_opt.\nRun again the script' );
    else
        rethrow(err)
    end
    
end


end