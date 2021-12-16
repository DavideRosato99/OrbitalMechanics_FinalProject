function [ DV, DV1, DV2, r_1, r_2, v_1, v_2, v_1_l, v_2_l ] = deltaV0( p1, p2, t_start, t_end)
% deltaV0 evaluate the Labert Problem
%
% PROTOTYPE:
%   [ DV, DV1, DV2, r_1, r_2, v_1, v_2, v_1_l, v_2_l ] = deltaV0( p1, p2, t_start, t_end)
%
% INPUT:
%	p1[1]          ID planet 1
%	p2[1]          ID planet 2
%                     1:   Mercury
%                     2:   Venus
%                     3:   Earth
%                     4:   Mars
%                     5:   Jupiter
%                     6:   Saturn
%                     7:   Uranus
%                     8:   Neptune
%                     9:   Pluto
%                     10:  Sun
%   t_start[1]      date on the departure orbit [mjd2000]
%   t_end[1]        date on the arrival orbit [mjd2000]
%
% OUTPUT:
% 	DV[1]           total DV [km/s]
% 	DV1[1]          departure DV [km/s]
% 	DV2[1]          arrival DV [km/s]
%   r_1[1x3]        position vector of the first planet at departure date [km]
%   r_2[1x3]        position vector of the second planet at arrival date [km]
%   v_1[1x3]        velocity vector of the first planet at departure date [km/s]
%   v_2[1x3]        velocity vector of the second planet at arrival date [km/s]
%   v_1_l[1x3]      velocity vector at beginning of the leg at departure date [km/s]
%   v_2_l[1x3]      velocity vector at the end of the leg at arrival date [km/s]
%
% CONTRIBUTORS:
%   Fabio Spada
%
% VERSIONS:
%   2021-02-11
%

% Physical parameters
muS = astroConstants(4);

%% Evaluate orbital parameters

    [kep_1] = uplanet( t_start, p1 );
    [ r_1, v_1 ] = kep2car( [kep_1, muS] );


    [kep_2] = uplanet( t_end, p2 );
    [ r_2, v_2 ] = kep2car( [kep_2, muS] );


        if( t_end > t_start )            % check if arrival date is more than departure date
            Dt = t_end - t_start; % [days]
            [~,~,~,~,v_1_l,v_2_l,~,~] = lambertMR( r_1, r_2, Dt*24*3600 , muS, 0, 0, 0, 1);
            v_1_l = v_1_l'; 
            v_2_l = v_2_l'; 
            DV1 = norm( v_1_l - v_1 );
            DV2 = norm( v_2_l - v_2 );
        else 
            DV1 = NaN; 
            DV2 = NaN;
        end

DV = DV1 + DV2;

end