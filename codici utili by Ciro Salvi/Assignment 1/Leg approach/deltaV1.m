function [ DV, T_1, T_2, DV1, DV2, r_1, r_2, v_1, v_2, v_1_l, v_2_l ] = deltaV1( p1, p2, t_start_j_1, t_end_j_1, t_start_j_2, t_end_j_2, s_size_1, s_size_2)
% deltaV1 evaluate the Labert Problem
%
% PROTOTYPE:
%   [ DV, T_1, T_2, DV1, DV2, r_1, r_2, v_1, v_2, v_1_l, v_2_l ] = deltaV1( p1, p2, T_1, t_start_j_2, t_end_j_2, s_size )
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
%   t_start_j_1[1]  lowest date on the departure orbit [mjd2000]
%   t_end_j_1[1]    highest date on the departure orbit [mjd2000]
%   t_start_j_2[1]  lowest date on the arrival orbit [mjd2000]
%   t_end_j_2[1]    highest date on the arrival orbit [mjd2000]
%   s_size_1[1]     size of the vector of dates on the departure orbit
%   s_size_2[1]     size of the vector of dates on the arrival orbit
%
% OUTPUT:
% 	DV[NxM]         matrix of total DV [km/s]
%   T_1[1xN]        vector of dates on the departure orbit [mjd2000]
%   T_2[1xM]        vector of dates on the arrival orbit [mjd2000]
% 	DV1[NxM]        matrix of departure DV [km/s]
% 	DV2[NxM]        matrix of arrival DV [km/s]
%   r_1[Nx3]        position vector of the first planet at each departure date [km]
%   r_2[Nx3]        position vector of the second planet at each arrival date [km]
%   v_1[Nx3]        velocity vector of the first planet at each departure date [km/s]
%   v_2[Mx3]        velocity vector of the second planet at each arrival date [km/s]
%   v_1_l[Nx3]      velocity vector at beginning of the leg at each departure date [km/s]
%   v_2_l[Mx3]      velocity vector at the end of the leg at each arrival date [km/s]
%
% CONTRIBUTORS:
%   Alessandro Staffolani
%   Fabio Spada
%   Suhailah Alkhawashke
%
% VERSIONS:
%   2021-02-11
%

if( nargin < 8 )
    s_size_2 = s_size_1;end
% Physical parameters
muS = astroConstants(4);
T_1 = linspace( t_start_j_1, t_end_j_1, s_size_1); % Date range on the departure orbit
T_2 = linspace( t_start_j_2, t_end_j_2, s_size_2); % Date range on the arrival orbit

%% Evaluate orbital parameters

for ii = 1 : length(T_1)
    [kep_1(ii,:)] = uplanet( T_1(ii), p1 ); % 1 Mercury - 2 Venus
    [ r_1(ii,:),v_1(ii,:) ] = kep2car( kep_1(ii,1), kep_1(ii,2), kep_1(ii,3), kep_1(ii,4), kep_1(ii,5), kep_1(ii,6), muS ,'rad');
end
for ii = 1 : length(T_2)
    [kep_2(ii,:)] = uplanet( T_2(ii), p2 ); %2 Venus - 5 Jupiter
    [ r_2(ii,:),v_2(ii,:) ] = kep2car( kep_2(ii,1), kep_2(ii,2), kep_2(ii,3), kep_2(ii,4), kep_2(ii,5), kep_2(ii,6), muS , 'rad');
end

% DV = zeros( length(T_1), length(T_2) );
DV1 = zeros( length(T_1), length(T_2) );
DV2 = zeros( length(T_1), length(T_2) );

for ii = 1 : length(T_1)
    for jj = 1 : length(T_2)
        if( T_2(jj) > T_1(ii) )            % check if arrival date is more than departure date
            Dt(ii,jj) = T_2(jj) - T_1(ii); % [days]
            [~,~,~,~,v_1_l(ii,jj,:),v_2_l(ii,jj,:),~,~] = lambertMR( r_1(ii,:), r_2(jj,:), Dt(ii,jj)*24*3600 , muS, 0, 0, 0, 1);
            DV1(ii,jj) = norm( [ v_1_l(ii,jj,1) v_1_l(ii,jj,2) v_1_l(ii,jj,3)] - v_1(ii,:) );
            DV2(ii,jj) = norm( v_2(jj,:) - [ v_2_l(ii,jj,1) v_2_l(ii,jj,2) v_2_l(ii,jj,3)] );
        end
    end
end
DV = DV1 + DV2;

end