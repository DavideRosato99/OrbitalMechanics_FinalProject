function [DeltaV_1, DeltaV_2, DeltaV_3, DeltaV_tot]  = interplTrajEvaluate(dep_ID, fb_ID, arr_ID, t1, TOF1, TOF2)
% PROTOTYPE:
%   [DeltaV_1, DeltaV_2, DeltaV_3, DeltaV_tot]  = interplTrajEvaluate(dep_ID, fb_ID, arr_ID, t1, TOF1, TOF2)
% 
% DESCRIPTION:
%  Returns singles and total cost of the interplanetary maneuver, given
%  planets of departure, flyby and arrival and corresponding times
% 
% INPUT:
% dep_ID, fb_ID, arr_ID: Integer number identifying the celestial body (< 10)
%                   1:   Mercury
%                   2:   Venus
%                   3:   Earth
%                   4:   Mars
%                   5:   Jupiter
%                   6:   Saturn
%                   7:   Uranus
%                   8:   Neptune
%                   9:   Pluto
%   t1, ToF1, ToF2:    time of departure, first and second time of flight
%   [Modified Julian Days 2000]
% 
% OUTPUT:
% DeltaV_1, DeltaV_2, DeltaV_3: cost of the departure maneuver,
% of the powered gravity assist and of the arrival insertion
% DeltaV_tot: total cost of the mission
% 
% CALLED FUNCTIONS:
%     astroConstants.m
%     uplanet.m
%     date2string.m
%     kep2car.m
%     lambertMR.m
%     poweredGA.m

mu_s = astroConstants(4); % Gravitational parameter of the Sun [km^3/s^2]

sun_ID = 10;

mu_fb = astroConstants(10 + fb_ID);  % Gravitational parameter of flyby planet     [km^3/s^2]


    kep_dep = uplanet(t1, dep_ID);
    [RR1, VV1] = kep2car(kep_dep(1), kep_dep(2), kep_dep(3), kep_dep(4), kep_dep(5), kep_dep(6), mu_s);
            
        t2 = t1 + TOF1; % Gravity Assist time
        ToF_seconds = TOF1 * (24*3600); % [s]
        kep_fb = uplanet(t2, fb_ID);
        [RR2, VV2] = kep2car(kep_fb(1), kep_fb(2), kep_fb(3), kep_fb(4), kep_fb(5), kep_fb(6), mu_s);
        
        % Lambert-solver: 1st call
        [a,p,e,ERROR,VVT1a,VVT2a,TPAR,theta] = lambertMR( RR1, RR2 , ToF_seconds, mu_s, 0, 0, 2 );
        VVT1a = VVT1a(:);
        VVT2a = VVT2a(:);
        
        % [a,e,i,Om,om,theta] = car2kep(RR1, VVT1a, mu_s)

        DeltaV_1 = norm(VVT1a-VV1);
        
        
            t3 = t2 + TOF2; % Arrival time
            ToF_seconds = TOF2 * (24*3600); % [s]
            kep_arr = uplanet(t3, arr_ID);
            [RR3, VV3] = kep2car(kep_arr(1), kep_arr(2), kep_arr(3), kep_arr(4), kep_arr(5), kep_arr(6), mu_s);


            % Lambert-solver: 2nd call
            [a,p,e,ERROR,VVT1b,VVT2b,TPAR,theta] = lambertMR( RR2, RR3 , ToF_seconds, mu_s, 0, 0, 2 );
            VVT1b = VVT1b(:);
            VVT2b = VVT2b(:);
            
            % [a,e,i,Om,om,theta] = car2kep(RR2, VVT1b, mu_s)
            
            
            [DeltaV_2, ~] = poweredGA(fb_ID, VV2, VVT2a, VVT1b);


            DeltaV_3 = norm(VV3-VVT2b);            
            
            DeltaV_tot = DeltaV_1 + DeltaV_2 + DeltaV_3;            
    


end

