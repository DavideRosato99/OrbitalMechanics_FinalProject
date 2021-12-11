function DeltaV_tot = fitnessFcnDeltaV(X)
% PROTOTYPE:
%   DeltaV_tot = fitnessFcnDeltaV(X)
% 
% DESCRIPTION:
%   Evaluates the costs for a given interplanetary trajectories, defined in
%   terms of departure dates and times of fligth
% 
% INPUT:
% X     [3] vector with components: departure time, ToF1, ToF2
% 
% OUTPUT:
% DeltaV_tot    total mission cost
% 
% CALLED FUNCTIONS:
%   astroConstants.m
%   uplanet.m
%   kep2car.m
%   lambertMR.m
%   poweredGA.m

t1 = X(1);
TOF1 = X(2);
TOF2 = X(3);

mu_s = astroConstants(4); % Gravitational parameter of the Sun [km^3/s^2]


dep_ID = 8;  % Neptune ID [-]
fb_ID = 5;   % Jupiter ID [-]
arr_ID = 4;  % Mars ID    [-]
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

        DeltaV_1 = norm(VVT1a-VV1);
        
        
            t3 = t2 + TOF2; % Arrival time
            ToF_seconds = TOF2 * (24*3600); % [s]
            kep_arr = uplanet(t3, arr_ID);
            [RR3, VV3] = kep2car(kep_arr(1), kep_arr(2), kep_arr(3), kep_arr(4), kep_arr(5), kep_arr(6), mu_s);


            % Lambert-solver: 2nd call
            [a,p,e,ERROR,VVT1b,VVT2b,TPAR,theta] = lambertMR( RR2, RR3 , ToF_seconds, mu_s, 0, 0, 2 );
            VVT1b = VVT1b(:);
            VVT2b = VVT2b(:);
            
            
            [DeltaV_2, ~] = poweredGA(fb_ID, VV2, VVT2a, VVT1b);


            DeltaV_3 = norm(VVT2b-VV3);            
            
            DeltaV_tot = DeltaV_1 + DeltaV_2 + DeltaV_3;            
            % DeltaV_tot = sqrt(TOF1+TOF2)*(DeltaV_1 + DeltaV_2 + DeltaV_3);


end

