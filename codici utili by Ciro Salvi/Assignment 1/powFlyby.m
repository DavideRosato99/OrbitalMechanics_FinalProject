function [rpSol, hpSol, deltavp, deltaVect, eVect] = powFlyby(vInfMin, vInfPlus, idAtt)
% powFlyby evaluates the powered flyby characteristic required for a given boundary
% problem around a given planet; if the solution passes at a height smaller than the 
% radius of the planet plus its atmosphere the solution does not exist and outputs are NaNs 
%
% PROTOTYPE:
%   [rpSol, hpSol, deltavp, deltaVect, eVect] = powFlyby(vInfMin, vInfPlus, idAtt)
%
% INPUT:
%   vInfMin[3x1]    Entry velocity of the flyby solution in ecliptic components [km/s]
%   vInfPlus[3x1]   Exit velocity of the flyby solution in ecliptic components [km/s]
%   idAtt[1]        Flyby attractor identifier [-]
%
% OUTPUT:
%   rpSol[1]        Radius of perigee of the flyby solution [km]
%   hpSol[1]        Height of perigee of the flyby solution [km]
%   deltavp[1]      Impulse given at the perigee of the flyby trajectory [km/s]
%   deltaVect[3x1]  Total deflection angle and total deflection angles of both 
%                   hyperbolic branches [rad]
%   eVect[2x1]      Eccentricities of both hyperbolic branches [-]
% 
% CONTRIBUTORS:
%   Fabio Spada
%   Ciro Salvi
%
% VERSIONS:
%   2021-02-11
%

mu             = astroConstants(10+idAtt);
rPlanet        = astroConstants(20+idAtt);

switch idAtt
    case 2
        hAtm = 100; % [km]  
end

vMin        = norm(vInfMin);
vPlus       = norm(vInfPlus);
delta       = acos(dot(vInfMin,vInfPlus)/(vMin*vPlus)); 
eGuess      = 1/sin(delta/2);


%%% option 1
rpGuess1 = (eGuess - 1)*mu/vMin^2;
%%% option 2
rpGuess2 = (eGuess - 1)*mu/vPlus^2;

eMin        = @(rp) 1 + rp*vMin^2/mu;
ePlus       = @(rp) 1 + rp*vPlus^2/mu;
deltaMin    = @(rp) 2*asin(1/eMin(rp));
deltaPlus   = @(rp) 2*asin(1/ePlus(rp));

fdelta      = @(rp) delta - deltaMin(rp)/2 - deltaPlus(rp)/2;

% try 
    rpSol       = fzero(fdelta, [rpGuess1 rpGuess2]);
    hpSol       = rpSol - rPlanet;
% catch
%     rpSol = NaN;
%     hpSol = NaN;
% end

if rpSol <=  rPlanet + hAtm
    rpSol = NaN;
    vpMin = NaN; 
    vpPlus = NaN; 
    deltavp = NaN; 
    
else
    vpMin       = sqrt(vMin^2 + 2*mu/rpSol);
    vpPlus      = sqrt(vPlus^2 + 2*mu/rpSol);
    deltavp     = norm(vpPlus - vpMin);
end

eVect = [eMin(rpSol); ePlus(rpSol)];
deltaVect = [delta; deltaMin(rpSol); deltaPlus(rpSol)];

end