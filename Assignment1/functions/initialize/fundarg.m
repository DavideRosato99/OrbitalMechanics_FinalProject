function [l, l1, f, d, omega] = fundarg(ttt)
%  fundarg - the function calulates the delauany variables and planetary 
%            values for iau 1980 theory
%
%  PROTOTYPE:
% [ l, l1, f, d, omega, ...
% lonmer, lonven, lonear, lonmar, lonjup, lonsat, lonurn, lonnep, precrate]... 
%  = fundarg( ttt, opt );
%
%  INPUTS:
%    ttt         - julian centuries of tt
%    opt         - method option                  '06', '02', '96', '80'
%
%  OUTPUTS:
%    l           - delaunay element               rad
%    ll          - delaunay element               rad
%    f           - delaunay element               rad
%    d           - delaunay element               rad
%    omega       - delaunay element               rad
%    planetary longitudes                         rad
%
%  LOCALS:
%    ttt2,ttt3,  - powers of ttt
%
%  CALLED FUNCTIONS:
%    none        
%
%  AUTHOR: 
%   David Vallado  719-573-2600   16 jul 2004
%
%  REFERENCES:
%    vallado       2004, 212-214
%
% CHANGELOG:
%     2021-12-17, 2021-2022 Assignments changes for practical uses
%
% CONTRIBUTORS:
%   Rosato Davide               10618468
%   Saba Mohammadi Yengeje      10789462
%   Spinelli Jason              10618465
%   Tagliati Alessia            10635119
%  ------------------------------------------------------------------------

l     = ((((0.064) * ttt + 31.310) * ttt + 1717915922.6330) * ttt) / 3600.0 + 134.96298139;
l1    = ((((-0.012) * ttt - 0.577) * ttt + 129596581.2240) * ttt) / 3600.0 + 357.52772333;
f     = ((((0.011) * ttt - 13.257) * ttt + 1739527263.1370) * ttt) / 3600.0 + 93.27191028;
d     = ((((0.019) * ttt - 6.891) * ttt + 1602961601.3280) * ttt) / 3600.0 + 297.85036306;
omega = ((((0.008) * ttt + 7.455) * ttt - 6962890.5390) * ttt) / 3600.0 + 125.04452222;

end




