function prec = precess(centjd)
 %  precess - the function calulates the transformation matrix that accounts 
 %            for the effects of precession. In the original function both 
 %            the 1980 and 2006 theories are handled. Here only the 2006
 %            theory is valid
 %
 %  PROTOTYPE:
 %     [prec,psia,wa,ea,xa] = precess ( ttt, opt );
 %  
 %  INPUTS:         description                    range / units
 %    ttt         - julian centuries of tt
 %    opt         - method option                  '01', '02', '96', '80'
 %
 %  OUTPUTS:
 %    prec        - transformation matrix for mod - j2000 (80 only)
 %    psia        - cannonical precession angle    rad    (00 only)
 %    wa          - cannonical precession angle    rad    (00 only)
 %    ea          - cannonical precession angle    rad    (00 only)
 %    xa          - cannonical precession angle    rad    (00 only)
 %
 %  LOCALS:
 %    ttt2        - ttt squared
 %    ttt3        - ttt cubed
 %    zeta        - precession angle               rad
 %    z           - precession angle               rad
 %    theta       - precession angle               rad
 %    oblo        - obliquity value at j2000 epoch "%
 %
 %  CALLED FUNCTION: [-]
 % 
 %  AUTHOR: 
 %    david vallado 719-573-2600   25 jun 2002
 %
 %  REVINSION:
 %    vallado  - consolidate with iau 2000   14 feb 2005
 %
 %  REFERENCES:
 %    vallado       2004, 214-216, 219-221
 %
 % CHANGELOG:
 %     2021-12-17, 2021-2022 Assignments changes for practical uses
 %
 % CONTRIBUTORS:
 %   Rosato Davide               10618468
 %   Saba Mohammadi Yengeje      10789462
 %   Spinelli Jason              10618465
 %   Tagliati Alessia            10635119
 % ------------------------------------------------------------------------
convrt = pi / (180.0*3600.0);
ttt2 = centjd * centjd;
ttt3 = ttt2 * centjd;

oblo =  84381.406; 
psia =  (((( -0.0000000951 * centjd + 0.000132851 ) * centjd - 0.00114045 ) * centjd - 1.0790069 ) * centjd + 5038.481507 ) * centjd; % "
wa   =  ((((  0.0000003337 * centjd - 0.000000467 ) * centjd - 0.00772503 ) * centjd + 0.0512623 ) * centjd -    0.025754 ) * centjd + oblo;
ea   =  (((( -0.0000000434 * centjd - 0.000000576 ) * centjd + 0.00200340 ) * centjd - 0.0001831 ) * centjd -   46.836769 ) * centjd + oblo;
xa   =  (((( -0.0000000560 * centjd + 0.000170663 ) * centjd - 0.00121197 ) * centjd - 2.3814292 ) * centjd +   10.556403 ) * centjd;

zeta =  (((( -0.0000003173 * centjd - 0.000005971 ) * centjd + 0.01801828 ) * centjd + 0.2988499 ) * centjd + 2306.083227 ) * centjd + 2.650545; % "
theta=  (((( -0.0000001274 * centjd - 0.000007089 ) * centjd - 0.04182264 ) * centjd - 0.4294934 ) * centjd + 2004.191903 ) * centjd;
z    =  ((((  0.0000002904 * centjd - 0.000028596 ) * centjd + 0.01826837 ) * centjd + 1.0927348 ) * centjd + 2306.077181 ) * centjd - 2.650545;


% convert units to rad
psia = psia  * convrt; % rad
wa   = wa    * convrt;
ea   = ea    * convrt;
xa   = xa    * convrt;

zeta = zeta  * convrt;
theta= theta * convrt;
z    = z     * convrt;

    
coszeta  = cos(zeta);
sinzeta  = sin(zeta);
costheta = cos(theta);
sintheta = sin(theta);
cosz     = cos(z);
sinz     = sin(z);

% ----------------- form matrix mod to j2000 -----------------
prec(1,1) =  coszeta * costheta * cosz - sinzeta * sinz;
prec(1,2) =  coszeta * costheta * sinz + sinzeta * cosz;
prec(1,3) =  coszeta * sintheta;
prec(2,1) = -sinzeta * costheta * cosz - coszeta * sinz;
prec(2,2) = -sinzeta * costheta * sinz + coszeta * cosz;
prec(2,3) = -sinzeta * sintheta;
prec(3,1) = -sintheta * cosz;
prec(3,2) = -sintheta * sinz;
prec(3,3) =  costheta;
