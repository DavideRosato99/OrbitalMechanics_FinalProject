function [satData] = sgp4Init(satData, epoch)
% sgp4Init - the function initializes satellite data from TLEs
%
%  INPUT:
%   satData     struct    [-]      satellites data from TLE 
%   epoch       double    [Nx1]    satellite TLE epoch(mjd2000) 
%
%  OUTPUT:
%   satData     struct    [-]      initialized satellites data from TLE 
%
% CALLED FUNCTIONS: 
%   wgs84data
%   initl
%   dpper
%   dsinit
%   dscom
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
%% INITIALIZATION
N = size(satData, 1);

% Set all near earth variables to zero
satData.isimp   = zeros(N, 1); satData.method = char('n'*ones(N, 1)); satData.aycof    = zeros(N, 1);
satData.con41   = zeros(N, 1); satData.cc1    = zeros(N, 1); satData.cc4      = zeros(N, 1);
satData.cc5     = zeros(N, 1); satData.d2     = zeros(N, 1); satData.d3       = zeros(N, 1);
satData.d4      = zeros(N, 1); satData.delmo  = zeros(N, 1); satData.eta      = zeros(N, 1);
satData.argpdot = zeros(N, 1); satData.omgcof = zeros(N, 1); satData.sinmao   = zeros(N, 1);
satData.t       = zeros(N, 1); satData.t2cof  = zeros(N, 1); satData.t3cof    = zeros(N, 1);
satData.t4cof   = zeros(N, 1); satData.t5cof  = zeros(N, 1); satData.x1mth2   = zeros(N, 1);
satData.x7thm1  = zeros(N, 1); satData.mdot   = zeros(N, 1); satData.nodedot  = zeros(N, 1);
satData.xlcof   = zeros(N, 1); satData.xmcof  = zeros(N, 1); satData.nodecf   = zeros(N, 1);

% Set all deep space variables to zero
satData.irez    = zeros(N, 1); satData.d2201  = zeros(N, 1); satData.d2211    = zeros(N, 1);
satData.d3210   = zeros(N, 1); satData.d3222  = zeros(N, 1); satData.d4410    = zeros(N, 1);
satData.d4422   = zeros(N, 1); satData.d5220  = zeros(N, 1); satData.d5232    = zeros(N, 1);
satData.d5421   = zeros(N, 1); satData.d5433  = zeros(N, 1); satData.dedt     = zeros(N, 1);
satData.del1    = zeros(N, 1); satData.del2   = zeros(N, 1); satData.del3     = zeros(N, 1);
satData.didt    = zeros(N, 1); satData.dmdt   = zeros(N, 1); satData.dnodt    = zeros(N, 1);
satData.domdt   = zeros(N, 1); satData.e3     = zeros(N, 1); satData.ee2      = zeros(N, 1);
satData.peo     = zeros(N, 1); satData.pgho   = zeros(N, 1); satData.pho      = zeros(N, 1);
satData.pinco   = zeros(N, 1); satData.plo    = zeros(N, 1); satData.se2      = zeros(N, 1);
satData.se3     = zeros(N, 1); satData.sgh2   = zeros(N, 1); satData.sgh3     = zeros(N, 1);
satData.sgh4    = zeros(N, 1); satData.sh2    = zeros(N, 1); satData.sh3      = zeros(N, 1);
satData.si2     = zeros(N, 1); satData.si3    = zeros(N, 1); satData.sl2      = zeros(N, 1);
satData.sl3     = zeros(N, 1); satData.sl4    = zeros(N, 1); satData.gsto     = zeros(N, 1);
satData.xfact   = zeros(N, 1); satData.xgh2   = zeros(N, 1); satData.xgh3     = zeros(N, 1);
satData.xgh4    = zeros(N, 1); satData.xh2    = zeros(N, 1); satData.xh3      = zeros(N, 1);
satData.xi2     = zeros(N, 1); satData.xi3    = zeros(N, 1); satData.xl2      = zeros(N, 1);
satData.xl3     = zeros(N, 1); satData.xl4    = zeros(N, 1); satData.xlamo    = zeros(N, 1);
satData.zmol    = zeros(N, 1); satData.zmos   = zeros(N, 1); satData.atime    = zeros(N, 1);
satData.xli     = zeros(N, 1); satData.xni    = zeros(N, 1);

% Single averaged mean elements
satData.am = zeros(N, 1); satData.em = zeros(N, 1); satData.im = zeros(N, 1);
satData.Om = zeros(N, 1); satData.mm = zeros(N, 1); satData.nm = zeros(N, 1);


%% WGS-84 EARTH CONSTANTS
[tumin, mu, Re, xke, j2, j3, j4, j3oj2] = wgs84data;

satData.tumin = tumin * ones(N, 1);
satData.mu = mu * ones(N, 1);
satData.radiusearthkm = Re * ones(N, 1);
satData.xke = xke * ones(N, 1);
satData.j2 = j2 * ones(N, 1);
satData.j3 = j3 * ones(N, 1);
satData.j4 = j4 * ones(N, 1);
satData.j3oj2 = j3oj2 * ones(N, 1);

ss          = 78/Re + 1;
qzms2t      = ((120.0 - 78.0) / Re)^4;
x2o3        = 2/3;
temp4       = 1.5e-12;


%% INITIALIZE SGP4/SDP4 PROPAGATOR
satData.init = char('y'*ones(N, 1));

ao     = zeros(1, 1);
con42  = zeros(1, 1);
cosio  = zeros(1, 1);
cosio2 = zeros(1, 1);
eccsq  = zeros(1, 1);
omeosq = zeros(1, 1);
posq   = zeros(1, 1);
rp     = zeros(1, 1);
rteosq = zeros(1, 1);
sinio  = zeros(1, 1);

sfour  = zeros(1, 1);
qzms24 = zeros(1, 1);
perige = zeros(1, 1);
pinvsq = zeros(1, 1);


for i = 1:N
    
[satData.method(i), ao(i), satData.con41(i), con42(i), cosio(i), ...
    cosio2(i), eccsq(i), omeosq(i), posq(i), rp(i), rteosq(i),...
    sinio(i), satData.gsto(i), satData.no(i)]...
    = initl(satData.xke(i), satData.j2(i), satData.ecco(i), epoch(i),...
    satData.inclo(i), satData.no_kozai(i));

    satData.a(i)    = (satData.no(i) * satData.tumin(i))^(-2.0/3.0);
    satData.alta(i) = satData.a(i) * (1 + satData.ecco(i)) - 1.0;
    satData.altp(i) = satData.a(i) * (1 - satData.ecco(i)) - 1.0;
    
    if (omeosq(i) >= 0) || satData.no(i) >= 0
        satData.isimp(i) = 0;
        
        if rp(i) < (220 / satData.radiusearthkm(i) + 1)
            satData.isimp(i) = 1;
        end
        
        sfour(i) = ss;
        qzms24(i) = qzms2t;
        perige(i) = (rp(i) - 1) * satData.radiusearthkm(i);
        
        % For perigees below 156 km, s and qoms2t are altered
        if perige(i) < 156
            sfour(i) = perige(i) - 78;
            if perige(i) < 98
                sfour(i) = 20;
            end
            qzms24(i) = ((120 - sfour(i)) / satData.radiusearthkm(i))^4;
            sfour(i) = sfour(i) / satData.radiusearthkm(i) + 1;
        end
        
        pinvsq(i) = 1 / posq(i);
        
        tsi(i)  = 1 / (ao(i) - sfour(i));
        satData.eta(i)  = ao(i) * satData.ecco(i) * tsi(i);
        etasq(i) = satData.eta(i) * satData.eta(i);
        eeta(i)  = satData.ecco(i) * satData.eta(i);
        psisq(i) = abs(1 - etasq(i));
        coef(i)  = qzms24(i) * tsi(i)^4.0;
        coef1(i) = coef(i) / psisq(i)^3.5;
        cc2(i)   = coef1(i) * satData.no(i) * (ao(i) * (1 + 1.5 * etasq(i)...
            + eeta(i) * (4 + etasq(i))) + 0.375 * satData.j2(i) *...
            tsi(i) / psisq(i) * satData.con41(i) * (8 + 3 * etasq(i)...
            * (8 + etasq(i))));
        satData.cc1(i) = satData.bstar(i) * cc2(i);
        cc3(i)   = 0;
        
        if (satData.ecco(i) > 1.0e-4)
            cc3(i) = -2 * coef(i) * tsi(i) * satData.j3oj2(i) * satData.no(i) * sinio(i) / satData.ecco(i);
        end
        satData.x1mth2(i) = 1 - cosio2(i);
        satData.cc4(i)    = 2.0* satData.no(i) * coef1(i) * ao(i) * omeosq(i) *...
           (satData.eta(i) * (2.0 + 0.5 * etasq(i)) + satData.ecco(i) *...
           (0.5 + 2.0 * etasq(i)) - satData.j2(i) * tsi(i) / (ao(i) * psisq(i)) *...
           (-3.0 * satData.con41(i) * (1.0 - 2.0 * eeta(i) + etasq(i) *...
           (1.5 - 0.5 * eeta(i))) + 0.75 * satData.x1mth2(i) *...
           (2.0 * etasq(i) - eeta(i) * (1.0 + etasq(i))) * cos(2.0 * satData.argpo(i))));
        satData.cc5(i) = 2.0 * coef1(i) * ao(i) * omeosq(i) * (1.0 + 2.75 *...
           (etasq(i) + eeta(i)) + eeta(i) * etasq(i));
        cosio4(i) = cosio2(i) * cosio2(i);
        temp1(i)  = 1.5 * satData.j2(i) * pinvsq(i) * satData.no(i);
        temp2(i)  = 0.5 * temp1(i) * satData.j2(i) * pinvsq(i);
        temp3(i)  = -0.46875 * satData.j4(i) * pinvsq(i) * pinvsq(i) * satData.no(i);
        satData.mdot(i)     = satData.no(i) + 0.5 * temp1(i) * rteosq(i) * satData.con41(i) +...
           0.0625 * temp2(i) * rteosq(i) * (13.0 - 78.0 * cosio2(i) + 137.0 * cosio4(i));
        satData.argpdot(i)  = -0.5 * temp1(i) * con42(i) + 0.0625 * temp2(i) *...
           (7.0 - 114.0 * cosio2(i) + 395.0 * cosio4(i)) +...
           temp3(i) * (3.0 - 36.0 * cosio2(i) + 49.0 * cosio4(i));
        xhdot1(i)            = -temp1(i) * cosio(i);
        satData.nodedot(i) = xhdot1(i) + (0.5 * temp2(i) * (4.0 - 19.0 * cosio2(i)) +...
           2.0 * temp3(i) * (3.0 - 7.0 * cosio2(i))) * cosio(i);
        xpidot(i)            =  satData.argpdot(i) + satData.nodedot(i);
        satData.omgcof(i)   = satData.bstar(i) * cc3(i) * cos(satData.argpo(i));
        satData.xmcof(i)    = 0.0;
        if (satData.ecco(i) > 1.0e-4)
           satData.xmcof(i) = -x2o3 * coef(i) * satData.bstar(i) / eeta(i);
        end
        satData.nodecf(i) = 3.5 * omeosq(i) * xhdot1(i) * satData.cc1(i);
        satData.t2cof(i)   = 1.5 * satData.cc1(i);
        
        % Sgp4fix for divide by zero with xinco = 180 deg
        if (abs(cosio(i)+1.0) > 1.5e-12)
            satData.xlcof(i)   = -0.25 * satData.j3oj2(i) * sinio(i) *...
                (3.0 + 5.0 * cosio(i)) / (1.0 + cosio(i));
        else
            satData.xlcof(i)   = -0.25 * satData.j3oj2(i) * sinio(i) *...
                (3.0 + 5.0 * cosio(i)) / temp4(i);
        end   
        satData.aycof(i)   = -0.5 * satData.j3oj2(i) * sinio(i);
        satData.delmo(i)   = (1.0 + satData.eta(i) * cos(satData.mo(i)))^3;
        satData.sinmao(i)  = sin(satData.mo(i));
        satData.x7thm1(i)  = 7.0 * cosio2(i) - 1.0;
        
        % Deep space initialization
        if ((2*pi / satData.no(i)) >= 225.0)
            satData.method(i) = 'd';
            satData.isimp(i)  = 1;
            tc(i)    =  0.0;
            inclm(i) = satData.inclo(i);

            [sinim(i),cosim(i),sinomm(i),cosomm(i),snodm(i),cnodm(i),day(i),satData.e3(i),satData.ee2(i),...
                em(i),emsq(i),gam(i),satData.peo(i),satData.pgho(i),satData.pho(i),satData.pinco(i),...
                satData.plo(i),rtemsq(i),satData.se2(i),satData.se3(i),satData.sgh2(i),...
                satData.sgh3(i),satData.sgh4(i),satData.sh2(i),satData.sh3(i),satData.si2(i),...
                satData.si3(i),satData.sl2(i),satData.sl3(i),satData.sl4(i),s1(i),s2(i),s3(i),s4(i),s5(i),...
                s6(i),s7(i),ss1(i),ss2(i),ss3(i),ss4(i),ss5(i),ss6(i),ss7(i),sz1(i),sz2(i),sz3(i),sz11(i),sz12(i),...
                sz13(i),sz21(i),sz22(i),sz23(i),sz31(i),sz32(i),sz33(i),satData.xgh2(i),satData.xgh3(i),...
                satData.xgh4(i),satData.xh2(i),satData.xh3(i),satData.xi2(i),satData.xi3(i),...
                satData.xl2(i),satData.xl3(i),satData.xl4(i),nm(i),z1(i),z2(i),z3(i),z11(i),z12(i),z13(i),...
                z21(i),z22(i),z23(i),z31(i),z32(i),z33(i),satData.zmol(i),satData.zmos(i)] ...
                = dscom(epoch(i),satData.ecco(i),satData.argpo(i),tc(i),satData.inclo(i),...
                 satData.nodeo(i),satData.no(i));

            [satData.ecco(i),satData.inclo(i),satData.nodeo(i),satData.argpo(i),satData.mo(i)]...
                = dpper(satData.e3(i),satData.ee2(i),satData.peo(i),satData.pgho(i),...
                satData.pho(i),satData.pinco(i),satData.plo(i),satData.se2(i),satData.se3(i),...
                satData.sgh2(i),satData.sgh3(i),satData.sgh4(i),satData.sh2(i),satData.sh3(i),...
                satData.si2(i),satData.si3(i),satData.sl2(i),satData.sl3(i),satData.sl4(i),...
                satData.t(i),satData.xgh2(i),satData.xgh3(i),satData.xgh4(i),satData.xh2(i),...
                satData.xh3(i),satData.xi2(i),satData.xi3(i),satData.xl2(i),satData.xl3(i),...
                satData.xl4(i),satData.zmol(i),satData.zmos(i),inclm(i),satData.init(i),...
                satData.ecco(i),satData.inclo(i),satData.nodeo(i),satData.argpo(i),satData.mo(i));

            argpm(i)  = 0.0;
            nodem(i)  = 0.0;
            mm(i)     = 0.0;

            [em(i),argpm(i),inclm(i),mm(i),nm(i),nodem(i),satData.irez(i),satData.atime(i),...
                satData.d2201(i),satData.d2211(i),satData.d3210(i),satData.d3222(i),...
                satData.d4410(i),satData.d4422(i),satData.d5220(i),satData.d5232(i),...
                satData.d5421(i),satData.d5433(i),satData.dedt(i),satData.didt(i),...
                satData.dmdt(i),dndt(i),satData.dnodt(i),satData.domdt(i),satData.del1(i),...
                satData.del2(i),satData.del3(i),...
                satData.xfact(i),satData.xlamo(i),satData.xli(i),satData.xni(i)] ...
                = dsinit(...
                satData.xke(i), cosim(i),emsq(i),satData.argpo(i),s1(i),s2(i),s3(i),s4(i),s5(i),sinim(i),ss1(i),ss2(i),ss3(i),...
                ss4(i),ss5(i),sz1(i),sz3(i),sz11(i),sz13(i),sz21(i),sz23(i),sz31(i),sz33(i),satData.t(i),tc(i),...
                satData.gsto(i),satData.mo(i),satData.mdot(i),satData.no(i),satData.nodeo(i),...
                satData.nodedot(i),xpidot(i),z1(i),z3(i),z11(i),z13(i),z21(i),z23(i),z31(i),z33(i),em(i),...
                argpm(i),inclm(i),mm(i),nm(i),nodem(i),satData.ecco(i),eccsq(i));
        end
        
        % Set variables if not deep space
        if (satData.isimp(i) ~= 1)
            cc1sq(i)          = satData.cc1(i) * satData.cc1(i);
            satData.d2(i)    = 4.0 * ao(i) * tsi(i) * cc1sq(i);
            temp(i)           = satData.d2(i) * tsi(i) * satData.cc1(i) / 3.0;
            satData.d3(i)    = (17.0 * ao(i) + sfour(i)) * temp(i);
            satData.d4(i)    = 0.5 * temp(i) * ao(i) * tsi(i) *...
               (221.0 * ao(i) + 31.0 * sfour(i)) * satData.cc1(i);
            satData.t3cof(i) = satData.d2(i) + 2.0 * cc1sq(i);
            satData.t4cof(i) = 0.25 * (3.0 * satData.d3(i) + satData.cc1(i) *...
               (12.0 * satData.d2(i) + 10.0 * cc1sq(i)));
            satData.t5cof(i) = 0.2 * (3.0 * satData.d4(i) +...
               12.0 * satData.cc1(i) * satData.d3(i) +...
               6.0 * satData.d2(i) * satData.d2(i) +...
               15.0 * cc1sq(i) * (2.0 * satData.d2(i) + cc1sq(i)));
        end
        
    end
    
end


%% FINAL PARSING OF DATA
names = fieldnames(satData);

for i = 1:length(names)
    if size(satData.(names{i}),1) == 1
        satData.(names{i}) = satData.(names{i})';
    end
end
    

end
