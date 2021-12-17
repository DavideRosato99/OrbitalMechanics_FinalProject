function [satData] = sgp4Init(satData, epoch)

% [satData] = sgp4Init(whichconst, opsmode, satData,...
%     satData.jdsatepoch+satData.jdsatepochf-2433281.5, satData.bstar, ...
%     satData.ndot, satData.nddot, satData.ecco, satData.argpo, ...
%     satData.inclo, satData.mo, satData.no_kozai, satData.nodeo);
%% INITIALIZATION
N = size(satData.Name, 2);

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



    

end


%% FUNCTIONS
function [tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2] = wgs84data

mu            = astroConstants(13);
radiusearthkm = astroConstants(23);
xke           = 60.0 / sqrt((radiusearthkm^3)/mu);
tumin         = 1.0 / xke;
j2            =   0.00108262998905;
j3            =  -0.00000253215306;
j4            =  -0.00000161098761;
j3oj2         =  j3 / j2;

end

function [method, ao, con41, con42, cosio, cosio2, eccsq, omeosq,...
    posq, rp, rteosq, sinio, gsto, no_unkozai] = initl(xke, j2, ecco,...
    epoch, inclo, no_kozai)

% WGS-84 earth constants
x2o3 = 2/3;

% Calculate auxillary epoch quantities
eccsq  = ecco * ecco;
omeosq = 1.0 - eccsq;
rteosq = sqrt(omeosq);
cosio  = cos(inclo);
cosio2 = cosio^2;

% Un-kozai the mean motion
ak           = (xke / no_kozai)^x2o3;
d1          = 0.75 * j2 * (3.0 * cosio2 - 1.0) / (rteosq * omeosq);
del         = d1 / (ak * ak);
adel        = ak * (1.0 - del * del - del * (1.0 / 3.0 + 134.0 * del * del / 81.0));
del         = d1/(adel * adel);
no_unkozai  = no_kozai / (1.0 + del);

ao     = (xke / no_unkozai)^x2o3;
sinio  = sin(inclo);
po     = ao * omeosq;
con42  = 1.0 - 5.0 * cosio2;
con41  = -con42-cosio2-cosio2;
posq   = po * po;
rp     = ao * (1.0 - ecco);
method = 'n';

% SGP4 modern approach to finding sidereal time
gsto = mjd20002gmst(epoch); 
end


function [ sinim,cosim,sinomm,cosomm,snodm,cnodm,day,e3,ee2,em,emsq,gam,...
           peo,pgho,pho,pinco,plo,rtemsq,se2,se3,sgh2,sgh3,sgh4,sh2,sh3,si2,...
           si3,sl2,sl3,sl4,s1,s2,s3,s4,s5,s6,s7,ss1,ss2,ss3,ss4,ss5,ss6,ss7,...
           sz1,sz2,sz3,sz11,sz12,sz13,sz21,sz22,sz23,sz31,sz32,sz33,xgh2,xgh3,...
           xgh4,xh2,xh3,xi2,xi3,xl2,xl3,xl4,nm,z1,z2,z3,z11,z12,z13,z21,z22,...
           z23,z31,z32,z33,zmol,zmos]...
         = dscom (epoch, ep, argpp, tc, inclp, nodep, np)

   % /* -------------------------- constants ------------------------- */
   zes     =  0.01675;
   zel     =  0.05490;
   c1ss    =  2.9864797e-6;
   c1l     =  4.7968065e-7;
   zsinis  =  0.39785416;
   zcosis  =  0.91744867;
   zcosgs  =  0.1945905;
   zsings  = -0.98088458;
   twopi   =  2.0 * pi;

   % /* --------------------- local variables ------------------------ */
   nm     = np;
   em     = ep;
   snodm  = sin(nodep);
   cnodm  = cos(nodep);
   sinomm = sin(argpp);
   cosomm = cos(argpp);
   sinim  = sin(inclp);
   cosim  = cos(inclp);
   emsq   = em * em;
   betasq = 1.0 - emsq;
   rtemsq = sqrt(betasq);

   % /* ----------------- initialize lunar solar terms --------------- */
   peo    = 0.0;
   pinco  = 0.0;
   plo    = 0.0;
   pgho   = 0.0;
   pho    = 0.0;
   day    = epoch + 18261.5 + tc / 1440.0;
   xnodce = rem(4.5236020 - 9.2422029e-4 * day, twopi);
   stem   = sin(xnodce);
   ctem   = cos(xnodce);
   zcosil = 0.91375164 - 0.03568096 * ctem;
   zsinil = sqrt(1.0 - zcosil * zcosil);
   zsinhl = 0.089683511 * stem / zsinil;
   zcoshl = sqrt(1.0 - zsinhl * zsinhl);
   gam    = 5.8351514 + 0.0019443680 * day;
   zx     = 0.39785416 * stem / zsinil;
   zy     = zcoshl * ctem + 0.91744867 * zsinhl * stem;
   zx     = atan2(zx, zy);
   zx     = gam + zx - xnodce;
   zcosgl = cos(zx);
   zsingl = sin(zx);

   % /* ------------------------- do solar terms --------------------- */
   zcosg = zcosgs;
   zsing = zsings;
   zcosi = zcosis;
   zsini = zsinis;
   zcosh = cnodm;
   zsinh = snodm;
   cc    = c1ss;
   xnoi  = 1.0 / nm;

   for (lsflg = 1:2)
       a1  =   zcosg * zcosh + zsing * zcosi * zsinh;
       a3  =  -zsing * zcosh + zcosg * zcosi * zsinh;
       a7  =  -zcosg * zsinh + zsing * zcosi * zcosh;
       a8  =   zsing * zsini;
       a9  =   zsing * zsinh + zcosg * zcosi * zcosh;
       a10 =   zcosg * zsini;
       a2  =   cosim * a7 + sinim * a8;
       a4  =   cosim * a9 + sinim * a10;
       a5  =  -sinim * a7 + cosim * a8;
       a6  =  -sinim * a9 + cosim * a10;

       x1  =  a1 * cosomm + a2 * sinomm;
       x2  =  a3 * cosomm + a4 * sinomm;
       x3  = -a1 * sinomm + a2 * cosomm;
       x4  = -a3 * sinomm + a4 * cosomm;
       x5  =  a5 * sinomm;
       x6  =  a6 * sinomm;
       x7  =  a5 * cosomm;
       x8  =  a6 * cosomm;

       z31 = 12.0 * x1 * x1 - 3.0 * x3 * x3;
       z32 = 24.0 * x1 * x2 - 6.0 * x3 * x4;
       z33 = 12.0 * x2 * x2 - 3.0 * x4 * x4;
       z1  =  3.0 *  (a1 * a1 + a2 * a2) + z31 * emsq;
       z2  =  6.0 *  (a1 * a3 + a2 * a4) + z32 * emsq;
       z3  =  3.0 *  (a3 * a3 + a4 * a4) + z33 * emsq;
       z11 = -6.0 * a1 * a5 + emsq *  (-24.0 * x1 * x7-6.0 * x3 * x5);
       z12 = -6.0 *  (a1 * a6 + a3 * a5) + emsq *...
           (-24.0 * (x2 * x7 + x1 * x8) - 6.0 * (x3 * x6 + x4 * x5));
       z13 = -6.0 * a3 * a6 + emsq * (-24.0 * x2 * x8 - 6.0 * x4 * x6);
       z21 =  6.0 * a2 * a5 + emsq * (24.0 * x1 * x5 - 6.0 * x3 * x7);
       z22 =  6.0 *  (a4 * a5 + a2 * a6) + emsq *...
           (24.0 * (x2 * x5 + x1 * x6) - 6.0 * (x4 * x7 + x3 * x8));
       z23 =  6.0 * a4 * a6 + emsq * (24.0 * x2 * x6 - 6.0 * x4 * x8);
       z1  = z1 + z1 + betasq * z31;
       z2  = z2 + z2 + betasq * z32;
       z3  = z3 + z3 + betasq * z33;
       s3  = cc * xnoi;
       s2  = -0.5 * s3 / rtemsq;
       s4  = s3 * rtemsq;
       s1  = -15.0 * em * s4;
       s5  = x1 * x3 + x2 * x4;
       s6  = x2 * x3 + x1 * x4;
       s7  = x2 * x4 - x1 * x3;

       % /* ----------------------- do lunar terms ------------------- */
       if (lsflg == 1)
           ss1   = s1;
           ss2   = s2;
           ss3   = s3;
           ss4   = s4;
           ss5   = s5;
           ss6   = s6;
           ss7   = s7;
           sz1   = z1;
           sz2   = z2;
           sz3   = z3;
           sz11  = z11;
           sz12  = z12;
           sz13  = z13;
           sz21  = z21;
           sz22  = z22;
           sz23  = z23;
           sz31  = z31;
           sz32  = z32;
           sz33  = z33;
           zcosg = zcosgl;
           zsing = zsingl;
           zcosi = zcosil;
           zsini = zsinil;
           zcosh = zcoshl * cnodm + zsinhl * snodm;
           zsinh = snodm * zcoshl - cnodm * zsinhl;
           cc    = c1l;
       end
   end

   zmol = rem(4.7199672 + 0.22997150  * day - gam, twopi);
   zmos = rem(6.2565837 + 0.017201977 * day, twopi);

   % /* ------------------------ do solar terms ---------------------- */
   se2  =   2.0 * ss1 * ss6;
   se3  =   2.0 * ss1 * ss7;
   si2  =   2.0 * ss2 * sz12;
   si3  =   2.0 * ss2 * (sz13 - sz11);
   sl2  =  -2.0 * ss3 * sz2;
   sl3  =  -2.0 * ss3 * (sz3 - sz1);
   sl4  =  -2.0 * ss3 * (-21.0 - 9.0 * emsq) * zes;
   sgh2 =   2.0 * ss4 * sz32;
   sgh3 =   2.0 * ss4 * (sz33 - sz31);
   sgh4 = -18.0 * ss4 * zes;
   sh2  =  -2.0 * ss2 * sz22;
   sh3  =  -2.0 * ss2 * (sz23 - sz21);

   % /* ------------------------ do lunar terms ---------------------- */
   ee2  =   2.0 * s1 * s6;
   e3   =   2.0 * s1 * s7;
   xi2  =   2.0 * s2 * z12;
   xi3  =   2.0 * s2 * (z13 - z11);
   xl2  =  -2.0 * s3 * z2;
   xl3  =  -2.0 * s3 * (z3 - z1);
   xl4  =  -2.0 * s3 * (-21.0 - 9.0 * emsq) * zel;
   xgh2 =   2.0 * s4 * z32;
   xgh3 =   2.0 * s4 * (z33 - z31);
   xgh4 = -18.0 * s4 * zel;
   xh2  =  -2.0 * s2 * z22;
   xh3  =  -2.0 * s2 * (z23 - z21);
end


function [  ep,     inclp,  nodep, argpp,  mp]...
          = dpper(...
            e3,     ee2,    peo,    pgho,   pho,    pinco,  plo,    se2,...
            se3,    sgh2,   sgh3,   sgh4,   sh2,    sh3,    si2,    si3,...
            sl2,    sl3,    sl4,    t,      xgh2,   xgh3,   xgh4,   xh2,...
            xh3,    xi2,    xi3,    xl2,    xl3,    xl4,    zmol,...
            zmos,   inclo,  init,   ep,     inclp,  nodep, argpp,  mp)

   % change to variable passed in
   % global opsmode

   % /* --------------------- local variables ------------------------ */
   twopi = 2.0 * pi;

   % /* ---------------------- constants ----------------------------- */
   zns   = 1.19459e-5;
   zes   = 0.01675;
   znl   = 1.5835218e-4;
   zel   = 0.05490;

   % /* --------------- calculate time varying periodics ----------- */
   zm    = zmos + zns * t;
   % // be sure that the initial call has time set to zero
   if (init == 'y')
       zm = zmos;
   end
   zf    = zm + 2.0 * zes * sin(zm);
   sinzf = sin(zf);
   f2    =  0.5 * sinzf * sinzf - 0.25;
   f3    = -0.5 * sinzf * cos(zf);
   ses   = se2* f2 + se3 * f3;
   sis   = si2 * f2 + si3 * f3;
   sls   = sl2 * f2 + sl3 * f3 + sl4 * sinzf;
   sghs  = sgh2 * f2 + sgh3 * f3 + sgh4 * sinzf;
   shs   = sh2 * f2 + sh3 * f3;
   zm    = zmol + znl * t;
   if (init == 'y')
       zm = zmol;
   end
   zf    = zm + 2.0 * zel * sin(zm);
   sinzf = sin(zf);
   f2    =  0.5 * sinzf * sinzf - 0.25;
   f3    = -0.5 * sinzf * cos(zf);
   sel   = ee2 * f2 + e3 * f3;
   sil   = xi2 * f2 + xi3 * f3;
   sll   = xl2 * f2 + xl3 * f3 + xl4 * sinzf;
   sghl  = xgh2 * f2 + xgh3 * f3 + xgh4 * sinzf;
   shll  = xh2 * f2 + xh3 * f3;
   pe    = ses + sel;
   pinc  = sis + sil;
   pl    = sls + sll;
   pgh   = sghs + sghl;
   ph    = shs + shll;

   if (init == 'n')
       %   //  0.2 rad = 11.45916 deg
       pe    = pe - peo;
       pinc  = pinc - pinco;
       pl    = pl - plo;
       pgh   = pgh - pgho;
       ph    = ph - pho;
       inclp = inclp + pinc;
       ep    = ep + pe;
       sinip = sin(inclp);
       cosip = cos(inclp);

       % /* ----------------- apply periodics directly ------------ */
       %   //  sgp4fix for lyddane choice
       %   //  strn3 used original inclination - this is technically feasible
       %   //  gsfc used perturbed inclination - also technically feasible
       %   //  probably best to readjust the 0.2 limit value and limit discontinuity
       %   //  use next line for original strn3 approach and original inclination
       %   //  if (inclo >= 0.2)
       %   //  use next line for gsfc version and perturbed inclination
       if (inclp >= 0.2)
           ph     = ph / sinip;
           pgh    = pgh - cosip * ph;
           argpp  = argpp + pgh;
           nodep  = nodep + ph;
           mp     = mp + pl;
       else
           % /* ---- apply periodics with lyddane modification ---- */
           sinop  = sin(nodep);
           cosop  = cos(nodep);
           alfdp  = sinip * sinop;
           betdp  = sinip * cosop;
           dalf   =  ph * cosop + pinc * cosip * sinop;
           dbet   = -ph * sinop + pinc * cosip * cosop;
           alfdp  = alfdp + dalf;
           betdp  = betdp + dbet;
           nodep  = rem(nodep, twopi);
           % sgp4fix for afspc written intrinsic functions 
           % nodep used without a trigonometric function ahead
           xls    = mp + argpp + cosip * nodep;
           dls    = pl + pgh - pinc * nodep * sinip;
           xls    = xls + dls;
           xnoh   = nodep;
           nodep  = atan2(alfdp, betdp);
           % sgp4fix for afspc written intrinsic functions 
           % nodep used without a trigonometric function ahead
           if (abs(xnoh - nodep) > pi)
               if (nodep < xnoh)
                   nodep = nodep + twopi;
               else
                   nodep = nodep - twopi;
               end
           end
           mp    = mp + pl;
           argpp = xls - mp - cosip * nodep;
       end
   end
end

function [  em,     argpm,  inclm,  mm,     nm,     nodem, irez,...
            atime,  d2201,  d2211,  d3210,  d3222,  d4410,  d4422,...
            d5220,  d5232,  d5421,  d5433,  dedt,   didt,   dmdt,...
            dndt,   dnodt,  domdt,  del1,   del2,   del3,   xfact,...
            xlamo,  xli,    xni]...
          = dsinit( ...
            xke,    cosim,  emsq,   argpo,  s1,     s2,     s3,     s4,...
            s5,     sinim,  ss1,    ss2,    ss3,    ss4,    ss5,...
            sz1,    sz3,    sz11,   sz13,   sz21,   sz23,   sz31,...
            sz33,   t,      tc,     gsto,   mo,     mdot,   no,...
            nodeo,  nodedot,xpidot, z1,     z3,     z11,...
            z13,    z21,    z23,    z31,    z33,    em,     argpm,...
            inclm,  mm,     nm,     nodem,  ecco,   eccsq)

   % /* --------------------- local variables ------------------------ */
   twopi = 2.0 * pi;
   aonv  = 0.0;
   q22    = 1.7891679e-6;
   q31    = 2.1460748e-6;
   q33    = 2.2123015e-7;
   root22 = 1.7891679e-6;
   root44 = 7.3636953e-9;
   root54 = 2.1765803e-9;
   rptim  = 4.37526908801129966e-3;
   root32 = 3.7393792e-7;
   root52 = 1.1428639e-7;
   x2o3   = 2.0 / 3.0;
   znl    = 1.5835218e-4;
   zns    = 1.19459e-5;

   %     // sgp4fix identify constants and allow alternate values
   % sgp4fix no longer needed, pass xke in
   %global tumin mu radiusearthkm xke j2 j3 j4 j3oj2  

   % /* -------------------- deep space initialization ------------ */
   irez = 0;
   if ((nm < 0.0052359877) && (nm > 0.0034906585))
       irez = 1;
   end
   if ((nm >= 8.26e-3) && (nm <= 9.24e-3) && (em >= 0.5))
       irez = 2;
   end
   d2201 = 0;
   d2211 = 0;
   d3210 = 0;
   d3222 = 0;
   d4410 = 0;
   d4422 = 0;
   d5220 = 0;
   d5232 = 0;
   d5421 = 0;
   d5433 = 0;
   del1  = 0;
   del2  = 0;
   del3  = 0;
   atime = 0;
   xfact = 0;
   xlamo = 0;
   xli   = 0;
   xni   = 0;

   % /* ------------------------ do solar terms ------------------- */
   ses  =  ss1 * zns * ss5;
   sis  =  ss2 * zns * (sz11 + sz13);
   sls  = -zns * ss3 * (sz1 + sz3 - 14.0 - 6.0 * emsq);
   sghs =  ss4 * zns * (sz31 + sz33 - 6.0);
   shs  = -zns * ss2 * (sz21 + sz23);
   %   // sgp4fix for 180 deg incl
   if ((inclm < 5.2359877e-2) || (inclm > pi - 5.2359877e-2))
       shs = 0.0;
   end
   if (sinim ~= 0.0)
       shs = shs / sinim;
   end
   sgs  = sghs - cosim * shs;

   % /* ------------------------- do lunar terms ------------------ */
   dedt = ses + s1 * znl * s5;
   didt = sis + s2 * znl * (z11 + z13);
   dmdt = sls - znl * s3 * (z1 + z3 - 14.0 - 6.0 * emsq);
   sghl = s4 * znl * (z31 + z33 - 6.0);
   shll = -znl * s2 * (z21 + z23);
   %   // sgp4fix for 180 deg incl
   if ((inclm < 5.2359877e-2) || (inclm > pi - 5.2359877e-2))
       shll = 0.0;
   end
   domdt = sgs + sghl;
   dnodt = shs;
   if (sinim ~= 0.0)
       domdt = domdt - cosim / sinim * shll;
       dnodt = dnodt + shll / sinim;
   end

   % /* ----------- calculate deep space resonance effects -------- */
   dndt   = 0.0;
   theta  = rem(gsto + tc * rptim, twopi);
   em     = em + dedt * t;
   inclm  = inclm + didt * t;
   argpm  = argpm + domdt * t;
   nodem  = nodem + dnodt * t;
   mm     = mm + dmdt * t;
   % //   sgp4fix for negative inclinations
   % //   the following if statement should be commented out
   % //if (inclm < 0.0)
   % //  {
   % //    inclm  = -inclm;
   % //    argpm  = argpm - pi;
   % //    nodem = nodem + pi;
   % //  }

   %  /* - update resonances : numerical (euler-maclaurin) integration - */
   %  /* ------------------------- epoch restart ----------------------  */
   %  //   sgp4fix for propagator problems
   %  //   the following integration works for negative time steps and periods
   %  //   the specific changes are unknown because the original code was so convoluted

   % /* -------------- initialize the resonance terms ------------- */
   if (irez ~= 0)
       aonv = (nm / xke)^x2o3;

       % /* ---------- geopotential resonance for 12 hour orbits ------ */
       if (irez == 2)
           cosisq = cosim * cosim;
           emo    = em;
           em     = ecco;
           emsqo  = emsq;
           emsq   = eccsq;
           eoc    = em * emsq;
           g201   = -0.306 - (em - 0.64) * 0.440;

           if (em <= 0.65)
               g211 =    3.616  -  13.2470 * em +  16.2900 * emsq;
               g310 =  -19.302  + 117.3900 * em - 228.4190 * emsq +  156.5910 * eoc;
               g322 =  -18.9068 + 109.7927 * em - 214.6334 * emsq +  146.5816 * eoc;
               g410 =  -41.122  + 242.6940 * em - 471.0940 * emsq +  313.9530 * eoc;
               g422 = -146.407  + 841.8800 * em - 1629.014 * emsq + 1083.4350 * eoc;
               g520 = -532.114  + 3017.977 * em - 5740.032 * emsq + 3708.2760 * eoc;
           else
               g211 =   -72.099 +   331.819 * em -   508.738 * emsq +   266.724 * eoc;
               g310 =  -346.844 +  1582.851 * em -  2415.925 * emsq +  1246.113 * eoc;
               g322 =  -342.585 +  1554.908 * em -  2366.899 * emsq +  1215.972 * eoc;
               g410 = -1052.797 +  4758.686 * em -  7193.992 * emsq +  3651.957 * eoc;
               g422 = -3581.690 + 16178.110 * em - 24462.770 * emsq + 12422.520 * eoc;
               if (em > 0.715)
                   g520 =-5149.66 + 29936.92 * em - 54087.36 * emsq + 31324.56 * eoc;
               else
                   g520 = 1464.74 -  4664.75 * em +  3763.64 * emsq;
               end
           end
           if (em < 0.7)
               g533 = -919.22770 + 4988.6100 * em - 9064.7700 * emsq + 5542.21  * eoc;
               g521 = -822.71072 + 4568.6173 * em - 8491.4146 * emsq + 5337.524 * eoc;
               g532 = -853.66600 + 4690.2500 * em - 8624.7700 * emsq + 5341.4  * eoc;
           else
               g533 =-37995.780 + 161616.52 * em - 229838.20 * emsq + 109377.94 * eoc;
               g521 =-51752.104 + 218913.95 * em - 309468.16 * emsq + 146349.42 * eoc;
               g532 =-40023.880 + 170470.89 * em - 242699.48 * emsq + 115605.82 * eoc;
           end

           sini2=  sinim * sinim;
           f220 =  0.75 * (1.0 + 2.0 * cosim+cosisq);
           f221 =  1.5 * sini2;
           f321 =  1.875 * sinim  *  (1.0 - 2.0 * cosim - 3.0 * cosisq);
           f322 = -1.875 * sinim  *  (1.0 + 2.0 * cosim - 3.0 * cosisq);
           f441 = 35.0 * sini2 * f220;
           f442 = 39.3750 * sini2 * sini2;
           f522 =  9.84375 * sinim * (sini2 * (1.0 - 2.0 * cosim- 5.0 * cosisq) +...
               0.33333333 * (-2.0 + 4.0 * cosim + 6.0 * cosisq) );
           f523 = sinim * (4.92187512 * sini2 * (-2.0 - 4.0 * cosim +...
               10.0 * cosisq) + 6.56250012 * (1.0+2.0 * cosim - 3.0 * cosisq));
           f542 = 29.53125 * sinim * (2.0 - 8.0 * cosim+cosisq *...
               (-12.0 + 8.0 * cosim + 10.0 * cosisq));
           f543 = 29.53125 * sinim * (-2.0 - 8.0 * cosim+cosisq *...
               (12.0 + 8.0 * cosim - 10.0 * cosisq));
           xno2  =  nm * nm;
           ainv2 =  aonv * aonv;
           temp1 =  3.0 * xno2 * ainv2;
           temp  =  temp1 * root22;
           d2201 =  temp * f220 * g201;
           d2211 =  temp * f221 * g211;
           temp1 =  temp1 * aonv;
           temp  =  temp1 * root32;
           d3210 =  temp * f321 * g310;
           d3222 =  temp * f322 * g322;
           temp1 =  temp1 * aonv;
           temp  =  2.0 * temp1 * root44;
           d4410 =  temp * f441 * g410;
           d4422 =  temp * f442 * g422;
           temp1 =  temp1 * aonv;
           temp  =  temp1 * root52;
           d5220 =  temp * f522 * g520;
           d5232 =  temp * f523 * g532;
           temp  =  2.0 * temp1 * root54;
           d5421 =  temp * f542 * g521;
           d5433 =  temp * f543 * g533;
           xlamo =  rem(mo + nodeo + nodeo-theta - theta, twopi);
           xfact =  mdot + dmdt + 2.0 * (nodedot + dnodt - rptim) - no;
           em    = emo;
           emsq  = emsqo;
       end

       % /* ---------------- synchronous resonance terms -------------- */
       if (irez == 1)
           g200  = 1.0 + emsq * (-2.5 + 0.8125 * emsq);
           g310  = 1.0 + 2.0 * emsq;
           g300  = 1.0 + emsq * (-6.0 + 6.60937 * emsq);
           f220  = 0.75 * (1.0 + cosim) * (1.0 + cosim);
           f311  = 0.9375 * sinim * sinim * (1.0 + 3.0 * cosim) - 0.75 * (1.0 + cosim);
           f330  = 1.0 + cosim;
           f330  = 1.875 * f330 * f330 * f330;
           del1  = 3.0 * nm * nm * aonv * aonv;
           del2  = 2.0 * del1 * f220 * g200 * q22;
           del3  = 3.0 * del1 * f330 * g300 * q33 * aonv;
           del1  = del1 * f311 * g310 * q31 * aonv;
           xlamo = rem(mo + nodeo + argpo - theta, twopi);
           xfact = mdot + xpidot - rptim + dmdt + domdt + dnodt - no;
       end

       % /* ------------ for sgp4, initialize the integrator ---------- */
       xli   = xlamo;
       xni   = no;
       atime = 0.0;
       nm    = no + dndt;
   end
end









