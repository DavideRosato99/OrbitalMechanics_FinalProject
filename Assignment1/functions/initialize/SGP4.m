function [rr, vv] = SGP4(satData, date)

%% INITIALIZE
N = size(satData.Name, 2);
rr = zeros(N, 3);
vv = zeros(N, 3);

%% SET MATHEMATICAL CONSTANTS
twopi = 2.0 * pi;
x2o3  = 2.0 / 3.0;

% Divide by zero check
temp4    =   1.5e-12;

%% TLEs loop
for i = 1:N
    vkmpersec = satData.radiusearthkm(i) * satData.xke(i) / 60;
    
    %%% T since
    Ytle = satData.epochyr(i) + 2000;
    Dtle = satData.epochdays(i);
    [mon, day, hr, minute, sec] = days2mdh(Ytle, Dtle);
    
    daysSince = datenum(date) - datenum([Ytle, mon, day, hr, minute, sec]);
    tsince = daysSince * 24 * 60;

    % Update for secular gravity and atmospheric drag
    xmdf    = satData.mo(i) + satData.mdot(i) * tsince;
    argpdf  = satData.argpo(i) + satData.argpdot(i) * tsince;
    nodedf  = satData.nodeo(i) + satData.nodedot(i) * tsince;
    argpm   = argpdf;
    mm      = xmdf;
    t2      = tsince * tsince;
    nodem   = nodedf + satData.nodecf(i) * t2;
    tempa   = 1.0 - satData.cc1(i) * tsince;
    tempe   = satData.bstar(i) * satData.cc4(i) * tsince;
    templ   = satData.t2cof(i) * t2;
    
    if (satData.isimp ~= 1)
        delomg = satData.omgcof(i) * tsince;
        delm   = satData.xmcof(i) * ((1.0 + satData.eta(i) * cos(xmdf))^3 ...
            - satData.delmo(i));
        temp   = delomg + delm;
        mm     = xmdf + temp;
        argpm  = argpdf - temp;
        t3     = t2 * tsince;
        t4     = t3 * tsince;
        tempa  = tempa - satData.d2(i) * t2 - satData.d3(i) * t3 -...
            satData.d4(i) * t4;
        tempe  = tempe + satData.bstar(i) * satData.cc5(i) * (sin(mm) -...
            satData.sinmao(i));
        templ  = templ + satData.t3cof(i) * t3 + t4 * (satData.t4cof(i) +...
            tsince * satData.t5cof(i));
    end
    
    nm    = satData.no(i);
    em    = satData.ecco(i);
    inclm = satData.inclo(i);
    
    if (satData.method(i) == 'd')
        tc = tsince;
        [satData.atime(i),em,argpm,inclm,satData.xli(i),mm,...
            satData.xni(i),nodem,dndt,nm] = dspace(...
            satData.d2201(i),satData.d2211(i),satData.d3210(i),...
            satData.d3222(i),satData.d4410(i),satData.d4422(i),...
            satData.d5220(i),satData.d5232(i),satData.d5421(i),...
            satData.d5433(i),satData.dedt(i),satData.del1(i),...
            satData.del2(i),satData.del3(i),satData.didt(i),...
            satData.dmdt(i),satData.dnodt(i),satData.domdt(i),...
            satData.irez(i),satData.argpo(i),satData.argpdot(i),tsince,...
            tc,satData.gsto(i),satData.xfact(i),satData.xlamo(i),satData.no(i),...
            satData.atime(i),em,argpm,inclm,satData.xli(i),mm,...
            satData.xni(i),nodem,nm);
    end
    
    if nm <= 0
        satData.error(i) = 2;
    end
    
    am = (satData.xke(i) / nm)^x2o3 * tempa * tempa;
    nm = satData.xke(i) / am^1.5;
    em = em - tempe;
    
    if ((em >= 1.0) || (em < -0.001) || (am < 0.95))
        satData.error(i) = 1;
    end 
    if (em < 1.0e-6)
        em  = 1.0e-6;
    end
    
    mm     = mm + satData.no(i) * templ;
    xlm    = mm + argpm + nodem;
    emsq   = em * em;
    temp   = 1.0 - emsq;
    nodem  = rem(nodem, twopi);
    argpm  = rem(argpm, twopi);
    xlm    = rem(xlm, twopi);
    mm     = rem(xlm - argpm - nodem, twopi);
    
    %%% COMPUTE EXTRA MEAN QUANTITIES
    sinim = sin(inclm);
    cosim = cos(inclm);
    
    %%% ADD LUNAR-SOLAR PERIODICS
    ep     = em;
    xincp  = inclm;
    argpp  = argpm;
    nodep  = nodem;
    mp     = mm;
    sinip  = sinim;
    cosip  = cosim;
    
    
    if (satData.method(i) == 'd')
        [ep,xincp,nodep,argpp,mp] = dpper(...
        satData.e3(i),satData.ee2(i),satData.peo(i),...
        satData.pgho(i),satData.pho(i),satData.pinco(i),...
        satData.plo(i),satData.se2(i),satData.se3(i),...
        satData.sgh2(i),satData.sgh3(i),satData.sgh4(i),...
        satData.sh2(i),satData.sh3(i),satData.si2(i),...
        satData.si3(i),satData.sl2(i),satData.sl3(i),...
        satData.sl4(i),tsince,satData.xgh2(i),...
        satData.xgh3(i),satData.xgh4(i),satData.xh2(i),...
        satData.xh3(i),satData.xi2(i),satData.xi3(i),...
        satData.xl2(i),satData.xl3(i),satData.xl4(i),...
        satData.zmol(i),satData.zmos(i),satData.inclo(i),...
        satData.init(i),ep,xincp,nodep,argpp,mp);
        if (xincp < 0.0)
            xincp  = -xincp;
            nodep = nodep + pi;
            argpp  = argpp - pi;
        end
        if ((ep < 0.0 ) || ( ep > 1.0))
            satData.error(i) = 3;
        end
    end
    
    %%% LONG PERIOD PERIODICS
    if (satData.method(i) == 'd')
        sinip =  sin(xincp);
        cosip =  cos(xincp);
        satData.aycof(i) = -0.5 * satData.j3oj2(i) * sinip;

        if (abs(cosip+1.0) > 1.5e-12)
            satData.xlcof(i) = -0.25 * satData.j3oj2(i) * sinip * (3.0 + 5.0 * cosip) / (1.0+cosip);
        else
            satData.xlcof(i) = -0.25 * satData.j3oj2(i) * sinip * (3.0 + 5.0 * cosip) / temp4;
        end
    end
    axnl = ep * cos(argpp);
    temp = 1.0 / (am * (1.0 - ep * ep));
    aynl = ep* sin(argpp) + temp * satData.aycof(i);
    xl   = mp + argpp + nodep + temp * satData.xlcof(i) * axnl;
    
    %%% SOLVE KEPLER EQUATION
    u    = rem(xl - nodep, twopi);
    eo1  = u;
    tem5 = 9999.9;
    ktr = 1;

    while (( abs(tem5) >= 1.0e-12) && (ktr <= 10) )
        sineo1 = sin(eo1);
        coseo1 = cos(eo1);
        tem5   = 1.0 - coseo1 * axnl - sineo1 * aynl;
        tem5   = (u - aynl * coseo1 + axnl * sineo1 - eo1) / tem5;
        if(abs(tem5) >= 0.95)
            if tem5 > 0.0
                tem5 = 0.95;
            else
                tem5 = -0.95;
            end
        end
        eo1    = eo1 + tem5;
        ktr = ktr + 1;
    end
    
    %%% SHORT PERIOD PRELIMINARY QUANTITIES
    ecose = axnl*coseo1 + aynl*sineo1;
    esine = axnl*sineo1 - aynl*coseo1;
    el2   = axnl*axnl + aynl*aynl;
    pl    = am*(1.0-el2);
    
    if pl < 0.0
        satData.error(i) = 4;
        rr(i, :) = [NaN NaN NaN];
        vv(i, :) = [NaN NaN NaN];
    else
        rl     = am * (1.0 - ecose);
        rdotl  = sqrt(am) * esine/rl;
        rvdotl = sqrt(pl) / rl;
        betal  = sqrt(1.0 - el2);
        temp   = esine / (1.0 + betal);
        sinu   = am / rl * (sineo1 - aynl - axnl * temp);
        cosu   = am / rl * (coseo1 - axnl + aynl * temp);
        su     = atan2(sinu, cosu);
        sin2u  = (cosu + cosu) * sinu;
        cos2u  = 1.0 - 2.0 * sinu * sinu;
        temp   = 1.0 / pl;
        temp1  = 0.5 * satData.j2(i) * temp;
        temp2  = temp1 * temp;

        % /* -------------- update for short period periodics ------------ */
        if (satData.method(i) == 'd')
            cosisq                 = cosip * cosip;
            satData.con41(i)  = 3.0*cosisq - 1.0;
            satData.x1mth2(i) = 1.0 - cosisq;
            satData.x7thm1(i) = 7.0*cosisq - 1.0;
        end
        mrt   = rl * (1.0 - 1.5 * temp2 * betal * satData.con41(i)) +...
        0.5 * temp1 * satData.x1mth2(i) * cos2u;
        su    = su - 0.25 * temp2 * satData.x7thm1(i) * sin2u;
        xnode = nodep + 1.5 * temp2 * cosip * sin2u;
        xinc  = xincp + 1.5 * temp2 * cosip * sinip * cos2u;
        mvt   = rdotl - nm * temp1 * satData.x1mth2(i) * sin2u / satData.xke(i);
        rvdot = rvdotl + nm * temp1 * (satData.x1mth2(i) * cos2u +...
        1.5 * satData.con41(i)) / satData.xke(i);

        % /* --------------------- orientation vectors ------------------- */
        sinsu =  sin(su);
        cossu =  cos(su);
        snod  =  sin(xnode);
        cnod  =  cos(xnode);
        sini  =  sin(xinc);
        cosi  =  cos(xinc);
        xmx   = -snod * cosi;
        xmy   =  cnod * cosi;
        ux    =  xmx * sinsu + cnod * cossu;
        uy    =  xmy * sinsu + snod * cossu;
        uz    =  sini * sinsu;
        vx    =  xmx * cossu - cnod * sinsu;
        vy    =  xmy * cossu - snod * sinsu;
        vz    =  sini * cossu;

        % /* --------- position and velocity (in km and km/sec) ---------- */
        rr(i, 1) = (mrt * ux)* satData.radiusearthkm(i);
        rr(i, 2) = (mrt * uy)* satData.radiusearthkm(i);
        rr(i, 3) = (mrt * uz)* satData.radiusearthkm(i);
        vv(i, 1) = (mvt * ux + rvdot * vx) * vkmpersec;
        vv(i, 2) = (mvt * uy + rvdot * vy) * vkmpersec;
        vv(i, 3) = (mvt * uz + rvdot * vz) * vkmpersec;
    end
    
    
end






