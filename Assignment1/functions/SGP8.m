function [pos, vel] = SGP8(year, fracDay, satdata)
% SGP8 - function to compute position and velocity vectors w.r.t. TEME
% reference frame (True Equator, Mean Equinox) starting from satellites
% TLEs
%
% PROTOTYPE
%   [pos, vel] = SGP8(year, fracDay, satdata)
%
% INPUT:
%   year     double [1x1]  Year at which calculations are performed  [-]
%   fracDay  double [1x1]  Fractional day in teh year                [days]
%   satData  table  [Nx16] TLEs table                                [-]
%
% OUTPUT:
%   pos      double [Nx3]  Position vector (TEME reference frame)    [km]
%   vel      double [Nx3]  Velocity vector (TEME reference frame)    [km/s]
%
% NOTE: Year must be greater or equal to 2000
%
% CONTRIBUTORS:
%   Rosato Davide               10618468
%   Saba Mohammadi Yengeje      10789462
%   Spinelli jason              10618465
%   Tagliati Alessia            10635119
%
% REFERENCES
%   [1] Celestrak documentation, NORAD,
%       URL: https://celestrak.com/NORAD/documentation/spacetrk.pdf
%
% VERSIONS
%   2021-10-21: Release
%
% -------------------------------------------------------------------------


%% PRE-CALCULATION ALLOCATION
N = size(satdata, 1);
pos = zeros(N, 3);
vel = zeros(N, 3);

%% CONSTANTS
ae     = 1;
xmnpda = 1440;
tothrd = 2/3;
XJ3    = -2.53881e-6;
e6a    = 1e-6;
xkmper = astroConstants(23);
ge     = astroConstants(13);
CK2    = (1.0826158e-3 / 2);
CK4    = (-3 * (-1.65597e-6/8) );

%% ACTUAL EPOCH
ActYear = mod(year, 100);

%% LOOP FOR EACH TLE
for ind = 1:N
    TLEyear = (satData.epoch(ind) - mod(satData.epoch(ind), 1000)) / 1000;
    yearSince = ActYear - TLEyear;
    TLEday = mod(satData.epoch(ind), 1000);
    daySince = fracDay - TLEday;
    tsince = yearSince*365*24*60 + daySince*24*60;
    
    flag = false;
    
    %%% LOOP CONSTANTS
    rho    = (2.461 * 0.00001 * xkmper);
    s      = ae + 78/xkmper;
    qo     = ae + 120/xkmper;
    xke    = sqrt((3600.0 * ge)/(xkmper^3));
    qoms2t = ((qo - s)^2)^2;
    temp2  = xke / (satdata.xn0(ind));
    a1     = (temp2^tothrd);
    cosio  = cos(satdata.xincl(ind));
    sinio  = sin(satdata.xincl(ind));
    delta1 = 1.5 * CK2 * (3.0 * cosio * cosio - 1.0) / ((a1^2) * ((1.0 - (satdata.e0(ind)^2))^1.5));
    a0     = a1 * (1.0 - (1.0/3.0) * delta1 - (delta1^2) - (134.0/81.0) * (delta1^2) * delta1);
    delta0 = 1.5 * CK2 * (3.0 * cosio * cosio - 1.0) / ((a0^2) * ((1.0 - (satdata.e0(ind)^2))^1.5));
    noii   = satdata.xn0(ind) / (1.0 + delta0);
    aoii   = a0 / (1.0 - delta0);

    %%% B TERM
    B = 2 * (satdata.bstar(ind))/rho;
    
    %%% LOOP CONSTANTS II
    beta = sqrt(1.0 - (satdata.e0(ind)^2));
    theta = cosio;
    temp = (noii * CK2) / ((aoii^2) * (beta^2) * beta);
    M1dot = -1.5 * temp * (1.0 - 3.0 * (theta^2));
    omega1dot = -1.5 * (temp / beta) * (1.0 - 5.0*(theta^2));
    xnode1dot = -3.0 * (temp / beta) * theta;
    M2dot = (3.0/16.0) * (temp * CK2 / ((aoii^2) * (beta^4.0))) * (13.0 - 78.0 * (theta^2) + 137.0 * (theta^4.0));
    omega2dot = (3.0/16.0) * (temp * CK2 / ((aoii^2) * (beta^5.0))) * (7.0 - 114.0 *(theta^2)+ 395.0 * (theta^4)) + 1.25 * (noii * CK4) / ((aoii^4.0) * (beta^8.0))* (3.0 - 36.0 * (theta^2) + 49.0 * (theta^4));
    xnode2dot = 1.5 * (temp * CK2 / ((aoii^2) * (beta^5.0))) * theta * (4.0 - 19.0 * (theta^2))+ 2.5 * (noii * CK4) / ((aoii^4.0) * (beta^8.0)) * theta * (3.0 - 7.0*(theta^2));
    ldot = noii + M1dot + M2dot;
    omegadot = omega1dot + omega2dot;
    xnodedot = xnode1dot + xnode2dot;
    tsi = 1.0/(aoii * (beta^2) - s);
    eta = satdata.e0(ind) * s * tsi;
    phi = sqrt(1.0 - (eta^2));
    alpha = sqrt(1.0 + (satdata.e0(ind)^2));
    C0 = 0.5 * B * rho * qoms2t * noii * aoii * (tsi^4.0) / (alpha * (phi^7.0));
    C1 = 1.5 * noii * (alpha^4.0) * C0;
    D1 = (tsi * (1.0 / (phi^2))) / (aoii * (beta^2));
    D2 = 12.0 + 36.0 * (eta^2) + 4.5 * (eta^4.0);
    D3 = 15.0 * (eta^2) + 2.5 * (eta^4.0);
    D4 = 5.0 * eta + (15.0 / 4.0) * (eta^3.0);
    D5 = tsi / (phi^2);
    B1 = -CK2 * (1.0 - 3.0 * (theta^2));
    B2 = -CK2 * (1.0 - (theta^2));
    A30 = -XJ3 * (ae^3.0) / CK2;
    B3 = A30 * sinio;
    C2 = D1 * D3 * B2;
    C3 = D4 * D5 * B3;
    n0dot = C1 * ((2.0 + (eta^2) * (3.0 + 34.0 * (satdata.e0(ind)^2)) + 5.0 * (satdata.e0(ind)) * (eta) * (4.0 + (eta^2)) + 8.5 * (satdata.e0(ind)^2)) + D1 * D2 * B1 + C2 * cos(2.0 * satdata.omega0(ind)) + C3 * sin(satdata.omega0(ind)));
    D6 = 30 * eta + 22.5 * (eta^3);
    D7 = 5.0 * eta + 12.5 * (eta^3);
    D8 = 1.0 + 6.75 * (eta^2) + (eta^4);
    C4 = D1 * D7 * B2;
    C5 = D5 * D8 * B3;
    e0dot = -C0 * (eta * (4.0 + (eta^2) + (satdata.e0(ind)^2) * (15.5 + 7.0 * (eta^2))) + satdata.e0(ind) * (5.0 + 15.0 * (eta^2)) + D1 * D6 * B1 + C4 * cos(2.0 * satdata.omega0(ind)) + C5 * sin(satdata.omega0(ind)));
    alphadottoalpha = satdata.e0(ind) * e0dot / ((alpha^2));
    C6 = (1.0/3.0) * (n0dot/noii);
    tsidottotsi = 2.0 * aoii * tsi * (C6 * (beta^2) + satdata.e0(ind) * e0dot);
    etadot = (e0dot + satdata.e0(ind) * tsidottotsi) * s * tsi;
    phidottophi = -eta * etadot / (phi^2);
    C0dottoC0 = C6 + 4.0 * tsidottotsi - alphadottoalpha - 7.0 * phidottophi;
    C1dottoC1 = n0dot / noii + 4.0 * alphadottoalpha + C0dottoC0;
    D9 = 6.0 * eta + 20.0 * satdata.e0(ind) + 15.0 * satdata.e0(ind) * (eta^2) + 68.0 * (satdata.e0(ind)^2) * eta;
    D10 = 20.0 * eta + 5.0 * (eta^3) + 17.0 * satdata.e0(ind) + 68 * (satdata.e0(ind)) * (eta^2);
    D11 = eta * (72.0 + 18.0 * (eta^2));
    D12 = eta * (30.0 + 10.0 * (eta^2));
    D13 = 5.0 + 11.25*(eta^2);
    D14 = tsidottotsi - 2.0 * phidottophi;
    D15 = 2.0 * (C6 + satdata.e0(ind) * e0dot / (beta^2));
    D1dot = D1 * (D14 + D15);
    D2dot = etadot * D11;
    D3dot = etadot * D12;
    D4dot = etadot * D13;
    D5dot = D5 * D14;
    C2dot = B2 * (D1dot * D3 + D1 * D3dot);
    C3dot = B3 * (D5dot * D4 + D5 * D4dot);
    D16 = D9 * etadot + D10 * e0dot + B1 * (D1dot * D2 + D1 * D2dot) + C2dot * cos(2.0 * satdata.omega0(ind)) + C3dot * sin(satdata.omega0(ind))+omega1dot * (C3 * cos(satdata.omega0(ind)) - 2 * C2 * sin(2.0 * satdata.omega0(ind)));
    n0ddot = n0dot * C1dottoC1 + C1 * D16;
    e0ddot = e0dot * C0dottoC0 - C0 * ((4.0 + 3.0 * (eta^2) + 30.0 * satdata.e0(ind) * eta + 15.5 * (satdata.e0(ind)^2) + 21.0 * (eta^2) * (satdata.e0(ind)^2)) * etadot + (5.0 + 15.0 * (eta^2) + 31.0 * satdata.e0(ind) * eta + 14.0 * satdata.e0(ind) * (eta^3)) * e0dot + B1 * (D1dot * D6 + D1 * etadot * (30.0 + 67.5 * (eta^2)))+ B2 * (D1dot * D7 + D1 * etadot * (5.0 + 37.5 * (eta^2))) * cos(2.0 * satdata.omega0(ind))+ B3 * (D5dot * D8 + D5 * eta * etadot * (13.5 + 4 * (eta^2))) * sin(satdata.omega0(ind)) + omega1dot * (C5 * cos(satdata.omega0(ind)) - 2.0 * C4 * sin(2.0 * satdata.omega0(ind))));
    D17 = n0ddot / noii - (n0dot/noii)^2;
    tsiddottotsi = 2.0 * (tsidottotsi - C6) *tsidottotsi + 2.0 * aoii * tsi * ((1.0/3.0) * D17 * (beta^2) - 2 * C6 * satdata.e0(ind) * e0dot + (e0dot^2) + satdata.e0(ind) * e0ddot);
    etaddot = (e0ddot + 2.0 * e0dot * tsidottotsi) * s * tsi + eta * tsiddottotsi;
    D18 = tsiddottotsi - (tsidottotsi^2);
    D19 = -(phidottophi^2) * (1.0 + 1.0 / (eta^2)) - eta * etaddot / (phi^2);
    D1ddot = D1dot * (D14+D15) + D1 * (D18 - 2.0 * D19 + (2.0/3.0) * D17 + 2.0 * (alpha^2) * (e0dot^2) / (beta^4.0) + 2.0 * satdata.e0(ind) * e0ddot / (beta^2));
    n0tdot = n0dot * (2.0 * (2.0/3.0) * D17 + 3.0 * ((e0dot^2) + satdata.e0(ind) * e0ddot) / (alpha^2) - 6.0 * ( satdata.e0(ind) * e0dot / (alpha^2))^2 + 4.0 * D18 - 7.0 * D19 ) + C1dottoC1 * n0ddot;
    n0tdot = n0tdot + C1 * (C1dottoC1 * D16 + D9 * etaddot + D10 * e0ddot + (etadot^2) * (6.0 + 30.0 * (satdata.e0(ind)) * eta + 68.0 * (satdata.e0(ind)^2)) + etadot * e0dot * (40.0 + 30.0 * (eta^2) + 272.0 * satdata.e0(ind) * eta) + (e0dot^2) * (17.0 + 68.0 * (eta^2)) + B1 * (D1ddot * D2+ 2.0 * D1dot * D2dot + D1 * (etaddot * D11 + (etadot^2) * (72. + 54. * (eta^2))))+ B2 * (D1ddot * D3 + 2.0 * D1dot * D3dot + D1 * (etaddot * D12 + (etadot^2) * (30.0 + 30.0 * (eta^2)))) * cos(2.0 * satdata.omega0(ind)) + B3 * ((D5dot * D14 + D5 * (D18 - 2.0 * D19)) * D4 + 2.0 * D4dot * D5dot + D5 * (etaddot * D13 + 22.5 * eta * (etadot^2))) * sin(satdata.omega0(ind)) + omega1dot * ((7.0 * (1.0/3.0) * (n0dot / noii) + 4.0 * satdata.e0(ind) * e0dot / (beta^2)) * (C3 * cos(satdata.omega0(ind)) - 2.0 * C2 * sin(2.0 * satdata.omega0(ind)))+ ((2.0 * C3dot * cos(satdata.omega0(ind)) - 4.0 * C2dot * sin(2.0 * satdata.omega0(ind)) - omega1dot * (C3 * sin(satdata.omega0(ind)) + 4 * C2 * cos(2.0 * satdata.omega0(ind)))))));
    smalln0ddot = n0ddot * 1.0E9;
    temp = (smalln0ddot^2) - n0dot * 1.0E18 * n0tdot;
    p = (temp + (smalln0ddot^2)) / temp;
    gamma = - n0tdot / (n0ddot * (p - 2));
    nd = n0dot / (p * gamma);
    q = 1.0 - e0ddot / (e0dot * gamma);
    ed = e0dot / (q * gamma);
    
    %%% ATMOSPHERIC DRAG AND GRAVITATION
    if ((abs( n0dot / noii * xmnpda)) < 0.00216)
        n = noii + n0dot * tsince;
        edot = -(2.0/3.0) * n0dot * (1.0 - satdata.e0(ind)) / noii;
        e = satdata.e0(ind) + edot * tsince;
        Z1 = 0.5 * n0dot * (tsince^2);
    else
        n = noii + nd * (1.0 - ((1.0 - gamma * tsince)^p));
        edot = e0dot;
        e = satdata.e0(ind) + ed * (1.0 - ((1.0 - gamma * tsince)^q));
        Z1 = gamma * tsince;
        Z1 = 1.0 - Z1;
        if Z1 < 0
            flag = true;
            pos(ind, 1:3) = nan(1, 3);
            vel(ind, 1:3) = nan(1, 3);
        else
            Z1 = (Z1^(p+1.0));
            Z1 = Z1 - 1.0;
            Z1 = Z1 / (gamma * (p+1));
            Z1 = Z1 + tsince;
            Z1 = Z1 * n0dot / (p * gamma);
        end
    end
    
    if not(flag)
        omega = satdata.omega0(ind) + omegadot * tsince + (7.0 * Z1) / (3.0 * noii) * omega1dot;
        xnode = satdata.xnode0(ind) + xnodedot * tsince + xnode1dot * (7.0 * Z1)/(3.0 * noii);
        M = satdata.xm0(ind) + tsince * ldot + Z1 + M1dot * (7.0 * Z1) / (3.0 * noii);

        %%% KEPLE R
        M    = fmod2p(M);
        temp = M + satdata.e0(ind) * sin(M) * (1 + e * cos(M));
        i    = 1;

        while 1
            deltaEw = (M + e * sin(temp) - temp) / (1.0 - e * cos(temp));
            Ew = deltaEw + temp;
            temp = Ew;
            i=i+1;

            if ( (i>10) || (abs(deltaEw) <= e6a) )
                break
            end
        end

        %%% SHORT-PERIOD PERIODICS
        a = ((xke / n)^tothrd);
        beta = sqrt(1.0 - (e^2));
        sinf = beta * sin(Ew) * (1.0 / (1.0 - e * cos(Ew)));
        cosf = (cos(Ew) - e) * (1.0 / (1.0 - e * cos(Ew)));
        f = actan(sinf, cosf);
        sinu = sinf * cos(omega) + cosf * sin(omega);
        cosu = cosf * cos(omega) - sinf * sin(omega);
        rii = (a * (1.0 - (e^2))) / (1.0 + (e * cosf));
        deltaR = (0.5 * CK2 * (1.0 / (a * (1.0 - (e^2))))) * ((1.0 - (theta^2)) * (2.0 * (cosu^2) - 1.0) - 3.0 * (3.0 * (theta^2) - 1.0)) - (0.25 * A30 * sinio) * sinu;
        temp = 3.0 * (0.5 * CK2 * (1.0 / (a * (1.0 - (e^2))))^2) * sinio * (2.0 * (cosu^2) - 1) - (0.25 * A30 * (1.0 / (a * (1 - (e^2))))) * e * sin(omega);
        deltaI = temp * cosio;
        deltaU = sin(satdata.xincl(ind) / 2) * ((0.5 * CK2 * (1.0 / (a * (1 - (e^2))))^2) * (0.5 * (1.0 - 7.0 * (theta^2)) * (2.0 * sinu * cosu)- 3.0 * (1.0 - 5.0 * (theta^2)) * (f - M + e * sinf)) - (0.25 * A30 * (1.0 / (a * (1.0 - (e^2))))) * sinio * cosu * (2.0 + (e*cosf))) - 0.5 * (0.25 * A30 * (1.0 / (a * (1 - (e^2))))) * (theta^2) * e * cos(omega) / cos(satdata.xincl(ind)/2);
        lambda = f + omega + xnode + (0.5 * CK2 * (1.0 / (a * (1.0 - (e^2))))^2) * (0.5 * (1.0 + 6.0 * cosio - 7.0 * (theta^2)) * (2.0 * sinu * cosu) - 3.0 * ((1.0 - 5.0 * (theta^2)) + 2.0 * cosio) * (f - M + e * sinf )) + (0.25 * A30 * (1.0 / (a *(1 - (e^2))))) * sinio * (cosio * e * cos(omega) / (1.0 + cosio) - (2.0 + (e*cosf)) * cosu);
        y4 = sin(satdata.xincl(ind) / 2.0) * sinu + cosu * deltaU + 0.5 * sinu * cos(satdata.xincl(ind) / 2.0)*deltaI;
        y5 = sin(satdata.xincl(ind) / 2.0) * cosu - sinu * deltaU + 0.5 * cosu * cos(satdata.xincl(ind) / 2.0)*deltaI;
        r = rii + deltaR;
        rdot = n * a * e * sinf / beta + (-n * (a / rii)^2) * (1.0 * CK2 * (1.0 / (a * (1.0 - (e^2)))) * (1.0 - (theta^2)) * (2.0 * sinu * cosu) + (0.25 * A30 * sinio) * cosu);
        rfdot = n * (a^2) * beta / rii + (-n * (a/rii)^2) * deltaR + a * (n * (a / rii)) * sinio * temp;


        %%% UNIT-ORIENTATION VECTOR
        cosI2 = sqrt(1.0 - (y4^2) - (y5^2));
        UVv(1) = 2.0 * y4 * (y5 * sin(lambda) - y4 * cos(lambda)) + cos(lambda);
        UVv(2) = -2.0 * y4 * (y5 * cos(lambda) + y4 * sin(lambda)) + sin(lambda);
        UVv(3) = 2.0 * y4 *cosI2;
        VVv(1) = 2.0 * y5 * (y5 * sin(lambda) - y4 * cos(lambda)) - sin(lambda);
        VVv(2) = -2.0 * y5 * (y5 * cos(lambda) + y4 * sin(lambda)) + cos(lambda);
        VVv(3) = 2.0 * y5 * cosI2;

        %%% POSITION AND VELOCITY
        posI(1:3) = r * UVv;
        velI(1:3) = rdot * UVv + rfdot * VVv;

        [pos(ind, 1:3), vel(ind, 1:3)] = ConvertSatState(posI, velI);
    end
    
    
end

end


function x = fmod2p(x)
    x = modulus(x, 2*pi);
end

function modu = modulus(arg1, arg2)
    modu = arg1 - floor(arg1/arg2) * arg2;
    
    if modu >= 0
        return
    else
        modu = modu + arg2;
        return
    end
end

function [p, v] = ConvertSatState(pos, vel)
    xkmper = astroConstants(23);
    
    p = pos * xkmper;
    v = vel * xkmper / 60;
end

function t = actan(y, x)
    t = atan2(y,x);
    if (t < 0)
        t = t + 2*pi;
    end
end