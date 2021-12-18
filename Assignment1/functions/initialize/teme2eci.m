function [rECI, vECI] = teme2eci(rTEME, vTEME, timemjd2000)

if size(rTEME, 1) == 1
    rTEME = rTEME';
end

if size(vTEME, 1) == 1
    vTEME = vTEME';
end

ddpsi = -0.054522 * pi / (180*3600);
ddeps = -0.006209 * pi / (180*3600);
centjd = (timemjd2000)/(36525);

prec = precess(centjd);

[deltapsi, meaneps, nut] = nutation(centjd, ddpsi, ddeps);

eqeg = deltapsi* cos(meaneps);

eqeg = rem (eqeg, 2.0*pi);

eqe(1,1) =  cos(eqeg);
eqe(1,2) =  sin(eqeg);
eqe(1,3) =  0.0;
eqe(2,1) = -sin(eqeg);
eqe(2,2) =  cos(eqeg);
eqe(2,3) =  0.0;
eqe(3,1) =  0.0;
eqe(3,2) =  0.0;
eqe(3,3) =  1.0;

tm = prec * nut * eqe';

rECI = tm * rTEME;
vECI = tm * vTEME;