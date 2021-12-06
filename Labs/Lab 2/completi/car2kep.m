function [a,e,ideg,Omdeg,omdeg,thetadeg] = car2kep(rr,vv,   mu)
% To obtain the results in radiants remove deg from outputs in the function

r = norm(rr);
v = norm(vv);
k = [0 0 1]';

a = -mu / (v^2 - 2*mu/r);

hh = cross(rr,vv);
h = norm(hh);

ee = cross(vv,hh)/mu - rr/r;
e = norm(ee);

i = acos(hh(3)/h);

N = cross(k,hh) / norm(cross(hh,k));

Om = acos(N(1));
if N(2) < 0
    Om = 2*pi - Om;
end

om = acos(dot(N,ee)/e);
if ee(3) < 0
    om = 2*pi - om;
end

theta = acos(dot(rr,ee)/(r*e));
if dot(rr,vv) < 0
    theta = 2*pi - theta;
end

ideg=rad2deg(i);
omdeg=rad2deg(om);
Omdeg=rad2deg(Om);
thetadeg=rad2deg(theta);


end

