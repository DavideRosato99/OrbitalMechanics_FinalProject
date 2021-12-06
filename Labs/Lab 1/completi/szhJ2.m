function acc=szhJ2(rv)

J2 = 0.00108263;
Re= 6378.137;
mu=astroConstants(13);
r=sqrt(rv(1)^2+rv(2)^2+rv(3)^2);
acc=3/2*(J2*mu*Re/r^4)*[rv(1)/r*(5*rv(3)^2/r^2-1) rv(2)/r*(5*rv(3)^2/r^2-1) rv(3)/r*(rv(3)^2/r^2 -3)];

end