tDepIn= date2mjd2000([2003 04 01 12 00 00]);
tDepFin= date2mjd2000([2003 08 01 12 00 00]);
tDepValues= linspace(tDepIn, tDepFin,1000);
muSun= astroConstants(14);
kepDep=uplanet( tDepValues(k1),3);
sDep=kep2car(kepDep,muSun);
kepArr=uplanet( tArrValues(k2),3);
sDep=kep2car(kepArr,muSun);
ToF= (tArrValues(k2)-tDepValues(k1))*24*3600;


[C,h]=contour(t1_vec_MJD2000+Dtime,t2_vec_MJD2000+Dtime,Dvtot',.. floor(minDV_GS)+(0:1:5));
clabel(C,h,floor(minDV_GS)+(0:1:5));
datethick ('x','yyyy mm dd' 'keeplimits'=;
Dtime= datenum([2000 01 01 12 00 00])