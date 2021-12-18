function [deltapsi, meaneps, nut] = nutation(ttt, ddpsi, ddeps)

        deg2rad = pi/180.0;

        [iar80, rar80] = iau80in;  % coeff in deg

        % ---- determine coefficients for iau 1980 nutation theory ----
        ttt2= ttt*ttt;
        ttt3= ttt2*ttt;

        meaneps = -46.8150 *ttt - 0.00059 *ttt2 + 0.001813 *ttt3 + 84381.448;
        meaneps = rem( meaneps/3600.0, 360.0 );
        meaneps = meaneps * deg2rad;

        [l, l1, f, d, omega] = fundarg(ttt);
%fprintf(1,'nut del arg %11.7f  %11.7f  %11.7f  %11.7f  %11.7f  \n',l*180/pi,l1*180/pi,f*180/pi,d*180/pi,omega*180/pi );

        deltapsi= 0.0;
        deltaeps= 0.0;
        for i= 106:-1: 1
            tempval= iar80(i,1)*l + iar80(i,2)*l1 + iar80(i,3)*f + ...
                     iar80(i,4)*d + iar80(i,5)*omega;
            deltapsi= deltapsi + (rar80(i,1)+rar80(i,2)*ttt) * sin( tempval );
            deltaeps= deltaeps + (rar80(i,3)+rar80(i,4)*ttt) * cos( tempval );
        end

        % --------------- find nutation parameters --------------------
        deltapsi = rem( deltapsi + ddpsi, 2.0 * pi );
        deltaeps = rem( deltaeps + ddeps, 2.0 * pi );
        trueeps  = meaneps + deltaeps;

%fprintf(1,'meaneps %11.7f dp  %11.7f de  %11.7f te  %11.7f  ttt  %11.7f \n',meaneps*180/pi,deltapsi*180/pi,deltaeps*180/pi,trueeps*180/pi, ttt );

        cospsi  = cos(deltapsi);
        sinpsi  = sin(deltapsi);
        coseps  = cos(meaneps);
        sineps  = sin(meaneps);
        costrueeps = cos(trueeps);
        sintrueeps = sin(trueeps);

        nut(1,1) =  cospsi;
        nut(1,2) =  costrueeps * sinpsi;
        nut(1,3) =  sintrueeps * sinpsi;
        nut(2,1) = -coseps * sinpsi;
        nut(2,2) =  costrueeps * coseps * cospsi + sintrueeps * sineps;
        nut(2,3) =  sintrueeps * coseps * cospsi - sineps * costrueeps;
        nut(3,1) = -sineps * sinpsi;
        nut(3,2) =  costrueeps * sineps * cospsi - sintrueeps * coseps;
        nut(3,3) =  sintrueeps * sineps * cospsi + costrueeps * coseps;
        
end


function [l, l1, f, d, omega] = fundarg(ttt)
l     = ((((0.064) * ttt + 31.310) * ttt + 1717915922.6330) * ttt) / 3600.0 + 134.96298139;
l1    = ((((-0.012) * ttt - 0.577) * ttt + 129596581.2240) * ttt) / 3600.0 + 357.52772333;
f     = ((((0.011) * ttt - 13.257) * ttt + 1739527263.1370) * ttt) / 3600.0 + 93.27191028;
d     = ((((0.019) * ttt - 6.891) * ttt + 1602961601.3280) * ttt) / 3600.0 + 297.85036306;
omega = ((((0.008) * ttt + 7.455) * ttt - 6962890.5390) * ttt) / 3600.0 + 125.04452222;

end