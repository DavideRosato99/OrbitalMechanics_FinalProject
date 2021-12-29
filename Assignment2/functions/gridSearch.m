function [data] = gridSearch(data, settings)

%%
minDate = data.timeWindows.depDate;
maxDate = data.timeWindows.arrDate;
depPlanID   = data.starting.depPlanID;
flyByPlanID = data.starting.flyByPlanID;
arrPlanID   = data.starting.arrPlanID;

%%
minDatemjd2000 = date2mjd2000(minDate);
maxDatemjd2000 = date2mjd2000(maxDate);
minHfl = 200;
ctr = 1;
for t = 0:30:10*365
    t/365
    for TOF1 = 0:50:500
        for TOF2 = 0:50:500
            [kepDep, ksun] = uplanet(minDatemjd2000 + t, depPlanID);
            [rrDep, vvDep] = par2car(kepDep, ksun);

            [kepFlyBy, ksun] = uplanet(minDatemjd2000 + TOF1 + t, flyByPlanID);
            [rrFlyBy, vvFlyBy] = par2car(kepFlyBy, ksun);

            [kepArr, ksun] = uplanet(minDatemjd2000 + TOF1 + TOF2 + t, arrPlanID);
            [rrArr, vvArr] = par2car(kepArr, ksun);

            [a1, p1, e1, error1, vi1, vf1, Tpar1, theta1, a2, p2, e2, ...
                error2, vi2, vf2, Tpar2, theta2, DV1(ctr), DV2(ctr), errorFB, r_p, h_ga, delta, ...
                delta_V_powFB(ctr), e_minus, e_plus, a_minus, a_plus] = deltaVtot(...
                rrDep, rrFlyBy, rrArr, vvDep, vvFlyBy, vvArr, TOF1*24*3600, TOF2*24*3600, minHfl, flyByPlanID);
            DVTOT(ctr) = DV1(ctr) + DV2(ctr) + delta_V_powFB(ctr);
            ctr = ctr + 1;
        end
    end
end

[minDVTOT, ind] = min(DVTOT)
DV1(ind)
DV2(ind)
delta_V_powFB(ind)









