function [VV_minus, VV_plus] = flybyPlotHelio(vv_inf_minus, VV_p, Delta, rr_s2p , mu_s, mu_p, inc_asymptote_angularPos, handle)


    % Preraring Figure Object
    if nargin<8
        HAXIS = gca;
    elseif ishandle(handle)==0
            msg = ['The figure handle is not valid'];
            error(msg)
    else
        try
            HAXIS=gca(handle);
        catch
            HAXIS=handle;  
        end
        hold on
    end

v_inf = norm(vv_inf_minus);

a = - mu_p/v_inf^2;
delta= 2*atan(-a/Delta);
e = 1/sin(delta/2);
r_p = a * (1-e);
Delta_V = 2*v_inf*sin(delta/2);


VV_minus = VV_p + vv_inf_minus;

norm2helioPlane1 = cross(rr_s2p, VV_minus);
norm2helioPlane1 = norm2helioPlane1/norm(norm2helioPlane1);

    switch inc_asymptote_angularPos
        case 'leading'
            rot_dir = -norm2helioPlane1;
        case 'trailing'
            rot_dir = +norm2helioPlane1;
        case 'under'
            rot_dir = - VV_p;
        case 'over'
            rot_dir = + VV_p;      
        otherwise
            theta = inc_asymptote_angularPos; % [rad]
            % counterclockwise angle that defines the location
            % of incoming asymptote, considering a rotation around
            % the position vector and s.t. theta=0 -> leading-side
            rot_dir = rotVecAroundVecByAngle(norm2helioPlane1, rr_s2p, -theta);  
    end
rot_dir = rot_dir/norm(rot_dir);

vv_inf_plus = rotVecAroundVecByAngle(vv_inf_minus, rot_dir, delta);
VV_plus = VV_p + vv_inf_plus;



norm2helioPlane2 = cross(rr_s2p, VV_plus);
norm2helioPlane2 = norm2helioPlane2/norm(norm2helioPlane2);

hold on
grid on
axis equal
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
dt = 10000;
tspan = [-200*86400:dt:0];
y0 = [rr_s2p; -VV_minus];
[time_vec, StateMat] = ode113( @(t,y) ode_2body(t,y, mu_s), tspan, y0, options );
X_arr = StateMat(:,1); Y_arr = StateMat(:,2); Z_arr = StateMat(:,3);
%VX = StateMat_arr(:,4); VY = StateMat_arr(:,5); VZ = StateMat_arr(:,6);
plot3(HAXIS, X_arr,Y_arr,Z_arr, 'LineWidth', 1.5, 'Color', 1/255*[66, 134, 244])

tspan = [0:dt:200*86400];
y0 = [rr_s2p; VV_plus];
[time_vec, StateMat] = ode113( @(t,y) ode_2body(t,y, mu_s), tspan, y0, options );
X_arr = StateMat(:,1); Y_arr = StateMat(:,2); Z_arr = StateMat(:,3);
%VX = StateMat_arr(:,4); VY = StateMat_arr(:,5); VZ = StateMat_arr(:,6);
plot3(HAXIS, X_arr,Y_arr,Z_arr, 'LineWidth', 1.5, 'Color', 1/255*[193, 0, 0])


sunGlobe = plotPlanet(10, [0,0,0], HAXIS, 30);






end

