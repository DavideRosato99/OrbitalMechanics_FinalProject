function dy = ode_2body( t, y, mu )



rr = y(1:3);
vv = y(4:6);

r = norm(rr);
v = norm(vv);

% Set the derivatives of the state
dy = [  vv(1)            ;
        vv(2)            ;
        vv(3)            ;
        -mu/r^3 *  rr(1) ;
        -mu/r^3 *  rr(2) ;
        -mu/r^3 *  rr(3)];

  
% Validity check:
%  epsilon = v^2/2 - mu/r;
%  hh = cross(rr, vv);                h = norm(hh);
%  ee = 1/mu * cross(vv, hh) - rr/r;  e = norm(ee);
% 
%  [epsilon h e]'
 
end

