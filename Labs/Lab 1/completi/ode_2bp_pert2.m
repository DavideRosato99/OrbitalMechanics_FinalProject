function dx=ode_2bp_pert2(t, x, mu)

    J2 = 0.00108263;
    Re= 6378.137;
    
    r=norm(x);
    acc=3/2*(J2*mu*(Re^2)/(r^4))*[x(1)/r*(5*x(3)^2/r^2-1)
                                  x(2)/r*(5*x(3)^2/r^2-1)
                                  x(3)/r*(5*x(3)^2/r^2-3)];
    
    dx=[x(4) x(5) x(6) -mu/(r^3)*x(1)+acc(1) -mu/(r^3)*x(2)+acc(2) -mu/(r^3)*x(3)+acc(3)]';
    
end