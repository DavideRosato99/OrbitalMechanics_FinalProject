function dx=ode_2bp(t, x, mu)

    %Model: Unperturbed restricted 2bp
    %dx -> State vector derivative - in column, components of velocity and
    %acceleration vectors
    %x -> State vector - in column, components of position and velocity
    %vectors
    %mu -> mu 
    %r -> Distance from celestial body

    r=norm(x);
    dx=[x(4:6); -mu/(r^3)*x(1:3)];
end