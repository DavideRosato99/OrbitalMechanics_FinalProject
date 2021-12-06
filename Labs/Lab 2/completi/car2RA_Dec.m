function [alpha, delta, l, m, n] = car2RA_Dec(rr_mat)

xx = rr_mat(:,1);  yy = rr_mat(:,2);  zz = rr_mat(:,3);
rr = vecnorm(rr_mat')'; % vecnorm() calculates column-wise norm, but we need it as row-wise

l = xx./rr;     m = yy./rr;    n = zz./rr;

delta = asin(n);

alpha = acos(l./cos(delta));
if m <= 0
    alpha = 2*pi - alpha;
end


end

