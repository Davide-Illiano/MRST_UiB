function thetafirst = getThetafirst(p, theta_R, theta_S, alpha, n)

    thetafirst_neg = -(alpha.*(-alpha.*p).^(n - 1).*(theta_R - theta_S).*(1./((-alpha.*p).^n + 1)).^((n - 1)./n - 1).*(n - 1))./((-alpha.*p).^n + 1).^2;
    thetafirst_pos = 0;
    
    neg = p <= 0;
    thetafirst = thetafirst_neg.*neg + thetafirst_pos.*~neg;
end