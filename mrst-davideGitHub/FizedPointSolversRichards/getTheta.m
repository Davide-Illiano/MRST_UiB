function theta = getTheta(p, theta_R, theta_S, alpha, n)
    theta_neg = theta_R + (theta_S - theta_R).*(1./(1 + (-alpha*p).^n)).^((n-1)/n);
    theta_pos = theta_S;
    
    neg = p <= 0;
    theta = theta_neg.*neg + theta_pos.*~neg;    

end