function theta = getThetaCoupled(p, c, theta_r, theta_s, alpha, nGM, ng, a, b)

    %theta_res+(theta_sat-theta_res)*(1./(1+abs(alpha*(1./(1-bGM*log(c/aGM+1)).*p).^nGM))).^ng;

    % Nicu's Phi
     Phi = (1./(1+real((-alpha*(1./(1-b*log(c/a+1)).*p)).^ng))).^nGM;
    %     (1./(1+real((-alpha*(1./(1-bGM*log(c/aGM+1)).*p)).^nGM))).^ng
    %Phi = (1./(1+real((-alpha*p).^ng))).^nGM;
    %-alpha p 
    
    theta_neg = theta_r + (theta_s - theta_r) .* Phi;
    theta_pos = theta_s;
    
    neg = p < 0;
    theta = theta_neg.*neg + theta_pos.*~neg;
end