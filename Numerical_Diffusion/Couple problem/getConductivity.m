function K = getConductivity(p,theta, theta_r, theta_s, Ks, ng)   % i added K_s

    % Ks.*((tht - theta_res)./(theta_sat-theta_res)).^0.5.*(1-(1-((tht - theta_res)./(theta_sat-theta_res)).^(1/ng)).^ng).^2; % 8,7 min. 
    ng = (ng-1)/ng;
    K_neg = Ks.*((theta - theta_r)./(theta_s-theta_r)).^0.5.*(1-(1-((theta - theta_r)./(theta_s-theta_r)).^(1/ng)).^ng).^2;
    K_pos = Ks;
    
    
    neg = p < 0;
    K = K_neg.*neg + K_pos.*~neg;
    %}

end