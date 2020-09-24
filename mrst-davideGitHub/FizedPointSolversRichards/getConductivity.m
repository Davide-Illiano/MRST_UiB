function K = getConductivity(p, theta, K_s, n)   % i added K_s
    K_pos = K_s.*theta.^(1/2).*(1 - (1 - theta.^(n./(n-1))).^((n-1)./n)).^2;
    %K_pos = K_s;
    K_neg = K_s;
    
    
    neg = p <= 0;
    K = K_neg.*neg + K_pos.*~neg;
   
    %K = theta;
%{    
theta_s = .3;
theta_r = .05;
n = 1.82;
m = 1-1/n;
alpha = .0154;
K_s = 0.000166666667*.06;
l=.153;
aa = 0.044;
bb = 0.4745;
m_k=1;
Phi = .377;
xi = 1;
D_0 = 6*10^(-4);
   % K = K_pos;
    %%Salinity
    THETA = (theta - theta_r )/ (theta_s - theta_r);
    K = K_s*THETA.^l.*(1-(1-THETA.^(1/m).^m)).^2;    %theta.^(1/2).*(1 - (1 - theta.^(n/(n-1))).^((n-1)/n)).^2; 
    %K = K_s*(( ((theta_s - theta_r)./(1+abs(alpha*(1-bb*log(c./aa+1))^(-1).*Psi).^n).^m +theta_r) - theta_r )/ (theta_s - theta_r)).^l ...
    %.*(1-(1-(( ((theta_s - theta_r)./(1+abs(alpha*(1-bb*log(c./aa+1))^(-1).*Psi).^n).^m +theta_r) - theta_r )/ (theta_s - theta_r)).^(1/m_k)).^m_k).^2;
%}
end