%% Try to compute external forces using manufactured solutions p_a = t*x*(1-x)*y*(1-y) and c_a = p_a
% we need also analytical expressions for theta(p,c) and K(theta)
syms x y t
p_a = t*y;
c_a = t*y;

theta = 1/2*p_a*c_a^2;
K = theta.^3;

f_1 = diff(theta,t) - divergence( K.* gradient(p_a, [x y]) + 1,[x y])
f_2 = diff(theta.*c_a,t) + divergence( - theta.*gradient(c_a, [x y]) + (K.* gradient(p_a, [x y]) + 1).*c_a ,[x y])

