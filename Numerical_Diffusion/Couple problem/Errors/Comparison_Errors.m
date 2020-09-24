%% Compute the errors between solutions
 clc
 clear all; close all;

I = 41; J=61; %I = 41; J=61; % 
a=0; b=2;
c=0; d=3;
dx = (b-a)/(I-1);
x = [a:dx:b];
x2 = (x(1:I-1)+x(2:I))/2;
dy=(d-c)/(J-1);
y=c:dy:d;
levels=12;

disp('Clay Errors Gauss')
load('Results_clay_41_61.mat')
load('sol_grw_Richy_2D_sat_clay41_61.mat')

disp('Error p')
norm(p-p_clay)/norm(p_clay)

disp('Error c')
norm(c-c_clay)/norm(c_clay)

disp('Error theta')
norm(tht-theta_clay)/norm(theta_clay)

disp('Error vx')
norm(Vx-vx_clay)/norm(vx_clay)

disp('Error vy')
norm(Vy-vy_clay)/norm(vy_clay)
% plots MRST results
plot_contours(I,J,x,y,p_clay,c_clay,theta_clay,vx_clay,vy_clay,levels);

% plots errors
% plot_contours(I,J,x,y,p_clay-p,c_clay-c,theta_clay-tht,vx_clay-Vx,vy_clay-Vy,levels);
%%
disp('Loam Errors Gauss')
load('Results_loam_41_61.mat')
load('sol_grw_Richy_2D_sat_loam41_61.mat')

disp('Error p')
norm(p-p_loam)/norm(p_loam)

disp('Error c')
norm(c-c_loam)/norm(c_loam)

disp('Error theta')
norm(tht-theta_loam)/norm(theta_loam)

disp('Error vx')
norm(Vx-vx_loam)/norm(vx_loam)

disp('Error vy')
norm(Vy-vy_loam)/norm(vy_loam)

% plots MRST results
plot_contours(I,J,x,y,p_loam,c_loam,theta_loam,vx_loam,vy_loam,levels);

% plots errors
% plot_contours(I,J,x,y,p_loam-p,c_loam-c,theta_loam-tht,vx_loam-Vx,vy_loam-Vy,levels);