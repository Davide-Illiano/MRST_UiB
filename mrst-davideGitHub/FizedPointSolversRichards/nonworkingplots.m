%%plot presentation
% complex case
%% Mesh size
close all

clear all
x = [10 20 40 80];
z = [10 20];
t = [10 20 40];
%% Only non working newton
%% Monolithic
% Lp_1 = .1    Lc_1 = .2
t_L_1 = [1.84 2.22 3.84 15];
n_L_1 = [110 120 132 149];

% Picard
t_P_1 = [2.34 3.05 5.73 22.7];
n_P_1 = [128 139 150 185];

% Newton
t_N_1 = [ 1.38 1.6];
n_N_1 = [ 40 40];
% Newton smaller timestep
t_N_1_small = [ 1.38 2.84 11.51];
n_N_1_small = [ 40 90 270];

%% Classical sequential
% Newton sequential
t_N_2 = [2.23 2.4];
n_N_2 = [74 84];

t_N_2_small = [2.23 4.1 17.2];
n_N_2_small = [74 164 672];

% Newton transport but L richards
t_mixed_2 = [7.6 9.16 14.11 37.28];
n_mixed_2 = [593 625 692 803];

% L sequential
% Lp_1 = .1    Lc_1 = .2
t_L_2 = [7.85 9.1 13.8 35];
n_L_2 = [538 604 673 769];

% Picard
t_P_2 = [2.21 2.65 ];
n_P_2 = [82 84];

t_P_2_small = [2.21 4.18 15.64];
n_P_2_small = [82 176 727];

%% Alternative sequential
% Lscheme
t_L_3 = [10.4 15.4 19.9 44.8];
n_L_3 = [62 84 96 108];

% Picard did not work

% Newton did not work

% Lscheme on richards newton on transport
t_mixed_3 = [11.3 15.1 20.1 46.8];
n_mixed_3 = [74 86 100 112];



%% plot also time just curiosity
%{
plot(x,t_L_1,'--o', x,t_L_2 ,'--+' ,...
    x,t_L_3,'-.', x, t_P_1,'--o', x, t_mixed_2,'--+', x, t_mixed_3, '-.',...
    z, t_N_1, '--o', z, t_P_2,'--+', z, t_N_2,'--+', ...
    t, t_N_1_small, t, t_N_2_small, t, t_P_2_small );


xlabel('1/h', 'FontUnits','points',...
'interpreter','latex',...
'FontSize',11,...
'FontName','Times'), ylabel('CPU Time', 'FontUnits','points',...
'interpreter','latex',...
'FontSize',11,...
'FontName','Times') ,...
legend({'Lp_1 = .1 Lc_1 = .2 (FIM)','Lp_1 = .1 Lc_1 = .2 (sequential)',...
'Lp_1 = .1 Lc_1 = .2 (alternative)', 'Picard (FIM)', ...
'Mixed approach (Seq)','Mixed approach (Altern)','Newton (FIM)',...
'Picard(sequential)','Newton(sequential)', 'Newton dt=dx(FMI)',...
'Newton dt=dx(sequential)','Newton dt=dx(alternative)'},...
      'FontUnits','points',...
'interpreter','latex',...
'FontSize',11,...
'FontName','Times',...
    'Location','NorthWest');
set(gca,...
'Units','normalized',...
'FontUnits','points',...
'FontWeight','normal',...
'FontSize',12,...
'FontName','Times')
%}

p = plot(x,n_L_1,'--o', x,n_L_2 ,'--o' ,...
    x,n_L_3,'--o', x, n_P_1,'-.*', ...
    z, n_N_1, ':s', z, n_N_2,':s', z, n_P_2,'-.*' );

%t, n_N_1_small, t, n_N_2_small,t, n_P_2_small
% 'Newton dt=dx(FMI)',...'Newton dt=dx(sequential)','Newton dt=dx(alternative)


xlabel('1/dx'), ylabel('Total Number Iterations') ,...
title('Total number of iterations'),...
legend({'L Scheme (mono)','L Scheme (seq)',...
'L Scheme (alt seq)','Picard (mono)',...
'Newton (mono)','Newton (seq)','Picard (seq)'});
set(gca,...
'Units','normalized',...
'FontUnits','points',...
'FontWeight','normal',...
'FontSize',17,...
'FontName','Times')
grid on
for i = 1:7
p(i).LineWidth = 1.5;
end

figure
plot(x,n_L_1,'--o',...
    x,n_L_3,'-.', x, n_P_1,'--o',  x, n_mixed_3, '--*',...
    z, n_N_1, '--o', z, n_P_2,'--+', z, n_N_2,'--+');
xlabel('1/h', 'FontUnits','points',...
'interpreter','latex',...
'FontSize',14,...
'FontName','Times'), ylabel('Number Iterations', 'FontUnits','points',...
'interpreter','latex',...
'FontSize',14,...
'FontName','Times') ,...
legend({'Lp_1 = .1 Lc_1 = .2 (FIM)',...
'Lp_1 = .1 Lc_1 = .2 (alternative)', 'Picard (FIM)', ...
'Mixed approach (Altern)','Newton (FIM)',...
'Picard(sequential)','Newton(sequential)'},...
      'FontUnits','points',...
'interpreter','latex',...
'FontSize',16,...
'FontName','Times',...
    'Location','NorthWest');
set(gca,...
'Units','normalized',...
'FontUnits','points',...
'FontWeight','normal',...
'FontSize',15,...
'FontName','Times')



%% Compare the diffrent schemes  
% example with inital and boundary pressure -2

%Comapre different L (max of derivatives was Theta_p = .15 Theta_c = .008)
% method 1 (monolithic)
% Lp_1 = .2    Lc_1 = .01
t_L_1_1 = [5.77 7.43 13.53 58.11];
n_L_1_1 = [59 66 74 85];

% Lp_2 = .1    Lc_2 = .005  
t_L_2_1 = [2.91 3.77 5.71 24.62];
n_L_2_1 = [27 31 34 39];



% method 2  (calssical approach sequential)
% Lp_1 = .2    Lc_1 = .01
t_L_1_2 = [12.24 15.61 24.15 62.34];
n_L_1_2 = [121 136 154 177];

% Lp_2 = .1    Lc_2 = .005
t_L_2_2 = [6.98 9.75 13.13 39.61];
n_L_2_2 = [68 76 86 98];



% method 3   (alternative sequental)
% Lp_1 = .2    Lc_1 = .01
t_L_1_3 = [16.9 21.23 31.10 73.47];
n_L_1_3 = [122 142 162 180];

% Lp_2 = .1    Lc_2 = .005
t_L_2_3 = [12.73 13.07 19.26 40.69];
n_L_2_3 = [66 78 88 98];



% Newton
t_N_1 = [ 1.48 1.62];
n_N_1 = [ 5 5];

t_N_2 = [2.48 2.77];
n_N_2 = [11 11];
n_N_ext_2 = [0 0 0 0 ];

t_N_3 = [0 0 0 0];
n_N_3 = [0 0 0 0];
n_N_ext_3 = [0 0 0 0 ];

% Picard
t_P_1 = [0 0 0 0];
n_P_1 = [0 0 0 0];

t_P_2 = [2.48 2.55 ];
n_P_2 = [11 11];
n_P_ext_2 = [0 0 ];

t_P_3 = [0 0 0 0];
n_P_3 = [0 0 0 0];
n_P_ext_3 = [0 0 0 0];

% mixed method L_1 = .1 L_2=.005 

figure('Units','inches',...
'PaperPositionMode','auto');
hold on
grid on
plot(x,n_L_1_1,'--o',x,n_L_2_1,'--+');

axis([0 90 0 250])

plot(x,n_L_1_2,'--o',x,n_L_2_2 ,'--+' ,...
    x,n_L_1_3,'--o', x, n_L_2_3,'--o', z, n_N_1, '+',...
    z, n_P_2, 'o', z, n_N_2, '+' );

plot(y, n_N_2 )



xlabel('1/h', 'FontUnits','points',...
'interpreter','latex',...
'FontSize',11,...
'FontName','Times'), ylabel('Number Iterations', 'FontUnits','points',...
'interpreter','latex',...
'FontSize',11,...
'FontName','Times') ,...
legend({'Lp_1 = .2 Lc_1 = .01 (FIM)','Lp_1 = .1 Lc_1 = .005 (FIM)','Lp_1 = .2 Lc_1 = .01 (sequential)',...
'Lp_1 = .1 Lc_1 = .005 (sequential)','Lp_1 = .2 Lc_1 = .01 (alternative)',...
'Lp_1 = .1 Lc_1 = .005 (alternative)','Newton (FIM)',...
'Picard(sequential)','Newton(alternative)'},...
      'FontUnits','points',...
'interpreter','latex',...
'FontSize',11,...
'FontName','Times',...
    'Location','NorthWest');
set(gca,...
'Units','normalized',...
'FontUnits','points',...
'FontWeight','normal',...
'FontSize',12,...
'FontName','Times')
 print -depsc2 iterationsdifferentschemesnonworkingcase.jpg

hold off
%%
%dx=dt
%Comapre different L (max of derivatives was Theta_p = .15 Theta_c = .008)
% method 1 (monolithic)
% Lp_1 = .2    Lc_1 = .01
t_L_1_1 = [];
n_L_1_1 = [];

% Lp_2 = .1    Lc_2 = .005  
t_L_2_1 = [2.87 ];
n_L_2_1 = [31 ];



% method 2  (calssical approach sequential)
% Lp_1 = .2    Lc_1 = .01
t_L_1_2 = [];
n_L_1_2 = [ ];

% Lp_2 = .1    Lc_2 = .005
t_L_2_2 = [  ];
n_L_2_2 = [ ];



% method 3   (alternative sequental)
% Lp_1 = .2    Lc_1 = .01
t_L_1_3 = [];
n_L_1_3 = [ ];

% Lp_2 = .1    Lc_2 = .005
t_L_2_3 = [];
n_L_2_3 = [ ];



% Newton
t_N_1 = [ ];
n_N_1 = [ ];

t_N_2 = [];
n_N_2 = [];
n_N_ext_2 = [ ];

t_N_3 = [0 0 0 0];
n_N_3 = [0 0 0 0];
n_N_ext_3 = [0 0 0 0 ];

% Picard
t_P_1 = [0 0 0 0];
n_P_1 = [0 0 0 0];

t_P_2 = [0 0 0 0];
n_P_2 = [0 0 0 0];
n_P_ext_2 = [0 0 0 0];

t_P_3 = [0 0 0 0];
n_P_3 = [0 0 0 0];
n_P_ext_3 = [0 0 0 0];

