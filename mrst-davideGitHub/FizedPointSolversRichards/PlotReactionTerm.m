%% Plot of different schemes
% here anayticl case, data obtained from List Radu + REACTION TERM! r(c) =  c.^2./(c+1)
% Computational times and number of iterations
close all
clear all

%% Mesh size
x = [10 20 40 80];
%% Compare the diffrent schemes  
% example with inital and boundary pressure -2

%Comapre different L (max of derivatives was Theta_p = .15 Theta_c = .008)
% method 1 (monolithic)
% Lp_1 = .2    Lc_1 = .01
t_L_1_1 = [3.442 4.480 8.705 38.210];
n_L_1_1 = [37 39 42 44];

% Lp_2 = .1    Lc_2 = .005   %not working
%t_L_2_1 = [3.2 4.1 7.6 31.5];
%n_L_2_1 = [44 48 51 54];

% Lp_3 = .5    Lc_3 = .1
t_L_3_1 = [9.732 12.706 26.116 115.218];
n_L_3_1 = [107 114 121 129];

% method 2  (calssical approach sequential)
% Lp_1 = .2    Lc_1 = .01
t_L_1_2 = [11.532 15.101 23.815 59.899];
n_L_1_2 = [131 147 153 210];

% Lp_2 = .1    Lc_2 = .005
t_L_2_2 = [7.432 9.838 13.715 36.141];
n_L_2_2 = [80 92 93 124];

% Lp_3 = .5    Lc_3 = .1
t_L_3_2 = [23.283 31.099 47.377 131.063];
n_L_3_2 = [282 313 330 437];

% method 3   (alternative sequental)
% Lp_1 = .2    Lc_1 = .01
t_L_1_3 = [32.054 35.852 46.878 93.741];
n_L_1_3 = [135 143 151 159];

% Lp_2 = .1    Lc_2 = .005
t_L_2_3 = [19.966 22.051 28.749 55.638];
n_L_2_3 = [57 83 89 93];

% Lp_3 = .5    Lc_3 = .1
t_L_3_3 = [66.126 76.938 101.39 204.855];
n_L_3_3 = [295 319 339 359];
%%
%{
figure, grid on, plot(x,t_L_1_1,'--o', x,t_L_3_1,'--o'), xlabel('1/h'), ylabel('Time (s)') , title('Computational Times different L (FIM)'), legend('Lp_1 = .2 Lc_1 = .01', 'Lp_1 = .5 Lc_1 = .1','Location','NorthWest');
figure, grid on, plot(x,n_L_1_1,'--o',x,n_L_3_1,'--o'), title('Number of Iterations different L (FIM)'), legend('Lp_1 = .2 Lc_1 = .01', 'Lp_1 = .5 Lc_1 = .1','Location','NorthWest'), xlabel('1/h'), ylabel('Iterations') ; 

figure, grid on, plot(x,n_L_1_1,'--o',x,n_L_3_1,'--o'), title('Number of Iterations different L (FIM)'), legend('Lp_1 = .2 Lc_1 = .01', 'Lp_1 = .5 Lc_1 = .1','Location','NorthWest'), xlabel('1/h'), ylabel('Iterations') ; 

figure, grid on, plot(x,t_L_1_2,'--o',x,t_L_2_2,'--o',x,t_L_3_2,'--o'), xlabel('1/h'), ylabel('Time (s)') , title('Computational Times different L (Sequential)'), legend('Lp_1 = .2 Lc_1 = .01','Lp_1 = .1 Lc_1 = .005', 'Lp_1 = .5 Lc_1 = .1','Location','NorthWest');
figure, grid on, plot(x,n_L_1_2,'--o',x,n_L_2_2,'--o',x,n_L_3_2,'--o'), title('Number of Iterations different L (Sequential)'), legend('Lp_1 = .2 Lc_1 = .01','Lp_1 = .1 Lc_1 = .005', 'Lp_1 = .5 Lc_1 = .1','Location','NorthWest'), xlabel('1/h'), ylabel('Iterations') ; 

figure, grid on, plot(x,t_L_1_3,'--o',x,t_L_2_3,'--o',x,t_L_3_3,'--o'), xlabel('1/h'), ylabel('Time (s)') , title('Computational Times different L (Alternative Sequential)'), legend('Lp_1 = .2 Lc_1 = .01','Lp_1 = .1 Lc_1 = .005', 'Lp_1 = .5 Lc_1 = .1','Location','NorthWest');
figure, grid on, plot(x,n_L_1_3,'--o',x,n_L_2_3,'--o',x,n_L_3_3,'--o'), title('Number of Iterations different L (Alternative Sequential)'), legend('Lp_1 = .2 Lc_1 = .01','Lp_1 = .1 Lc_1 = .005', 'Lp_1 = .5 Lc_1 = .1','Location','NorthWest'), xlabel('1/h'), ylabel('Iterations') ; 
%}

% Newton
t_N_1 = [1.163 1.352 1.993 6.677];
n_N_1 = [4 4 4 4];

t_N_2 = [1.653 1.924 2.734 6.218];
n_N_2 = [8 12 13 13];
n_N_ext_2 = [2 3 3 3];

t_N_3 = [3.162 3.5 4.404 8.429];
n_N_3 = [9 11 11 11];
n_N_ext_3 = [5 6 6 6];

% Picard
t_P_1 = [.999 1.092 1.517 4.105];
n_P_1 = [4 4 4 4];

t_P_2 = [1.585 1.765 2.522 5.853];
n_P_2 = [8 12 13 17];
n_P_ext_2 = [2 3 3 4];

t_P_3 = [3.466 3.796 4.950 10.304];
n_P_3 = [9 11 11 11];
n_P_ext_3 = [5 6 6 6];
%{
figure, grid on, plot(x,t_L_1_1,'--o',x,t_N_1,'--*',x,t_P_1,'--o'), title('Computational Times different Schemes (FIM)'), legend('L Scheme', 'Newton Method', 'Modified Picard','Location','NorthWest'), xlabel('1/h'), ylabel('Time (s)');
figure, grid on, plot(x,n_L_1_1,'--o',x,n_N_1,'--*',x,n_P_1,'--o'), title('Number of Iterations different Schemes (FIM)'), legend('L Scheme', 'Newton Method', 'Modified Picard','Location','NorthWest'), xlabel('1/h'), ylabel('Iterations') ;

figure, grid on, plot(x,t_L_2_2,'--o',x,t_N_2,'--*',x,t_P_2,'--o'), title('Computational Times different Schemes (Sequential)'), legend('L Scheme', 'Newton Method', 'Modified Picard','Location','NorthWest'), xlabel('1/h'), ylabel('Time (s)');
figure, plot(x,n_L_2_2,'--o',x,n_N_2,'--*',x,n_P_2,'--o'), title('Number of Iterations different Schemes (Sequential)'), legend('L Scheme', 'Newton Method', 'Modified Picard','Location','NorthWest'), xlabel('1/h'), ylabel('Iterations') ;

figure, grid on, plot(x,t_L_1_3,'--o',x,t_N_3,'--*',x,t_P_3,'--o'), title('Computational Times different Schemes (Alternative Sequential)'), legend('L Scheme', 'Newton Method', 'Modified Picard','Location','NorthWest'), xlabel('1/h'), ylabel('Time (s)');
figure, grid on, plot(x,n_L_1_3,'--o',x,n_N_3,'--*',x,n_P_3,'--o'), title('Number of Iterations different Schemes (Alternative Sequential)'), legend('L Scheme', 'Newton Method', 'Modified Picard','Location','NorthWest'), xlabel('1/h'), ylabel('Iterations') ;
%}
%%
figure, grid on, plot(x,t_L_1_1,'--o', x,t_L_3_1,'--o', x,t_L_1_2,'--o', x,t_L_2_2,'--o',x,t_L_3_2,'--o', x,t_L_1_3,'--o', x,t_L_2_3,'--o',x,t_L_3_3,'--o'), xlabel('1/h'), ylabel('Time (s)') , title('Computational Times different L Schemes'), legend('Lp_1 = .2 Lc_1 = .01 (FIM)', 'Lp_3 = .5 Lc_3 = .1 (FIM)',...
    'Lp_1 = .2 Lc_1 = .01 (Sequential)', 'Lp_2 = .1 Lc_2 = .005 (Sequential)', 'Lp_3 = .5 Lc_3 = .1 (Sequential)', ...
    'Lp_1 = .2 Lc_1 = .01 (Alternative Seq)', 'Lp_2 = .1 Lc_2 = .005 (Alternative Seq)', 'Lp_3 = .5 Lc_3 = .1 (Alternative Seq)', ...
    'Location','NorthWest');

figure, grid on, plot(x,n_L_1_1,'--o', x,n_L_3_1,'--o', x,n_L_1_2,'--o', x,n_L_2_2,'--o',x,n_L_3_2,'--o', x,n_L_1_3,'--o', x,n_L_2_3,'--o',x,n_L_3_3,'--o'), xlabel('1/h'), ylabel('Number Iterations') , title('Number Iterations different L Schemes'), legend('Lp_1 = .2 Lc_1 = .01 (FIM)', 'Lp_3 = .5 Lc_3 = .1 (FIM)',...
    'Lp_1 = .2 Lc_1 = .01 (Sequential)', 'Lp_2 = .2 Lc_2 = .005 (Sequential)', 'Lp_3 = .5 Lc_3 = .1 (Sequential)', ...
    'Lp_1 = .2 Lc_1 = .01 (Alternative Seq)', 'Lp_2 = .1 Lc_2 = .005 (Alternative Seq)', 'Lp_3 = .5 Lc_3 = .1 (Alternative Seq)', ...
    'Location','NorthWest');

figure, grid on, plot(x,t_L_1_1,'--o', x,t_N_1,'--o', x,t_P_1,'--o',...
    x,t_L_2_2,'--o',x,t_N_2,'--o', x,t_P_2,'--o',...
    x,t_L_2_3,'--o',x,t_N_3,'--o', x,t_P_3,'--o'), xlabel('1/h'), ylabel('Time (s)') , title('Computational Times different Schemes and Approaches'),...
    legend('Lp_1 = .2 Lc_1 = .01 (FIM)', 'Newton (FIM)', 'Picard (FIM)',...
    'Lp_2 = .1 Lc_2 = .005 (Sequential)', 'Newton (Sequential)', 'Picard (Sequential)',...
    'Lp_2 = .1 Lc_2 = .005 (Alternative Seq)', 'Newton (Alternative Seq)', 'Picard (Alternative Seq)',...
    'Location','NorthWest');



figure('Units','inches',...
'PaperPositionMode','auto');
grid on, plot(x,n_L_1_1,'--o', x,n_N_1,'-.o', x,n_P_1,'--*',...
    x,n_L_2_2,'-.o',x,n_N_2,'-.o', x,n_P_2,'-.*',...
    x,n_L_2_3,'--o',x,n_N_3,'-.o', x,n_P_3,'--*'), xlabel('1/h', 'FontUnits','points',...
'interpreter','latex',...
'FontSize',11,...
'FontName','Times'), ylabel('Number Iterations', 'FontUnits','points',...
'interpreter','latex',...
'FontSize',11,...
'FontName','Times') ,...
    legend('Lp_1 = .2 Lc_1 = .01 (FIM)', 'Newton (FIM)', 'Picard (FIM)',...
    'Lp_2 = .1 Lc_2 = .005 (Sequential)', 'Newton (Sequential)', 'Picard (Sequential)',...
    'Lp_2 = .1 Lc_2 = .005 (Alternative Seq)', 'Newton (Alternative Seq)', 'Picard (Alternative Seq)',...
    'Location','NorthWest');
set(gca,...
'Units','normalized',...
'FontUnits','points',...
'FontWeight','normal',...
'FontSize',12,...
'FontName','Times')
 print -depsc2 ﻿NumberIterationsDifferentschemesReaction2.eps

%%

figure, grid on, plot( x,t_N_1,'--o', x,t_P_1,'--o',...
    x,t_N_2,'--o', x,t_P_2,'--o',...
   x,t_N_3,'-.o', x,t_P_3,'--o'), xlabel('1/h'), ylabel('Time (s)') , title('Computational Times different Schemes and Approaches'),...
    legend( 'Newton (FIM)', 'Picard (FIM)',...
     'Newton (Sequential)', 'Picard (Sequential)',...
    'Newton (Alternative Seq)', 'Picard (Alternative Seq)',...
    'Location','NorthWest');

figure('Units','inches',...
'PaperPositionMode','auto'); grid on, plot( ...
    x,n_N_2,'-.*', x,n_P_2,'-.*',...
    x,n_N_3,'-.+', x,n_P_3,'--o'), xlabel('1/h', 'FontUnits','points',...
'interpreter','latex',...
'FontSize',11,...
'FontName','Times') , ylabel('Number of iterations', 'FontUnits','points',...
'interpreter','latex',...
'FontSize',11,...
'FontName','Times')  ,
    legend( ...
     'Newton (Sequential)', 'Picard (Sequential)',...
     'Newton (Alternative Seq)', 'Picard (Alternative Seq)',...
    'Location','NorthWest'), axis([10 80 3 21]);
 set(gca,...
'Units','normalized',...
'FontUnits','points',...
'FontWeight','normal',...
'FontSize',12,...
'FontName','Times')



figure('Units','inches',...
'PaperPositionMode','auto');
grid on, plot( x,n_N_ext_2,'--*', x,n_N_ext_3,'-.+',...
    x,n_P_ext_2,'--o', x,n_P_ext_3,'-.o'), xlabel('1/h', 'FontUnits','points',...
'interpreter','latex',...
'FontSize',11,...
'FontName','Times') , ylabel('Number of iterations', 'FontUnits','points',...
'interpreter','latex',...
'FontSize',11,...
'FontName','Times') , ...
    legend( 'Newton (Sequential)','Newton (Alternative Seq)',...
     'Picard (Sequential)', 'Picard (Alternative Seq)',...
     'Location','SouthEast'), axis([10 80 2 6.5]);
 set(gca,...
'Units','normalized',...
'FontUnits','points',...
'FontWeight','normal',...
'FontSize',12,...
'FontName','Times')

 
  %% Try with extra lage times dx = [10 20 40 80]  dt = dx
 % newton sequential
 t_N_seq = [2.583 5.778 16.5 79.427]; % case dx=160 dt=80
 n_N_seq = [10 11 21 20];
 n_N_ex = [2 2 3 2];
 
 t_N_2seq = [3.455 8.356 29.653 158.723];
 n_N_2seq = [9 15 29 25];
 n_N_2ex = [5 8 15 13];
 
 t_P_seq = [2.146 4.49 14.558 65.452];
 n_P_seq = [8 11 21 20];
 n_P_ex = [2 2 3 2];
 
 t_P_2seq = [3.413 9.191 35.53 ]; % no convergence
 n_P_2seq = [9 15 29 ];
 n_P_2ex = [5 8 15 ]; 

 %L_p = .2  L_c= .01
 t_L1_seq = [19.309];
 n_L1_seq = [133];
 n_L1_ex = [15 ];
 
 t_L1_2seq = [36.05];
 n_L1_2seq = [135 ];
 n_L1_2ex = [68 ]; 
 
 %L_p = .1  L_c= .005
 t_L1_seq = [11.908]; 
 n_L1_seq = [80];
 n_L1_ex = [10 ];
 
 t_L1_2seq = [22.529];
 n_L1_2seq = [79 ];
 n_L1_2ex = [40 ]; 
  