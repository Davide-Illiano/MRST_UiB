%% Surfactant PLOTS
close all
x = [50 100 150 200];

% Monolithic
n_N_1 = [100 100 100 100];

n_L_1 = [2100 2100 2100 2100];

n_P_1 = [465 470 480 490] ;

% Sequential classical
n_N_2 = [350 350 350 350];

n_L_2 = [3700 3800 3900 4100];

n_P_2 = [400 450 500 500 ] ;

% Sequential alternative
n_N_3 = [200 200 200 200];

n_L_3 = [2400 2700 3100 3200];

n_P_3 = [200 200 200 200] ;

figure
grid on, 
p = plot(x,n_L_1,'--o', x,n_N_1,':s', x,n_P_1,'-.*',...
    x,n_L_2,'--o',x,n_N_2,':s', x,n_P_2,'-.*',...
    x,n_L_3,'--o',x,n_N_3,':s', x,n_P_3,'-.*'),...
   xlabel('1/dx'), ylabel('Number of iterations') ,...
title('Total number of iterations'),...
    legend({'Mon-LS', 'Mon-Newton', 'Mon-Picard',...
    'NonLinS L-Scheme', 'NonLinS Newton', 'NonLinS Picard',...
    'AltS-LS', 'AltS-Newton', 'AltS-Picard'});
set(gca,...
'Units','normalized',...
'FontUnits','points',...
'FontWeight','normal',...
'FontSize',17,...
'FontName','Times')
set(gca, 'YScale', 'log')
%print -depsc2 cmwr2018.eps
for i = 1:9
p(i).LineWidth = 1.5;
end


%% Variably saturared domain
x = [10 20 40]
y = [10 20]
n_L_Mono = [10000 26000];
n_L_clasic = [4000 16000 44000];
n_L_alt = [3000 6000 12000];

figure('Units','inches',...
'PaperPositionMode','auto');
grid on, 
p = plot(y,n_L_Mono,'--o', x,n_L_clasic,'--o', x,n_L_alt,'--o'),...
   xlabel('1/dx'), ylabel('Number of iterations') ,...
title('Total number of iterations'),...
    legend({'Mon-LS', 'NonLinS L-Scheme','AltS-LS'});
set(gca,...
'Units','normalized',...
'FontUnits','points',...
'FontWeight','normal',...
'FontSize',17,...
'FontName','Times')
%set(gca, 'YScale', 'log')
%print -depsc2 cmwr2018.eps
for i = 1:9
p(i).LineWidth = 1.5;
end
