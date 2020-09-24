clear all
close all
x = [10 20 40 80];

%% Number iterations
L_mono = [41 40 40 45 ];
L_seq = [20 20 20 31];
L_alt = [10 10 10 15];

figure
plot(x, L_mono, '*-', x, L_seq, '*-', x, L_alt, '*-')

axis([10 80 5 60])
title('Numbers of iterations')
legend('L Scheme (mono) L_1=.008', 'L Scheme (seq) L_1=.0001', 'L Scheme (alt) L_1=.0001')


%% Residual plot
o = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17];
x = [1 2 3 4];
y= [1 2 3];
z = [1 2]; 

L_mono_P = [3.02e-03 1.79e-03 1.06e-03 6.30e-04 3.74e-04 2.22e-04 1.31e-04 7.80e-05 4.63e-05 2.74e-05 1.63e-05 9.66e-06 ...
    5.73e-06 3.40e-06 2.02e-06 1.20e-06 7.10e-07];

L_seq_P = [3.02e-03 1.35e-05 1.45e-08];

L_alt_P = [3.02e-03 1.35e-05 1.45e-08];

figure
grid on
plot(o, log(L_mono_P), '*-', y, log(L_seq_P), '*-', y, log(L_alt_P), 'o-')
title('Pressure residuals')
legend('L Scheme (mono) L_1=.008', 'L Scheme (seq) L_1=.0001', 'L Scheme (alt) L_1=.0001')
xlabel('Number of iterations')
ylabel('Log(residuals)')



pause
clear all
close all
x = [10 20 40 80];

%% Also hysteresis
%% Number iterations
N_mono = [30 30 30 30];
N_seq = [30 30 30 30];
N_alt = [20 20 20 20];

L_mono = [50 40 50 70];
L_seq = [70 60 60 60 ];
L_alt = [50 40 40 40];


plot(x, N_mono, 'o-', x, N_seq, 'o-', x, N_alt, 'o-',...
    x, L_mono, '*-', x, L_seq, '*-', x, L_alt, '*-')

axis([10 80 10 80])
title('Numbers of iterations')
legend('Newton (mono)', 'Newton (seq)', 'Newton (alt)', ...
    'L Scheme (mono)', 'L Scheme (seq)', 'L Scheme (alt)')


%% Residual plot
o = [1 2 3 4 5];
x = [1 2 3 4];
y= [1 2 3];
z = [1 2]; 

N_mono_P = [2.15e-02 1.43e-04 1.03e-09];
N_mono_c = [1.79e-02 4.79e-05 8.58e-10];
N_mono_t = [1.68e-02 2.41e-06];


N_seq_P = [2.04e-02 7.64e-06];
N_seq_c = [1.69e-02 2.45e-06];
N_seq_t = [1.96e-02 1.94e-06];


N_alt_P = [2.04e-02 7.64e-06];
N_alt_c = [1.96e-02 3.72e-07];
N_alt_t = [1.69e-02 2.45e-06];


L_mono_P = [2.15e-02 9.59e-05 1.26e-05 1.18e-06];
L_mono_c = [1.63e-02 1.15e-04 5.53e-07];
L_mono_t = [1.68e-02 1.69e-03 1.65e-04 1.63e-05 1.62e-06];


L_seq_P = [2.04e-02 4.88e-05 1.12e-06];
L_seq_c = [2.24e-02 5.49e-06];
L_seq_t = [1.69e-02 1.71e-03 1.70e-04 1.69e-05  1.68e-06];


L_alt_P = [2.04e-02 4.88e-05 1.12e-06];
L_alt_c = [2.24e-02 1.35e-05  1.27e-06];
L_alt_t = [1.69e-02 1.71e-03 1.70e-04 1.69e-05 1.68e-06];

figure
grid on
plot(y, log(N_mono_P), 'o-', z,log( N_seq_P), 'o-', z,log( N_alt_P), 'o-',...
    x, log(L_mono_P), '*-', y, log(L_seq_P), '*-', y, log(L_alt_P), '*-')
title('Pressure residuals')
legend('Newton (mono)', 'Newton (seq)', 'Newton (alt)', ...
    'L Scheme (mono)', 'L Scheme (seq)', 'L Scheme (alt)')
xlabel('Number of iterations')
ylabel('Log(residuals)')


figure
grid on
plot(y, log(N_mono_c), 'o-', z,log( N_seq_c), 'o-', z,log( N_alt_c), 'o-',...
    y, log(L_mono_c), '*-', z, log(L_seq_c), '*-', y, log(L_alt_c), '*-')
title('Concentration residuals')
legend('Newton (mono)', 'Newton (seq)', 'Newton (alt)', ...
    'L Scheme (mono)', 'L Scheme (seq)', 'L Scheme (alt)')
xlabel('Number of iterations')
ylabel('Log(residuals)')


figure
grid on
plot(z, log(N_mono_t), 'o-', z,log( N_seq_t), 'o-', z,log( N_alt_t), 'o-',...
    o, log(L_mono_t), '*-', o, log(L_seq_t), '*-', o, log(L_alt_t), '*-')
title('Theta residuals')
legend('Newton (mono)', 'Newton (seq)', 'Newton (alt)', ...
    'L Scheme (mono)', 'L Scheme (seq)', 'L Scheme (alt)')
xlabel('Number of iterations')
ylabel('Log(residuals)')
