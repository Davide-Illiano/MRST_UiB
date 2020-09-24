% Construct 3D grid with 50 cells in the x-direction
mrstModule add ad-blackoil ad-core ad-props mrst-gui blackoil-sequential
clear all
close all
clc

gravity reset on

G = cartGrid([40, 40], [100, 100]*meter); %cartGrid([40, 1, 40], [100, 100, 10]*meter);
G = computeGeometry(G);

% Homogenous rock properties
rock = struct('perm', 1000*darcy*ones(G.cells.num, 1), ...
              'poro', 0.5*ones(G.cells.num, 1));

% Default oil-water fluid with unit values
fluid = initSimpleADIFluid('phases', 'W');
fluid.theta = @(p, c) getThetaCoupled(p, c, 0.05, 0.3, 0.0551, .6, .07, 2.43);  % (p, c, theta_R, theta_S, alpha, zeta, sigma, n)
fluid.Kmult = @(p, theta) getConductivity(p, theta, 2.9);   %.8 is K_s probably unrealistic number, required in def of K
%fluid.theta_first = @(p) getThetafirst(p, 0.026, 0.42, 0.95, 2.9);

% Set up model and initial state.
model = TransportEquationModel(G, rock, fluid);
state0 = initResSol(G, -1); %-1
state0.c = ones(G.cells.num,1);

% Set up drive mechanism: constant rate at x=0, constant pressure at x=L
time = 365*day;
n = 80;
dT = time/n;

 %pv = poreVolume(G, rock);
 %injRate = -sum(pv)/(time);

bc = [];
%bc = fluxside(bc, G, 'xmin', -injRate, 'sat', [1]);
bc = pside(bc, G, 'ymax', 20, 'sat', [1]);%pside(bc, G, 'zmax', 20, 'sat', [1]); %10
bc.c = 4.*ones(size(bc.sat,1), 1);

schedule = simpleSchedule(repmat(dT,1,n), 'bc', bc);


model.nonlinearTolerance = 1e-4;
%%
[~,states, report] = simulateScheduleAD(state0, model, schedule);

%% L scheme
modelFPLScheme = RichardsTransportEquationModelLScheme(G, rock, fluid, 'L_p', .07, 'L_c', .8);  
modelFPLScheme.nonlinearTolerance = model.nonlinearTolerance;
nls = NonLinearSolver('maxIterations', 10000);
[~,statesFPLScheme, reportFPLScheme] = simulateScheduleAD(state0, modelFPLScheme, schedule, 'nonlinearsolver', nls);

%% Picard scheme
modelFPPicard = RichardsTransportEquationModelPicardScheme(G, rock, fluid);
modelFPPicard.nonlinearTolerance = model.nonlinearTolerance;
nls = NonLinearSolver('maxIterations', 10000);
[~,statesFPPicard, reportFPPicard] = simulateScheduleAD(state0, modelFPPicard, schedule, 'nonlinearsolver', nls);


%% Plot the result
%% Plot the result Newton
mrstModule add mrst-gui
close all
plotToolbar(G, states)
title('Newton')
colorbar

h = colorbar;
set(h, 'ylim', [1 2]) 
%% Plot results modified Picard
figure;
plotToolbar(G, statesFPPicard)
title('Modified Picard')

%% Plot results L scheme
figure;
plotToolbar(G, statesFPLScheme)
title('L-scheme')

%%
figure;
plot([report.Iterations, reportFPPicard.Iterations, reportFPLScheme.Iterations])
legend('Newton','ModifiedPicard', 'L-scheme')
set(gca, 'YScale', 'log')