%% Version from SourceTree
%Construct 3D grid with 50 cells in the x-direction
mrstModule add ad-core ad-props
gravity reset on
mrstVerbose off

close all

G = cartGrid([20, 1, 20], [100, 100, 10]*meter);
G = computeGeometry(G);

% Homogenous rock properties
rock = struct('perm', 1000*darcy*ones(G.cells.num, 1), ...
              'poro', 0.5*ones(G.cells.num, 1));

% Default oil-water fluid with unit values
fluid = initSimpleADIFluid('phases', 'W');
fluid.theta = @(p) getTheta(p, 0.026, 0.42, 0.95, 2.9);
fluid.Kmult = @(p, theta) getConductivity(theta, 2.9);   %.8 is K_s probably unrealistic number, required in def of K
fluid.theta_first = @(p) getThetafirst(p, 0.026, 0.42, 0.95, 2.9);

% Set up model and initial state.
model = RichardsEquationModel(G, rock, fluid);
state0 = initResSol(G, -1);

% Set up drive mechanism: constant rate at x=0, constant pressure at x=L
time = 365*day;
n = 50;
dT = time/n;

% pv = poreVolume(G, rock);
% injRate = -sum(pv)/(time);

bc = [];
% bc = fluxside(bc, G, 'xmin', -injRate, 'sat', [1]);
bc = pside(bc, G, 'zmax', 10, 'sat', [1]);

schedule = simpleSchedule(repmat(dT,1,n), 'bc', bc);



model.nonlinearTolerance = 1e-4;
%% Newton method
[~,states, report] = simulateScheduleAD(state0, model, schedule);

%% Picard method
modelFP = RichardsEquationModelPicardScheme(G, rock, fluid);
modelFP.nonlinearTolerance = model.nonlinearTolerance;

nls = NonLinearSolver('maxIterations', 10000);
[~,statesFP, reportFP] = simulateScheduleAD(state0, modelFP, schedule, 'nonlinearsolver', nls);

%% L scheme
modelFPL = RichardsEquationModelLScheme(G, rock, fluid, 'L', .5);
modelFPL.nonlinearTolerance = model.nonlinearTolerance;

nls = NonLinearSolver('maxIterations', 10000);
[~,statesFPL, reportFPL] = simulateScheduleAD(state0, modelFPL, schedule, 'nonlinearsolver', nls);


%% Plot the result Newton
mrstModule add mrst-gui
close all
plotToolbar(G, states)
title('Newton')

%% Plot results modified Picard
figure;
plotToolbar(G, statesFP)
title('Modified Picard')

%% Plot results L scheme
figure;
plotToolbar(G, statesFP)
title('L-scheme')

%%
figure;
plot([report.Iterations, reportFP.Iterations, reportFPL.Iterations])
legend('Newton','ModifiedPicard', 'L-scheme')
set(gca, 'YScale', 'log')
