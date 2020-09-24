%% Version from SourceTree
%Construct 3D grid with 50 cells in the x-direction
mrstModule add ad-core ad-props
gravity reset on

clea all
close all

mrstVerbose off
o=1;

for X = [ 50 ]
G = cartGrid([X, 1, X], [100, 100, 10]*meter);
G = computeGeometry(G);

% Homogenous rock properties
rock = struct('perm', 1000*darcy*ones(G.cells.num, 1), ...
              'poro', 0.5*ones(G.cells.num, 1));

% Default oil-water fluid with unit values
fluid = initSimpleADIFluid('phases', 'W');
fluid.theta = @(p) getTheta(p, 0.026, 0.42, 0.95, 2.9);
fluid.Kmult = @(p, theta) getConductivity(p, theta, 2.9, .8);   %.8 is K_s probably unrealistic number, required in def of K
fluid.theta_first = @(p) getThetafirst(p, 0.026, 0.42, 0.95, 2.9);

% Set up model and initial state.
model = RichardsEquationModel(G, rock, fluid);
state0 = initResSol(G, -1);

% Set up drive mechanism: constant rate at x=0, constant pressure at x=L
time = 365*day;
n = X; %100
dT = time/n;

% pv = poreVolume(G, rock);
% injRate = -sum(pv)/(time);

bc = [];
% bc = fluxside(bc, G, 'xmin', -injRate, 'sat', [1]);
bc = pside(bc, G, 'zmax', 10, 'sat', [1]);

schedule = simpleSchedule(repmat(dT,1,n), 'bc', bc);



model.nonlinearTolerance = 1e-4;
%%
[~,states, report] = simulateScheduleAD(state0, model, schedule);

%% L scheme
modelFPLScheme = RichardsEquationModelLScheme(G, rock, fluid, 'L', 5);  
modelFPLScheme.nonlinearTolerance = model.nonlinearTolerance;

nls = NonLinearSolver('maxIterations', 10000);
[~,statesFPLScheme, reportFPLScheme] = simulateScheduleAD(state0, modelFPLScheme, schedule, 'nonlinearsolver', nls);

%% Picard scheme
modelFPPicard = RichardsEquationModelPicardScheme(G, rock, fluid);
modelFPPicard.nonlinearTolerance = model.nonlinearTolerance;

nls = NonLinearSolver('maxIterations', 10000);
[~,statesFPPicard, reportFPPicard] = simulateScheduleAD(state0, modelFPPicard, schedule, 'nonlinearsolver', nls);

%% Plotting section
%
%% Plot the result Newton
mrstModule add mrst-gui
close all
plotToolbar(G, states)
title('Newton')

%% Plot results L scheme
figure;
plotToolbar(G, statesFPLScheme)
title('L-scheme')

%% Plot results Picard scheme
figure;
plotToolbar(G, statesFPPicard)
title('Picard-scheme')

%% Plot of number iterations
figure;
plot([report.Iterations, reportFPLScheme.Iterations, reportFPPicard.Iterations])
title('Logarithmic plot of number iteration')
legend('Newton', 'L-scheme', 'Picard-scheme')
set(gca, 'YScale', 'log')

%% Plot number iteration Newton and Picard
figure;
plot([report.Iterations, reportFPPicard.Iterations])
legend('Newton', 'Picard-scheme')
%}
%% study of the errors
%{
e_L(o) = norm(states{X}.pressure- statesFPLScheme{X}.pressure);
e_P(o) = norm(states{X}.pressure- statesFPPicard{X}.pressure);

o=o+1;
%}
end
%% plot of errors
%{
figure;
plot([e_L, e_P])
legend('Newton', 'Picard-scheme')
%}