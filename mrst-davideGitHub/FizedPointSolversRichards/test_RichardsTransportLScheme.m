%% Version from SourceTree
%Construct 3D grid with 50 cells in the x-direction
mrstModule add ad-core ad-props
gravity reset on
clear all
close all
clc

o=1;
%for X = [10 20 40 80]

G = cartGrid([50, 1, 50], [100, 100, 10]*meter);
G = computeGeometry(G);

% Homogenous rock properties
rock = struct('perm', 1000*darcy*ones(G.cells.num, 1), ...
              'poro', 0.5*ones(G.cells.num, 1));

% Default oil-water fluid with unit values
fluid = initSimpleADIFluid('phases', 'W');
fluid.theta = @(p) getTheta(p, 0.026, 0.42, 0.95, 2.9);
fluid.Kmult = @(p, theta) getConductivity(p, theta, 2.9, .8);

% Set up model and initial state.
model = RichardsTransportEquationModel(G, rock, fluid);
state0 = initResSol(G, -1);
state0.c = zeros(G.cells.num,1);

% Set up drive mechanism: constant rate at x=0, constant pressure at x=L
time = 365*day;
n = 100;  % it was 100
dT = time/n;

% pv = poreVolume(G, rock);
% injRate = -sum(pv)/(time);

bc = [];
% bc = fluxside(bc, G, 'xmin', -injRate, 'sat', [1]);
bc = pside(bc, G, 'zmax', 10, 'sat', [1]);
bc.c = 4.*ones(size(bc.sat,1), 1);


schedule = simpleSchedule(repmat(dT,1,n), 'bc', bc);

model.nonlinearTolerance = 1e-4;
%%
[~,states, report] = simulateScheduleAD(state0, model, schedule);

%%
modelFP = RichardsTransportEquationModelLScheme(G, rock, fluid, 'L', 10);
modelFP.nonlinearTolerance = model.nonlinearTolerance;

nls = NonLinearSolver('maxIterations', 10000);
[~,statesFP, reportFP] = simulateScheduleAD(state0, modelFP, schedule, 'nonlinearsolver', nls);

%{
%% Plot the result the result Newton method
mrstModule add mrst-gui
close all
plotToolbar(G, states)
title('Newton')

%% Plot the result L scheme
figure;
plotToolbar(G, statesFP)
title('L-scheme')


%%
figure;
plot([report.Iterations, reportFP.Iterations])
legend('Newton', 'L-scheme')
set(gca, 'YScale', 'log')
%}

% computation of fidderence between results
e(o) = norm(states{X}.pressure - statesFP{X}.pressure);
o=o+1;

%end
