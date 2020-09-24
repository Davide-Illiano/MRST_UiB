%% Classic Buckley-Leverett Problem: 1D Horizontal Displacement
% permeabilities that follow a standard quadratic Corey model
%% Example built to implement a fixed point solver
clear all
mrstModule add ad-blackoil ad-core ad-props mrst-gui blackoil-sequential
mrstVerbose off
tic

%% Set up model Easy example
%{
% Construct 3D grid with 50 cells in the x-direction
ncells = 50;
G = cartGrid([ncells, 1, 1], [ncells, 1, 1]*meter);
G = computeGeometry(G);

% Homogenous rock properties
rock = struct('perm', ones(G.cells.num, 1), ...
              'poro', ones(G.cells.num, 1));

% Default oil-water fluid with unit values
fluid = initSimpleADIFluid('phases', 'WO', 'n', [2 2]);

% Set up model and initial state.
model = TwoPhaseOilWaterModel(G, rock, fluid);
state0 = initResSol(G, 0, [0.5, 0.5]);
state0.wellSol = initWellSolAD([], model, state0);

% Set up drive mechanism: constant rate at x=0, constant pressure at x=L
% time = 1;
n = 50;
time = n;

pv = poreVolume(G, rock);
injRate = -sum(pv)/time;
bc = fluxside([], G, 'xmin', -injRate, 'sat', [1, 0]);
bc = pside(bc, G, 'xmax', 0*barsa, 'sat', [0, 1]);

dT = repmat(time/n, n, 1);
%}

%% Set up BuckleyL
%
% Construct 3D grid with 50 cells in the x-direction
G = cartGrid([50, 1, 1], [1000, 10, 10]*meter);
G = computeGeometry(G);

% Homogenous rock properties
rock = struct('perm', darcy*ones(G.cells.num, 1), ...
              'poro', .3*ones(G.cells.num, 1));

% Default oil-water fluid with unit values
fluid = initSimpleADIFluid('phases', 'WO', 'n', [2 2]);

% Set up model and initial state.
model = TwoPhaseOilWaterModel(G, rock, fluid);
state0 = initResSol(G, 50*barsa, [0, 1]);
state0.wellSol = initWellSolAD([], model, state0);

% Set up drive mechanism: constant rate at x=0, constant pressure at x=L
pv = poreVolume(G, rock);
injRate = -sum(pv)/(500*day);
bc = fluxside([], G, 'xmin', -injRate, 'sat', [1, 0]);
bc = pside(bc, G, 'xmax', 0*barsa, 'sat', [0, 1]);

% Time schedule
dT = 20*day;
n = 25;
%}

%% Simulation with general solver
%schedule = simpleSchedule(dT, 'bc', bc);
%% BuckleyL
schedule = simpleSchedule(repmat(dT,1,25), 'bc', bc);

%% Call the Fixed Point Model appositely defined 
fmodel = FixedPointModelLScheme(model);

%% Use sequential solvers
seqModel = getSequentialModelFromFI(model);

seqModel.transportModel.conserveWater = true;
seqModel.transportModel.conserveOil = ~seqModel.transportModel.conserveWater;

%% Implement L Scheme
%seqModel.transportModel = FixedPointModelLScheme(seqModel.transportModel);

%% Implement Picard Method
seqModel.transportModel = FixedPointModelPicardScheme(seqModel.transportModel);

nls = NonLinearSolver('maxIterations', 100);
[~, fsstates] = simulateScheduleAD(state0, seqModel, schedule, 'nonlinearsolver', nls);

% Computational time to run the code
time = toc
%% Plot the results obtained
close all
plotToolbar(G, fsstates, 'field', 's:1','lockCaxis',true), 
caxis([0 1]), view(10,10)
colorbar
