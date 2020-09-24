%% Classic Buckley-Leverett Problem: 1D Horizontal Displacement
% permeabilities that follow a standard quadratic Corey model
clear all
mrstModule add ad-blackoil ad-core ad-props mrst-gui deckformat ad-core ad-blackoil ad-eor ad-fi ad-props mrst-gui
mrstVerbose off
%% Set up model

% Construct 3D grid with 50 cells in the x-direction
ncells = 50;
G = cartGrid([ncells, 1, 1], [ncells, 1, 1]*meter);
G = computeGeometry(G);

% Homogenous rock properties
rock = struct('perm', ones(G.cells.num, 1), ...
              'poro', ones(G.cells.num, 1));

% Default oil-water fluid with unit values
%fluid = initSimpleADIFluid('phases', 'WO', 'n', [2 2]);

fluid = struct('phases', 'WO', ... %'mu', [1, 10]*centi*poise, ...
                           'pcOW'      ,@(sw) 1-sw.^2, ...  
                           'krW'       ,@(sw) (1-sw).^(1/2).*(1-sw.^(1/(1/2))), ...   % from Van Genuchten
                           'krOW'      ,@(sw) (sw).^(1.2).*(1-(1-sw.^(2)).^(1/2)).^2, ...
                           'sWcon'     ,.15, ...
                           'rhoOS'     ,500, ...
                           'rhoWS'     ,1000, ...
                           'cW'        ,4.28e-10, ...
                           'muWr'      ,4.8e-04, ...
                           'bW'        ,@(pw) pw.*0 + 1,...   % multiplier density water
                           'muW'       ,@(pw) pw.*0 + 1*centi*poise,...
                           'cO'        ,6.65e-10, ...
                           'bO'        ,@(po) po.*0 + 1, ...  % multiplier oil density
                           'BOxmuO'    ,@(po) po.*0 + 1, ...  %not sure what is this
                           'muO'       ,@(po) po.*0 + 10*centi*poise,...
                           'cR'        ,3e-10, ...
                           'pvMultR'   ,@(p) (p.*0 + 1), ...  % multiplier due to pressure on pores
                           'muWMult'   ,@(c) c.*0+1, ...      % used in multiplier due to polymer
                           'dps'       ,0, ...                % used in eq. for conservation of polymer
                           'rrf'       ,1.1, ...              % PROBLEM NOT ABLE TO SET TO 1
                           'rhoR'      ,2000, ...
                           'adsInx'    ,1, ...                % if set equal 2 adsorption is function of c_max
                           'adsMax'    ,5e-5, ...             % PROBLEM, NOT ABLE TO SET TO ZERO
                           'ads'       ,@(c) c.*0 + 1, ...    % adsorption
                           'mixPar'    ,1, ...                % being equal 1 multiplier due to polymer are equal 1
                           'cmax'      ,4 ...                 % set equal to 4 in the example
                            ); 

% Set up model and initial state.
%model = TwoPhaseOilWaterModel(G, rock, fluid);
model = OilWaterPolymerModel(G, rock, fluid);

state0 = initResSol(G, 0, [0.5, 0.5]);
% if polymer model
state0.c    = zeros(G.cells.num,1);
state0.cmax = zeros(G.cells.num,1);

state0.wellSol = initWellSolAD([], model, state0);

% Set up drive mechanism: constant rate at x=0, constant pressure at x=L
% time = 1;
n = 50;
time = n;

pv = poreVolume(G, rock);
injRate = -sum(pv)/time;
bc = fluxside([], G, 'xmin', -injRate, 'sat', [1, 0]);
bc = pside(bc, G, 'xmax', 0*barsa, 'sat', [0, 1]);
%if polymer
bc.c = 4.*ones(size(bc.sat,1), 1); 

dT = repmat(time/n, n, 1);
%% Simulation with general solver
schedule = simpleSchedule(dT, 'bc', bc);

%[~,sstates] = simulateScheduleAD(state0, model, schedule);
%%
%close all
%plotToolbar(G, sstates, 'field', 's:1','lockCaxis',true), 
%caxis([0 1]), view(10,10)
%colorbar
%%
%
%fmodel = FixedPointModelPicardScheme(model);
%fmodel = FixedPointModelLScheme(model);
fmodel = FixedPointModelNewton(model);

%fmodel.parentModel.dsMaxAbs = inf;
%nls = NonLinearSolver('maxIterations', 100);
%tic
%[~, fstates] = simulateScheduleAD(state0, fmodel, schedule, 'nonlinearsolver', nls);
%toc
%}
%%

mrstModule add blackoil-sequential
tic
seqModel = getSequentialModelFromFI(model);

seqModel.transportModel.conserveWater = true;
seqModel.transportModel.conserveOil = ~seqModel.transportModel.conserveWater;

seqModel.transportModel = FixedPointModelNewton(seqModel.transportModel);
seqModel.transportNonLinearSolver.maxIterations = 500;

%seqModel.pressureModel = FixedPointModel(seqModel.pressureModel);
%seqModel.pressureNonLinearSolver.maxIterations = 5000;

[~, fsstates] = simulateScheduleAD(state0, seqModel, schedule, 'nonlinearsolver', nls);
toc
close all
plotToolbar(G, fsstates, 'field', 's:1','lockCaxis',true), 
caxis([0 1]), view(10,10)
colorbar