% Construct 3D grid with 50 cells in the x-direction
mrstModule add ad-blackoil ad-core ad-props mrst-gui blackoil-sequential
clear all
close all
clc
mrstVerbose on

gravity reset on
X = 10;
G = cartGrid([X X], [1 1]*meter); % 50, 50 100 100
G = computeGeometry(G);

% Homogenous rock properties
rock = struct('perm', 1*ones(G.cells.num, 1), ...
              'poro', .5*ones(G.cells.num, 1)); % perm 1000*darcy*ones(G.cells.num, 1) poro 0.5*ones(G.cells.num, 1)

% Default oil-water fluid with unit values
fluid = initSimpleADIFluid('phases', 'W');
fluid.theta = @(p, c) getThetaCoupled(p, c, 0.026, 0.42, 0.551, 2.9, .44, .0046);  % (p, c, theta_R, theta_S, alpha, zeta, sigma, n)
fluid.Kmult = @(p, theta) getConductivity(p, theta, .12, 2.9);   %.8 is K_s probably unrealistic number, required in def of K


%% Compute external forces, assume f(x,y) = xy 
x = G.cells.centroids(:,1);
y = G.cells.centroids(:,2);

% Set up model and initial state.c
%model = RichardsTransportEquationModel(G, rock, fluid);
[ii, jj] = gridLogicalIndices(G);

lower = jj < X*1/4;
p = -2*ones(G.cells.num, 1);
p(lower) = - y(lower) - 1/4;
%s = [sW, 1 - sW];

state0 = initResSol(G, p); % - y - 3/4   -1
state0.c = 0*ones(G.cells.num,1);


% Set up drive mechanism: constant rate at x=0, constant pressure at x=L
time = 1; %365*day;
n = 10;
dT = time/n;


%f = x.*y;
% Plot force f on domain G
%plotCellData(G,f);
%clear f;
%% Introduce dependency from the time, f(x,y,t) = xyt

%plotCellData(G,f);

%% Boundary conditions
bc = [];
%bc = fluxside(bc, G, 'xmin', -injRate, 'sat', [1]);
bc = pside(bc, G, 'ymax', -2, 'sat', [1]);%10  %in paper xmin -3 or -2
bc.c = 4.*ones(size(bc.sat,1), 1); %4.*ones(size(bc.sat,1), 1);

schedule = simpleSchedule(repmat(dT,1,n), 'bc', bc);


model.nonlinearTolerance = 1e-6;
nls = NonLinearSolver('maxIterations', 10000);


% Set up model and initial state.
Richardsmodel = RichardsEquationModelSequential(G, rock, fluid, 'Newton', .1);
Transportmodel = SequentialTransportEquationModel(G, rock, fluid, 'Newton', .005);
seqmodel = SequentialRichardsTransportModel(Richardsmodel, Transportmodel);

[seqmodel.forces] = getValidDrivingForces(seqmodel);

for i=1:n
    force1 = .006.*cos(-4/3*pi.*y).*sin(x);
    force1(lower) = 0;
   seqmodel.forces(i).analyticalForce1 = force1;
   seqmodel.forces(i).analyticalForce2 = force1;
   
end

seqmodel.nonlinearTolerance = model.nonlinearTolerance;

tic
[~,seqstates, seqreport] = simulateScheduleAD(state0, seqmodel, schedule, 'nonlinearsolver', nls);
toc
disp('Total Number iterations:')
sum(seqreport.Iterations)

%
%% Plot the result Newton
%{
mrstModule add mrst-gui
close all
figure
plotToolbar(G, seqstates)
title('Newton')
colorbar

%}