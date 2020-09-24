% Construct 3D grid with 50 cells in the x-direction
mrstModule add ad-blackoil ad-core ad-props mrst-gui %blackoil-sequential diagnostics
clear all
close all
clc
mrstVerbose on

gravity on
X = 10;
nx = X;
ny = X;
G = cartGrid([X X], [1 1]*meter); % 50, 50 100 100
G = computeGeometry(G);

% Homogenous rock properties
rock = struct('perm', ones(G.cells.num, 1), ...
              'poro', .5*ones(G.cells.num, 1)); % perm 1000*darcy*ones(G.cells.num, 1) poro 0.5*ones(G.cells.num, 1)

% Default oil-water fluid with unit values
fluid = initSimpleADIFluid('phases', 'W');
fluid.theta = @(p, c) getThetaCoupled(p, c, 0.026, 0.42, 0.551, 2.9, .44, .0046);  % (p, c, theta_R, theta_S, alpha,  n, a,b)
fluid.Kmult = @(p, theta) getConductivity(p, theta, .12, 2.9);   %.8 is K_s probably unrealistic number, required in def of K
%fluid.theta_first = @(p) getThetafirst(p, 0.026, 0.42, 0.95, 2.9);

%% Fixed point solver
modelFP = RichardsTransportEquationFixedPointSchemesnewBC(G, rock, fluid, 'Newton', 1 );    %'Mixed', 1, 'L_p', .2, 'L_c', .01 

if modelFP.Mixed == 1    
                prompt = 'Number of minimum iteration with L scheme ';
                modelFP.n = input(prompt);   %suggested .04
                %prompt = 'Introduce L_p: '
                %modelFP.L_p = input(prompt); %suggested .07
                modelFP.L_p = .15;
                %prompt = 'Introduce L_c: '
                %modelFP.L_c = input(prompt);
                modelFP.L_c = .005;
end

if modelFP.LScheme == 1    
                %prompt = 'Introduce L_p: '
                %modelFP.L_p = input(prompt); %suggested .07
                modelFP.L_p = .1;
                %prompt = 'Introduce L_c: '
                %modelFP.L_c = input(prompt);
                modelFP.L_c = .1;
end

%% Compute external forces, assume f(x,y) = xy 
x = G.cells.centroids(:,1);
y = G.cells.centroids(:,2);

% Set up model and initial state.
%model = RichardsTransportEquationModel(G, rock, fluid);
[ii, jj] = gridLogicalIndices(G);

lower = jj <= X*1/4;  % vadose zone
higher = jj > X*1/4;
p = -2*ones(G.cells.num, 1); %-2
p(lower) = -y(lower) - 1/4;  % correct is -y(lower) + 1/4


state0 = initResSol(G, p); 
state0.c = 1*ones(G.cells.num,1); %0*
%state0.wellSol  = initWellSolAD([], modelFP, state0);

% Set up drive mechanism: constant rate at x=0, constant pressure at x=L
time = 1;
n = 10;
dT = time/n;
%inspectFluidModel(modelFP)

%f = x.*y;
% Plot force f on domain G
%plotCellData(G,f);
%clear f;
%% Introduce dependency from the time, f(x,y,t) = xyt

%plotCellData(G,f);

%% Boundary conditions
bc = [];
bc = pside(bc, G, 'ymax', -2, 'sat', [1]);%-3
bc.c = 4*ones(size(bc.sat,1), 1);%4.*ones(size(bc.sat,1), 1);

schedule = simpleSchedule(repmat(dT,1,n),'bc', bc)

model.nonlinearTolerance = 1e-6; %-4
%%
%{
tic
[~,states, report] = simulateScheduleAD(state0, model, schedule);
toc
%% L scheme
modelFPLScheme = RichardsTransportEquationModelLScheme(G, rock, fluid, 'L_p', .07, 'L_c', .8);  
modelFPLScheme.nonlinearTolerance = model.nonlinearTolerance;
nls = NonLinearSolver('maxIterations', 10000);
[~,statesFPLScheme, reportFPLScheme] = simulateScheduleAD(state0, modelFPLScheme, schedule, 'nonlinearsolver', nls);
%}


[modelFP.forces] = getValidDrivingForces(modelFP);
modelFP.nonlinearTolerance = 1e-6;
%% Introduce external forces
%
for i=1:n
   analyticalForce1 = .006.*cos(-4/3*pi.*y).*sin(x);   %we want external forces in vadose zone
   analyticalForce1(lower) = 0;       % lower = jj > X*1/4;
   
   analyticalForce2 = .006.*cos(-4/3*pi.*y).*sin(x);
   analyticalForce2(lower) = 0;       % lower = jj > X*1/4;
   
   modelFP.forces(i).analyticalForce1 = analyticalForce1; % original: .006.*cos(4/3*pi.*y).*sin(x);
   modelFP.forces(i).analyticalForce2 = analyticalForce2;
   
end
%}

nls = NonLinearSolver('maxIterations', 10000);
tic
[~,statesFP, reportFP] = simulateScheduleAD(state0, modelFP, schedule, 'nonlinearsolver', nls);
toc

disp('Total number iterations:')
sum(reportFP.Iterations)
%{
%% Plot the result
%% Plot the result Newton
mrstModule add mrst-gui
close all
plotToolbar(G, states)
title('Newton')
colorbar
%}

%% Plot results modified Picard
%
figure;
plotToolbar(G, statesFP)
title('Newton method')
colorbar;
%}

%{
hold on
for i = 1:15
   l=figure(1), plotCellData(G,statesFP{i}.pressure),title('Pressure'), colorbar
   saveas(l,sprintf('FIGpressure%d.png',i)); 
   h=figure(2), plotCellData(G,statesFP{i}.c),title('Concentration'), colorbar
   saveas(h,sprintf('FIGconcentration%d.png',i)); 
   k=figure(3), plotCellData(G,statesFP{i}.theta),title('Water content'), colorbar
   saveas(k,sprintf('FIGtheta%d.png',i)); 
end
%}
   %figure(3), plotCellData(G,statesFP{i}.c), title('Concentration'),colorbar

%{
statescoupled = statesFP;
load('notcoupled.mat')

figure
title('Difference between thetas final time step (Couple -non coupled)')
for i = 1:10
    hold on
    plotToolbar(G, statescoupled{i}.theta - statesFP{i}.theta)
end
%}
%{
%% Plot results L scheme
figure;
plotToolbar(G, statesFPLScheme)
title('L-scheme')

%%
figure;
plot([report.Iterations, reportFPPicard.Iterations, reportFPLScheme.Iterations])
legend('Newton','ModifiedPicard', 'L-scheme')
set(gca, 'YScale', 'log')
%}