% Construct 3D grid with 50 cells in the x-direction
mrstModule add ad-blackoil ad-core ad-props mrst-gui blackoil-sequential
clear all
close all
clc
mrstVerbose on

gravity on
X=10;
nx = X;
ny = X;
G = cartGrid([X X], [1 1]*meter); % 50, 50 100 100
G = computeGeometry(G);

% Homogenous rock properties
rock = struct('perm', 1*darcy*ones(G.cells.num, 1), ...
              'poro', .5*ones(G.cells.num, 1)); % perm 1000*darcy*ones(G.cells.num, 1) poro 0.5*ones(G.cells.num, 1)

% Default oil-water fluid with unit values
fluid = initSimpleADIFluid('phases', 'W');
fluid.theta = @(p, c) getThetaCoupled(p, c, 0.026, 0.42, 0.95, 2.49, 73, 2.9);  % (p, c, theta_R, theta_S, alpha, zeta, sigma, n)
fluid.Kmult = @(p, theta) getConductivity(p, theta, .12, 2.9);   %.8 is K_s probably unrealistic number, required in def of K
%fluid.theta_first = @(p) getThetafirst(p, 0.026, 0.42, 0.95, 2.9);

%% Fixed point solver
modelFP = RichardsTransportEquationFixedPointSchemes(G, rock, fluid, 'Newton', 1 );    %'Mixed', 1, 'L_p', .2, 'L_c', .01 

if modelFP.Mixed == 1    
                prompt = 'Number of minimum iteration with L scheme ';
                modelFP.n = input(prompt);
                prompt = 'Introduce L_p: '
                modelFP.L_p = input(prompt);
                prompt = 'Introduce L_c: '
                modelFP.L_c = input(prompt);
end

if modelFP.LScheme == 1    
                prompt = 'Introduce L_p: '
                modelFP.L_p = input(prompt);
                prompt = 'Introduce L_c: '
                modelFP.L_c = input(prompt);
end


%% Compute external forces, assume f(x,y) = xy 
x = G.cells.centroids(:,1);
y = G.cells.centroids(:,2);

% Set up model and initial state.
%model = RichardsTransportEquationModel(G, rock, fluid);
[ii, jj] = gridLogicalIndices(G);

lower = jj < X*1/4;  %ground zone above is vadose zone
p = -3*ones(G.cells.num, 1);
p(lower) = - y(lower) - 3/4;



%{
rate = 0.01*meter^3/second;
bhp  = 100*barsa;
radius = 1/10*1/X;
W = [];
W = verticalWell(W, G, rock, nx - 5, ny - 5, [],     ...
            'Type', 'bhp' , 'Val', bhp, ...
            'Radius', radius, 'InnerProduct', 'ip_tpf', ...
            'Comp_i', [0, 1]);
        %disp('Well #1: '); display(W(1));
%}


state0 = initResSol(G, 0); %- y - 3/4    -1
state0.c = 0*ones(G.cells.num,1);

%state0.wellSol  = initWellSolAD([], modelFP, state0);

% Set up drive mechanism: constant rate at x=0, constant pressure at x=L
time = 1;
n = 20;
dT = time/n;


%f = x.*y;
% Plot force f on domain G
%plotCellData(G,f);
%clear f;
%% Introduce dependency from the time, f(x,y,t) = xyt

%plotCellData(G,f);

%% Boundary conditions
bc = [];
bc = pside(bc, G, 'ymin', 0, 'sat', [1]);
bc = pside(bc, G, 'xmin', 0, 'sat', [1]);%-3
bc.c = 0*ones(size(bc.sat,1), 1);%4.*ones(size(bc.sat,1), 1);


schedule = simpleSchedule(repmat(dT,1,n),'bc', bc)


model.nonlinearTolerance = 1e-7; %-4
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
modelFP.nonlinearTolerance = 1e-7;
for i=1:n
   
   t = i*dT;
 
   analyticalForce1 = 2*t.*x.^2.*y.^2 - 2*t^3.*x.*y.^3 - 2*t^3.*x.^3.*y - 2*t^2.*x.^2.*y 
   
   analyticalForce2 = 3*t^2.*x.^3.*y.^3 - 2*t^3.*x.*y.^3 -2*t^3.*x.^3.*y
   
   modelFP.forces(i).analyticalForce1 = analyticalForce1; % original: .006.*cos(4/3*pi.*y).*sin(x);
   modelFP.forces(i).analyticalForce2 = analyticalForce2;
   
   %solution_a = t.*x;
   
   %schedule.control(i).forcesAnal = i.*dT*x.*y;
end

nls = NonLinearSolver('maxIterations', 10000);
tic
[~,statesFP, reportFP] = simulateScheduleAD(state0, modelFP, schedule, 'nonlinearsolver', nls);
toc



%{
%% Plot the result
%% Plot the result Newton
mrstModule add mrst-gui
close all
plotToolbar(G, states)
title('Newton')
colorbar
%}
%h = colorbar;
%set(h, 'ylim', [1 3]) 
%% Plot results modified Picard
figure;
plotToolbar(G, statesFP)
title('Newton not coupled theta')
colorbar

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