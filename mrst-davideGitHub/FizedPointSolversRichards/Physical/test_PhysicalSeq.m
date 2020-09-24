% Construct 3D grid with 50 cells in the x-direction
mrstModule add ad-blackoil ad-core ad-props mrst-gui blackoil-sequential
clear all
close all
clc
mrstVerbose on

gravity reset on
gravity y
X = 40;

nx = X;
ny = X;
G = cartGrid([X X], [2 3]*meter); % 50, 50 100 100
G = computeGeometry(G);

% Homogenous rock properties

% Default oil-water fluid with unit values
fluid = initSimpleADIFluid('phases', 'W');
K_s = 8.2e-4/day;
fluid.theta = @(p, c) getThetaCoupled(p, c, 0, 0.446, 0.152, 1.17);  % (p, c, theta_R, theta_S, alpha,  n, a,b)
fluid.Kmult = @(p, theta) getConductivity(p, theta, K_s, 1.17);   %.8 is K_s probably unrealistic number, required in def of K

rock = struct('perm', K_s*ones(G.cells.num, 1), ... % should be 8.2e-4
              'poro', .5*ones(G.cells.num, 1)); 
          
          
% Set up model and initial state.
Richardsmodel = RichardsEquationModelSequential(G, rock, fluid, 'L_p', .00001);
Transportmodel = SequentialTransportEquationModel(G, rock, fluid, 'L_c', .05);
seqmodel = SequentialRichardsTransportModel(Richardsmodel, Transportmodel);

x = G.cells.centroids(:,1);
y = G.cells.centroids(:,2);

% Set up model and initial state.
%model = RichardsTransportEquationModel(G, rock, fluid);
[ii, jj] = gridLogicalIndices(G);

p = 1-y;  

state0 = initResSol(G, p); 
state0.c = 0*ones(G.cells.num,1); 
%state0.wellSol  = initWellSolAD([], modelFP, state0);

% Set up drive mechanism: constant rate at x=0, constant pressure at x=L
time = 3*day;
n = 10;
dT = time/n;
%inspectFluidModel(modelFP)

%f = x.*y;
% Plot force f on domain G
%plotCellData(G,f);
%clear f;

%% Boundary conditions
z = linspace(0,1,floor(X/3));
p_z = 1-z;

bfaces  = find(any(G.faces.neighbors==0,2));
%{
figure
hold on
plotFaces(G, bfaces(1:2:2*X), 'Linewidth', 3); %left
plotFaces(G, bfaces(2:2:2*X+1),'r', 'Linewidth', 3); %right
plotFaces(G, bfaces(2:2:2*floor(X/3)),'r', 'Linewidth', 8); %right lower
plotFaces(G, bfaces(2*X+1:3*X), 'y', 'Linewidth', 3); %down
plotFaces(G, bfaces(3*X+1:3*X+floor(X/2)), 'y', 'Linewidth', 8); %top left
plotFaces(G, bfaces(3*X+floor(X/2)+1:4*X), 'y', 'Linewidth',3); %top right
plotFaces(G, bfaces(3*X+1:4*X),'g', 'Linewidth', 3); %up

grid on
title('Boundaries');
%}

z = linspace(0,1,size(bfaces(2:2:2*floor(X/3)-1),1));
p_z = 1-z;

for j = 1:n
bc = [];
%{
% no flow conditions left
bc = addBC(bc, bfaces(1:2:2*X), 'flux', 0, 'sat', [1]);
bc.c = 0*ones(size(bc.sat,1), 1);
% no flow conditions down
bc = addBC(bc, bfaces(2*X+1:3*X), 'flux', 0, 'sat', [1]);
bc.c = 0*ones(size(bc.sat,1), 1);


% Mixed conditions on up and right sides
bc = addBC(bc, bfaces(3*X+floor(X/2)+1: 4*X), 'flux', 0, 'sat', [1]); % up no flux
bc.c = 0*ones(size(bc.sat,1), 1);
%}
%bc = pside(bc, G, 'ymax', -2, 'sat', [1]);
bc = addBC(bc, bfaces(2:2:2*floor(X/3)-1), 'pressure', p_z, 'sat', [1]); % pressure right
bc.c = 0*ones(size(bc.sat,1), 1);

if j*dT <= 1*day
    bc = addBC(bc, bfaces(3*X+1:3*X+floor(X/2)), 'pressure', (-2 + 2.2*(j*dT/day)), 'sat', [1]); % up pressure
    bc.c = 0*ones(size(bc.sat,1), 1);
else
    bc = addBC(bc, bfaces(3*X+1:3*X+floor(X/2)), 'pressure', .2, 'sat', [1]); % up pressure
    bc.c = 0*ones(size(bc.sat,1), 1);  
end
%{
bc = addBC(bc, bfaces(2*floor(X/3):2: 2*X+1), 'flux', 0, 'sat', [1]); % right no flux
bc.c = 0*ones(size(bc.sat,1), 1);
%}

BC(j) = bc;
end

schedule = simpleSchedule(repmat(dT,1,n),'bc', BC)

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


[seqmodel.forces] = getValidDrivingForces(seqmodel);
seqmodel.nonlinearTolerance = 1e-6;

nls = NonLinearSolver('maxIterations', 10000);
tic
[~,statesFP, reportseq] = simulateScheduleAD(state0, seqmodel, schedule, 'nonlinearsolver', nls);
toc

disp('Total number iterations:')
sum(reportseq.Iterations)
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
figure;
plotToolbar(G, statesFP)
title('Richards + Transport')
colorbar;

%{
figure('Units','inches',...
'PaperPositionMode','auto');
grid on, 
plotCellData(G,statesFP{10}.c),...
title('Concentration');
set(gca,...
'Units','normalized',...
'FontUnits','points',...
'FontWeight','normal',...
'FontSize',13,...
'FontName','Times')
colorbar

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