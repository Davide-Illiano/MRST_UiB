% Construct 3D grid with 50 cells in the x-direction
mrstModule add ad-blackoil ad-core ad-props mrst-gui %blackoil-sequential diagnostics
clear all
close all
clc
mrstVerbose on
gravity off

% Create grid
G = cartGrid([50 20]); % 50, 50 100 100
G = computeGeometry(G);

% Resize rock
poro = gaussianField(G.cartDims, [.2 .6],[11 3], 3.5);
K = poro.^3.*(1e-5)./(.81*72*(1-poro).^2);
rock = makeRock(G, K(:),poro(:));


%% inspect
%
% show poro
figure;
plotCellData(G,rock.poro,'EdgeColor', 'none');
colorbar('horiz'); axis equal tight;
title('Rock porosity')

%show perm
figure;
plotCellData(G,rock.perm);
colorbar('horiz'); axis equal tight;
title('Rock permeability');
%}

% Default oil-water fluid with unit values
fluid = initSimpleADIFluid('phases', 'W');
fluid.theta = @(p,c) getThetaCoupled(p, c, 0.026, 0.42, 1/(50*barsa).*.8, 2.9, .44, .0046);  % (p, c, theta_R, theta_S, alpha,  n, a,b)
fluid.Kmult = @(p, theta) 1+0.*theta+0.*p; %getConductivity(p, theta, .12, 2.9);   %.8 is K_s probably unrealistic number, required in def of K


%% Fixed point solver
modelFP = RichardsTransportEquationFixedPointSchemes(G, rock, fluid, 'Newton', 1 );    %'Mixed', 1, 'L_p', .2, 'L_c', .01 

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
                modelFP.L_p = 1.5;
                %prompt = 'Introduce L_c: '
                %modelFP.L_c = input(prompt);
                modelFP.L_c = 0;
end


state0 = initResSol(G, -100*barsa, [1]); 
state0.c = 0*ones(G.cells.num,1); %0*
%state0.wellSol  = initWellSolAD([], modelFP, state0);

% Set up drive mechanism: constant rate at x=0, constant pressure at x=L
time = 50*day;
n = 50;
dT = time/n;
%inspectFluidModel(modelFP)

%% Introduce dependency from the time, f(x,y,t) = xyt

%plotCellData(G,f);

%% Boundary conditions
pv = sum(poreVolume(G,rock));
src = addSource([], [1],1e-2,'sat',1);
src.c = 0;
src = addSource(src, [1000], -1e-2, 'sat',1);
src.c = 0; 
%bc = [];
%bc = pside(bc, G, 'xmin', -1, 'sat', [1]);%-3
%bc.c = 0*ones(size(bc.sat,1), 1);

schedule = simpleSchedule(repmat(dT,1,n),'src', src);



[modelFP.forces] = getValidDrivingForces(modelFP);
modelFP.nonlinearTolerance = 1e-4;


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