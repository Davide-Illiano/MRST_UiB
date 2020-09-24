% Construct 3D grid with 50 cells in the x-direction
mrstModule add ad-blackoil ad-core ad-props mrst-gui %blackoil-sequential diagnostics
clear all
close all

% clc
mrstVerbose off

% gravity off
gravity on
gravity([0,-1,0])


X = 101;
Y = 151;
I = X;
J = Y;
a=0; b=2;
c=0; d=3;
x1=a; x2=b;
y1=c; y2=d;


G = cartGrid([X Y], [2 3]);
G = computeGeometry(G);

%% Initial condition
dx = (b-a)/(I-1);
x = [a:dx:b];
x2 = (x(1:I-1)+x(2:I))/2;
dy=(d-c)/(J-1);
Lx=I-2; Ly=J-2;

if gravity == [0,-1,0]
        Umean = -0.0331;
else
        Umean = 0.0165;
end
% Homogenous rock properties
rock = struct('perm', ones(G.cells.num, 1), ...
              'poro', ones(G.cells.num, 1)); % perm 1000*darcy*ones(G.cells.num, 1) poro 0.5*ones(G.cells.num, 1)
clay = 0; % clay = 0 then loam
if  clay == 1       % L_p=10 does not converge at 65 
    K_s = 8.2*1e-4;
    T = 3;
    DeltaT = T/3;
    n = 9;  %3*day * 1/ (1/3)
    theta_r = 0;
    theta_s = .446;
    alpha = .152;
    nGM = 1.17;
    ng = (nGM-1)/nGM; 
else                % L_10 
    K_s = 4.96*1e-2;
    T = 3;       %% change to 3 days!
    DeltaT = T/3; 
    n = 10; % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    theta_s = 1;
    theta_r = .131;
    alpha = .423;
    nGM = 2.06; 
    ng = (nGM-1)/nGM;
end
A=0.44;
B=0.0046; 
    
% Default oil-water fluid with unit values
fluid = initSimpleADIFluid('phases', 'W');
fluid.theta = @(p,c)  theta_s + 0.*p + 0.*c;% getThetaCoupled(p, c, theta_r, theta_s, alpha, ng, nGM, A, B);     %theta_s + 0.*p+ 0.*c;%     ./ (1 - p - 1/10.*c);   % theta is a function of pressure and concentration
fluid.Kmult = @(p,theta) K_s + 0.*p + 0.*theta;% getConductivity(p,theta, theta_r, theta_s, K_s, nGM); %K_s + 0.*p + 0.*theta;


%% Fixed point solver
modelFP = RichardsTransportEquationFixedPointSchemes(G, rock, fluid, 'LScheme', 1 );    %'Mixed', 1, 'L_p', .2, 'L_c', .01 

if modelFP.LScheme == 1    
                %prompt = 'Introduce L_p: '
                %modelFP.L_p = input(prompt); %suggested .07
                modelFP.L_p = .1;
                %prompt = 'Introduce L_c: '
                %modelFP.L_c = input(prompt);
                modelFP.L_c = .1;
end

y = G.cells.centroids(:,2);
x = G.cells.centroids(:,1);

p0 = 1 - y/3; % constant-parameters case!!!
c0 = y/3;

state0 = initResSol(G, p0);
state0.c = c0;


% Set up drive mechanism: constant rate at x=0, constant pressure at x=L
dT = T/n;


%% Boundary conditions
bc = [];
bc = pside(bc,G,'YMax', 0 , 'sat', [1]);
bc.c = ones(size(bc.sat,1),1);
bc = pside(bc,G,'YMin', 1, 'sat', [1]);  % 1
bc.c = [bc.c' 0.*bc.c']';

bfaces  = find(any(G.faces.neighbors==0,2));
%{
figure
hold on
% plotFaces(G, bfaces(1:2:2*Y), 'Linewidth', 1); %left
% plotFaces(G, bfaces(2:2:2*Y+1),'r', 'Linewidth', 1); %right
% plotFaces(G, bfaces(2:2:2*floor(Y/3)),'r', 'Linewidth', 2); %right lower
plotFaces(G, bfaces(2*Y+1:2*Y+X), 'y', 'Linewidth', 3); %down
% plotFaces(G, bfaces(84:93), 'y', 'Linewidth', 2); %top left
%plotFaces(G, bfaces(3*X+floor(X/2)+1:4*X), 'y', 'Linewidth',3); %top right
plotFaces(G, bfaces(2*Y+X+1:2*X+2*Y),'g', 'Linewidth', 3); %up

grid on
title('Boundaries');
return
%}

% volume = mean(G.cells.volumes);
% bc = [];
% bc = addBC(bc,bfaces(2*Y+X+1:2*X+2*Y), 'flux', -Umean*dx, 'sat', [1] ); %TOP
% bc.c = ones(size(bc.sat,1),1);
% bc = addBC(bc,bfaces(2*Y+1:2*Y+X),  'flux', Umean*dx, 'sat', [1]); %BOTTOM
% bc.c = [bc.c' 0.*bc.c']';

schedule = simpleSchedule(repmat(dT,1,n),'bc', bc);


[modelFP.forces] = getValidDrivingForces(modelFP);
modelFP.nonlinearTolerance = 1e-6;

nls = NonLinearSolver('maxIterations', 10000);

tic
[~,statesFP, reportFP] = simulateScheduleAD(state0, modelFP, schedule, 'nonlinearsolver', nls);
toc

disp('Total number iterations:')
sum(reportFP.Iterations)

% create vector with all the condition numbers
o=1;
%%
%}
%% Plot the result
%% Plot the result Newton
%{
mrstModule add mrst-gui
%close all
plotToolbar(G, statesFP)
title('Monolithic')
set(gca,...
'Units','normalized',...
'FontUnits','points',...
'FontWeight','normal',...
'FontSize',19,...
'FontName','Times')
colorbar
%}

%% Plot results modified Picard
%
% figure;
% plotToolbar(G, statesFP)
% title('Final time step')
% colorbar;

states = statesFP;

p = zeros(J,I);
v = faceFlux2cellVelocity(G,states{n}.flux);
theta = p;
c=p;

for l = 1:J
    p_mrst(l,:) = states{n}.pressure(((l-1)*I+1:l*I));
    theta(l,:) = states{n}.theta(((l-1)*I+1:l*I));
    c_mrst(l,:) = states{n}.c(((l-1)*I+1:l*I));
    vx(l,:) = v(((l-1)*I+1:l*I),1);
    vy(l,:) = v(((l-1)*I+1:l*I),2);
end

figure
plot(G.cells.centroids(:, 2), states{10}.c, '.')
title('Final concentration')
hold on 
plot(y,y/3)

figure
plot(G.cells.centroids(:, 2), states{10}.pressure, '.')
title('Final pressure')
hold on 
plot(y,1-y/3)

figure
plot(states{10}.flux)
title('Fluxes')

I = X; J=Y; %I = 41; J=61; % 
a=0; b=2;
c=0; d=3;
dx = (b-a)/(I-1)
x = [a:dx:b];
x2 = (x(1:I-1)+x(2:I))/2;
dy=(d-c)/(J-1)
y=c:dy:d;
levels=12;
plot_contours(I,J,x,y,p_mrst,c_mrst,theta,vx,vy,levels);

disp('Norm errror of mean vy')
norm(vy-Umean)

disp('Norm errror of mean vy')
norm(vy-mean(mean(vy)))



%{
figure('Units','inches',...
'PaperPositionMode','auto');
grid on, 
plotCellData(G,statesFP{n}.pressure),...
title('Pressure unsaturated domain');
set(gca,...
'Units','normalized',...
'FontUnits','points',...
'FontWeight','normal',...
'FontSize',19,...
'FontName','Times')
colorbar


figure('Units','inches',...
'PaperPositionMode','auto');
grid on, 
plotCellData(G,statesFP{n}.c),...
title('Concentration unsaturated domain');
set(gca,...
'Units','normalized',...
'FontUnits','points',...
'FontWeight','normal',...
'FontSize',19,...
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