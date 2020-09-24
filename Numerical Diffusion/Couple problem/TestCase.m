%% First test case
% No reaction term, easy theta and K expressions
%% IMPO: run mrst-bitbucket run code

% global state;
% initstate; state;

mrstModule add ad-blackoil ad-core ad-props mrst-gui blackoil-sequential
clear all
close all
clc
mrstVerbose off
gravity off
% gravity on
% gravity ([0,1,0])

X = 21;
Y = 31;
I = X;
J = Y;
a=0; b=2;
c=0; d=3;
x1=a; x2=b;
y1=c; y2=d;

%% Initial condition
dx = (b-a)/(I-1);
x = [a:dx:b];
x2 = (x(1:I-1)+x(2:I))/2;
dy=(d-c)/(J-1);
Lx=I-2; Ly=J-2;


G = cartGrid([X,Y], [2, 3]);
G = computeGeometry(G);

% Homogenous rock properties
rock = struct('perm', ones(G.cells.num, 1), ...
              'poro', ones(G.cells.num, 1));
          

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

icorr=0;
irand=0; %1=random Ksat; % irand=0; %0=const. Ksat %    L_p = .1;
if irand==1
    a=0; b=2;
    c=0; d=3;
    dx = (b-a)/(I-1);
    xx = [a:dx:b];
    x2 = (xx(1:I-1)+xx(2:I))/2;
    dy=(d-c)/(J-1);
    yy=c:dy:d;
    Nmod = 100;
    ZC1 = .1; 
    ZC2 = .01; 
    varK = 0.5;
    [XX,YY] = meshgrid(xx,yy);
    if icorr==1
        [wavenum, phi] = Kraichnan_Exp_param_short(Nmod,ZC1,ZC2);
    else
        [wavenum, phi] = Kraichnan_Gauss_param_short(Nmod,ZC1,ZC2);
    end
    C1 = wavenum(:,1);
    C2 = wavenum(:,2);
    Ks= K_r(XX(:)',YY(:)',Nmod,K_s,varK,C1,C2,phi);
    K_s=reshape(Ks,J,I);
    K_s = K_s';
    K_s = K_s(:);
else
    K_s = K_s*ones(J,I);
    K_s = K_s;
    K_s = K_s(:);
end

A=0.44;
B=0.0046;              

% Default oil-water fluid with unit values
fluid = initSimpleADIFluid('phases', 'W');
fluid.theta = @(p,c)   theta_s + 0.*p + 0.*c; %getThetaCoupled(p, c, theta_r, theta_s, alpha, ng, nGM,A,B);     %theta_s + 0.*p+ 0.*c;%     ./ (1 - p - 1/10.*c);   % theta is a function of pressure and concentration
fluid.Kmult = @(p,theta) K_s + 0.*theta + 0.*p; %getConductivity(p,theta, theta_r, theta_s, K_s, nGM); %K_s + 0.*p + 0.*theta;

% diffusion/dispersion coeff
D = 0.001;%0.001;  
D1 = D;
D2 = D;

L = .1; 

Richardsmodel = RichardsEquationModelSequential(G, rock, fluid,  'L_p', L);   % in the splitting solver we define the L_p and L_c separately
Transportmodel = SequentialTransportEquationModelnewBC(G, rock, fluid, 'L_c', 2*L, 'D', D);  % 
model = SequentialRichardsTransportModel(Richardsmodel, Transportmodel,'D', D);

y = G.cells.centroids(:,2);
x = G.cells.centroids(:,1);

p0 = 1 - y/3; % constant-parameters case!!!
c0 = (y-y1)/(y2-y1)/1.2;
% c0(size(c0,1)-X+1 : size(c0,1)-(X-1)/2) = 1;     % upper left side c = 1;
% c0(X:X:floor(X*Y/3)-1) = 0;  %% Change here the concentration on that side!

state0 = initResSol(G, p0);
state0.c = c0;

%{
figure
plotToolbar(G,state0);
colorbar;
%}


%% Time domain [0,1]
time = T;  % max time
dT = time/n; 

%% Boundary conditions, Dirichlet, pressure equal to zero on each side of the domain
nls = NonLinearSolver('maxIterations', 1000);
%% Boundary conditions, Dirichlet, pressure equal to zero on each side of the domain
% p_z = p0(1:X:floor(size(p0,1)/3)); %1 - z; 
% p_z = p_z(1:(Y-1)/3);
% 
% c_z = c0(1:X:floor(size(p0,1)/3));
% c_z = c_z(1:(Y-1)/3);
% 
% bfaces  = find(any(G.faces.neighbors==0,2));
%{
figure
hold on
plotFaces(G, bfaces(1:2:2*Y), 'Linewidth', 1); %left
plotFaces(G, bfaces(2:2:2*Y+1),'r', 'Linewidth', 1); %right
%plotFaces(G, bfaces(2:2:2*floor(Y/3)),'r', 'Linewidth', 2); %right lower
% plotFaces(G, bfaces((X+Y)*2-floor(X)+1 : (X+Y)*2-floor(X/2)-1),'r', 'Linewidth', 2);   %top left
%plotFaces(G, bfaces(2*X+1:3*X), 'y', 'Linewidth', 3); %down
% plotFaces(G, bfaces((X+Y)*2-floor(X)+1 : (X+Y)*2-floor(X/2)-1), 'y', 'Linewidth', 2); %top left
% plotFaces(G, bfaces(3*X+floor(X/2)+1:4*X), 'y', 'Linewidth',3); %top right
% plotFaces(G, bfaces(3*X:4*X-1),'g', 'Linewidth', 3); %down
plotFaces(G, bfaces(4*X:5*X-1),'g', 'Linewidth', 3); %down

grid on
title('Boundaries');
return
%}

%{
for j = 1:n
bc = [];
%
% no flow conditions left
%bc = addBC(bc, 'xmin', 'flux', 0, 'sat', [1]);

% no flow conditions down
%bc = addBC(bc, 'xmax', 'flux', 0, 'sat', [1]);

% Mixed conditions on up and right sides
%bc = addBC(bc, bfaces(3*X+floor(X/2)+1: 4*X), 'flux', 0, 'sat', [1]); % up no flux

bc = addBC(bc, bfaces(2:2:2*floor(Y/3)), 'pressure', p_z, 'sat', [1]); % pressure right
bc.c = 0.*ones(size(bc.sat,1), 1); 
%bc = addBC(bc, , 'flux', 0, 'sat', [1]); %right no flow

if j*dT <= DeltaT
    bc = addBC(bc, bfaces((X+Y)*2-floor(X)+1 : (X+Y)*2-floor(X/2)-1), 'pressure', (-2 + 2.2*(j*dT/DeltaT)), 'sat', [1]); % up pressure
    bc.c = [bc.c' 1 + 0.*bc.c']';
else
    bc = addBC(bc, bfaces((X+Y)*2-floor(X)+1 : (X+Y)*2-floor(X/2)-1), 'pressure', .2, 'sat', [1]); % up pressure 
    bc.c = [bc.c' 1 + 0.*bc.c']';
end

BC(j) = bc;
end

nls = NonLinearSolver('maxIterations', 1000);

schedule = simpleSchedule(repmat(dT,1,n),'bc', BC);
%}

%{
bc = [];
bc = addBC(bc, bfaces(2:2:2*floor(Y/3)), 'pressure', p_z, 'sat', [1]); % pressure right
bc.c = 0.*ones(size(bc.sat,1), 1); 
bc = pside(bc,G,'YMax', 0, 'sat', [1]);
bc.c = horzcat(bc.c',ones(size(bc.sat,1)-size(bc.c,1), 1)')';
bc = pside(bc,G,'YMin', 1, 'sat', [1]);
bc.c = horzcat(bc.c', 0*ones(size(bc.face,1)-size(bc.c,1), 1)')';
%}
%
bc = [];
bc = pside(bc,G,'YMax', 0 , 'sat', [1]);
bc.c = ones(size(bc.sat,1),1);
bc = pside(bc,G,'YMin', 1, 'sat', [1]);
bc.c = [bc.c' 0.*bc.c']';
%}

%{
for j = 1:n
bc = [];
%
% no flow conditions left
%bc = addBC(bc, 'xmin', 'flux', 0, 'sat', [1]);

% no flow conditions down
%bc = addBC(bc, 'xmax', 'flux', 0, 'sat', [1]);

% Mixed conditions on up and right sides
%bc = addBC(bc, bfaces(3*X+floor(X/2)+1: 4*X), 'flux', 0, 'sat', [1]); % up no flux

bc = addBC(bc, bfaces(2:2:2*floor(Y/3)), 'pressure', p_z, 'sat', [1]); % pressure right
bc.c = c_z.*ones(size(bc.sat,1), 1); 
%bc = addBC(bc, , 'flux', 0, 'sat', [1]); %right no flow

if j*dT <= DeltaT
    bc = addBC(bc, bfaces((X+Y)*2-floor(X)+1 : (X+Y)*2-floor(X/2)-1), 'pressure', (-2 + 2.2*(j*dT/DeltaT)), 'sat', [1]); % up pressure
    bc.c = [bc.c' 1 + 0.*bc.c']';
else
    bc = addBC(bc, bfaces((X+Y)*2-floor(X)+1 : (X+Y)*2-floor(X/2)-1), 'pressure', .2, 'sat', [1]); % up pressure 
    bc.c = [bc.c' 1 + 0.*bc.c']';
end

BC(j) = bc;
end
%}
schedule = simpleSchedule(repmat(dT,1,n),'bc', bc);

tic
[~,states, report] = simulateScheduleAD(state0, model, schedule, 'nonlinearsolver', nls);
toc

% figure
% mesh(xx(1:I),yy(1:J),c)
% title('Final numerical c') 
% print -depsc2 FinalNumericalC.eps
% 
% err_c = norm(states{n}.c-C_gss)./norm(C_gss);
% message = ['Relative c error: ', num2str(err_c)];
% disp(message)
% 
% figure
% mesh(xx(1:I),yy(1:J),c-c_gss')
% title('Numerical error c') 
% print -depsc2 NumericalErrorC.eps
% 
% 
% figure
% mesh(xx(1:I),yy(1:J),c_gss')
% title('Final analythical c')
% print -depsc2 FinalAnalythicalC.eps
%%
p = zeros(J,I);
v = faceFlux2cellVelocity(G,states{n}.flux);
% v = cellFlux2faceFlux(G, sum(states{n}.flux, 2));
theta = p;
c=p;

for l = 1:J
    p_mrst(l,:) = states{n}.pressure(((l-1)*I+1:l*I));
    theta(l,:) = states{n}.theta(((l-1)*I+1:l*I));
    c_mrst(l,:) = states{n}.c(((l-1)*I+1:l*I));
    vx(l,:) = v(((l-1)*I+1:l*I),1);
    vy(l,:) = v(((l-1)*I+1:l*I),2);
end

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
%}

% figure
% plotToolbar(G,states)
% title('Final states')
% colorbar

%% Save the results
vx_clay = vx;
vy_clay = vy;
c_clay = c_mrst;
p_clay = p_mrst;
theta_clay = theta;
mean_vx_clay = mean(mean(vx_clay));
mean_vy_clay = mean(mean(vy_clay));
% save('Results_clay_21_31.mat','vx_clay', 'vy_clay', 'c_clay', 'p_clay', 'theta_clay','mean_vx_clay', 'mean_vy_clay')

