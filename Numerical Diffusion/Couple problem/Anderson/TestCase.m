%% First test case
% No reaction term, easy theta and K expressions
%% IMPO: run mrst-bitbucket run code

%global state;
%initstate; state;

mrstModule add ad-blackoil ad-core ad-props mrst-gui blackoil-sequential
clear all
close all
clc
mrstVerbose off
gravity on

X = 21;
Y = 31;
I = X;
J = Y;
a=0; b=2;
c=0; d=3;

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
          
T = 3;
DeltaT = T/3;

itest=1; % 1: U=-0.008; 2: U=-0.08; 3: U=0.8 
if itest==1
    K_s = 4.96*10^-2; % by mofifying K_sat ===> diferent Pecelt numbers
    %U_MEAN = - 0.0331;    %-0.033 x=401 y=601
    n = 10;
    theta_s = .396;
    theta_r = .131;
    alpha = .423;
    nGM = 2.06;
    ng = (nGM-1)/nGM;
end

theta_s = .396;
% Default oil-water fluid with unit values
fluid = initSimpleADIFluid('phases', 'W');
fluid.theta = @(p,c)   theta_s + 0.*p+ 0.*c; %getThetaCoupled(p, c, theta_r, theta_s, alpha, ng, nGM);     %theta_s + 0.*p+ 0.*c;%     ./ (1 - p - 1/10.*c);   % theta is a function of pressure and concentration
fluid.Kmult = @(p,theta) getConductivity(p,theta, theta_r, theta_s, K_s, nGM); %K_s + 0.*p + 0.*theta;

% diffusion/dispersion coeff
D = 0.001;%0.001;  
D1 = D;
D2 = D;

Richardsmodel = RichardsEquationModelSequential(G, rock, fluid,  'L_p', 1);   % in the splitting solver we define the L_p and L_c separately
Transportmodel = SequentialTransportEquationModelnewBC(G, rock, fluid, 'L_c', 1, 'D', D);  % 
model = SequentialRichardsTransportModel(Richardsmodel, Transportmodel,'D', D);

y = G.cells.centroids(:,2);
p0 = 1 - y;  
% how is 'y' defined? is 'p0' the same initial pressure as in the GRW codes? !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c0 = y/3/1.2;
c0(size(c0,1)-X+1 : size(c0,1)-(X-1)/2) = 1;     % upper left side c = 1;
c0(X:X:X*Y/3-1) = 0;  %% Change here the concentration on that side!

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
nls = NonLinearSolver('maxIterations', 10000);
%% Boundary conditions, Dirichlet, pressure equal to zero on each side of the domain
p_z = p0(1:X:floor(size(p0,1)/3)); %1 - z; 
p_z = p_z(1:(Y-1)/3);

bfaces  = find(any(G.faces.neighbors==0,2));

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
p = zeros(J,I);
v = faceFlux2cellVelocity(G,states{n}.flux);
theta = p;
c=p;

for l = 1:J
    p_mrst(l,:) = states{n}.pressure(((l-1)*I+1:l*I));
    theta(l,:) = p_mrst(l,:).*0 + theta_s; %states{n}.theta(((l-1)*I+1:l*I));
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

%% Save the results
vx_loam = vx;
vy_loam = vy;
c_loam = c_mrst;
p_loam = p_mrst;
theta_loam = theta;
mean_vx_loam = mean(mean(vx_loam));
mean_vy_loam = mean(mean(vy_loam));
%  save('Results_loam_41_61.mat','vx_loam', 'vy_loam', 'c_loam', 'p_loam', 'theta_loam','mean_vx_loam', 'mean_vy_loam')

