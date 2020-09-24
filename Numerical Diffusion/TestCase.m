%% numerical diffusion test case
% No reaction term, easy theta and K expressions
%% IMPO: run mrst-bitbucket run code

%global state;
%initstate; state;

mrstModule add ad-blackoil ad-core ad-props mrst-gui blackoil-sequential
clear all
close all
clc
mrstVerbose off
gravity off

gravity on
gravity([0,-1,0])

%%   Grid Initialization
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

itest=1; % 1: U=-0.008; 2: U=-0.08; 3: U=0.8 
if itest==1
    y0=(Ly+1)*0.7*dy; x0=(Lx+1)/2*dx;
    K_s = 4.96*10^-2; % by mofifying K_sat ===> diferent Pecelt numbers
    U_MEAN = - 0.0331;    %-0.033 x=401 y=601
    n = 10;
elseif itest==2
    y0=(Ly+1)*0.7*dy; x0=(Lx+1)/2*dx;
    K_s = 0.12; 
    U_MEAN = - 0.08;
    n = 50;
else
    y0=(Ly+1)*0.99*dy; x0=(Lx+1)/2*dx;
    K_s = 1.2;
    U_MEAN = - 0.8;
    n = 15;
end
theta_s = 1;
          
% Default oil-water fluid with unit values
fluid = initSimpleADIFluid('phases', 'W');
fluid.theta = @(p,c) theta_s + 0.*p + 0.*c;%getThetaCoupled(p, c, theta_r, theta_s, alpha, ng, nGM, A, B);     %1 ./ (1 - p - 1/10.*c);   % theta is a function of pressure and concentration
fluid.Kmult = @(p,theta)  K_s + 0.*p + 0.*theta;%getConductivity(p,theta, theta_r, theta_s, K_s, nGM); 

% diffusion/dispersion coeff
D = 0.001;%0.001;  
D1 = D;
D2 = D;

Richardsmodel = RichardsEquationModelSequential(G, rock, fluid,  'L_p', .01);   % in the splitting solver we define the L_p and L_c separately
Transportmodel = SequentialTransportEquationModelnewBC(G, rock, fluid, 'L_c', .001, 'D', D);  % 
model = SequentialRichardsTransportModel(Richardsmodel, Transportmodel,'D', D);

xx = G.cells.centroids(:,1);
p0 = 1 - ((xx-a)/(b-a));

% how is 'y' defined? is 'p0' the same initial pressure as in the GRW codes? !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ti=1;
gss=Gauss_IC(ti,dx,dy,x0,y0,Lx,Ly,U_MEAN,D1,D2);
c0 = gss;%/sum(sum(gss)); % thin Gaussian distribution IC for MRST !!!

%% initial moments
x1=0; x2=2;
y1=0; y2=3;
xx=x1:dx:x2; yy=y1:dy:y2;
X=xx(2:I-1); Y=yy(2:J-1);

xr=(X-x0); yr=(Y-y0);

figure
mesh(xx(2:I-1),yy(2:J-1),c0')
mesh(X,Y,c0')
title('Initial c')

%% Extend initial concentration
c0 = [c0; zeros((size(c0,2)),1)'];
c0 = [zeros((size(c0,1)),1) c0];
c0 = [c0 zeros((size(c0,1)),1)];
c0 = [zeros((size(c0,2)),1)';  c0];

Mx0=0; Varx0=0;
My0=0; Vary0=0; 
x = G.cells.centroids(:,1);
y = G.cells.centroids(:,2);

    x = x(1:I);
    y = y(1:I:size(y,1));
    
        for XX = 1 : size(x,1)
               for YY = 1 : size(y,1)
                    Mx0 = Mx0 + x(XX) .* c0(XX,YY);%./(size(x,1)*size(y,1)) ;
                    Varx0 = Varx0 +  x(XX).^2 .* c0(XX,YY);
                    My0 = My0 + y(YY) .* c0(XX,YY); 
                    Vary0 = Vary0 + y(YY).^2 .*c0(XX,YY); 
               end
        end
       
 S_xx0 = (Varx0 - Mx0.^2);   % they should be the same because it is symmetrical
 S_yy0 = (Vary0 - My0.^2);
    
    
c0 = c0(:); 

state0 = initResSol(G, p0);
state0.c = c0;


%% Time domain [0,1]
time = T;  % max time
dT = time/n; 

%% Boundary conditions, Dirichlet, pressure equal to zero on each side of the domain
nls = NonLinearSolver('maxIterations', 10000);

bc = [];
bc = pside(bc,G,'YMax', 0, 'sat', [1]);
bc = pside(bc,G,'YMin', 1, 'sat', [1]);
bc.c = 0.*ones(size(bc.sat,1), 1);

schedule = simpleSchedule(repmat(dT,1,n),'bc', bc);

tic
[~,states, report] = simulateScheduleAD(state0, model, schedule, 'nonlinearsolver', nls);
toc

tt=T+ti;
gss=Gauss_IC(tt,dx,dy,x0,y0,Lx,Ly,U_MEAN,D1,D2);
c_gss=gss/sum(sum(gss)); %c=c_gss;
c_gss = [c_gss; zeros((size(c_gss,2)),1)'];
c_gss = [zeros((size(c_gss,1)),1) c_gss];
c_gss = [c_gss zeros((size(c_gss,1)),1)];
c_gss = [zeros((size(c_gss,2)),1)';  c_gss];

C_gss = c_gss(:); 

c = zeros(J,I);
    for l = 1:J
        c(l,:) = states{n}.c(((l-1)*I+1:l*I));
    end
    
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
states{n}.flux = -1.*abs(states{n}.flux);
v = faceFlux2cellVelocity(G,states{n}.flux);
%%
for l = 1:J
    vx(l,:) = v(((l-1)*I+1:l*I),1);
    vy(l,:) = v(((l-1)*I+1:l*I),2);
end

%% compute now numerical diffusion/dispersion coeff by 2011 List&al articel
% m_x(t) = int(int(x*c(t,x,y),dx),dy)
% m_y(t) = int(int(y*c(t,x,y),dx),dy)

c = zeros(J,I);
TR = n;
Mx=zeros(1,TR); Varx=zeros(1,TR);
My=zeros(1,TR); Vary=zeros(1,TR);
figure
hold on
for N = 1:n

    for l = 1:J
        c(l,:) = states{N}.c(((l-1)*I+1:l*I));
    end
    
    C = c'; 
    
x = G.cells.centroids(:,1);
y = G.cells.centroids(:,2);

    x = x(1:I);
    y = y(1:I:size(y,1));
    
        for XX = 1 : size(x,1)
               for YY = 1 : size(y,1)
                   
                    Mx(N) = Mx(N) + x(XX) .* C(XX,YY);
                    Varx(N) = Varx(N) + x(XX).^2 .* C(XX,YY); % for MRST !!!,
                    My(N) = My(N) + y(YY) .*C(XX,YY); 
                    Vary(N) = Vary(N) + y(YY).^2 .*C(XX,YY); % where c(x,y)=concentration                    
              
               end
    end
    %}
    
N = N + 1;
end

xx=x1:dx:x2; yy=y1:dy:y2;
X=xx(1:I); Y=yy(1:J);
figure
mesh(xx(1:I),yy(1:J),c)
% in the GRW function 'testBGRW_const.m' it was mesh(X,Y,c'), where X=xx(2:I-1); Y=yy(2:J-1); !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
title('Final c')

ntot=1;
dt = T/TR;

for t=1:TR
    %{
% in GRW 'ntot' is used to define the concentration as normalized nr. of particles; in MRST you have the concentrations !!!!!!!!!!!!!!!    
    Varx(t)=Varx(t)/ntot-Mx(t)*Mx(t)-Varx0; 
    Vary(t)=Vary(t)/ntot-My(t)*My(t)-Vary0;  % add initlal variance
    
    Dx(t)=Varx(t)/(2*t*dt); 
    Dy(t)=Vary(t)/(2*t*dt);
    
    Mx(t) = (Mx(t) - Mx0)/(t*dt); 
    My(t) = (My(t) - My0)/(t*dt);
    %}
    
    S_xx(t) = Varx(t) - Mx(t).^2;
    S_yy(t) = Vary(t) - My(t).^2;
    
    Dx(t) = (S_xx(t) - S_xx0)/(2*(t*dt));
    Dy(t) = (S_yy(t) - S_yy0)/(2*(t*dt));
    
    %Dx(t)=Varx(t)/(2*t*dt); 
    %Dy(t)=Vary(t)/(2*t*dt);
    
    Mx(t) = (Mx(t) - Mx0)/(t*dt); 
    My(t) = (My(t) - My0)/(t*dt);
end
%
%%
%
t=1:TR;

figure;
subplot(1,2,1)
plot(t*dt,Mx,t*dt,My); 
legend('Mx(t)','My(t)','Location','best'); legend('boxoff');
xlabel('t'); 
subplot(1,2,2)
plot(t*dt,Dx,t*dt,Dy); 
legend('Dx(t)','Dy(t)','Location','best'); legend('boxoff');
xlabel('t'); 
%print -depsc2 MxMy_DxDy_plots.eps

figure; hold all
subplot(1,2,1)
semilogy(t*dt,norm(Mx-0)/norm(mean(mean(vx))),'.',t*dt,norm(My-U_MEAN)/norm(U_MEAN),'.');
legend('|V_{x}(t) - V_{0x}| /|V_{0x}|','|V_{y}(t) - V_{0y}| /|V_{0y}|','Location','best');
subplot(1,2,2)
semilogy(t*dt,norm(Dx-D1)/D1,'.',t*dt,norm(Dy-D2)/D2,'.');
legend('|D_{x}(t) - D_{0x}| / D_{0x}','|D_{y}(t) - D_{0y}| / D_{0y}','Location','best');
%print -depsc2 VxVy_D0xD0y_plots.eps
%}


%{
t=1:TR;
figure(5);
subplot(1,2,1)
plot(t*dt,Mx,t*dt,My); 
legend('Mx(t)','My(t)','Location','best'); legend('boxoff');
xlabel('t'); 

subplot(1,2,2)
plot(t*dt,Dx,t*dt,Dy); 
legend('Dx(t)','Dy(t)','Location','best'); legend('boxoff');
xlabel('t'); 
print -depsc2 MxMy_DxDy_plots.eps

figure(6); hold all
subplot(1,2,1)
semilogy(t*dt,abs(Mx-0)/abs(U_MEAN),'.',t*dt,abs(My-U_MEAN)/abs(U_MEAN),'.');
legend('|V_{x}(t) - V_{0x}| /|V_{0x}|','|V_{y}(t) - V_{0y}| /|V_{0y}|','Location','best');
subplot(1,2,2)
semilogy(t*dt,abs(Dx-D1)/D1,'.',t*dt,abs(Dy-D2)/D2,'.');
legend('|D_{x}(t) - D_{0x}| / D_{0x}','|D_{y}(t) - D_{0y}| / D_{0y}','Location','best');
print -depsc2 VxVy_D0xD0y_plots.eps
%}
err_c = norm(states{n}.c-C_gss)./norm(C_gss)
eps_D1=sqrt(dt)*norm(Dx-D1)/D1/T
eps_D2=sqrt(dt)*norm(Dy-D2)/D2/T
time_steps = TR
U_mean = mean(mean(vy))
U_MEAN = -0.033;
E_mean_error = (U_mean - U_MEAN)



%% Compare concentration
%{
load('c_nicu');
c_Nicu = c_nicu;

c_Nicu = [c_Nicu; zeros((size(c_Nicu,2)),1)'];
c_Nicu = [zeros((size(c_Nicu,1)),1) c_Nicu];
c_Nicu = [c_Nicu zeros((size(c_Nicu,1)),1)];
c_Nicu = [zeros((size(c_Nicu,2)),1)';  c_Nicu];

c_Nicu = c_Nicu(:);
c_mrst = states{n}.c;

disp('Error: ')
norm(c_Nicu-c_mrst)/norm(c_mrst)
%}
%{
figure
plotToolbar(G,states)
colorbar
title('Numerical solution')
%}
%{
p = vec2mat(states{n}.pressure,X);
figure
contourf(p,12)
colormap(flipud(parula)); colorbar;
set(gca,...
'Units','normalized',...
'FontUnits','points',...
'FontWeight','normal',...
'FontSize',17,...
'FontName','Times')
xlabel('$x$','Interpreter','latex'); ylabel('$z$','Interpreter','latex'); 
title('$\Psi(x,z,t)$  loam, rand=1/Gauss','Interpreter','latex'); 
print -depsc2 Psi_loam_gauss.eps

c = vec2mat(states{n}.c,X);
figure
contourf(c,12)
colormap(flipud(parula)); colorbar;
set(gca,...
'Units','normalized',...
'FontUnits','points',...
'FontWeight','normal',...
'FontSize',17,...
'FontName','Times')
xlabel('$x$','Interpreter','latex'); ylabel('$z$','Interpreter','latex'); 
title('$c(x,z,t)$ loam, rand=1/Gauss','Interpreter','latex'); 
print -depsc2 c_loam_gauss.eps

theta = vec2mat(states{n}.theta,X);
figure
contourf(theta,12)
colormap(flipud(parula)); colorbar;
set(gca,...
'Units','normalized',...
'FontUnits','points',...
'FontWeight','normal',...
'FontSize',17,...
'FontName','Times')
xlabel('$x$','Interpreter','latex'); ylabel('$z$','Interpreter','latex'); 
title('$\theta(\Psi,c)$ loam, rand=1/Gauss','Interpreter','latex'); 
print -depsc2 theta_loam_gauss.eps


v = faceFlux2cellVelocity(G,states{n}.flux);

for l = 1:J
    vx(l,:) = v(((l-1)*I+1:l*I),1);
    vy(l,:) = v(((l-1)*I+1:l*I),2);
end

vx = vec2mat(vx,X);
figure
contourf(vx,12)
colormap(flipud(parula)); colorbar;
set(gca,...
'Units','normalized',...
'FontUnits','points',...
'FontWeight','normal',...
'FontSize',17,...
'FontName','Times')
xlabel('$x$','Interpreter','latex'); ylabel('$z$','Interpreter','latex'); 
title('$vx(x,z,t)$ loam, rand=1/Gauss','Interpreter','latex'); 
print -depsc2 vx_loam_gauss.eps

vy = vec2mat(vy,X);
figure
contourf(vy,12)
colormap(flipud(parula)); colorbar;
set(gca,...
'Units','normalized',...
'FontUnits','points',...
'FontWeight','normal',...
'FontSize',17,...
'FontName','Times')
xlabel('$x$','Interpreter','latex'); ylabel('$z$','Interpreter','latex'); 
title('$vy(x,z,t)$ loam, rand=1/Gauss','Interpreter','latex'); 
print -depsc2 vy_loam_gauss.eps

%}
%% Plot rate of convergence of linearization scheme plotting residual 
% for pressure and concentration at the final time step
%{

for i = 1:size(report.ControlstepReports{n}.StepReports{1}.NonlinearReport,1)
    res_p(i) = report.ControlstepReports{n}.StepReports{1}.NonlinearReport{i}.PressureSolver.StepReports{1}.NonlinearReport{1}.Residuals(1);
end

for i = 1:size(report.ControlstepReports{n}.StepReports{1}.NonlinearReport,1)
    res_c(i) = report.ControlstepReports{n}.StepReports{1}.NonlinearReport{i}.TransportSolver.StepReports{1}.NonlinearReport{1}.Residuals(1);
end


%figure
hold on
plot(1:size(res_p,2),log10(res_p), 'LineWidth', 1.5);
plot(1:size(res_c,2),log10(res_c), 'LineWidth', 1.5);
legend('Pressure residuals', 'Concentration residuals')

ylabel('log(Residuals)')
xlabel('Iterations')

title('Rate fo convergence')
%}