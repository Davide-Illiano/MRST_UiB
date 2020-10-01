%% numerical diffusion test case
% No reaction term, easy theta and K expressions
%% IMPO: run mrst-bitbucket run code
% maxNumCompThreads(1)
%global state; 
%initstate; state;

mrstModule add ad-blackoil ad-core ad-props mrst-gui blackoil-sequential
clear all
close all
clc
mrstVerbose off
gravity off

tic
for test = [1 2]

    tic
    clear Vx Vy Mx My Dx Dy p_mrst
    numbRealiz = 100;
    Nmod = 10; %10^2 ;
    varK= 0.1 ;
    ZC1 = 1.0;
    ZC2 = 1.0;
    KMean = 15;
    
    I= test * 130 + 1; %801; %401; %
    J= test * 70 + 1;  %401; %201; %
    a=0; b=13;  %original 20 and 10 we use half
    c=0; d=7;
    dx=(b-a)/(I-1);
    x=a:dx:b;
    x2=(x(1:I-1)+x(2:I))/2;
    dy=(d-c)/(J-1);
    y=c:dy:d;
    y2=(y(1:J-1)+y(2:J))/2;
    
% diffusion/dispersion coeff
    D = 0.01; %0.01;
    D1 = D;
    D2 = D;

    U_mean = 0.7134;
    Pe = U_mean * dx/D; 
    
    Lx=I-2; Ly=J-2;
    x0 = round(Lx*dx/8);  %x0=round(Lx*dx/10);
    y0 = round(Ly*dy/2);
    
    G = cartGrid([I,J], [b, d]);  % original 20, 10 we use hald
    G = computeGeometry(G);
    
    % Homogenous rock properties
    rock = struct('perm', ones(G.cells.num, 1), ...
        'poro', ones(G.cells.num, 1));
    
    T = 10;   % original 10 we use half
    
    itest=1; % 1: U=-0.008; 2: U=-0.08; 3: U=0.8
    if itest==1
        K_s = 15; % by mofifying K_sat ===> diferent Pecelt numbers
        U_MEAN = 0.7134;    %-0.033 x=401 y=601
        n = 20 * test;
    end
    
    theta_s = 1;
    
    % Default oil-water fluid with unit values
    fluid = initSimpleADIFluid('phases', 'W');
    fluid.theta = @(p,c) theta_s + 0.*p + 0.*c;%getThetaCoupled(p, c, theta_r, theta_s, alpha, ng, nGM, A, B);     %1 ./ (1 - p - 1/10.*c);   % theta is a function of pressure and concentration
    fluid.Kmult = @(p,theta)  1 + 0.*p + 0.*theta;%getConductivity(p,theta, theta_r, theta_s, K_s, nGM);
    

    
    xx = G.cells.centroids(:,1);
    p0 = 1 - ((xx-a)/(b-a));
    
    ti = 0.01;
    gss=Gauss_IC(ti,dx,dy,x0,y0,Lx,Ly,U_MEAN,D);
    c0 = gss/sum(sum(gss)); %/sum(sum(gss));
    
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
    
    
    
    state0 = initResSol(G, p0);
    state0.c = c0(:);
    
    %{
figure
plotToolbar(G,state0);
colorbar;
    %}
    TR = n;
    Mx=zeros(1,TR); Varx=zeros(1,TR);
    My=zeros(1,TR); Vary=zeros(1,TR);
    
    for numb = 1:numbRealiz
        %% Time domain [0,1]
        time = T;  % max time
        dT = time/n;
        
        %% Boundary conditions, Dirichlet, pressure equal to zero on each side of the domain
        nls = NonLinearSolver('maxIterations', 10000);
        
        bc = [];
        bc = pside(bc,G,'XMin', 0.6183 , 'sat', [1]);
        bc = pside(bc,G,'XMax', 0 , 'sat', [1]);
        bc.c = 0.*ones(size(bc.sat,1), 1);
        
        schedule = simpleSchedule(repmat(dT,1,n),'bc', bc);
        
        [wavenum, phi] = Kraichnan_Gauss_param_short(Nmod,ZC1,ZC2);
        
        %% Assemble the conductivity
        C1 = wavenum(:,1);
        C2 = wavenum(:,2);
        [X,Y] = meshgrid(x(2:I-1),y(2:J-1));
        d = K_repmat(X(:)',Y(:)',Nmod,KMean,varK,C1,C2,phi);
        d = reshape(d,I-2,J-2); % it was J-2 I-2
        
        % Need to extend size of d to J and I
        d = [d; zeros((size(d,2)),1)'];
        d = [zeros((size(d,1)),1) d];
        d = [d zeros((size(d,1)),1)];
        d = [zeros((size(d,2)),1)';  d];
        
        % set the values
        d(1,:) = d(2,:);
        d(size(d,1),:) = d(size(d,1)-1,:);
        d(:,1) = d(:,2);
        d(:,size(d,2)) = d(:,size(d,2)-1);
        %d = 0.* d + 15;
        
        model.K = d;
        
        %% Fixed point solver
        model = RichardsTransportEquationFixedPointSchemes(G, rock, fluid, 'Newton', 1, 'D', D , 'K', d);    %'Mixed', 1, 'L_p', .2, 'L_c', .01
        model.diffusion = 1;
        
        
        [~,states, report] = simulateScheduleAD(state0, model, schedule, 'nonlinearsolver', nls);
        
               
        v = faceFlux2cellVelocity(G,states{n}.flux);
        %%
        for l = 1:J
            Vx(l,:) = v(((l-1)*I+1:l*I),1);
            Vy(l,:) = v(((l-1)*I+1:l*I),2);
            p_mrst(l,:) = states{n}.pressure(((l-1)*I+1:l*I));
        end
        
        c = p_mrst;
        %     save(['varK_GAUSS_10_MRST_Mesh(',num2str(I),',',num2str(J),')_'num2str(numb),'.mat'],'c','Vx','Vy') ;
        
      
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
                    %                   xr=(x*dx-x0); yr=(y*dy-y0); % from 'testBGRW_const.m'
                    %xr=(xxx*dx-x0); yr=(yyy*dy-y0); % modified by me !!!!!!!!!!!!!!!!!!!
                    
                    % change c(yy,xx) by c(xx,yy) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    % if MRST gives c(nr. lines, nr. columns) please transpose c in the following !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    Mx(N) = Mx(N) + x(XX) .* C(XX,YY);
                    Varx(N) = Varx(N) + x(XX).^2 .* C(XX,YY); % for MRST !!!,
                    My(N) = My(N) + y(YY) .*C(XX,YY);
                    Vary(N) = Vary(N) + y(YY).^2 .*C(XX,YY); % where c(x,y)=concentration
                    
                end
            end
            %}
            
            N = N + 1;
        end
        numb
    end
    %%
    %     save(['varK_GAUSS_10_MRST_Mesh(',num2str(I),',',num2str(J),')_',num2str(numb),'.mat'],'c','Vx','Vy') ;
    
    
    dt = T/TR;
    
    for t=1:TR
        
        Mx(t) = Mx(t)/numbRealiz;
        My(t) = My(t)/numbRealiz;
        Varx(t) = Varx(t)/numbRealiz;
        Vary(t) = Vary(t)/numbRealiz;
        
        S_xx(t) = Varx(t) - Mx(t).^2;
        S_yy(t) = Vary(t) - My(t).^2;
        
        Dx(t) = (S_xx(t) - S_xx0)/(2*(t*dt));
        Dy(t) = (S_yy(t) - S_yy0)/(2*(t*dt));
        
        Mx(t) = (Mx(t) - Mx0)/(t*dt);
        My(t) = (My(t) - My0)/(t*dt);
        %}
        %      for t=1:n
        %         Mx(t)=Mx(t)/numbRealiz; My(t)=My(t)/numbRealiz;
        %         Varx(t)=Varx(t)/numbRealiz-Mx(t)^2; Vary(t)=Vary(t)/numbRealiz-My(t)^2;
        %         Dx(t)=Varx(t)/(2.0*t*dt); Dy(t)=Vary(t)/(2.0*t*dt);
        %         Mx(t)= Mx(t)/(t*dt); My(t) = My(t)/(t*dt);
        %     end
        %
        %      for t=1:n
        %         Mx(t)=Mx(t)/numbRealiz; My(t)=My(t)/numbRealiz;
        %         Varx(t)=Varx(t)/numbRealiz-Mx(t)*Mx(t)-Varx0; Vary(t)=Vary(t)/numbRealiz-My(t)*My(t)-Vary0;
        %
        %         S_xx(t) = Varx(t) - Mx(t).^2;
        %         S_yy(t) = Vary(t) - My(t).^2;
        %
        %         Dx(t) = (S_xx(t) - S_xx0)/(2*(t*dt));
        %         Dy(t) = (S_yy(t) - S_yy0)/(2*(t*dt));
        %
        %
        %         %Dx(t)=Varx(t)/(2.0*t*dt); Dy(t)=Vary(t)/(2.0*t*dt);
        %         Mx(t)=(Mx(t)-Mx0)/(t*dt); My(t)=(My(t)-My0)/(t*dt);
        %     end
        
    end
 %save(['ConstantK_MXEtc_MRST_Mesh(',num2str(I),',',num2str(J),')_n(,',num2str(n),').mat'], 'Mx', 'Dx', 'My', 'Dy');

toc
        
    %%
    %{
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
    %}
    % figure; hold all
    % subplot(1,2,1)
    % semilogy(t*dt,norm(Mx-0)/norm(mean(mean(vx))),'.',t*dt,norm(My-U_MEAN)/norm(U_MEAN),'.');
    % legend('|V_{x}(t) - V_{0x}| /|V_{0x}|','|V_{y}(t) - V_{0y}| /|V_{0y}|','Location','best');
    % subplot(1,2,2)
    % semilogy(t*dt,norm(Dx-D1)/D1,'.',t*dt,norm(Dy-D2)/D2,'.');
    % legend('|D_{x}(t) - D_{0x}| / D_{0x}','|D_{y}(t) - D_{0y}| / D_{0y}','Location','best');
    %print -depsc2 VxVy_D0xD0y_plots.eps
    %}
    
    
    %{
t=1:TR;
figure(5);
subplot(1,2,1)
plot(t*dt,Mx,t*dt,My);
legend('Mx(t)','My(t)','Location','best'); legend('boxoff');
xlabel('t');
%
subplot(1,2,2)
%
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
    
    eps_D1=sqrt(dt)*norm(Dx-D1)/D1/T
    eps_D2=sqrt(dt)*norm(Dy-D2)/D2/T
    
    dx
    dy
    Pe
    
    time_steps = TR;

    %Rand = randi(100000,1)
%save(['ConstK_MXEtc_MRST_Mesh(',num2str(I),',',num2str(J),')_n(,',num2str(n),').mat'], 'dx', 'Pe', 'Mx', 'Dx', 'My', 'Dy', 'eps_D1', 'eps_D2');


     Rand = randi(1000000,1)
save(['NewProb_RandK_MRST_Mesh(',num2str(dx),',',num2str(dy),')_n(',num2str(n),')',num2str(Rand),'.mat'], 'dx', 'Pe', 'Mx', 'Dx', 'My', 'Dy', 'eps_D1', 'eps_D2');

% save(['Test20','/tentative(',num2str(I),',',num2str(J),')_n(,',num2str(n),').mat'], 'dx', 'Pe', 'Mx', 'Dx', 'My', 'Dy', 'eps_D1', 'eps_D2');


end


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
=======
     Rand = randi(100000,1)
save(['NewProb_RandK_MRST_Mesh(',num2str(I),',',num2str(J),')_n(',num2str(n),')',num2str(Rand),'.mat'], 'dx', 'Pe', 'Mx', 'Dx', 'My', 'Dy', 'eps_D1', 'eps_D2');
>>>>>>> f77ee09fec118f0064c76c76946f698ee856a2f7

% save(['NewProb_constantK_MRST_Mesh(',num2str(I),',',num2str(J),')_n(',num2str(n),').mat'], 'dx', 'Pe', 'Mx', 'Dx', 'My', 'Dy', 'eps_D1', 'eps_D2');


end


<<<<<<< HEAD
title('Rate fo convergence')
%}
toc
