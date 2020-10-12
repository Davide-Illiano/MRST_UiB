%% numerical diffusion test case
% No reaction term, easy theta and K expressions
%% IMPO: run mrst-bitbucket run code
maxNumCompThreads(1)
%global state; 
%initstate; state;

mrstModule add ad-blackoil ad-core ad-props mrst-gui blackoil-sequential
clear all
close all
clc
mrstVerbose off
gravity off

tic
for test = [10]
     clear Vx Vy Mx My Dx Dy p_mrst
    numbRealiz = 10;
    Nmod = 10; %10^2 ;
    varK= 0.1 ;
    ZC1 = 0.1;
    ZC2 = 0.1;
    KMean = 15;
    
    I= test * 200 + 1; %801; %401; %
    J= test * 100 + 1;  %401; %201; %
    a=0; b=2;  %original 20 and 10 we use half
    c=0; d=1;
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
    x0 = 0.2;%round(Lx*dx/8);  %x0=round(Lx*dx/10);
    y0 = 0.5;%round(Ly*dy/2);
    
    G = cartGrid([I,J], [b, d]);  % original 20, 10 we use hald
    G = computeGeometry(G);
    
    % Homogenous rock properties
    rock = struct('perm', ones(G.cells.num, 1), ...
        'poro', ones(G.cells.num, 1));
    
    T = 0.5;   % original 10 we use half
    
    itest=1; % 1: U=-0.008; 2: U=-0.08; 3: U=0.8
    if itest==1
        K_s = 15; % by mofifying K_sat ===> diferent Pecelt numbers
        U_MEAN = 0.7134;    %-0.033 x=401 y=601
        n = test * 20;
        n_flow = 1;
    end
    
    theta_s = 1;
    
    % Default oil-water fluid with unit values
    fluid = initSimpleADIFluid('phases', 'W');
    fluid.theta = @(p) theta_s + 0.*p + 0.*c;%getThetaCoupled(p, c, theta_r, theta_s, alpha, ng, nGM, A, B);     %1 ./ (1 - p - 1/10.*c);   % theta is a function of pressure and concentration    

    
    xx = G.cells.centroids(:,1);
    p0 = 0.1 - 0.1*((xx-a)/(b-a));  % not really need
    
    ti = 0.0001;
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
    

    TR = n;
    Mx=zeros(1,TR); Varx=zeros(1,TR);
    My=zeros(1,TR); Vary=zeros(1,TR);
    
    
    for numb = 1:numbRealiz
        tic
        %% Time domain [0,1]
        time = T;  % max time
        dT = time/n;
        dT_flow = time/n_flow;
        
        %% Boundary conditions, Dirichlet, pressure equal to zero on each side of the domain
        nls = NonLinearSolver('maxIterations', 10000);
        
        bc_flow = [];
        bc_flow = pside(bc_flow,G,'XMin', 0.1 , 'sat', [1]);
        bc_flow = pside(bc_flow,G,'XMax', 0   , 'sat', [1]);
        bc = bc_flow;
        bc.c = 0.*ones(size(bc.sat,1), 1);
        
        % create 2 schedueles: one for flow one for transport
        schedule_flow = simpleSchedule(repmat(dT_flow,1,n_flow),'bc', bc_flow);
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
%         d = 0.* d + 15;
        
        model.K = d;
        
        %% Fixed point solver
        % model only flow:
        model_flow = RichardsEquationModel(G, rock, fluid, 'K', d);
        
        
        % solve flow
        tic
        [~,states_flow, report_flow] = simulateScheduleAD(state0, model_flow, schedule_flow, 'nonlinearsolver', nls);
        time_flow = toc;
        
        v = faceFlux2cellVelocity(G,states_flow{n_flow}.flux);
        
        for l = 1:J
            Vx(l,:) = v(((l-1)*I+1:l*I),1);
            Vy(l,:) = v(((l-1)*I+1:l*I),2);
            p_mrst(l,:) = states_flow{n_flow}.pressure(((l-1)*I+1:l*I));
        end
        
        
        
        Richardsmodel = RichardsEquationModelSequential(G, rock, fluid,  'Newton', 1, 'K', d);   % in the splitting solver we define the L_p and L_c separately
        Transportmodel = SequentialTransportEquationModelnewBC(G, rock, fluid, 'Newton', 1,'K', d, 'D', D);%
           % 'p_flow', states_flow{n_flow}.pressure, 'vw', states_flow{n_flow}.vw);  
        model = SequentialRichardsTransportModel(states_flow{n_flow}, Richardsmodel, Transportmodel,'D', D, 'K', d);
        state0.pressure = states_flow{n_flow}.pressure;
        
        % model transport
%         model = SequentialTransportEquationModelnewBC(G, rock, fluid, 'Newton', 1,...
%             'D', D , 'K', d,'p_flow', states_flow{n_flow}.pressure, 'vw', states_flow{n_flow}.vw);    
%         model.diffusion = 1;
        tic 
        [~,states, report] = simulateScheduleAD(state0, model, schedule, 'nonlinearsolver', nls);
        time_transp = toc;
               
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
                    %xr=(x*dx-x0); yr=(y*dy-y0); % from 'testBGRW_const.m'
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
       save(['Partial_SplittingTentative0810_RandK_MRST_Mesh(',num2str(dx),',',num2str(dy),')_n(',num2str(n),')_3_num',num2str(numb),'.mat'], 'dx', 'Pe', 'Mx', 'My','Varx', 'Vary', 'time_flow', 'time_transp');
    time_realiz = toc
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

        
    end


 %save(['ConstantK_MXEtc_MRST_Mesh(',num2str(I),',',num2str(J),')_n(,',num2str(n),').mat'], 'Mx', 'Dx', 'My', 'Dy');

        
    %%

    
    eps_D1=sqrt(dt)*norm(Dx-D1)/D1/T
    eps_D2=sqrt(dt)*norm(Dy-D2)/D2/T
    
    dx
    dy
    Pe
    
    time_steps = TR;

    %Rand = randi(100000,1)
%save(['ConstK_MXEtc_MRST_Mesh(',num2str(I),',',num2str(J),')_n(,',num2str(n),').mat'], 'dx', 'Pe', 'Mx', 'Dx', 'My', 'Dy', 'eps_D1', 'eps_D2');


     %Rand = randi(1000000,1)
%save(['SplittingTentative0810_Reduced_RandK_MRST_Mesh(',num2str(dx),',',num2str(dy),')_n(',num2str(n),').mat'], 'dx', 'Pe', 'Mx', 'Dx', 'My', 'Dy', 'eps_D1', 'eps_D2');

% save(['Test20','/tentative(',num2str(I),',',num2str(J),')_n(,',num2str(n),').mat'], 'dx', 'Pe', 'Mx', 'Dx', 'My', 'Dy', 'eps_D1', 'eps_D2');


end

toc
