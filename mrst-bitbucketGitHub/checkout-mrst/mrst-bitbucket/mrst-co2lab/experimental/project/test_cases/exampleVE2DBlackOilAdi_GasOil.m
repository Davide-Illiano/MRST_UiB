%% VE simulation in a standard black-oil solver
%  In this example we show how to set up a standard format black-oil
%  model that can be used to simulate a VE model. For the actual
%  simulation,  we use the fully-implicit solver in MRST from the 'ad-fi'
%  module, which is based on automatic differentiation. 

try
   require deckformat ad-fi
catch %#ok<CTCH>
   mrstModule add deckformat ad-fi
end

%% Parameters for the simulation
gravity on
close all;
[nx,ny,nz] = deal(40, 1, 1);     % Cells in Cartsian grid
[Lx,Ly,H]  = deal(2000,1000,15); % Physical dimensions of reservoir
total_time = 5*year;             % Total simulation time
nsteps     = 10;                 % Number of time steps in simulation
dt         = total_time/nsteps;  % Time step length
perm       = 100;                % Permeability in milli darcies
phi        = 0.1;                % Porosity
depth      = 1000;               % Initial depth
ipress     = 240;                % Initial pressure

%% Create input deck and construct grid
% Create an input deck that can be used together with the fully-implicit
% solver from the 'ad-fi' module. Since the grid is constructed as part of
% setting up the input deck, we obtain it directly. 
[deck, G] = sinusDeckAdi_GasOil([nx ny nz], [Lx Ly H], nsteps, dt, ...
                         -.1*pi/180, depth, phi, perm, ...
                         (H*phi*Lx*Ly)*0.2*day/year, ipress);

% Alternatively, we could read deck from file and construct the grid
% deck = readEclipseDeck( ...
%    fullfile(VEROOTDIR,'data','decks','sinusDeckAdi.DATA');
% G = initEclipseGrid(deck);

figure, plotGrid(G),view([0 -1 0]), box on
                      
                 
%% Initialize data structures
% First, we convert the input deck to SI units, which is the unit system
% used by MRST. Second, we initialize the rock parameters from the deck;
% the resulting data structure may have to be post-processed to remove
% inactive cells. Then we set up the fluid object and tell the ad-fi solver
% that that we are working with an oil-gas system.
deck  = convertDeckUnits(deck);
rock  = initEclipseRock(deck);
rock  = compressRock(rock, G.cells.indexMap);
fluid = initDeckADIFluid(deck);
% set the capillary pressure and the VE relperms explicitely
Gt = topSurfaceGrid(G);
Gt.parent=G;
rock3D=rock;
rock = averageRock(rock3D, Gt);
rock.parent=rock;
clear G;

l_fac=0e-6;
fluid.bW=@(po,varargin) (po-200*barsa)*l_fac+1;
fluid.BW=@(po,varargin) 1./fluid.bW(po,varargin);
fluid.bG=@(pg,varargin) pg*0.0+1;
fluid.BG=@(pg,varargin) pg*0.0+1;
fluid.muW=@(po,varargin) 0.4e-3*(po*0+1);
fluid.muG=@(pg,varargin) 1e-4*(pg*0+1);

fluid_case='hystersis';
fluid_case='simple';
fluid_case='cap_linear';
fluid_case='cap_1D_table_P';
%fluid_case='cap_1D_table_kscaled' 
%fluid_case='integrate'
%fluid_case='sharp_interface' 
res_gas = 0.1;
res_water= 0.2;
switch fluid_case
    case 'simple'
       fluid.krG=@(sg,varargin) sg;
       fluid.krW=@(so,varargin) so;
       fluid.pcWG=@(sg, p, varargin) norm(gravity)*(fluid.rhoWS.*fluid.bW(p)-fluid.rhoGS.*fluid.bG(p)).*(sg).*Gt.cells.H;
       fluid=rmfield(fluid,'relPerm');
       res_gas=0;
    case 'sharp_interface'    
       fluid = addVERelperm(fluid, Gt, 
                            'res_water',res_water,...
                            'res_gas',res_gas);
    case 'integrate'
        res_gas = 0.0;
        fluid = addVERelpermCap(fluid,'alpha',1,'beta',1,'cap_scale',100*barsa,'H',Gt.cells.H,'kr_pressure',false);
    case 'cap_linear'
        fluid = addVERelpermCapLinear(fluid,...
                                      'res_gas',res_gas,...
                                      'res_water',res_water,...
                                      'beta',4,...
                                      'cap_scale',0.5*max(Gt.cells.H)*10*(fluid.rhoWS-fluid.rhoGS),...
                                      'H',Gt.cells.H,'kr_pressure',true);
    case 'cap_1D_table_P'
        drho=400;
        C=max(Gt.cells.H)*1*drho*norm(gravity);
        alpha=(1/2);
        beta = 4;
        samples=100;
        table_co2_1d = makeVEtables('invPc3D', @(p) (C./(p+C)).^(1/alpha),...
            'is_kscaled', false,'kr3D', @(s) s.^beta,...
            'drho', drho,...
            'Gt', Gt, 'samples', samples);
        S_tab=linspace(0,1,10)';
        S_tab_w=S_tab;
        kr_tab_w=S_tab_w;
        table_water_1d=struct('S',1-S_tab_w,'kr',kr_tab_w,'h',[]);
        fluid=addVERelperm1DTablesPressure(fluid,...
                                           'res_water', res_water,...
                                           'res_gas', res_gas,...
                                           'height',Gt.cells.H,...
                                           'table_co2',table_co2_1d,...
                                           'table_water',table_water_1d,...
                                           'kr_pressure',true);
  case 'cap_1D_table_kscaled'        
        fluid.surface_tension=100;  
        kscale=sqrt(0.1/(100*milli*darcy))*fluid.surface_tension;
        drho=400;
        C=max(Gt.cells.H)*0.4*drho*norm(gravity)/kscale;
        alpha=0.5;
        beta = 3;
        samples=100;
        table_co2_1d = makeVEtables('invPc3D', @(p) (C./(p+C)).^(1/alpha),...
                                    'is_kscaled', true,....
                                    'kr3D', @(s) s.^beta,...
                                    'drho', drho,...
                                    'Gt', Gt,...
                                    'samples', samples,'kscale',kscale);
        S_tab=linspace(0,1,10)';
        S_tab_w=S_tab;
        kr_tab_w=S_tab_w;
        table_water_1d=struct('S',1-S_tab_w,'kr',kr_tab_w,'h',[]);
        fluid=addVERelperm1DTablesPressure(fluid,...
                                           'res_water', res_water,...
                                           'res_gas', res_water,...
                                           'height', Gt.cells.H,...
                                           'table_co2',table_co2_1d,...
                                           'table_water',table_water_1d,...
                                           'rock',rock,...
                                           'kr_pressure',true);                                      
                               
    case 'hystersis'
        
        fluid = addVERelperm(fluid, Gt, 'res_water',0,'res_gas',res_gas);
    otherwise
       disp('Use deck as fluid')
end

s=setupSimCompVe(Gt,rock);
systemOG = initADISystemVE({'Oil', 'Gas'}, Gt, rock, fluid,'simComponents',s,'VE',true);


%% Run the schedule setup in the file
% Before we can run the schedule, we make sure that we have an initial
% hydrostatic pressure distribution. Then we pick the schedule from the
% input deck and start the simulator.
%x0 = initEclipseState(G, deck, initEclipseFluid(deck));
z  = Gt.cells.z;
%x0.flux=zeros(Gt.faces.num,1);
clear x0;
x0.pressure = ipress*barsa +(z(:)-z(end))*norm(gravity)*deck.PROPS.DENSITY(2);
x0.s(:,1)=deck.SOLUTION.SOIL;
x0.s(:,2)=deck.SOLUTION.SGAS;
x0.smax=x0.s;
x0.smin=x0.s;

%x0.s=x0.s(:,[2,1]);
%x0.s=x0.z(:,[2,1]);


[wellSols, states] = runScheduleADIwithVE(x0, Gt, rock, systemOG, deck.SCHEDULE);

%% Plot results
%figure

xc = Gt.cells.centroids(:,1);
zt = Gt.cells.z;
zb = zt + Gt.cells.H;
for nn=1:numel(states)
    figure(1)
    clf
    state=states{nn};
    
    subplot(2,2,1),cla
    title('pressure')
    plotCellData(Gt,state.pressure/barsa);colorbar('horiz'), caxis([100 200])
    
    subplot(2,2,2),cla
    title('saturation')
    plotCellData(Gt,state.s(:,2));colorbar('horiz'), caxis([0 1])
    
    % plot as VE
    subplot(2,2,3),cla
    plot(xc,state.pressure/barsa); set(gca,'YLim',[100 200]);

    subplot(2,2,4),cla,hold on
    patch(xc([1 1:end end]), [zt(end)-10; zt; zt(end)-10],.7*[1 1 1]);
    patch(xc([1 1:end end]), [zb(end)+10; zb; zb(end)+10],.7*[1 1 1]);
    z_free=zt + Gt.cells.H.*(state.s(:,2)-(state.smax(:,2)-state.s(:,2))*res_gas);
    z_res=zt + Gt.cells.H.*state.smax(:,2);
    patch(xc([1:end end:-1:1]), ...
      [z_free; zt(end:-1:1)], 'r')
    patch(xc([1:end end:-1:1]), ...
        [z_res; z_free(end:-1:1)], 'g')
    patch(xc([1:end end:-1:1]), ...
      [z_res; zb(end:-1:1)], 'b')
    set(gca,'YDir','reverse'), axis tight

    %%
    
    %%
    drawnow;
    pause(0.01)
end
%%
figure(2),clf;
z=linspace(0,H,1000);
for nn=1:numel(states)
    state=states{nn};
    p=state.pressure;
    pc=fluid.pcWG(state.s(:,2), p,'sGmax',state.smax(:,2));
    pcmax=fluid.pcWG(state.smax(:,2), p,'sGmax',state.smax(:,2));
    drho=norm(gravity)*(fluid.rhoWS.*fluid.bW(p)-fluid.rhoGS.*fluid.bG(p));
    h=pc./drho;
    h_max=pcmax./drho;   
    cells=1:ceil(nx/9):nx;
    k=0;
    
    for cell=cells;
        k=k+1;
        %subplot(numel(cells),1,k)
        subplot(1,numel(cells),k),cla
        dp=-drho(cell).*(z-h(cell));
        s3D_free=fluid.invPc3D(dp);
        s3D_free(dp<0)=1;
        dp_max=-drho(cell).*(z-h_max(cell));
        s3D_max=fluid.invPc3D(dp_max);
        s3D_max(dp_max<0)=1;
        if(false)
        plot(1-s3D_max,z,1-s3D_free,z);
        axis([0 1.1 0 H])
        else
          xval=[(1-s3D_free)';0];
          yval=[z';H*10];
          patch(xval([1:end end:-1:1]), ...
          [zeros(size(yval)); yval(end:-1:1)], 'r')
          xxval=[(1-s3D_free)'+((1-s3D_max)-(1-s3D_free))'*fluid.res_gas/(1-fluid.res_water);0];
          %xxval=[(1-s3D_max)';0];
          %yyval=[z';H*10];
          patch([xval([1:end]);xxval([end:-1:1])], ...
          [yval(1:end); yval(end:-1:1)], 'g')
          xxxval=[(1-s3D_max)';0];
          yyval=[z';H*10];      
          %xxval=[(1-s3D_max)';0];
          %yyval=[z';H*10];
          patch([xxval([1:end]);xxxval([end:-1:1])], ...
          [yval(1:end); yval(end:-1:1)], [0 1 1])
      
          patch([xxxval([1:end]); ones(size(xval))], ...
          [yval(1:end); yval(end:-1:1)], 'b') 
          axis([0 1 0 H])
          box on
          set(gca,'LineWidth',2)
        end
        
        set(gca,'YDir','reverse')
    end
    %drawnow;
    pause(0.1)
end