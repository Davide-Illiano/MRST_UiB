clc; clear; close all;
%% Preliminaries
gravity reset
gravity on
T          = 40*year();
stopInject = 10*year();
dT         = .1*year();%.1*year();
dTplot     = 1*dT;
injectTop = false;

rate = 1.4e2*meter^3/day;

% Well radius
r = 0.1;
% Fluid data at p = 300 bar
muw = 0.30860;  rhow = 975.86; sw    = 0.1;
muc = 0.056641; rhoc = 686.54; srco2 = 0.2;
kwm = [0.2142 0.85];
kwm = [1 1];

v = [-30 20];
dims = [30, 1, 3];
physdims=[1000, 1000,40];
% dims = [10, 10, 10];
%mrstModule add internal/vertical-equil

%% Create grids
% Create a logically Cartesian grid and add in geometry data
grdecl = simpleGrdecl(dims, 0, 'undisturbed', true, 'physDims', physdims);
G_full = processGRDECL(grdecl);

clear ijk;
[ijk{1:3}] = ind2sub(G_full.cartDims, G_full.cells.indexMap(:));
ijk        = [ijk{:}];

mi = max(ijk(:,1));
% % Define a 3D region in the middle of the domain
r3D = ismember(ijk(:,1), round(mi/3):round(2*mi/3) );
% Define a 2D region in the middle of the domain
% r3D = ~ismember(ijk(:,1), round(mi/3):round(2*mi/3) );
clf;
plotGrid(G_full, 'facea', 0)
plotGrid(G_full, r3D)

% G = cartGrid(dims, dims*100);

% Perturb the z coordinates based on the x coordinates
perturb = @(x) (.025*(max(x) - min(x))*sin(x)  + .25*(x-min(x)))*0.0;


G_full.nodes.coords(:,3) = G_full.nodes.coords(:,3); + perturb(G_full.nodes.coords(:,1));

G_full = computeGeometry(G_full);
% Create top surface grid
[G_top, G_full] = topSurfaceGrid(G_full);

[G_coupled, region3D,grdecl_coupled] = make2D3Dgrid(grdecl, r3D, 'g3D', G_full,  'g_top', G_top);
G_coupled.nodes.coords(:,3) = G_coupled.nodes.coords(:,3) + perturb(G_coupled.nodes.coords(:,1));

%% MSFV
% p = partitionUI(G_top, [4,1]);
% CG = generateCoarseGrid(G_top, p);
% CG = coarsenGeometry(CG);
% DG = partitionUIdual(CG, [4, 1]);
% clf;
% plotDual(G_top, DG)
% outlineCoarseGrid(G_top, p, 'red')
% view(0,90), axis tight
%%
clf;
subplot(3,1,1)
plotGrid(G_full);
subplot(3,1,2)
plotGrid(G_top);
subplot(3,1,3)
plotGrid(G_coupled);
%% Show grids
h = figure(1); clf;
% Show the full 3D grid
subplot(2,1,1)
plotGrid(G_full);
title('Full grid')
axis equal tight; view(v)
% set(gca, 'ZDir', 'normal')
% Show the top grid
subplot(2,1,2)
plotGrid(G_top, 'FaceColor', 'red')
title('Surface grid')
axis equal tight; view(v)
% set(gca, 'ZDir', 'normal')
%% Fluid and petrophysical properties
rock.perm = 25*milli*darcy*ones(G_full.cells.num,1);
% rock.perm(:,2) = rock.perm(:,2)./1000;
% rock.perm = (1 + 24*rand(G.cells.num,1)).*milli*darcy;
rock.poro = 0.3*ones(G_full.cells.num,1);

rock2D    =    averageRock(rock, G_top);

rock_coupled = rockCoupled(G_coupled, rock, rock2D, region3D);



S_2D_h       = computeMimeticIP(G_top, rock2D, 'Innerproduct','ip_simple');
rock2D_s=rock2D;
rock2D_s.perm=rock2D.perm.*G_top.cells.H;
S_2D_s       = computeMimeticIP(G_top, rock2D_s, 'Innerproduct','ip_simple');
S          = computeMimeticIP(G_full, rock, 'Innerproduct','ip_simple');
%T_coupled = computeTransECLIPSE(G_coupled, rock_coupled,struct('GRID',grdecl_coupled));
T_coupled = computeCoupledTrans(G_coupled, rock_coupled);
T_msfv = computeTrans(G_top, rock2D);
preComp = initTransportVE(G_top, rock2D);

%%
clear ijk
[ijk{1:3}] = ind2sub(G_full.cartDims, G_full.cells.indexMap(:));
y_middle = double(median(double(ijk{2})));
x_min    = double(min(double(ijk{1})));
x_max    = double(max(double(ijk{1})));
% Injector in the left side of the domain
wI = ijk{1} == x_min + 3 & ijk{2} == y_middle;
if injectTop
    wI = wI & ijk{3} == max(ijk{3});
end
% Producer in the right side of the domain
wP       = ijk{1} == x_max - 3 & ijk{2} == y_middle;

% Injector
% W = addWell([], G_full, rock, find(wI),...
%    'Type', 'rate', 'Val', rate, 'Radius', r, 'comp_i', [1,0], 'name', 'I');
% % Producer
% W = addWell(W, G_full, rock, find(wP),...
%    'Type', 'bhp', 'Val', 0, 'Radius', r, 'comp_i', [1,0], 'name', 'P');

W = verticalWell([], G_full, rock, x_min, y_middle, [],...
   'Type', 'rate', 'Val', rate, 'Radius', r, 'comp_i', [1,0], 'name', 'I');
% Producer
W = verticalWell(W, G_full, rock, x_max, y_middle, [],...
   'Type', 'bhp', 'Val', 0, 'Radius', r, 'comp_i', [1,0], 'name', 'P');



% Convert the wells to 2D equivialent
W_2D_h = convertwellsVE(W, G_full, G_top, rock2D,'ip_simple');
W_2D_s = convertwellsVE_s(W, G_full, G_top, rock2D,'ip_simple');

% Create the same wells for the coupled grid. Note that while essentially
% similar to the 3D grid, it has fewer cells and its own ijk-indices.

clear ijk
[ijk{1:3}] = ind2sub(G_full.cartDims, G_coupled.cells.indexMap(:));
y_middle = double(median(double(ijk{2})));
x_min    = double(min(ijk{1}));
x_max    = double(max(ijk{1}));

W_coupled = verticalWell([],        G_coupled, rock_coupled, x_min, y_middle, [],...
   'Type', 'rate', 'Val', rate, 'Radius', r, 'comp_i', [1,0], 'name', 'I','InnerProduct','ip_simple');

W_coupled = verticalWell(W_coupled, G_coupled, rock_coupled, x_max, y_middle, [],...
   'Type', 'bhp', 'Val', 0, 'Radius', r, 'comp_i', [1,0], 'name', 'P','InnerProduct','ip_simple');


clf; 
plotGrid(G_full, 'edgea', .1, 'facea', .1);
% plotWell(G, W, 'radius', 10);
plotGrid(G_full, vertcat(W.cells), 'facec', 'red');
axis equal tight; view(v)
%%
mu = [muc muw] .* centi*poise; rho = [rhoc rhow] .* kilogram/meter^3;

fluid_H = initVEFluidHForm(G_top, 'mu' , mu, ...
                                 'rho', rho, ...
                                 'sr', srco2, 'sw', sw, 'kwm', kwm);

fluid_S = initSimpleVEFluid_s('mu' , mu , 'rho', rho, ...
                                'height'  , G_top.cells.H,...
                                'sr', [srco2, sw],'kwm',kwm);

fluid_coupled = initCoupledVEFluid('mu' , mu , 'rho', rho, ...
                               'sr', [srco2, sw], 'region3D', region3D, ...
                               'n', [1, 1], 'g', G_coupled);  
                       
                       
fluid_full = initCoreyFluid('mu' , mu , 'rho', rho, 'sr', [srco2, sw], 'kwm', [1 1], 'n', [1 1]);  
%%
% H formulation
sol_H = initResSolVE(G_top, 0);
sol_H.wellSol = initWellSol(W, 0);
% S formulation
sol_S = initResSolVE(G_top, 0);
sol_S.wellSol = initWellSol(W, 0);
% S formulation (MsFV)
% sol_msfv = initResSolVE(G_top, 0);
% sol_msfv.wellSol = initWellSol(W, 0);
% Coupled formulation
sol_coupled  = initResSolVE(G_coupled, 0);
sol_coupled.wellSol = initWellSol(W, 0);
% Full grid
sol_full = initResSol(G_full, 0);
sol_full.wellSol = initWellSol(W, 0);


% % %% Prepare plot
% % 
% opts = {'slice', double([x_min + 3 y_middle]), 'Saxis', [0 1-fluid_H.sw], 'maxH', max(G.nodes.coords(:,3)), ...
%         'view', v, 'wireH', true, 'wireS', true};
% % % plotPanelVE(G, G_2D, W, sol_S, 0.0, [0 0 0 0], opts{:});
% plotPanelVE(G, G_2D, W, sol_H, 0.0, [0 0 0 0], opts{:});
%% Main loop
% Run the simulation using a sequential splitting with pressure and
% transport computed in separate steps. 
t = 0;
fprintf(1,'\nSimulating %d years on injection',convertTo(stopInject,year));
fprintf(1,' and %d years of migration\n', convertTo(T-stopInject,year));
fprintf(1,'Time: %4d years', convertTo(t,year));
tic;

while t<T
    
   % Advance solution:
   % First we compute pressures
   sol_H =    solveIncompFlowVE  (sol_H, G_top, S_2D_h, rock, fluid_H, 'wells', W_2D_h);
   sol_S =    solveIncompFlow(sol_S, G_top, S_2D_s,       fluid_S, 'wells', W_2D_s);
   sol_S.pressure = sol_S.pressure;
   sol_coupled = incompTPFAVE_coupled(sol_coupled, G_coupled, T_coupled, fluid_coupled, ...
                            'wells', W_coupled, 'region3D', region3D);
   sol_full = solveIncompFlow(sol_full , G_full, S, fluid_full, 'wells', W);
%    sol_msfv = solveMSFV_TPFA_Incomp(sol_msfv, G_top, CG, T_msfv, fluid_S, ...
%                                    'wells', W_2D, 'Dual', DG...
%                                     );
   
   
   % which is then used to advance saturations
   sol_H = explicitTransportVE(sol_H, G_top, dT, rock, fluid_H, ...
                              'wells', W_2D_h, ...
                              'preComp', preComp, 'Verbose', false, 'computeDt', true);
                          
   
   sol_S       = implicitTransport(sol_S, G_top, dT, rock2D, fluid_S, 'wells', W_2D_s,'Verbose',true);
%    sol_msfv    = implicitTransport(sol_msfv, G_top, dT, rock2D, fluid_S, 'wells', W_2D);
   sol_full    = implicitTransport(sol_full, G_full, dT, rock, fluid_full, 'wells', W,'Verbose',true);
   %
   % Here we use options only in used in the jacobian.
   % To avoid warning in mearge_options in implicitTransport which do not know this options
   % we trun it off
   warning('off','implicitTransport:Option:Unsupported')
   sol_coupled = implicitTransport(sol_coupled, G_coupled, dT, rock_coupled,...
                            fluid_coupled,'vert_avrg',true, ...
                            'wells',W_coupled,'region3D', region3D, 'verbose', true);
   warning('on','implicitTransport:Option:Unsupported')                   

   
   % Reconstruct 'saturation' defined as s=h/H, where h is the height of
   % the CO2 plume and H is the total height of the formation
   sol_H.s = height2Sat(sol_H, G_top, fluid_H);
   % For the s formulation, construct height from the saturation values
   sol_S.h = fluid_S.sat2height(sol_S);
   %... monotone?
   sol_S.h_max = max(sol_S.h, sol_S.h_max);
   
   %
   sol_full.h = sat2height(sol_full.s, G_top, rock);
   
   
   assert( max(sol_H.s(:,1))<1+eps && min(sol_H.s(:,1))>-eps );
   
   
   
   
   
   t = t + dT;

   % Check if we are to stop injecting
   if t>= stopInject
      W_2D  = []; bcVE = []; dT = 5*year(); dTplot = dT;
   end

   % Compute trapped and free volumes of CO2
   fprintf(1,'\b\b\b\b\b\b\b\b\b\b%4d years', int32(convertTo(t,year)));
%    fprintf(1, '%d years\n', int32(convertTo(t, year)));
   vols = volumesVE(G_top, sol_H, rock2D, fluid_H);
   
%    figure(1); clf; hold on
%    x = G_top.cells.centroids(:,1);
%    plot(x, [sol_S.h sol_H.h, sol_full.h])
%    legend({'S formulation', 'H formulation', 'Full 3d'})
%    axis([0 max(G_full.nodes.coords(:,1)) 0 max(G_top.cells.H)])
%    drawnow
   figure(2);
   subplot(2,2,1)
   
   [s h] = normalizeValuesVE(G_top, sol_S, fluid_S);
   plotCellData(G_full, s)
%    plotCellData(G_full, height2Sat(sol_S, G_top, fluid_S))
   axis square tight; view([0 0]);colorbar
   title('S formulation')
   ca2=caxis();
   
%    subplot(2,2,2)
%    [s h] = normalizeValuesVE(G_top, sol_msfv, fluid_S);
%    plotCellData(G_full, s)
% %    plotCellData(G_full, height2Sat(sol_S, G_top, fluid_S))
%    axis square tight; view([0 0])
%    title('MsFV S formulation')
   
   subplot(2,2,2)
   [s h] = normalizeValuesVE(G_top, sol_H, fluid_H);
   plotCellData(G_full, s)
   plotCellData(G_full, height2Sat(sol_H, G_top, fluid_H))
   axis square tight; view([0 0]);colorbar
   title('H formulation');
   caxis(ca2);
   
   subplot(2,2,3)
   plotCellData(G_full, sol_full.s)
   title('Full 3D')
   axis square tight; view([0 0]);colorbar
   caxis(ca2);
   
   subplot(2,2,4)
   [s h] = normalizeValuesVE(G_top, sol_coupled, fluid_coupled,...
                            'CoupledGrid', G_coupled);
   plotCellData(G_full, s)
%    plotCellData(G_coupled, sol_coupled.s)
   title('Coupled 3D TPFA')
   axis square tight; view([0 0]) ; colorbar
   caxis(ca2);
   drawnow
   
   
   figure(3);
   subplot(2,2,1)
   plotCellData(G_top, sol_S.pressure);colorbar
   axis square tight; view(0,90)
   title('S formulation')
   ca3=caxis();
   
   subplot(2,2,2)
   plotCellData(G_top, sol_H.pressure);colorbar
   axis square tight; view(0,90)
   title('H formulation')
   caxis(ca3);
   
   subplot(2,2,3)
   plotCellData(G_full, sol_full.pressure);colorbar
   title('Full 3D')
   axis square tight; view(0,90)
    caxis(ca3);
   subplot(2,2,4)
   plotCellData(G_coupled, sol_coupled.pressure);colorbar
   title('Coupled 3D TPFA')
   axis square tight; view(0,90)
   caxis(ca3);
   drawnow;
   %{
   figure(4);
   clf;
   subplot(1,2,1);
   plotCellData(G_coupled, sol_coupled.mob(:,1))
   plotGrid(G_coupled, 'facea', 0);
   axis square tight; view([0 0])
   subplot(1,2,2);
   plotCellData(G_coupled, sol_coupled.mob(:,2))
   plotGrid(G_coupled, 'facea', 0);
   axis square tight; view([0 0])
   %}
   drawnow
   
%    fprintf('Pressure_%s Min: %s \t Max:%s\n', 'S', min(sol_S.pressure), max(sol_S.pressure))
%    fprintf('Pressure_%s Min: %s \t Max:%s\n', 'H', min(sol_H.pressure), max(sol_H.pressure))
%    fprintf('Pressure_%s Min: %s \t Max:%s\n', 'F', min(sol_full.pressure), max(sol_full.pressure))
%    Plotting
%    if mod(t,dTplot)~= 0 && t<T,
%       continue
%    else
% %       plotPanelVE(G, G_2D, W, sol_S, t, [vols sum(vols)], opts{:});
% %       plotPanelVE(G, G_2D, W, sol_H, t, [vols sum(vols)], opts{:});
%       plotPanelVE(G, G_2D, W, sol_full, t, [vols sum(vols)], opts{:});
%       drawnow
%    end
% pause
%    break
end
fprintf(1,'\n\n');
%% Check fluxes etc
bndf = G_coupled.cells.faces(G_coupled.facesBnd.cellFace3D);
% sol_coupled.flux(bndf)
sum(sol_coupled.flux(bndf).*G_coupled.faces.areas(bndf))
mean(sol_coupled.flux(sol_coupled.flux>0))

% cellNo = rldecode(1:G_coupled.cells.num, double(diff(G_coupled.cells.facePos)), 2).';

figure(5); clf;
% f = boundaryFaces(G);
% f = find(G_coupled.faces.tag > 0);
f = setdiff(find(G_coupled.faces.normals(:,3) == 0), boundaryFaces(G_coupled) );
flux = sol_coupled.flux./G_coupled.faces.areas;
% flux = T_coupled;
% flux = flux(G_coupled.cells.faces(:,1));
plotFaces(G_coupled, f, flux(f), 'facea', .2);
colorbar;
std(flux(f))

figure(6); clf;
% f = boundaryFaces(G);
% f = find(G_coupled.faces.tag > 0);
flux = sol_full.flux;
f = setdiff(find(G_full.faces.normals(:,3) == 0), boundaryFaces(G_full) );
% flux = computeTrans(G_full, rock);
% flux = flux(G_full.cells.faces(:,1));
plotFaces(G_full, f, flux(f), 'facea', .2);
colorbar;
std(flux(f))
%%
tarea = T_coupled./G_coupled.faces.areas(G_coupled.cells.faces(:,1));
tt = unique(tarea);
figure(7); clf;

plotFaces(G_coupled, G_coupled.cells.faces((tarea == tt(1)), 1), 'FaceC', 'r')
% plotFaces(G_coupled, G_coupled.cells.faces((tarea == tt(2)), 1), 'FaceC', 'g')
% plotFaces(G_coupled, G_coupled.cells.faces((tarea == tt(3)), 1), 'FaceC', 'b')
