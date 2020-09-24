function [eqs, fluxes] = eqsfiOGExplicitWellsCon(state0, state, dt, G, W, s, f, varargin)

opt = struct('Verbose', mrstVerbose, ...
             'reverseMode', false,...
             'scaling', [],...
             'resOnly', false,...
             'fluxes',false,...
             'bc',[]);
                      %'history', [],...

opt = merge_options(opt, varargin{:});

if ~isempty(opt.scaling)
    scalFacs = opt.scaling;
else
    scalFacs.rate = 1; scalFacs.pressure = 1;
end

%hst = opt.history;

% current variables: ------------------------------------------------------
p    = state.pressure;
sG   = state.s(:,2);
%sGmax=state0.smax(:,2);
if(isfield(state0,'smax'))
    sGmax=max(state.s(:,2),state0.smax(:,2));
else
    sGmax=[];
end
pBHP = vertcat(state.wellSol.bhp);
qGs  = vertcat(state.wellSol.qGs);
qOs  = vertcat(state.wellSol.qOs);

% previous variables ------------------------------------------------------
p0  = state0.pressure;
sG0 = state0.s(:,2);
%--------------------------------------------------------------------------


%Initialization of independent variables ----------------------------------

zw = zeros(size(pBHP));
if ~opt.resOnly
   % ADI variables needed since we are not only computing residuals.

   if ~opt.reverseMode,

      [p, sG, qGs, qOs, pBHP] = initVariablesADI(p, sG, qGs, qOs, pBHP);

   else

      [p0, sG0, zw, zw, zw] = ...
         initVariablesADI(p0, sG0,          ...
                          zeros(size(qGs)), ...
                          zeros(size(qOs)), ...
                          zeros(size(pBHP)));                          %#ok

   end
end

g  = norm(gravity);
[Tw, dzw, Rw, wc, perf2well, pInx, iInxG, iInxO] = getWellStuffOG(W);

%--------------------
%check for p-dependent tran mult:
trMult = 1;
if isfield(f, 'tranMultR'), trMult = f.tranMultR(p); end

%check for p-dependent porv mult:
pvMult = 1; pvMult0 = 1;
if isfield(f, 'pvMultR')
    pvMult =  f.pvMultR(p);
    pvMult0 = f.pvMultR(p0);
end



% -------------------------------------------------------------------------
% [krW, krO] = f.relPerm(sW);
pcOG = 0;
if isfield(f, 'pcOG')
    pcOG  = f.pcOG(sG,'sGmax',sGmax);
else
    pcOG = zeros(G.cells.num,1);
end
pcOGw = pcOG(wc);
bG     = f.bG(p+pcOG);
muG     = f.muG(p+pcOG);
rhoG   = bG.*f.rhoGS;
% rhoW on face, avarge of neighboring cells (E100, not E300)
rhoGf  = s.faceAvg(rhoG);
if(~any(strcmp(G.type,'topSurfaceGrid')))
    dpG     = s.grad(p+pcOG) - g*(rhoGf.*s.grad(G.cells.centroids(:,3)));
else
    dpG     = s.grad(p+pcOG) - g*(rhoGf.*s.grad(G.cells.z));
end

if(isfield(f,'con_fluid'))
    upc = (double(dpG)>=0);
    sGu= s.faceUpstr(upc, sG);
    sGmaxu= s.faceUpstr(upc, sGmax);
    mubGu= s.faceUpstr(upc, bG./muG);
    con_krG = f.con_krG(sGu,'sGmax',sGmaxu);
    bGvG = mubGu.*con_krG.*s.T.*dpG;
else
    krG = f.krG(sG,'sGmax',sGmax);
    krO = f.krOG(1-sG,'sGmax',sGmax);
    mobG   = trMult.*krG./f.muG(p+pcOG);
    % water upstream-index
    upc = (double(dpG)>=0);
    bGvG = s.faceUpstr(upc, bG.*mobG).*s.T.*dpG;
end

% oil props
bO     = f.bO(p);
muO    = f.muO(p);
rhoO   = bO.*f.rhoOS;
rhoOf  = s.faceAvg(rhoO);


if(~any(strcmp(G.type,'topSurfaceGrid')))
    dpO    = s.grad(p) - g*(rhoOf.*s.grad(G.cells.centroids(:,3)));
else
    dpO    = s.grad(p) - g*(rhoOf.*s.grad(G.cells.z));
end

if(isfield(f,'con_fluid'))
    upc = (double(dpO)>=0);
    sGu= s.faceUpstr(upc, sG);
    sGmaxu= s.faceUpstr(upc, sGmax);
    mubOu= s.faceUpstr(upc, bO./muO);
    con_krO = f.con_krO(1-sGu,'sGmax',sGmaxu);
    bOvO = mubOu.*con_krO.*s.T.*dpO;
else
    %mobO   = trMult.*krO./f.BOxmuO(p);
    mobO   = trMult.*krO./f.muO(p);
    % oil upstream-index
    upc = (double(dpO)>=0);
    bOvO   = s.faceUpstr(upc, bO.*mobO).*s.T.*dpO;
end

%WELLS ----------------------------------------------------------------
bGw     = bG(wc);
bOw     = bO(wc);
if(isfield(f,'con_fluid'))
    sGw    = sG(wc);
    muGw   =muG(wc);
    muOw   =muO(wc);
    mobGw  = f.con_krGw(sGw)./muGw;
    mobOw  = f.con_krOw(1-sGw)./muOw;
else
    mobGw  = mobG(wc);
    mobOw  = mobO(wc);
end
%producer mobility
bGmobGw  = bGw.*mobGw;
bOmobOw  = bOw.*mobOw;

%set water injector mobility: mobw = mobw+mobo+mobg, mobo = 0;
if(isfield(f,'con_fluid'))
    % this is different since this is taken intoaccunt in the con_kr?w
    % need to fix the posibility of zero mobility
    %bGmobGw(iInxG) = bGw(iInxG).*(mobGw(iInxG)+mobOw(iInxG).*(double(mobGw(iInxG))<0.5));
    bGmobGw(iInxG) = bGw(iInxG).*(mobGw(iInxG)+mobOw(iInxG));
    %bGmobGw(iInxG) = bGw(iInxG).*1./(1./mobGw(iInxG) + 1./mobOw(iInxG));
    bOmobOw(iInxG) = 0;
    bGmobGw(iInxO) = 0;
    %bOmobOw(iInxO) = bOw(iInxO).*(mobOw(iInxO)+mobGw(iInxO).*(double(mobOw(iInxO))<0.5) );
    bOmobOw(iInxO) = bOw(iInxO).*(mobOw(iInxO)+mobGw(iInxO) );
else
if(true)
    % this sum may have wrong number in the end point
    bGmobGw(iInxG) = bGw(iInxG).*(mobGw(iInxG) + mobOw(iInxG));
    %bGmobGw(iInxG) = bGw(iInxG).*1./(1./mobGw(iInxG) + 1./mobOw(iInxG));
    bOmobOw(iInxG) = 0;
    bGmobGw(iInxO) = 0;
    bOmobOw(iInxO) = bOw(iInxO).*(mobGw(iInxO) + mobOw(iInxO));
else
   bGmobGw(iInxG) = bGw(iInxG)./1e-3;
   bOmobOw(iInxG) = 0;
   bGmobGw(iInxO) = 0;
   bOmobOw(iInxO) = bOw(iInxO)./1e-3;
end
end
pw  = p(wc);

bGqG  = -bGmobGw.*Tw.*(pBHP(perf2well) - pw + 0.0*pcOGw + g*dzw.*rhoG(wc));
bOqO  = -bOmobOw.*Tw.*(pBHP(perf2well) - pw + g*dzw.*rhoO(wc));

if(~isempty(opt.bc))
    assert(all(strcmp(opt.bc.type,'pressure')));
    Tbc=s.T_all(opt.bc.face);
    bc_cell=sum(G.faces.neighbors(opt.bc.face,:),2);
    %bc_cell_loc=false(G.cells.num,1);
    %bc_cell_loc(bc_cell_nr)=true;
    % assume cappillary pressure zero outside
    dzbc=(G.cells.z(bc_cell)-G.faces.z(opt.bc.face));
    pObc=p(bc_cell);rhoObc=rhoO(bc_cell);
    dpbc_o=opt.bc.value-pObc+g*(rhoObc.*dzbc);
    pGbc= pObc+pcOG(bc_cell);rhoGbc=rhoG(bc_cell);
    dpbc_g=opt.bc.value-(pGbc)+g*(rhoGbc.*dzbc);
    bGmobGbc = bG(bc_cell).*mobG(bc_cell);
    bOmobObc = bO(bc_cell).*mobO(bc_cell);
    bGmobGbc(dpbc_g>0)=0;
    %bGmobGbc(dpbc_o>0)=0;
    %bOmobObc(dpbc_g>0)=0;
    if(any(dpbc_o>0))
        bOmobObc(dpbc_o>0)=bO(bc_cell(dpbc_o>0)).*(mobO(bc_cell(dpbc_o>0))+mobG(bc_cell(dpbc_o>0)));
    end
    bGqGbc  = -bGmobGbc.*Tbc.*(dpbc_g);
    bOqObc  = -bOmobObc.*Tbc.*(dpbc_o);
else
    bc_cell=[];
end

% EQUATIONS ---------------------------------------------------------------


% oil:
eqs{1} = (s.pv/dt).*( pvMult.*bO.*(1-sG) - pvMult0.*f.bO(p0).*(1-sG0) ) + s.div(bOvO);
%eqs{1}(wc) = eqs{1}(wc) + bOqO;

% water:
eqs{2} = (s.pv/dt).*( pvMult.*bG.*sG - pvMult0.*f.bG(p0).*sG0 ) + s.div(bGvG);
%eqs{2}(wc) = eqs{2}(wc) + bGqG;
if(~isempty(W))
    eqs{1}(wc) = eqs{1}(wc) + bOqO;
    eqs{2}(wc) = eqs{2}(wc) + bGqG;
end
if(~isempty(opt.bc))
    eqs{1}(bc_cell)  = eqs{1}(bc_cell) + bOqObc;
    eqs{2}(bc_cell)  = eqs{2}(bc_cell) + bGqGbc;
    assert(sum(double(bGqGbc))>=0)
end
% well equations
zeroW = 0*zw;

eqs{3} = Rw'*bGqG + qGs + zeroW;
eqs{4} = Rw'*bOqO + qOs + zeroW;

% Last eq: boundary cond
eqs{5} = handleBC(W, pBHP, [] , qOs, qGs, scalFacs) + zeroW;
if(opt.fluxes)
    upc = (double(dpO)>=0);
    bOf   = s.faceUpstr(upc, bO);
    upc = (double(dpG)>=0);
    bGf   = s.faceUpstr(upc, bG);
    fluxes=struct('bGvG', double(bGvG),...
                 'vG' , double(bGvG./bGf),...
                'bOvO', double(bOvO),...
                 'vO' , double(bOvO./bOf),...
                'bGqG', double(bGqG),...
                 'qG' , double(bGqG./bGw),...
                'bOqO', double(bOqO),...
                 'qO' , double(bOqO./bOw));

else
    fluxes=[];
end

end
%--------------------------------------------------------------------------










