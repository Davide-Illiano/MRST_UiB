function [problem, state] = equationsRichardsTransport(state0, state, model, dt, drivingForces, varargin)
opt = struct('Verbose', mrstVerbose, ...
             'reverseMode', false,...
             'resOnly', false,...
             'iteration', -1);

opt = merge_options(opt, varargin{:});

W = drivingForces.W;
s = model.operators;
fluid = model.fluid;

% Properties at current timestep
[p, c, wellSol] = model.getProps(state, 'pressure', 'concentration', 'wellsol');
% Properties at previous timestep
[p0, c0, wellSol0] = model.getProps(state0, 'pressure', 'concentration',  'wellSol');

[wellVars, wellVarNames, wellMap] = model.FacilityModel.getAllPrimaryVariables(wellSol);

% Initialize independent variables.
if ~opt.resOnly
    % ADI variables needed since we are not only computing residuals.
    if ~opt.reverseMode
        [p, c,  wellVars{:}] = initVariablesADI(p, c, wellVars{:});
    else
        wellVars0 = model.FacilityModel.getAllPrimaryVariables(wellSol0);
        [p0, c, wellVars0{:}] = initVariablesADI(p0, c, wellVars0{:}); %#ok
    end
end
primaryVars = {'pressure', 'concentration',  wellVarNames{:}};

% Gravity contribution
gdz = model.getGravityGradient();


bW     = fluid.bW(p);
bW0     = fluid.bW(p0);


rhoW   = bW.*fluid.rhoWS;
rhoWf  = s.faceAvg(rhoW);

dp = s.Grad(p) - rhoWf.*gdz;  
upc  = (double(dp)<=0);

%% Obtain the saturation sW as fraction between the water content theta
%  and the porosity. How to use rock.poro
% Compute transmissibility
theta = fluid.theta;
rock = model.rock;
sW = theta(p,c)./rock.poro;
sW0 = theta(p0,c0)./rock.poro;

% Evaluate relative permeabilities
[krW] = model.evaluateRelPerm({sW});

% Multipliers for properties
[pvMult, transMult, mobMult, pvMult0] = getMultipliers(model.fluid, p, p0);

% Modifiy relperm by mobility multiplier (if any)
krW = mobMult.*krW; 

Kmult = model.K(:); %fluid.Kmult(p, theta(p,c));   
T = s.T.*s.faceUpstr(upc, Kmult);


% Evaluate water properties
[vW, bW, mobW, rhoW, pW, upcw] = getFluxAndPropsWater_BO(model, p, sW, krW, T, gdz);
mobW = 0*mobW + 1;

bW0 = model.fluid.bW(p0);

%
vW = -T.*dp;

bWvW = s.faceUpstr(upc, bW).*vW;


varK=0.1;
K_MEAN=15;
Jslope=0.05;
cvel=exp(varK/2)/(Jslope*K_MEAN); 

VcW = s.faceUpstr(upc,c).*cvel.*bWvW;

VcD = -model.D.*s.Grad(c); %-s.Grad(c).*s.faceAvg(theta);

Vc = VcW + VcD;


% Conservation of mass for water
water = (s.pv/dt).*( bW.*theta(p,c) - bW0.*theta(p0,c0) ) + s.Div(bWvW); %
transport = (s.pv./dt).*(bW.*theta(p,c).*c - bW0.*theta(p0,c0).* c0) + s.Div(Vc);% + c.*R;

eqs = {water, transport};
names = {'water', 'concentration'};
types = {'cell', 'cell'};


% Add in any fluxes / source terms prescribed as boundary conditions.
rho = {rhoW};
mob = {mobW};
sat = {sW};


if model.outputFluxes
    state = model.storeFluxes(state, vW, [], []);
end


[eqs, state] = addBoundaryConditionsAndSources(model, eqs, names, types, state, ...
                                                                 {p}, sat, mob, rho, ...
                                                                 {}, {c}, ...
                                                                 drivingForces);
% Finally, add in and setup well equations
%
[eqs, names, types, state.wellSol] = model.insertWellEquations( eqs, names, ...
                                                     types, wellSol0, wellSol, ...
                                                     wellVars, wellMap, ...
                                                     p, mob, rho, ...
                                                     {}, {c}, ...
                                                     dt, opt);
                                                 
                                                 

if model.extraStateOutput
    isActive = model.getActivePhases();

    state.bfactor = zeros(model.G.cells.num, sum(isActive));
    b = {double(bW)};
    state = model.setPhaseData(state, b, 'bfactor');  
    
    state.mob = zeros(model.G.cells.num, sum(isActive));
    mob = {double(mobW)};
    state = model.setPhaseData(state, mob, 'mob');
    
    nInterfaces = size(model.operators.N, 1);
    state.upstreamFlag = false(nInterfaces, sum(isActive));
    mob = {upcw};
    state = model.setPhaseData(state, mob, 'upstreamFlag');
    
    state.rho = zeros(model.G.cells.num, sum(isActive));
    rho = {double(rhoW)};
    state = model.setPhaseData(state, rho, 'rho');    
end
%}                                                 
                                                 
%}
% Apply scaling
eqs{1} = eqs{1}.*(dt./s.pv);
eqs{2} = eqs{2}.*(dt./s.pv);

problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);

state.theta = double(sW).*rock.poro;
end


%{
Copyright 2009-2017 SINTEF ICT, Applied Mathematics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
