function [problem, state] = equationsRichardsTransportLScheme(state0, state, model, dt, drivingForces, varargin)
opt = struct('Verbose', mrstVerbose, ...
             'reverseMode', false,...
             'resOnly', false,...
             'iteration', -1);

opt = merge_options(opt, varargin{:});

W = drivingForces.W;
s = model.operators;
fluid = model.fluid;

% Properties at current timestep
[p_prev, c_prev, wellSol] = model.getProps(state, 'pressure', 'concentration', 'wellsol');
% Properties at previous timestep
[p0, c0, wellSol0] = model.getProps(state0, 'pressure', 'concentration', 'wellSol');

[wellVars, wellVarNames, wellMap] = model.FacilityModel.getAllPrimaryVariables(wellSol);

p = p_prev;
c = c_prev;
% Initialize independent variables.
if ~opt.resOnly
    % ADI variables needed since we are not only computing residuals.
    if ~opt.reverseMode
        [p, c, wellVars{:}] = initVariablesADI(p, c, wellVars{:});
    else
        wellVars0 = model.FacilityModel.getAllPrimaryVariables(wellSol0);
        [p0, c0, wellVars0{:}] = initVariablesADI(p0, c0, wellVars0{:}); %#ok
    end
end
primaryVars = {'pressure','concentration', wellVarNames{:}};

% Gravity contribution
gdz = model.getGravityGradient();


bW     = fluid.bW(p_prev);
bW0     = fluid.bW(p0);

rhoW   = bW.*fluid.rhoWS;
rhoWf  = s.faceAvg(rhoW);

dp = s.Grad(p) - rhoWf.*gdz;
upc  = (double(dp)<=0);


% Compute transmissibility
theta = fluid.theta;

sW = theta(p_prev,c_prev)./model.rock.poro;
sW0 = theta(p0,c0)./model.rock.poro;

% Evaluate water properties
%[vW, bW, mobW, rhoW, pW, upcw] = getFluxAndPropsWater_BO(model, p_prev, sW, krW, T, gdz);
%mobW = 0*mobW + 1;


% p_face = s.faceAvg(p);
Kmult = fluid.Kmult(p_prev, c_prev);
T = s.T.*s.faceUpstr(upc, Kmult);
vW = -T.*dp;
bWvW = s.faceUpstr(upc, bW).*vW;


dpc = s.Grad(p0) - rhoWf.*gdz;
Kmultc = fluid.Kmult(p0, c0);
Tc = s.T.*s.faceUpstr(upc, Kmultc);
vWc = -Tc.*dpc;
bWvWc = s.faceUpstr(upc, bW).*vWc;

VcW = s.faceUpstr(upc,c).*bWvWc;
VcD = -s.Grad(c).*s.faceAvg(sW);

Vc = VcW + VcD;


% Conservation of mass for water
V = model.G.cells.volumes;
water = (V./dt).*( theta(p_prev,c_prev) - theta(p0,c0) + 5e-7.*(p - p_prev)) + s.Div(bWvW);

% Transport equation
transport = (V./dt).* (theta(p_prev,c_prev).*c - theta(p0,c0).* c0 + 5e-7.*(c - c_prev)) + s.Div(Vc)  ;

eqs = {water,transport};
names = {'water', 'concentration'};
types = {'cell', 'cell'};

% Add in any fluxes / source terms prescribed as boundary conditions.
rho = {rhoW};
mob = {sW};
sat = {sW};

[eqs, state] = addBoundaryConditionsAndSources(model, eqs, names, types, state, ...
                                                                 {p}, sat, mob, rho, ...
                                                                 {}, {c}, ...
                                                                 drivingForces);
% Finally, add in and setup well equations
[eqs, names, types, state.wellSol] = model.insertWellEquations(eqs, names, types, wellSol0, wellSol, wellVars, wellMap, p, mob, rho, {}, {c}, dt, opt);
% Apply scaling
eqs{1} = eqs{1}.*(dt./V);
eqs{2} = eqs{2}.*(dt./V);

problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);

if model.outputFluxes
    state = model.storeFluxes(state, vW, [], []);
end
%{
if model.extraStateOutput
    isActive = model.getActivePhases();

    state.bfactor = zeros(model.G.cells.num, sum(isActive));
    b = {double(bW)};
    state = model.setPhaseData(state, b, 'bfactor');
    
    
    state.mob = zeros(model.G.cells.num, sum(isActive));
    mob = {double(mob)};
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
state.theta = double(sW).*model.rock.poro;
end

%{
Copyright 2009-2017 SINTEF ICT, Applied Mathematics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any latezr version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
