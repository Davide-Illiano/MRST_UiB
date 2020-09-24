function [problem, state] = equationsRichards(state0, state, model, dt, drivingForces, varargin)
opt = struct('Verbose', mrstVerbose, ...
             'reverseMode', false,...
             'resOnly', false,...
             'iteration', -1);

opt = merge_options(opt, varargin{:});

W = drivingForces.W;
s = model.operators;
fluid = model.fluid;

% Properties at current timestep
[p, sW, wellSol] = model.getProps(state, 'pressure', 'water', 'wellsol');
% Properties at previous timestep
[p0, sW0, wellSol0] = model.getProps(state0, 'pressure', 'water', 'wellSol');

[wellVars, wellVarNames, wellMap] = model.FacilityModel.getAllPrimaryVariables(wellSol);

% Initialize independent variables.
if ~opt.resOnly
    % ADI variables needed since we are not only computing residuals.
    if ~opt.reverseMode
        [p, sW, wellVars{:}] = initVariablesADI(p, sW, wellVars{:});
    else
        wellVars0 = model.FacilityModel.getAllPrimaryVariables(wellSol0);
        [p0, sW0 ,wellVars0{:}] = initVariablesADI(p0, sW0, wellVars0{:}); %#ok
    end
end
primaryVars = {'pressure','sW',  wellVarNames{:}};

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
sW = theta(p)./rock.poro;
sW0 = theta(p0)./rock.poro;
sO = 1 - sW;
sO0 = 1 - sW0;

% Evaluate relative permeabilities
[krW] = model.evaluateRelPerm({sW});

% Multipliers for properties
[pvMult, transMult, mobMult, pvMult0] = getMultipliers(model.fluid, p, p0);

% Modifiy relperm by mobility multiplier (if any)
krW = mobMult.*krW; 

Kmult = fluid.Kmult(p, theta(p));   
T = s.T.*s.faceUpstr(upc, Kmult);

vW = -T.*dp;


% Evaluate water properties
[vW, bW, mobW, rhoW, pW, upcw] = getFluxAndPropsWater_BO(model, p, sW, krW, T, gdz);
bW0 = model.fluid.bW(p0);

bWvW = s.faceUpstr(upc, bW).*vW;

% Conservation of mass for water
V = model.G.cells.volumes;
water = (V./dt).*( bW.*theta(p) - bW0.*theta(p0) ) + s.Div(bWvW); %(V./dt).*( bW.*theta(p) - bW0.*theta(p0) ) + s.Div(bWvW); %

eqs = {water};
names = {'water'};
types = {'cell'};

% Add in any fluxes / source terms prescribed as boundary conditions.
rho = {rhoW};
mob = {mobW};
sat = {sW};
pressures = {p};

[eqs, state] = addBoundaryConditionsAndSources(model, eqs, names, types, state, ...
                                                                 pressures, sat, mob, rho, ...
                                                                 {}, {}, ...
                                                                 drivingForces);
% Finally, add in and setup well equations
[eqs, names, types, state.wellSol] = model.insertWellEquations(eqs, names, types, wellSol0, wellSol, wellVars, wellMap, p, mob, rho, {}, {}, dt, opt);
% Apply scaling
eqs{1} = eqs{1}.*(dt./V); 

problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);

if model.outputFluxes
    %state = model.storeFluxes(state, vW, [], []);
end


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
