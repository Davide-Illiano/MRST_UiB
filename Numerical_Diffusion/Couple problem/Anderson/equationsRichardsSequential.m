function [problem, state] = equationsRichardsSequential(state0, state, model, dt, drivingForces, varargin)
opt = struct('Verbose', mrstVerbose, ...
             'reverseMode', false,...
             'resOnly', false,...
             'iteration', -1);

opt = merge_options(opt, varargin{:});

W = drivingForces.W;
s = model.operators;
fluid = model.fluid;

% Properties at current timestep
[p, c, wellSol] = model.getProps(state, 'pressure', 'concentration','wellsol');
% Properties at previous timestep
[p0, c0, wellSol0] = model.getProps(state0, 'pressure', 'concentration', 'wellSol');

[wellVars, wellVarNames, wellMap] = model.FacilityModel.getAllPrimaryVariables(wellSol);

% Initialize independent variables.
if ~opt.resOnly
    % ADI variables needed since we are not only computing residuals.
    if ~opt.reverseMode
        [p, wellVars{:}] = initVariablesADI(p, wellVars{:});
    else
        wellVars0 = model.FacilityModel.getAllPrimaryVariables(wellSol0);
        [p0,wellVars0{:}] = initVariablesADI(p0, wellVars0{:}); %#ok
    end
end
primaryVars = {'pressure', wellVarNames{:}};

rhoW   = fluid.rhoWS;
rhoWf  = s.faceAvg(rhoW);

ez = model.operators.Grad(model.G.cells.centroids) * [0 1]';
dp = s.Grad(p) +  ez; 
upc  = (double(dp)<=0);

% Compute transmissibility
theta = fluid.theta;
rock = model.rock;
sW = theta(p,c)./rock.poro;
sW0 = theta(p0,c0)./rock.poro;


% p_face = s.faceAvg(p);
Kmult = fluid.Kmult(p,theta(p,c));  
T = s.T.*s.faceUpstr(upc, Kmult);

vW = -T.*dp;


% Conservation of mass for water
V = model.G.cells.volumes;

water = (V./dt).*( theta(p,c) - theta(p0,c0) ) + s.Div(vW); %(V./dt).*( bW.*theta(p) - bW0.*theta(p0) ) + s.Div(bWvW); %

eqs = {water};
names = {'water'};
types = {'cell'};

% Add in any fluxes / source terms prescribed as boundary conditions.
rho = {rhoW};
mob = {sW};
sat = {sW};

[eqs, state] = addBoundaryConditionsAndSources(model, eqs, names, types, state, ...
                                                                 {p}, sat, mob, rho, ...
                                                                 {}, {}, ...
                                                                 drivingForces);
% Finally, add in and setup well equations
[eqs, names, types, state.wellSol] = model.insertWellEquations( eqs, names, ...
                                                     types, wellSol0, wellSol, ...
                                                     wellVars, wellMap, ...
                                                     p, mob, rho, ...
                                                     {}, {}, ...
                                                     dt, opt);
% Apply scaling
eqs{1} = eqs{1}.*(dt./V);  

problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);

if model.outputFluxes
    state = model.storeFluxes(state, vW, [], []);
end


state.theta = double(theta(p,c));
%state.thetafirst = fluid.thetafirst(p);
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
