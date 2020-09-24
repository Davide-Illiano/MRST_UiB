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
[p0, c0, wellSol0] = model.getProps(state0, 'pressure','concentration', 'wellSol');

[wellVars, wellVarNames, wellMap] = model.FacilityModel.getAllPrimaryVariables(wellSol);

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
primaryVars = {'pressure', 'concentration',  wellVarNames{:}};
x = model.G.cells.centroids(:,1);
y = model.G.cells.centroids(:,2);


% Gravity contribution
gdz = model.getGravityGradient();

bW     = fluid.bW(p);
bW0     = fluid.bW(p0);

rhoW   = bW.*fluid.rhoWS;
rhoWf  = s.faceAvg(rhoW);

dp = s.Grad(p) - rhoWf.*gdz;  
upc  = (double(dp)<=0);

% Compute transmissibility
theta = fluid.theta;
rock = model.rock;
sW = theta(p,c)./rock.poro;
sW0 = theta(p0,c0)./rock.poro;

% p_face = s.faceAvg(p);
Kmult = fluid.Kmult(p, sW);
T = s.T.*s.faceUpstr(upc, Kmult);
model.rock.perm = Kmult;

vW = -T.*dp;

bWvW = s.faceUpstr(upc, bW).*vW;

VcW = s.faceUpstr(upc,c).*bWvW;
VcD = -s.Grad(c).*s.faceAvg(sW);

Vc = VcW + VcD;

% Conservation of mass for water
V = model.G.cells.volumes;

water = (V./dt).*( theta(p,c) - theta(p0,c0)) + s.Div(vW);   % s.Div(bWvW);
transport = (V./dt).*(theta(p,c).*c - theta(p0,c0).* c0)  - s.Div(Vc);% + c.*R;

eqs = {water, transport};
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
%{
[eqs, names, types, state.wellSol] = model.insertWellEquations( eqs, names, ...
                                                     types, wellSol0, wellSol, ...
                                                     wellVars, wellMap, ...
                                                     p, mob, rho, ...
                                                     {}, {c}, ...
                                                     dt, opt);
%}
% Apply scaling
eqs{1} = eqs{1}.*(dt./V);
eqs{2} = eqs{2}.*(dt./V);

problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);

if model.outputFluxes
    state = model.storeFluxes(state, vW, [], []);
end

%state.theta = double(sW.*rock.poro);
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
