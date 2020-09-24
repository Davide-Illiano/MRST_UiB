function [problem, state] = equationsTransportSequentialLScheme(state0, state, model, dt, drivingForces, varargin)
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
[p0, c0, wellSol0] = model.getProps(state0, 'pressure','concentration', 'wellSol');

[wellVars, wellVarNames, wellMap] = model.FacilityModel.getAllPrimaryVariables(wellSol);

p = p_prev;
c = c_prev;

% Initialize independent variables.
if ~opt.resOnly
    % ADI variables needed since we are not only computing residuals.
    if ~opt.reverseMode
        [c, wellVars{:}] = initVariablesADI( c, wellVars{:});
    else
        wellVars0 = model.FacilityModel.getAllPrimaryVariables(wellSol0);
        [c0, wellVars0{:}] = initVariablesADI(c0, wellVars0{:}); %#ok
    end
end
primaryVars = { 'concentration',  wellVarNames{:}};

%ez = model.operators.Grad(model.G.cells.centroids) * [0 1]';
ez = model.G.cells.centroids(:,2);
dp = s.Grad(p ) ;%; + ez

% Compute transmissibility
theta = fluid.theta;
rock = model.rock;
sW = theta(p_prev,c_prev)./rock.poro;
sW0 = theta(p0,c0)./rock.poro;


% p_face = s.faceAvg(p);
V = model.G.cells.volumes;
Kmult = fluid.Kmult(p_prev, theta(p_prev,c_prev));  
T = s.T.*s.faceAvg(Kmult);

vW =  -T.*dp ;%vW =  -T.*dp - T ;%.* [zeros(2450, 1)' ones(2450,1)']';%;

VcW1 = s.faceAvg(c).*vW;


VcD = -0.*1e-3.*s.Grad(c);%.*s.faceAvg(sW);


% % Flux for transport
% mobSft = 1.*c;
% upc  = (double(dp)<=0);
% vSft = -s.faceUpstr(upc, mobSft).*T.*dp;
% VcW = s.faceUpstr(upc, 0.*mobSft + 1).*vSft;

Vc = VcW1 + VcD;
    
% Conservation of mass for water
transport = (V./dt).*( theta(p_prev,c_prev).*c - theta(p0,c0).*c0 + model.L_c.*(c - c_prev)) ...
     + s.Div(Vc) -0.*1e-3.*c./(c_prev+1);


eqs = { transport};
names = { 'concentration'};
types = {'cell'};

% Add in any fluxes / source terms prescribed as boundary conditions.
rho = {sW};
mob = {sW};
sat = {sW};

boundary_info.T = s.T.*s.faceAvg(Kmult);
boundary_info.T_all = s.T_all.*mean(Kmult);
boundary_info.vW = state.vw;
boundary_info.p = p;

[eqs, state] = addBoundaryConditionsAndSourcesNEW(model, eqs, names, types, state, ...
                                                                 {p}, sat, mob, rho, ...
                                                                 {}, {c}, ...
                                                                 drivingForces, boundary_info);
                                                             
% [eqs, state] = addBoundaryConditionsAndSources(model, eqs, names, types, state, ...
%                                                                  {p}, sat, mob, rho, ...
%                                                                  {}, {c}, ...
%                                                                  drivingForces);                                                             
% Finally, add in and setup well equations
[eqs, names, types, state.wellSol] = model.insertWellEquations(eqs, names, types, wellSol0, wellSol, wellVars, wellMap, p, mob, rho, {}, {c}, dt, opt);
% Apply scaling
eqs{1} = eqs{1}.*(dt./V);

problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);

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
