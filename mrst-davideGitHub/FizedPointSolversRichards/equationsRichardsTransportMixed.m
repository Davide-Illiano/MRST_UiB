%% We want a mixed scheme, we will use few l scheme iterations and then we move to newton.
function [problem, state] = equationsRichardsTransportMixed(state0, state, model, dt, drivingForces, varargin)
opt = struct('Verbose', mrstVerbose, ...
             'reverseMode', false,...
             'resOnly', false,...
             'iteration', -1);
opt = merge_options(opt, varargin{:});

%% L scheme
if opt.iteration < model.n
    


W = drivingForces.W;
s = model.operators;
fluid = model.fluid;
% Properties at current timestep
[p_prev, c_prev, theta_prev] = model.getProps(state, 'pressure', 'concentration', 'theta');
% Properties at previous timestep
[p0, c0, theta0] = model.getProps(state0, 'pressure','concentration', 'theta');

%[wellVars, wellVarNames, wellMap] = model.FacilityModel.getAllPrimaryVariables(wellSol);

% Initialize independent variables.
p = p_prev;
c = c_prev;
theta = theta_prev;

if ~opt.resOnly
    % ADI variables needed since we are not only computing residuals.
    if ~opt.reverseMode
        [p, c, theta] = initVariablesADI(p, c,theta);
    else
        %wellVars0 = model.FacilityModel.getAllPrimaryVariables(wellSol0);
        [p0, c0, theta0] = initVariablesADI(p0, c0, theta); %#ok
    end
end
primaryVars = {'pressure', 'concentration', 'theta'};

% Gravity contribution
gdz = model.getGravityGradient();

bW     = fluid.bW(p_prev);
bW0     = fluid.bW(p0);

rhoW   = bW.*fluid.rhoWS;
rhoWf  = s.faceAvg(rhoW);

dp = s.Grad(p) - rhoWf.*gdz;  
upc  = (double(dp)<=0);

% Compute transmissibility
rock = model.rock;
%sW = theta(p,c)./rock.poro;
%sW0 = theta(p0,c0)./rock.poro;


% p_face = s.faceAvg(p);
Kmult = fluid.Kmult(theta_prev);
T = s.T.*s.faceUpstr(upc, Kmult);


vW = -T.*dp;

bWvW = s.faceUpstr(upc, bW).*vW;
VcW = s.faceUpstr(upc,c).*bWvW;

VcD = -s.Grad(c); %-s.Grad(c).*s.faceAvg(theta);

Vc = VcW - VcD;

% Conservation of mass for water
V = model.G.cells.volumes;

dTheta = theta - theta0;
pos = dTheta >= fluid.delta;
neg = dTheta <= -fluid.delta;
zero = abs(dTheta) < fluid.delta;

%{
F = 1/fluid.tau.*(p + fluid.p_cap(theta_prev) - fluid.gamma) .* pos + ...
     1/fluid.tau.*(p + fluid.p_cap(theta_prev) + fluid.gamma) .* neg + ...
     fluid.delta/(fluid.tau*fluid.delta + fluid.gamma).* (p + fluid.p_cap(theta_prev)) .* zero;

F_first = full(diag(F.jac{1}));
F_second = full(diag(F.jac{3}));
one = sign(F_first)
two = sign(F_second);
%}
%{
water = (V).*( 1/fluid.tau.*(p + fluid.p_cap(theta_prev) - fluid.gamma) .* pos + ...
     1/fluid.tau.*(p + fluid.p_cap(theta_prev) + fluid.gamma) .* neg + ...
     fluid.delta/(fluid.tau*fluid.delta + fluid.gamma).* (p + fluid.p_cap(theta_prev)) .* zero )  + ...
     s.Div(vW) + model.L_p.*(p - p_prev) ;  
%}
%
water = (V./dt).*( theta - theta0 )+ model.L_p.*(p - p_prev) + ...
      s.Div(vW)  ;  
%}
transport = (V./dt).*(theta.*c - theta0.*c0) + s.Div(Vc) + V.*c./(c_prev + 1) +...
    model.L_c.*(c - c_prev)  ;
%{
three =  p + fluid.p_cap(theta_prev) - fluid.gamma .* fluid.getGamma(theta_prev, theta0)...
         - fluid.tau .*1/dt.*(theta - theta0) + model.L_t.*(theta - theta_prev) ;
%}
%
three = (fluid.tau + dt.*model.L_t + fluid.gamma./fluid.delta).*theta - (fluid.tau + fluid.gamma./fluid.tau.*zero).*theta0 - dt.*...
    (p + fluid.p_cap(theta_prev) + model.L_t.*theta_prev + fluid.gamma.*neg - fluid.gamma.*pos);
%}
eqs = {water, transport, three};
names = {'water', 'concentration', 'theta'};
types = {'cell', 'cell', 'cell'};

% Add in any fluxes / source terms prescribed as boundary conditions.
rho = {rhoW};
mob = {theta./rock.poro};
sat = {theta./rock.poro};

[eqs, state] = addBoundaryConditionsAndSources(model, eqs, names, types, state, ...
                                                                 {p}, sat, mob, rho, ...
                                                                 {}, {c,theta}, ...
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
eqs{3} = eqs{3}.*(dt);

problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);

if model.outputFluxes
    state = model.storeFluxes(state, vW, [], []);
end
state.p_prev = p_prev;

%% Newton method
else

    
W = drivingForces.W;
s = model.operators;
fluid = model.fluid;

% Properties at current timestep
[p, c, theta] = model.getProps(state, 'pressure', 'concentration', 'theta');
% Properties at previous timestep
[p0, c0, theta0] = model.getProps(state0, 'pressure','concentration', 'theta');

%[wellVars, wellVarNames, wellMap] = model.FacilityModel.getAllPrimaryVariables(wellSol);

% Initialize independent variables.

if ~opt.resOnly
    % ADI variables needed since we are not only computing residuals.
    if ~opt.reverseMode
        [p, c, theta ] = initVariablesADI(p, c,theta);
    else
        %wellVars0 = model.FacilityModel.getAllPrimaryVariables(wellSol0);
        [p0, c0, theta0] = initVariablesADI(p0, c0, theta); %#ok
    end
end
primaryVars = {'pressure', 'concentration', 'theta'};

% Gravity contribution
gdz = model.getGravityGradient();

bW     = fluid.bW(p);
bW0     = fluid.bW(p0);

rhoW   = bW.*fluid.rhoWS;
rhoWf  = s.faceAvg(rhoW);

dp = s.Grad(p) - rhoWf.*gdz;  
upc  = (double(dp)<=0);

% Compute transmissibility
rock = model.rock;

% p_face = s.faceAvg(p);
Kmult = fluid.Kmult(theta);
T = s.T.*s.faceUpstr(upc, Kmult);


vW = -T.*dp;

bWvW = s.faceUpstr(upc, bW).*vW;
VcW = s.faceUpstr(upc,c).*bWvW;

VcD = -s.Grad(c); %-s.Grad(c).*s.faceAvg(theta);

Vc = VcW - VcD;

% Conservation of mass for water
V = model.G.cells.volumes;

dTheta = theta - theta0;
pos = dTheta./dt >= fluid.delta;
neg = dTheta./dt <= -fluid.delta;
zero = abs(dTheta./dt) < fluid.delta;

%{
water = (V).*( 1/fluid.tau.*(p + fluid.p_cap(theta) - fluid.gamma) .* pos + ...
     1/fluid.tau.*(p + fluid.p_cap(theta) + fluid.gamma) .* neg + ...
     fluid.delta/(fluid.tau*fluid.delta + fluid.gamma).* (p + fluid.p_cap(theta)) .* zero )  + ...
     s.Div(vW);  
%}
%
water = (V./dt).*( theta- theta0 )  + ...
     s.Div(vW);
%} 
transport = (V./dt).*(theta.*c - theta0.*c0 ) + s.Div(Vc) + V.*c./(c+1); 

three =  p + fluid.p_cap(theta) - fluid.gamma.* fluid.getGamma(theta,theta0)  - fluid.tau.*1/dt.*(theta - theta0);


%noWater = double(theta./rock.poro) == 0;
%transport(noWater) = c(noWater);

eqs = {water, transport, three};
names = {'water', 'concentration', 'theta'};
types = {'cell', 'cell', 'cell'};

% Add in any fluxes / source terms prescribed as boundary conditions.
rho = {rhoW};
mob = {theta./rock.poro};
sat = {theta./rock.poro};

[eqs, state] = addBoundaryConditionsAndSources(model, eqs, names, types, state, ...
                                                                 {p}, sat, mob, rho, ...
                                                                 {}, {c,theta}, ...
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
eqs{3} = eqs{3}.*(dt);

problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);

if model.outputFluxes
    state = model.storeFluxes(state, vW, [], []);
end
end


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
