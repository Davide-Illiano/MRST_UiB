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
[p_prev, c_prev, wellSol] = model.getProps(state, 'pressure', 'concentration','wellsol');
% Properties at previous timestep
[p0, c0,wellSol0] = model.getProps(state0, 'pressure', 'concentration', 'wellSol');

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

dp = s.Grad(p) - gdz;
upc  = (double(dp)<=0);

% Compute transmissibility
theta = fluid.theta;

sW = theta(p_prev,c_prev);
sW0 = theta(p0,c0);


% p_face = s.faceAvg(p);
Kmult = fluid.Kmult(p_prev, sW);
T = s.T.*s.faceUpstr(upc, Kmult); %    faceUpstr(upc, Kmult);


vW = -T.*dp;

bWvW = vW;

varK=0.1;
K_MEAN=15;
Jslope=0.05;
cvel=exp(varK/2)/(Jslope*K_MEAN); 

VcW = s.faceUpstr(upc, c).*cvel.*bWvW;
VcD = -1e-2.*s.Grad(c);

Vc = VcW + VcD;

% Conservation of mass for water
V = model.G.cells.volumes;
water = (V./dt).*( theta(p_prev,c_prev) - theta(p0,c0)  + model.L_p.*(p - p_prev)) ...
    + s.Div(bWvW);


% Transport equation
transport = V./dt.*( theta(p_prev,c_prev).*c - theta(p0,c0).*c0 +  model.L_c.*(c - c_prev) ) ...
     + s.Div(Vc);

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
    state1 = model.storeFluxes(state, vW, [], []);    
                                                    %% not vW you need to use qw!!!!
                                                    %% the boundary flux is canceled out be careful!
    v = find(state.flux);
    state1.flux(v) = state.flux(v);
    state.flux = state1.flux;
                                                    
%     X = model.G.cartDims(1);
%     q = state.flux(size(state.flux,1)-2*X);%mean(state.flux(floor(size(state.flux,1)/2)+2*X:size(state.flux,1)-2*X));   % second half is vertical velocity which is now constnat
%     q = ones(2*X,1).*q;
%     q = [-1*q(1:X)' q(X+1:2*X)']';      % first half must be negative sign
%     state = model.storeBoundaryFluxes(state, q, [], [], drivingForces);
%     v = faceFlux2cellVelocity(model.G,sum(state.flux,2));
%      state.flux(size(state.flux,1)-X+1:size(state.flux,1)) = state.flux(size(state.flux,1)-X)*ones(X,1);
%     mean(v(:,2))
%             act = model.getActivePhases();
%             tmp = cell(numel(act), 1);
%             tmp(act) = bnd_cond.phaseVolume;
%             state = model.storeBoundaryFluxes(state, tmp{:}, forces);
end


state.theta = double(theta(p_prev,c_prev));
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
