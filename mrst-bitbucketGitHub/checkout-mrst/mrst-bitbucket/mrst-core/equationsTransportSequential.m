function [problem, state] = equationsTransportSequential(state0, state, model, dt, drivingForces, varargin)
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

% p = model.p_flow;
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

% Gravity contribution
gdz = model.getGravityGradient();


bW     = fluid.bW(p);

rhoW   = bW.*fluid.rhoWS;
rhoWf  = s.faceAvg(rhoW);

dp = s.Grad(p) - rhoWf.*gdz;
upc  = (double(dp)<=0);

Kmult = model.K(:); %fluid.Kmult(p, theta(p,c));   
T = s.T.*s.faceUpstr(upc, Kmult);

vW = -T.*dp;

VcW = s.faceUpstr(upc,c).*vW;  
vW = vW .* 0;
vW(1:size(vW,1)/2) = 0.075/100;

% VcW = s.faceUpstr(upc,c).*vW;   
% VcD = -model.D.*s.Grad(c);
% 
% Vc = VcW + VcD;

VcW = s.faceAvg(c).*vW;   
VcD = -model.D.*s.Grad(c);

Vc = VcW + VcD;

% Conservation of mass for water
V = model.G.cells.volumes;

sW = 1;

transport = (V./dt).*(c - c0) + s.Div(Vc);


eqs = { transport};
names = {'concentration'};
types = {'cell'};

% Add in any fluxes / source terms prescribed as boundary conditions.
rho = {rhoW};
mob = {sW};
sat = {sW};

[eqs, state] = addBoundaryConditionsAndSources(model, eqs, names, types, state, ...
                                                                 {p}, sat, mob, rho, ...
                                                                 {}, {c}, ...
                                                                 drivingForces);
                                                             
% Finally, add in and setup well equations
[eqs, names, types, state.wellSol] = model.insertWellEquations( eqs, names, ...
                                                     types, wellSol0, wellSol, ...
                                                     wellVars, wellMap, ...
                                                     p, mob, rho, ...
                                                     {}, {c}, ...
                                                     dt, opt);
% Apply scaling
eqs{1} = eqs{1}.*(dt./V);

problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);

end
