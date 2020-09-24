function obj = CO2props(varargin); %rhofile, hfile, noassert)

  opt.rhofile              = 'rho_extensive';
  opt.mufile               = 'mu_extensive';
  opt.hfile                = 'h_small';
  opt.assert               = true;
  opt.sharp_phase_boundary = true;
  opt = merge_options(opt, varargin{:});
    
%% ESTABLISHING STATIC DATA
    
  rho = load(opt.rhofile); % load density data
  h   = load(opt.hfile);   % load enthalpy data
  mu  = load(opt.mufile);  % load viscosity data
  assert_valid_range = opt.assert; % do optional checking whether P and T are
                                   % within covered range (if not checked,
                                   % NaN values will be returned in such
                                   % cases instead)

  % make separate matrices by phase, and compute derivative tables
  rho = establish_data_tables(rho, opt.sharp_phase_boundary);
  h   = establish_data_tables(h  , opt.sharp_phase_boundary);
  mu  = establish_data_tables(mu , opt.sharp_phase_boundary);
  
      
  %% DECLARING MEMBER FUNCTIONS

  % return CO2 phase at (P, T) (0: supercritical; 1: liquid; 2: gas)
  obj.phaseOf     = @(P, T) phase_of(P, T, opt.sharp_phase_boundary);

  % CO2 density and its derivatives
  obj.rhoDPPT = @(P, T) calcMulti(P, T, @(P, T) calcDPPT(rho, P, T), @noder,       @noder);
  obj.rhoDPTT = @(P, T) calcMulti(P, T, @(P, T) calcDPTT(rho, P, T), @noder,       @noder);
  obj.rhoDTTT = @(P, T) calcMulti(P, T, @(P, T) calcDTTT(rho, P, T), @noder,       @noder);
  obj.rhoDPPP = @(P, T) calcMulti(P, T, @(P, T) calcDPPP(rho, P, T), @noder,       @noder);
  obj.rhoDPT  = @(P, T) calcMulti(P, T, @(P, T) calcDPT (rho, P, T), @obj.rhoDPPT, @obj.rhoDPTT);
  obj.rhoDTT  = @(P, T) calcMulti(P, T, @(P, T) calcDTT (rho, P, T), @obj.rhoDPTT, @obj.rhoDTTT);
  obj.rhoDPP  = @(P, T) calcMulti(P, T, @(P, T) calcDPP (rho, P, T), @obj.rhoDPPP, @obj.rhoDPPT);
  obj.rhoDT   = @(P, T) calcMulti(P, T, @(P, T) calcDT  (rho, P, T), @obj.rhoDPT,  @obj.rhoDTT);
  obj.rhoDP   = @(P, T) calcMulti(P, T, @(P, T) calcDP  (rho, P, T), @obj.rhoDPP,  @obj.rhoDPT);
  obj.rho     = @(P, T) calcMulti(P, T, @(P, T) calcVal (rho, P, T), @obj.rhoDP,   @obj.rhoDT);
  obj.rhoInfo = @() printInfo(rho);
  
  
  % CO2 enthalpy and its derivatives
  obj.h     = @(P, T) calcMulti(P, T, @(P, T) calcVal  (h, P, T), @obj.hDP,   @obj.hDT);
  obj.hDP   = @(P, T) calcMulti(P, T, @(P, T) calcDP   (h, P, T), @obj.hDPP,  @obj.hDPT);
  obj.hDPP  = @(P, T) calcMulti(P, T, @(P, T) calcDPP  (h, P, T), @obj.hDPPP, @noder);
  obj.hDPPP = @(P, T) calcMulti(P, T, @(P, T) calcDPPP (h, P, T), @noder,     @noder);
  obj.hDT   = @(P, T) calcMulti(P, T, @(P, T) calcDT   (h, P, T), @obj.hDPT,  @obj.hDTT);
  obj.hDTT  = @(P, T) calcMulti(P, T, @(P, T) calcDTT  (h, P, T), @noder,     @obj.hDTTT);
  obj.hDTTT = @(P, T) calcMulti(P, T, @(P, T) calcDTTT (h, P, T), @noder,     @noder);
  obj.hDPT  = @(P, T) calcMulti(P, T, @(P, T) calcDPT  (h, P, T), @noder,     @noder);
  obj.hInfo = @() printInfo(h);
  
  % CO2 dynamic viscosity and its derivatives
  obj.muDPT  = @(P, T) calcMulti(P, T, @(P, T) calcDPT  (mu, P, T), @noder,     @noder);
  obj.muDTTT = @(P, T) calcMulti(P, T, @(P, T) calcDTTT (mu, P, T), @noder,     @noder);
  obj.muDPPP = @(P, T) calcMulti(P, T, @(P, T) calcDPPP (mu, P, T), @noder,     @noder);
  obj.muDTT  = @(P, T) calcMulti(P, T, @(P, T) calcDTT  (mu, P, T), @noder,     @obj.muDTTT);
  obj.muDPP  = @(P, T) calcMulti(P, T, @(P, T) calcDPP  (mu, P, T), @obj.muDPPP, @noder);
  obj.muDT   = @(P, T) calcMulti(P, T, @(P, T) calcDT   (mu, P, T), @obj.muDPT,  @obj.muDTT);
  obj.muDP   = @(P, T) calcMulti(P, T, @(P, T) calcDP   (mu, P, T), @obj.muDPP,  @obj.muDPT);
  obj.mu     = @(P, T) calcMulti(P, T, @(P, T) calcVal  (mu, P, T), @obj.muDP,   @obj.muDT);
  obj.muInfo = @() printInfo(mu);
  
  % compressibility factors
  obj.beta  = @(P, T)   obj.rhoDP(P, T) ./ obj.rho(P, T);  % Compressibility coef.
  obj.gamma = @(P, T) - obj.rhoDT(P, T) ./ obj.rho(P, T);  % Coef. of thermal expansion
  
  % compressibility factor cross derivative (equal for both Pbeta and Tbeta)
  obj.betaDPDT = @(P,T) (obj.rhoDPT(P, T) + obj.rhoDP(P, T) .* obj.rhoDT(P, T)) ./ obj.rho(P, T);
  
  obj.dispose  = @dispose;        % Destructor - should always be called when
                                  % object has reached end of its life cycle

  
%% MEMBER METHODS IMPLEMENTATION

function m = printInfo(m)
    m.P
    m.T
    size(m.vals)
end

function noder(P, T)
% Function used as placeholder where no derivative function is available/implemented.
% Only used to throw an error
    error('The requested derivative is not implemented')
end

function res = calcVal(m, P, T)
    res = extractValues(P, T, m.P.span, m.T.span, m.vals, m.vals, m.vals);
end

function res = calcDP(m, P, T)
% First partial derivative in pressure
    pspan = shrink_span(m.P.span, m.P.stepsize);
    res = extractValues(P, T, pspan, m.T.span, m.hyp.dp, m.liq.dp, m.gas.dp);
end

function res = calcDT(m, P, T)
% First partial derivative in temperature
    tspan = shrink_span(m.T.span, m.T.stepsize);
    res = extractValues(P, T, m.P.span, tspan, m.hyp.dt, m.liq.dt, m.gas.dt);
end

function res = calcDPP(m, P, T)
% second partial derivative in pressure
    pspan = shrink_span(m.P.span, 2 * m.P.stepsize);
    res = extractValues(P, T, pspan, m.T.span, m.hyp.dpp, m.liq.dpp, m.gas.dpp);
end

function res = calcDTT(m, P, T)
% second partial derivative in temperature
    tspan = shrink_span(m.T.span, 2 * m.T.stepsize);
    res = extractValues(P, T, m.P.span, tspan, m.hyp.dtt, m.liq.dtt, m.gas.dtt);
end

function res = calcDPT(m, P, T)
% Cross derivative in pressure and temperature
    pspan = shrink_span(m.P.span, m.P.stepsize);
    tspan = shrink_span(m.T.span, m.T.stepsize);
    res = extractValues(P, T, pspan, tspan, m.hyp.dpt, m.liq.dpt, m.gas.dpt);
end

function res = calcDPPT(m, P, T)
% Mixed derivative; twice in pressure, once in temperature
    pspan = shrink_span(m.P.span, 2 * m.P.stepsize);
    tspan = shrink_span(m.T.span, 1 * m.T.stepsize);
    res = extractValues(P, T, pspan, tspan, m.hyp.dppt, m.liq.dppt, m.gas.dppt);
end

function res = calcDPTT(m, P, T)
% Mixed derivative; once in pressure, twice in temperature
    pspan = shrink_span(m.P.span, 1 * m.P.stepsize);
    tspan = shrink_span(m.T.span, 2 * m.T.stepsize);
    res = extractValues(P, T, pspan, tspan, m.hyp.dptt, m.liq.dptt, m.gas.dptt);
end

function res = calcDPPP(m, P, T)
% third partial derivative in pressure
    pspan = shrink_span(m.P.span, 3 * m.P.stepsize);
    res = extractValues(P, T, pspan, m.T.span, m.hyp.dppp, m.liq.dppp, m.gas.dppp);
end

function res = calcDTTT(m, P, T)
% third partial derivative in temperature
    tspan = shrink_span(m.T.span, 3 * m.T.stepsize);
    res = extractValues(P, T, m.P.span, tspan, m.hyp.dttt, m.liq.dttt, m.gas.dttt);
end

% @@ THIS FUNCTION TAKEN FROM ADI.M, WHERE IT WAS A PRIVATE FUNCTION.  SO
% IT IS A BIT OF A HACK TO PASTE IT HERE.
function J = lMultDiag(d, J1)
    n = numel(d);
    D = sparse((1:n)', (1:n)', d, n, n);
    J = cell(1, numel(J1));
    for k = 1:numel(J)
        J{k} = D*J1{k};
    end
end


function res = calcMulti(P, T, fun, fun_dp, fun_dt)
% function used as wrapper around the functions used to compute rho, enthalpy, and 
% their derivatives.  The reason for this wrapper is to add support for the ADI
% framework when applicable.
    
    assert(isa(T,'ADI') | isvector(T));
    assert(isa(P,'ADI') | isvector(P));
    
    if(isa(P,'ADI'))
        pval=P.val;
    else
        pval=P;
    end
    
    if (isa(T,'ADI'))
        tval = T.val;
    else
        tval = T;
        if(numel(tval)==1)
           tval=ones(numel(pval),1)*T; 
        else
           tval = T;  
        end
    end

    res = fun(pval(:), tval(:));
    
    if(isa(P,'ADI'))
        % Check if the function has an available derivative function
        dres_dp = fun_dp(pval, tval);
        if(isa(T,'ADI'))
            dres_dt = fun_dt(pval, tval);
            assert(numel(P.jac)==numel(T.jac));
            
            %res = ADI(res, lMultDiag(dres_dp, P.jac) + lMultDiag(dres_dt, T.jac));
            tmp1 = lMultDiag(dres_dp, P.jac);
            tmp2 = lMultDiag(dres_dt, T.jac);
            res = ADI(res, cell_adder(tmp1, tmp2));%{tmp1{1} + tmp2{1}, tmp1{2} + tmp2{2}});
        else
            res = ADI(res, lMultDiag(dres_dp, P.jac));
       end
    else
        if(isa(T,'ADI'))
            dres_dt = fun_dt(pval, tval);
            res = ADI(res, lMultDiag(dres_dt, T.jac));
        end 
    end
end
  
function c = cell_adder(c1, c2)
    for i = 1:numel(c1)
        c{i} = c1{i} + c2{i};
    end
end

function res = extractValues(P, T, spanP, spanT, supsamples, liqsamples, gassamples)
    
    if assert_valid_range
        assert_in_range(P, spanP);
        assert_in_range(T, spanT);
    else
        [P, T]   = truncate_PT_vectorized(P, T, spanP, spanT);
    end
    
    phase    = phase_of(P, T, opt.sharp_phase_boundary);
    sph      = (phase == 0);
    lph      = (phase == 1);
    gph      = (phase == 2);
    
    res      = nan * ones(numel(P), 1);
    res(sph) = extract_val_vectorized(supsamples, spanP, spanT, P(sph), T(sph));
    res(lph) = extract_val_vectorized(liqsamples, spanP, spanT, P(lph), T(lph));
    res(gph) = extract_val_vectorized(gassamples, spanP, spanT, P(gph), T(gph));
    
end

function dispose()
    clear rho;
    clear h;
    fprintf('Member data deleted.\n');
end

  
end
  
%% HELPER FUNCTIONS IMPLEMENTATION

function span = shrink_span(span, ds)
    % shrink span by amount ds, centered
    span = [span(1) + ds/2, span(2) - ds/2];
end

function tables = establish_data_tables(tables, sharp_phase_boundary)
% Based on the orignal table, generate additional, separate tables for each
% phase and for derivatives in P and T
    vals = tables.vals;
    Psteplen = tables.P.stepsize;
    Tsteplen = tables.T.stepsize;
    
    tables.hyp = deriv_tables(vals, Psteplen, Tsteplen);
            
    if sharp_phase_boundary
        % Establish separate tables for hypercritical, liquid and gas regions
        [liq, gas] = prepare_by_phase(vals, tables.P.span, tables.T.span);

        % compute table of derivatives
        tables.liq = deriv_tables(liq , Psteplen, Tsteplen);
        tables.gas = deriv_tables(gas , Psteplen, Tsteplen);
    else
        % This is the easy way.  A more comprehensive approach would ensure
        % that these tables were not needed at all if we do not enforce the
        % sharp boundary.
        tables.liq = tables.hyp;
        tables.gas = tables.hyp;
    end
end

function tables = deriv_tables(base, delta_p, delta_t)
    
    % computing derivatives in pressure
    tables.dp     = diff(base,        1, 1) / delta_p;
    tables.dpp    = diff(tables.dp,   1, 1) / delta_p;
    tables.dppp   = diff(tables.dpp,  1, 1) / delta_p;
    
    % computing derivatives in temperature
    tables.dt     = diff(base,        1, 2) / delta_t;
    tables.dtt    = diff(tables.dt,   1, 2) / delta_t;
    tables.dttt   = diff(tables.dtt,  1, 2) / delta_t;

    % computing cross derivative
    tables.dpt    = diff(tables.dp,   1, 2) / delta_t;
    
    % computer higher mixed derivatives
    tables.dppt   = diff(tables.dpt,  1, 1) / delta_p;
    tables.dptt   = diff(tables.dpt,  1, 2) / delta_t;
end


function phase = phase_of(P, T, sharp_phase_boundary)
% Determine phase of CO2 at pressure P and temperature T.  (Vectorized)
% Flags are:
% * 0 designate "supercritical"
% * 1 designate "liquid"
% * 2 designate "gas"    
    if ~sharp_phase_boundary
        % We do not consider the sharp boundary, and can therefore make a
        % shortcut by considering everything as a single phase
        phase = zeros(numel(P), 1); % everything considered "supercritical"
        return;
    end
    
    phase = nan * ones(numel(P), 1);
  
    [Pc, Tc] = CO2CriticalPoint();
    
    Tsuper = (T>=Tc);
    Psuper = (P>=Pc);
    
    phase( Tsuper &  Psuper) = 0; % hypercritical
    phase( Tsuper & ~Psuper) = 2; % gas
    phase(~Tsuper &  Psuper) = 1; % liquid
    
    % For the remaining entries, we are neither beyond critical temperature
    % nor critical pressure.  We will have to do a more careful analysis to
    % determine at which side of the vapor-liquid boundary (P,T) is located.
    
    remains = find(isnan(phase) & ~isnan(P) & ~isnan(T));
  
    if numel(remains) > 0
        Tinf = T(remains); % vector of temperature of the remaining cases
        Pvap = CO2VaporPressure(Tinf);
        phase_rem = ones(numel(remains), 1); % set to liquid by default
        phase_rem(P(remains) < Pvap) = 2;    % when pressure is lower than
                                             % vaporation pressure, the phase is gas
        phase(remains) = phase_rem;
    end
end

function [liq, gas] = prepare_by_phase(base, span_p, span_t)
% Prepare separate density matrices for the liquid and gas regions ('rliq', 
%'rgas'). The returned matrices will have the same size as 'base',
% but will only contain numeric values for the samples corresponding to the
% liquid (or gas) phase, as well as a few extrapolated values beyond the
% discontinuity in order to faciltiate estimation of derivatives in this
% region.  

    % This function should only be called if we use a sharp phase boundary,
    % hence the 'true' argument passed along to 'compute_phase' here.
    phase = compute_phase(size(base, 1), size(base, 2), span_p, span_t, true);

    liq = NaN * ones(size(base));
    gas = NaN * ones(size(base));

    lphase = (phase ~= 2); % everything that is not gas (i.e. including supercritical)
    gphase = (phase ~= 1); % everything that is not liquid (i.e. including supercritical)

    % Copying relevant values from the original sample matrix to the new,
    % phase-specific matrices
    liq(lphase(:)) = base(lphase(:));
    gas(gphase(:)) = base(gphase(:));

    % Determining largest column number below critical temperature (i.e. the
    % largest column for which extrapolation of values across discontinuity
    % is needed below)
    maxcol = max(find((cellfun(@(x)any(isnan(x)), num2cell(gas(end,:), 1)))));
    %dc = 1:maxcol+1;
    dc = 1:maxcol+3; %@ Strictly speaking, we should extrapolate along
                     %T-direction in order to get the additional values
                     %across boundary necessary for computing derivative.
                     %Currently, extrapolation is only in T.  Fix when time.

    % Determining row indices along discontinuity (column indices is given by (1:maxcol))
    dr = cellfun(@(x)min(find(~isnan(x)))-1, num2cell(liq(:, dc), 1));

    % Removing an extra row of values on each side of discontinuity, to
    % ensure that no value ended on the 'wrong side' (a bit of a hack  @@)
    drl = dr+1;
    drg = dr-1;
    liq(sub2ind(size(liq), drl, 1:maxcol+3)) = NaN;
    gas(sub2ind(size(gas), dr , 1:maxcol+3)) = NaN;
    
    % Extrapolating values along the boundary
    % @@ Perhaps extrapolation should also be based on temperature derivative?  
    % @@ As of now, the extrapolation takes place along isotherms, using pressure derivatives 
    % @@ only.
    ev = 5; % number of rows of extrapolated values
    for i = 0:ev-1
        % quadratic extrapolation of liquid density across boundary
        liq(sub2ind(size(liq),drl-i, dc)) = 3 * (liq(sub2ind(size(liq), drl-i+1, dc)) - ...
                                                liq(sub2ind(size(liq), drl-i+2, dc))) + ...
                                                liq(sub2ind(size(liq),drl-i+3, dc));
        
        % quadratic extrapolation of gas density across boundary
        gas(sub2ind(size(gas), drg+i+1, dc)) = 3 * (gas(sub2ind(size(gas), drg+i, dc)) - ...
                                                   gas(sub2ind(size(gas), drg+i-1, dc))) + ...
                                                   gas(sub2ind(size(gas), drg+i-2, dc));
    end
end

function phase = compute_phase(rows, cols, span_p, span_t, sharp_phase_boundary)

    Pvals = linspace(span_p(1), span_p(2), rows);
    Tvals = linspace(span_t(1), span_t(2), cols);
    [T, P] = meshgrid(Tvals, Pvals);
    phase = reshape(phase_of(P(:), T(:), sharp_phase_boundary), rows, cols);
end
  
function assert_in_range(val, range)
    assert(all(val >= range(1)))
    assert(all(val <= range(2)))
end

function [P, T] = truncate_PT_vectorized(P, T, Pspan, Tspan)
  
    P(P < Pspan(1) | P > Pspan(2)) = NaN;
    T(T < Tspan(1) | T > Tspan(2)) = NaN;
end
  
  
function [p_ix, t_ix, p, t, nans] = ix_and_local_par(grid, spanP, spanT, P, T)
    
    nans = isnan(P) | isnan(T);
    
    Psteps = size(grid, 1) - 1;
    Tsteps = size(grid, 2) - 1;
    
    dSpanP = diff(spanP);
    dSpanT = diff(spanT);
    Dp = dSpanP/Psteps;
    Dt = dSpanT/Tsteps;
  
    p_ix = 1 + floor(((P(~nans) - spanP(1)) * Psteps)/dSpanP);    
    t_ix = 1 + floor(((T(~nans) - spanT(1)) * Tsteps)/dSpanT);
    
    p_ix = min(p_ix, Psteps);
    t_ix = min(t_ix, Tsteps);
    
    Psampled = spanP(1) + (p_ix - 1) * Dp;
    Tsampled = spanT(1) + (t_ix - 1) * Dt;
    
    p = (P(~nans) - Psampled) / Dp; % should be between 0 and 1
    t = (T(~nans) - Tsampled) / Dt; % should be between 0 and 1
    
end
  
function res = extract_val_vectorized(grid, spanP, spanT, P, T)
% Vectorized version of 'extract_val', where P and T are 
% allowed to be vectors
    [p_ix, t_ix, p, t, nans] = ix_and_local_par(grid, spanP, spanT, P, T);
    
    
    lines = size(grid, 1);
    corners = [grid(p_ix   + lines * (t_ix - 1)), ...
               grid(p_ix+1 + lines * (t_ix - 1)), ...
               grid(p_ix   + lines * (t_ix    )), ...
               grid(p_ix+1 + lines * (t_ix    ))];
  
    weight = [(1-p).*(1-t), p.*(1-t), (1-p).*t, p.*t];
                   
    res = NaN * ones(numel(P), 1);
    res(~nans) = sum(corners .* weight, 2);
    
end
