function fluid = initCoreyFluidROC(varargin)
%Initialize incompressible two-phase fluid model compatible with reordering
%
% SYNOPSIS:
%   fluid = initCoreyFluidROC('pn1', pv1, ...)
%
% PARAMETERS:
%   'pn'/pv - List of 'key'/value pairs defining specific fluid
%             characteristics.  The following parameters must be defined
%             with one value for each of the two fluid phases:
%               - mu  -- Phase viscosities [muw,  muo] in units of Pa*s.
%               - rho -- Phase densities   [rhow, rhoo] in units of
%                        kilogram/meter^3.
%               - n   -- #reg x 2 array of phase relative permeability
%                        exponents.
%               - sr  -- #reg x 2 array of residual phase saturation.
%               - kwm -- #reg x 2 array of phase relative permeability at
%                        residual saturation.
%
%             Optionally, one may specify the following parameter
%               - reg -- Saturation regions, one value per cell in the
%                        model. NB! This option will only work if the
%                        fluid model is evaluated for each cell in the
%                        model at once.
%
% RETURNS:
%   fluid - Fluid data structure as described in 'fluid_structure'
%           representing the current state of the fluids within the
%           reservoir model.
%
% EXAMPLE:
%   fluid = initCoreyFluid('mu' , [   1,  10]*centi*poise     , ...
%                          'rho', [1014, 859]*kilogram/meter^3, ...
%                          'n'  , [   2,   2]                 , ...
%                          'sr' , [ 0.2, 0.2]                 , ...
%                          'kwm', [   1,   1]);
%
%   s = linspace(0, 1, 1001).'; kr = fluid.relperm(s);
%   plot(s, kr), legend('kr_1', 'kr_2')
%
% SEE ALSO:
%   `fluid_structure`, `initSimpleFluid`, `solveIncompFlow`.

%{
Copyright 2009-2018 SINTEF ICT, Applied Mathematics.

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


   opt = struct('mu', [], 'rho', [], 'n', [], 'sr', [], 'kwm', [],...
                'reg', []);
   opt = merge_options(opt, varargin{:});

   if ~isempty(opt.reg)
      m = 2*max(opt.reg);
      assert((numel(opt.n) >= m) & (numel(opt.sr) >= m) ...
           & (numel(opt.kwm) >= m));
      assert(numel(opt.mu) == 2);
      assert(numel(opt.rho) == 2);
   end

   prop = @(  varargin) properties(opt, varargin{:});
   kr   = @(s,varargin) relperm(s, opt, varargin{:});

   % Check fluid options:
   n = size(opt.n, 1);
   assert(numel(opt.mu)     == 2,...
       'Dimensions mismatch in ''mu'' option');
   assert(numel(opt.rho)    == 2,...
       'Dimensions mismatch in ''rho'' option');
   assert(all(size(opt.sr)  == [n, 2]),...
       'Dimensions mismatch in ''sr'' option');
   assert(all(size(opt.n)   == [n, 2]),...
       'Dimensions mismatch in ''n'' option');
   assert(all(size(opt.kwm) == [n, 2]),...
       'Dimensions mismatch in ''kwm'' option');

   fluidopts = struct ('viscw', opt.mu (:,1), 'visco', opt.mu (:,2), ...
                       'srw',   opt.sr (:,1), 'sro',   opt.sr (:,2), ...
                       'nw',    opt.n  (:,1), 'no',    opt.n  (:,2),  ...
                       'krmw',  opt.kwm(:,1), 'krmo',  opt.kwm(:,2),...
                       'satnum',opt.reg);

   fluid    = struct ('properties', prop             , ...
                      'saturation', @(x,varargin) x.s, ...
                      'relperm'   , kr               , ...
                      'param'     , fluidopts);
end

%--------------------------------------------------------------------------
% Helpers follow.
%--------------------------------------------------------------------------

function varargout = properties(opt, varargin)
   varargout{1}                 = opt.mu ;
   if nargout > 1, varargout{2} = opt.rho; end
   if nargout > 2, varargout{3} = []     ; end
end

%--------------------------------------------------------------------------

function varargout = relperm(s, opt, varargin)
   [s1, s2, den] = modified_saturations(s, opt);

   if isempty(opt.reg)
      n   = opt.n;
      kwm = opt.kwm;
   else
      assert(size(s,1)==numel(opt.reg));
      I = [2*opt.reg-1 2*opt.reg];
      n = opt.n(I);
      kwm = opt.kwm(I);
   end
   varargout{1}    = [ kwm(:,1) .* s1 .^ n(:,1), kwm(:,2) .* s2 .^ n(:,2)];

   if nargout > 1,
      null = zeros(size(s1));
      varargout{2} = [ kwm(:,1) .* n(:,1) .* s1 .^ (n(:,1) - 1), ...
                       null, null                       , ...
                       kwm(:,2) .* n(:,2) .* s2 .^ (n(:,2) - 1)] ./ den;
   end

   if nargout > 2,
      a = n .* (n - 1);
      varargout{3} = [ kwm(:,1) .* a(:,1) .* s1 .^ (n(:,1) - 2), ...
                       kwm(:,2) .* a(:,2) .* s2 .^ (n(:,2) - 2)] ./ den;
   end
end

%--------------------------------------------------------------------------

function [s1, s2, den] = modified_saturations(s, opt)
   if isempty(opt.reg)
      sr  = opt.sr;
   else
      assert(size(s,1)==numel(opt.reg));
      sr = opt.sr([2*opt.reg-1 2*opt.reg]);
   end

   den = 1 - sum(sr,2);
   s1  = (    s(:,1) - sr(:,1)) ./ den;  s1(s1 < 0) = 0;  s1(s1 > 1) = 1;
   s2  = (1 - s(:,1) - sr(:,2)) ./ den;  s2(s2 < 0) = 0;  s2(s2 > 1) = 1;
end
