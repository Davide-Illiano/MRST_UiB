function [state,dt,report, sreport] = impesTPFA_URC_time(state0, G, T, fluid, dt, pv, varargin)
% [state,dt,report, sreport] = impesTPFA_URC_time(state0, G, T, fluid, dt,pv, varargin)
% URC formulation with time step search. Used impesTPFA_URC.
      opt = struct('bc', [], 'src', [], 'wells', [], ...
                'ATol', 5.0e-7,...
                'ntol', 1e-5,...
                'URC_p',true,...
                'URC_s','linear',...
                'URC_opt',struct('maxnewt',10,'lstrials',3),...% only used if URC_s='nonlinear'
                'report', [],...
                'cfl_factor',10,...
                'max_time_cut',10,...
                'verbose',mrstVerbose,...
                'rock',[]);
     opt = merge_options(opt, varargin{:});
     ok=false;
     report0=opt.report;
     it=0;
     init_state=state0;
     while ~ok && it < opt.max_time_cut
        [state,dt,report,sreport] = impesTPFA_URC(state0, G, T, fluid, dt, pv,  ...
           'wells', opt.wells,'src',opt.src,'bc',opt.bc,...
           'ATol', opt.ATol, ...
           'RTol', 5.0e-100, ...
           'ntol',opt.ntol,...
           'EstimateTimeStep', true,...
           'DynamicMobility',false, ...
           'face_update',false,...
           'cfl_factor',opt.cfl_factor,...
           'URC_p', opt.URC_p,...
           'URC_s', opt.URC_s,...
           'LineSearch',false,...
           'report', report0,'rock',opt.rock,...
           'init_state',init_state);
           if(sreport.success==false)
              if(sreport.upwind_changed)
                 dispif(opt.verbose,['Redo timestep due to upwind change in impes URC to', num2str(dt/day), ' days'])
                 init_state=sreport.init_state;
              else
                 dt=dt/2.0;
              end
              ok=false;
              it=it+1;
              dispif(opt.verbose,['Cutting timestep in impes URC to', num2str(dt/day), ' days'])
              if(opt.verbose)
                 disp(sreport)
              end
           else
              dispif(opt.verbose,['Successful step with time step in impes URC to', num2str(dt/day), ' days'])
              ok=true;
           end
     end
     if(~ok)
        error('Did find a valid time step')
     end

end

