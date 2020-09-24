function lfluid= blackOilLinearTemp(c_p,h_t,fluid,p_ref,T_ref,varargin)
    opt=struct('only_enthalpy',true);
    [p,T]=initVariablesADI(p_ref,T_ref);
    if(isfield(fluid,'bG'))
        names = {'W', 'O', 'G'};
    else
        names = {'W', 'O'};
    end
   for i = 1:numel(names)
       n = names{i};       
       rho=fluid.(['rho', n, 'S']).*fluid.(['b', n])(p);
       cp_rho=rho.jac{1}(1,1);ct_rho=rho.jac{2}(1,1);
       rho=double(rho);
       assert(cp_rho>=0);
       mu=fluid.(['mu', n])(p);
       cp_mu=mu.jac{1}(1,1);ct_mu=mu.jac{2}(1,1);
       mu=double(mu);
       ct_h=c_p(i);
       h=h_t(i);       
       obj{i}=linearFluid(cp_rho,ct_rho,ct_h,cp_mu,ct_mu,rho,h,mu,T_ref,p_ref)
       if(opt.only_enthalpy)
        fluid.(['h', n]) =@(p,T) obj{i}.enthalpy(p,T);
        fluid.(['u', n]) =@(p,T) obj{i}.enthalpy(p,T)-p./obj{i}.density(p,T);
       end
   end
   if(opt.only_enthalpy)
      lfluid=fluid;
   else
      lfluid=fluid2BlackOil(fluid_o,fluid_pvt,opt.p_ref,opt.T_ref);
   end
end
