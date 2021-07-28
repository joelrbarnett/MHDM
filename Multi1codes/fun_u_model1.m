function y=fun_u_model1(uk,lambda,f,alpha,pn,tau,oldun,options)
y=(1./uk-f./(uk).^2)+alpha*(1./f-1./sqrt(f.*uk))-lambda*div(pn(:,:,1),pn(:,:,2),options)+(1/tau)*(uk-oldun);