function y=fun_u_model1_blur(uk,lambda,f,alpha,pn,tau,oldun,options,h,Auk)
temp=sqrt(Auk);
y=blurA((Auk-f)./((Auk).^2)+alpha*(temp-sqrt(f))./(f.*temp),h,'s')-lambda*div(pn(:,:,1),pn(:,:,2),options)+(1/tau)*(uk-oldun);