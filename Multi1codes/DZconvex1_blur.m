function [un,fun_value,diff,it]=DZconvex1_blur(f,sigma,lambda,alpha,options,tau,maxit,h,x,y)

pn=0*grad(f)+0.1;
bun=f;
un=bun;

Aun=blurA(un,h);
grad_u=grad(un,options);
tv_u=sqrt(grad_u(:,:,1).^2+grad_u(:,:,2).^2);
fun_value=zeros(1,maxit);
fun_value(1)=mean2(log(Aun)+f./Aun)+alpha*mean2((sqrt(Aun./f)-1).^2)+lambda*mean2(tv_u);

diff=1; it=1;
while diff(it)>=1e-4 && it<=maxit

   if mod(it,10)==0
        [it,fun_value(it),psnr(x,min((mean2(y)/mean2(un))*un,255))]
        figure(100),imshow(uint8((mean2(y)/mean2(un))*un))
   end
   
   %update pn
   tmp=lambda*sigma*grad(bun,options)+pn;
   t1=tmp(:,:,1); t2=tmp(:,:,2); ntmp=max(sqrt(t1.^2+t2.^2),1);
   tmp(:,:,1)=t1./ntmp;tmp(:,:,2)=t2./ntmp;
   pn=tmp;
      
   %update un
   oldun=un;
   %Newton method
   uk=un;
   for k=1:12
       Auk=blurA(uk,h);
       gw=blurA(-1./(Auk.^2)+2*(f./(Auk.^3)),h,'s')+(alpha/2).*blurA(1./sqrt(f.*(Auk.^3)),h,'s')+1/tau;
       gpw=fun_u_model1_blur(uk,lambda,f,alpha,pn,tau,oldun,options,h,Auk);
       uk=uk-gw.\gpw;
       resk(k)=norm(fun_u_model1_blur(uk,lambda,f,alpha,pn,tau,oldun,options,h,Auk),'fro');
       if resk(k)/resk(1)<1e-3
           break;
       end
   end
   un=max(uk,1);
   bun=2*un-oldun;
   
   %find energy
   Aun=blurA(un,h);
   grad_u=grad(un,options);
   tv_u=sqrt(grad_u(:,:,1).^2+grad_u(:,:,2).^2);
   fun_value(it+1)=mean2(log(Aun)+f./Aun)+alpha*mean2((sqrt(Aun./f)-1).^2)+lambda*mean2(tv_u);
    if it<=20
        diff(it+1)=1;
    else
        diff(it+1)=(fun_value(it)-fun_value(it+1))/fun_value(it);
    end
    
    
    it=it+1;
end
   
un=(mean2(y)/mean2(un))*un;