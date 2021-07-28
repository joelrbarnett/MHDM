function [un,fun_value,diff,it]=DZconvex1(f,sigma,lambda,alpha,options,tau,maxit,x,y)

pn=0*grad(f)+0.01;
bun=f;
un=bun;

grad_u=grad(un,options);
tv_u=sqrt(grad_u(:,:,1).^2+grad_u(:,:,2).^2);
fun_value=zeros(1,maxit);
fun_value(1)=mean2(log(un)+f./un)+alpha*mean2((sqrt(un./f)-1).^2)+lambda*mean2(tv_u);

diff=1; it=1;
while diff(it)>=4e-4 && it<=maxit

    if mod(it,50)==0
        [it, fun_value(it), psnr(x,min((mean2(y)/mean2(un))*un,255))]
    end

    %updating pn
    tmp=sigma*lambda*grad(bun,options)+pn;
    t1=tmp(:,:,1); t2=tmp(:,:,2); ntmp=max(sqrt(t1.^2+t2.^2),1);
    tmp(:,:,1)=t1./ntmp;tmp(:,:,2)=t2./ntmp;
    pn=tmp;

    %updating un
    oldun=un;
    %Newton method
    uk=un;
    for k=1:10
        gw=(-1./(uk.^2)+2*(f./(uk.^3)))+alpha./(2*sqrt(f.*(uk.^3)))+1/tau;
        gpw=fun_u_model1(uk,lambda,f,alpha,pn,tau,oldun,options);
        uk=uk-gw.\gpw;
        resk(k)=norm(fun_u_model1(uk,lambda,f,alpha,pn,tau,oldun,options),'fro');
        if resk(k)/resk(1)<1e-6
            break;
        end
    end
    un=max(uk,1);
    bun=2*un-oldun;

    grad_u=grad(un,options);
    tv_u=sqrt(grad_u(:,:,1).^2+grad_u(:,:,2).^2);
    fun_value(it+1)=mean2(log(un)+f./un)+alpha*mean2((sqrt(un./f)-1).^2)+lambda*mean2(tv_u);
    if it<=20
        diff(it+1)=1;
    else
        diff(it+1)=(fun_value(it)-fun_value(it+1))/fun_value(it);
    end
    
    
    it=it+1;

end

un=(mean2(y)/mean2(un))*un;
