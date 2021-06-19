%Multiplicative denoising and deblurring  (Semi-implicit)
%Author: Joel Barnett

function [u]=AAlog_blur_refined(f,xk,dt,lambda,ak,T,epsilon, maxIter)
f = double(f);
[Image_h,Image_w]=size(f);

%initialization: u^0=f/Txk Here, I'm epsilon regularizing to prevent /0
init=sum(sum(f./(epsilon+imfilter(xk,T,'symmetric','same'))))/(Image_h*Image_w); 
u=init.*ones(size(f));
phi=zeros(size(f));
for n=1:maxIter
    w=u.*xk;
    %These aren't square matrices
    DXF=dxf(u); %forward diff in x, for all y indexes
    DXB=dxb(u); % back diff in x, for all y indexes
    DYF=dyf(u); %forward diff in y, for all x indexes
    DYB=dyb(u); %back diff in y, for all x indexes
    
    DXFphi=dxf(phi); %forward diff in x, for all y indexes
    DXBphi=dxb(phi); % back diff in x, for all y indexes
    DYFphi=dyf(phi); %forward diff in y, for all x indexes
    DYBphi=dyb(phi); %back diff in y, for all x indexes
    
    DXFw=dxf(w); %forward diff in x, for all y indexes
    DXBw=dxb(w); % back diff in x, for all y indexes
    DYFw=dyf(w); %forward diff in y, for all x indexes
    DYBw=dyb(w); %
% EPSILON REGUALRIZATION of div terms
    %u terms
    c1= 1./sqrt(epsilon^2 + DXF(:,2:end-1).^2 + DYF(2:end-1,:).^2);
    c2= 1./sqrt(epsilon^2 + DXB(:,2:end-1).^2 + DYF(1:end-2,:).^2);
    c3= c1; % 1./sqrt(epsilon^2 + DXF(:,2:end-1).^2 + DYF(2:end-1,:).^2);
    c4= 1./sqrt(epsilon^2 + DXF(:,1:end-2).^2 + DYB(2:end-1,:).^2);
    
    %phi terms
    cPhi1= 1./sqrt(epsilon^2 + DXFphi(:,2:end-1).^2 + DYFphi(2:end-1,:).^2);
    cPhi2= 1./sqrt(epsilon^2 + DXBphi(:,2:end-1).^2 + DYFphi(1:end-2,:).^2);
    cPhi3= cPhi1; % 1./sqrt(epsilon^2 + DXF(:,2:end-1).^2 + DYF(2:end-1,:).^2);
    cPhi4= 1./sqrt(epsilon^2 + DXFphi(:,1:end-2).^2 + DYBphi(2:end-1,:).^2);
    
    %w terms
    cw1= 1./sqrt(epsilon^2 + DXFw(:,2:end-1).^2 + DYFw(2:end-1,:).^2);
    cw2= 1./sqrt(epsilon^2 + DXBw(:,2:end-1).^2 + DYFw(1:end-2,:).^2);
    cw3= c1; % 1./sqrt(epsilon^2 + DXF(:,2:end-1).^2 + DYF(2:end-1,:).^2);
    cw4= 1./sqrt(epsilon^2 + DXFw(:,1:end-2).^2 + DYBw(2:end-1,:).^2);
    
% Compute phi update
    if n==1 %after first iteration, gradPhi is already updated 
        %DXCphi=dxc(phi); DYCphi=dyc(phi);
        gradPhi = sqrt(epsilon^2 + DXFphi(:,2:end-1).^2+DYFphi(2:end-1,:).^2);
    end
    %S_n=integral( log(u)*phi)/integral(|gradPhi|)
    S_n=sum(sum(log(u(2:end-1,2:end-1)+epsilon).*phi(2:end-1,2:end-1)./gradPhi));
    
    phi(2:end-1,2:end-1) = 1./(1+S_n*dt*(cPhi1+cPhi2+cPhi3+cPhi4)).*(phi(2:end-1,2:end-1)+...
        dt.*(log(u(2:end-1,2:end-1)+epsilon) + S_n*(cPhi1.*phi(3:end,2:end-1)+...
        cPhi2.*phi(1:end-2,2:end-1)+cPhi3.*phi(2:end-1,3:end)+cPhi4.*phi(2:end-1,1:end-2))) );
    %now, phi=phi^{n+1}, so we can recompute the gradient for u update
    DXCphi=dxc(phi); DYCphi=dyc(phi);
    gradPhi = sqrt(epsilon^2 + DXCphi(:,2:end-1).^2+DYCphi(2:end-1,:).^2);
    
% for div(grad(u)/|u||grad u|), we need to divide the ci terms by |u(i,j)|
% appropriately regularized
    %u terms
    c1=c1./abs(sqrt(epsilon^2+u(2:end-1,2:end-1).^2)); %divide by |u(i,j)|
    c2=c2./abs(sqrt(epsilon^2+u(1:end-2,2:end-1).^2)); %divide by |u(i-1,j)|
    c3= c1; %c3./abs(sqrt(epsilon^2+u(2:end-1,2:end-1).^2)); %divide by |u(i,j)|
    c4=c4./abs(sqrt(epsilon^2+u(2:end-1,1:end-2).^2)); %divide by |u(i,j-1)|
    %w terms
    cw1=cw1./abs(sqrt(epsilon^2+w(2:end-1,2:end-1).^2)); %divide by |w(i,j)|
    cw2=cw2./abs(sqrt(epsilon^2+w(1:end-2,2:end-1).^2)); %divide by |w(i-1,j)|
    cw3= cw1; %c3./abs(sqrt(epsilon^2+u(2:end-1,2:end-1).^2)); %divide by |w(i,j)|
    cw4=cw4./abs(sqrt(epsilon^2+w(2:end-1,1:end-2).^2)); %divide by |w(i,j-1)|
% form |grad u|/u|u| term. Here, I tried centered differences for |grad(u)|
    DXC = dxc(u); DYC=dyc(u); DXCw=dxc(w); DYCw=dyc(w);
    %u terms
    UabsU= (u(2:end-1,2:end-1)+epsilon).*sqrt(epsilon^2+u(2:end-1,2:end-1).^2);
    gradU_UabsU = sqrt(DXC(:,2:end-1).^2 + DYC(2:end-1,:).^2)./UabsU;
    %w terms
    WabsW= (w(2:end-1,2:end-1)+epsilon).*sqrt(epsilon^2+w(2:end-1,2:end-1).^2);
    gradW_WabsW = sqrt(DXCw(:,2:end-1).^2 + DYCw(2:end-1,:).^2)./WabsW;

% fidelity term
    Tuxk=imfilter(u.*xk,T,'symmetric','same');
    fidelity = (1./(epsilon+Tuxk(2:end-1,2:end-1)) - f(2:end-1,2:end-1)./(epsilon+Tuxk(2:end-1,2:end-1).^2));
    fidelity = imfilter(fidelity,T,'symmetric','same');
    fidelity = lambda.*fidelity.*xk(2:end-1,2:end-1);
%update interior grid points of u            
    u(2:end-1,2:end-1)=1./(1+dt.*(c1+c2+c3+c4)+dt*lambda*ak.*(cw1+cw2+cw3+cw4).*xk(2:end-1,2:end-1).^2).*...
        (u(2:end-1,2:end-1) - dt.*fidelity + dt.* (gradU_UabsU+ak*lambda.*gradW_WabsW.*xk(2:end-1,2:end-1))+...
        dt.*(c1.*u(3:end,2:end-1)+c2.*u(1:end-2,2:end-1)+c3.*u(2:end-1,3:end)+c4.*u(2:end-1,1:end-2))+...
        dt*ak*lambda.*xk(2:end-1,2:end-1).*(cw1.*w(3:end,2:end-1)+cw2.*w(1:end-2,2:end-1)+cw3.*w(2:end-1,3:end)+cw4.*w(2:end-1,1:end-2)) -phi(2:end-1,2:end-1)./gradPhi );
            
%Update boundary conditions for u and phi
    %Left and Right Boundary Conditions
    u(1,2:Image_w-1)=u(2,2:Image_w-1); 
    u(Image_h,2:Image_w-1)=u(Image_h-1,2:Image_w-1); 
    
    %Top and Bottom boundary conditions
    u(2:Image_h-1,1)=u(2:Image_h-1,2); 
    u(2:Image_h-1,Image_w)=u(2:Image_h-1,Image_w-1);
    
    %Corner Boundary Conditions:
    u(1,1)=u(2,2); u(1,Image_w) = u(2,Image_w-1); 
    u(Image_h,1) =u(Image_h-1,2); 
    u(Image_h,Image_w)=u(Image_h-1,Image_w-1);
    
    %Left and Right Boundary Conditions
    phi(1,2:Image_w-1)=phi(2,2:Image_w-1); 
    phi(Image_h,2:Image_w-1)=phi(Image_h-1,2:Image_w-1); 
    
    %Top and Bottom boundary conditions
    phi(2:Image_h-1,1)=phi(2:Image_h-1,2); 
    u(2:Image_h-1,Image_w)=u(2:Image_h-1,Image_w-1);
    
    %Corner Boundary Conditions:
    phi(1,1)=phi(2,2); phi(1,Image_w) = phi(2,Image_w-1); 
    phi(Image_h,1) =phi(Image_h-1,2); 
    phi(Image_h,Image_w)=phi(Image_h-1,Image_w-1);
end
end

%function dxf performs the forward difference in x, the first index-position
%of matrix M, and returns DXF, a matrix of size size(M,1)-2 x size(M,2)
%representing the forward-difference approximations to the x-derviative of
%the matrix M at interior grid points (2:end-1,1:end)
function DXF=dxf(M)
    DXF =M(3:end,:) -M(2:end-1,:);
end

%backwards difference x-derivative approximation at interior grid points
% (2:end-1,1:end)
function DXB=dxb(M) 
    DXB=M(2:end-1,:) - M(1:end-2,:);
end

%forwards difference y-derivative approximation at interior grid points
% (1:end, 2:end-1)
function DYF=dyf(M) 
    DYF=M(:,3:end) - M(:,2:end-1);
end

%backward difference y-derivative approximation at interior grid points
% (1:end, 2:end-1)
function DYB=dyb(M)
    DYB=M(:,2:end-1) - M(:,1:end-2);
end

%Centered difference x and y derivatives at interior grid points
function DXC=dxc(M)
    DXC = M(3:end,:) - M(1:end-2,:);
end
function DYC= dyc(M)
    DYC = M(:,3:end) - M(:,1:end-2);
end






