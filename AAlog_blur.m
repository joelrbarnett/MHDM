%Multiplicative denoising and deblurring  (Semi-implicit)
%Author: Joel Barnett

function [u]=AAlog_blur(f,xk,dt,lambda,T,epsilon, maxIter)

f = double(f);
[Image_h,Image_w]=size(f);

%initialization: u^0=f/Txk Here, I'm epsilon regularizing to prevent /0
init=sum(sum(f./(epsilon+imfilter(xk,T,'symmetric','same'))))/(Image_h*Image_w); 
u=init.*ones(size(f));
for n=1:maxIter
        %These aren't square matrices
        DXF=dxf(u); %forward diff in x, for all y indexes
        DXB=dxb(u); % back diff in x, for all y indexes
        DYF=dyf(u); %forward diff in y, for all x indexes
        DYB=dyb(u); %back diff in y, for all x indexes
% EPSILON REGUALRIZATION,
        c1= 1./sqrt(epsilon^2 + DXF(:,2:end-1).^2 + DYF(2:end-1,:).^2);
        c2= 1./sqrt(epsilon^2 + DXB(:,2:end-1).^2 + DYF(1:end-2,:).^2);
        c3= c1; % 1./sqrt(epsilon^2 + DXF(:,2:end-1).^2 + DYF(2:end-1,:).^2);
        c4= 1./sqrt(epsilon^2 + DXF(:,1:end-2).^2 + DYB(2:end-1,:).^2);
% for div(grad(u)/|u||grad u|), we need to divide the ci terms by |u(i,j)|
% appropriately regularized
        c1=c1./abs(sqrt(epsilon^2+u(2:end-1,2:end-1).^2)); %divide by |u(i,j)|
        c2=c2./abs(sqrt(epsilon^2+u(1:end-2,2:end-1).^2)); %divide by |u(i-1,j)|
        c3= c1; %c3./abs(sqrt(epsilon^2+u(2:end-1,2:end-1).^2)); %divide by |u(i,j)|
        c4=c4./abs(sqrt(epsilon^2+u(2:end-1,1:end-2).^2)); %divide by |u(i,j-1)|
% form |grad u|/u|u| term. Here, I tried centered differences for |grad(u)|
    DXC = dxc(u); DYC=dyc(u);
    UabsU= (u(2:end-1,2:end-1)+epsilon).*sqrt(epsilon^2+u(2:end-1,2:end-1).^2);
    gradU_UabsU = sqrt(DXC(:,2:end-1).^2 + DYC(2:end-1,:).^2)./UabsU;

% fidelity term
    Tuxk=imfilter(u.*xk,T,'symmetric','same');
    fidelity = (1./(epsilon+Tuxk(2:end-1,2:end-1)) - f(2:end-1,2:end-1)./(epsilon+Tuxk(2:end-1,2:end-1).^2));
    fidelity = imfilter(fidelity,T,'symmetric','same');
    fidelity = lambda.*fidelity.*xk(2:end-1,2:end-1);
%update interior grid points            
        u(2:end-1,2:end-1)=1./(1+dt.*(c1+c2+c3+c4)).*...
            (u(2:end-1,2:end-1) - dt.*fidelity + dt.* gradU_UabsU+...
            dt.*(c1.*u(3:end,2:end-1)+c2.*u(1:end-2,2:end-1)+c3.*u(2:end-1,3:end)+c4.*u(2:end-1,1:end-2)));
            
%Update boundary conditions
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

%Cetnered difference x and y derivatives at interior grid points
function DXC=dxc(M)
    DXC = M(3:end,:) - M(1:end-2,:);
end
function DYC= dyc(M)
    DYC = M(:,3:end) - M(:,1:end-2);
end





