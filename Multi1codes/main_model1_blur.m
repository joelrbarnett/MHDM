clear all

%x=255*phantom(256);
%x=double(imread('cameraman.tif'));
x=double(imread('parrot.png'));
n=size(x,1);

sigma = 1;
K = 10;
noise = reshape(mean(exponential_rnd(1/sigma,n*n,K),2),n,n);

%h = fspecial('motion',5,30); %
h = fspecial('gaussian',7,2);
hx = blurA(x,h);

y = hx .* noise;
f=max(y,1);
load parrot_gaussK10;

disp('====================================')

%% Our method

lambda=0.07; sigma=3; tau=3; alpha=16;
options.bound = 'per';
[lambda,sigma,tau,alpha]
maxit=500;

tic;
%un=DZconvex1_blur_cpu(f,sigma,lambda,alpha,options,tau,293,h,y);
[un,fun_value,diff,it]=DZconvex1_blur(f,sigma,lambda,alpha,options,tau,maxit,h,x,y);
t=toc

figure,imshow(uint8(un));
psnr(x,min(un,255))
figure,plot(fun_value(2:it)),xlabel('iterations'),ylabel('function values')  
 
