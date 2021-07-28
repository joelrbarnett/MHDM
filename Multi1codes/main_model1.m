 clear all
% 
% %x=255*phantom(256);
% %x=max(min(x,255),1);
% %x=double(imread('cameraman.tif'));
% x=double(imread('parrot.png'));
% n=size(x,1);
% 
% sigma = 1;
% K = 10;
% noise = reshape(mean(exponential_rnd(1/sigma,n*n,K),2),n,n);
% y = x .* noise;
% f=max(y,1);
% figure,imshow(f,[])
% %load parrot_noiseK6;

disp('====================================')

%% Our method
folder_path="../Test_Images_plus1/"; %read images with no zero values
fileNames=["barbara","cameraman","pollen","mandril","circles","geometry","disc_square"]; 
images=["barbara.png","cameraman.tif","pollen.tif","mandril_gray.tif","circles.tif","geometry.tif","disc_square.png"];
imagesPNG=["barbara.png","cameraman.png","pollen.png","mandril.png","circles.png","geometry.png","disc_square.png"];

f_orig=double(imread(char(folder_path+imagesPNG(4))));
%%%Gamma noise %%%
rng(10); %fix seed across all runs
a=25; %gamma noise with mean 1, standard deviation 0.2. 
GamNoise=gamrnd(a,1/a,size(f_orig));
f_data=f_orig.*GamNoise;

lambda=0.06; sigma=3; tau=3; alpha=16; %20*sqrt(6)/9;
[lambda,sigma,tau,alpha]
options.bound = 'per';
maxit=1000;

%Find the best lambda to use
% lambdaArray=0.05:0.01:0.2;
% snrArray=zeros(size(lambdaArray));
% iter=1;
% for lambda = lambdaArray
%     tic;
%     %un=DZconvex1_cpu(f,sigma,lambda,alpha,options,tau,179,x,y);
%     [un,fun_value,diff,it]=DZconvex1(f_data,sigma,lambda,alpha,options,tau,maxit,f_orig,f_data);
%     t=toc
%     snrArray(iter)=20*log10(norm(un,'fro')/norm(un-f_orig,'fro'));
%     iter=iter+1;
% end

%Get the best snr value, and use that lambda for final run
% [~,maxIt]=max(snrArray);
% lambda=lambdaArray(maxIt);



tic;
%un=DZconvex1_cpu(f,sigma,lambda,alpha,options,tau,179,x,y);
[un,fun_value,diff,it]=DZconvex1(f_data,sigma,lambda,alpha,options,tau,maxit,f_orig,f_data);
t=toc
figure,imshow(uint8(un));
%psnr(f_orig,min(un,255))
%real snr 
20*log10(norm(un,'fro')/norm(un-f_orig,'fro'))
figure,plot((2:it)),xlabel('iterations'),ylabel('function values')

