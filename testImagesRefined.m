% is the original image. 
% model: minimize int[lambda*( log(Tuxk) + f/Tuxk] + ak lambda TV(log(uxk))+ TV(log(u)) over u. 
% Multiscale: u = u0*u1*...*uk*..., decomposition of u. w We refer to
% partial prodcut xk=u0*u1*...*uk

clear all

%for saving
filePrefix="cameraman_noise_refined/";
figPrefix="cameraman_";

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%choose images
images=["barbara.png","cameraman.tif","pollen.tif","lena_gray_512.tif","peppers_gray.tif","mandril_gray.tif","circles.tif","geometry.tif"];
folder_path='../Test Images/';

%F_orig=imread([folder_path,'cameraman.tif']); %cameraman
%F_orig=imread([folder_path,'barbara.png']);
F_orig=imread(char(folder_path+images(2))); 

F_orig=double(F_orig(:,:));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%setup parameters
[n,m]=size(F_orig);

maxIters=50; %time iterations in solving for wk
numScales=16; %number of decompositions
lambda0=0.01; %0.02;
ak=1; %for tighter version...aditional scaling factor
epsilon= 0.01; %for regularizing TV
dt=0.01; %0.025; %timestep
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Form noisy image: switch T to identity for only noise. Switch GamNoise to
%ones(size(F_orig)) for only blurring. 
%%% Blur and Gamma noise %%%
a=25; %gamma noise with mean 1, standard deviation 0.2. 
GamNoise=gamrnd(a,1/a,size(F_orig));
%GamNoise=ones(size(F_orig)); %For only deblurring, no noise

%T=fspecial('gaussian',[3 3],sqrt(2)); %blurring component/operator
T=fspecial('average',[1 1]); %identity, for no blur
F_blur=imfilter(F_orig,T,'symmetric','same'); %create blurred image
F_data=F_blur.*GamNoise; %multiply noise into blurred image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set up arrays to hold metrics of fit performance.
% Bregman distance between F_data and F_blur
delta= sum(sum(F_data./F_blur + log(F_blur) - log(F_data) -ones(size(F_orig))));
%delta=norm(F_blur-F,'fro')/norm(F_blur,'fro'); %rel reisdual

rmse_final=zeros(numScales,1);
xk_f_norm2=zeros(numScales,1);
D_f_data_Txk_D_f_Data_f_orig=zeros(numScales,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Image Storage Arrays:
ukArray=zeros([[m n 1], numScales]);
xkArray=zeros([[m n 1], numScales]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Run decomposition
xk=ones(size(F_data));
qk=2;
q=3;
lambda=lambda0;
for k=1:numScales    
    %get decomposed piece uk.
    uk=AAlog_blur_refined(F_data,xk,dt,lambda,ak,T, epsilon, maxIters);
    %uk=AAlog_blur_tight(F_data,xk,dt,lambda,ak,T, epsilon, maxIters);
    %uk=AAlog_blur(F_data,xk,dt,lambda,T,epsilon,maxIters);
    %update xk and lambda_k
    xk=uk.*xk;
    lambda=lambda* q; %alternatively, use *qk for adaptive lambda
    ak=1/(k^(3/2)); %update ak

    %Store images:
    ukArray(:,:,1,k)=uk; %image piece, single scale
    xkArray(:,:,1,k)=xk; %updated multiscale image
        
    %capture errors
    xk_f_norm2(k)=norm(F_orig-xk,'fro'); %difference from true image
    rmse_final(k)=xk_f_norm2(k)/sqrt(m*n);   %RMSE error from true image
    %Bregman distance from F_data to Txk (restored, blurred image)
    D_f_data_Txk=sum(sum(F_data./imfilter(xk,T,'symmetric','same') + log(imfilter(xk,T,'symmetric','same')) - log(F_data)-ones(size(F_orig))  ));
    D_f_data_Txk_D_f_Data_f_orig(k) = D_f_data_Txk^2/(delta^2);%ratio of bregman distances
    
    % For adaptive lambda updates, specifies how to choose qk
    val=D_f_data_Txk_D_f_Data_f_orig(k); %take last bregman ratio
    val0=D_f_data_Txk_D_f_Data_f_orig(1);%initial bregman ratio
    qk = 1.5/(1+5*exp(-(val-1)))^10+1;  %sigmoid qk update
    %qk=1+log(val); %log qk update
    %qk=2.5/(val0-1)*(val-1)+1; %linear fit between (1,1), (val0,2)
    qk=min(qk,2.0); %choose bounds of qk
    qk=max(qk,1.05);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find min error and stopping criteria
%Min error and index of that error. Divide by (m*n) to get RMSE
[minVal,mink]=min(xk_f_norm2)
%k_star = max_k D(F_data,Txk)^2/ D(F_data,Tu)^2 \geq tau, with tau>1.
k_star=min(find((D_f_data_Txk_D_f_Data_f_orig<=1)==1))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot results
noisyRMSE=norm(F_orig-F_data,'fro')/sqrt(m*n); %original RMSE error

%Plot the original, degraded, and recovered images
figure()
image(F_orig); axis image; axis off; colormap(gray(256));
title('Original Image','FontSize',16)
% 
figure()
image(F_data); axis image; axis off; colormap(gray(256));
title(['Noisy, blurred image. RMSE=',num2str(noisyRMSE)],'FontSize',16)
figName=char(filePrefix+figPrefix+"noisy.png");
saveas(gcf,figName)
%
figure()
image(xkArray(:,:,1,mink)); axis image; axis off; colormap(gray(256));
subtitle="\lambda_k ="+lambda0+ "*"+q +"^k, \epsilon= "+epsilon+", \Delta t= "+ dt+", maxIters= "+maxIters;
titleStr="Restored image: image(xk), k="+num2str(mink)+" RMSE="+num2str(minVal/sqrt(m*n));
title([titleStr,subtitle],'FontSize',16)
figName=char(filePrefix+figPrefix+"restored.png");
saveas(gcf,figName)

%Plot montage of multiscales near optimal
figure()
numPlots =11;
montage(xkArray(:,:,:,end-numPlots:end),gray(256))

titleStr = "Multiscales k = "+num2str(numScales -numPlots)+" through k= "+num2str(numScales);
title(titleStr,'FontSize',16)
figName=char(filePrefix+figPrefix+"multiscales.png");
saveas(gcf,figName)

%Plot the RMSE and Total error, as well as the stopping criterion plot
figure()
subplot(1,2,1)
yyaxis left
plot(1:numScales,rmse_final)
xlabel('Multiscales','FontSize',16)
ylabel('RMSE','FontSize',16)
title('RMSE vs multiscale-decompositions','FontSize',16)
yyaxis right
plot(1:numScales, xk_f_norm2)
ylabel('|xk-f|','FontSize',16)

subplot(1,2,2)
semilogy(1:numScales, D_f_data_Txk_D_f_Data_f_orig,1:numScales,ones(numScales,1))
xlabel('k','FontSize',16)
title('D(F_{data},Tx_k)^2/D(F_{data},Tu)^2','FontSize',16)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
