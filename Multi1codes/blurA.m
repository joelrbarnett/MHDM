function hu = blurA(u,h,star);
% This is to calculate Au or A^* u;
% hu = blurA(u,h)-----------Au;
% hsu = blurA(u,h,'s')-------A^*u;

center=round((size(h)+1)/2);
[Ny,Nx]=size(u);
trans=@(X) 1/sqrt(Ny*Nx)*fft2(X);
itrans=@(X) sqrt(Ny*Nx)*ifft2(X);

blur_A=zeros(size(u)); blur_A(1:size(h,1),1:size(h,2))=h;
blur_matrix_trans=fft2(circshift(blur_A,1-center));

if nargin == 2
  hu=real(itrans(trans(u).*blur_matrix_trans));
else   
  hsu=real(itrans(trans(u).*conj(blur_matrix_trans)));
  hu=hsu;
end

end