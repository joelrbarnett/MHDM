function p = psnr(x,y)

% psnr - compute the Peack Signal to Noise Ratio, defined by :
%       PSNR(x,y) = 10*log10( max(max(x),max(y))^2 / |x-y|^2 ).
%
%   p = psnr(x,y);
%

d = mean( mean( (x(:)-y(:)).^2 ) );
m1 = max( abs(x(:)) );
m2 = max( abs(y(:)) );
m = m1;

p = 10*log10( m^2/d );

%
% Part of CurvL1MultDenoise Version:100
% Written by: Jalal Fadili, GREYC CNRS-ENSICAEN-Univ. Caen
%      	      Mila Nikolova, CMLA CNRS-ENS Cachan
%	      Sylvain Durand, MAP5 CNRS-Univ. Paris 5
% Created March 2009
% This is Copyrighted Material
% E-mail: Jalal.Fadili@greyc.ensicaen.fr
%
% Modified Wed Dec 16 18:41:37 CET 2009
%

