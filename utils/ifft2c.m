function res = ifft2c(x)
%
%
% res = ifft2c(x)
% 
% orthonormal centered 2D ifft
%
% (c) Michael Lustig 2005

% res = sqrt(length(x(:)))*ifftshift(ifft2(fftshift(x)));

% modified by Florian Knoll for compatibility with odd matrices, according to:
% https://www.mathworks.com/matlabcentral/newsreader/view_thread/285244
res = sqrt(length(x(:)))*fftshift(ifft2(ifftshift(x)));


