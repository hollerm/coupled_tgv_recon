function res = fft2c(x)
% res = fft2c(x)
% orthonormal forward 2D FFT

res = 1/sqrt(length(x(:)))*fftshift(fft2(ifftshift(x)));

