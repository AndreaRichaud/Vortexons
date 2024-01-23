function [Laplacian] = FFT_Laplacian(psi,Mat_k_squared)

psi_hat = fftshift(fft2(psi));

% Fourier Space Laplacian
Laplacian = ifft2(ifftshift(- Mat_k_squared .* psi_hat));

end

