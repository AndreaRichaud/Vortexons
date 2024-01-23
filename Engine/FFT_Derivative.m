function [Gradient_x, Gradient_y, Laplacian] = FFT_Derivative(psi,Mat_k_x,Mat_k_y)


psi_hat = fftshift(fft2(psi));

% Fourier Space Gradient
Gradient_x = ifft2(ifftshift(  1i*Mat_k_x  .* psi_hat));
Gradient_y = ifft2(ifftshift(  1i*Mat_k_y  .* psi_hat));


% Fourier Space Laplacian
Laplacian = ifft2(ifftshift(-( Mat_k_x.^2 + Mat_k_y.^2 ) .* psi_hat));

end

