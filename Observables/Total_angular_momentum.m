function [L_z_tot] = Total_angular_momentum(N_a,N_b,psi_a,psi_b,Mat_x,Mat_y,Mat_k_x,Mat_k_y,hbar,dx,dy)
[Gradient_x_a,Gradient_y_a,~] = FFT_Derivative(psi_a,Mat_k_x,Mat_k_y);
[Gradient_x_b,Gradient_y_b,~] = FFT_Derivative(psi_b,Mat_k_x,Mat_k_y);

L_z_a = sum(sum(N_a*dx*dy*real(-1i*hbar * conj(psi_a).*(Mat_x.*Gradient_y_a - Mat_y.*Gradient_x_a))));
L_z_b = sum(sum(N_b*dx*dy*real(-1i*hbar * conj(psi_b).*(Mat_x.*Gradient_y_b - Mat_y.*Gradient_x_b))));

L_z_tot = L_z_a+L_z_b;
end

