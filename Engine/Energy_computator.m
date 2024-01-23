function [Ene] = Energy_computator(psi_a,psi_b,m_a,m_b,g_a,g_b,g_ab,N_a,N_b,L_z,Mat_V_a,Mat_V_b,Mat_k_squared,hbar,dx,dy)
Laplacian_a = FFT_Laplacian(psi_a,Mat_k_squared);
Laplacian_b = FFT_Laplacian(psi_b,Mat_k_squared);

Kinetic_a =  -N_a*hbar^2/(2*m_a)*sum(sum(Laplacian_a.*conj(psi_a)));
Kinetic_b =  -N_b*hbar^2/(2*m_b)*sum(sum(Laplacian_b.*conj(psi_b)));

Pot_a = N_a*sum(sum(Mat_V_a .* psi_a .* conj(psi_a)));
Pot_b = N_b*sum(sum(Mat_V_b .* psi_b .* conj(psi_b)));

Intra_a = g_a*N_a^2/(2*L_z)* sum(sum(psi_a .* conj(psi_a) .* psi_a .* conj(psi_a)));
Intra_b = g_b*N_b^2/(2*L_z)* sum(sum(psi_b .* conj(psi_b) .* psi_b .* conj(psi_b)));

Inter_ab = g_ab*N_a*N_b/L_z*sum(sum(psi_a .*conj(psi_a) .* psi_b .* conj(psi_b)));


Ene = dx*dy*real(Kinetic_a + Kinetic_b + Pot_a + Pot_b + Intra_a + Intra_b + Inter_ab);

end

