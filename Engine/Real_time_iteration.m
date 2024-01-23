function [Hpsi_a,Hpsi_b] = Real_time_iteration(psi_a,psi_b,m_a,m_b,g_a,g_b,g_ab,N_a,N_b,L_z,Mat_V_a,Mat_V_b,Mat_k_squared,hbar)
Laplacian_a = FFT_Laplacian(psi_a,Mat_k_squared);
Laplacian_b = FFT_Laplacian(psi_b,Mat_k_squared);

Intra_a = g_a*N_a/L_z* psi_a.*conj(psi_a) .* psi_a;
Intra_b = g_b*N_b/L_z* psi_b.*conj(psi_b) .* psi_b;

Inter_ab = g_ab*N_b/L_z* psi_b.*conj(psi_b) .* psi_a;
Inter_ba = g_ab*N_a/L_z* psi_a.*conj(psi_a) .* psi_b;

Pot_a = Mat_V_a .* psi_a;
Pot_b = Mat_V_b .* psi_b;


Hpsi_a = -hbar^2/(2*m_a)*Laplacian_a + Intra_a + Inter_ab + Pot_a;
Hpsi_b = -hbar^2/(2*m_b)*Laplacian_b + Intra_b + Inter_ba + Pot_b;


end

