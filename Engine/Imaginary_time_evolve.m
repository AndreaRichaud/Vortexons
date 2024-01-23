function [psi_a,psi_b,Hpsi_a,Hpsi_b,L_z_psi_a,L_z_psi_b]=Imaginary_time_evolve(psi_a,psi_b,m_a,m_b,g_a,g_b,g_ab,N_a,N_b,L_z,Mat_V_a,Mat_V_b,Mat_x,Mat_y,dx,dy,Mat_k_x,Mat_k_y,hbar,dt,Omega)

[Hpsi_a,Hpsi_b,L_z_psi_a,L_z_psi_b] = Imaginary_time_iteration(psi_a,psi_b,m_a,m_b,g_a,g_b,g_ab,N_a,N_b,L_z,Mat_V_a,Mat_V_b,Mat_x,Mat_y,Mat_k_x,Mat_k_y,hbar);

k_1_a=-dt/hbar*(Hpsi_a-Omega*L_z_psi_a);
k_1_b=-dt/hbar*(Hpsi_b-Omega*L_z_psi_b);

[Hpsi_a,Hpsi_b,L_z_psi_a,L_z_psi_b] = Imaginary_time_iteration(psi_a+k_1_a/2,psi_b+k_1_b/2,m_a,m_b,g_a,g_b,g_ab,N_a,N_b,L_z,Mat_V_a,Mat_V_b,Mat_x,Mat_y,Mat_k_x,Mat_k_y,hbar);

k_2_a=-dt/hbar*(Hpsi_a-Omega*L_z_psi_a);
k_2_b=-dt/hbar*(Hpsi_b-Omega*L_z_psi_b);

[Hpsi_a,Hpsi_b,L_z_psi_a,L_z_psi_b] = Imaginary_time_iteration(psi_a+k_2_a/2,psi_b+k_2_b/2,m_a,m_b,g_a,g_b,g_ab,N_a,N_b,L_z,Mat_V_a,Mat_V_b,Mat_x,Mat_y,Mat_k_x,Mat_k_y,hbar);

k_3_a=-dt/hbar*(Hpsi_a-Omega*L_z_psi_a);
k_3_b=-dt/hbar*(Hpsi_b-Omega*L_z_psi_b);

[Hpsi_a,Hpsi_b,L_z_psi_a,L_z_psi_b] = Imaginary_time_iteration(psi_a+k_3_a,psi_b+k_3_b,m_a,m_b,g_a,g_b,g_ab,N_a,N_b,L_z,Mat_V_a,Mat_V_b,Mat_x,Mat_y,Mat_k_x,Mat_k_y,hbar);

k_4_a=-dt/hbar*(Hpsi_a-Omega*L_z_psi_a);
k_4_b=-dt/hbar*(Hpsi_b-Omega*L_z_psi_b);

psi_a = psi_a+ 1/6*(k_1_a+2*k_2_a+2*k_3_a+k_4_a);
psi_b = psi_b+ 1/6*(k_1_b+2*k_2_b+2*k_3_b+k_4_b);

norm_a = dx*dy*sum(sum(psi_a .* conj(psi_a)));
norm_b = dx*dy*sum(sum(psi_b .* conj(psi_b)));

psi_a=psi_a/sqrt(norm_a);
psi_b=psi_b/sqrt(norm_b);

end

