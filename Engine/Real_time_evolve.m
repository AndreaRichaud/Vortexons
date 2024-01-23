function [psi_a,psi_b]=Real_time_evolve(psi_a,psi_b,m_a,m_b,g_a,g_b,g_ab,N_a,N_b,L_z,Mat_V_a,Mat_V_b,dx,dy,Mat_k_squared,hbar,dt)

[Hpsi_a,Hpsi_b] = Real_time_iteration(psi_a,psi_b,m_a,m_b,g_a,g_b,g_ab,N_a,N_b,L_z,Mat_V_a,Mat_V_b,Mat_k_squared,hbar);

k_1_a=-1i*dt/hbar*(Hpsi_a);
k_1_b=-1i*dt/hbar*(Hpsi_b);

[Hpsi_a,Hpsi_b] = Real_time_iteration(psi_a+k_1_a/2,psi_b+k_1_b/2,m_a,m_b,g_a,g_b,g_ab,N_a,N_b,L_z,Mat_V_a,Mat_V_b,Mat_k_squared,hbar);

k_2_a=-1i*dt/hbar*(Hpsi_a);
k_2_b=-1i*dt/hbar*(Hpsi_b);

[Hpsi_a,Hpsi_b] = Real_time_iteration(psi_a+k_2_a/2,psi_b+k_2_b/2,m_a,m_b,g_a,g_b,g_ab,N_a,N_b,L_z,Mat_V_a,Mat_V_b,Mat_k_squared,hbar);

k_3_a=-1i*dt/hbar*(Hpsi_a);
k_3_b=-1i*dt/hbar*(Hpsi_b);

[Hpsi_a,Hpsi_b] = Real_time_iteration(psi_a+k_3_a,psi_b+k_3_b,m_a,m_b,g_a,g_b,g_ab,N_a,N_b,L_z,Mat_V_a,Mat_V_b,Mat_k_squared,hbar);

k_4_a=-1i*dt/hbar*(Hpsi_a);
k_4_b=-1i*dt/hbar*(Hpsi_b);

psi_a = psi_a+ 1/6*(k_1_a+2*k_2_a+2*k_3_a+k_4_a);
psi_b = psi_b+ 1/6*(k_1_b+2*k_2_b+2*k_3_b+k_4_b);

norm_a = dx*dy*sum(sum(psi_a .* conj(psi_a)));
norm_b = dx*dy*sum(sum(psi_b .* conj(psi_b)));

psi_a=psi_a/sqrt(norm_a);
psi_b=psi_b/sqrt(norm_b);

end

