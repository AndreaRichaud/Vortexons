function psi_b = Kick_psi_b(psi_b,v_x,v_y,Mat_x,Mat_y,hbar,m_b,dx,dy)


alpha_x = m_b*v_x/hbar;
alpha_y = m_b*v_y/hbar;

psi_b = psi_b .* exp(1i*(Mat_x*alpha_x + Mat_y*alpha_y));


norm_b = dx*dy*real(sum(sum(psi_b .* conj(psi_b))));
psi_b  = psi_b/sqrt(norm_b);


end

