function Omega = Omega_minus(hbar,N_a,N_b,m_a,m_b,disp_x,disp_y,r_Border)
%This function ouptuts, in S.I. units Eq. (33) of our PRA 2021
mu=N_b*m_b/(N_a*m_a);
r0=sqrt(disp_x^2+disp_y^2)/r_Border;
Omega=hbar/(m_a*r_Border^2)*(2/(1-r0^2))/(1+sqrt(1-2*mu/(1-r0^2)));
end

