function psi_out = Horizontal_translation_psi(psi,N_horiz_points,N_verti_points,offset_x)

psi_out=zeros(N_verti_points,N_horiz_points);

psi_out(:,offset_x:end)=psi(:,1:end-offset_x+1);

end

