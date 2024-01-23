function [psi_b] = Quadrupolar_perturbation(psi_b,x0,y0,alpha,N_verti_points,N_horiz_points,vec_x,vec_y,dx,dy)

Mat_perturbation=zeros(N_verti_points,N_horiz_points);

for i=1:N_verti_points
    for j=1:N_horiz_points
        Mat_perturbation(i,j)=exp(1i/alpha^2 *((vec_x(j)-x0)^2 - (vec_y(i)-y0)^2));
    end
end

psi_b=psi_b.*Mat_perturbation;
norm_b = dx*dy*real(sum(sum(psi_b .* conj(psi_b))));
psi_b  = psi_b/sqrt(norm_b);

end

