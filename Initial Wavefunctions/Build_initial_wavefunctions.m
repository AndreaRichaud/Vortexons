function [psi_a,psi_b] = Build_initial_wavefunctions(r_Border,disp_x,disp_y,vortex_charge,N_horiz_points,N_verti_points,vec_x,vec_y,dx,dy)

psi_a=zeros(N_verti_points,N_horiz_points);
psi_b=zeros(N_verti_points,N_horiz_points);

sigma_0_psi_b=1*10^-6;

if disp_x==0 && disp_y==0 %If the vortex is at the center, then no need of the image vortex
    for i=1:N_verti_points
        for j=1:N_horiz_points
            if (vec_x(j))^2+(vec_y(i))^2<(r_Border)^2
                psi_a(i,j)=exp(1i*(+vortex_charge*atan2((vec_y(i)-disp_y),(vec_x(j)-disp_x))  ));
                psi_b(i,j)=exp( -((vec_x(j)-disp_x)^2 + (vec_y(i)-disp_y)^2)/(2*(sigma_0_psi_b)^2));
            else
                psi_a(i,j)=0;
                psi_b(i,j)=0;
            end
        end
    end
else
    % It he vortex is not at the center, than it's better to include the
    % presence of an image vortex
    disp_x_prime=r_Border^2/(disp_x^2+disp_y^2)*disp_x;
    disp_y_prime=r_Border^2/(disp_x^2+disp_y^2)*disp_y;
    for i=1:N_verti_points
        for j=1:N_horiz_points
            if (vec_x(j))^2+(vec_y(i))^2<(r_Border)^2
                psi_a(i,j)=exp(1i*( +vortex_charge*atan2((vec_y(i)-disp_y),(vec_x(j)-disp_x)) - vortex_charge*atan2((vec_y(i)-disp_y_prime),(vec_x(j)-disp_x_prime))   ));
                psi_b(i,j)=exp( -((vec_x(j)-disp_x)^2 + (vec_y(i)-disp_y)^2)/(2*(sigma_0_psi_b)^2));
            else
                psi_a(i,j)=0;
                psi_b(i,j)=0;
            end
        end
    end
end

norm_a = dx*dy*sum(sum(psi_a .* conj(psi_a)));
norm_b = dx*dy*sum(sum(psi_b .* conj(psi_b)));

psi_a=psi_a/sqrt(norm_a);
psi_b=psi_b/sqrt(norm_b);
end

