function [x_b,y_b] = Core_b_position(psi_b,Mat_x,Mat_y,dx,dy)

rho_b=real(psi_b.*conj(psi_b));

x_b=dx*dy*sum(sum(Mat_x.*rho_b))/(dx*dy*sum(sum(rho_b)));
y_b=dx*dy*sum(sum(Mat_y.*rho_b))/(dx*dy*sum(sum(rho_b)));

end

