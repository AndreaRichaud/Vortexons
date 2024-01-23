function [Mat_V_a, Mat_V_b]=Build_confining_potential_slope(r_Border,disp_x,disp_y,N_horiz_points,N_verti_points,vec_x,vec_y,V_0,time)

assert(disp_x^2+disp_y^2<r_Border^2,'Error: the position of the vortex lies outside the circular trap')

pendenza_V_b=10^-28/0.2;

% Build the flat confining potential (which is of course 2D):
Mat_V_a=zeros(N_verti_points,N_horiz_points);
Mat_V_b=zeros(N_verti_points,N_horiz_points);


for i=1:N_verti_points
    for j=1:N_horiz_points
        if (vec_x(j))^2+(vec_y(i))^2<(r_Border)^2
            Mat_V_a(i,j)=0;
            Mat_V_b(i,j)=vec_y(i)*pendenza_V_b*time;
        else
            Mat_V_a(i,j)=V_0;
            Mat_V_b(i,j)=vec_y(i)*pendenza_V_b*time;
        end
    end
end




end

