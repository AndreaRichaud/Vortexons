function [Mat_V_a, Mat_V_b]=Build_confining_potential(r_Border,disp_x,disp_y,N_horiz_points,N_verti_points,vec_x,vec_y,V_0,Pinning_flag,V_pinning,sigma_pinning)

assert(disp_x^2+disp_y^2<r_Border^2,'Error: the position of the vortex lies outside the circular trap')

% Build the flat confining potential (which is of course 2D):
Mat_V_a=zeros(N_verti_points,N_horiz_points);
Mat_V_b=zeros(N_verti_points,N_horiz_points);

if Pinning_flag==0
    
    for i=1:N_verti_points
        for j=1:N_horiz_points
            if (vec_x(j))^2+(vec_y(i))^2<(r_Border)^2
                Mat_V_a(i,j)=0;
                Mat_V_b(i,j)=0;
            else
                Mat_V_a(i,j)=V_0;
                Mat_V_b(i,j)=V_0;
            end
        end
    end
    
elseif Pinning_flag==1
    
    for i=1:N_verti_points
        for j=1:N_horiz_points
            if (vec_x(j))^2+(vec_y(i))^2<(r_Border)^2
                Mat_V_a(i,j)=V_pinning*exp(- ((vec_x(j)-disp_x)^2 + (vec_y(i)-disp_y)^2)/(2*(sigma_pinning)^2));
                Mat_V_b(i,j)=0;
            else
                Mat_V_a(i,j)=V_0;
                Mat_V_b(i,j)=V_0;
            end
        end
    end
    
    
end

