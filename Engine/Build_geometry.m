function [dx,dy,vec_x,vec_y,Mat_x,Mat_y,vec_k_x,vec_k_y,Mat_k_x,Mat_k_y,Mat_k_squared] = Build_geometry(L_x,L_y,N_horiz_points,N_verti_points)

% Spatial steps:
dx=L_x/(N_horiz_points-1);
dy=L_y/(N_verti_points-1);

% Vector with the x-coordinates (cenetered in x=0)
vec_x = linspace(-L_x/2,L_x/2,N_horiz_points);
% Vector with the y-coordinates (cenetered in y=0)
vec_y = linspace(-L_y/2,L_y/2,N_verti_points);

% Matrix with the x-coordinates (each row is vector vec_x)
Mat_x = zeros(N_verti_points,N_horiz_points);
% Matrix with the y-coordinates (each column is vector vec_y)
Mat_y = zeros(N_verti_points,N_horiz_points);

for i=1:N_verti_points
    Mat_x(i,:)=vec_x;
end

for j=1:N_horiz_points
    Mat_y(:,j)=vec_y;
end

% Fourier space vectors
vec_k_x = 2*pi/L_x*(-N_horiz_points/2:N_horiz_points/2-1);
vec_k_y = 2*pi/L_y*(-N_verti_points/2:N_verti_points/2-1);

[Mat_k_x, Mat_k_y] = meshgrid(vec_k_x, vec_k_y);

Mat_k_squared = Mat_k_x.^2 + Mat_k_y.^2;

end

