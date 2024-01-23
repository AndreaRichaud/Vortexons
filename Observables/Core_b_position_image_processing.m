function [x_b,y_b] = Core_b_position_image_processing(psi_b,r_Border,vec_x,vec_y,dx,dy)
% Normalized height at which one recognizes the vortex core
threshold_intensity=0.2;

% Compute the density:
rho_b=real(psi_b.*conj(psi_b));
% Normalize the density in such a way that the maximum value is one, so it
% can be processed by the image-processing tool.
rho_b=rho_b/(max(max(rho_b)));

% Radius of the circumference used to crop the density plot:
coefficiente_ritaglio=0.9;
rho_mask=zeros(size(rho_b,1),size(rho_b,2));
for i=1:size(rho_mask,1)
    for j=1:size(rho_mask,2)
        if sqrt((vec_x(j))^2+(vec_y(i))^2)<coefficiente_ritaglio*r_Border
            rho_mask(i,j)=true;
        else
            rho_mask(i,j)=false;
        end
    end
end

% The target region (i.e. the vortex core) must have a low intensity AND
% lie within the (cropped) circular region.
rho_b_reduced=(rho_b>threshold_intensity)&logical(rho_mask);

s = regionprops(rho_b_reduced,rho_b,{'Centroid','WeightedCentroid'});

x_b_raw=s(1).WeightedCentroid(1);
y_b_raw=s(1).WeightedCentroid(2);

x_b=vec_x(floor(x_b_raw))+dx*(x_b_raw-floor(x_b_raw));
y_b=vec_y(floor(y_b_raw))+dy*(y_b_raw-floor(y_b_raw));


end

