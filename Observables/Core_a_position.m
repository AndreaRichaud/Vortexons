function [x_a,y_a] = Core_a_position(psi_a,r_Border,vec_x,vec_y,dx,dy)
% Normalized height at which one recognizes the vortex core
threshold_intensity=0.7;

% Compute the density:
rho_a=real(psi_a.*conj(psi_a));
% Normalize the density in such a way that the maximum value is one, so it
% can be processed by the image-processing tool.
rho_a=rho_a/(max(max(rho_a)));

rho_a=1-rho_a;

% Radius of the circumference used to crop the density plot:
coefficiente_ritaglio=0.9;
rho_mask=zeros(size(rho_a,1),size(rho_a,2));
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
rho_a_reduced=(rho_a>threshold_intensity)&logical(rho_mask);

s = regionprops(rho_a_reduced,rho_a,{'Centroid','WeightedCentroid'});

x_a_raw=s(1).WeightedCentroid(1);
y_a_raw=s(1).WeightedCentroid(2);

x_a=vec_x(floor(x_a_raw))+dx*(x_a_raw-floor(x_a_raw));
y_a=vec_y(floor(y_a_raw))+dy*(y_a_raw-floor(y_a_raw));


end

