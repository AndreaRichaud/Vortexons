function [vec_Ene,vec_L_z,vec_position_a,vec_position_b] = Initialize_Output(n_iterations,sample_frequency)
n_elements=floor(n_iterations/sample_frequency);

vec_Ene=zeros(1,n_elements);
vec_L_z=zeros(1,n_elements);
vec_position_a=zeros(n_elements,2);
vec_position_b=zeros(n_elements,2);

end

