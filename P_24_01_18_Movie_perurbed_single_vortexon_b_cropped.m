clear all
close all
clc

tic

% Determine where your m-file's folder is.
folder = fileparts(which(mfilename)); 
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));

mkdir Temp_real

load('Single_vortexon.mat')
filename=sprintf('Vortexon_dynamics_%s',datestr(now,'yyyy_mm_dd_HH_MM_SS'));


% Do you want a pinning potential?
Pinning_flag=0;
[Mat_V_a, Mat_V_b]=Build_confining_potential_b_cropped(r_Border,disp_x,disp_y,N_horiz_points,N_verti_points,vec_x,vec_y,V_0,Pinning_flag,V_pinning,sigma_pinning);

c_p=1;

n_iterations=2*10^5;
sample_frequency=1000;

[vec_Ene,vec_L_z,vec_position_a,vec_position_b]=Initialize_Output(n_iterations,sample_frequency);

for i=1:n_iterations
   
    
    [psi_a,psi_b]=Real_time_evolve(psi_a,psi_b,m_a,m_b,g_a,g_b,g_ab,N_a,N_b,L_z,Mat_V_a,Mat_V_b,dx,dy,Mat_k_squared,hbar,dt);
    
    if mod(i,sample_frequency)==1
        
        fprintf('Avanzamento = %.2f %% in %.0f s \n', i/n_iterations*100,toc)
        
        %Fig 1:
        fig=Plot_frame_off(psi_a,psi_b,vec_x,vec_y,Mat_x,Mat_y,Mat_k_x,Mat_k_y,hbar,r_Border,dt,i);
        fig_name = sprintf('Temp_real/%04d.jpg',c_p);
        saveas(fig, fig_name);
        close
       
        vec_Ene(c_p)=Energy_computator(psi_a,psi_b,m_a,m_b,g_a,g_b,g_ab,N_a,N_b,L_z,Mat_V_a,Mat_V_b,Mat_k_squared,hbar,dx,dy);
        vec_L_z(c_p)=Total_angular_momentum(N_a,N_b,psi_a,psi_b,Mat_x,Mat_y,Mat_k_x,Mat_k_y,hbar,dx,dy);
        [vec_position_a(c_p,1),vec_position_a(c_p,2)] = Core_a_position(psi_a,r_Border,vec_x,vec_y,dx,dy);
        [vec_position_b(c_p,1),vec_position_b(c_p,2)] = Core_b_position(psi_b,Mat_x,Mat_y,dx,dy);
        
        c_p=c_p+1;
        
        
    end
    
end

c_p=c_p-1;
clf
%Fig 2:
Plot_trajecory(vec_position_a,vec_position_b,c_p,r_Border,filename,2);
%Fig 3:
Plot_trajecory_split(vec_position_a,vec_position_b,c_p,r_Border,filename,3);
%Fig 4:
Plot_Ene_and_L_z(vec_Ene,vec_L_z,c_p,filename,4);

Produce_animation(filename)


save(sprintf('%s.mat',filename),'vec_position_a', 'vec_position_b')

