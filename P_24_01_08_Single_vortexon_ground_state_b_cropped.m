clear all
close all
clc

tic

% Determine where your m-file's folder is.
folder = fileparts(which(mfilename)); 
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));

mkdir Temp_img

%Planck's constant
hbar = 1.0545718*10^-34;
%Atomic mass unit
amu = 1.66054*10^-27;
%Bohr radius:
a0 = 5.291777*10^-11;

%Sodium-23 mass:
m_a = 23*amu;
%Potassium-39 mass:
m_b = 39*amu;

% Number of atoms in species-a
N_a = 90000;
% Number of atoms in species-b
N_b = 2000;

% Radius of the circular trap
r_Border=50*10^-6;

% Size of the simulated area:
L_x = 120*10^-6;
L_y = 120*10^-6;
% Direction z is assumed to be compactified:
L_z = 2*10^-6;

% Intra- and Inter-species interactions
a_aa = 52.0*a0;
a_bb =  7.6*a0;
%a_ab = 14.0*a0; 
 
% Intra- and Inter-species interactions
g_a  = 4*pi*hbar^2*a_aa/m_a;
g_b  = 4*pi*hbar^2*a_bb/m_b;
%m_ab=(1/m_a+1/m_b)^-1;
%g_ab=2*pi*hbar^2*a_ab/m_ab;
g_ab = 2*sqrt(g_a*g_b);

% Number of iterations:
n_iterations = 1*10^4;
% Time step:
dt = 10^-5;

% Number of points in the spatial grid (make an integer power of 2 in order to
% optimize the FFT algorithms):
N_horiz_points = 256;
N_verti_points = N_horiz_points;


%% Build geometry:
[dx,dy,vec_x,vec_y,Mat_x,Mat_y,vec_k_x,vec_k_y,Mat_k_x,Mat_k_y,Mat_k_squared] = Build_geometry(L_x,L_y,N_horiz_points,N_verti_points);

%% Vortex properties:
% Vortex charge:
vortex_charge=+1;
% Vortex center:
disp_x=0.57*r_Border;
disp_y=0;
% Do you want a pinning potential?
Pinning_flag=1;

% Rotation frequency of the trap:
Omega=0;

%% Hard_wall_confining_potential properties
V_0=10^-30;
V_pinning = V_0;
sigma_pinning = 0.1*10^-6;

%% Build confining potential
[Mat_V_a, Mat_V_b]=Build_confining_potential_b_cropped(r_Border,disp_x,disp_y,N_horiz_points,N_verti_points,vec_x,vec_y,V_0,Pinning_flag,V_pinning,sigma_pinning);

%% Build initial wavefunctions
[psi_a,psi_b] = Build_initial_wavefunctions(r_Border,disp_x,disp_y,vortex_charge,N_horiz_points,N_verti_points,vec_x,vec_y,dx,dy);

% Define sample frequency:
sample_frequency=1000;

%Initialize some internal paramters
vec_mu_a=zeros(floor(n_iterations/sample_frequency),1);
vec_mu_b=zeros(floor(n_iterations/sample_frequency),1);
vec_Ene=zeros(floor(n_iterations/sample_frequency),1);

i_sampling=1;
for i=1:n_iterations
    
    [psi_a,psi_b,Hpsi_a,Hpsi_b,L_z_psi_a,L_z_psi_b]=Imaginary_time_evolve(psi_a,psi_b,m_a,m_b,g_a,g_b,g_ab,N_a,N_b,L_z,Mat_V_a,Mat_V_b,Mat_x,Mat_y,dx,dy,Mat_k_x,Mat_k_y,hbar,dt,Omega);
    
    if mod(i,sample_frequency)==1

        fprintf('Avanzamento = %.2f %% in %.0f s \n', i/n_iterations*100,toc)

        vec_mu_a(i_sampling)=real(dx*dy*sum(sum(conj(psi_a).*(Hpsi_a-Omega*L_z_psi_a))));
        vec_mu_b(i_sampling)=real(dx*dy*sum(sum(conj(psi_b).*(Hpsi_b-Omega*L_z_psi_b))));
        vec_Ene(i_sampling)=Energy_computator(psi_a,psi_b,m_a,m_b,g_a,g_b,g_ab,N_a,N_b,L_z,Mat_V_a,Mat_V_b,Mat_k_squared,hbar,dx,dy);

        fig=Plot_frame_off(psi_a,psi_b,vec_x,vec_y,Mat_x,Mat_y,Mat_k_x,Mat_k_y,hbar,r_Border,dt,i);
        fig_name = sprintf('Temp_img/%d.jpg',i_sampling);
        saveas(fig, fig_name);
        close

        i_sampling=i_sampling+1;
        
    end
    
end

label_fig=2;
Plot_convergence_imaginary_time(vec_mu_a,vec_mu_b,vec_Ene,i_sampling-1,sample_frequency,label_fig)        

save('Single_vortexon.mat')


