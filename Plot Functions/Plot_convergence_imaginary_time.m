function [] = Plot_convergence_imaginary_time(vec_mu_a,vec_mu_b,vec_Ene,i_sampling,sample_frequency,label_fig)

figure(label_fig)

subplot(1,3,1)
plot(vec_mu_a(1:i_sampling));
xlabel('n. iteration / sample frequency','interpreter','latex')
ylabel('$\mu_a$','interpreter','latex')
%
subplot(1,3,2)
plot(vec_mu_b(1:i_sampling));
xlabel('n. iteration / sample frequency','interpreter','latex')
ylabel('$\mu_b$','interpreter','latex')
%
subplot(1,3,3)
plot(vec_Ene(1:i_sampling));
xlabel('n. iteration / sample frequency','interpreter','latex')
ylabel('$E_{tot}$','interpreter','latex')

pause(0.1)
end

