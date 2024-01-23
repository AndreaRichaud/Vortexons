function [] = Plot_Ene_and_L_z(vec_Ene,vec_L_z,c_p,filename,n_fig)
fig_Ene_L_z=figure(n_fig);
subplot(1,2,1)
plot(vec_Ene(1:c_p))
xlabel('n. iterations / sample frequency','interpreter','latex')
ylabel('$E_{tot}$','interpreter','latex')
subplot(1,2,2)
plot(vec_L_z(1:c_p))
xlabel('n. iterations / sample frequency','interpreter','latex')
ylabel('$L_{z,tot}$','interpreter','latex')
saveas(fig_Ene_L_z,sprintf('%s_Ene_L_z',filename),'png');
end

