function [] = Plot_trajecory(vec_position_a,vec_position_b,c_p,r_Border,filename,n_fig)

fig_Trajectory=figure(n_fig);
hold on
scatter(vec_position_a(1:c_p,1)/r_Border,vec_position_a(1:c_p,2)/r_Border,80,'.b')
scatter(vec_position_b(1:c_p,1)/r_Border,vec_position_b(1:c_p,2)/r_Border,50,'.r')
hold off
xlabel('$x/R$','interpreter','latex')
ylabel('$y/R$','interpreter','latex')
legend('A','B','Location','Best');
xlim([-1.1,1.1])
ylim([-1.1,1.1])
viscircles([0,0],1,'Color','k');
daspect([1 1 1])
grid on
grid minor

saveas(fig_Trajectory,sprintf('%s_Trajectory',filename),'png');

end

