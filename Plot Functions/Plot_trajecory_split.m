function [] = Plot_trajecory_split(vec_position_a,vec_position_b,c_p,r_Border,filename,n_fig)

fig_Trajectory=figure(n_fig);
subplot(1,2,1)
scatter(vec_position_a(1:c_p,1)/r_Border,vec_position_a(1:c_p,2)/r_Border,10,'.b')
xlabel('$x/R$','interpreter','latex')
ylabel('$y/R$','interpreter','latex')
xlim([-1.1,1.1])
ylim([-1.1,1.1])
viscircles([0,0],1,'Color','k');
daspect([1 1 1])
grid on
grid minor
%
subplot(1,2,2)
scatter(vec_position_b(1:c_p,1)/r_Border,vec_position_b(1:c_p,2)/r_Border,10,'.r')
xlabel('$x/R$','interpreter','latex')
ylabel('$y/R$','interpreter','latex')
xlim([-1.1,1.1])
ylim([-1.1,1.1])
viscircles([0,0],1,'Color','k');
daspect([1 1 1])
grid on
grid minor


saveas(fig_Trajectory,sprintf('%s_Trajectory_split',filename),'png');

end

