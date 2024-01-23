function [fig] = Plot_frame_off(psi_a,psi_b,vec_x,vec_y,Mat_x,Mat_y,Mat_k_x,Mat_k_y,hbar,r_Border,dt,i)

[Gradient_x_a,Gradient_y_a,]=FFT_Derivative(psi_a,Mat_k_x,Mat_k_y);
[Gradient_x_b,Gradient_y_b,]=FFT_Derivative(psi_b,Mat_k_x,Mat_k_y);
Momentum_x_a = real(-1i*hbar*conj(psi_a).*Gradient_x_a);
Momentum_y_a = real(-1i*hbar*conj(psi_a).*Gradient_y_a);

Momentum_x_b = real(-1i*hbar*conj(psi_b).*Gradient_x_b);
Momentum_y_b = real(-1i*hbar*conj(psi_b).*Gradient_y_b);

fig=figure('Visible', 'off', 'Position', [100, 100, 1920, 1080]);
ax(1)=subplot(2,3,1);
surf(Mat_x/r_Border,Mat_y/r_Border,real(psi_a.*conj(psi_a)),'EdgeColor','none')
colormap(ax(1),gray)
view(2)
xlabel('$x/R$','interpreter','latex')
ylabel('$y/R$','interpreter','latex')
title('$|\psi_a|^2$','interpreter','latex')
pbaspect([1 1 1])
set(gca,'TickLabelInterpreter','latex','FontSize',14)
xlim([min(vec_x)/r_Border max(vec_x)/r_Border])
ylim([min(vec_y)/r_Border max(vec_y)/r_Border])
%
ax(2)=subplot(2,3,4);
surf(Mat_x/r_Border,Mat_y/r_Border,real(psi_b.*conj(psi_b)),'EdgeColor','none')
colormap(ax(2),gray)
view(2)
xlabel('$x/R$','interpreter','latex')
ylabel('$y/R$','interpreter','latex')
title('$|\psi_b|^2$','interpreter','latex')
pbaspect([1 1 1])
set(gca,'TickLabelInterpreter','latex','FontSize',14)
xlim([min(vec_x)/r_Border max(vec_x)/r_Border])
ylim([min(vec_y)/r_Border max(vec_y)/r_Border])


subplot(2,3,2)
quiver(Mat_x(1:7:end,1:7:end)/r_Border,Mat_y(1:7:end,1:7:end)/r_Border,Momentum_x_a(1:7:end,1:7:end),Momentum_y_a(1:7:end,1:7:end))
xlim([min(vec_x)/r_Border max(vec_x)/r_Border])
ylim([min(vec_y)/r_Border max(vec_y)/r_Border])
xlabel('$x/R$','interpreter','latex')
ylabel('$y/R$','interpreter','latex')
title('$\vec{J}_a$','interpreter','latex')
pbaspect([1 1 1])
set(gca,'TickLabelInterpreter','latex','FontSize',14)
subplot(2,3,5)
quiver(Mat_x(1:7:end,1:7:end)/r_Border,Mat_y(1:7:end,1:7:end)/r_Border,Momentum_x_b(1:7:end,1:7:end),Momentum_y_b(1:7:end,1:7:end))
xlim([min(vec_x)/r_Border max(vec_x)/r_Border])
ylim([min(vec_y)/r_Border max(vec_y)/r_Border])
xlabel('$x/R$','interpreter','latex')
ylabel('$y/R$','interpreter','latex')
title('$\vec{J}_b$','interpreter','latex')
pbaspect([1 1 1])
set(gca,'TickLabelInterpreter','latex','FontSize',14)


ax(5)=subplot(2,3,3);
hold on
surf(Mat_x/r_Border,Mat_y/r_Border,0*psi_a,angle(psi_a),'EdgeColor','none')
contour(Mat_x/r_Border,Mat_y/r_Border,angle(psi_a),[-pi:2*pi/10:pi],'k')
hold off
colormap(ax(5),hsv)
view(2)
xlabel('$x/R$','interpreter','latex')
ylabel('$y/R$','interpreter','latex')
title('$\theta_a$','interpreter','latex')
pbaspect([1 1 1])
set(gca,'TickLabelInterpreter','latex','FontSize',14)
xlim([min(vec_x)/r_Border max(vec_x)/r_Border])
ylim([min(vec_y)/r_Border max(vec_y)/r_Border])
%
ax(6)=subplot(2,3,6);
surf(Mat_x/r_Border,Mat_y/r_Border,angle(psi_b),'EdgeColor','none')
colormap(ax(6),hsv)
view(2)
xlabel('$x/R$','interpreter','latex')
ylabel('$y/R$','interpreter','latex')
title('$\theta_b$','interpreter','latex')
pbaspect([1 1 1])
set(gca,'TickLabelInterpreter','latex','FontSize',14)
xlim([min(vec_x)/r_Border max(vec_x)/r_Border])
ylim([min(vec_y)/r_Border max(vec_y)/r_Border])

sgtitle(sprintf('t = %.3f s', dt*i))

pause(0.1)
end

