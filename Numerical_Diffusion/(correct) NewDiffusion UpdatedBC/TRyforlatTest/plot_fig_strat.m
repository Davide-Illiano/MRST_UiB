function plot_fig_strat(I,J,x,y,p,c,tht,Vx,Vy)

figure;
mesh(x,y,p); %,'mo--','LineWidth',1.2,'MarkerFaceColor','yellow','MarkerSize',7.2) ;
xlabel('$x$','Interpreter','latex'); ylabel('$z$','Interpreter','latex');
zlabel('$\psi(x,z,t)$','Interpreter','latex'); view(115,15); %(20,50);
grid on ;

figure;
mesh(x,y,c); %,'mo--','LineWidth',1.2,'MarkerFaceColor','yellow','MarkerSize',7.2) ;
xlabel('$x$','Interpreter','latex'); ylabel('$z$','Interpreter','latex');
zlabel('$c(x,z,t)$','Interpreter','latex'); view(115,15); %(20,50);
grid on ;
figure;
mesh(x,y,tht); %,'mo--','LineWidth',1.2,'MarkerFaceColor','yellow','MarkerSize',7.2) ;
xlabel('$x$','Interpreter','latex'); ylabel('$z$','Interpreter','latex');
zlabel('$\theta(x,z,t)$','Interpreter','latex'); view(115,15); %view(-45,25); %(20,50);
grid on ;
figure
mesh(x,y,Vx)
xlabel('$x$','Interpreter','latex'); ylabel('$z$','Interpreter','latex'); 
zlabel('$q_x(x,z)$','Interpreter','latex'); view(115,15);
figure
mesh(x,y,Vy)
xlabel('$x$','Interpreter','latex'); ylabel('$z$','Interpreter','latex'); 
zlabel('$q_z(x,z)$','Interpreter','latex'); view(115,15);


