function plot_contours(I,J,x,y,p,c,tht,Vx,Vy,levels,name)

figure;
contourf(x,y,p,levels); colorbar; colormap(flipud(parula)); 
xlabel('$x$','Interpreter','latex'); ylabel('$z$','Interpreter','latex');
title('$\psi(x,z,t)$','Interpreter','latex'); 
xlabel('$\mathbf{x}$','Interpreter','latex'); ylabel('$\mathbf{z}$','Interpreter','latex');
title('$\mathbf{\Psi(x,z,t)}$','Interpreter','latex'); 

figure;
contourf(x,y,c,levels); colorbar; colormap(flipud(parula)); 
xlabel('$x$','Interpreter','latex'); ylabel('$z$','Interpreter','latex');
title('$c(x,z,t)$','Interpreter','latex'); 
xlabel('$\mathbf{x}$','Interpreter','latex'); ylabel('$\mathbf{z}$','Interpreter','latex');
title('$\mathbf{c(x,z,t)}$','Interpreter','latex'); 

figure;
contourf(x,y,tht,levels); colorbar; colormap(flipud(parula)); 
xlabel('$x$','Interpreter','latex'); ylabel('$z$','Interpreter','latex');
title('$\theta(x,z,t)$','Interpreter','latex'); 

figure
contourf(x,y,Vx,levels); colorbar; colormap(flipud(parula)); 
xlabel('$x$','Interpreter','latex'); ylabel('$z$','Interpreter','latex'); 
title('$q_x(x,z)$','Interpreter','latex'); 

figure
contourf(x,y,Vy,levels); colorbar; colormap(flipud(parula)); 
xlabel('$x$','Interpreter','latex'); ylabel('$z$','Interpreter','latex'); 
title('$q_z(x,z)$','Interpreter','latex'); 
