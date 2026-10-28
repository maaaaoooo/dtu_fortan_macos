 % Plotting Element Values
 close all
 clear all
 cant_5x1_plotelements_data_0002;
 % Determine colorscale
 colormap('jet')
 % Make plot
 title('Shear stress')
 hold on
 patch("Faces", IX(:,1:4), "Vertices", X(:,1:2), "CData", plotval, "FaceColor", "flat", "EdgeColor", "k");
 axis equal;
 axis off;
 colorbar;
 hold off
 set(gcf,'color',[ 1  1 1]);
