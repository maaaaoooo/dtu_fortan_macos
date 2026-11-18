 % Plotting Element Values
 close all
 clear all
 example1_plotelements_data_0002;
 % Determine colorscale
 colormap('jet')
 % Make plot
 title('Stresses')
 hold on
 cmap = jet(256);
 relval=(plotval-min(plotval))/(max(plotval)-min(plotval));
 XX=X(:,1);
 XY=X(:,2);
 h=plot(XX(IX(:,1:2))',XY(IX(:,1:2))',"LineWidth",1.5);
 for i = 1:length(h)
     set(h(i),"color",cmap(floor(relval(i)*255+1),:));
 end
 clim([min(plotval) max(plotval)]);
 axis equal;
 axis off;
 colorbar;
 hold off
 set(gcf,'color',[ 1  1 1]);
