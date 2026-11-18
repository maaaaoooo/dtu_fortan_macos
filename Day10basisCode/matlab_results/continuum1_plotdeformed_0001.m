 % Plotting Un-Deformed and Deformed Structure
 close all
 clear all
 continuum1_plotdeformed_data0001;
 % Make plot
 figure
 set(gcf,'Name','Deformed')
 xx=zeros(size(IX))';yy=xx;
 subplot(2,1,2)
 patch("Faces", IX(:,1:4), "Vertices", X(:,1:2)+[D(1:2:end)'; D(2:2:end)';]', "FaceColor", [1 1 0], "EdgeColor", "k");
 title('Deformed')
 axis equal;
 xaxes = get(gca,'xlim');
 yaxes = get(gca,'ylim');
 axis off;
 subplot(2,1,1)
 patch("Faces", IX(:,1:4), "Vertices", X(:,1:2), "FaceColor", [1 1 0], "EdgeColor", "k");
 title('Undeformed')
 axis([min(xaxes(1),min(X(:,1))) max(xaxes(2),max(X(:,1)))...
  min(yaxes(1),min(X(:,2))) max(yaxes(2),max(X(:,2))) ]);
 axis equal;
 axis off;
 hold off
 set(gcf,'color',[ 1  1 1]);
