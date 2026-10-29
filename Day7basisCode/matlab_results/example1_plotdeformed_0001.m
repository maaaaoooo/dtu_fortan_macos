 % Plotting Un-Deformed and Deformed Structure
 close all
 clear all
 example1_plotdeformed_data0001;
 % Make plot
 figure
 set(gcf,'Name','Deformed')
 subplot(2,1,2)
     XX=X(:,1)+D(1:2:end);
     XY=X(:,2)+D(2:2:end);
     plot(XX(IX(:,1:2))',XY(IX(:,1:2))',"LineWidth",1.5,"color","b")
 title('Deformed')
 axis equal
 xaxes = get(gca,'xlim');
 yaxes = get(gca,'ylim');
 axis off;
 subplot(2,1,1)
     XX=X(:,1);
     XY=X(:,2);
     plot(XX(IX(:,1:2))',XY(IX(:,1:2))',"LineWidth",1.5,"color","b")
 title('Undeformed')
 axis([min(xaxes(1),min(X(:,1))) max(xaxes(2),max(X(:,1)))...
  min(yaxes(1),min(X(:,2))) max(yaxes(2),max(X(:,2))) ]);
 axis off;
 hold off
 set(gcf,'color',[ 1  1 1]);
