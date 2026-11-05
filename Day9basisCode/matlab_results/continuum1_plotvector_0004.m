 % Plotting Vector Field, i.e. principle stresses
 close all
 clear all
 continuum1_plotvector_data0004;
 % Make plot
 ne = size(IX,1);
 % User scale parameter: See fedata - scale_vec
 scale_vec = 1;
 % --- Characteristic length (max corner-to-corner extent in x and y) ---
 % max_{i,j} |x_i - x_j| over the 4 corners equals max(x) - min(x)
 xe  = reshape(X(IX(:,1:4),1), ne, 4);
 ye  = reshape(X(IX(:,1:4),2), ne, 4);
 clx = max(max(xe,[],2) - min(xe,[],2));
 cly = max(max(ye,[],2) - min(ye,[],2));
 clmax = max(clx, cly);
 scal  = max(max(abs(vect(:,1:2))))*sqrt(10)/clmax / scale_vec;
 % --- Element centres (average of all element nodes) ---
 xc = mean(reshape(X(IX,1), size(IX)), 2);
 yc = mean(reshape(X(IX,2), size(IX)), 2);
 % --- Principal directions ---
 th = -vect(:,3);
 d1 = [ cos(th)  sin(th)];      % direction of vect_1
 d2 = [-sin(th)  cos(th)];      % direction of vect_2 (rotated +pi/2)
 U1 = (abs(vect(:,1))/scal) .* d1;
 U2 = (abs(vect(:,2))/scal) .* d2;
 % Stack all four arrows per element: +v1, -v1, +v2, -v2
 P   = repmat([xc yc], 4, 1);
 U   = [U1; -U1; U2; -U2];
 sgn = [vect(:,1); vect(:,1); vect(:,2); vect(:,2)];
 pos = sgn >= 0;
 % Make plot
 figure
 hold on
 % Mesh: one patch object for all elements
 patch("Faces", IX(:,1:4), "Vertices", X(:,1:2), ...
 "FaceColor", "w", "EdgeColor", "k");
 % Arrows: one quiver object per colour (blue = tension, red = compression)
 quiver(P(pos,1),  P(pos,2),  U(pos,1),  U(pos,2),  0, "b","ShowArrowHead","off");
 quiver(P(~pos,1), P(~pos,2), U(~pos,1), U(~pos,2), 0, "r","ShowArrowHead","off");
 title( 'Principal stresses')
 axis equal;  axis off;  hold off
 set(gcf,'color',[ 1  1  1]);
