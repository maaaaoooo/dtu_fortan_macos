 % Plotting Eigenmodes
 close all
 clear all
 beam_modal_ploteig_data_0002;
 freq =   0.23575573060796046      ;
 timeint =   0.88837505539832540      ;
 timetot =    53.302503323899522      ;
 % Find max window size.
 lxmin = min(X(:,1));        lxmax = max(X(:,1));
 lymin = min(X(:,2));        lymax = max(X(:,2));
 dxmin = min(D(1:2:end));    dxmax = max(D(1:2:end));
 dymin = min(D(2:2:end));    dymax = max(D(2:2:end));
 lxmin = lxmin - max(abs(dxmin),abs(dxmax))*1.05;
 lxmax = lxmax + max(abs(dxmin),abs(dxmax))*1.05;
 lymin = lymin - max(abs(dymin),abs(dymax))*1.05;
 lymax = lymax + max(abs(dymin),abs(dymax))*1.05;
 % Make plot
 figure
 set(gcf,'color',[ 1  1 1]);
 times = 0:timeint:timetot;
 for i = 1:length(times)
     tfact = sin(freq*times(i));
     clf;
     patch("Faces", IX(:,1:4), "Vertices", X(:,1:2)+tfact*[D(1:2:end)'; D(2:2:end)';]', "FaceColor", [1 1 0], "EdgeColor", "k");
     axis([lxmin lxmax lymin lymax])
     axis off
     title( 'Mode')
     pause(0.01)
 end
