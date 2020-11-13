
%--------------------------------------------------------------------------
%--------------------------------Plotting----------------------------------
%--------------------------------------------------------------------------

%COLORING:
%Specify that the first six colors used will be red, blue, green, black,
%yellow, and magenta.
colorMatrix(1:6, 1:3)=[
   0, 1, 0;
   1, 0, 1;
   1, 0, 0;
   1, 1, 0;
   0, 0, 0;
   0, 0, 1
];
%For maore than six bodies the plot will be a mess anyway, but fill in the
%rest of the colors with the 'hot' command
colorMatrix(7:N+2, 1:3)=hot((N+2-6));

%% Plotting observables
if exist('I5', 'var')
    plot(t(1:timeMax),I5(1:timeMax))
    grid on
    title('Moment of Inertia vs. time')
    xlabel('time axis'), ylabel ('I')
end

if exist('E', 'var')
    figure
    plot(t(1:timeMax),E(1:timeMax))
    grid on
    title('Energy vs. Time')
    xlabel('time axis'), ylabel ('Energy')
end

if exist('CM1', 'var')
    figure
    plot3(CM1(1:timeMax), CM2(1:timeMax), CM3(1:timeMax))
    grid on
    title('Path of the center of Mass')
    xlabel('x coor CM'), ylabel('y coor CM'), zlabel('z coor CM')
end

if exist('deltaE', 'var')
    figure
    plot(t(1:timeMax), deltaE(1:timeMax))
    grid on
    title('delta energy as a function of time')
    xlabel('time'), ylabel('delta energy')
end


%% PLOTS OF TRAJECTORIES
theta = [0:pi/180:2*pi];
earth.x = R_earth_sun.*cos(theta)/ud;
earth.y = R_earth_sun.*sin(theta)/ud;
earth.z = theta*0/ud;
mars.x = R_mars_sun.*cos(theta)/ud;
mars.y = R_mars_sun.*sin(theta)/ud;
mars.z = theta*0/ud;
venus.x = R_venus_sun.*cos(theta)/ud;
venus.y = R_venus_sun.*sin(theta)/ud;
venus.z = theta*0/ud;

if exist('y', 'var')
    figure
    hold on
    for k=1:N
        plot3(y(:, 3*k-2), y(:, 3*k-1), y(:, 3*k), 'Color', colorMatrix(k,:))
        plot3(y(1, 3*k-2), y(1, 3*k-1), y(1, 3*k), 'o' ,'Color', colorMatrix(k,:))
        plot3(y(timeMax, 3*k-2), y(timeMax, 3*k-1), y(timeMax, 3*k), '*','Color', colorMatrix(k,:))
    end
    grid on
    title('The trajectories of the bodies (o is start, * is end)')
    xlabel('x axis'), ylabel('y axis'), zlabel('z axis');
    plot3(earth.x, earth.y, earth.z, '--', 'Color', 'blue');
    plot3(mars.x, mars.y, mars.z, '--', 'Color', 'red');
    plot3(venus.x, venus.y, venus.z, '--', 'Color', 'black');
    axis equal;
end

%%
if exist('CM1', 'var')
    figure
    CM = [CM1' CM2' CM3'];
    hold on
    for k=1:N
        plot3(y(:, 3*k-2)-CM(:,1), y(:, 3*k-1)-CM(:,2), y(:, 3*k)-CM(:,3),'Color', colorMatrix(k,:))
        plot3(y(1, 3*k-2)-CM(1,1), y(1, 3*k-1)-CM(1,2), y(1, 3*k)-CM(1,3), 'o' ,'Color', colorMatrix(k,:))
        plot3(y(timeMax, 3*k-2)-CM(timeMax,1), y(timeMax, 3*k-1)-CM(timeMax,2), y(timeMax, 3*k)-CM(timeMax,3), '*' ,'Color', colorMatrix(k,:))
    end
    grid on
    title('The trajectories of the bodies in center of mass coordinates (o is start, * is end)')
    xlabel('x axis'), ylabel('y axis'), zlabel('z axis');
    plot3(earth.x, earth.y, earth.z, 'Color', 'blue');
    plot3(mars.x, mars.y, mars.z, 'Color', 'red');
    plot3(earth.x, earth.y, earth.z, '--', 'Color', 'blue');
    plot3(mars.x, mars.y, mars.z, '--', 'Color', 'red');
    plot3(venus.x, venus.y, venus.z, '--', 'Color', 'black');
    axis equal;
end