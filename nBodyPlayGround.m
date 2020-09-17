%N body simulation in inerital frame.

clear   %clears workspace to get ready for run

%--------------------------------------------------------------------------
%----------------------NOTES-----------------------------------------------
%--------------------------------------------------------------------------

%USE:
%User specifies the number of bodies 'N', the 'simulationTime', and the
%gravitationa constant 'G'.

%The program generates random masses and initial conditions for the bodies,
%integrates them and plots the trajectories and some analysis.

%OUTPUT:
%The output is several plots.

%Two three dimensional plots are given. One is a plot of all the
%trajectories in the given reference frame, and the other in coordinates
%with the center of mass at the orgin.  Both are plots in R^3 and both sets
%of coordinates are inertial.  Each trajectory has it's own color although
%for large N these will get harder to tell apart.  Each trajectory is
%marked with a circle at time zero and an astric at the final time to mark
%the beginning and end of the trajectory.

%Several two dimensional plots are output as well.  The most important is
%the plot of the change in energy as a function of time.  This is a good
%indicator of how much integration errof is creeping into the problem.  The
%energy will usually jump when the system passes near a collision, but as
%long as the scale of the axis is less than the acceptable error tolerence
%you have in mind, then the results are likely good (at least for short
%times).  For example if you requair nine signigfigant figures, and the
%drift in the energy is on the order of 10^(-11), this is a good indicator
%that your requirements have been met.

%PURPOSE:
%The program is fun to play with, but more importantly it tests the file
%'nBodyWpar' (which is the N-Body vector field with parameters for G, N
%and masses).  The present program demonstrates that 'nBodyWpar' properly
%adapts to any number of masses, with any values of the mass and
%gravitational constant.  It also demonstrates how to use it in a program.

%The code to compute the constants of motion could be put into their own
%m-files.  This would greatly reduce the number of lined of code in theis
%program but I wanted to keep the number of files down.


%--------------------------------------------------------------------------
%----------------------USER SPECIFIED--------------------------------------
%------------------------PARAMETERS----------------------------------------
%--------------------------------------------------------------------------

%As the user you specify these.  The Program does the rest.
%Remember that the the system of ODEs is a 6*N system of first order
%scalar ODEs and expect matlab to bog down accordingly.


%If G=1 then integration show interesting behavior for 'simulationTime'
%between 2 and 10.  For such values the program run in a resonable amount
%of time for 'N' between 2 and 12.  (Run times often one minute or much
%less, however if the randomly chosen initial conditions pass the system
%near a singularity then the program will run slower even in this range).

%I've run the program with as many as 40 bodies on a home desk top and
%obtained results in an hour.  I've had out of memory errore for more than
%45 bodies (with simulation times over 10).

N=12                 %number of bodies
simulationTime=4    %How long to run the simulation
G=1                  %gravitational constant

%The value of the gravitation constant in the metric system is
%
%G_metric=6.67300 × 10-11 (m3 kg-1 s-2)
%
%which is very small

%--------------------------------------------------------------------------
%----------PROGRAM:  mess with the following at your own risk...-----------
%--------------------------------------------------------------------------


%----------------------DATA------------------------------------------------

%Random Masses
Mass =(1/G)*rand(1, N);       %must be a row vector
                              %rand gives uniformly distributed random
                              %values for the masses between 0 and 1.
                              %if graviational constant is small (like
                              %in the metric system) then multiplying by
                              %1/G makes the masses ''big'' enough to
                              %interact strongly.

%Random Initial conditions
randIc=((1)/(2*G))*randn(1, 6*N); %must be a row vector
                                  %'randn' gives a Gaussian distribution of
                                  %random values centered around zero.
                                  %Then positions and velocities may be
                                  %positive or negative.
                                  %The 1/2G gives ''big'' positions
                                  %and velocities if G is small.  There are
                                  %6*N as each body has 3 position
                                  %components and 3 velocity components and
                                  %there are N such bodies

%------------------INTEGRATION---------------------------------------------


%Integrate the System
tspan = [0 simulationTime];
y0=randIc;
options=odeset('RelTol',1e-12,'AbsTol',1e-12);
[t,y] = ode113('nBodyWpar',tspan,y0,options,flag,N,G,Mass);


integrationDone=1

%--------Quantities neede for the constants of motion loops----------------


A=size(y)                 %vector containing (rows,columns)
timeMax=A(1,1)             %number of rows of y is number of time steps

M=0;                       %initialize M
for i=1:N
    M=M+Mass(1,i);   %sum of the masses is needed to
                     %compute the constantsclo of motion
end

%--------------------------------------------------------------------------
%-----------------Analysis of Integrals of Motion--------------------------
%--------------------------------------------------------------------------


%The next section computes the 10 classical integrals of motion.

%First the components of the velociety of the center
%of mass.  Then the components of the position of the
%center of mass (which are not constant, but linear in
%time.)  Then the appropriate difference which is a
%constant.

%Next the three components of angular momentum are computed.

%And finally energy, moment of inertia, and it's derivatives.

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%Linear Momentum


%first integral; component one of V_cm
I1=0;
for i=1:timeMax
    S=0;              %Initialize S
        for j=1:N
            S=S+Mass(1,j)*y(i,3*N+(3*j-2));
        end
    I1(i)=S/M;
end

%second integral; component two of V_cm
I2=0;
for i=1:timeMax
    S=0;              %Initialize S
        for j=1:N
            S=S+Mass(1,j)*y(i,3*N+(3*j-1));
        end
    I2(i)=S/M;
end

%third integral; component three of V_cm
I3=0;
for i=1:timeMax
    S=0;              %Initialize S
        for j=1:N
            S=S+Mass(1,j)*y(i,3*N+(3*j));
        end
    I3(i)=S/M;
end

%Center of Mass as a function of time; component one of R_cm(t)
CM1=0;
for i=1:timeMax
    S=0;              %Initialize S
        for j=1:N
            S=S+Mass(1,j)*y(i,3*j-2);
        end
    CM1(i)=S/M;
end

%Center of Mass as a function of time; component two of R_cm(t)
CM2=0;
for i=1:timeMax
    S=0;              %Initialize S
        for j=1:N
            S=S+Mass(1,j)*y(i,3*j-1);
        end
    CM2(i)=S/M;
end

%Center of Mass as a function of time; component three of R_cm(t)
CM3=0;
for i=1:timeMax
    S=0;              %Initialize S
        for j=1:N
            S=S+Mass(1,j)*y(i,3*j);
        end
    CM3(i)=S/M;
end

%Fourth integral: component one of R_cm(t)-[V_cm(t)]*t
I4=0;
for i=1:timeMax
    S=0;              %Initialize S
        for j=1:N
            S=CM1(i)-I1(i)*t(i);
        end
    I4(i)=S;
end

%Fifth integral: component two of R_cm(t)-[V_cm(t)]*t
I5=0;
for i=1:timeMax
    S=0;              %Initialize S
        for j=1:N
            S=CM2(i)-I2(i)*t(i);
        end
    I5(i)=S;
end

%Sixth integral: component three of R_cm(t)-[V_cm(t)]*t
I6=0;
for i=1:timeMax
    S=0;              %Initialize S
        for j=1:N
            S=CM3(i)-I3(i)*t(i);
        end
    I6(i)=S;
end


%Angular Momentum

%seventh integral; component one of h
I7=0;
for i=1:timeMax
    S=0;              %Initialize S
        for j=1:N
            S=S+Mass(1,j)*[y(i,3*j-1)*y(i,3*N+3*j)-y(i,3*j)*y(i,3*N+3*j-1)];
        end
    I7(i)=S/M;
end

%eighth integral; component two of h
I8=0;
for i=1:timeMax
    S=0;              %Initialize S
        for j=1:N
            S=S+Mass(1,j)*[y(i,3*j)*y(i,3*N+3*j-2)-y(i,3*j-2)*y(i,3*N+3*j)];
        end
    I8(i)=S/M;
end

%nineth integral; component three of h
I9=0;
for i=1:timeMax
    S=0;              %Initialize S
        for j=1:N
            
S=S+Mass(1,j)*[y(i,3*j-2)*y(i,3*N+3*j-1)-y(i,3*j-1)*y(i,3*N+3*j-2)];
        end
    I9(i)=S/M;
end

%Energy: (T-U)(t)

%Potential Energy:  U(t)
U=0;
for i=1:timeMax
    S=0;              %Initialize S
        for j=1:N
            s=0;
            for k=1:N
                
Rjk=(y(i,3*j-2)-y(i,3*k-2))^2+(y(i,3*j-1)-y(i,3*k-1))^2+(y(i,3*j)-y(i,3*k))^2;
                if k~=j
                    s=s+Mass(1,k)/[sqrt(Rjk)];
                else
                   s=s+0;
                end
            end
            S=S+Mass(1,j)*s;
        end
    U(i)=(G/2)*S;
end


%Kenetic Energy;  T(t)
T=0;
for i=1:timeMax
    S=0;              %Initialize S
    for j=1:N
        
S=S+Mass(1,j)*[(y(i,3*N+(3*j-2)))^2+(y(i,3*N+(3*j-1)))^2+(y(i,3*N+(3*j)))^2];
    end
    T(i)=S/2;
end

%Tenth integral; Energy   E(t)
E=0;
for i=1:timeMax
    E(i)=T(i)-U(i);
end


%Aux Quanties


%I or moment of inerita
I=0;
for i=1:timeMax
    S=0;              %Initialize S
    for j=1:N
        S=S+Mass(1,j)*[(y(i,3*j-2))^2+(y(i,3*j-1))^2+(y(i,3*j))^2];
    end
    I(i)=S;
end

%d/dt(I)
Idot=0;
for i=1:timeMax
    S=0;              %Initialize S
    for j=1:N
        
S=S+Mass(1,j)*[y(i,3*N+(3*j-2))*y(i,3*j-2)+y(i,3*N+(3*j-1))*y(i,3*j-1)+y(i,3*N+(3*j))*y(i,3*j)];
    end
    Idot(i)=S*2;
end

%d^2/(dt)^2(I)
Idotdot=0;
for i=1:timeMax
    Idotdot(i)=2*T(i)+2*E(i);
end

%Change in Energy.  If integration is perfect then this is always zero.
%Its deviation from zero measures the solutions deviation from the real
%one.

%deltaE
deltaE=0;
deltaE(1)=0;
for i=2:timeMax
    deltaE(i)=E(i)-E(1);
end


%--------------------------------------------------------------------------
%--------------------------------Plotting----------------------------------
%--------------------------------------------------------------------------

%COLORING:
%Specify that the first six colors used will be red, blue, green, black,
%yellow, and magenta.
colorMatrix(1:6, 1:3)=[1, 0, 0;
                       0, 0, 1;
                       0, 1, 0;
                       0, 0, 0;
                       1, 1, 0;
                       1, 0, 1];
%For maore than six bodies the plot will be a mess anyway, but fill in the
%rest of the colors with the 'hot' command
colorMatrix(7:N+2, 1:3)=hot((N+2-6));


%Plotting observables

plot(t(1:timeMax),I5(1:timeMax))
grid on
title('Moment of Inertia vs. time')
xlabel('time axis'), ylabel ('I')

figure
plot(t(1:timeMax),E(1:timeMax))
grid on
title('Energy vs. Time')
xlabel('time axis'), ylabel ('Energy')



figure
plot3(CM1(1:timeMax), CM2(1:timeMax), CM3(1:timeMax))
grid on
title('Path of the center of Mass')
xlabel('x coor CM'), ylabel('y coor CM'), zlabel('z coor CM')

figure
plot(t(1:timeMax), deltaE(1:timeMax))
grid on
title('delta energy as a function of time')
xlabel('time'), ylabel('delta energy')



%PLOTS OF TRAJECTORIES
figure
hold on
    for k=1:N
        plot3(y(:, 3*k-2), y(:, 3*k-1), y(:, 3*k), 'Color', colorMatrix(k,:))
        plot3(y(1, 3*k-2), y(1, 3*k-1), y(1, 3*k), 'o' ,'Color', colorMatrix(k,:))
        plot3(y(timeMax, 3*k-2), y(timeMax, 3*k-1), y(timeMax, 3*k), '*','Color', colorMatrix(k,:))
    end

grid on
title('The trajectories of the bodies (o is start, * is end)')
xlabel('x axis'), ylabel('y axis'), zlabel('z axis')



figure
CM= [CM1' CM2' CM3'];
hold on
 for k=1:N
  plot3(y(:, 3*k-2)-CM(:,1), y(:, 3*k-1)-CM(:,2), y(:, 3*k)-CM(:,3),'Color', colorMatrix(k,:))
  plot3(y(1, 3*k-2)-CM(1,1), y(1, 3*k-1)-CM(1,2), y(1, 3*k)-CM(1,3), 'o' ,'Color', colorMatrix(k,:))
  plot3(y(timeMax, 3*k-2)-CM(timeMax,1), y(timeMax, 3*k-1)-CM(timeMax,2), y(timeMax, 3*k)-CM(timeMax,3), '*' ,'Color', colorMatrix(k,:))
 end
grid on
title('The trajectories of the bodies in center of mass coordinates (o is start, * is end)')
xlabel('x axis'), ylabel('y axis'), zlabel('z axis')




