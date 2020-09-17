
%--------------------------------------------------------------------------
%----------------------NOTES-----------------------------------------------
%--------------------------------------------------------------------------

%Dependencies:
%This program calls 'nBodyWpar.m' in the integration.
%It calls 'rstTOijk.m' when transforming coordinates of the initial
%conditions.

%Use:
%Just type 'NbodyOne' from the matlab command line with the program  and
%all it's dependencies in the search path.

%The values of the masses and the initial conditions can easily be tinkered
%with if you want to experiment with the system.  More drastic changes
%would be needed to adjust the number of bodies (because of the way the
%initial conditions are entered) but this program could be used as a modle
%for others with other numbers of bodies.

%Purpose:
%N body simulation of 5 bodies with symmetric initial data.
%This program and it's results are described in the introduction set of
%notes for the N-Body problem.  The special component of this code is the
%way that the intial conditions are specified (essentially in spherical
%coordinates.

%Output:
%The output is a plot of the trajectories in R^3.




N=5;                         %number of bodies
G=1;                       %Gravational Constant
Mass=[1.0 1.0 1.0 1.0 1.0];  %their masses

%NOTE:  The masses are currently all equal.  I the problem described in the
%notes, the output of this program with the present masses is compaired to
%the output when the value of the first mass is changed to 2.

%Coordinates given in special system have to be converted to Rectangular.

PI=3.1415926535898;
c=PI/180;

%Initial positions in r alpha and beta
r01=1.0;
r02=1.0;
r03=1.0;
r04=1.0;
r05=1.0;
ar1=0.0*c;
ar2=72.0*c;
ar3=144.0*c;
ar4=216.0*c;
ar5=288.0*c;
br1=0.0;
br2=0.0;
br3=0.0;
br4=0.0;
br5=0.0;

%initial velocities in v alpha and beta.

v01=0.6;
v02=0.6;
v03=0.6;
v04=0.6;
v05=0.6;
av1=90.0*c;
av2=90.0*c;
av3=90.0*c;
av4=90.0*c;
av5=90.0*c;
bv1=0.0;
bv2=0.0;
bv3=0.0;
bv4=0.0;
bv5=0.0;

Rec=0;
%convert positions to cartesean frame

Rec(1:3,1)= [r01*cos(ar1)*cos(br1);r01*sin(ar1)*cos(br1);r01*sin(br1)];
Rec(4:6,1)= [r02*cos(ar2)*cos(br2);r02*sin(ar2)*cos(br2);r01*sin(br2)];
Rec(7:9,1)= [r03*cos(ar3)*cos(br3);r03*sin(ar3)*cos(br3);r03*sin(br3)];
Rec(10:12,1)= [r04*cos(ar4)*cos(br4);r04*sin(ar4)*cos(br4);r04*sin(br4)];
Rec(13:15,1)= [r05*cos(ar5)*cos(br5);r05*sin(ar5)*cos(br5);r05*sin(br5)];


A=0;
%convert spherical velocities to rst frame
A(1:3,1)=[v01*cos(av1)*cos(bv1);v01*sin(av1)*cos(bv1);v01*sin(bv1)];
A(4:6,1)=[v02*cos(av2)*cos(bv2);v02*sin(av2)*cos(bv2);v02*sin(bv2)];
A(7:9,1)=[v03*cos(av3)*cos(bv3);v01*sin(av3)*cos(bv3);v03*sin(bv3)];
A(10:12,1)=[v04*cos(av4)*cos(bv4);v04*sin(av4)*cos(bv4);v04*sin(bv4)];
A(13:15,1)=[v05*cos(av5)*cos(bv5);v05*sin(av5)*cos(bv5);v05*sin(bv5)];

%convert rst velocities to ijk.  Each velociety has a different matrix.

Rec(16:18,1)= rstTOijk(Rec(1:3,1))*A(1:3,1);
Rec(19:21,1)= rstTOijk(Rec(4:6,1))*A(4:6,1);
Rec(22:24,1)= rstTOijk(Rec(7:9,1))*A(7:9,1);
Rec(25:27,1)= rstTOijk(Rec(10:12,1))*A(10:12,1);
Rec(28:30,1)= rstTOijk(Rec(13:15,1))*A(13:15,1);




%Integrate the System
tspan = [0 10];
y0=Rec

options=odeset('RelTol',1e-10,'AbsTol',1e-10);
[t,y] = ode113('nBodyWpar',tspan,y0,options,flag,N,G,Mass);




%--------------Below Should not be modified----------------------------
%------------The program should be fit to the problem------------------
%------------data by adjusting the above and 'Nbody'-------------------
%-----------The only exception is the plotting-------------------------
A=size(y)                 %vector containing (rows,columns)
timeMax=A(1,1)             %number of rows of y is number of time steps

M=0;                       %initialize M
for i=1:N
    M=M+Mass(1,i);
end


%The next section computes the 10 classical integrals of motion.

%First the components of the velociety of the center
%of mass.  Then the components of the position of the
%center of mass (which are not constant, but linear in
%time.)  Then the appropriate difference which is a
%constant.

%Next the three components of angular momentum are computed.

%And finally energy.


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


%plotting

%Three dimensional plotting

%Plot all bodies in physical space
%plot3(y(:,1), y(:,2), y(:,3))
%plot3(y(:,4),y(:,5),y(:,6))
plot3(y(:,1), y(:,2), y(:,3),'r', y(:,4), y(:,5), y(:,6), 'b', y(:,7),y(:,8), y(:,9), 'g', y(:,10), y(:,11), y(:,12), 'k', y(:,13), y(:,14),y(:,15), 'y')
grid on
xlabel('x axis'), ylabel('y axis'), zlabel('z axis')

%Plotting the integrals

%plot(t(1:timeMax),I5(1:timeMax))




