

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

