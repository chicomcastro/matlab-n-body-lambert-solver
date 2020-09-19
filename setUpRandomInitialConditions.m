
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

y0 = randIc;