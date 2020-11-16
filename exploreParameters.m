
% i
passo_fracao = 0.01;
n = length(passo_fracao:passo_fracao:(1-passo_fracao));

% j
passo_magnitude = 0.001;
m = length(0:passo_magnitude:2);

custo2 = NaN(m,n);

clear fracao_impulso magnitude_impulso;
tic
for j=1:m
    disp(j/m*100 + "%");
    magnitude_impulso=(j-1)*passo_magnitude;
    for i=1:n
        fracao_impulso=i*passo_fracao;
        pseudocode;
        custo2(j,i) = delta_v_total_2B;
    end
end
toc
%%
figure;
y = passo_fracao:passo_fracao:(1-passo_fracao);
x = 0:passo_magnitude:2;
[X, Y] = meshgrid(x, y);
Z = custo2(1:m,1:n)';
f = @(x) x > 10;
Z(f(Z)) = NaN;
contourf(X, Y, Z);
pcolor(X, Y, Z);
shading interp
title("Custo total em função do empuxo durante swing by por Vênus");
ylabel("Fração do tempo de swing by");
xlabel("Magnitude do empuxo [km/s]");
colorbar

%%
clear min;
min.coust = min(custo(:));
for i_=1:size(custo,1)
    for j_=1:size(custo,2)
        if abs(custo(i_,j_) - min.coust) < 1e-4
            min.i = i_;
            min.j = j_;
            break;
        end
    end
end
fracao_impulso=min.i*passo_fracao;
magnitude_impulso=(min.j-1)*passo_magnitude;
disp("Min coust of " + min.coust + " km/s for delta_v = " + magnitude_impulso + " km/s at " + fracao_impulso);