
% i
passo_fracao = 0.005;
n = length(passo_fracao:passo_fracao:(1-passo_fracao));

% j
passo_magnitude = 0.005;
m = length(0:passo_magnitude:2);

%%
resultados_custo = NaN(m,n);
resultados_t_voo = NaN(m,n);
resultados_r_p = NaN(m,n);
resultados_phase_mars = NaN(m,n);
clear fracao_impulso magnitude_impulso min;
tic
for j=1:m
    disp(j/m*100 + "%");
    magnitude_impulso=(j-1)*passo_magnitude;
    for i=1:n
        fracao_impulso=i*passo_fracao;
        pseudocode;
        resultados_custo(j,i) = delta_v_total_2B;
        if ~isnan(delta_v_total_2B)
            resultados_t_voo(j,i) = t_voo_total;
            resultados_r_p(j,i) = r_p;
            resultados_phase_mars(j,i) = mars_final_phase;
        end
    end
end
toc
%%
y = passo_fracao:passo_fracao:(1-passo_fracao);
x = 0:passo_magnitude:2;
[X, Y] = meshgrid(x, y);
%%
figure;
Z = resultados_custo(1:m,1:n)';
%f = @(x) x > 7.6;
%Z(f(Z)) = NaN;
contourf(X, Y, Z);
pcolor(X, Y, Z);
shading interp
title("Custo total em função do empuxo durante swing by por Vênus");
ylabel("Fração do tempo de swing by");
xlabel("Magnitude do impulso [km/s]");
colorbar
%%
figure;
Z = resultados_custo(1:m,1:n)';
%f = @(x) x > 7.6;
%Z(f(Z)) = NaN;
contourf(X, Y, Z);
pcolor(X, Y, Z);
shading interp
title("Custo total em função do empuxo durante swing by por Vênus");
ylabel("Fração do tempo de swing by");
xlabel("Magnitude do impulso [km/s]");
colormap(parula(5))
colorbar
%%
figure;
Z = resultados_custo(1:m,1:n)';
f = @(x) x > 8;
Z(f(Z)) = NaN;
contourf(X, Y, Z);
pcolor(X, Y, Z);
shading interp
title("Custo total filtrado em função do empuxo durante swing by por Vênus");
ylabel("Fração do tempo de swing by");
xlabel("Magnitude do impulso [km/s]");
colorbar
%%
figure;
Z = resultados_t_voo(1:m,1:n)';
%f = @(x) x < 100;
%Z(f(Z)) = NaN;
contourf(X, Y, Z);
pcolor(X, Y, Z);
shading interp
title("Tempo de voo em função do empuxo durante swing by por Vênus");
ylabel("Fração do tempo de swing by");
xlabel("Magnitude do impulso [km/s]");
colormap(parula(5))
colorbar
%%
figure;
Z = resultados_r_p(1:m,1:n)';
%f = @(x) x < 100;
%Z(f(Z)) = NaN;
contourf(X, Y, Z);
pcolor(X, Y, Z);
shading interp
title("R_p em função do empuxo durante swing by por Vênus");
ylabel("Fração do tempo de swing by");
xlabel("Magnitude do impulso [km/s]");
colormap(parula(5))
colorbar
%%
figure;
Z = resultados_phase_mars(1:m,1:n)';
%f = @(x) x < 100;
%Z(f(Z)) = NaN;
contourf(X, Y, Z);
pcolor(X, Y, Z);
shading interp
title("Fase final de Marte em graus em função do empuxo");
ylabel("Fração do tempo de swing by");
xlabel("Magnitude do impulso [km/s]");
colormap(parula(5))
colorbar

%%
clear min;
resultados_custo = custo;
min.coust = min(resultados_custo(:));
for i_=1:size(resultados_custo,1)
    for j_=1:size(resultados_custo,2)
        if abs(resultados_custo(i_,j_) - min.coust) < 1e-4
            min.i = i_;
            min.j = j_;
            break;
        end
    end
end
fracao_impulso=min.i*passo_fracao;
magnitude_impulso=(min.j-1)*passo_magnitude;
disp("Min coust of " + min.coust + " km/s for delta_v = " + magnitude_impulso + " km/s at " + fracao_impulso);

