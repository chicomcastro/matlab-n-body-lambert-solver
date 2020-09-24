clear;
opts = optimset(...
    'TolFun', 1e-8,...
    'TolX', 1e-8,...
    'MaxFunEvals', 1e0,...
    'MaxIter', 1e2...
);
tic
itr = 0;
global shouldPlot;
shouldPlot = 0;
x = [ ...
    2*pi,  ...  % rad
    50.460 ...  % km/s
];%random_uniform(lower_boundary, upper_boundary);
while itr < 20
    itr = itr + 1;
    [y, v] = coust(x);

    disp("---");
    disp("Iteração: " + itr);
    disp("Distância a Marte: " + y + " UA");
    disp("Velocidade de saída: " + x(2) + " km/s");
    disp("Custo: " + v + " km/s");
    
    x(2) = x(2) + 1e0/1e3;
end
toc