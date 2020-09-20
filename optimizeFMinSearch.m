
opts = optimset(...
    'TolFun', 1e-14,...
    'TolX', 1e-14,...
    'MaxFunEvals', 1e4,...
    'MaxIter', 1e4...
);

itr = 0;
while itr < 1
    itr = itr + 1;
    x0 = [2*pi 4.9775e+04];%random_uniform(lower_boundary, upper_boundary);
    x = fminsearch(@coust, x0, opts);

    disp("Iteração: " + itr);
    disp(x);
    disp("Custo: " + custo(x));
end