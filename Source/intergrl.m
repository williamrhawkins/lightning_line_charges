function [E_z] = intergrl(k, P, data)
%INTEGRAL Summary of this function goes here
%   Detailed explanation goes here

meanX = data{k,1};
dirVect = data{k,2};
t = data{k,3};
length = data{k,5};

fun = @(x) ((P(3) - dirVect(3).*x - meanX(3)))./(sqrt((P(1)-(dirVect(1).*x + meanX(1))).^2 + (P(2)-(dirVect(2).*x + meanX(2))).^2 + (P(3)-(dirVect(3).*x + meanX(3))).^2).^3);
int = integral(fun, t(1), t(2));

E_z = (9e9 ./ length) .* int;



end

