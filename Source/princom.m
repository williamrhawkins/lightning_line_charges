function [meanX, dirVect, t, endpts, len] = princom(X)

%performs PCA curve fitting on a 3x3 x, y, z matrix.
%   Detailed explanation goes here

meanX = mean(X, 1);
[coeff, score, ~] = pca(X);
dirVect = coeff(:,1);
t = [min(score(:,1)), max(score(:,1))];                     
endpts = [meanX + t(1)*dirVect'; meanX + t(2)*dirVect'];
len = ((endpts(1,1) - endpts(2,1)).^2 + (endpts(1,2) - endpts(2,2)).^2 + (endpts(1,3) - endpts(2,3)).^2 ).^(1/2);

end

