% clc;

% disp("Loading [A] and {b}")
% tic
% S = tripletToSparse('./big.matrix');
% force = readVector('./big.force');
% toc

disp("Solving...")
tic
R = chol(S);
toc