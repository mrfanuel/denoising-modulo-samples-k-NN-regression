% by courtesy of Prof. Mihai Cucuringu
%% UCQP: returns the solution of the unconstrained "least squares with smooth regularization" problem
%%
%%
%% Input: 
%% complex vector z on product manifold of unit magnitude complex numbers
%% L: Laplacian of input graph
%% reg_param: regularizer parameter "gamma"

function gest = UCQP_denoise(z, L, reg_param,n)

A = eye(n) + reg_param*L;
gest = A\z;