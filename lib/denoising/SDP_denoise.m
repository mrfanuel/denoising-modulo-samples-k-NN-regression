%% Input: 
%% complex vector z on product manifold of unit magnitude complex numbers
%% L: Laplacian of input graph
%% reg_param: regularizer parameter "gamma"

function gest = SDP_denoise(z, L, reg_param,n)

T = zeros(n+1);
T = [reg_param*L  -z; -z' 0]; % T is a (n+1)x(n+1) matrix

% run SDP
cvx_begin sdp

variable W(n+1,n+1) hermitian

minimize trace(T*W)
subject to 
diag(W) == 1
W == semidefinite(n+1,n+1)
cvx_end

% Extract the first n entries of the last column of W
gest = W(1:n,n+1);
