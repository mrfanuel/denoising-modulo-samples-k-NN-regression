%% Input: 
%% complex vector z on product manifold of unit magnitude complex numbers
%% real vector x of the node grid
%% k: number of neighbours
%% Output:
%% kNN estimate (before rounding on the product of circle manifold) 


function g_est = kNN_denoise(z,x,k)


[index_nn,D] = knnsearch(x,x,'K',k);

g_est = zeros(size(z));

for i=1:length(x)

    nb_i = index_nn(i,:); 
    h_i = (1/k)*sum(z(nb_i)); 
    g_est(i) = h_i;

end
