%% Input: x is a complex vector 
%% Output: entry wise projection of x onto the product manifold of unit circles

function out_proj = project_manifold(x)

x(x==0) = 1; % entries of x which are 0 are set to 1
x = x./abs(x); % other entries are divided by their absolute value
out_proj = x; 