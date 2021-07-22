% by courtesy of Prof. Mihai Cucuringu
function [ x , glob_shift ] = optimalShift(x, f_orig)

x = reshape(x, 1, length(x));
f_orig = reshape(f_orig, 1, length(f_orig));

difs =  x - f_orig;

[s, intervals_cuts] = hist(difs, 20);
    %  intervals_cuts : also returns the position of the bin centers
[mval, m_ind ] = max(s);
glob_shift = intervals_cuts(m_ind);
x = x - glob_shift ;

end