%% Input: x is a complex vector with each entry of unit magnitude
%% Output: extract angle (modulo 2*pi) of each entry of x, and divide it by "2*pi"

function mod_output = extract_modulo(x) 

temp = angle(x);
mod_output = mod(temp,2*pi);
mod_output = mod_output/(2*pi);