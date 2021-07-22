%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input: 
%% real vector f_mod1_denoised on a 1D grid 
%% Output:
%% real vector f_unwrapped (unwrapped signal on a 1D grid)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_unwrapped = unwrap_1D(f_mod1_denoised)

f_unwrapped = zeros(size(f_mod1_denoised));

f_unwrapped(1) = f_mod1_denoised(1);

for i=2:length(f_mod1_denoised)
    d = f_mod1_denoised(i) - f_mod1_denoised(i-1);
    if abs(d)<0.5
        f_unwrapped(i) = f_unwrapped(i-1) + d;
    elseif d<-0.5
        f_unwrapped(i) = f_unwrapped(i-1) + d + 1;
    elseif d>0.5
        f_unwrapped(i) = f_unwrapped(i-1) + d - 1;
    end
end