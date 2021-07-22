%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input: 
%% real matrix f_mod1_denoised on a square 2D grid 
%% Output:
%% real matrix f_unwrapped (unwrapped signal on a 2D grid)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_unwrapped = unwrap_2D(f_mod1_denoised)

%% pre-allocate the output
f_unwrapped = zeros(size(f_mod1_denoised));

%% unwrapping the first column
f_unwrapped(:,1) = unwrap_1D(f_mod1_denoised(:,1));

%% number of rows
m = size(f_mod1_denoised,1);

%% unwrapping each row
for i=1:m
    for j=2:m
        d = f_mod1_denoised(i,j) - f_mod1_denoised(i,j-1);
        if abs(d)<0.5
            f_unwrapped(i,j) = f_unwrapped(i,j-1) + d;
        elseif d<-0.5
            f_unwrapped(i,j) = f_unwrapped(i,j-1) + d + 1;
        elseif d>0.5
            f_unwrapped(i,j) = f_unwrapped(i,j-1) + d - 1;
        end
    end
end
	
	