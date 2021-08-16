%% Input: vectors x and y
%% Output: square root of mean squared error between x and y (root-mean-square error: RMSE)
function res = MS_error(x,y)

	[x_shifted , glob_shift ] = optimalShift(x, y);
	res = sqrt( sum( (x_shifted'-y) .^ 2) / length(y) );
    
end
