function res = MS_wrap_around_error(x,y)
	res = zeros(size(x));
	for i=1:length(x)
		diff = abs(x(i)-y(i));
		res(i) = min(diff, 1-diff);
	end
	res = sqrt( sum( (res) .^ 2) / length(x) );
		
end