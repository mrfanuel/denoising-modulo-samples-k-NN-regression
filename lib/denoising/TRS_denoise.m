% by courtesy of Prof. Mihai Cucuringu

function [ gest_trs ] = TRS_denoise(z,L,reg_param,n)

zbar = [real(z); imag(z)]; % 2n x 1 vector
Hbar = [reg_param*L zeros(n,n); zeros(n,n) reg_param*L];  % 2n x 2n matrix

%--------------------------------------------------------------------------
% Call solver to return solution to QCQP:
%
% min gbar^T Hbar gbar - gbar^T zbar
% s.t gbar^T gbar = n
% Note: gbar is a real valued vector of length 2n; gbar = [real(g) imag(g)]
% where g is a complex valued vector of length n,
%--------------------------------------------------------------------------

[gbar_sol,lam1] = TRSgep(2*Hbar,(-2)*zbar,speye(2*n),sqrt(n),1);

gest_trs  = gbar_sol(1:n,1) + 1i* gbar_sol(n+1:2*n,1); % solution of TRS as a complex vector of length n

% g = temp./abs(temp); % Normalize each coordintate to have unit magnitude
% [theta,r] = cart2pol(real(g), imag(g)); % Cartesian to Polar coordinates
% theta1 = mod(theta,2*pi); 
% y_denoise = theta1/(2*pi); % Extract the final mod 1 denoised value
% 
% y_denoise = y_denoise';  % so that all the returned values are of size 1 x n
end
