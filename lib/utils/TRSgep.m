% by courtesy of Prof. Mihai Cucuringu

function [x,lam1] = TRSgep(A,a,B,Del,estimategap)
% Solves the trust-region subproblem by a generalized eigenproblem without iterations
% Compute 2 eigenvalues of 2n*2n matrix to estimate gap
% optional input: estimategap = 1 estimates the gap between largest
% eigenvalue and rest, to distinguish hard case
% 
% minimize (x^TAx)/2+ ax
% subject to x^TBx <= Del^2
%
% A: nxn symmetric, a: nx1 vector
% B: nxn symmetric positive definite
% 
% Yuji Nakatsukasa, 2015

if nargin<5
    estimategap = 1; % use gap estimate for hardcase
end

n = size(A,1);
if issparse(B)
MM1 = [sparse(n,n) B;B sparse(n,n)];
else
MM1 = [zeros(n) B;B zeros(n)];    
end

%%! Commented this out since we only check for solution on boundary!!

% [p1,flag] = pcg(A,-a,1e-12); % possible interior solution
% if norm(A*p1+a)/norm(a)<1e-5,
% if p1'*B*p1>=Del^2, p1 = nan;
% end
% else
%     p1 = nan;
% end

%!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!

% This is the core of the code
opts.tol = 1e-10;
opts.issym = 1;
opts.maxit = 5000;
if issparse(A)
    if norm(B(2:end,1))>0
%[V,lam1] = eigs(@(x)MM0timesxswap(A,B,a,Del,x),2*n,blkdiag(B,B),1,'lr',opts); 
[V,lam1] = eigs(@(x)MM0timesxswap(A,B,a,Del,-x),2*n,-blkdiag(B,B),1,'lr',opts);
%V = [V(n+1:end);V(1:n)];
    else
[V,lam1] = eigs(@(x)MM0timesx(A,B,a,Del,x),2*n,-MM1,1,'lr',opts); 
    end
else
[V,lam1] = eigs([-B A;A -a*a'/(Del^2)],-MM1,1,'lr',opts);     
%[V,lam1] = eigs([A -a*a'/(Del^2);-B A;],1,'sr',opts);
end

    if norm(real(V(:,1))) < 1e-3 %sometimes complex
        V = imag(V);    else        V = real(V);
    end

    lam1 = real(lam1(1));
    x = V(1:length(A),1); % this is parallel to soln
    normx = sqrt(x'*(B*x));         
    x = x/normx*Del; % in the easy case, this naive normalization improves accuracy
    if x'*a>0, x = -x; end % take correct sign

    if estimategap % estimate gap
    w = randn(2*n,1); VMM1 = V'*MM1; 
    w = w-VMM1*w/sqrt(VMM1*V); % random vector B-orthogonal to V
    [w,flag] = minres(@(x)MM0timesx(A,B,a,Del,x)-lam1*MM1*x,w,1e-6,2); % shifted inverse iteration with crude MINRES solver
    Mw = MM0timesx(A,B,a,Del,w);
    M1w = MM1*w; [~,ix] = max(abs(Mw));
    lam2 = Mw(ix)/M1w(ix); gap = abs(lam1-lam2); % gap estimate, factor 10 used since estimate is rough
    tolhardcase = sqrt(eps*sqrt(norm(A,'fro')^2+norm(B,'fro')+Del^2)*sqrt(n)/gap);  % tolerance for hard-case    
    else
    tolhardcase = 1e-4; % tolerance for hard-case    
    end


if normx < tolhardcase % hard case
disp(['hard case!',num2str(normx)])
x1 = V(length(A)+1:end,1);

%[xA,lamA] = eigs(A,B,1,'SA');
[x1,alpha1] = eigs(A,B,1,'SA');
%keyboard
lam1 = -alpha1;

Pvect = x1;  %first try only k=1, often enough
%Pvect = Pvect/chol(Pvect'*B*Pvect);

[x2,flag] = minres(@(x)pcgforAtilde(A,B,lam1,Pvect,alpha1,x),-a,1e-12,3000);

%x2 = pcg(@(x)pcgforAtilde(A,B,lam1,Pvect,alpha1,x),-a,1e-12,1000);
%{
[Pvect,DD] = eigs(A,B,3,'sr');
Pvect = Pvect/chol(Pvect'*B*Pvect);
x2 = pcg(@(x)pcgforAtilde(A,B,lam1,Pvect,alpha1,x),-a,1e-8,500);    
%}
if norm((A+lam1*B)*x2+a)/norm(a) > tolhardcase, % large residual, repeat
    for ii = 3*[1:3]
    if ii>length(A), break,end
    [Pvect,DD] = eigs(A,B,ii,'sr');
    %x2 = pcg(@(x)pcgforAtilde(A,B,lam1,Pvect,alpha1,x),-a,1e-12,3000);
    [x2,flag] = minres(@(x)pcgforAtilde(A,B,lam1,Pvect,alpha1,x),-a,1e-12,3000);
    if norm((A+lam1*B)*x2+a)/norm(a) < tolhardcase, break, end
    end
end

Bx = B*x1; Bx2 = B*x2; aa = x1'*(Bx); bb = 2*x2'*Bx; cc = (x2'*Bx2-Del^2); 
alp = (-bb+sqrt(bb^2-4*aa*cc))/(2*aa); %norm(x2+alp*x)-Delta
if ((x2+alp*x1)'*A*(x2+alp*x1))/2+a'*(x2-alp*x1) < ((x2-alp*x1)'*A*(x2-alp*x1))/2+(x2-alp*x1)'*p1,
x = x2+alp*x1;
else
x = x2-alp*x1;    
end
end

end



function [y] = MM0timesx(A,B,g,Delta,x)
% MM0 = [-B A;
%         A -g*g'/Delta^2];
n = size(A,1); 
x1 = x(1:n); x2 = x(n+1:end);
y1 = -B*x1 + A*x2;
y2 = A*x1-g*(g'*x2)/Delta^2;
y = [y1;y2];
end

function [y] = MM0timesxswap(A,B,g,Delta,x)
% MM0 = [-B A;
%         A -g*g'/Delta^2];
n = size(A,1); 
x1 = x(1:n); x2 = x(n+1:end);
y1 = -B*x1 + A*x2;
y2 = A*x1-g*(g'*x2)/Delta^2;
y = -[y2;y1];
end

function [y] = pcgforAtilde(A,B,lamA,Pvect,alpha1,x)

[n,m] = size(Pvect);
y = A*x+lamA*(B*x);

for i=1:m
    y = y+(alpha1*(x'*(B*Pvect(:,i))))*(B*Pvect(:,i));
end
end
