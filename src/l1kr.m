function x_kr = l1kr(A, y)
%  x = l1kr(A, y)
% 
%  INPUTS 
%  @A: measurement matrix  
%  @y: obeserved vector 
%
%  OUTPUTS
%  @x: solution to the compressed system
%
%  AUTHORS
%    Belhassan Bayer, Nidhal Bouynaya, Gregory Ditzler and Roman 
%
%  MAINTAINER
%    Gregory Ditzler (gregory.ditzler@gmail.com)
%
%  LICENSE
%    MIT
X = null(A);
s = size(X,2); % "s=dim(ker(A))"
n = size(A,2);

cvx_begin quiet
  variable x(n,1)
  minimize(norm(x,1))
  subject to
    A*x == y; 
cvx_end

[~, i] = sort(abs(x));
j = setdiff(1:n, i(1:s));
xhat = inv(A(:, j))*y;
x_kr = zeros(n,1);
x_kr(j) = xhat;