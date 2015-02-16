function [x_kr, x_l1] = l1_approximate_reconstruction(A, y)
%  x = l1_approximate_reconstruction(A, y)
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

x_l1 = x;
x_kr = x;
[~, i] = sort(abs(x));

j = setdiff(1:n, i(1:s));
xhat = inv(A(:, j))*y;
x_kr = zeros(n,1);
x_kr(j) = xhat;
