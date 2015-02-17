function [x_kr, x_l1] = l1_approximate_reconstruction(A, y, noise, delta)
%  [x_kr, x_l1] = l1_approximate_reconstruction(A, y)
% 
%  INPUTS 
%  @A: measurement matrix  
%  @y: obeserved vector 
%  @noise: use noise implementation [boolean] [optional]
%  @delta: 
%
%  OUTPUTS
%  @x_kr: solution to the system using the l1-kr approx
%  @x_l1: solution to the system using l1
%
%  AUTHORS
%    Belhassan Bayer, Nidhal Bouynaya, Gregory Ditzler and Roman 
%
%  MAINTAINER
%    Gregory Ditzler (gregory.ditzler@gmail.com)
%
%  LICENSE
%    MIT
if nargin == 2 
  noise = false;
  delta = 0;
end
if nargin == 3
  delta = .001/sqrt(size(A,2));
end
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

if noise == false
  [~, i] = sort(abs(x));
  j = setdiff(1:n, i(1:s));
  xhat = inv(A(:, j))*y;
  x_kr = zeros(n,1);
  x_kr(j) = xhat;
else
  x_tmp = x;
  x_tmp(abs(x_tmp) <= delta/sqrt(n)) = 0;
  
  [~, i] = sort(abs(x));
  j = setdiff(1:n, i(1:s));
  xhat = inv(A(:, j))*y;
  x_kr = zeros(n,1);
  x_kr(j) = xhat;
end