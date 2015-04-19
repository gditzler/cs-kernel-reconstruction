function [x_kr, x_l1] = l1kr(A, y, shift)
%  [x_kr, x_l1] = l1_approximate_reconstruction(A, y)
% 
%  INPUTS 
%  @A: measurement matrix  
%  @y: obeserved vector 
%  @shift: 
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
  shift = 3;
end
X = null(A);
s = size(X,2); % "s=dim(ker(A))"
n = size(A,2);

% minimize the l1-norm
cvx_begin quiet
  variable x(n,1)
  minimize(norm(x,1))
  subject to
    A*x == y; 
cvx_end

% save the solution to the l1 problem before approximating the kernel
% reconstruction. 
x_l1 = x;

[~, i] = sort(abs(x));

smallest = i(1:(s+shift));     % find the s+delta smallest entriries 
combrows = combnk(smallest, s);  % generate combinations of the s+delta indices

% loop over the possibilites of the s+delta entries that could be tested
% for being a `zero` entry. 
sp = zeros(size(combrows, 1), 1);
parfor r = 1:size(combrows, 1)
  j = setdiff(1:n, combrows(r, :));

  xhat = A(:, j)\y;
  x_kr = zeros(n, 1);
  x_kr(j) = xhat;
  sp(r) = sum(abs(x_kr) > 10e-15);  % sparisty 
end
[~, i] = sort(sp);

% solve for the sparest solution again
j = setdiff(1:n, combrows(i(1), :));
xhat = A(:, j)\y;
x_kr = zeros(n,1);
x_kr(j) = xhat;