function x = l1_approximate_reconstruction(A, y)
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
%    Belhassan Bayer, Nidhal Bouynaya, and Gregory Ditzler
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
x(i(1:s)) = 0;