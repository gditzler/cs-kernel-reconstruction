function [A, x, y] = cs_model(m, n, k)
% [A, x, y] = cs_model(m, n, k)
% 
%  INPUTS 
%  @m: size(A,1):   
%  @n: size(A,2): 
%  @k: the k in k-sparse  
%
%  OUTPUTS
%  @A: measurement matrix sampled from Normal(0,1)
%  @x: sparsest solution  
%  @y: observation vector 
%
%  AUTHORS
%    Belhassan Bayer, Nidhal Bouynaya, and Gregory Ditzler
%
%  MAINTAINER
%    Gregory Ditzler (gregory.ditzler@gmail.com)
%
%  LICENSE
%    MIT

A = randn(m,n); 
while m ~= rank(A)
  clear A;
  A = randn(m,n);
end
x = zeros(n,1);
p = randn(k,1);   % generate "k" non-zero elements of "x"
rp = randperm(n); % generate random permutations
x(rp(1:k)) = p;   % place the "k" non-zero elements in random positions in "x"
y = A*x;          % generate the linear model: "y=Ax"