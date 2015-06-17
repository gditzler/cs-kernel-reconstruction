function [A, x, y] = cs_model(m, n, k, type)
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
%
%  MAINTAINER
%    Gregory Ditzler (gregory.ditzler@gmail.com)
%
%  LICENSE
%    MIT

if strcmp(type, 'Gaussian')
  A = randn(m,n); 
  while m ~= rank(A)
    A = randn(m,n);
  end
  x = zeros(n,1);
  p = randn(k,1);   % generate "k" non-zero elements of "x"
  rp = randperm(n); % generate random permutations
  x(rp(1:k)) = p;   % place the "k" non-zero elements in random positions in "x"
  y = A*x;          % generate the linear model: "y=Ax"
elseif strcmp(type, 'GaussianNoise')
  A = randn(m,n); 
  while m ~= rank(A)
    A = randn(m,n);
  end
  x = zeros(n,1);
  p = randn(k,1);   % generate "k" non-zero elements of "x"
  rp = randperm(n); % generate random permutations
  x(rp(1:k)) = p;   % place the "k" non-zero elements in random positions in "x"
  
  lvl = .005;
  e = randn(m,1);
  e = rand*lvl*e./norm(e);
  y = A*x+e; 
elseif strcmp(type, 'Uni')
  A = rand(m,n); 
  while m ~= rank(A)
    A = rand(m,n);
  end
  x = zeros(n,1);
  p = randn(k,1);   % generate "k" non-zero elements of "x"
  rp = randperm(n); % generate random permutations
  x(rp(1:k)) = p;   % place the "k" non-zero elements in random positions in "x"
  y = A*x;          % generate the linear model: "y=Ax"
elseif strcmp(type, 'Fourier')
  [Y,X] = meshgrid(0:n-1,0:n-1);
  A = 1/sqrt(n) * exp( -1i/n *X.*Y);
  
  x = zeros(n,1);
  p = randn(k,1);   % generate "k" non-zero elements of "x"
  rp = randperm(n); % generate random permutations
  x(rp(1:k)) = p;   % place the "k" non-zero elements in random positions in "x"
  y = A*x;
elseif  strcmp(type, 'Bernoulli')
  A = randn(m,n);
  A(A < 0) = 0;
  A(A ~= 0) = 1;
  x = zeros(n,1);
  p = randn(k,1);   % generate "k" non-zero elements of "x"
  rp = randperm(n); % generate random permutations
  x(rp(1:k)) = p;   % place the "k" non-zero elements in random positions in "x"
  y = A*x;          % generate the linear model: "y=Ax"
end