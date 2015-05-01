function [x_kr, x_l1] = l1kr_noise(A, y, epsilon)
%  [x_kr, x_l1] = l1kr_noise(A, y, shift)
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

X = null(A);
s = size(X, 2); % "s=dim(ker(A))"
n = size(A, 2);

% minimize the l1-norm
cvx_begin quiet
  variable x(n,1)
  minimize(norm(x, 1))
  subject to
    norm(A*x - y, 2) <= epsilon;
cvx_end

% save the solution to the l1 problem before approximating the kernel
% reconstruction. 
x_tmp = x;
x_l1 = x;
sgn = sign(x);

[~, i] = sort(abs(x));

smallest = i(1:s);     % find the s+delta smallest entriries 
largest = setdiff(1:n, smallest);
x_tmp = x_tmp(largest);
lgst_smallest = max(abs(x(smallest)));
sgn_largest = sign(x_tmp);

x_tmp(sgn_largest == 1) =  x_tmp(sgn_largest == 1) - lgst_smallest;
x_tmp(sgn_largest == -1) = x_tmp(sgn_largest == -1) + lgst_smallest;
x_tmp2 = zeros(n, 1);
x_tmp2(largest) = x_tmp;

sup_set = find(abs(x_tmp2) <= epsilon);
smallest = union(smallest, sup_set);

combrows = combnk(smallest, s);  % generate combinations of the s+delta indices

% loop over the possibilites of the s+delta entries that could be tested
% for being a `zero` entry. 
sp = zeros(size(combrows, 1), 1);
for r = 1:size(combrows, 1)
  j = setdiff(1:n, combrows(r, :));

  xhat = A(:, j)\y;
  x_kr = zeros(n, 1);
  x_kr(j) = xhat;
  [~, q] = sort(abs(x_kr));
  x_kr_final = x_kr;
  
  for k = 1:n
    if norm(x_kr(q(1:k)), 2) < epsilon
      x_kr_final(q(k)) = 0;
    else
      break;
    end
  end
  
  sp(r) = sum(abs(x_kr_final) > 10e-15);  % sparisty 
end
[~, i] = sort(sp);


% solve for the sparest solution again
j = setdiff(1:n, combrows(i(1), :));
xhat = A(:, j)\y;
x_kr = zeros(n, 1);
x_kr(j) = xhat;
q = sort(abs(x_kr));
x_kr_final = x_kr;
for k = 1:n
  if norm(x_kr(q(1:k)), 2) < epsilon
    x_kr_final(q(k)) = 0;
  else
    break;
  end
end
x_kr = x_kr_final;
