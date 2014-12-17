function A=zero_out_smalls(A, small_cri)
%% Zero Out Smalls 
% To zero out the small elements in an input vector or matrix whose absolute values 
% are smaller than a given criterion, specified by small_cri. By default, 
% small_cri=eps. Due to operations on floating point numbers, many quantities 
% should be zero mathematically but are actually not. If these "fuzzy" quantities
% are carried over in further calculations, erroneous results may arise.
% They will also unnecessarily increase the number of non-zero elements in
% sparse matrices and slow down the calculations.
%
% Inputs: 
%   A, the vector or matrix of which the small elements are desired to be zeroed out.
%   small_cri is an optional input as a criterion for being small, its
%   default value is EPS.
% 
% Output: 
%    A, the same vector or matrix as the input but with the small elements zeroed out.     
%

% Zhigang Xu, May 2008, xuz@dfo-mpo.gc.ca, updated on July 5, 2008 for the help message. 

%% Default criterion for smalls
if nargin < 2
    small_cri=eps;
end

%% Validate the inputs
if small_cri < 0
    error('The small_cri should be a non negtive number.')
end

%% Replace the smalls with zeros for a full or a sparse matrix
% When A is non-sparse, to replace the small elements with zero is easy.
% When A is sparse, we have to be cautious for the replacement. It is necessary to avoid 
% indexing on A on the left side of the equation sign, otherwise Matlab will get very slow 
% when the size of A is very large! (cf
% http://blogs.mathworks.com/loren/2007/03/01/creating-sparse-finite-element-matrices-in-matlab/#1. )

if ~issparse(A)
    id=abs(A)<=small_cri;
    A(id)=0;
else
   B=abs(A)>small_cri;
   A=A.*B;
end

end