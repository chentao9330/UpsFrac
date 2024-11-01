function varargout = initVariablesADI(varargin)
% Initialize a set of automatic differentiation variables
%
% SYNOPSIS:
%  a            = initVariablesADI(a);
%  [a, b, c, d] = initVariablesADI(a, b, c, d);
%
% PARAMETERS:
%   Any number of variables in either column vector format or as scalars.
%   These variables will be instansiate as adi objects containing both a
%   .val field and a .jac jacobian. These variables will start with
%   identity jacobians with regards to themselves and zero jacobians with
%   regards to the other variables (implicitly defined by the ordering of
%   input and output). 
%
%   These variables can then be used to create more complex expressions,
%   resulting in automatic compuation of the first order derivatives
%   leading to easy implementation of Newton-like nonlinear solvers.
%   
% EXAMPLE:
%        x = 1;
%        y = 5;
%        [x, y] = initVariablesADI(x, y)
%
%        This gives x.jac ->  {[1]  [0]} and y.jac ->  {[0]  [1]}.
%
%        If we compute z = x.*y.^2 we get
%
%        z.val = 25 (as is expected),
%        z.jac{1} = d(x*y^2)/dx = y^2 = 5^2 = 25
%        z.jac{2} = d(x*y^2)/dy = 2*x*y = 2*1*5 = 10;
%
%        Note that as this is meant for vector operations, the
%        element-wise operations should be used (.* instead of *) even when
%        dealing with scalars.
%       
% RETURNS:
%   The same variables as inputted, as ADI objects.
%
% SEE ALSO:
%   ADI.m

%{
#COPYRIGHT#
%}

% $Date: $
% $Revision: $
n = nargin;
assert(nargin == nargout, 'Number of output variables must equal the number of input variables!');

jacList = cell(zeros(1,n));
numvals = cellfun(@numel, varargin);
varargout = cell(1,n);

for k = 1:n
    jac = jacList;
    for k1 = 1:n
        sz = [numvals(k), numvals(k1)];
        if k==k1
            jac{k1} = sparse((1:sz(1))',(1:sz(1)), 1, sz(1), sz(1));
        else
            jac{k1} = sparse(sz(1),sz(2));
        end
    end
    varargout(k) = {ADI(varargin{k}, jac)};
end
end

