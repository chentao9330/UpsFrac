classdef ADI
% ADI class: simple implementation of automatic differentiation for easy construction of jacobian matrices.
% 
% SYNOPSIS:
%   x = ADI(value, jacobian);   
%
% PARAMETERS:
%   value - The numerical value of the object
%
%   jacobian - The jacobian of the object.
%
%   
%
% RETURNS:
%   ADI object.
%
% COMMENTS:
%  This class is typically instansiated for a set of different variables
%  using initVariablesADI. The file contains a worked example demonstrating
%  the usage for several variables.
%
% SEE ALSO:
%   initVariablesADI

%{
#COPYRIGHT#
%}

% $Date: $
% $Revision: $
   properties
      val  %function value
      jac  %list of sparse jacobian matrices
   end
   methods
      function obj = ADI(a,b)
         %ADI class constructor
         if nargin == 0 % empty constructor
            obj.val     = [];
            obj.jac     = {};
         elseif nargin == 1 % 
             if isa(a, 'ADI')
                 obj = a;
             else
                 error('Contructor requires 2 inputs')
             end
         elseif nargin == 2 % values + jacobians
             obj.val = a; % value
             if ~iscell(b)
                 b = {b};
             end
             obj.jac = b; % jacobian or list of jacobians
         else
             error('Input to constructor not valid')
         end
      end

      %--------------------------------------------------------------------
      function h = numval(u)
          h = numel(u.val);
      end
      
      %--------------------------------------------------------------------
      function h = double(u)
          h = u.val;
      end
      
      %--------------------------------------------------------------------

      function h = ge(u, v)
          h = ge(double(u), double(v));
      end
      
      %--------------------------------------------------------------------

      function h = uplus(u) 
          h = u;
      end
      
      %--------------------------------------------------------------------
      
      function h = uminus(u)
         h = ADI(-u.val, uminusJac(u.jac));
      end
      
      %--------------------------------------------------------------------
      
      function h = plus(u,v)
         if ~isa(u,'ADI')       %u is a vector
            h = ADI(u+v.val, v.jac);
         elseif ~isa(v,'ADI')   %v is a vector
            h = ADI(u.val+v, u.jac);
         else
            h = ADI(u.val+v.val, plusJac(u.jac, v.jac) ); 
         end
      end
      
      %--------------------------------------------------------------------
      
      function h = minus(u,v) 
         h = plus(u, uminus(v)); 
      end
      
      %--------------------------------------------------------------------
      
      function h = mtimes(u,v)% '*'
          if ~isa(u,'ADI') %u is a scalar/matrix
              h = ADI(u*v.val, mtimesJac(u, v.jac));
          elseif and(~isa(v,'ADI'), numel(v)==1) %v is a scalar
              h = mtimes(v,u);
          else
              error('Operation not supported')
          end
      end
      
      %--------------------------------------------------------------------
      
      function h = times(u,v)% '.*'
         if ~isa(u,'ADI') %u is a scalar/vector
             if numel(u)==1 % u is scalar
                 h = mtimes(u,v);
             else %u is vector
                 h = ADI(u.*v.val, lMultDiag(u, v.jac));
             end
         elseif ~isa(v,'ADI') %v is a scalar/vector
             h = times(v,u);
         else
             h = ADI(u.val.*v.val, timesJac(u.val, v.val, u.jac, v.jac));
         end
      end
      
      %--------------------------------------------------------------------
      
      function h = mrdivide(u,v)% '/'
         if ~isa(v,'ADI') %v is a scalar
            h = mtimes(u, 1/v);
         else
            error('Operation not supported');
         end
      end
      
      %--------------------------------------------------------------------
      
      function h = mldivide(u,v)% '\'
          if ~isa(u,'ADI') %u is a scalar/matrix
              h = ADI(u\v.val, mldivideJac(u, v.jac));
          else
              error('Operation not supported');
          end
      end
      
      %--------------------------------------------------------------------
      
      function h = power(u,v)% '.^'
          if ~isa(v,'ADI') %v is a scalar
              h = ADI(u.val.^v, lMultDiag(v.*u.val.^(v-1), u.jac));
          else
              error('Operation not supported')
          end
      end
      
      %--------------------------------------------------------------------
      
      function h = rdivide(u,v)% './'
          h = times(u, power(v, -1));
      end
      
      %--------------------------------------------------------------------
      
      function h = ldivide(u,v)% '.\'
          h = rdivide(v,u);
      end
      
      %--------------------------------------------------------------------
      
      function h = subsref(u,s)
          switch s(1).type
              case '.'
                  h = builtin('subsref',u,s);
              case '()'
                  subs  = s.subs{:};
                  if subs == ':' 
                      h = u; 
                  else
                      if islogical(subs), subs = find(subs); end
                      h = ADI(u.val(subs), subsrefJac(u.jac, subs)); 
                  end
              case '{}'     
                  error('Operation not supported');
          end
      end
      
      %--------------------------------------------------------------------
      
      function u = subsasgn(u,s,v)
          switch s(1).type
              case '.'
                  u = builtin('subsasgn',u,s,v);
              case '()'
                  subs  = s.subs{:};
                  if ~isa(u, 'ADI') % u is a vector
                      u = double2AD(u, v.jac);
                  end
                  u.val(subs) = v;
                  if ~isa(v, 'ADI') % v is a constant vector
                      u.jac = subsasgnJac(u.jac, subs); % set rows to zero
                  else
                      u.jac = subsasgnJac(u.jac, subs, v.jac);
                  end
              case '{}'
                  error('Operation not supported');
          end
      end
      
      %--------------------------------------------------------------------
      
      function h = exp(u)
          eu = exp(u.val);
          h  = ADI(eu, lMultDiag(eu, u.jac));
      end
      
      %--------------------------------------------------------------------
      
      function h = max(u,v) % this function should be expanded
          if ~isa(u,'ADI') %u is a vector
              [value, inx] = max([u v.val], [], 2);
              h  = ADI(value, lMultDiag(inx==2, v.jac));
          elseif ~isa(v,'ADI') %v is a vector
              h = max(v,u);
          else
              error('Not yet implemented ...');
          end
      end
      
      %--------------------------------------------------------------------
      
      function h = vertcat(varargin)
          nv    = numel(varargin);
          nj    = numel(varargin{1}.jac);
          vals  = cell(1,nv);
          jacs  = cell(1,nv);
          sjacs = cell(1, nj);
          for k = 1:nv
              vals{k} = varargin{k}.val;
          end
          for k = 1:nj
              for k1 = 1:nv
                  sjacs{k1} = varargin{k1}.jac{k};
              end
              jacs{k} = vertcatJac(sjacs{:});
          end
          h = ADI(vertcat(vals{:}), jacs);
      end
      
      
      %--------------------------------------------------------------------
      
      function h = cat(varargin)
          h = vertcat(varargin{:});
          h = ADI(h.val, horzcatJac(h.jac{:}));
      end
      
      %--------------------------------------------------------------------
      
      function horzcat(varargin)
          error('horzcat doesn not make sense for class ADI')
      end
      
      %--------------------------------------------------------------------
      
      function h = interpReg(T, u, reginx)
          [y, dydu] = interpReg(T, u.val, reginx);
          h = ADI(y, lMultDiag(dydu, u.jac));
      end
      
      %--------------------------------------------------------------------
      
      function h = interpRegPVT(T, x, v, flag, reginx)
          [y, dydx, dydv] = interpRegPVT(T, x.val, v.val, flag, reginx);
          h = ADI(y, timesJac(dydx, dydv, v.jac, x.jac)); %note order of input
      end
      
      %--------------------------------------------------------------------
      
%       function u = addToVals(u, inx, v)
%           % adds v to u(inx)
%           assert(numel(inx)==numel(v.val));
%           u.val(inx) = u.val(inx) + v.val;
%           for k = 1:numel(u.jac)
%               u.jac{k} = addToRows(u.jac{k}, inx, v.jac{k});
%           end
%       end
      %--------------------------------------------------------------------
   end
end

%**************************************************************************
%-------- Helper functions involving Jacobians  ---------------------------
%**************************************************************************

function J = uminusJac(J1)
J = cellfun(@uminus, J1, 'UniformOutput', false);
end

%--------------------------------------------------------------------------

function J = plusJac(J1, J2)
J = cellfun(@plus, J1, J2, 'UniformOutput', false);
end

%--------------------------------------------------------------------------

function J = mtimesJac(M, J1)
J = cell(1, numel(J1));
for k = 1:numel(J)
    J{k} = M*J1{k};
end
end

%--------------------------------------------------------------------------

function J = lMultDiag(d, J1)
n = numel(d);
D = sparse((1:n)', (1:n)', d, n, n);
J = cell(1, numel(J1));
for k = 1:numel(J)
    J{k} = D*J1{k};
end
end

%--------------------------------------------------------------------------

function J = timesJac(v1, v2, J1, J2)
n  = numel(v1);
D1 = sparse((1:n)', (1:n)', v1, n, n);
D2 = sparse((1:n)', (1:n)', v2, n, n);
J = cell(1, numel(J1));
for k = 1:numel(J)
    J{k} = D1*J2{k} + D2*J1{k};
end
end

%--------------------------------------------------------------------------

function J = mldivideJac(M, J1)
J = cell(1, numel(J1));
for k = 1:numel(J)
    J{k} = M\J1{k};
end
end

%--------------------------------------------------------------------------

function J = subsrefJac(J1, subs)
J = cell(1, numel(J1));
for k = 1:numel(J)
    J{k} = J1{k}(subs,:);
end
end

%--------------------------------------------------------------------------

function u = double2AD(u, J1)
% u is vector, J reference jacobian
nr = numel(u);
J  = cell(1, numel(J1));
for k = 1:numel(J)
    nc   = size(J1{k}, 2);
    J{k} = sparse(nr, nc);
end
u = ADI(u, J);
end

%--------------------------------------------------------------------------

function J = subsasgnJac(J, subs, J1)
if nargin == 3
    for k = 1:numel(J)
        J{k}(subs,:) = J1{k};
    end
else
    for k = 1:numel(J)
        J{k}(subs,:) = 0;
    end
end
end

%--------------------------------------------------------------------------

function J = vertcatJac(varargin)
J = vertcat(varargin{:});
end

%--------------------------------------------------------------------------

function J = horzcatJac(varargin)
J = horzcat(varargin{:});
end         
         
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
