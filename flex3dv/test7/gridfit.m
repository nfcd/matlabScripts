% Copyright (c) 2006, John D'Errico
% All rights reserved.

% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.


function [zgrid,xgrid,ygrid] = gridfit(x,y,z,xnodes,ynodes,varargin)
% gridfit: estimates a surface on a 2d grid, based on scattered data
%          Replicates are allowed. All methods extrapolate to the grid
%          boundaries. Gridfit uses a modified ridge estimator to
%          generate the surface, where the bias is toward smoothness.
%
%          Gridfit is not an interpolant. Its goal is a smooth surface
%          that approximates your data, but allows you to control the
%          amount of smoothing.
%
% usage #1: zgrid = gridfit(x,y,z,xnodes,ynodes);
% usage #2: [zgrid,xgrid,ygrid] = gridfit(x,y,z,xnodes,ynodes);
% usage #3: zgrid = gridfit(x,y,z,xnodes,ynodes,prop,val,prop,val,...);
%
% Arguments: (input)
%  x,y,z - vectors of equal lengths, containing arbitrary scattered data
%          The only constraint on x and y is they cannot ALL fall on a
%          single line in the x-y plane. Replicate points will be treated
%          in a least squares sense.
%
%          ANY points containing a NaN are ignored in the estimation
%
%  xnodes - vector defining the nodes in the grid in the independent
%          variable (x). xnodes need not be equally spaced. xnodes
%          must completely span the data. If they do not, then the
%          'extend' property is applied, adjusting the first and last
%          nodes to be extended as necessary. See below for a complete
%          description of the 'extend' property.
%
%          If xnodes is a scalar integer, then it specifies the number
%          of equally spaced nodes between the min and max of the data.
%
%  ynodes - vector defining the nodes in the grid in the independent
%          variable (y). ynodes need not be equally spaced.
%
%          If ynodes is a scalar integer, then it specifies the number
%          of equally spaced nodes between the min and max of the data.
%
%          Also see the extend property.
%
%  Additional arguments follow in the form of property/value pairs.
%  Valid properties are:
%    'smoothness', 'interp', 'regularizer', 'solver', 'maxiter'
%    'extend'
%
%  Any UNAMBIGUOUS shortening (even down to a single letter) is
%  valid for property names. All properties have default values,
%  chosen (I hope) to give a reasonable result out of the box.
%
%   'smoothness' - scalar - determines the eventual smoothness of the
%          estimated surface. A larger value here means the surface
%          will be smoother. Smoothness must be a non-negative real
%          number.
%
%          Note: the problem is normalized in advance so that a
%          smoothness of 1 MAY generate reasonable results. If you
%          find the result is too smooth, then use a smaller value
%          for this parameter. Likewise, bumpy surfaces suggest use
%          of a larger value. (Sometimes, use of an iterative solver
%          with too small a limit on the maximum number of iterations 
%          will result in non-convergence.)
%
%          DEFAULT: 1
%
%
%   'interp' - character, denotes the interpolation scheme used
%          to interpolate the data.
%
%          DEFAULT: 'triangle'
%          
%          'bilinear' - use bilinear interpolation within the grid
%                     (also known as tensor product linear interpolation)
%
%          'triangle' - split each cell in the grid into a triangle,
%                     then linear interpolation inside each triangle
%
%          'nearest' - nearest neighbor interpolation. This will
%                     rarely be a good choice, but I included it
%                     as an option for completeness.
%
%
%   'regularizer' - character flag, denotes the regularization
%          paradignm to be used. There are currently three options.
%
%          DEFAULT: 'gradient'
%
%          'diffusion' or 'laplacian' - uses a finite difference
%              approximation to the Laplacian operator (i.e, del^2).
%              
%              We can think of the surface as a plate, wherein the
%              bending rigidity of the plate is specified by the user
%              as a number relative to the importance of fidelity to
%              the data. A stiffer plate will result in a smoother
%              surface overall, but fit the data less well. I've
%              modeled a simple plate using the Laplacian, del^2. (A
%              projected enhancement is to do a better job with the
%              plate equations.)
%              
%              We can also view the regularizer as a diffusion problem,
%              where the relative thermal conductivity is supplied.
%              Here interpolation is seen as a problem of finding the
%              steady temperature profile in an object, given a set of
%              points held at a fixed temperature. Extrapolation will
%              be linear. Both paradigms are appropriate for a Laplacian
%              regularizer.
%              
%          'gradient' - attempts to ensure the gradient is as smooth
%              as possible everywhere. Its subtly different from the
%              'diffusion' option, in that here the directional
%              derivatives are biased to be smooth across cell
%              boundaries in the grid.
% 
%              The gradient option uncouples the terms in the Laplacian.
%              Think of it as two coupled PDEs instead of one PDE. Why
%              are they different at all? The terms in the Laplacian
%              can balance each other.
%
%          'springs' - uses a spring model connecting nodes to each 
%              other, as well as connecting data points to the nodes
%              in the grid. This choice will cause any extrapolation
%              to be as constant as possible.
%              
%              Here the smoothing parameter is the relative stiffness
%              of the springs connecting the nodes to each other compared
%              to the stiffness of a spting connecting the lattice to
%              each data point. Since all springs have a rest length
%              of zero, any extrapolation will be minimized.
%
%          Note: I don't terribly like the 'springs' strategy.
%          It tends to drag the surface towards the mean of all
%          the data. Its been left in only because the paradigm
%          interests me.
%
%
%   'solver' - character flag - denotes the solver used for the
%          resulting linear system. Different solvers will have
%          different solution times depending upon the specific
%          problem to be solved. Up to a certain size grid, the
%          direct \ solver will often be speedy, until memory
%          swaps causes problems.
%
%          What solver should you use? Problems with a significant
%          amount of extrapolation should avoid lsqr. \ may be
%          best numerically for small smoothnesss parameters and
%          high extents of extrapolation.
%
%          Large numbers of points will slow down the direct
%          \, but when applied to the normal equations, \ can be
%          quite fast. Since the equations generated by these
%          methods will tend to be well conditioned, the normal
%          equations are not a bad choice of method to use. Beware
%          when a small smoothing paramter is used, since this will
%          make the equations less well conditioned.
% 
%          DEFAULT: 'normal'
%
%          '\' - uses matlab's backslash operator to solve the sparse
%                     system. 'backslash' is an alternate name.
%
%          'symmlq' - uses matlab's iterative symmlq solver
%
%          'lsqr' - uses matlab's iterative lsqr solver
%
%          'normal' - uses \ to solve the normal equations.
%
%
%   'maxiter' - only applies to iterative solvers - defines the
%          maximum number of iterations for an iterative solver
%
%          DEFAULT: min(10000,length(xnodes)*length(ynodes))
%
%
%   'extend' - character flag - controls whether the first and last
%          nodes in each dimension are allowed to be adjusted to
%          bound the data, and whether the user will be warned if
%          this was deemed necessary to happen.
%
%          DEFAULT: 'warning'
%
%          'warning' - Adjust the first and/or last node in
%                     x or y if the nodes do not FULLY contain
%                     the data. Issue a warning message to this
%                     effect, telling the amount of adjustment
%                     applied.
%
%          'never'  - Issue an error message when the nodes do
%                     not absolutely contain the data.
%
%          'always' - automatically adjust the first and last
%                     nodes in each dimension if necessary.
%                     No warning is given when this option is set.
%
%
% Arguments: (output)
%  zgrid   - (nx,ny) array containing the fitted surface
%
%  xgrid, ygrid - as returned by meshgrid(xnodes,ynodes)
%
%
% Speed considerations:
%  Remember that gridfit must solve a LARGE system of linear
%  equations. There will be as many unknowns as the total
%  number of nodes in the final lattice. While these equations
%  may be sparse, solving a system of 10000 equations may take
%  a second or so. Very large problems will benefit from the
%  iterative solvers.
%
%
% Example usage:
%
%  x = rand(100,1);
%  y = rand(100,1);
%  z = exp(x+2*y);
%  xnodes = 0:.1:1;
%  ynodes = 0:.1:1;
%
%  g = gridfit(x,y,z,xnodes,ynodes);
%
% Note: this is equivalent to the following call:
%
%  g = gridfit(x,y,z,xnodes,ynodes,'smooth',1, ...
%              'interp','triangle','solver','normal', ...
%              'regularizer','gradient','extend','warning');

% set defaults
params.smoothness = 1;
params.interp = 'triangle';
params.regularizer = 'gradient';
params.solver = 'normal';
params.maxiter = [];
params.extend = 'warning';

% and check for any overrides
params = parse_pv_pairs(params,varargin);

% check the parameters for acceptability
% smoothness == 1 by default
if isempty(params.smoothness)
  params.smoothness = 1;
else
  if (params.smoothness<=0)
    error 'Smoothness must be real, finite, and positive.'
  end
end
% regularizer  - must be one of 4 options - the second and
% third are actually synonyms.
valid = {'springs', 'diffusion', 'laplacian', 'gradient'};
if isempty(params.regularizer)
  params.regularizer = 'diffusion';
end
ind = strmatch(lower(params.regularizer),valid);
if (length(ind)==1)
  params.regularizer = valid{ind};
else
  error(['Invalid regularization method: ',params.regularizer])
end

% interp must be one of:
%    'bilinear', 'nearest', or 'triangle'
% but accept any shortening thereof.
valid = {'bilinear', 'nearest', 'triangle'};
if isempty(params.interp)
  params.interp = 'triangle';
end
ind = strmatch(lower(params.interp),valid);
if (length(ind)==1)
  params.interp = valid{ind};
else
  error(['Invalid interpolation method: ',params.interp])
end

% solver must be one of:
%    'backslash', '\', 'symmlq', 'lsqr', or 'normal'
% but accept any shortening thereof.
valid = {'backslash', '\', 'symmlq', 'lsqr', 'normal'};
if isempty(params.solver)
  params.solver = '\';
end
ind = strmatch(lower(params.solver),valid);
if (length(ind)==1)
  params.solver = valid{ind};
else
  error(['Invalid solver option: ',params.solver])
end

% extend must be one of:
%    'never', 'warning', 'always'
% but accept any shortening thereof.
valid = {'never', 'warning', 'always'};
if isempty(params.extend)
  params.extend = 'warning';
end
ind = strmatch(lower(params.extend),valid);
if (length(ind)==1)
  params.extend = valid{ind};
else
  error(['Invalid extend option: ',params.extend])
end

% ensure all of x,y,z,xnodes,ynodes are column vectors,
% also drop any NaN data
x=x(:);
y=y(:);
z=z(:);
k = isnan(x) | isnan(y) | isnan(z);
if any(k)
  x(k)=[];
  y(k)=[];
  z(k)=[];
end
xmin = min(x);
xmax = max(x);
ymin = min(y);
ymax = max(y);

% did they supply a scalar for the nodes?
if length(xnodes)==1
  xnodes = linspace(xmin,xmax,xnodes)';
  xnodes(end) = xmax; % make sure it hits the max
end
if length(ynodes)==1
  ynodes = linspace(ymin,ymax,ynodes)';
  ynodes(end) = ymax; % make sure it hits the max
end

xnodes=xnodes(:);
ynodes=ynodes(:);
dx = diff(xnodes);
dy = diff(ynodes);
nx = length(xnodes);
ny = length(ynodes);
ngrid = nx*ny;

% default for maxiter?
if isempty(params.maxiter)
  params.maxiter = min(10000,nx*ny);
end

% check lengths of the data
n = length(x);
if (length(y)~=n)|(length(z)~=n)
  error 'Data vectors are incompatible in size.'
end
if n<3
  error 'Insufficient data for surface estimation.'
end

% verify the nodes are distinct
if any(diff(xnodes)<=0)|any(diff(ynodes)<=0)
  error 'xnodes and ynodes must be monotone increasing'
end

% do we need to tweak the first or last node in x or y?
if xmin<xnodes(1)
  switch params.extend
    case 'always'
      xnodes(1) = xmin;
    case 'warning'
      warning(['xnodes(1) was decreased by: ',num2str(xnodes(1)-xmin),', new node = ',num2str(xmin)])
      xnodes(1) = xmin;
    case 'never'
      error(['Some x (',num2str(xmin),') falls below xnodes(1) by: ',num2str(xnodes(1)-xmin)])
  end
end
if xmax>xnodes(end)
  switch params.extend
    case 'always'
      xnodes(end) = xmax;
    case 'warning'
      warning(['xnodes(end) was increased by: ',num2str(xmax-xnodes(end)),', new node = ',num2str(xmax)])
      xnodes(end) = xmax;
    case 'never'
      error(['Some x (',num2str(xmax),') falls above xnodes(end) by: ',num2str(xmax-xnodes(end))])
  end
end
if ymin<ynodes(1)
  switch params.extend
    case 'always'
      ynodes(1) = ymin;
    case 'warning'
      warning(['ynodes(1) was decreased by: ',num2str(ynodes(1)-ymin),', new node = ',num2str(ymin)])
      ynodes(1) = ymin;
    case 'never'
      error(['Some y (',num2str(ymin),') falls below ynodes(1) by: ',num2str(ynodes(1)-ymin)])
  end
end
if ymax>ynodes(end)
  switch params.extend
    case 'always'
      ynodes(end) = ymax;
    case 'warning'
      warning(['ynodes(end) was increased by: ',num2str(ymax-ynodes(end)),', new node = ',num2str(ymax)])
      ynodes(end) = ymax;
    case 'never'
      error(['Some y (',num2str(ymax),') falls above ynodes(end) by: ',num2str(ymax-ynodes(end))])
  end
end

% only generate xgrid and ygrid if requested.
if nargout>1
  [xgrid,ygrid]=meshgrid(xnodes,ynodes);
end

% determine which cell in the array each point lies in
[junk,indx] = histc(x,xnodes);
[junk,indy] = histc(y,ynodes);
% any point falling at the last node is taken to be
% inside the last cell in x or y.
k=(indx==nx);
indx(k)=indx(k)-1;
k=(indy==ny);
indy(k)=indy(k)-1;

% interpolation equations for each point
tx = min(1,max(0,(x - xnodes(indx))./dx(indx)));
ty = min(1,max(0,(y - ynodes(indy))./dy(indy)));
ind = indy + ny*(indx-1);
% Future enhancement: add cubic interpolant
switch params.interp
  case 'triangle'
    % linear interpolation inside each triangle
    k = (tx > ty);
    L = ones(n,1);
    L(k) = ny;
    
    t1 = min(tx,ty);
    t2 = max(tx,ty);
    A = sparse(repmat((1:n)',1,3),[ind,ind+ny+1,ind+L], ...
       [1-t2,t1,t2-t1],n,ngrid);
    
  case 'nearest'
    % nearest neighbor interpolation in a cell
    k = round(1-ty) + round(1-tx)*ny;
    A = sparse((1:n)',ind+k,ones(n,1),n,ngrid);
    
  case 'bilinear'
    % bilinear interpolation in a cell
    A = sparse(repmat((1:n)',1,4),[ind,ind+1,ind+ny,ind+ny+1], ...
       [(1-tx).*(1-ty), (1-tx).*ty, tx.*(1-ty), tx.*ty], ...
       n,ngrid);
    
end
rhs = z;

% Build regularizer. Add del^4 regularizer one day.
switch params.regularizer
  case 'springs'
    % zero "rest length" springs
    [i,j] = meshgrid(1:nx,1:(ny-1));
    ind = j(:) + ny*(i(:)-1);
    m = nx*(ny-1);
    stiffness = 1./dy;
    Areg = sparse(repmat((1:m)',1,2),[ind,ind+1], ...
       stiffness(j(:))*[-1 1],m,ngrid);
    
    [i,j] = meshgrid(1:(nx-1),1:ny);
    ind = j(:) + ny*(i(:)-1);
    m = (nx-1)*ny;
    stiffness = 1./dx;
    Areg = [Areg;sparse(repmat((1:m)',1,2),[ind,ind+ny], ...
       stiffness(i(:))*[-1 1],m,ngrid)];
    
    [i,j] = meshgrid(1:(nx-1),1:(ny-1));
    ind = j(:) + ny*(i(:)-1);
    m = (nx-1)*(ny-1);
    stiffness = 1./sqrt(dx(i(:)).^2 + dy(j(:)).^2);
    Areg = [Areg;sparse(repmat((1:m)',1,2),[ind,ind+ny+1], ...
       stiffness*[-1 1],m,ngrid)];
    
    Areg = [Areg;sparse(repmat((1:m)',1,2),[ind+1,ind+ny], ...
       stiffness*[-1 1],m,ngrid)];
    
  case {'diffusion' 'laplacian'}
    % thermal diffusion using Laplacian (del^2)
    [i,j] = meshgrid(1:nx,2:(ny-1));
    ind = j(:) + ny*(i(:)-1);
    dy1 = dy(j(:)-1);
    dy2 = dy(j(:));
    
    Areg = sparse(repmat(ind,1,3),[ind-1,ind,ind+1], ...
      [-2./(dy1.*(dy1+dy2)), 2./(dy1.*dy2), ...
       -2./(dy2.*(dy1+dy2))],ngrid,ngrid);
    
    [i,j] = meshgrid(2:(nx-1),1:ny);
    ind = j(:) + ny*(i(:)-1);
    dx1 = dx(i(:)-1);
    dx2 = dx(i(:));
    
    Areg = Areg + sparse(repmat(ind,1,3),[ind-ny,ind,ind+ny], ...
      [-2./(dx1.*(dx1+dx2)), 2./(dx1.*dx2), ...
       -2./(dx2.*(dx1+dx2))],ngrid,ngrid);
    
  case 'gradient'
    % Subtly different from the Laplacian. A point for future
    % enhancement is to do it better for the triangle interpolation
    % case.
    [i,j] = meshgrid(1:nx,2:(ny-1));
    ind = j(:) + ny*(i(:)-1);
    dy1 = dy(j(:)-1);
    dy2 = dy(j(:));

    Areg = sparse(repmat(ind,1,3),[ind-1,ind,ind+1], ...
      [-2./(dy1.*(dy1+dy2)), 2./(dy1.*dy2), ...
      -2./(dy2.*(dy1+dy2))],ngrid,ngrid);

    [i,j] = meshgrid(2:(nx-1),1:ny);
    ind = j(:) + ny*(i(:)-1);
    dx1 = dx(i(:)-1);
    dx2 = dx(i(:));

    Areg = [Areg;sparse(repmat(ind,1,3),[ind-ny,ind,ind+ny], ...
      [-2./(dx1.*(dx1+dx2)), 2./(dx1.*dx2), ...
      -2./(dx2.*(dx1+dx2))],ngrid,ngrid)];

end
nreg = size(Areg,1);

% Append the regularizer to the interpolation equations,
% scaling the problem first. Use the 1-norm for speed.
NA = norm(A,1);
NR = norm(Areg,1);
A = [A;Areg*(params.smoothness*NA/NR)];
rhs = [rhs;zeros(nreg,1)];

% solve the full system, with regularizer attached
switch params.solver
  case {'\' 'backslash'}
    % permute for minimum fill in for R (in the QR)
    p = colamd(A);
    zgrid=zeros(ny,nx);
    zgrid(p) = A(:,p)\rhs;
    
  case 'normal'
    % The normal equations, solved with \. Can be fast
    % for huge numbers of data points.
    
    % Permute for minimum fill-in for \ (in chol)
    APA = A'*A;
    p = symamd(APA);
    zgrid=zeros(ny,nx);
    zgrid(p) = APA(p,p)\(A(:,p)'*rhs);
    
  case 'symmlq'
    % iterative solver - symmlq - requires a symmetric matrix,
    % so use it to solve the normal equations. No preconditioner.
    tol = abs(max(z)-min(z))*1.e-13;
    [zgrid,flag] = symmlq(A'*A,A'*rhs,tol,params.maxiter);
    zgrid = reshape(zgrid,ny,nx);
    
    % display a warning if convergence problems
    switch flag
      case 0
        % no problems with convergence
      case 1
        % SYMMLQ iterated MAXIT times but did not converge.
        warning(['Symmlq performed ',num2str(params.maxiter), ...
          ' iterations but did not converge.'])
      case 3
        % SYMMLQ stagnated, successive iterates were the same
        warning 'Symmlq stagnated without apparent convergence.'
      otherwise
        warning(['One of the scalar quantities calculated in',...
          ' symmlq was too small or too large to continue computing.'])
    end
    
  case 'lsqr'
    % iterative solver - lsqr. No preconditioner here.
    tol = abs(max(z)-min(z))*1.e-13;
    [zgrid,flag] = lsqr(A,rhs,tol,params.maxiter);
    zgrid = reshape(zgrid,ny,nx);
    
    % display a warning if convergence problems
    switch flag
      case 0
        % no problems with convergence
      case 1
        % lsqr iterated MAXIT times but did not converge.
        warning(['Lsqr performed ',num2str(params.maxiter), ...
          ' iterations but did not converge.'])
      case 3
        % lsqr stagnated, successive iterates were the same
        warning 'Lsqr stagnated without apparent convergence.'
      case 4
        warning(['One of the scalar quantities calculated in',...
          ' LSQR was too small or too large to continue computing.'])
    end
    
end

function params=parse_pv_pairs(params,pv_pairs)
% parse_pv_pairs: parses sets of property value pairs
% usage: params=parse_pv_pairs(default_params,pv_pairs)
%
% arguments: (input)
%  default_params - structure, with one field for every potential
%             property/value pair. Each field will contain the default
%             value for that property. If no default is supplied for a
%             given property, then that field must be empty.
%
%  pv_array - cell array of property/value pairs.
%             Case is ignored when comparing properties to the list
%             of field names. Also, any unambiguous shortening of a
%             field/property name is allowed.
%
% arguments: (output)
%  params   - parameter struct that reflects any updated property/value
%             pairs in the pv_array.

npv = length(pv_pairs);
n = npv/2;

if n~=floor(n)
  error 'Property/value pairs must come in PAIRS.'
end
if n<=0
  % just return the defaults
  return
end

if ~isstruct(params)
  error 'No structure for defaults was supplied'
end

% there was at least one pv pair. process any supplied
propnames = fieldnames(params);
lpropnames = lower(propnames);
for i=1:n
  pi = lower(pv_pairs{2*i-1});
  vi = pv_pairs{2*i};
  
  ind = strmatch(pi,lpropnames,'exact');
  if isempty(ind)
    ind = strmatch(pi,lpropnames);
    if isempty(ind)
      error(['No matching property found for: ',pv_pairs{2*i-1}])
    elseif length(ind)>1
      error(['Ambiguous property name: ',pv_pairs{2*i-1}])
    end
  end
  pi = propnames{ind};
  
  % override the corresponding default in params
  params = setfield(params,pi,vi);
  
end








