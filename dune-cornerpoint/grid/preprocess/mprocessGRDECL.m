function G = mprocessGRDECL(grdecl, varargin)
%Compute grid topology and geometry from pillar grid description.
%
% SYNOPSIS:
%   G = mprocessGRDECL(grdecl)
%   G = mprocessGRDECL(grdecl, 'pn1', pv1, ...)
%
% PARAMETERS:
%   grdecl - Raw pillar grid structure, as defined by function
%            'readGRDECL', with fields COORDS, ZCORN and, possibly, ACTNUM.
%
%  'pn'/pv - List of 'key'/value pairs defining optional parameters.  The
%            supported options are:
%
%              Verbose -- Whether or not to display progress information
%                         Logical.  Default value: Verbose = false.
%
%              Tolerance --
%                         Minimum distinguishing vertical distance for
%                         points along a pillar.  Specifically, two points
%                         (x1,y1,z1) and (x2,y2,z2) are considered separate
%                         only if ABS(z2 - z1) > Tolerance.
%                         Non-negative scalar.
%                         Default value: Tolerance = 0.0 (distinguish all
%                         points along a pillar whose z coordinate differ
%                         even slightly).
%
%              CheckGrid --
%                         Whether or not to perform basic consistency
%                         checks on the resulting grid.
%                         Logical.  Default value: CheckGrid = true.
%
%              SplitDisconnected --
%                         Whether or not to split disconnected grid
%                         components into separate grids/reservoirs.
%                         Logical.  Default value: SplitDisconnected=true.
%
% RETURNS:
%   G      - Valid grid definition containing connectivity, cell
%            geometry, face geometry and unique nodes.
%
% EXAMPLE:
%   G = processGRDECL(readGRDECL('small.grdecl'));
%   plotGrid(G); view(10,45);
%
% SEE ALSO:
%   readGRDECL, deactivateZeroPoro, removeCells, processGRDECL

%{
#COPYRIGHT#
%}

% $Date: 2012-01-27 11:03:23 +0100 (Fri, 27 Jan 2012) $
% $Revision: 950 $
   opt = struct('Verbose', mrstVerbose, 'Tolerance', 0.0, 'CheckGrid', true, ...
             'SplitDisconnected', true);
   opt = merge_options(opt, varargin{:});
   
   if isfield(grdecl, 'ACTNUM'),
         grdecl.ACTNUM  = int32(grdecl.ACTNUM);      
   else
       grdecl.ACTNUM = int32(ones(prod(grdecl.cartDims),1));
   end    
   G = processgrid_mex(grdecl,opt.Tolerance);
   G.griddim = 3;
   if(opt.SplitDisconnected)
    G = splitDisconnectedGrid(G, false);   
   end
   
   if opt.CheckGrid,
      assert(all(diff(G.cells.facePos)>3));
      assert(all(diff(G.faces.nodePos)>2));
      assert(all(all(~isinf(G.nodes.coords))));
   end
 %{  
   if isfield(grdecl, 'MAPAXES'),
      for i = 1 : numel(G),
         G(i).nodes.coords(:,1:2) = ...
            mapAxes(G(i).nodes.coords(:,1:2), grdecl.MAPAXES);
      end
   end
 %}
end

function G = splitDisconnectedGrid(G, verbose)
   % Check if grid is connected
   ix = all(G.faces.neighbors~=0, 2);
   I  = [G.faces.neighbors(ix,1);G.faces.neighbors(ix,2)];
   J  = [G.faces.neighbors(ix,2);G.faces.neighbors(ix,1)];
   N  = double(max(G.faces.neighbors(:)));
   A  = sparse(double(I),double(J),1,N,N)+speye(N);
   clear ix I J
   [a,b,c,d]=dmperm(A);                                %#ok
   clear A b d
   if numel(c) > 2,
      dispif(verbose, '\nGrid has %d disconnected components\n', ...
         numel(c)-  1);
      % Partition grid into connected subgrids
      for i = 1:numel(c) - 1,
         g(i)  = extractSubgrid(G, a(c(i):c(i+1)-1));  %#ok
         sz(i) = g(i).cells.num;                       %#ok
         g(i).cartDims = G.cartDims;       %#ok
      end
      
      % Return largest (in number of cells) grid first
      [i,i] = sort(-sz);                                                 %#ok
      G     = g(i);
   end
end
