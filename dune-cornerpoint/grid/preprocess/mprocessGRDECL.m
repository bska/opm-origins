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

% $Date$
% $Revision$

   opt = struct('Verbose', mrstVerbose, 'Tolerance', 0.0, ...
                'CheckGrid', true, 'SplitDisconnected', true);
   opt = merge_options(opt, varargin{:});

   if isfield(grdecl, 'ACTNUM') && ...
         ~isa(grdecl.ACTNUM, 'int32'),

      grdecl.ACTNUM = int32(grdecl.ACTNUM);
   end

   require mex/libgeometry

   G = processgrid_mex(grdecl,opt.Tolerance);

   if opt.CheckGrid,
      assert(all(diff(G.cells.facePos)>3));
      assert(all(diff(G.faces.nodePos)>2));
      assert(all(all(~isinf(G.nodes.coords))));
   end

   if(opt.SplitDisconnected)
    G = splitDisconnectedGrid(G, false);
   end

   [ G(:).type    ] = deal({ mfilename });
   [ G(:).griddim ] = deal(3);

%{
   if isfield(grdecl, 'MAPAXES'),
      for i = 1 : numel(G),
         G(i).nodes.coords(:,1:2) = ...
            mapAxes(G(i).nodes.coords(:,1:2), grdecl.MAPAXES);
      end
   end
%}
end

%--------------------------------------------------------------------------

function G = splitDisconnectedGrid(G, verbose)
   % Check if grid is connected
   [a, c, c] = dmperm(adjacency(G));                                   %#ok

   ncomp = numel(c) - 1;
   if ncomp > 1,
      dispif(verbose, '\nGrid has %d disconnected components\n', ncomp);

      % Partition grid into connected subgrids
      g = arrayfun(@(i) extractSubgrid(G, a(c(i) : c(i + 1) - 1)), ...
                   1 : ncomp, 'UniformOutput', false);

      % Return grids in order of decreasing number of cells.
      cartDims = G.cartDims;
      [i, i]   = sort(- cellfun(@(g) g.cells.num, g));                 %#ok
      G        = [ g{i} ];

      [ G(:).cartDims ] = deal(cartDims);
   end
end

%--------------------------------------------------------------------------

function A = adjacency(G)
   N = double(G.faces.neighbors(~ any(G.faces.neighbors == 0, 2), :));
   I = [ N(:,1) ; N(:,2) ; (1 : G.cells.num) .' ];
   J = [ N(:,2) ; N(:,1) ; (1 : G.cells.num) .' ];
   A = sparse(I, J, 1, G.cells.num, G.cells.num);
end
