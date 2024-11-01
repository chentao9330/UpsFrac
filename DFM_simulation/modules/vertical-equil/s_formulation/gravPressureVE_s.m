function ff = gravPressureVE_s(g, omega)
% computes innerproduct cf (face_centroid - cell_centroid) * g for each face

   g_vec = gravity();   % Must be a 1-by-3 row vector for subsequent code.

   if norm(g_vec) > 0,       
      dim = size(g.nodes.coords,2);

      assert (1 < dim && dim < 4);
      assert (all(size(g_vec) == [1,3]));
      cellno = rldecode(1 : g.cells.num, diff(g.cells.facePos), 2) .';     
      
      if (any(strcmp(g.type, 'topSurfaceGrid')))
         % VE - 2D (i.e. g is created using topSurfaceGrid)
         cvec   = (g.faces.z(g.cells.faces(:,1), :) - ...
            g.cells.z(cellno            , :));
         ff     = omega(cellno) .* cvec*norm(gravity);
      else
         % VE - 3D
         cvec   = g.faces.centroids(g.cells.faces(:,1), :) - ...
                  g.cells.centroids(cellno            , :);
         cellF2D = (1:size(g.cells.faces,1))';
         
         
         hftb = any(bsxfun(@eq,g.cells.faces(cellF2D,2),[5,6]),2);
         % remove contribution from top and bottom, but only for 2D (VE)
         % cells i.e. not for real 3D cells (in the coupled version)         
         cvec(cellF2D(hftb),:) = 0;
         ff     = omega(cellno) .* (cvec * g_vec(1:dim).');
         
      end
   else
      ff     = zeros([size(g.cells.faces,1), 1]);
   end
end
