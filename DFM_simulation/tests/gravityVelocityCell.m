function g = gravityVelocityCell(G,rock)
if(norm(gravity())>0)
  if( ( size(rock.perm,2)==1 || size(rock.perm,2)==size(G.nodes.coords,2) ))
      g     =  bsxfun(@times,rock.perm,gravity);      % Recall: 'grav' is a *ROW* vector.
   else
      error('gravity not jet implemented for non diagonal perm tensor')
  end
else
   g=0;
end
