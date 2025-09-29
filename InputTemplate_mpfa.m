%% DFM simulation template file modified from the MRST code //  Part of package: UpsFrac
% A somewhat advanced example, where a fracture network is gridded
% with triangles using a fairly simple approach. The flow equations are
% then discretized with multi-point flux expressions, and a
% fluid flow problem is solved. The linear boundary condition is applied.

% Create an initial distribution of points
% N     = 100;
dx    = 100;
dy    = 100;
[X Y] = meshgrid(0:dx/20:dx, 0:dy/20:dy);%dx/20 and dy/20 are matrix grid size
X     = X;
Y     = Y;
p     = [X(:), Y(:)];

% Estimate of the grid size
h = .000001;

% Define endpoints of the fractures.
% If a nertex here is also contained in p, one of them should be removed to
% avoid mapping errors in the subsequent Delaunay triangulation
vertices = [ 
    
];






constraints = reshape(1 : size(vertices,1),2,[])' ;
constraints = [constraints, (1 : size(constraints,1))'];
box = [0 0; dx dy];
perimiter = h ;
p = remove_closepoints(vertices,constraints,p,perimiter);
args = struct('precision',0.00001);
[vertices, constraints] = removeFractureIntersections(vertices,constraints,box,args);
numOrdPt = size(p,1);
p  = [p ; vertices];
constraints(:,1:2) = constraints(:,1:2) + numOrdPt;
[p, constraints, map] = partition_edges(p,constraints,dx/30,box,args);;%dx/30 is fracture size
tags = constraints(:,3);
constraints =  constraints(:,1:2);
delTri     = DelaunayTri(p, constraints);
G     = triangleGrid(delTri.X, delTri.Triangulation);
G = computeGeometry(G);

G.faces.tags = zeros(G.faces.num,1);
faceNodes = sort(reshape(G.faces.nodes,2,[])',2);
constraints = sort(delTri.Constraints,2);
for iter = 1 : size(constraints)
    fracFace = find(ismember(faceNodes,constraints(iter,:),'rows'));
    G.faces.tags(fracFace) = tags(iter);
end
aperture = zeros(G.faces.num,1);

% Set different apertures for each
for iter = unique(tags)'
    hit = G.faces.tags == iter;
    aperture(hit) = leng_aper(iter,5);%aperture(hit) = 1.2*10^-4; % This choice is random in more than one way..
end

G = addhybrid(G,G.faces.tags > 0,aperture);
% 
% figure
% plotGrid_DFM(G)
% plotFractures(G)




% text(G.faces.centroids(:,1)-0.02, G.faces.centroids(:,2)-0.01, num2str((1:G.faces.num)'),'FontSize',8);
% 
% 
% 
% 
% 
% 
% title('subregion')
% axis equal, axis off



strsite = 'hybrid';
%mkdir RESULT;
filename = strcat('./RESULT/',strsite,subtxt);

fid=fopen(filename,'wt');

fprintf(fid,'%g\n',G.faces.hybrid);
fclose(fid)



strsite = 'centroids';
filename = strcat('./RESULT/',strsite,subtxt);
fid=fopen(filename,'wt');
[m,n]=size(G.faces.centroids);
for i=1:1:m
    for j=1:1:n   
        if j==n      
            fprintf(fid,'%g\n',G.faces.centroids(i,j));    
        else
            fprintf(fid,'%g\t',G.faces.centroids(i,j));
        end
    end
end
fclose(fid);










%% Set parameters

% Find indices of hybrid cells
hybridInd = find(G.cells.hybrid);
nCells = G.cells.num;

% Define permeability and porosity
rock.perm = 1*milli * darcy * ones(nCells,2);
rock.poro = 0*0.01 * ones(nCells,1);
rock.perm(hybridInd,:) = aperture(G.cells.tags(hybridInd)).^2/12 * [1 1];
rock.poro(hybridInd) = 0.5;

% Create fluid object. Non-linear rel perms, resident fluid has the higher
% viscosity
fluid = initSingleFluid('mu', 1*centi*poise, 'rho', 1014*kilogram/meter^3);

B=[G.faces.centroids]
% % display(B)
% 
l=find(abs(B(:,1)-0)<1e-6) %find the first(1) column: x, the value is 0
bl=l'
r=find(abs(B(:,1)-dx)<1e-6)
br=r'

d=find(abs(B(:,2)-0)<1e-6)
bd=d'
u=find(abs(B(:,2)-dy)<1e-6)
bu=u'




for i=1:2

prob=i
    
    
if prob==1

bc = addBC([],bl, 'pressure', dx);
bc = addBC(bc,br, 'pressure', 0);

bc = addBC(bc,bu, 'pressure', dx-G.faces.centroids(bu,1));
bc = addBC(bc,bd, 'pressure', dx-G.faces.centroids(bd,1));

state = initResSol(G,0,1);
T_MPFA = computeMultiPointTrans_DFM(G,rock,'hybrid',true);
[G,T_MPFA] = computeHybridMPTrans(G,T_MPFA);
state_MPFA = incompMPFA_DFM(state,G,T_MPFA,fluid,'bc',bc,'cellConnections',true);

% figure
% plotCellData_DFM(G,state_MPFA.pressure)
% plotFractures(G,hybridInd,state_MPFA.pressure)
% title('Problem 2')
% axis equal, axis off



strsite = 'ffluxpro1';
filename = strcat('./RESULT/',strsite,subtxt);

fid=fopen(filename,'wt');

fprintf(fid,'%30.25g\n',state_MPFA.flux);
fclose(fid)





else

bc = addBC([],bl, 'pressure', dy-G.faces.centroids(bl,2));
bc = addBC(bc,br, 'pressure', dy-G.faces.centroids(br,2));

bc = addBC(bc,bu, 'pressure', 0);
bc = addBC(bc,bd, 'pressure', dy);


state = initResSol(G,0,1);
T_MPFA = computeMultiPointTrans_DFM(G,rock,'hybrid',true);
[G,T_MPFA] = computeHybridMPTrans(G,T_MPFA);
state_MPFA = incompMPFA_DFM(state,G,T_MPFA,fluid,'bc',bc,'cellConnections',true);



% figure
% plotCellData_DFM(G,state_MPFA.pressure)
% plotFractures(G,hybridInd,state_MPFA.pressure)
% title('Problem 2')
% axis equal, axis off



strsite = 'ffluxpro2';
filename = strcat('./RESULT/',strsite,subtxt);

fid=fopen(filename,'wt');

fprintf(fid,'%30.25g\n',state_MPFA.flux);
fclose(fid)





end

end