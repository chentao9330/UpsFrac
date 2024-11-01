function grid_case = makeAnalyticCases()
grid_case={};
for k=1:2
G = cartGrid([30 30 1], [100 100 10]);
if(k==2)
   G.nodes.coords = twister(G.nodes.coords);
   grid_case{k}.name='cartgrid10x10x1twister';  %#ok
else
  grid_case{k}.name='cartgrid10x10x1';  %#ok
end
G = computeGeometry(G, 'Verbose', true);
grid_case{k}.G= G;  %#ok
rock.perm=ones(G.cells.num,1)*milli*darcy;
rock.poro=ones(G.cells.num,1)*0.2;
grid_case{k}.rock= rock;  %#ok
fluid = initSimpleFluid('mu' , [   1,  1]*centi*poise     , ...
                            'rho', [1014, 859]*kilogram/meter^3, ...
                            'n'  , [   2,   2]);
                         
grid_case{k}.fluid = fluid;  %#ok


   
%% constant velocity with pressure boudary conditon
%
veldir= [1 2 3]*barsa;
apressure =@(pos) sum((100e5+bsxfun(@times,veldir,pos)),2);
grid_case{k}.sim_cases={};  %#ok
initsat=1;
grid_case{k}.sim_cases{end+1}.sol =initResSol(G, initsat);  %#ok
grid_case{k}.sim_cases{end}.name = 'linear pressure at boundary';  %#ok
grid_case{k}.sim_cases{end}.sol.pressure=apressure(G.cells.centroids);  %#ok
grid_case{k}.sim_cases{end}.apressure=apressure;  %#ok
bfaces=find(sum(G.faces.neighbors==0,2)==1);
bc=addBC([],bfaces,'pressure',apressure(G.faces.centroids(bfaces,:)),'sat',ones(initsat,1));
W=[];
grid_case{k}.sim_cases{end}.bc=bc;  %#ok
grid_case{k}.sim_cases{end}.W=W;  %#ok
grid_case{k}.sim_cases{end}.reltol=1e-13;  %#ok
%% single well solution
%
W = verticalWell(W, G, rock,  floor(G.cartDims(1)/2),   floor(G.cartDims(1)/2), (1:G.cartDims(3)),     ...
                     'Type', 'rate', 'Val', 1e-5, ...
                     'Radius', 0.125, 'Name', 'P1');
%% find xy faces                  
bfx = find(G.cells.faces(:,2) < 5);
bfx = G.cells.faces(bfx,1);  %#ok
bf = bfx(find(any(G.faces.neighbors(bfx,:) == 0, 2)));  %#ok
cbfxy = find(G.cells.faces(:,2) < 5);
ci = find(any(G.faces.neighbors(G.cells.faces(cbfxy,1),:) == 0,2));
bfxy = G.cells.faces(cbfxy(ci),1);                 %#ok

[kr] = fluid.relperm(1);
mu              = fluid.properties();
tmob =sum(kr./mu,2);
lamperm = rock.perm(1)*tmob;

%% analytic solution for wells
clear pos;
W.dZ=10;
W.compi=[1 0];
apressure =@(pos) analytic2dSolutionWell(W,G,pos,lamperm);
grid_case{k}.sim_cases{end+1}.sol =initResSol(G, initsat);  %#ok
grid_case{k}.sim_cases{end}.name = 'one vertical well';  %#ok
grid_case{k}.sim_cases{end}.sol.pressure=apressure(G.cells.centroids);  %#ok
grid_case{k}.sim_cases{end}.apressure=apressure;  %#ok
bc=addBC([],bfxy,'pressure',apressure(G.faces.centroids(bfxy,:)),'sat',ones(initsat,1));
%bc=addBC([],bfxy,'pressure',repmat(100e5,numel(bfxy),1),'sat',ones(initsat,1));
grid_case{k}.sim_cases{end}.bc=bc;  %#ok
grid_case{k}.sim_cases{end}.W=W;  %#ok
grid_case{k}.sim_cases{end}.reltol=1e-4;  %#ok
end
end

function pressure = analytic2dSolutionWell(W,G,xyz,lamperm)
pressure = zeros(size(xyz,1),1);
nxyz=sqrt(xyz.^2);  %#ok
for i=1:numel(W)
    if(strcmp(W(i).type,'rate') &&  ~isempty(W(i).cells))
        rate = W(i).val/W(i).dZ;
        cell = W(i).cells(1);
        tmp = (xyz-repmat(G.cells.centroids(cell,:),size(xyz,1),1));
        nxyz=sqrt(sum(tmp.^2,2));
        pressure = pressure-(rate/(lamperm))*log(nxyz)/(2*pi)+100e5;       
    elseif(strcmp(W(i).type,'bhp'))
        error('hei');
    end
    
end
end
