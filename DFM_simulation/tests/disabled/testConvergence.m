function testConvergence(mycase,refdim)
if(nargin<2)
   mycase='qfs'
   refdim=[2,3,2];
end
%c = onCleanup(@()exportWorkspace);
%init fluid
fluid = initSimpleFluid('mu' , [1.0, 1.0] .* centi*poise    , ...
   'rho', [1000, 1000] .* kilogram/meter^3, ...
   'n'  , [1, 1]);
gravity off
% defining grdecl for setup
nx=4;ny=4;nz=1;
grdecl = simpleGrdecl([nx,ny,nz],0,'undisturbed',true)
sjakk=true;
switch mycase
   case 'qfs'
      rate=2+03;
      bhp=230
      %% set up quater five spot
      %middle producer
      compdat(5,:)={'P1',[ceil(nx/2)]    [ceil(ny/2)]    [1]    [nz]    'OPEN'    [-1]    [-1]    [0.2160]    [-1]    [0]    'Default'    'Z'    [-1]};
      wellspeck(5,:)={'P1'    'MANI-F'    [ceil(nx/2)]    [ceil(ny)]    'Default'    'OIL'    [0]    'STD'    'SHUT'    'NO'    [0]    'SEG'    [0]};
      wconprod(1,:) = {'P1'    'OPEN'    'BHP'     [Inf]    [Inf]    [Inf]    [ Inf]    [Inf]    int2str(bhp)    [0]    [0]    [0]};
      % injectors
      compdat(1,:)={'I1' [1]    [1]    [1]    [nz]    'OPEN'    [-1]    [-1]    [0.2160]    [-1]    [0]    'Default'    'Z'    [-1]};
      wellspeck(1,:)={'I1'    'MANI-F'    [-1]    [-1]    'Default'    'OIL'    [0]    'STD'    'SHUT'    'NO'    [0]    'SEG'    [0]};
      wconinje(1,:) = {'I1'    'WATER'    'OPEN'    'RESV'    [Inf]    [rate/4]    'Default'    [Inf]    [0]    [0]};
      %
      compdat(2,:)={'I2' [nx]    [1]    [1]    [nz]    'OPEN'    [-1]    [-1]    [0.2160]    [-1]    [0]    'Default'    'Z'    [-1]};
      wellspeck(2,:)={'I2'    'MANI-F'    [-1]    [-1]    'Default'    'OIL'    [0]    'STD'    'SHUT'    'NO'    [0]    'SEG'    [0]};
      wconinje(2,:) = {'I2'    'WATER'    'OPEN'    'RESV'    [Inf]    [rate/4]    'Default'    [Inf]    [0]    [0]};
      %
      compdat(3,:)={'I3' [1]    [ny]    [1]    [nz]    'OPEN'    [-1]    [-1]    [0.2160]    [-1]    [0]    'Default'    'Z'    [-1]};
      wellspeck(3,:)={'I3'    'MANI-F'    [-1]    [-1]    'Default'    'OIL'    [0]    'STD'    'SHUT'    'NO'    [0]    'SEG'    [0]};
      wconinje(3,:) = {'I3'    'WATER'    'OPEN'    'RESV'    [Inf]    [rate/4]    'Default'    [Inf]    [0]    [0]};
      %
      compdat(4,:)={'I4',[nx]    [ny]    [1]    [nz]    'OPEN'    [-1]    [-1]    [0.2160]    [-1]    [0]    'Default'    'Z'    [-1]};
      wellspeck(4,:)={'I4'    'MANI-F'    [ceil(nx/2)]    [ceil(ny)]    'Default'    'OIL'    [0]    'STD'    'SHUT'    'NO'    [0]    'SEG'    [0]};
      wconinje(4,:) = {'I4'    'WATER'    'OPEN'    'RESV'    [Inf]    [rate/4]    'Default'    [Inf]    [0]    [0]};
      %
   case 'one_well'
      rate=1e3
      compdat(1,:)={'I1',[ceil(nx/2)]    [ceil(ny/2)]    [1]    [nz]    'OPEN'    [0]    [0]    [0.2160]    [0]    [0]    'Default'    'Z'    [-1]};
      wellspeck(1,:)={'I1'    'MANI-F'    [ceil(nx/2)]    [ceil(ny/2)]    'Default'    'OIL'    [0]    'STD'    'SHUT'    'NO'    [0]    'SEG'    [0]};
      wconinje(1,:) = {'I1'    'WATER'    'OPEN'    'RESV'    [Inf]    [rate]    'Default'    [Inf]    [0]    [0]};
      %compdat(2,:)={'I2' [nx]    [1]    [1]    [nz]    'OPEN'    [-1]    [-1]    [0.2160]    [-1]    [0]    'Default'    'Z'    [-1]}
      %wellspeck(2,:)={'I2'    'MANI-F'    [-1]    [-1]    'Default'    'OIL'    [0]    'STD'    'SHUT'    'NO'    [0]    'SEG'    [0]}
      %wconinje(1,:) = {'I2'    'WATER'    'OPEN'    'RESV'    [Inf]    [rate/4]    'Default'    [Inf]    [0]    [0]}
      wconprod=[];
   case 'linear'
      wconprod=[];
      wellspeck=[];
      compdat=[];
      wconinje=[];
   otherwise
      error('no case defined')
end
if(~isempty(wellspeck))
   grdecl.SCHEDULE=struct('control',[]);
   grdecl.SCHEDULE.control.COMPDAT=compdat;
   grdecl.SCHEDULE.control.WELSPECS=wellspeck;
   grdecl.SCHEDULE.control.WCONINJE=wconinje;
   grdecl.SCHEDULE.control.WCONPROD=wconprod;
end
grdecl.METRIC=1;
grdecl.PORO=ones(nx,ny,nz)*0.1;
grdecl.PORO=grdecl.PORO(:);
grdecl.PERMX=ones(nx,ny,nz)*100*milli*darcy;
if(sjakk)
   grdecl.PERMX(1:2:end,1:2:end,1:2:end)=0.1*grdecl.PERMX(1:2:end,1:2:end,1:2:end);
   grdecl.PERMX(2:2:end,2:2:end,1:2:end)=0.1*grdecl.PERMX(2:2:end,2:2:end,1:2:end);
end
grdecl.PERMX=grdecl.PERMX(:);
grdecl.PERMY=grdecl.PERMX;
grdecl.PERMZ=0.1*grdecl.PERMX;

grdecl_old=grdecl;
%%
%%
grdecl_old=refineGrdecl(grdecl,refdim);% init well and grid
test_upscaling=true;
if(test_upscaling)
   grdecl_old=refineGrdecl(grdecl,refdim);% init well and grid
   g_old=processGRDECL(grdecl_old);   
   g_old=computeGeometry(g_old);
   rock_old=grdecl2Rock(grdecl_old);
   %grdecl_coarse=coarseGrdecl(grdecl_old,grdecl.cartDims,'perm_up_scaling'
   %,'flowbased')
   grdecl_coarse=coarseGrdecl(grdecl_old,refdim.*[2 2 1],'perm_up_scaling','flowbased')
   grdecl_coarse=coarseGrdecl(grdecl_old,refdim.*[2 2 1],'perm_up_scaling','simple')
   g_coarse=processGRDECL(grdecl_coarse);
   g_coarse=computeGeometry(g_coarse);
   rock_coarse=grdecl2Rock(grdecl_coarse);
   figure(1),clf
   subplot(2,2,1),
   plotGrid(g_old);view([1 1 1])
   subplot(2,2,2),
   plotGrid(g_coarse);view([1 1 1])
   subplot(2,2,3),
   plotCellData(g_old,log(rock_old.perm(:,1)));view([1 1 1])
   ca=caxis();
   subplot(2,2,4),
   plotCellData(g_coarse,log(rock_coarse.perm(:,1)));view([1 1 1]);caxis(ca);colorbar
   %%   
   disp(['Coarse rock perm X min ',num2str(min(rock_coarse.perm(:,1)/(milli*darcy))),' max ', num2str(max(rock_coarse.perm(:,1)/(milli*darcy)))])
   disp(['Coarse rock perm Y min ',num2str(min(rock_coarse.perm(:,2)/(milli*darcy))),' max ', num2str(max(rock_coarse.perm(:,2)/(milli*darcy)))])
   disp(['Coarse rock perm Z min ',num2str(min(rock_coarse.perm(:,3)/(milli*darcy))),' max ', num2str(max(rock_coarse.perm(:,3)/(milli*darcy)))])
   %keyboard()
   return   
end


g=processGrdecl(grdecl);
g=computeGeometry(g);
rock=grdecl2Rock(grdecl);


inner_product='ip_simple';
if(isfield(grdecl,'SCHEDULE'))
   W_mim=processWells(g,rock,grdecl.SCHEDULE.control,grdecl.METRIC.','InnerProduct',inner_product,'verbose',true);
   W_tpf=processWells(g,rock,grdecl.SCHEDULE.control,grdecl.METRIC.','InnerProduct','ip_tpf','verbose',true);
   for i=1:numel(W_mim)
      W_mim(i).compi=W_mim(i).compi(1:2);
      W_tpf(i).compi=W_tpf(i).compi(1:2);
   end
   
else
   W_mim=[];
   W_tpf=[];
end
W_mpfa=W_tpf;

state = initState(g, W_mim, 0, [0, 1]);

%init rock
rock=grdecl2Rock(grdecl);



%%plot configuration
%
figure(1),clf
plotCellData(g,rock.perm(:,1))
if(~isempty(W_mim))
   W=W_mim;
   plotWell(g, W, 'height', 200);
   nIW=find([[W(:).val]>0 ] & [strcmp({W(:).type},'rate')]);
   nPW=find(strcmp({W(:).type},'bho'));
   
   plotGrid(g, vertcat(W(nIW).cells), 'FaceColor', 'b');
   plotGrid(g, vertcat(W(nPW).cells), 'FaceColor', 'r');
end

%plotFaults(grdecl, G, 'height', 100, 'fontsize', 1, 'color', 'g')
%%
% preprocessing for TPFA,MPFA, mimetic
verbose=true;
Trans_mpfa   = computeMultiPointTrans(g, rock);
S    = computeMimeticIP(g, rock, 'Verbose', verbose,'InnerProduct',inner_product);
Trans    = computeTrans(g, rock, 'Verbose', verbose,'grdecl',grdecl);
if(strcmp(mycase,'one_well'))
   bfx = find(g.cells.faces(:,2) < 5);
   bfx = g.cells.faces(bfx,1);
   bf = bfx(find(any(g.faces.neighbors(bfx,:) == 0, 2)));
   cbfxy = find(g.cells.faces(:,2) < 5);
   ci = find(any(g.faces.neighbors(g.cells.faces(cbfxy,1),:) == 0,2));
   bfxy = g.cells.faces(cbfxy(ci),1);
   kr=fluid.relperm(0);
   mu=fluid.properties(state);
   lt=sum(kr/mu);
   lamperm = rock.perm(1,1)*lt*2;
   bfc = g.faces.centroids(bfxy,:);
   bpress =  analytic2dSolutionWell(W_mim,rock,g,bfc,lamperm);
   shift=mean(bpress)-100*barsa;
   bpress=bpress-shift;
   refpress = analytic2dSolutionWell(W_mim,rock,g,g.cells.centroids,lamperm)-shift;
   ref_well= analytic2dSolutionWell(W_mim,rock,g,g.cells.centroids(W_mim(:).cells,:)+[W_mim(:).r,0,0],lamperm)-shift;
   if(sjakk)
      bc = addBC([],bfxy,'pressure',100*barsa,'sat',1)';
   else      
      bc = addBC([],bfxy,'pressure',bpress,'sat',1)';
   end
elseif(strcmp(mycase,'linear'))
   bc=pside([],g,'XMin',200*barsa);
   bc=pside(bc,g,'XMax',100*barsa);
else
   bc=[];
end
condition_number=false;
disp('Calculating mimetic')
tic,state_mim  = solveIncompFlow(state, g, S, fluid,'wells',W_mim,'bc',bc,'condition_number',condition_number);toc
disp('Calculating tpfa')
tic,state_tpfa = incompTPFA(state, g, Trans, fluid,'wells',W_tpf,'bc',bc,'condition_number',condition_number);toc
disp('Calculating mpfa')
tic,state_mpfa = incompMPFA(state, g, Trans_mpfa, fluid,'wells',W_mpfa,'bc',bc,'condition_number',condition_number);toc
%    [W,wells_ok] = fixWells(W,wSol)
%wells_ok = true

%%plot results
%%
figure(2),clf
subplot(1,4,1)
plotCellData(g,state_mim.pressure/barsa);colorbar
ca=caxis();
subplot(1,4,2)
plotCellData(g,state_tpfa.pressure/barsa);colorbar;caxis(ca)
subplot(1,4,3)
plotCellData(g,state_mpfa.pressure/barsa);colorbar;;caxis(ca)
%subplot(1,4,4)
%plotCellData(g,refpress/barsa);colorbar;caxis(ca)
%if(strcmp(mycase,'one_well'))
   refpress=state_tpfa.pressure;
   %refpress(W_mim(:).cells)=nan;
   figure(3),clf
   subplot(1,3,1)
   plotCellData(g,(state_mim.pressure-refpress)/barsa);colorbar
   ca=caxis();
   subplot(1,3,2)
   plotCellData(g,(state_tpfa.pressure-refpress)/barsa);colorbar;caxis(ca)
   subplot(1,3,3)
   plotCellData(g,(state_mpfa.pressure-refpress)/barsa);colorbar;;caxis(ca)
if(strcmp(mycase,'one_well'))
   disp(['Ref well ',num2str(ref_well/barsa)])
   disp(['BHP mim ',num2str(state_mim.wellSol.pressure/barsa)])
   disp(['BHP mpfa ',num2str(state_mpfa.wellSol.pressure/barsa)])
   disp(['BHP tpfa ',num2str(state_tpfa.wellSol.pressure/barsa)])
elseif(strcmp(mycase,'linear'))
   disp(['Flux mim ', num2str(sum(state_mim.flux(bc.face).*(g.faces.neighbors(bc.face,1)==0)))])
   disp(['Flux tpfa ', num2str(sum(state_tpfa.flux(bc.face).*(g.faces.neighbors(bc.face,1)==0)))])
   disp(['Flux mpfa ', num2str(sum(state_mpfa.flux(bc.face).*(g.faces.neighbors(bc.face,1)==0)))])
end
%end
%%
end
function pressure = analytic2dSolutionWell(W,rock,G,xyz,lamperm)
pressure = zeros(size(xyz,1),1);
nxyz=sqrt(xyz.^2);
for i=1:numel(W)
    if(strcmp(W(i).type,'rate') &  ~isempty(W(i).cells))
        rate = W(i).val/(max(G.faces.centroids(:,3))-min(G.faces.centroids(:,3)));
        cell = W(i).cells(1)
        tmp = (xyz-repmat(G.cells.centroids(cell,:),size(xyz,1),1));
        nxyz=sqrt(sum(tmp.^2,2));
        pressure = pressure-(rate/(lamperm))*log(nxyz)/(2*pi);       
    elseif(strcmp(W(i).type,'bhp'))
        error('hei');
    end    
end
end
