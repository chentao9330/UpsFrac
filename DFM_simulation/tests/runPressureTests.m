function [G,solall,anasol,rell_err]=runPressureTests()
%% run test of all the test cases for the trid cases with different
%% pressure solvers
grid_cases=makeAnalyticCases();
%verbose=mrstVerbose();
%grid_cases=ioNorneAnalyticCases();


MODS = mrstModule;
mrstModule add mpfa


for i=1:numel(grid_cases)
   gcase=grid_cases{i};
   psolvers=makePressureSolvers(gcase.G,gcase.rock,gcase.fluid,[],'nonwetting');
   for j=1:numel(gcase.sim_cases)
      sim_case=gcase.sim_cases{j};
      solall=cell(numel(psolvers),1);
      rell_err=zeros(numel(psolvers),1);
      for k=1:numel(psolvers)
         solver=psolvers{k}.solver;
         sol = solver(sim_case.sol,sim_case.bc,sim_case.W);
         solall{k} = solver(sim_case.sol,sim_case.bc,sim_case.W);
         solall{k}.name=psolvers{k}.name;
         
         disp(['Solver ',psolvers{k}.name,' grid_case ', gcase.name, ' with  sim case ', sim_case.name ]);
         noninf=find(isfinite(sim_case.sol.pressure));
         diff = sol.pressure(noninf)-sim_case.sol.pressure(noninf);         
         relnorm = norm(diff)/norm(sim_case.sol.pressure(noninf));
         
         disp(['reltive error in pressure ', num2str(relnorm)])
         disp(['max diff ', num2str(max(abs(diff))./barsa) , ' max in reservoir ', num2str((max(sol.pressure)-min(sol.pressure))/barsa)])
         rell_err(k)=max(abs(diff))/(max(sol.pressure)-min(sol.pressure));
         if(relnorm>sim_case.reltol || ~isfinite(relnorm))
            disp('*****************************')
            disp(['Error for solver ',psolvers{k}.name,' grid_case ', gcase.name, ' with  sim case ', sim_case.name ])
            disp('*****************************')
            if(~(strcmp(psolvers{k}.name,'TPF') && strcmp(gcase.name,'cartgrid10x10x1twister')))
%               error('No linear accuracy of pressure solver')
            end
         end
      end
   end
end
%%
G=gcase.G;
anasol=sim_case.sol;

figure(1)
subplot(2,1,1)
cla,plotCellData(gcase.G,sol.pressure/barsa);ca=caxis();colorbar
subplot(2,1,2)
cla,plotCellData(gcase.G,(sim_case.sol.pressure)/barsa);caxis(ca);colorbar
figure(2)
cla,plotCellData(gcase.G,(sol.pressure(noninf)-sim_case.sol.pressure(noninf))/barsa,noninf);colorbar


mrstModule clear
mrstModule('add', MODS{:});

end

%--------------------------------------------------------------------------

function psolvs = makePressureSolvers(G, rock, fluid, src, pc_form)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
innerProducts = {'ip_simple', 'ip_quasirt'};
psolvs = {};
for i = 1 : numel(innerProducts),
   Sh = computeMimeticIP(G, rock, 'Verbose', mrstVerbose(), ...
                         'InnerProduct', innerProducts{i});

   psolvs{i}.solver = ...
      @(sol, bc, W) solveIncompFlow(sol, G, Sh, fluid, 'bc', bc,   ...
                                    'src', src, 'wells', W,        ...
                                    'MatrixOutput', mrstVerbose(), ...
                                    'pc_form', pc_form);  %#ok

   psolvs{i}.name = ['Hybrid_', innerProducts{i}];  %#ok
end

%inn_param=[1:1:10]
inn_param=[1,4,6];%,10,100]
%inn_param=[6]
for i = 1 : numel(inn_param),
   t  = inn_param(i);
   Sh = computeMimeticIP(G, rock, 'Verbose', mrstVerbose(), ...
                        'InnerProduct', 'ip_qfamily', 'qparam', t);

   psolvs{end+1}.solver = ...
      @(sol, bc, W) solveIncompFlow(sol, G, Sh, fluid, 'bc', bc,   ...
                                    'src', src, 'wells', W,        ...
                                    'MatrixOutput', mrstVerbose(), ...
                                    'pc_form', pc_form);  %#ok

   psolvs{end}.name=['Hybrid_', '_ip_dimfamily_', num2str(t)];  %#ok
end

%mixed solver
%%{
for i = 1 : numel(innerProducts),
   Sm = computeMimeticIP(G, rock, 'Verbose', mrstVerbose(), ...
                         'InnerProduct', innerProducts{i},  ...
                         'Type', 'mixed');

   psolvs{end+1}.solver = ...
      @(sol,bc,W) solveIncompFlow(sol, G, Sm, fluid, 'bc', bc, ...
                                  'src', src, 'wells', W,      ...
                                  'Solver', 'mixed',           ...
                                  'MatrixOutput', true);  %#ok

   psolvs{end}.name = ['Mixed_', innerProducts{i}];  %#ok
end

%inn_param=[1:1:10]
for i = 1 : numel(inn_param),
   t  = inn_param(i);
   Sm = computeMimeticIP(G, rock, 'Verbose', mrstVerbose(), ...
                         'InnerProduct', 'ip_qfamily',      ...
                         'Type', 'mixed', 'qparam', t);

   psolvs{end+1}.solver = ...
      @(sol, bc, W) solveIncompFlow(sol, G, Sm, fluid,   ...
                                    'bc', bc, 'src', src,  ...
                                    'wells', W, ...
                                    'Solver', 'mixed',     ...
                                    'MatrixOutput', true);  %#ok

   psolvs{end}.name = ['Mixed_', 'ip_dimfamily_', num2str(t)];  %#ok
end
%}

%tpf solver
Trans = computeTrans(G, rock, 'Verbose', mrstVerbose());

psolvs{end+1}.solver = ...
   @(sol, bc, W) incompTPFA(sol, G, Trans, fluid, 'bc', bc, 'src', src, ...
                            'wells', W, 'MatrixOutput', true, ...
                            'pc_form', pc_form);
psolvs{end}.name = 'TPF';

%mpfa solver
Trans_mpfa = computeMultiPointTrans(G, rock, 'Verbose', mrstVerbose());

psolvs{end+1}.solver = ...
   @(sol, bc, W) incompMPFA(sol, G, Trans_mpfa, fluid,  'bc', bc, ...
                            'src', src, 'wells', W, 'MatrixOutput', true);
psolvs{end}.name = 'MPFA';
end
