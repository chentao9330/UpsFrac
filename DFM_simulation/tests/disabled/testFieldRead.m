% testFieldRead -- Test of reading Eclipse output fields

% Restart file. ASCII
filename = '../../branches/ioNorne/meisam-data/BC0407_IO.F0009';

if exist(filename, 'file') == 2,
   [fid, msg] = fopen(filename,'r');
   if fid < 0, error(msg); end

   %          readEclipseOutput
   disp('   readEclipseOutput')
   ticif(1)
   celldata1 = readEclipseOutput(filename);
   tocif(1)
   numcells = numel(celldata1.PRESSURE.values)  %#ok
   fprintf('PRESSURE: %f  %f', ...
           celldata1.PRESSURE.values(1), ...
           celldata1.PRESSURE.values(end))
   fprintf('SWAT: %f  %f', ...
           celldata1.SWAT.values(1), ...
           celldata1.SWAT.values(end))

   %          readSelectedEclipseOutput
   disp('readSelectedEclipseOutput')
   fieldnames = {'PRESSURE','SWAT'};
   ticif(1)
   celldata2 = readSelectedEclipseOutput(filename, fieldnames);
   tocif(1)

   numcells = numel(celldata2.PRESSURE.values) %#ok
   sprintf('PRESSURE: %f  %f', ...
           celldata2.PRESSURE.values(1), ...
           celldata2.PRESSURE.values(end))

   sprintf('SWAT: %f  %f', ...
           celldata2.SWAT.values(1), ...
           celldata2.SWAT.values(end))
else
   fprintf('File ''%s'' does not exist.  Test ignored\n', filename);
end
