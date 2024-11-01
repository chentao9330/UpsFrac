% testBinRestartWrite -- Test of writing Eclipse Restart files in binary format.
function testBinRestartWrite()

format compact;

% Read ASCII restart files.
disp('READING ORIGINAL ASCII RESTART FILES')
prefix = '../../branches/testCase/FORMAT/ASCI/QFIVE_ASC';
for i = 1 : 2,
   filename = sprintf('%s.F%04d', prefix, i-1)
   restart{i} = readEclipseOutput(filename);
end

%{
% Read binary restart files.
disp('READING ORIGINAL BINARY RESTART FILES')
prefix = '../../branches/testCase/FORMAT/BINARY/QFIVE_BI';
for i = 1 : 2,
   filename = sprintf('%s.X%04d', prefix, i-1)
   restart{i} = readBinaryEclipseFile(filename);
end
%}


% Write new binary restart files.
disp('WRITING NEW BINARY FILES')
new_prefix = 'NEW_BINARY';
writeBinaryEclipseRestartData(new_prefix, restart);


% Read new binary restart files.
disp('READING NEW BINARY RESTART FILES')
for i = 1 : numel(restart),
  filename = sprintf('%s.X%04d', new_prefix, i-1)
  new_restart{i} = readBinaryEclipseFile(filename);
  disp_ecl(new_restart{i});
  %disp_ecl_all(new_restart{i});
end
end

% Display first and last values in the fields
function disp_ecl(ecl_fields)
disp('ecl_fields :');
fnames = fieldnames(ecl_fields);
 for i = 1 : numel(fnames),
    name = fnames{i}
    if numel(ecl_fields.(name).values) > 0,
        type  = ecl_fields.(name).type
        first = ecl_fields.(name).values(1)
        last  = ecl_fields.(name).values(end)
    end
 end
end

% Display all values in the fields
function disp_ecl_all(ecl_fields)
disp('ecl_fields :');
fnames = fieldnames(ecl_fields);
 for i = 1 : numel(fnames),
    name = fnames{i}
    if numel(ecl_fields.(name).values) > 0,
      ecl_fields.(name).type
      ecl_fields.(name).values'
   end
 end
end

