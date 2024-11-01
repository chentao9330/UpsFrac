function testReadGrdec()
% intended to test the functionality of readGrdecl
%
verbose=mrstVerbose();
readGRDECL('/data/ioNorne/norne2d.grdecl','verbose',verbose)
readGRDECL('/data/ioNorne/BC0407_IO.DATA','verbose',verbose)
end

