function [qWs, qOs, qGs, bhp] = wellSolToVector(wellsols)
% Helper function which makes cell arrays of well solutions easier to plot
    nt = numel(wellsols);
    nw = numel(wellsols{1});
    ws = vertcat(wellsols{:});
    
    fix = @(v) reshape(v, nt, nw);
    sgn = fix([ws.sign]);
    bhp = fix([ws.pressure]);
    qWs = sgn.*fix([ws.qWs]);
    qGs = sgn.*fix([ws.qGs]);
    qOs = sgn.*fix([ws.qOs]);

end
