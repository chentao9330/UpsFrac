function isOscillating = detectNewtonOscillations(history, primary, current, tol)
    if current < 3
        isOscillating = false;
        return
    end
    
    tmp = history(current-2:current, primary);
    
    oscillate =  relChange(tmp(1,:), tmp(3,:)) < tol & ...
                 relChange(tmp(2,:), tmp(3,:)) > tol;
    
    isOscillating = sum(oscillate) > 1;
end

function v = relChange(a,b)
    v = abs((a-b)./b);
end
