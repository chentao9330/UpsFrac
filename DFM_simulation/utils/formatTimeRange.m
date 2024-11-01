function s = formatTimeRange(time)
% Small utility which returns a human readable string from seconds.
    s = '';
    
    timesys = {year, day, hour, second};
    timen = {'Years', 'Days', 'Hours', 'Seconds'};
    added = false;
    for i = 1:numel(timesys)
        a = floor(time/timesys{i});
        if a > 0
            if added
                space = ', ';
            else
                space = '';
            end
            s = [s sprintf([space '%d ' timen{i}], a) ]; %#ok
            time = time - a*timesys{i};
            added = true;
        end
    end
end
