function fluid = initDeckADIFluid(deck)
reg = handleRegions(deck);
fluid = [];
%props
props = deck.PROPS;
fns = fieldnames(props);
for k = 1:numel(fns)
    fn   = fns{k}; 
    asgn = str2func(['assign',fn]);
    try
        fluid = asgn(fluid, props.(fn), reg);
    catch
        warning(['Could not assign property ', fn])
    end
end
fluid = assignRelPerm(fluid);
end    

