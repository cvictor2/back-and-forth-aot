function mask = observed_mask(m,c,t)
%OBSERVED_MASK Component windows are half-open in physical time.
% Return the coefficients/components receiving feedback at physical time t.
% For PDEs this is simply m.obs, the static Fourier mask built in model.m.
% For windowed Lorenz, m.obs names eligible components and this routine
% additionally turns their feedback on/off during their allocated time windows.
% It does not alter any state array. "Observed" is not a state projection.
mask=m.obs;
if strcmp(m.name,'lorenz') && strcmp(c.observe.schedule,'thirds')
    mask(:)=false;
    for j=c.observe.components
        % Component j gets the j-th third of the physical window. Its fraction
        % selects an initial part of that third. A zero fraction gives no observations
        % under the half-open rule used in the paper.
        left=(j-1)*c.window/3; right=left+c.observe.fractions(j)*c.window/3;
        % legacyClosed selects closed-right-end windows. The standard branch
        % includes the left endpoint and excludes the right endpoint.
        if isfield(c.lorenz,'windowRule') && strcmp(c.lorenz.windowRule,'legacyClosed')
            if j==1,leftTest=t>=left;else,leftTest=t>left;end
            if leftTest && t<=right,mask(j)=true;end
        elseif t>=left && t<right, mask(j)=true; end
    end
end
end
