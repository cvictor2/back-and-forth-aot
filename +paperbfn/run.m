function [result,runDir] = run(c)
%RUN Run complete forward/backward cycles and persist enough data for replotting.
% Execute ONE configured experiment. Suites call this repeatedly for parameters.
% Call chain for NSE and 1D models
%   config/suite -> model -> reference -> step -> nonlinear -> errors
% Lorenz assimilation replaces repeated step calls with lorenz_leg.m.
% KdV uses the separate kdv_run.m driver.
%
% truth(:,j+1) is the reference at physical time j*dt. Each column is flattened
% for storage, then reshaped to m.shape before a PDE step. For NSE this is an
% N-by-N Fourier streamfunction. v is the current assimilated state.
% One cycle contains a complete forward pass and a complete backward pass.
% The estimate is carried between passes/cycles, not restarted each time.
% Validate settings once and derive stepsPerPass=window/dt. Invalid grids or
% an unstable isolated explicit feedback step are rejected before integration.
c=paperbfn.validate(c);
% Resume requires identical settings, current source text and input-file bytes.
% The pathname alone cannot identify a restart that may have been overwritten.
% A cached result must include an input manifest to establish data identity.
currentSources=paperbfn.source_manifest();
currentInputs=paperbfn.input_manifest(c);
if isfield(c.output,'resume') && c.output.resume && isfolder(c.output.root)
    saved=dir(fullfile(c.output.root,'*','result.mat'));
    current=c;current.output=rmfield(current.output,'resume');
    for file=saved'
        a=load(fullfile(file.folder,file.name),'result');previous=a.result.config;
        if isfield(previous.output,'resume'),previous.output=rmfield(previous.output,'resume');end
        sameSources=isfield(a.result,'sourceManifest') && isequaln(a.result.sourceManifest,currentSources);
        sameInputs=isfield(a.result,'inputManifest') && isequaln(a.result.inputManifest,currentInputs);
        finished=ismember(a.result.status,{'completed','tolerance_reached'}) && a.result.completedCycles==c.cycles;
        if isequaln(previous,current) && sameSources && sameInputs && finished
            result=a.result;runDir=file.folder;
            fprintf('Reusing completed run %s\n',c.name);return
        end
    end
end
% Build geometry/operators and evolve an unnudged reference over one window.
% reference also performs the configured physical spinup before defining t=0.
% The same reference window is reused on every forward/backward cycle.
m=paperbfn.model(c); [truth,reference,m]=paperbfn.reference(c,m);
u0=reshape(truth(:,1),m.shape);
% Choose the assimilated state at t=0. NSE defaults to the observed projection
% of the reference. That is an INITIAL guess only, not a projection after steps.
% "exact" is useful when testing backward instability from exact initial data.
switch c.initial.guessType
    case 'zero', v=zeros(m.shape);
    case 'projection', v=u0.*m.obs;
    case 'exact', v=u0;
    case 'explicit', v=c.initial.guess(:);
end
% reconstruction stores the last COMPLETED backward endpoint at physical t=0.
% v can later belong to a partial failed pass, so keep the two distinct.
initialGuess=v; reconstruction=v;
% Diagnostics are sampled every diagnosticStride steps, plus pass endpoints
% and failures. Sampling affects saved rows, not the number of solver steps.
% Columns 1:5 describe time/pass bookkeeping. Columns 6:11 are six norms.
% The initial row compares the initial guess with the reference at t=0.
maxRows=2*c.cycles*(ceil(c.stepsPerPass/c.output.diagnosticStride)+1)+1;
diagnostics=nan(maxRows,11); row=1;
% cycle coordinate, physical time, cycle number, leg (+1/-1), step, six errors
initialError=paperbfn.errors(v,u0,m);
diagnostics(row,:)=[0,0,0,0,0,initialError];
cycleErrors=nan(c.cycles,7); % cycle plus six errors
% Preallocate per-cycle endpoints. Only completed cycles enter these arrays.
% Never treat the last finite point of an incomplete backward pass as t=0.
cycleStates=cell(c.cycles,1); status='completed'; failure=''; completed=0;
% Outer BFN loop. The forward pass follows t=0 to window with h=+dt.
% The backward pass follows window to 0 with h=-dt. Both use old-time data.
for cycle=1:c.cycles
    for leg=1:2
        back=(leg==2); h=c.dt*(1-2*back);
        % Reference index 1 means physical t=0, index stepsPerPass+1 means t=window.
        % indices lists OLD state times. Each step advances to idx+sign(h).
        if back, indices=c.stepsPerPass+1:-1:2; else, indices=1:c.stepsPerPass; end
        % The optimized Lorenz path uses the same physical RK4 and explicit feedback,
        % with its time loop in lorenz_leg.m so large parameter sweeps stay practical.
        if strcmp(c.model,'lorenz')
            [v,legRows,nextIdx,bad]=paperbfn.lorenz_leg(v,h,truth,indices,m,c,cycle,leg);
            diagnostics(row+1:row+size(legRows,1),:)=legRows;row=row+size(legRows,1);
            if bad
                status='failed';failure=sprintf('Cycle %d, leg %d: nonfinite state/error or error limit.',cycle,leg);
            end
        else
        % Actual PDE time loop. j counts elapsed steps in this pass, not physical time.
        % For NSE each iteration updates an entire Fourier streamfunction array.
        for j=1:c.stepsPerPass
            idx=indices(j); nextIdx=idx+sign(h); physical=(idx-1)*c.dt;
            % u is the old-time observation source used for explicit feedback.
            % unext is the reference at the NEW physical time used to measure error.
            % Passing both makes the time convention visible. step.m ignores unext.
            u=reshape(truth(:,idx),m.shape); unext=reshape(truth(:,nextIdx),m.shape);
            % The complete numerical update happens here. Open step.m to see the
            % integrating factor, physical nonlinear update and separate nudging term.
            v=paperbfn.step(v,h,physical,u,unext,m,c,back,false);
            % Always compare against the matching physical time, including backward steps.
            % For NSE, errors converts streamfunction coefficients to velocity norms.
            % The stop guard is checked at every step, even if diagnostics are sampled.
            err=paperbfn.errors(v,unext,m);
            bad=any(~isfinite(v(:))) || any(~isfinite(err)) || max(err)>c.stop.maxNorm;
            if mod(j,c.output.diagnosticStride)==0 || j==c.stepsPerPass || bad
                row=row+1;
                % Plot time is completed-cycle coordinate, not physical PDE time. Each pass
                % occupies 1/2 unit, so forward endpoints are half-integers and backward
                % endpoints are integers. Physical time is saved separately in column 2.
                coordinate=(cycle-1)+(leg-1)/2+j/(2*c.stepsPerPass);
                diagnostics(row,:)=[coordinate,(nextIdx-1)*c.dt,cycle,sign(h),j,err];
            end
            if bad
                status='failed'; failure=sprintf('Cycle %d, leg %d, step %d: nonfinite state/error or error limit.',cycle,leg,j); break
            end
        end
        end
        if strcmp(status,'failed'), break; end
    end
    if strcmp(status,'failed'), break; end
    % Both passes succeeded. We are now back at physical t=0. Save this estimate
    % and its error against the reference initial state for endpoint tables.
    completed=cycle; reconstruction=v;
    cycleErrors(cycle,:)=[cycle,paperbfn.errors(v,u0,m)];
    if c.output.saveCycleStates, cycleStates{cycle}=v; end
    fprintf('%s: cycle %d/%d, initial-state L2 error %.6g\n',c.name,cycle,c.cycles,cycleErrors(cycle,2));
    % Optional early stop uses total L2 (Euclidean for Lorenz) at a completed cycle.
    % Without a tolerance, run the full requested number of cycles unless it fails.
    if ~isempty(c.stop.tolerance) && cycleErrors(cycle,2)<=c.stop.tolerance
        status='tolerance_reached'; break
    end
end
% Persist settings, source hashes, time histories and endpoint states so the
% result can be inspected/replotted without rerunning. finalState may be partway
% through a failed pass. reconstruction remains the last valid t=0 endpoint.
result.sourceManifest=currentSources; result.inputManifest=currentInputs;
result.config=c; result.reference=reference; result.status=status; result.failure=failure;
result.completedCycles=completed; result.diagnostics=diagnostics(1:row,:);
result.diagnosticColumns={'cycleCoordinate','physicalTime','cycle','direction','step', ...
    'L2total','L2observed','L2unobserved','H1total','H1observed','H1unobserved'};
result.cycleColumns={'cycle','L2total','L2observed','L2unobserved','H1total','H1observed','H1unobserved'};
result.cycleErrors=cycleErrors(1:completed,:); result.cycleStates=cycleStates(1:completed);
result.reconstruction=reconstruction; result.initialReference=u0; result.initialGuess=initialGuess; result.finalState=v;
result.finalReference=reshape(truth(:,nextIdx),m.shape);
result.initialError=initialError;
result.runtime=struct('matlab',version,'computer',computer,'created',char(datetime('now','Format','yyyyMMdd''T''HHmmss')));
result.normConvention='absolute integral L2; velocity L2/H1 for NSE; Euclidean for Lorenz';
result.cycleConvention='one cycle = one forward pass and one backward pass';
% Use an isolated output directory. MAT files retain numerical arrays, CSVs
% expose the same error histories and tables, and JSON records all settings.
% A failed status and its reason are written rather than hidden by plotting.
if ~isfolder(c.output.root), mkdir(c.output.root); end
safe=regexprep(c.name,'[^A-Za-z0-9_-]','_');
[~,token]=fileparts(tempname(c.output.root));
runDir=fullfile(c.output.root,[safe '_' char(datetime('now','Format','yyyyMMdd_HHmmss')) '_' token]);mkdir(runDir);
save(fullfile(runDir,'result.mat'),'result','-v7.3');
fid=fopen(fullfile(runDir,'config.json'),'w'); assert(fid>=0); cleaner=onCleanup(@()fclose(fid));
fwrite(fid,jsonencode(c),'char'); clear cleaner
writetable(array2table(result.diagnostics,'VariableNames',result.diagnosticColumns),fullfile(runDir,'diagnostics.csv'));
writetable(array2table(result.cycleErrors,'VariableNames',result.cycleColumns),fullfile(runDir,'cycle_errors.csv'));
if c.output.plot, paperbfn.plot_result(result,runDir); end
fprintf('Saved %s (%s)\n',runDir,status);
end
