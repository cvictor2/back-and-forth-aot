function verify_resume
%VERIFY_RESUME Confirm that reuse requires current settings, input bytes and code.
% Runs tiny NSE cases in an isolated temporary directory. No paper data or
% release sources are edited. An unchanged run should reuse its output folder.
% Replacing either initial data or forcing at the same filename must recompute.
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);
folder=tempname;mkdir(folder);closer=onCleanup(@()rmdir(folder,'s'));
restart=fullfile(folder,'restart.mat');[x,y]=ndgrid((0:31)*2*pi/32);
psi_hat_initial=fft2(cos(x)+.1*cos(6*y));f_hat=zeros(32);
save(restart,'psi_hat_initial','f_hat');
c=paperbfn.config('nse');c.name='resume_check';c.N=32;c.nonlinearCutoff=32/3;
c.observe.cutoff=5;c.initial.highMode=6;c.initial.restartFile=restart;
c.window=.002;c.cycles=1;c.output.plot=false;c.output.resume=true;c.output.root=folder;
[r1,p1]=paperbfn.run(c);[~,p2]=paperbfn.run(c);
assert(strcmp(p1,p2),'Identical inputs and source should reuse the completed run.');
% Change initial coefficients without changing the configured path.
psi_hat_initial=2*psi_hat_initial;save(restart,'psi_hat_initial','f_hat');
[r2,p3]=paperbfn.run(c);
assert(~strcmp(p1,p3) && isequal(r2.initialReference,psi_hat_initial));
assert(~strcmp(r1.inputManifest.sha256,r2.inputManifest.sha256));
% A forcing-only change also invalidates reuse.
f_hat=fft2(sin(x));save(restart,'psi_hat_initial','f_hat');
[r3,p4]=paperbfn.run(c);assert(~strcmp(p3,p4));
assert(~strcmp(r2.inputManifest.sha256,r3.inputManifest.sha256));
% Reject a saved result with an outdated solver identity. Earlier results
% have different input hashes, so none is a valid substitute for this case.
result=r3;result.sourceManifest(1).sha256=repmat('0',1,64);
save(fullfile(p4,'result.mat'),'result','-v7.3');
[r4,p5]=paperbfn.run(c);assert(~strcmp(p4,p5));
% A cached result without input hashes cannot establish data identity.
result=rmfield(r4,'inputManifest');save(fullfile(p5,'result.mat'),'result','-v7.3');
[~,p6]=paperbfn.run(c);assert(~strcmp(p5,p6));
% Missing data must raise an error instead of silently returning cached data.
delete(restart);missingRejected=false;
try,paperbfn.run(c);catch exception,missingRejected=strcmp(exception.identifier,'paperbfn:InputFile');end
assert(missingRejected);
% Binary SHA-256 must preserve values above 127 when passing bytes to Java.
file=fullfile(folder,'bytes.bin');fid=fopen(file,'wb');fwrite(fid,uint8([0 127 128 255]),'uint8');fclose(fid);
c.initial.restartFile=file;manifest=paperbfn.input_manifest(c);
assert(strcmp(manifest.sha256,'89273d2f70b93285bb7ddb4bcee86a5347ca7159352e3cbdd20c23e9d1e507d3'));
fprintf('Resume checks passed for unchanged input, changed initial data/forcing, source identity, incomplete manifests and missing files.\n');
end
