function manifest=source_manifest()
%SOURCE_MANIFEST SHA-256 of each MATLAB source used to identify a run's code.
% Record which MATLAB source text produced a new run. The paths are relative
% to the repository so results remain interpretable after moving the folder.
% This is provenance, not a numerical accuracy test or a reference trajectory.
% Comment edits intentionally change the hash even when arithmetic is identical.
root=fileparts(fileparts(mfilename('fullpath')));
% Fingerprint public root runners and shared solver/suite helpers. Result,
% figure and data files are separate and are not MATLAB source code here.
files=[dir(fullfile(root,'*.m'));dir(fullfile(root,'+paperbfn','*.m'))];
manifest=struct('file',{},'sha256',{});
for j=1:numel(files)
    path=fullfile(files(j).folder,files(j).name);
    % SHA-256 is computed from file text, including comments. The digest is stored
    % as a lowercase hexadecimal string alongside each source filename.
    md=java.security.MessageDigest.getInstance('SHA-256');
    md.update(uint8(fileread(path)));
    hash=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
    manifest(j).file=strrep(path,[root filesep],'');manifest(j).sha256=hash;
end
end
