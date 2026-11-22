function manifest=input_manifest(c)
%INPUT_MANIFEST SHA-256 identity of files actually loaded by this experiment.
% NSE loads its initial streamfunction and forcing from one restart MAT file.
% Other shared models construct data from settings, already compared by run.m.
% Hash the entire input file, not its name, size or modification timestamp.
% Replacing data at the same path therefore invalidates resume. Results
% require this manifest for reuse. Reference caching additionally hashes
% the loaded numerical arrays inside reference.m.
manifest=struct('file',{},'sha256',{});
if ~strcmp(c.model,'nse'),return;end
path=c.initial.restartFile;
fid=fopen(path,'rb');
assert(fid>=0,'paperbfn:InputFile','Cannot read NSE input file %s.',path);
closer=onCleanup(@()fclose(fid));
md=java.security.MessageDigest.getInstance('SHA-256');
while true
    bytes=fread(fid,1024*1024,'*uint8');
    if isempty(bytes),break;end
    % Reinterpret all eight bits when passing MATLAB bytes to Java signed bytes.
    md.update(typecast(bytes,'int8'));
end
hash=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
manifest=struct('file',path,'sha256',hash);
end
