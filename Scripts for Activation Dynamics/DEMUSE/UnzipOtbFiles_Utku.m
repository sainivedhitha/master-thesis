clear all
[name, path]=uigetfile;
list=dir([path, '*.otb']);
for i=1:length(list),
    unzip([path, list(i).name],path);
    abs=xmlread([path, 'abstract.xml']);
    abstract=xmlwrite([path, 'abstract.xml']);
    ind1=strfind(abstract,'<protocol_code>');
    ind2=strfind(abstract,'</protocol_code>');
    newfilename=abstract(ind1+length(['<protocol_code>']):ind2-1);
    newfilename=[newfilename,'.sig'];
    if ~isempty(dir([path, newfilename])),
        delete([path, newfilename]);
    end
    n=list(i).name(1:end-4);
    movefile([path,n,'.sig'],[path, newfilename]);
    delete([path, n,'.xml'],[path, 'abstract.xml']);
end