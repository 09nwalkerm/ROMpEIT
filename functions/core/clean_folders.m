function clean_folders(direct,folder,tag)

% folder = e.g. '/Inverse_result_trad'

if isempty(tag)
    return
else
    cd([direct folder])
    delete([tag '.*'])
    cd ..
    delete([tag '.*'])
    cd ..
    delete([tag '.*'])
    cd ..
    delete([tag '.*'])
end