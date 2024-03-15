function output = loadJsonFile(pathJsonFile)
    %LOADJSONFILE load json file and save into a structure
    if exist(pathJsonFile,'file')
        json = jsondecode(fileread(pathJsonFile));
    else
        error("missing speed map analysis parameters file");
    end

    for fn = fieldnames(json)'
        try
            output.(fn{1}) = json.(fn{1});
        catch
            warning('Could not copy field %s', fn{1});
        end
    end   
end

