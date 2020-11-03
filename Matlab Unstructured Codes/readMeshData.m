%% Reading the Mesh File

fid = fopen(meshfile_name,'r');

for k=1:4   % Ignoring starting lines
    cur_line = fgets(fid);
end

results         = sscanf(cur_line,'(%g (%g %g %x %g %g)(')';    % Reading the number nodes line
number_nodes    = results(4)-results(3)+1 ;

for k=1:2
    cur_line = fgets(fid);
end

vertices   	= fscanf(fid,'%f',[2 number_nodes])';

for k=1:3
    cur_line = fgets(fid);
end

results   	= sscanf(cur_line,'(%g (%g %g %x %g %g)(')';        % Reading the Total number of cells line

num_cells = results(4)-results(3)+1 ;

for k=1:2
    cur_line = fgets(fid);
end

results   	= sscanf(cur_line,'(%g (%g %g %x %g %g)(')';        % Reading the Total number of faces line
num_faces   = results(4)-results(3)+1 ;

for k=1:2
    cur_line = fgets(fid);
end

results   	= sscanf(cur_line,'(%g (%x %g %x %g %g)(')';        % Reading the number of Internal faces line i,e without boundaries

interior_faces  = results(4)-results(3)+1 ;

int_faces = fscanf(fid,'%x',[4 interior_faces])';

num_bcFaces = zeros(4,1) ;

% Reading Boundary Conditions

for num_bc = 1:4
    
    for k=1:4
        cur_line = fgets(fid);
    end
    result = strsplit(cur_line,{')','"','(',' ','"' });
    
    for k=1:1
        cur_line = fgets(fid);
    end
    
    results   	= sscanf(cur_line,'(%g (%x %x %x %g %g)(')';        % Reading number of boundary faces
    
    num_bcFaces(num_bc) = results(4)-results(3)+1 ;
    
    assignin('base',char(result(6)),fscanf(fid,'%x',[4 num_bcFaces(num_bc)])') ;
    
end

fclose(fid) ;       % Reading Mesh completed