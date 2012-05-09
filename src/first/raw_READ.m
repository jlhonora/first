%**********************************************************
% Raw data reader. 
%
% Shaihan Malik 15-03-05
%
% Raw reader for .list .data files. This reader is not as fast as it could
% be because it reads data only one kx vector at a time. This is deliberate
% since it was designed with kt under-sampled data in mind. The matrix is
% constrained to have only 5 dimensions at present and data with more will
% be rejected. 
%
% Have tested on some multi slice multi coil data and some volume data and
% all appears to be working correctly. Further work: average extra
% dimensions. Phase correction not implemented. After reading data must be
% ifft'd in kx,ky (and kz if applicable) and then shifted by half a FOV in
% x (not sure why this is).
%
% Modified by JLH Nov. 2011 to include the reading of FRX, NOI and PHX.
% Their usage still needs to be cleared. 
%
%***********************************************************

function [kspace, petable, frx, noi, phx] = raw_READ(filename, varargin)
type = [];
NoConsoleOutput = false;
try
    if(length(varargin)>0)
        if(strcmp(lower(varargin{1}),'noconsoleoutput')==1);
            NoConsoleOutput = true;
        end
    end
catch
end
% [kspace, petable] = raw_READ_STD(filename,type);
% [frx, petable2] = raw_READ_FRX(filename,type);
% phx = [];
% %[phx, petable3] = raw_READ_PHX(filename,type);
% [noi, petable4] = raw_READ_NOI(filename,type);

[kspace, petable] = raw_READ_SOME(filename,type, 'STD', NoConsoleOutput);
[frx, petable2] = raw_READ_SOME(filename,type, 'FRX', NoConsoleOutput);
phx = [];
%[phx, petable3] = raw_READ_SOME(filename,type, 'PHX', NoConsoleOutput);
[noi, petable4] = raw_READ_SOME(filename,type, 'NOI', NoConsoleOutput);

end

function [kspace, petable] = raw_READ_SOME(filename,type, DataString, NoConsoleOutput)

headername = [filename '.list'];
dataname = [filename '.data'];

% first look in the list file to find the number of dynamics, dimensions
% etc. do this by looking for lines starting with a '.' At present this
% will not work for the kx,ky,kz sizes as they have 2 entries and all
% others have one. Make a cell array to hold the variables. There are 8
% main variables (not including number of mixes), and if there are 
% multiple coils then there will be n_loc more lines giving coil combinations per slice

list_id = fopen(headername,'r','l');
if (list_id==-1)
    disp('Error: file not found or something');
    return
end

line_counter=1;
while (line_counter<=9)
    line = fgetl(list_id);
    if (strcmp(line(1),'.'))
        s = sscanf(line,'%*s %*s %*s %*s %s %s %s');
        % get the string with colon and use this to divide it
        div = strfind(s,':');
        % get the string above and store in cell array
        attributes{line_counter,1} = s(1:div-1);
        attributes{line_counter,2} = str2num(s(div+1:length(s)));
        line_counter = line_counter + 1;
    end
end

% attributes row num.   variable
%       1               mixes
%       2               spatial dimensions
%       3               dynamics
%       4               cardiac phases
%       5               echoes
%       6               locations (slices)
%       7               extra 1
%       8               extra 2
%       9               averages

% Now see if there are any extra coils .. if there are then the next line
% starting with '.' should be for receiver coils
done = false;
while (~done)
    line = fgetl(list_id);
    if (strcmp(line(1),'.'))
        % attempt to find the string 'coil_channel_combination'
        found = strfind(line,'coil_channel_combination');
        if (length(found)>0)
            multiple_coils = true;
        else
            multiple_coils = false;
        end
        done = true;
    end
end

% If there are multiple coils, which ones? Assume that the coils used are
% the same for each slice.
if multiple_coils
    div = strfind(line,':');
    % assume there are 10 coils (0:9) and there is a 1 if they are used and
    % a zero otherwise
    coils = sscanf(line((div+1):length(line)),'%i %i %i %i %i %i %i %i %i %i');
end

% There isn't really any need right now for the cell array... Store the
% variables in the order above in a numeric array
for i=1:9,
    data_dimensions(i) = cell2mat(attributes(i,2));
end

% complain if number of mixes is not 1
if (data_dimensions(1) ~= 1)
    disp('Error: number of mixes is greater than one, case not handled by code');
    return
end

% now read the .list file to get the phase encode table.
done=false;
rowcount=1;
while (~done)
    line = fgetl(list_id);
    if(line == -1)
        done=true;
        continue;
    end

    if ((line(1)=='#')||(line(1)=='.'))
        continue
    end

    % now we want to loop through lines beginning STD....
    linelength = numel(line);
    tag = strtok(line);
    % to read in NOI or PHX for example replace following with a case switch
    if (tag == DataString) %%%%%%%%%%%%%% !!!!!!!!!!!!!!!!!!!! %%%%%%%%%%%
        row = sscanf(line(6:linelength)','%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f' );
        petable(rowcount,:) = row;
        rowcount = rowcount+1;
    end
end
fclose(list_id);

% the 19 columns that are read in are:
% mix   dyn   card  echo  loca  chan  extr1 extr2 ky    kz    aver  sign
% rf    grad  enc   rtop  rr    size   offset
% Not sure how to deal with all possible eventualities. For now I'll
% try to do it generally. Matrix dimension can't be arbitrarily high so
% this is a constraint.

% don't want to cause more than a 5D matrix (inc x, 4D+coils)
% N = number of dimensions to read not including x (inc coils)
N = numel(find(data_dimensions>1));
% account for 3D case
if (data_dimensions(2)>2)
    N=N+1;
end

if multiple_coils
    if(~NoConsoleOutput)
        disp('Found multiple channel data')
    end
    N = N+1;
end


if (N>5)
    disp('Error too many dimensions, case not handled by code');
    return
else
    if(~NoConsoleOutput)
        disp(['Data matrix will have ' num2str(N+1) ' dims']);
    end
end


% store the contributing dimension column numbers in a vector. 
% this is a bit inelegant but never mind
nonzero = find(data_dimensions>1);
if(~NoConsoleOutput)
    nonzero
end
counter =1;
if(~NoConsoleOutput)
    disp('dimension 1 = kx');
end
for i=1:length(nonzero),
    switch nonzero(i)
        case 1
            disp('Error number of mixes greater than 1');
        case 2
            % this should always be one
            dims(counter)=9; %ky entry
            if(~NoConsoleOutput)
                disp(['dimension ',num2str(counter+1),' = ky']);
            end
            counter = counter+1;
            if (data_dimensions(2)>2)
                %have 2 PE dirs, so also read kz
                dims(counter)=10;
                if(~NoConsoleOutput)
                    disp(['dimension ',num2str(counter+1),' = kz']);
                end
                counter=counter+1;
            end
        case 3
            %dynamic data
             dims(counter)=2;
             if(~NoConsoleOutput)
                disp(['dimension ',num2str(counter+1),' = dynamic']);
             end
             counter = counter +1;
        case 4
            %multiple cardiac phases
            dims(counter)=3;
            if(~NoConsoleOutput)
                disp(['dimension ',num2str(counter+1),' = cardiac phase']);
            end
            counter=counter+1;
        case 5
            %multiple echoes
            dims(counter)=4;
            if(~NoConsoleOutput)
                disp(['dimension ',num2str(counter+1),' = echo number']);
            end
            counter=counter+1;
        case 6
            %multi slice data
            dims(counter)=5;
            if(~NoConsoleOutput)
                disp(['dimension ',num2str(counter+1),' = slice']);
            end
            counter = counter+1;
        case 7
            %extra 1
            dims(counter)=7;
            if(~NoConsoleOutput)
                disp(['dimension ',num2str(counter+1),' = extra 1']);
            end
            counter=counter+1;
        case 8
            %extra 2
            dims(counter)=8;
            if(~NoConsoleOutput)
                disp(['dimension ',num2str(counter+1),' = extra 2']);
            end
            counter=counter+1;
        case 9
            %multi averages
            dims(counter)=11;
            if(~NoConsoleOutput)
                disp(['dimension ',num2str(counter+1),' = averages']);
            end
            counter=counter+1;
    end
end
if multiple_coils
    dims(N)=6;
    if(~NoConsoleOutput)
        disp(['dimension ',num2str(counter+1),' = coil channels']);
    end
end

if (length(dims)>N)
    % must be a problem
    display('Error with dimensions ....');
    return
end
    

% Have N dimensions plus x to worry about.
recsize=4;
xsize = petable(1,18)/(2*recsize);
% hopefully all kx-vectors are the same length
for i=1:N,
    dim_min(i) = min(petable(:,dims(i)));
    dim_max(i) = max(petable(:,dims(i)));
    dim_size(i)= dim_max(i)+abs(dim_min(i))+1;
end

% the above estimation doesn't work for multiple coil data so correct it
if multiple_coils
    dim_size(N) = dim_max(N);
    % doing it this way bypasses the need to use the info taken from the
    % header file
end

kspace = complex(zeros([xsize dim_size],'single'));


% now open the data file and start reading
data_id = fopen(dataname,'r','l');
num_entries = numel(petable(:,1));

if(~NoConsoleOutput)
    disp('Beginning binary read ... ');
end
for i=1:num_entries,

    offset = petable(i,19);
    % read in the x vector
    fseek(data_id,offset,'bof');
    xr = fread(data_id,xsize,'*float',recsize);
    % now read again but the imaginary part
    fseek(data_id,offset+recsize,'bof'); 
    xi = fread(data_id,xsize,'*float',recsize);
    % combine the vectors
    xc = complex(xr,xi);
    
    % now find out where it goes to put it back
    for j=1:N,
        dim_index(j) = petable(i,dims(j)) + abs(dim_min(j)) + 1;
    end
 
    % fix for multi coils
    if multiple_coils
        dim_index(N) = petable(i,dims(N));
    end
    
    % Can't think of a more succinct method for this switch?
    switch N
        case 1
            kspace(:,dim_index(1)) = xc;
        case 2
            kspace(:,dim_index(1),dim_index(2)) = xc;
        case 3
            kspace(:,dim_index(1),dim_index(2),dim_index(3)) = xc;
        case 4
            kspace(:,dim_index(1),dim_index(2),dim_index(3),dim_index(4)) = xc;
        %suribe
        case 5
            kspace(:,dim_index(1),dim_index(2),dim_index(3),dim_index(4),dim_index(5)) = xc;
        otherwise
            disp('Error with matrix size')
            return
    end
    
end

fclose(data_id);
if(~NoConsoleOutput)
    disp(' ........ done');
end
end

function [kspace, petable] = raw_READ_STD(filename,type)

headername = [filename '.list'];
dataname = [filename '.data'];

% first look in the list file to find the number of dynamics, dimensions
% etc. do this by looking for lines starting with a '.' At present this
% will not work for the kx,ky,kz sizes as they have 2 entries and all
% others have one. Make a cell array to hold the variables. There are 8
% main variables (not including number of mixes), and if there are 
% multiple coils then there will be n_loc more lines giving coil combinations per slice

list_id = fopen(headername,'r','l');
if (list_id==-1)
    disp('Error: file not found or something');
    return
end

line_counter=1;
while (line_counter<=9)
    line = fgetl(list_id);
    if (strcmp(line(1),'.'))
        s = sscanf(line,'%*s %*s %*s %*s %s %s %s');
        % get the string with colon and use this to divide it
        div = strfind(s,':');
        % get the string above and store in cell array
        attributes{line_counter,1} = s(1:div-1);
        attributes{line_counter,2} = str2num(s(div+1:length(s)));
        line_counter = line_counter + 1;
    end
end

% attributes row num.   variable
%       1               mixes
%       2               spatial dimensions
%       3               dynamics
%       4               cardiac phases
%       5               echoes
%       6               locations (slices)
%       7               extra 1
%       8               extra 2
%       9               averages

% Now see if there are any extra coils .. if there are then the next line
% starting with '.' should be for receiver coils
done = false;
while (~done)
    line = fgetl(list_id);
    if (strcmp(line(1),'.'))
        % attempt to find the string 'coil_channel_combination'
        found = strfind(line,'coil_channel_combination');
        if (length(found)>0)
            multiple_coils = true;
        else
            multiple_coils = false;
        end
        done = true;
    end
end

% If there are multiple coils, which ones? Assume that the coils used are
% the same for each slice.
if multiple_coils
    div = strfind(line,':');
    % assume there are 10 coils (0:9) and there is a 1 if they are used and
    % a zero otherwise
    coils = sscanf(line((div+1):length(line)),'%i %i %i %i %i %i %i %i %i %i');
end

% There isn't really any need right now for the cell array... Store the
% variables in the order above in a numeric array
for i=1:9,
    data_dimensions(i) = cell2mat(attributes(i,2));
end

% complain if number of mixes is not 1
if (data_dimensions(1) ~= 1)
    disp('Error: number of mixes is greater than one, case not handled by code');
    return
end

% now read the .list file to get the phase encode table.
done=false;
rowcount=1;
while (~done)
    line = fgetl(list_id);
    if(line == -1)
        done=true;
        continue;
    end

    if ((line(1)=='#')||(line(1)=='.'))
        continue
    end

    % now we want to loop through lines beginning STD....
    linelength = numel(line);
    tag = strtok(line);
    % to read in NOI or PHX for example replace following with a case switch
    if (tag == 'STD')
        row = sscanf(line(6:linelength)','%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f' );
        petable(rowcount,:) = row;
        rowcount = rowcount+1;
    end
end
fclose(list_id);

% the 19 columns that are read in are:
% mix   dyn   card  echo  loca  chan  extr1 extr2 ky    kz    aver  sign
% rf    grad  enc   rtop  rr    size   offset
% Not sure how to deal with all possible eventualities. For now I'll
% try to do it generally. Matrix dimension can't be arbitrarily high so
% this is a constraint.

% don't want to cause more than a 5D matrix (inc x, 4D+coils)
% N = number of dimensions to read not including x (inc coils)
N = numel(find(data_dimensions>1));
% account for 3D case
if (data_dimensions(2)>2)
    N=N+1;
end

if multiple_coils
    disp('Found multiple channel data')
    N = N+1;
end


if (N>5)
    disp('Error too many dimensions, case not handled by code');
    return
else
    disp(['Data matrix will have ' num2str(N+1) ' dims']);
end


% store the contributing dimension column numbers in a vector. 
% this is a bit inelegant but never mind
nonzero = find(data_dimensions>1)
counter =1;
disp('dimension 1 = kx');
for i=1:length(nonzero),
    switch nonzero(i)
        case 1
            disp('Error number of mixes greater than 1');
        case 2
            % this should always be one
            dims(counter)=9; %ky entry
            disp(['dimension ',num2str(counter+1),' = ky']);
            counter = counter+1;
            if (data_dimensions(2)>2)
                %have 2 PE dirs, so also read kz
                dims(counter)=10;
                disp(['dimension ',num2str(counter+1),' = kz']);
                counter=counter+1;
            end
        case 3
            %dynamic data
             dims(counter)=2;
             disp(['dimension ',num2str(counter+1),' = dynamic']);
             counter = counter +1;
        case 4
            %multiple cardiac phases
            dims(counter)=3;
            disp(['dimension ',num2str(counter+1),' = cardiac phase']);
            counter=counter+1;
        case 5
            %multiple echoes
            dims(counter)=4;
            disp(['dimension ',num2str(counter+1),' = echo number']);
            counter=counter+1;
        case 6
            %multi slice data
            dims(counter)=5;
            disp(['dimension ',num2str(counter+1),' = slice']);
            counter = counter+1;
        case 7
            %extra 1
            dims(counter)=7;
            disp(['dimension ',num2str(counter+1),' = extra 1']);
            counter=counter+1;
        case 8
            %extra 2
            dims(counter)=8;
            disp(['dimension ',num2str(counter+1),' = extra 2']);;
            counter=counter+1;
        case 9
            %multi averages
            dims(counter)=11;
            disp(['dimension ',num2str(counter+1),' = averages']);
            counter=counter+1;
    end
end
if multiple_coils
    dims(N)=6;
    disp(['dimension ',num2str(counter+1),' = coil channels']);
end

if (length(dims)>N)
    % must be a problem
    display('Error with dimensions ....');
    return
end
    

% Have N dimensions plus x to worry about.
recsize=4;
xsize = petable(1,18)/(2*recsize);
% hopefully all kx-vectors are the same length
for i=1:N,
    dim_min(i) = min(petable(:,dims(i)));
    dim_max(i) = max(petable(:,dims(i)));
    dim_size(i)= dim_max(i)+abs(dim_min(i))+1;
end

% the above estimation doesn't work for multiple coil data so correct it
if multiple_coils
    dim_size(N) = dim_max(N);
    % doing it this way bypasses the need to use the info taken from the
    % header file
end

kspace = complex(zeros([xsize dim_size],'single'));


% now open the data file and start reading
data_id = fopen(dataname,'r','l');
num_entries = numel(petable(:,1));

disp('Beginning binary read ... ');
for i=1:num_entries,

    offset = petable(i,19);
    % read in the x vector
    fseek(data_id,offset,'bof');
    xr = fread(data_id,xsize,'*float',recsize);
    % now read again but the imaginary part
    fseek(data_id,offset+recsize,'bof'); 
    xi = fread(data_id,xsize,'*float',recsize);
    % combine the vectors
    xc = complex(xr,xi);
    
    % now find out where it goes to put it back
    for j=1:N,
        dim_index(j) = petable(i,dims(j)) + abs(dim_min(j)) + 1;
    end
 
    % fix for multi coils
    if multiple_coils
        dim_index(N) = petable(i,dims(N));
    end
    
    % Can't think of a more succinct method for this switch?
    switch N
        case 1
            kspace(:,dim_index(1)) = xc;
        case 2
            kspace(:,dim_index(1),dim_index(2)) = xc;
        case 3
            kspace(:,dim_index(1),dim_index(2),dim_index(3)) = xc;
        case 4
            kspace(:,dim_index(1),dim_index(2),dim_index(3),dim_index(4)) = xc;
        %suribe
        case 5
            kspace(:,dim_index(1),dim_index(2),dim_index(3),dim_index(4),dim_index(5)) = xc;
        otherwise
            disp('Error with matrix size')
            return
    end
    
end

fclose(data_id);
disp(' ........ done');
end
 
function [kspace, petable] = raw_READ_FRX(filename,type)

headername = [filename '.list'];
dataname = [filename '.data'];

% first look in the list file to find the number of dynamics, dimensions
% etc. do this by looking for lines starting with a '.' At present this
% will not work for the kx,ky,kz sizes as they have 2 entries and all
% others have one. Make a cell array to hold the variables. There are 8
% main variables (not including number of mixes), and if there are 
% multiple coils then there will be n_loc more lines giving coil combinations per slice

list_id = fopen(headername,'r','l');
if (list_id==-1)
    disp('Error: file not found or something');
    return
end

line_counter=1;
while (line_counter<=9)
    line = fgetl(list_id);
    if (strcmp(line(1),'.'))
        s = sscanf(line,'%*s %*s %*s %*s %s %s %s');
        % get the string with colon and use this to divide it
        div = strfind(s,':');
        % get the string above and store in cell array
        attributes{line_counter,1} = s(1:div-1);
        attributes{line_counter,2} = str2num(s(div+1:length(s)));
        line_counter = line_counter + 1;
    end
end

% attributes row num.   variable
%       1               mixes
%       2               spatial dimensions
%       3               dynamics
%       4               cardiac phases
%       5               echoes
%       6               locations (slices)
%       7               extra 1
%       8               extra 2
%       9               averages

% Now see if there are any extra coils .. if there are then the next line
% starting with '.' should be for receiver coils
done = false;
while (~done)
    line = fgetl(list_id);
    if (strcmp(line(1),'.'))
        % attempt to find the string 'coil_channel_combination'
        found = strfind(line,'coil_channel_combination');
        if (length(found)>0)
            multiple_coils = true;
        else
            multiple_coils = false;
        end
        done = true;
    end
end

% If there are multiple coils, which ones? Assume that the coils used are
% the same for each slice.
if multiple_coils
    div = strfind(line,':');
    % assume there are 10 coils (0:9) and there is a 1 if they are used and
    % a zero otherwise
    coils = sscanf(line((div+1):length(line)),'%i %i %i %i %i %i %i %i %i %i');
end

% There isn't really any need right now for the cell array... Store the
% variables in the order above in a numeric array
for i=1:9,
    data_dimensions(i) = cell2mat(attributes(i,2));
end

% complain if number of mixes is not 1
if (data_dimensions(1) ~= 1)
    disp('Error: number of mixes is greater than one, case not handled by code');
    return
end

% now read the .list file to get the phase encode table.
done=false;
rowcount=1;
while (~done)
    line = fgetl(list_id);
    if(line == -1)
        done=true;
        continue;
    end

    if ((line(1)=='#')||(line(1)=='.'))
        continue
    end

    % now we want to loop through lines beginning STD....
    linelength = numel(line);
    tag = strtok(line);
    % to read in NOI or PHX for example replace following with a case switch
    if (tag == 'FRX')
        row = sscanf(line(6:linelength)','%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f' );
        petable(rowcount,:) = row;
        rowcount = rowcount+1;
    end
end
fclose(list_id);

% the 19 columns that are read in are:
% mix   dyn   card  echo  loca  chan  extr1 extr2 ky    kz    aver  sign
% rf    grad  enc   rtop  rr    size   offset
% Not sure how to deal with all possible eventualities. For now I'll
% try to do it generally. Matrix dimension can't be arbitrarily high so
% this is a constraint.

% don't want to cause more than a 5D matrix (inc x, 4D+coils)
% N = number of dimensions to read not including x (inc coils)
N = numel(find(data_dimensions>1));
% account for 3D case
if (data_dimensions(2)>2)
    N=N+1;
end

if multiple_coils
    disp('Found multiple channel data')
    N = N+1;
end


if (N>5)
    disp('Error too many dimensions, case not handled by code');
    return
else
    disp(['Data matrix will have ' num2str(N+1) ' dims']);
end


% store the contributing dimension column numbers in a vector. 
% this is a bit inelegant but never mind
nonzero = find(data_dimensions>1)
counter =1;
disp('dimension 1 = kx');
for i=1:length(nonzero),
    switch nonzero(i)
        case 1
            disp('Error number of mixes greater than 1');
        case 2
            % this should always be one
            dims(counter)=9; %ky entry
            disp(['dimension ',num2str(counter+1),' = ky']);
            counter = counter+1;
            if (data_dimensions(2)>2)
                %have 2 PE dirs, so also read kz
                dims(counter)=10;
                disp(['dimension ',num2str(counter+1),' = kz']);
                counter=counter+1;
            end
        case 3
            %dynamic data
             dims(counter)=2;
             disp(['dimension ',num2str(counter+1),' = dynamic']);
             counter = counter +1;
        case 4
            %multiple cardiac phases
            dims(counter)=3;
            disp(['dimension ',num2str(counter+1),' = cardiac phase']);
            counter=counter+1;
        case 5
            %multiple echoes
            dims(counter)=4;
            disp(['dimension ',num2str(counter+1),' = echo number']);
            counter=counter+1;
        case 6
            %multi slice data
            dims(counter)=5;
            disp(['dimension ',num2str(counter+1),' = slice']);
            counter = counter+1;
        case 7
            %extra 1
            dims(counter)=7;
            disp(['dimension ',num2str(counter+1),' = extra 1']);
            counter=counter+1;
        case 8
            %extra 2
            dims(counter)=8;
            disp(['dimension ',num2str(counter+1),' = extra 2']);;
            counter=counter+1;
        case 9
            %multi averages
            dims(counter)=11;
            disp(['dimension ',num2str(counter+1),' = averages']);
            counter=counter+1;
    end
end
if multiple_coils
    dims(N)=6;
    disp(['dimension ',num2str(counter+1),' = coil channels']);
end

if (length(dims)>N)
    % must be a problem
    display('Error with dimensions ....');
    return
end
    

% Have N dimensions plus x to worry about.
recsize=4;
xsize = petable(1,18)/(2*recsize);
% hopefully all kx-vectors are the same length
for i=1:N,
    dim_min(i) = min(petable(:,dims(i)));
    dim_max(i) = max(petable(:,dims(i)));
    dim_size(i)= dim_max(i)+abs(dim_min(i))+1;
end

% the above estimation doesn't work for multiple coil data so correct it
if multiple_coils
    dim_size(N) = dim_max(N);
    % doing it this way bypasses the need to use the info taken from the
    % header file
end

kspace = complex(zeros([xsize dim_size],'single'));


% now open the data file and start reading
data_id = fopen(dataname,'r','l');
num_entries = numel(petable(:,1));

disp('Beginning binary read ... ');
for i=1:num_entries,

    offset = petable(i,19);
    % read in the x vector
    fseek(data_id,offset,'bof');
    xr = fread(data_id,xsize,'*float',recsize);
    % now read again but the imaginary part
    fseek(data_id,offset+recsize,'bof'); 
    xi = fread(data_id,xsize,'*float',recsize);
    % combine the vectors
    xc = complex(xr,xi);
    
    % now find out where it goes to put it back
    for j=1:N,
        dim_index(j) = petable(i,dims(j)) + abs(dim_min(j)) + 1;
    end
 
    % fix for multi coils
    if multiple_coils
        dim_index(N) = petable(i,dims(N));
    end
    
    % Can't think of a more succinct method for this switch?
    switch N
        case 1
            kspace(:,dim_index(1)) = xc;
        case 2
            kspace(:,dim_index(1),dim_index(2)) = xc;
        case 3
            kspace(:,dim_index(1),dim_index(2),dim_index(3)) = xc;
        case 4
            kspace(:,dim_index(1),dim_index(2),dim_index(3),dim_index(4)) = xc;
        %suribe
        case 5
            kspace(:,dim_index(1),dim_index(2),dim_index(3),dim_index(4),dim_index(5)) = xc;
        otherwise
            disp('Error with matrix size')
            return
    end
    
end

fclose(data_id);
disp(' ........ done');
end

function [kspace, petable] = raw_READ_PHX(filename,type)

headername = [filename '.list'];
dataname = [filename '.data'];

% first look in the list file to find the number of dynamics, dimensions
% etc. do this by looking for lines starting with a '.' At present this
% will not work for the kx,ky,kz sizes as they have 2 entries and all
% others have one. Make a cell array to hold the variables. There are 8
% main variables (not including number of mixes), and if there are 
% multiple coils then there will be n_loc more lines giving coil combinations per slice

list_id = fopen(headername,'r','l');
if (list_id==-1)
    disp('Error: file not found or something');
    return
end

line_counter=1;
while (line_counter<=9)
    line = fgetl(list_id);
    if (strcmp(line(1),'.'))
        s = sscanf(line,'%*s %*s %*s %*s %s %s %s');
        % get the string with colon and use this to divide it
        div = strfind(s,':');
        % get the string above and store in cell array
        attributes{line_counter,1} = s(1:div-1);
        attributes{line_counter,2} = str2num(s(div+1:length(s)));
        line_counter = line_counter + 1;
    end
end

% attributes row num.   variable
%       1               mixes
%       2               spatial dimensions
%       3               dynamics
%       4               cardiac phases
%       5               echoes
%       6               locations (slices)
%       7               extra 1
%       8               extra 2
%       9               averages

% Now see if there are any extra coils .. if there are then the next line
% starting with '.' should be for receiver coils
done = false;
while (~done)
    line = fgetl(list_id);
    if (strcmp(line(1),'.'))
        % attempt to find the string 'coil_channel_combination'
        found = strfind(line,'coil_channel_combination');
        if (length(found)>0)
            multiple_coils = true;
        else
            multiple_coils = false;
        end
        done = true;
    end
end

% If there are multiple coils, which ones? Assume that the coils used are
% the same for each slice.
if multiple_coils
    div = strfind(line,':');
    % assume there are 10 coils (0:9) and there is a 1 if they are used and
    % a zero otherwise
    coils = sscanf(line((div+1):length(line)),'%i %i %i %i %i %i %i %i %i %i');
end

% There isn't really any need right now for the cell array... Store the
% variables in the order above in a numeric array
for i=1:9,
    data_dimensions(i) = cell2mat(attributes(i,2));
end

% complain if number of mixes is not 1
if (data_dimensions(1) ~= 1)
    disp('Error: number of mixes is greater than one, case not handled by code');
    return
end

% now read the .list file to get the phase encode table.
done=false;
rowcount=1;
while (~done)
    line = fgetl(list_id);
    if(line == -1)
        done=true;
        continue;
    end

    if ((line(1)=='#')||(line(1)=='.'))
        continue
    end

    % now we want to loop through lines beginning STD....
    linelength = numel(line);
    tag = strtok(line);
    % to read in NOI or PHX for example replace following with a case switch
    if (tag == 'PHX')
        row = sscanf(line(6:linelength)','%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f' );
        petable(rowcount,:) = row;
        rowcount = rowcount+1;
    end
end
fclose(list_id);

% the 19 columns that are read in are:
% mix   dyn   card  echo  loca  chan  extr1 extr2 ky    kz    aver  sign
% rf    grad  enc   rtop  rr    size   offset
% Not sure how to deal with all possible eventualities. For now I'll
% try to do it generally. Matrix dimension can't be arbitrarily high so
% this is a constraint.

% don't want to cause more than a 5D matrix (inc x, 4D+coils)
% N = number of dimensions to read not including x (inc coils)
N = numel(find(data_dimensions>1));
% account for 3D case
if (data_dimensions(2)>2)
    N=N+1;
end

if multiple_coils
    disp('Found multiple channel data')
    N = N+1;
end


if (N>5)
    disp('Error too many dimensions, case not handled by code');
    return
else
    disp(['Data matrix will have ' num2str(N+1) ' dims']);
end


% store the contributing dimension column numbers in a vector. 
% this is a bit inelegant but never mind
nonzero = find(data_dimensions>1)
counter =1;
disp('dimension 1 = kx');
for i=1:length(nonzero),
    switch nonzero(i)
        case 1
            disp('Error number of mixes greater than 1');
        case 2
            % this should always be one
            dims(counter)=9; %ky entry
            disp(['dimension ',num2str(counter+1),' = ky']);
            counter = counter+1;
            if (data_dimensions(2)>2)
                %have 2 PE dirs, so also read kz
                dims(counter)=10;
                disp(['dimension ',num2str(counter+1),' = kz']);
                counter=counter+1;
            end
        case 3
            %dynamic data
             dims(counter)=2;
             disp(['dimension ',num2str(counter+1),' = dynamic']);
             counter = counter +1;
        case 4
            %multiple cardiac phases
            dims(counter)=3;
            disp(['dimension ',num2str(counter+1),' = cardiac phase']);
            counter=counter+1;
        case 5
            %multiple echoes
            dims(counter)=4;
            disp(['dimension ',num2str(counter+1),' = echo number']);
            counter=counter+1;
        case 6
            %multi slice data
            dims(counter)=5;
            disp(['dimension ',num2str(counter+1),' = slice']);
            counter = counter+1;
        case 7
            %extra 1
            dims(counter)=7;
            disp(['dimension ',num2str(counter+1),' = extra 1']);
            counter=counter+1;
        case 8
            %extra 2
            dims(counter)=8;
            disp(['dimension ',num2str(counter+1),' = extra 2']);;
            counter=counter+1;
        case 9
            %multi averages
            dims(counter)=11;
            disp(['dimension ',num2str(counter+1),' = averages']);
            counter=counter+1;
    end
end
if multiple_coils
    dims(N)=6;
    disp(['dimension ',num2str(counter+1),' = coil channels']);
end

if (length(dims)>N)
    % must be a problem
    display('Error with dimensions ....');
    return
end
    

% Have N dimensions plus x to worry about.
recsize=4;
xsize = petable(1,18)/(2*recsize);
% hopefully all kx-vectors are the same length
for i=1:N,
    dim_min(i) = min(petable(:,dims(i)));
    dim_max(i) = max(petable(:,dims(i)));
    dim_size(i)= dim_max(i)+abs(dim_min(i))+1;
end

% the above estimation doesn't work for multiple coil data so correct it
if multiple_coils
    dim_size(N) = dim_max(N);
    % doing it this way bypasses the need to use the info taken from the
    % header file
end

kspace = complex(zeros([xsize dim_size],'single'));


% now open the data file and start reading
data_id = fopen(dataname,'r','l');
num_entries = numel(petable(:,1));

disp('Beginning binary read ... ');
for i=1:num_entries,

    offset = petable(i,19);
    % read in the x vector
    fseek(data_id,offset,'bof');
    xr = fread(data_id,xsize,'*float',recsize);
    % now read again but the imaginary part
    fseek(data_id,offset+recsize,'bof'); 
    xi = fread(data_id,xsize,'*float',recsize);
    % combine the vectors
    xc = complex(xr,xi);
    
    % now find out where it goes to put it back
    for j=1:N,
        dim_index(j) = petable(i,dims(j)) + abs(dim_min(j)) + 1;
    end
 
    % fix for multi coils
    if multiple_coils
        dim_index(N) = petable(i,dims(N));
    end
    
    % Can't think of a more succinct method for this switch?
    switch N
        case 1
            kspace(:,dim_index(1)) = xc;
        case 2
            kspace(:,dim_index(1),dim_index(2)) = xc;
        case 3
            kspace(:,dim_index(1),dim_index(2),dim_index(3)) = xc;
        case 4
            kspace(:,dim_index(1),dim_index(2),dim_index(3),dim_index(4)) = xc;
        %suribe
        case 5
            kspace(:,dim_index(1),dim_index(2),dim_index(3),dim_index(4),dim_index(5)) = xc;
        otherwise
            disp('Error with matrix size')
            return
    end
    
end

fclose(data_id);
disp(' ........ done');
end

function [kspace, petable] = raw_READ_NOI(filename,type)

headername = [filename '.list'];
dataname = [filename '.data'];

% first look in the list file to find the number of dynamics, dimensions
% etc. do this by looking for lines starting with a '.' At present this
% will not work for the kx,ky,kz sizes as they have 2 entries and all
% others have one. Make a cell array to hold the variables. There are 8
% main variables (not including number of mixes), and if there are 
% multiple coils then there will be n_loc more lines giving coil combinations per slice

list_id = fopen(headername,'r','l');
if (list_id==-1)
    disp('Error: file not found or something');
    return
end

line_counter=1;
while (line_counter<=9)
    line = fgetl(list_id);
    if (strcmp(line(1),'.'))
        s = sscanf(line,'%*s %*s %*s %*s %s %s %s');
        % get the string with colon and use this to divide it
        div = strfind(s,':');
        % get the string above and store in cell array
        attributes{line_counter,1} = s(1:div-1);
        attributes{line_counter,2} = str2num(s(div+1:length(s)));
        line_counter = line_counter + 1;
    end
end

% attributes row num.   variable
%       1               mixes
%       2               spatial dimensions
%       3               dynamics
%       4               cardiac phases
%       5               echoes
%       6               locations (slices)
%       7               extra 1
%       8               extra 2
%       9               averages

% Now see if there are any extra coils .. if there are then the next line
% starting with '.' should be for receiver coils
done = false;
while (~done)
    line = fgetl(list_id);
    if (strcmp(line(1),'.'))
        % attempt to find the string 'coil_channel_combination'
        found = strfind(line,'coil_channel_combination');
        if (length(found)>0)
            multiple_coils = true;
        else
            multiple_coils = false;
        end
        done = true;
    end
end

% If there are multiple coils, which ones? Assume that the coils used are
% the same for each slice.
if multiple_coils
    div = strfind(line,':');
    % assume there are 10 coils (0:9) and there is a 1 if they are used and
    % a zero otherwise
    coils = sscanf(line((div+1):length(line)),'%i %i %i %i %i %i %i %i %i %i');
end

% There isn't really any need right now for the cell array... Store the
% variables in the order above in a numeric array
for i=1:9,
    data_dimensions(i) = cell2mat(attributes(i,2));
end

% complain if number of mixes is not 1
if (data_dimensions(1) ~= 1)
    disp('Error: number of mixes is greater than one, case not handled by code');
    return
end

% now read the .list file to get the phase encode table.
done=false;
rowcount=1;
while (~done)
    line = fgetl(list_id);
    if(line == -1)
        done=true;
        continue;
    end

    if ((line(1)=='#')||(line(1)=='.'))
        continue
    end

    % now we want to loop through lines beginning STD....
    linelength = numel(line);
    tag = strtok(line);
    % to read in NOI or PHX for example replace following with a case switch
    if (tag == 'NOI')
        row = sscanf(line(6:linelength)','%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f' );
        petable(rowcount,:) = row;
        rowcount = rowcount+1;
    end
end
fclose(list_id);

% the 19 columns that are read in are:
% mix   dyn   card  echo  loca  chan  extr1 extr2 ky    kz    aver  sign
% rf    grad  enc   rtop  rr    size   offset
% Not sure how to deal with all possible eventualities. For now I'll
% try to do it generally. Matrix dimension can't be arbitrarily high so
% this is a constraint.

% don't want to cause more than a 5D matrix (inc x, 4D+coils)
% N = number of dimensions to read not including x (inc coils)
N = numel(find(data_dimensions>1));
% account for 3D case
if (data_dimensions(2)>2)
    N=N+1;
end

if multiple_coils
    disp('Found multiple channel data')
    N = N+1;
end


if (N>5)
    disp('Error too many dimensions, case not handled by code');
    return
else
    disp(['Data matrix will have ' num2str(N+1) ' dims']);
end


% store the contributing dimension column numbers in a vector. 
% this is a bit inelegant but never mind
nonzero = find(data_dimensions>1)
counter =1;
disp('dimension 1 = kx');
for i=1:length(nonzero),
    switch nonzero(i)
        case 1
            disp('Error number of mixes greater than 1');
        case 2
            % this should always be one
            dims(counter)=9; %ky entry
            disp(['dimension ',num2str(counter+1),' = ky']);
            counter = counter+1;
            if (data_dimensions(2)>2)
                %have 2 PE dirs, so also read kz
                dims(counter)=10;
                disp(['dimension ',num2str(counter+1),' = kz']);
                counter=counter+1;
            end
        case 3
            %dynamic data
             dims(counter)=2;
             disp(['dimension ',num2str(counter+1),' = dynamic']);
             counter = counter +1;
        case 4
            %multiple cardiac phases
            dims(counter)=3;
            disp(['dimension ',num2str(counter+1),' = cardiac phase']);
            counter=counter+1;
        case 5
            %multiple echoes
            dims(counter)=4;
            disp(['dimension ',num2str(counter+1),' = echo number']);
            counter=counter+1;
        case 6
            %multi slice data
            dims(counter)=5;
            disp(['dimension ',num2str(counter+1),' = slice']);
            counter = counter+1;
        case 7
            %extra 1
            dims(counter)=7;
            disp(['dimension ',num2str(counter+1),' = extra 1']);
            counter=counter+1;
        case 8
            %extra 2
            dims(counter)=8;
            disp(['dimension ',num2str(counter+1),' = extra 2']);;
            counter=counter+1;
        case 9
            %multi averages
            dims(counter)=11;
            disp(['dimension ',num2str(counter+1),' = averages']);
            counter=counter+1;
    end
end
if multiple_coils
    dims(N)=6;
    disp(['dimension ',num2str(counter+1),' = coil channels']);
end

if (length(dims)>N)
    % must be a problem
    display('Error with dimensions ....');
    return
end
    

% Have N dimensions plus x to worry about.
recsize=4;
xsize = petable(1,18)/(2*recsize);
% hopefully all kx-vectors are the same length
for i=1:N,
    dim_min(i) = min(petable(:,dims(i)));
    dim_max(i) = max(petable(:,dims(i)));
    dim_size(i)= dim_max(i)+abs(dim_min(i))+1;
end

% the above estimation doesn't work for multiple coil data so correct it
if multiple_coils
    dim_size(N) = dim_max(N);
    % doing it this way bypasses the need to use the info taken from the
    % header file
end

kspace = complex(zeros([xsize dim_size],'single'));


% now open the data file and start reading
data_id = fopen(dataname,'r','l');
num_entries = numel(petable(:,1));

disp('Beginning binary read ... ');
for i=1:num_entries,

    offset = petable(i,19);
    % read in the x vector
    fseek(data_id,offset,'bof');
    xr = fread(data_id,xsize,'*float',recsize);
    % now read again but the imaginary part
    fseek(data_id,offset+recsize,'bof'); 
    xi = fread(data_id,xsize,'*float',recsize);
    % combine the vectors
    xc = complex(xr,xi);
    
    % now find out where it goes to put it back
    for j=1:N,
        dim_index(j) = petable(i,dims(j)) + abs(dim_min(j)) + 1;
    end
 
    % fix for multi coils
    if multiple_coils
        dim_index(N) = petable(i,dims(N));
    end
    
    % Can't think of a more succinct method for this switch?
    switch N
        case 1
            kspace(:,dim_index(1)) = xc;
        case 2
            kspace(:,dim_index(1),dim_index(2)) = xc;
        case 3
            kspace(:,dim_index(1),dim_index(2),dim_index(3)) = xc;
        case 4
            kspace(:,dim_index(1),dim_index(2),dim_index(3),dim_index(4)) = xc;
        %suribe
        case 5
            kspace(:,dim_index(1),dim_index(2),dim_index(3),dim_index(4),dim_index(5)) = xc;
        otherwise
            disp('Error with matrix size')
            return
    end
    
end

fclose(data_id);
disp(' ........ done');
end
