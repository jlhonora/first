function rawinfo=raw_read_info(filename,table)
% Collects information from LIST/RAW files:
% PHX lines
% FRX lines
% x and y resolution
% kx and ky range
% kx and ky oversample factor
% x and y range

% MV - 08-10-2007

% gets file names
if exist('filename')
    headername = [filename '.list'];   
    dataname = [filename '.data'];      
else
    [file, path] = uigetfile('*.list');
    headername = [path file];
    [path, filename, ext] = fileparts(headername);
    dataname = [path filesep filename '.data'];
end

% opens .list and .data files
lis = fopen(headername,'r','l');
if (lis == -1)
    error('Problem opening .list file');
end
dat = fopen(dataname,'r','l');
if (dat == -1)
    error('Problem opening .data file');
end

% reads info, PHX and FRX lines from list file
i=1;
szline=4;
while i<=szline
    linha=fgetl(lis);
    if ~isempty(strfind(linha,'X-resolution'))
       xres=sscanf(linha,'%*55c %f',1);
    end
    if ~isempty(strfind(linha,'Y-resolution'))
       yres=sscanf(linha,'%*55c %f',1);
    end  
    if ~isempty(strfind(linha,'kx_range'))
       kx_range=sscanf(linha,'%*55c %f %f',[1,2]);
    end
    if ~isempty(strfind(linha,'ky_range'))
       ky_range=sscanf(linha,'%*55c %f %f',[1,2]);
    end
    if ~isempty(strfind(linha,'X_range'))
       x_range=sscanf(linha,'%*55c %f %f',[1,2]);
    end
    if ~isempty(strfind(linha,'Y_range'))
       y_range=sscanf(linha,'%*55c %f %f',[1,2]);
    end
    if ~isempty(strfind(linha,'kx_oversample_factor'))
       kx_over=sscanf(linha,'%*55c %f',1);
    end
    if ~isempty(strfind(linha,'ky_oversample_factor'))
       ky_over=sscanf(linha,'%*55c %f',1);
    end
    if  ((length(linha)>5)&(linha(1)~='#')&(linha(3:5)=='PHX'))
        line(i,:)=linha;
        i=i+1;
    end
    if  ((length(linha)>5)&(linha(1)~='#')&(linha(3:5)=='FRX'))
        line(i,:)=linha;
        i=i+1;
    end
end

% gets the size, starting position and signal of the PHX and FRX data
for i=1:szline
    info(i,:)=sscanf(line(i,:),'%*110c %f %f',[1,2]); % NNU
  % info(i,:)=sscanf(line(i,:),'%*115c %f %f',[1,2]); % 3T
end
offs=info(:,end);
sz=info(:,end-1)/4; % because each info has 4 bits
     
% opens data file in the selected position
stat=fseek(dat,offs(1),'bof'); 
if stat==-1
    error('Problem reading .data file');
end

% reads the FRX lines

echoes = max(table(:,4))+1
for ec=1:echoes
    for j=1:2
        freqdata(j,:,ec)=fread(dat,sz(j),'*float32'); 
        for i=1:sz(j)/2
            refreq(j,i,ec)=freqdata(j,2*i-1,ec);
            imfreq(j,i,ec)=freqdata(j,2*i,ec);
        end
        freq(j,:,ec)=complex(refreq(j,:,ec),imfreq(j,:,ec));
    end
end

% reads the PHX lines 
for ec=1:echoes
    for j=1:2 
        phxdata(j,:,ec)=fread(dat,sz(j+2),'*float32'); 
        for i=1:sz(j)/2
            rephx(j,i,ec)=phxdata(j,2*i-1,ec);
            imphx(j,i,ec)=phxdata(j,2*i,ec);
        end
        phx(j,:,ec)=complex(rephx(j,:,ec),imphx(j,:,ec));
    end
end
fclose('all');

% Constructs the output
rawinfo.phx=phx;
rawinfo.freq=freq;
rawinfo.xres=xres
rawinfo.yres=yres
rawinfo.kx_range=kx_range;
rawinfo.ky_range=ky_range;
rawinfo.x_range=x_range;
rawinfo.y_range=y_range;

