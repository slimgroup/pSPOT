%% Create a memory map to a file
    m=memmapfile('/users/slic/troberso/jSPOT/random_nums_file.mat', 'format', 'double')

%% Set writable to true
    m.writable=true;

%% Read the first value of the file, and change it to a random number on the file.
    m.Data(1)
    m.Data(1)=(rand(1)*100);
    m.Data(1)

%% Apply a SPOT operator to the data
    path (path, [pwd, '/ConocoPhillipsCode/spotbox-1.0p']);
    m.Data(1)
    D=opEye(size(m.Data,1))*2;
    m.Data = D*m.Data;
    m.Data(1)
    
    
    
    % memmap does not support complex numbers, so opDFT doesn't work
    % D = opDFT(size(m.Data,1), 1);
    % m.Data = complex((D*m.Data));

%% JavaSeis memmap
m=memmapfile('/users/slic/troberso/demo_CP_code/test.js/TraceFile', 'format', 'double')
m.Data(1)
%The format type is not correct for a .js file. 