% Daniel Wagenaar, Anat Kahan , Cell Reports 2021
function x = load2p(fn)
% LOAD2P - Load data from 2-photon microscope
%   x = LOAD2P(fn) loads the named file, which must have been produced
%   by the Gradinaru Lab 2-photon LabView system ("GLab2P").
%   The result is a structure with several fields:
%      dat - YxXxT PMT data
%      tt - Tx1 timestamps for each frame (in ms since the first frame)
%      scope - scope metainformation from the file
%      frame - frame metainformation from the file
%      seq - sequence metainformation from the file (or [] if n/a)
%   x = LOAD2P without a filename argument prompts the user
%   LOAD2P(fn) or LOAD2P saves the results in a tiff file, based on
%   user input.

if nargin==0
    [fn,pn] = uigetfile('*.dat', 'Select a Glab-2P data file to load', ...
        'C:\Users\Gradinaru Lab\Desktop\Two Photon Data');
    if isempty(fn)
        disp('No selection made - canceling operation');
        return
    end
    fn = [pn filesep fn];
end

fd = fopen(fn, 'rb');
if fd<0
    error('Could not open the named file');
end

str = fread(fd, [1 4], '*char');
if strcmp(str, 'G2P1')
    x = load2p_singleframe(fd);
    fclose(fd);
elseif strcmp(str, 'G2P2')
    x = load2p_sequence(fd);
    fclose(fd);
else
    fclose(fd);
    error('Not a GLab-2P file');
end

if nargout==0
    load2p_saveastiff(x);
end

function x = load2p_sequence(fd)
N = fread(fd, [1 1], 'int32');
str = fread(fd, [1 N], '*char');
x.seq = jsondecode(str);
N = fread(fd, [1 1], 'int32');
str = fread(fd, [1 N], '*char');
x.scope = jsondecode(str);
N = fread(fd, [1 1], 'int32');
str = fread(fd, [1 N], '*char');
x.frame = jsondecode(str);
x.tt = zeros(0,1);
x.dat = [];
while 1
    S = fread(fd, [1 2], 'int32');
    if isempty(S)
        break;
    end
    frm = fread(fd, fliplr(S), '*uint16')';
    if isempty(x.dat)
        x.dat = frm;
    else
        x.dat(:,:,end+1) = frm;
    end
    t = fread(fd, [1 1], 'int32');
    x.tt(end+1) = t;
end
x.tt = x.tt(:) - x.seq.StartTime;


function x = load2p_singleframe(fd)
N = fread(fd, [1 1], 'int32');
str = fread(fd, [1 N], '*char');
x.scope = jsondecode(str);
N = fread(fd, [1 1], 'int32');
str = fread(fd, [1 N], '*char');
x.frame = jsondecode(str);
S = fread(fd, [1 2], 'int32');
x.dat = fread(fd, fliplr(S), '*uint16')';
x.tt = [0];
x.seq = [];


function load2p_saveastiff(x)
[fn,pn] = uiputfile('*.tif', 'Specify name of TIF file to save', ...
    'C:\Users\Gradinaru Lab\Desktop');
if isempty(fn)
    disp('No selection made - canceling operation');
    return
end
fn = [pn filesep fn];

[Y X T] = size(x.dat);

mx = double(max(x.dat(:)));
compr = 'LZW';
imwrite(uint16(65535*double(x.dat(:,:,1))/mx), fn, 'TIFF', 'Compression', compr);
for t=2:T
    imwrite(uint16(65535*double(x.dat(:,:,t))/mx), fn, 'TIFF', 'Compression', compr, ...
        'WriteMode', 'append');
end
