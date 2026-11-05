function [bits, info] = file_to_bits(filePath)
%FILE_TO_BITS Read a binary file and convert it to a left-msb bitstream.

fid = fopen(filePath, 'rb');
if fid < 0
    error('Cannot open input file: %s', filePath);
end
cleaner = onCleanup(@() fclose(fid));
bytes = fread(fid, Inf, 'uint8=>uint8');

bits = zeros(numel(bytes) * 8, 1);
for k = 1:8
    bits(k:8:end) = double(bitget(bytes, 9 - k));
end

[folder, name, ext] = fileparts(filePath);
info.path = filePath;
info.folder = folder;
info.name = [name ext];
info.numBytes = numel(bytes);
info.numBits = numel(bits);
end
