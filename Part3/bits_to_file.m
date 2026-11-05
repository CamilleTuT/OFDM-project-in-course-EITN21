function bytes = bits_to_file(bits, outputFile, numBytes)
%BITS_TO_FILE Convert a left-msb bitstream back to bytes and write a file.

bits = bits(:) ~= 0;
requiredBits = numBytes * 8;
if numel(bits) < requiredBits
    error('Not enough bits to write %d bytes.', numBytes);
end
bits = bits(1:requiredBits);

bytes = zeros(numBytes, 1, 'uint8');
for k = 1:8
    bytes = bitor(bytes, uint8(bits(k:8:end)) .* uint8(2^(8-k)));
end

outDir = fileparts(outputFile);
if ~isempty(outDir) && ~exist(outDir, 'dir')
    mkdir(outDir);
end

fid = fopen(outputFile, 'wb');
if fid < 0
    error('Cannot open output file: %s', outputFile);
end
cleaner = onCleanup(@() fclose(fid));
fwrite(fid, bytes, 'uint8');
end
