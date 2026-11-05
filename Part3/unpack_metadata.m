function meta = unpack_metadata(bits, cfg)
%UNPACK_METADATA Parse the fixed-length overhead transport block.

bits = bits(:);
if numel(bits) < cfg.metadataBitLength
    error('Metadata block is too short.');
end
bits = bits(1:cfg.metadataBitLength);

pos = 1;
[magic, pos] = take_uint(bits, pos, 16);
[version, pos] = take_uint(bits, pos, 8);
[modId, pos] = take_uint(bits, pos, 4);
[fileBitLength, pos] = take_uint(bits, pos, 32);
[fileByteLength, pos] = take_uint(bits, pos, 32);
[payloadOfdmSymbols, pos] = take_uint(bits, pos, 16);
[payloadEncodedBitCount, ~] = take_uint(bits, pos, 32);

if uint32(magic) ~= cfg.metadataMagic
    error('Metadata magic mismatch. Packet may be decoded at the wrong position.');
end
if uint8(version) ~= cfg.metadataVersion
    error('Unsupported metadata version: %d.', double(version));
end
if double(modId) ~= 1
    error('Unsupported metadata modulation id: %d.', double(modId));
end

meta.magic = uint32(magic);
meta.version = uint8(version);
meta.fileBitLength = double(fileBitLength);
meta.fileByteLength = double(fileByteLength);
meta.payloadOfdmSymbols = double(payloadOfdmSymbols);
meta.payloadEncodedBitCount = double(payloadEncodedBitCount);
end

function [value, nextPos] = take_uint(bits, pos, width)
slice = bits(pos:pos + width - 1);
value = bits_to_uint(slice);
nextPos = pos + width;
end
