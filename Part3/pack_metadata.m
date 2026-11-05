function bits = pack_metadata(meta, cfg)
%PACK_METADATA Build the fixed-length overhead transport block.

bits = [
    uint_to_bits(cfg.metadataMagic, 16)
    uint_to_bits(cfg.metadataVersion, 8)
    uint_to_bits(1, 4)
    uint_to_bits(meta.fileBitLength, 32)
    uint_to_bits(meta.fileByteLength, 32)
    uint_to_bits(meta.payloadOfdmSymbols, 16)
    uint_to_bits(meta.payloadEncodedBitCount, 32)
];

if numel(bits) ~= cfg.metadataBitLength
    error('Metadata length mismatch: expected %d bits, got %d bits.', ...
        cfg.metadataBitLength, numel(bits));
end
end
