function codedBits = ldpc_encode_bits(bits, codeRate, modulation, encodedBitCount, rv)
%LDPC_ENCODE_BITS 5G NR LDPC encode and rate-match a transport block.

if nargin < 5
    rv = 0;
end

bits = double(bits(:) ~= 0);
encodedBitCount = double(encodedBitCount);
cbsInfo = nrDLSCHInfo(numel(bits), codeRate);

bitsWithCrc = nrCRCEncode(bits, cbsInfo.CRC);
codeBlocks = nrCodeBlockSegmentLDPC(bitsWithCrc, cbsInfo.BGN);
encodedBlocks = nrLDPCEncode(codeBlocks, cbsInfo.BGN);
codedBits = nrRateMatchLDPC(encodedBlocks, encodedBitCount, rv, modulation, 1);
codedBits = double(codedBits(:));
end
