function [bits, crcError, info] = ldpc_decode_bits(llrs, transportBlockLength, ...
    codeRate, modulation, maxIterations, rv)
%LDPC_DECODE_BITS Rate recover, LDPC decode, desegment, and CRC check.

if nargin < 6
    rv = 0;
end

llrs = double(llrs(:));
cbsInfo = nrDLSCHInfo(transportBlockLength, codeRate);
rateRecovered = nrRateRecoverLDPC(llrs, transportBlockLength, codeRate, ...
    rv, modulation, 1, cbsInfo.C);
decodedBlocks = nrLDPCDecode(rateRecovered, cbsInfo.BGN, maxIterations);
[transportBlockWithCrc, segmentError] = nrCodeBlockDesegmentLDPC( ...
    decodedBlocks, cbsInfo.BGN, transportBlockLength + cbsInfo.L);
[bits, crcError] = nrCRCDecode(transportBlockWithCrc, cbsInfo.CRC);

bits = double(bits(1:transportBlockLength));
info.cbsInfo = cbsInfo;
info.segmentError = segmentError;
end
