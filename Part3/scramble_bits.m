function scrambledBits = scramble_bits(bits, seed)
%SCRAMBLE_BITS XOR input bits with a deterministic PRBS sequence.

bits = bits(:) ~= 0;
seq = prbs_bits(numel(bits), seed) ~= 0;
scrambledBits = double(xor(bits, seq));
end
