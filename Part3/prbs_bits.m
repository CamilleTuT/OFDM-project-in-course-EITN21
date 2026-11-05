function seq = prbs_bits(nBits, seed)
%PRBS_BITS Generate a deterministic x^15 + x^14 + 1 PRBS sequence.

state = false(15, 1);
seed = uint32(seed);
for k = 1:15
    state(k) = bitget(seed, 16 - k) ~= 0;
end
if ~any(state)
    state(1) = true;
end

seq = false(nBits, 1);
for k = 1:nBits
    seq(k) = state(end);
    feedback = xor(state(15), state(14));
    state = [feedback; state(1:end-1)];
end

seq = double(seq);
end
