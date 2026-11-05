function value = bits_to_uint(bits)
%BITS_TO_UINT Convert a left-msb bit vector to uint64.

bits = bits(:) ~= 0;
value = uint64(0);
for k = 1:numel(bits)
    value = bitshift(value, 1) + uint64(bits(k));
end
end
