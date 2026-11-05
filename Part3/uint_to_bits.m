function bits = uint_to_bits(value, width)
%UINT_TO_BITS Convert an unsigned integer to a left-msb column bit vector.

value = uint64(value);
bits = zeros(width, 1);
for k = 1:width
    bits(k) = double(bitget(value, width - k + 1));
end
end
