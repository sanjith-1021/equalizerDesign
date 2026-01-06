function bits = pnBits(nBits, taps, seed)
%PNBITS Generate a PN bit sequence using an LFSR.
if nargin < 2 || isempty(taps)
    taps = [10 7];
end
if nargin < 3 || isempty(seed)
    seed = 1;
end

m = max(taps);
state = de2bi(seed, m, 'left-msb');
if ~any(state)
    state(end) = 1;
end

bits = zeros(nBits, 1);
for k = 1:nBits
    bits(k) = state(end);
    feedback = mod(sum(state(taps)), 2);
    state = [feedback state(1:end-1)];
end
end
