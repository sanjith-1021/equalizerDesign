% =========================================================
% Local function: Convolutional encoder (rate 1/2)
% =========================================================
function coded = conv_encode_rate12(u, generators)
% CONV_ENCODE_RATE12
%   Simple rate-1/2 convolutional encoder with configurable generators.
%
%   INPUTS
%     u          : Column or row vector of input bits (0/1).
%                  Length = N (number of information bits).
%
%     generators : 2 x K matrix of binary generator coefficients (0/1).
%                  - 2 rows because this is a rate 1/2 encoder
%                    (2 output bits for each input bit).
%                  - K columns, where K = constraint length.
%                    Each row j is: [g_j0 g_j1 ... g_j(K-1)]
%                    g_j0 multiplies current input bit,
%                    g_j1..g_j(K-1) multiply memory bits.
%
%                  Example (K = 3):
%                    generators = [1 1 1;
%                                  1 0 1];
%
%   OUTPUT
%     coded      : Column vector of encoded bits (0/1).
%                  Length = 2 * N (2 output bits per input bit).
%
%   NOTES
%     - This version does NOT add flush bits (no zero-termination).
%       The encoder simply processes the N bits in u.
%     - Memory is initialized to all zeros at the start.
%     - XOR is used for mod-2 addition.

    % ------------------------------
    % Basic input formatting/checks
    % ------------------------------
    u = u(:);                      % force column vector
    [n, K] = size(generators);     % n = number of output streams, K = constraint length

    % This implementation supports only rate 1/2 (two output bits per input bit)
    if n ~= 2
        error('This version supports only rate 1/2 (2 rows in generators).');
    end

    mem = K - 1;                   % number of memory elements
    N   = length(u);               % number of information bits

    % ------------------------------------------------------
    % Memory bits [m1 ... m_mem], initialized to all zeros.
    % At any time k, the shift register contents are:
    %   reg = [current_input  m1  m2  ...  m_mem]
    % ------------------------------------------------------
    mem_bits = zeros(1, mem);

    % Preallocate output: 2 bits per input bit
    coded   = zeros(2*N, 1);
    out_idx = 1;                   % index into 'coded' vector

    % ===========================
    % Encoding loop over all bits
    % ===========================
    for k = 1:N
        % Current input bit at time k
        in = u(k);

        % Full shift register content at this time:
        %   reg(1) = current input bit
        %   reg(2:end) = previous memory bits
        reg = [in, mem_bits];

        % ----------------------------------------------
        % Compute the two output bits for this time k.
        % For each output j (1 or 2):
        %   v(j) = sum(generators(j,:) .* reg) mod 2
        % Here we implement it explicitly using XOR.
        % ----------------------------------------------
        v = zeros(1,2);            % [v1 v2]
        for j = 1:2
            acc = 0;               % accumulator for mod-2 sum
            for t = 1:K
                if generators(j,t) == 1
                    % XOR is addition in GF(2)
                    acc = xor(acc, reg(t));
                end
            end
            v(j) = acc;            % store computed output bit
        end
        
        % Store the two output bits in the coded sequence
        coded(out_idx)   = v(1);
        coded(out_idx+1) = v(2);
        out_idx = out_idx + 2;

        % ---------------------------------------------
        % Update memory for next time step:
        %   - New memory is the first 'mem' bits of reg.
        %   - Equivalent to shifting reg right and dropping
        %     the oldest bit.
        %
        % Example (mem = 2):
        %   Before: mem_bits = [m1 m2]
        %   reg    = [in m1 m2]
        %   After:  mem_bits = [in m1]
        % ---------------------------------------------
        mem_bits = reg(1:mem);
    end
end
