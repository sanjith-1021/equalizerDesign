% =========================================================
% Local function: Build trellis for rate 1/2 code
% =========================================================
function [next_state, out_bits] = build_trellis_rate12(generators)
% BUILD_TRELLIS_RATE12
%   Build trellis tables for a rate-1/2 convolutional code.
%
%   INPUT
%     generators : 2 x K matrix of 0/1 (two generator polynomials)
%                  - Row j is g_j = [g_j0, g_j1, ..., g_j(K-1)]
%                  - g_j0 multiplies current input u(k)
%                  - g_j1..g_j(K-1) multiply memory bits.
%
%                  Example (K = 3):
%                    generators = [1 1 1;
%                                  1 0 1];
%
%   OUTPUTS
%     next_state : [numStates x 2] matrix
%                  next_state(s+1, u+1) = next state index (0..numStates-1)
%                  when current state = s and input bit = u.
%
%     out_bits   : [numStates x 2] matrix
%                  out_bits(s+1, u+1) = encoded output label for that branch,
%                  stored as a decimal number in {0,1,2,3} representing the
%                  2-bit output [v0 v1] as:
%                      0 -> 00
%                      1 -> 01
%                      2 -> 10
%                      3 -> 11
%
%   PURPOSE
%     These tables are used by the Viterbi decoder to:
%       - know which next state is reached from (state, input)
%       - know what 2 output bits are expected on that branch
%     so it can compute branch metrics and do path extension.

    % ------------------------------
    % Basic checks / parameter setup
    % ------------------------------
    [n, K] = size(generators);
    if n ~= 2
        error('This version supports only rate 1/2 (2 rows in generators).');
    end

    mem = K - 1;                  % number of memory elements in the shift register
    numStates = 2^mem;            % number of possible memory contents (states)

    % Allocate trellis tables:
    %   - next_state(s+1, u+1): next state index for input u from state s
    %   - out_bits(s+1, u+1)  : packed output bits for that branch (0..3)
    next_state = zeros(numStates, 2);  % inputs u = 0,1
    out_bits   = zeros(numStates, 2);  % packed output [v0 v1] as decimal

    % ======================================================
    % Loop over all possible states (memory configurations)
    % ======================================================
    % States are numbered from 0 to numStates-1.
    % Each state corresponds to a pattern of memory bits [m1 ... m_mem].
    % We need to:
    %   1) Decode s -> [m1 ... m_mem]
    %   2) For each input bit u in {0,1}, compute:
    %        - output bits [v0 v1]
    %        - next memory bits -> next state index
    for s = 0:numStates-1
        % ---------------------------------------------
        % Decode state index s into memory bits [m1..m_mem]
        % ---------------------------------------------
        % We interpret s as a binary number:
        %   s = b_(mem-1)*2^(mem-1) + ... + b_1*2 + b_0
        % and map these bits into mem_bits.
        mem_bits = zeros(1, mem);
        tmp = s;
        for i = mem:-1:1
            mem_bits(i) = bitand(tmp, 1);   % LSB of tmp
            tmp = bitshift(tmp, -1);        % shift right
        end
        % After this loop:
        %   mem_bits = [m1 m2 ... m_mem] (oldest to newest)

        % ==========================================
        % For each possible input bit u = 0 or 1
        % ==========================================
        for u = 0:1
            % Full shift register contents at this transition:
            %   reg(1)     = current input u
            %   reg(2:end) = current memory bits [m1..m_mem]
            reg = [u, mem_bits];

            % ----------------------------------------------
            % Compute the two output bits for this branch.
            % For each generator row j:
            %   v(j) = (sum over t of generators(j,t) * reg(t)) mod 2
            % ----------------------------------------------
            v = zeros(1,2);  % [v0 v1]
            for j = 1:2
                acc = 0;     % accumulator for mod-2 sum (XOR)
                for k = 1:K
                    if generators(j,k) == 1
                        acc = xor(acc, reg(k));
                    end
                end
                v(j) = acc;
            end

            % ----------------------------------------------
            % Compute new memory contents for the next state
            % ----------------------------------------------
            % After shifting:
            %   new_mem_bits = [u, m1, ..., m_(mem-1)]
            new_mem_bits = reg(1:mem);

            % Convert new_mem_bits back into an integer state index ns
            % in the range 0..numStates-1.
            ns = 0;
            for i = 1:mem
                ns = bitor(bitshift(ns,1), new_mem_bits(i));
            end

            % Store in trellis tables:
            %   - next_state: the next state index reached
            %   - out_bits  : the output label as a decimal code
            next_state(s+1, u+1) = ns;
            out_bits(s+1, u+1)   = v(1) + 2*v(2);  % pack [v0 v1] into 0..3
        end
    end
end
