% =========================================================
% Local function: Hard-decision Viterbi decoder (rate 1/2)
% =========================================================
function u_hat = viterbi_decode_rate12(rx_bits, next_state, out_bits)
% VITERBI_DECODE_RATE12
%   Hard-decision Viterbi decoder for a rate-1/2 convolutional code.
%
%   INPUTS
%     rx_bits   : Column or row vector of received coded bits (0/1).
%                 Length must be even, since we have 2 coded bits
%                 per information bit (rate 1/2).
%
%     next_state: Trellis next-state table.
%                 Size: [numStates x 2]
%                 next_state(s+1, u+1) = next state index (0..numStates-1)
%                 when current state = s and input bit u = 0 or 1.
%
%     out_bits  : Trellis output table.
%                 Size: [numStates x 2]
%                 out_bits(s+1, u+1) is a decimal number in {0,1,2,3}
%                 encoding the 2 output bits for branch (s, u).
%                 Mapping: 0->00, 1->01, 2->10, 3->11.
%
%   OUTPUT
%     u_hat     : Column vector of decoded information bits (0/1).
%                 Length = number of trellis steps (N).
%
%   NOTES
%     - This implementation uses HARD decisions and Hamming distance
%       as the branch metric.
%     - It does NOT assume zero-termination; final state is chosen as
%       the state with minimum path metric at time N.
%     - Path metrics are stored for all time steps to make the logic
%       easier to understand (not memory-optimized).

    % ------------------------------
    % Basic input formatting/checks
    % ------------------------------
    rx_bits = rx_bits(:);               % force column vector
    L = length(rx_bits);               % total number of received bits

    % For rate 1/2, we must have 2 bits per time step
    if mod(L,2) ~= 0
        error('Length of rx_bits must be even (2 bits per symbol).');
    end

    N = L / 2;                         % number of trellis steps (time instants)
    numStates = size(next_state, 1);   % number of states in the trellis

    % -------------------------------------------
    % Group received bits into N x 2 matrix r(k)
    % -------------------------------------------
    % r(k, :) = [r_k0 r_k1] are the two received bits at time k.
    r = zeros(N, 2);
    idx = 1;
    for k = 1:N
        r(k,1) = rx_bits(idx);
        r(k,2) = rx_bits(idx+1);
        idx = idx + 2;
    end

    % -------------------------------------------------
    % Allocate data structures for the Viterbi algorithm
    % -------------------------------------------------
    INF = 1e9;                         % "infinity" for impossible paths

    % PM(k, s+1) = path metric (cumulative Hamming distance) of the
    % best path that ends in state s at time step (k-1).
    %
    % We use N+1 rows:
    %   - PM(1,:)   : metrics BEFORE processing any symbols (time 0)
    %   - PM(k+1,:) : metrics AFTER processing symbol k
    %   - PM(N+1,:) : final metrics after N steps
    PM = INF * ones(N+1, numStates);
    PM(1,1) = 0;                       % start in state 0 with metric 0
                                       % (all other states start as impossible)

    % DEC(k, s+1)  : input bit (0 or 1) chosen to reach state s at time k
    % PREV(k, s+1) : previous state index (0..numStates-1) on the best path
    DEC  = zeros(N, numStates);
    PREV = zeros(N, numStates);

    % ===========================
    % Forward recursion (trellis)
    % ===========================
    % For each time step k, we:
    %   - take the current received bits r(k,:)
    %   - extend all survivor paths from time k-1 by inputs u in {0,1}
    %   - update the best path metric into each next state
    for k = 1:N  % time step (1..N)
        % Received bits at time k
        r0 = r(k,1);
        r1 = r(k,2);

        % Temporary array to hold new path metrics for time k
        PM_next = INF * ones(1, numStates);

        % Loop over all current states at time k-1
        for s = 0:numStates - 1
            % Best metric for state s at time k-1
            prev_metric = PM(k, s+1);

            % If this state has INF metric, no valid path reaches it,
            % so we skip it.
            if prev_metric >= INF
                continue;
            end

            % Try both possible input bits u = 0,1 from this state s.
            for u = 0:1
                % Trellis:
                %   s  : current state (0..numStates-1)
                %   u  : input bit (0/1)
                %   ns : next state
                ns = next_state(s+1, u+1);  % index 0..numStates-1

                % out is the encoded output label on this branch,
                % packed as two bits in a number 0..3.
                out = out_bits(s+1, u+1);   % 0..3 => [v0 v1]

                % Decode packed output bits:
                %   v0 = LSB of out
                %   v1 = next bit of out
                v0 = bitand(out, 1);
                v1 = bitand(bitshift(out, -1), 1);

                % ---------------------------------------------
                % Branch metric = Hamming distance between:
                %   received bits [r0 r1]
                %   expected bits [v0 v1]
                % ---------------------------------------------
                bm = (r0 ~= v0) + (r1 ~= v1);

                % New candidate path metric if we go from state s
                % to ns using input bit u
                metric = prev_metric + bm;

                % If this candidate is better (smaller metric)
                % than the current best for state ns, update it.
                if metric < PM_next(ns+1)
                    PM_next(ns+1) = metric;  % new best path metric into ns
                    DEC(k, ns+1)  = u;       % store chosen input bit
                    PREV(k, ns+1) = s;       % store previous state
                end
            end
        end

        % After considering all transitions at time k,
        % commit the new metrics as the metrics for time k.
        PM(k+1,:) = PM_next;
    end

    % =================
    % Traceback (output)
    % =================
    % At the end of the forward recursion, PM(N+1,:) holds the best
    % path metric for each possible final state.
    %
    % Since this implementation does not enforce a known final state
    % (like zero-termination), we choose the state with the smallest
    % final metric as the winning end state.
    [~, best_idx] = min(PM(N+1,:));
    state = best_idx - 1;  % convert from 1-based index to 0-based state

    % Recover the input bits by walking backwards through DEC/PREV
    u_hat = zeros(N,1);    % decoded bits (time 1..N)
    for k = N:-1:1
        % DEC(k, state+1) is the input bit that led into 'state'
        % at time k on the best path.
        u_hat(k) = DEC(k, state+1);

        % PREV(k, state+1) is the previous state on that best path.
        state    = PREV(k, state+1);
    end
end
