classdef channelModel < matlab.System
    % channelModel  HF channel for Table XX cases (<= 2400 bps)
    %
    % Usage:
    %   cm = channelModel(cfg, caseId);
    %   [y, chanObj] = cm(x);

    properties (Nontunable, SetAccess = private)
        cfg
        caseId
    end

    properties (Access = private)
        caseConfig
        sampleRate
    end

    methods
        function obj = channelModel(cfg, caseId)
            obj.cfg = cfg;
            obj.caseId = caseId;
        end

        function [y, chanObj] = apply(obj, x)
            [y, chanObj] = obj(x);
        end
    end

    methods (Access = protected)
        function setupImpl(obj)
            cases = obj.getCases();
            obj.caseConfig = cases(obj.caseId);
            obj.sampleRate = obj.cfg.sampRate;
        end

        function [y, chanObj] = stepImpl(obj, x)
            c = obj.caseConfig;
            if ~c.isFading
                y = awgn(x, c.snrDb, 'measured');
                if nargout > 1
                    chanObj = [];
                end
                return;
            end

            [y, chanObj] = obj.applyFading(x, c);
        end
    end

    methods (Access = private)
        function [y, chanObj] = applyFading(obj, x, c)
            xIn = x(:);
            nSig = numel(xIn);

            f3dbNorm = 0.01;
            firLen = 257;
            fadeRate = 100 * c.fadingBwHz;

            h0 = obj.makeFadingFirNorm(obj.sampleRate, fadeRate, f3dbNorm, firLen, nSig, c.seed);
            h1 = obj.makeFadingFirNorm(obj.sampleRate, fadeRate, f3dbNorm, firLen, nSig, c.seed + 1);

            scale = 1 / sqrt(2);
            h0 = h0 * scale;
            h1 = h1 * scale;

            delaySamples = round(c.multipathMs * 1e-3 * obj.sampleRate);
            if delaySamples > 0
                xDelay = [zeros(delaySamples, 1); xIn(1:end - delaySamples)];
            else
                xDelay = xIn;
            end

            y = h0 .* xIn + h1 .* xDelay;
            y = awgn(y, c.snrDb, 'measured');

            if size(x, 1) == 1
                y = y.';
            end

            if nargout > 1
                chanObj = struct( ...
                    'pathDelays', [0, c.multipathMs * 1e-3], ...
                    'delaySamples', delaySamples, ...
                    'pathGains', [h0, h1], ...
                    'fadingBwHz', c.fadingBwHz, ...
                    'seed', c.seed, ...
                    'firLen', firLen, ...
                    'f3dbNorm', f3dbNorm);
            end
        end

        function hSig = makeFadingFirNorm(obj, FsSig, Ffade, f3dbNorm, L, Nsig, seed)
            stream = RandStream('mt19937ar', 'Seed', seed);
            b = obj.gaussFirNormF3db(f3dbNorm, L);
            Nfade = ceil(Nsig * Ffade / FsSig) + 2;
            w = (randn(stream, 1, Nfade) + 1j * randn(stream, 1, Nfade)) / sqrt(2);
            hFade = conv(w, b, 'same');
            hFade = hFade / sqrt(mean(abs(hFade).^2));

            tFade = (0:Nfade-1) / Ffade;
            tSig = (0:Nsig-1) / FsSig;
            hSig = interp1(tFade, hFade, tSig, 'linear', 'extrap').';
        end

        function b = gaussFirNormF3db(~, f3dbNorm, L)
            sigma = f3dbNorm / sqrt(log(2));
            Nfft = 2^nextpow2(16 * L);
            k = 0:Nfft-1;
            f = (k - floor(Nfft/2)) / Nfft;
            H = exp(-(f.^2) / (2 * sigma^2));
            H = ifftshift(H);
            hFull = real(ifft(H));
            hFull = fftshift(hFull);
            mid = floor(Nfft/2) + 1;
            idx = (mid - floor(L/2)) : (mid + floor(L/2));
            b = hFull(idx);
            b = b / sum(b);
        end
    end

    methods (Access = private, Static)
        function cases = getCases()
            add = @(isF, mp, bw, snr, seed) struct( ...
                'isFading', isF, ...
                'multipathMs', mp, ...
                'fadingBwHz', bw, ...
                'snrDb', snr, ...
                'seed', seed);

            cases = [ ...
                add(false, 0, 0, 10, 1001); ...         % 1: 2400, 1 fixed,  SNR = 10 dB          
                add(true,  2, 1, 30, 1002); ...         % 2: 2400, 2 fading, multipath = 2 ms, BW = 1 Hz,  SNR = 18 dB
                add(true,  2, 5, 30, 1003); ...         % 3: 2400, 2 fading, multipath = 2 ms, BW = 5 Hz,  SNR = 30 dB
                add(true,  5, 1, 30, 1004); ...         % 4: 2400, 2 fading, multipath = 5 ms, BW = 1 Hz,  SNR = 30 dB
                add(true,  2, 1, 11, 1005); ...         % 5: 1200, 2 fading, multipath = 2 ms, BW = 1 Hz,  SNR = 11 dB
                add(true,  2, 1, 7,  1006); ...         % 6: 600,  2 fading, multipath = 2 ms, BW = 1 Hz,  SNR = 7  dB
                add(true,  5, 5, 7,  1007); ...         % 7: 300,  2 fading, multipath = 5 ms, BW = 5 Hz,  SNR = 7  dB
                add(true,  5, 5, 5,  1008); ...         % 8: 150,  2 fading, multipath = 5 ms, BW = 5 Hz,  SNR = 5  dB
                add(true,  5, 5, 2,  1009)];            % 9: 75,   2 fading, multipath = 5 ms, BW = 5 Hz,  SNR = 2  dB
        end
    end
end
