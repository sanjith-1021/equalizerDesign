function [channelModel, chanObj] = channelModel(cfg, mode, varargin)
    % channelModel  HF channel factory for MIL-STD-188-110B-style tests
    % 
    % Usage examples:
    %   [cm, obj] = channelModel(cfg, 'table', caseId);  % Table-XX cases <= 2400 bps
    %   [cm, obj] = channelModel(cfg, 'lm', snrDb);      % ITU-R HFLM Watterson channel
    %   [cm, obj] = channelModel(cfg, 'toy', snrDb);     % simple static FIR toy channel
    %   [cm, obj] = channelModel(cfg, 'toy', snrDb, nChanTaps);
    %
    %   y = cm(x);   % apply selected channel model to waveform x

    switch lower(mode)
        case 'table'
            caseId = varargin{1};
            [channelModel, chanObj] = tableChannel(cfg, caseId);

        case 'lm'
            snrDb = varargin{1};
            [channelModel, chanObj] = lmChannel(cfg, snrDb);

        case 'toy'
            snrDb = varargin{1};
            [channelModel, chanObj] = toyChannel(cfg, snrDb);

        otherwise
            error('Unknown mode "%s". Use ''table'', ''lm'', or ''toy''.', mode);
    end

end

% -------------------------------------------------------------------------
function [channelModel, chanObj] = tableChannel(cfg, caseId)

    cases = getCases();
    c = cases(caseId);

    if ~c.isFading
        chanObj      = [];
        channelModel = @(x) awgn(x, c.snrDb, 'measured');
        return;
    end

    pathDelays   = [0, c.multipathMs * 1e-3];
    avgPathGains = [0 0];

    fdMax     = c.fadingBwHz;   % fading BW = 2*sigma, fdMax = BW => sigmaNorm = 0.5
    sigmaNorm = 0.5;
    doppSpec  = doppler('Gaussian', sigmaNorm);

    chanObj = comm.RayleighChannel( ...
        'SampleRate',          cfg.sampRate, ...
        'PathDelays',          pathDelays, ...
        'AveragePathGains',    avgPathGains, ...
        'NormalizePathGains',  true, ...
        'FadingTechnique',     'Filtered Gaussian noise', ...
        'MaximumDopplerShift', fdMax, ...
        'DopplerSpectrum',     doppSpec, ...
        'RandomStream',        'mt19937ar with seed', ...
        'Seed',                c.seed, ...
        'PathGainsOutputPort', false);

    channelModel = @(x) awgn(chanObj(x), c.snrDb, 'measured');

end

% -------------------------------------------------------------------------
function [channelModel, chanObj] = lmChannel(cfg, snrDb)

    fdMax  = 1;  % Hz
    chanObj = stdchan('iturHFLM', cfg.sampRate, fdMax);
    chanObj.RandomStream        = 'mt19937ar with seed';
    chanObj.Seed                = 9999;
    chanObj.PathGainsOutputPort = false;

    channelModel = @(x) awgn(chanObj(x), snrDb, 'measured');

end

% -------------------------------------------------------------------------
function [channelModel, chanObj] = toyChannel(~, snrDb)
    nChanTaps = 4;
    chanCoeffs = randn(nChanTaps, 1) + 1j * randn(nChanTaps, 1);
    chanCoeffs = chanCoeffs / norm(chanCoeffs);

    chanObj = chanCoeffs;
    channelModel = @(x) conv(awgn(x, snrDb, 'measured'), chanCoeffs, 'same');

end

% -------------------------------------------------------------------------
function cases = getCases()

    add = @(isF, mp, bw, snr, seed) struct( ...
        'isFading',    isF, ...
        'multipathMs', mp, ...
        'fadingBwHz',  bw, ...
        'snrDb',       snr, ...
        'seed',        seed);

    % Table XX rows for user bit rates <= 2400 bps
    % 1: 2400, 1 fixed,  SNR = 10 dB
    % 2: 2400, 2 fading, multipath = 2 ms, BW = 1 Hz,  SNR = 18 dB
    % 3: 2400, 2 fading, multipath = 2 ms, BW = 5 Hz,  SNR = 30 dB
    % 4: 2400, 2 fading, multipath = 5 ms, BW = 1 Hz,  SNR = 30 dB
    % 5: 1200, 2 fading, multipath = 2 ms, BW = 1 Hz,  SNR = 11 dB
    % 6: 600,  2 fading, multipath = 2 ms, BW = 1 Hz,  SNR = 7  dB
    % 7: 300,  2 fading, multipath = 5 ms, BW = 5 Hz,  SNR = 7  dB
    % 8: 150,  2 fading, multipath = 5 ms, BW = 5 Hz,  SNR = 5  dB
    % 9: 75,   2 fading, multipath = 5 ms, BW = 5 Hz,  SNR = 2  dB

    cases = [ ...
        add(false, 0, 0, 10, 1001); ...
        add(true,  2, 1, 18, 1002); ...
        add(true,  2, 5, 30, 1003); ...
        add(true,  5, 1, 30, 1004); ...
        add(true,  2, 1, 11, 1005); ...
        add(true,  2, 1, 7,  1006); ...
        add(true,  5, 5, 7,  1007); ...
        add(true,  5, 5, 5,  1008); ...
        add(true,  5, 5, 2,  1009)];

end
