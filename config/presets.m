function cfg = presets()

    cfg.M = 8;
    cfg.sampsPerSymb = 2;
    cfg.rolloff = 0.25;
    cfg.filterSpan = 20;
    cfg.sampRate = 9600;
    cfg.nHops = 240;
    cfg.nSlots = 10;
    cfg.nBlankSymbs = 32;
    cfg.nTrngSymbs01 = 360;
    cfg.nTrngSymbs02 = 16;
    cfg.nDataSymbs = 32;
    cfg.seed = 101;
    cfg.generators = [1 1 1 1 0 0 1;  % G0 = 171
                    1 0 1 1 0 1 1]; % G1 = 133
    cfg.eq = struct( ...
        'nFeedforwardTaps', 32, ...
        'nFeedbackTaps', 24, ...
        'train', struct( ...
            'initInvCorr', 1e-2, ...
            'processNoise', 1e-3, ...
            'measurementNoise', 1e-1 ...
            ), ...
        'data', struct( ...
            'processNoise', 1e-3, ...
            'measurementNoise', 4e-1) ...
    );

    cfg.fec = struct();
    cfg.fec.enable = true;
    cfg.fec.rate = 1/2;
    cfg.fec.constraintLength = 7;
    cfg.fec.generatorOctal = [171 133];
    cfg.fec.trellis = poly2trellis(cfg.fec.constraintLength, cfg.fec.generatorOctal);
    cfg.fec.tblen = 5 * (cfg.fec.constraintLength - 1);
    cfg.fec.useSoft = true;     % true=soft demap+soft viterbi, false=hard
    cfg.fec.nSoft = 3;          % soft metric bits for Viterbi
    cfg.fec.llrClip = 8;  

    cfg.nHops = 240;

end
