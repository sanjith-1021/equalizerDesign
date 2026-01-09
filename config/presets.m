function cfg = presets(presetName)
%PRESETS Return configuration struct for HF modem demos.
if nargin < 1
    presetName = 'ber-sweep';
end

switch lower(presetName)
    case {'ber-sweep'}
        cfg = baseCfg();

        cfg.snrDb = 0:4:30;
        cfg.nSlots = 1000;


    case {'tap-anim'}
        cfg = baseCfg();
        cfg.snrDb = 30;
        cfg.nSlots = 20;
        cfg.animatePause = 0.01;
        

    case {'basic-run'}
        cfg = baseCfg();
        cfg.snrDb = 30;

    case {'intrlvd-run'}
        cfg = baseCfg();
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



    otherwise
        error('Unknown preset: %s', presetName);
end
end

function cfg = baseCfg()
    cfg.M = 8;
    cfg.sampsPerSymb = 2;
    cfg.rolloff = 0.25;
    cfg.filterSpan = 12;
    cfg.sampRate = 9600;
    cfg.nHops = 1024;
    cfg.nBlankSymbs = 32;
    cfg.nTrngSymbs01 = 256;
    cfg.nTrngSymbs02 = 16;
    cfg.nDataSymbs = 32;
    cfg.seed = 101;
    cfg.generators = [1 1 1 1 0 0 1;  % G0 = 171
                    1 0 1 1 0 1 1]; % G1 = 133
    cfg.eq = struct( ...
        'nFeedforwardTaps', 48, ...
        'nFeedbackTaps', 12, ...
        'train', struct( ...
            'stepSize', 1e-5, ...
            'initInvCorr', 1e-2, ...
            'forgetFactor', 0.990, ...
            'measurementNoise', 1e-5, ...
            'processNoise', 5e-3 ...
            ), ...
        'data', struct( ...
            'stepSize', 2e-5, ...
            'forgetFactor', 0.995, ...
            'measurementNoise', 1e0, ...
            'processNoise', 1e-4) ...
        );
    cfg.fec = struct();

end
