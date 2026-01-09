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
        cfg.animatePause = 0.2;
        

    case {'basic-run'}
        cfg = baseCfg();
        cfg.snrDb = 30;


    otherwise
        error('Unknown preset: %s', presetName);
end
end

function cfg = baseCfg()
    cfg.M = 8;
    cfg.sampsPerSymb = 2;
    cfg.rolloff = 0.35;
    cfg.filterSpan = 12;
    cfg.sampRate = 9600;
    cfg.nHops = 256;
    cfg.nBlankSymbs = 32;
    cfg.nTrngSymbs01 = 512;
    cfg.nTrngSymbs02 = 20;
    cfg.nDataSymbs = 20;
    cfg.seed = 101;
    cfg.generators = [1 1 1 1 0 0 1;  % G0 = 171
                    1 0 1 1 0 1 1]; % G1 = 133
    cfg.eq = struct( ...
        'nFeedforwardTaps', 40, ...
        'nFeedbackTaps', 12, ...
        'train', struct( ...
            'stepSize', 1e-2, ...
            'initInvCorr', 5e-2, ...
            'forgetFactor', 0.965, ...
            'measurementNoise', 1e-4, ...
            'processNoise', 1e-2 ...
            ), ...
        'data', struct( ...
            'stepSize', 2e-4, ...
            'forgetFactor', 0.985, ...
            'measurementNoise', 1e-2, ...
            'processNoise', 1e-4) ...
        );
    

end
