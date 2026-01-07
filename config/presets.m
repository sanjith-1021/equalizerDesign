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
        cfg.nSlots = 100;
        cfg.animatePause = 0.2;
        
        cfg.eq.forgetFactor = 0.99;
        cfg.eq.initInvCorr = 1e-2;



    case {'basic-run'}
        cfg = baseCfg();
        cfg.snrDb = 30;


    otherwise
        error('Unknown preset: %s', presetName);
end
end

function cfg = baseCfg()
    cfg.M = 8;
    cfg.sampsPerSymb = 4;
    cfg.rolloff = 0.25;
    cfg.filterSpan = 8;
    cfg.sampRate = 9600;
    cfg.nHops = 30;
    cfg.nBlankSymbs = 8;
    cfg.nTrngSymbs01 = 256;
    cfg.nTrngSymbs02 = 16;
    cfg.nDataSymbs = 32;
    cfg.seed = 101;
    cfg.generators = [1 1 1 1 0 0 1;  % G0 = 171
                    1 0 1 1 0 1 1]; % G1 = 133
    cfg.eq = struct( ...
        'nFeedforwardTaps', 32, ...
        'nFeedbackTaps', 8, ...
        'stepSize', 1e-3, ...
        'forgetFactor', 0.95, ...
        'initInvCorr', 5e-2, ...
        'measurementNoise', 1e-3, ...
        'processNoise', 1e-4...
        );
    

end
