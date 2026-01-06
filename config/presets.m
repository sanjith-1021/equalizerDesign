function cfg = presets(presetName)
%PRESETS Return configuration struct for HF modem demos.
if nargin < 1
    presetName = 'ber-sweep';
end

switch lower(presetName)
    case {'ber-sweep'}
        cfg = baseCfg();
        cfg.nTrngSymbs02 = 16;
        cfg.nDataSymbs = 32;
        cfg.snrDb = 0:4:30;
        cfg.nSlots = 10;
        cfg.M = 8;

    case {'tap-anim'}
        cfg = baseCfg();
        cfg.nTrngSymbs02 = 20;
        cfg.nDataSymbs = 20;
        cfg.snrDb = 30;
        cfg.nSlots = 100;
        cfg.animatePause = 0.2;
        cfg.M = 8;

    case {'basic-run'}
        cfg = baseCfg();
        cfg.nTrngSymbs02 = 16;
        cfg.nDataSymbs = 32;
        cfg.snrDb = 30;
        cfg.nSlots = 2;
        cfg.M = 8;

    otherwise
        error('Unknown preset: %s', presetName);
end
end

function cfg = baseCfg()
cfg.nBlankSymbs = 8;
cfg.nTrngSymbs01 = 256;
cfg.nHops = 30;
cfg.sampsPerSymb = 4;
cfg.rolloff = 0.25;
cfg.filterSpan = 8;
cfg.sampRate = 9600;
cfg.fMax = 1;
cfg.seed = 101;
cfg.generators = [1 1 1 1 0 0 1;  % G0 = 171
                  1 0 1 1 0 1 1]; % G1 = 133
end
