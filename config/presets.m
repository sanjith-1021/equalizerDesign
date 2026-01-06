function cfg = presets(presetName)
%PRESETS Return configuration struct for HF modem demos.
if nargin < 1
    presetName = 'ber-sweep';
end

switch lower(presetName)
    case {'ber-sweep', 'hfmodemtop', 'top'}
        cfg = baseCfg();
        cfg.numDataSymbols = 16;
        cfg.numBlankSymbols = 64;
        cfg.snrDb = 0:4:30;
        cfg.numFrames = 100;
    case {'tap-anim', 'hfmodemanimtop', 'anim'}
        cfg = baseCfg();
        cfg.numDataSymbols = 64;
        cfg.numBlankSymbols = 8;
        cfg.snrDb = 25;
        cfg.numFrames = 100;
        cfg.animatePause = 0.2;
    otherwise
        error('Unknown preset: %s', presetName);
end
end

function cfg = baseCfg()
cfg.M = 8;
cfg.numPilotSymbols = 32;
cfg.samplesPerSymbol = 4;
cfg.rolloff = 0.25;
cfg.filterSpan = 8;
cfg.sampRate = 9600;
cfg.fMax = 1;
cfg.seed = 101;
cfg.nChanTaps = 5;
cfg.generators = [1 1 1 1 0 0 1;  % G0 = 171
                  1 0 1 1 0 1 1]; % G1 = 133
end
