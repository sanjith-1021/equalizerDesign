% HF modem tap animation with a static custom channel.
close all; clc; clear;

scriptDir = fileparts(mfilename('fullpath'));
addpath(fullfile(scriptDir, 'config'));
cfg = presets('tap-anim');

rng(cfg.seed);
pilotBits = randi([0 1], cfg.numPilotSymbols * log2(cfg.M), 1);
pilotInts = bi2de(reshape(pilotBits, log2(cfg.M), []).', 'left-msb');
cfg.PilotSymbols = pskmod(pilotInts, cfg.M, pi / cfg.M);

rrc = rcosdesign(cfg.rolloff, cfg.filterSpan, cfg.samplesPerSymbol, 'sqrt').';

% rng(cfg.seed);
% chanLM = stdchan('iturHFLM', cfg.sampRate, cfg.fMax);
% chanLM.RandomStream = 'mt19937ar with seed';
% chanLM.Seed = 9999;
% chanLM.PathGainsOutputPort = false;
% reset(chanLM);
% channelModel = @(waveform) chanLM(awgn(waveform, cfg.snrDb, 'measured'));

chanCoeffs = randn(cfg.nChanTaps, 1) + 1j * randn(cfg.nChanTaps, 1);
chanCoeffs = chanCoeffs / norm(chanCoeffs);
channelModel = @(waveform) conv(awgn(waveform, cfg.snrDb, 'measured'), chanCoeffs, 'same');

modulator = TxModulator('Cfg', cfg);
demodLms = RxDemodulator('Cfg', cfg, 'EqualizerAlgorithm', 'LMS');
demodRls = RxDemodulator('Cfg', cfg, 'EqualizerAlgorithm', 'RLS');
demodKalman = RxDemodulator('Cfg', cfg, 'EqualizerAlgorithm', 'Kalman');

berLms = zeros(cfg.numFrames, 1);
berRls = zeros(cfg.numFrames, 1);
berKalman = zeros(cfg.numFrames, 1);


tapFig = figure('Name', 'Error History', 'NumberTitle', 'off');
errCnt = cfg.numBlankSymbols * 2 + cfg.numPilotSymbols + cfg.numDataSymbols;
errAx = gobjects(3, 1);
errLineChan = gobjects(3, 1);
errLineEst = gobjects(3, 1);
algs = {'LMS', 'RLS', 'Kalman'};
errAxis = 0:errCnt-1;
for k = 1:3
    errAx(k) = subplot(3, 1, k);
    errLineEst(k) = plot(errAx(k), errAxis, zeros(1, errCnt), 'r.-', 'LineWidth', 1.1, 'MarkerSize', 8);
    grid(errAx(k), 'on');
    ylabel(errAx(k), '|err|');
    title(errAx(k), algs{k});
    if k == 3
        xlabel(errAx(k), 'Symbol Index');
    end
end

for frameIdx = 1:cfg.numFrames
    dataBits = randi([0 1], cfg.numDataSymbols * log2(cfg.M)/2, 1);
    codedBits = conv_encode_rate12(dataBits, cfg.generators);
    txWaveform = modulator(codedBits);

    rxWaveform = channelModel(txWaveform);

    [rcvdBitsLms, ~, errHistLms] = demodLms(rxWaveform);
    [rcvdBitsRls, ~, errHistRls] = demodRls(rxWaveform);
    [rcvdBitsKalman, ~, errHistKalman] = demodKalman(rxWaveform);

    berLms(frameIdx) = mean(rcvdBitsLms ~= codedBits);
    berRls(frameIdx) = mean(rcvdBitsRls ~= codedBits);
    berKalman(frameIdx) = mean(rcvdBitsKalman ~= codedBits);

    errHists = {errHistLms, errHistRls, errHistKalman};
    for k = 1:3
        set(errLineEst(k), 'YData', abs(errHists{k}));
        maxVal = max([abs(errHists{k}); 1e-3]);
        ylim(errAx(k), [0 maxVal * 1.1]);
    end

    tapFig.Name = sprintf('Error History | SNR %.1f dB | Frame %d/%d', ...
        cfg.snrDb, frameIdx, cfg.numFrames);
    drawnow;
    if cfg.animatePause > 0
        pause(cfg.animatePause);
    end
end
