% HF modem tap animation with a static custom channel.
close all; clc; clear;

scriptDir = fileparts(mfilename('fullpath'));
addpath(fullfile(scriptDir, 'config'));
cfg = presets('tap-anim');

rng(cfg.seed);
rng(cfg.seed);
trngBits01 = randi([0 1], cfg.nTrngSymbs01 * log2(cfg.M), 1);
trngBits02 = randi([0 1], cfg.nTrngSymbs02 * log2(cfg.M), 1);
trngInts01 = bi2de(reshape(trngBits01, log2(cfg.M), []).', 'left-msb');
trngInts02 = bi2de(reshape(trngBits02, log2(cfg.M), []).', 'left-msb');
cfg.TrngSeq01 = pskmod(trngInts01, cfg.M, pi / cfg.M);
cfg.TrngSeq02 = pskmod(trngInts02, cfg.M, pi / cfg.M);


rrc = rcosdesign(cfg.rolloff, cfg.filterSpan, cfg.sampsPerSymb, 'sqrt').';

chan = channelModel(cfg, 2);

modulator = TxModulator('Cfg', cfg);
demodLms = RxDemodulator('Cfg', cfg, 'EqualizerAlgorithm', 'LMS');
demodRls = RxDemodulator('Cfg', cfg, 'EqualizerAlgorithm', 'RLS');
demodKalman = RxDemodulator('Cfg', cfg, 'EqualizerAlgorithm', 'Kalman');

berLms = zeros(cfg.nSlots, 1);
berRls = zeros(cfg.nSlots, 1);
berKalman = zeros(cfg.nSlots, 1);


tapFig = figure('Name', 'Equalizer View', 'NumberTitle', 'off');
errCnt = numel(modulator.SlotSymbType);
errAx = gobjects(3, 1);
errLineEst = gobjects(3, 1);
scatAx = gobjects(3, 1);
scatLine = gobjects(3, 1);
algs = {'LMS', 'RLS', 'Kalman'};
errAxis = 0:errCnt-1;
for k = 1:3
    errAx(k) = subplot(3, 2, (k-1)*2 + 1);
    errLineEst(k) = plot(errAx(k), errAxis, zeros(1, errCnt), 'r.-', 'LineWidth', 1.1, 'MarkerSize', 8);
    grid(errAx(k), 'on');
    ylabel(errAx(k), '|err|');
    title(errAx(k), sprintf('%s | Error', algs{k}));
    xlabel(errAx(k), 'Symb Index');

    scatAx(k) = subplot(3, 2, (k-1)*2 + 2);
    scatLine(k) = plot(scatAx(k), 0, 0, 'b.');
    axis(scatAx(k), 'equal');
    grid(scatAx(k), 'on');
    title(scatAx(k), sprintf('%s | Scatter', algs{k}));
    xlabel(scatAx(k), 'I');
    ylabel(scatAx(k), 'Q');
end

for frameIdx = 1:cfg.nSlots
    dataBits = randi([0 1], cfg.nHops * cfg.nDataSymbs * log2(cfg.M)/2, 1);
    codedBits = conv_encode_rate12(dataBits, cfg.generators);
    txWaveform = modulator(codedBits);

    [rxWaveform,~] = chan(txWaveform);

    [rcvdBitsLms, eqSymbsLms, errHistLms] = demodLms(rxWaveform);
    [rcvdBitsRls, eqSymbsRls, errHistRls] = demodRls(rxWaveform);
    [rcvdBitsKalman, eqSymbsKalman, errHistKalman] = demodKalman(rxWaveform);

    berLms(frameIdx) = mean(rcvdBitsLms ~= codedBits);
    berRls(frameIdx) = mean(rcvdBitsRls ~= codedBits);
    berKalman(frameIdx) = mean(rcvdBitsKalman ~= codedBits);

    errHists = {errHistLms, errHistRls, errHistKalman};
    dataSymbs = { ...
        eqSymbsLms(demodLms.slotSymbType == 3), ...
        eqSymbsRls(demodRls.slotSymbType == 3), ...
        eqSymbsKalman(demodKalman.slotSymbType == 3)};
    for k = 1:3
        set(errLineEst(k), 'YData', abs(errHists{k}));
        maxVal = max([abs(errHists{k}); 1e-3]);
        ylim(errAx(k), [0 maxVal * 1.1]);
        set(scatLine(k), 'XData', real(dataSymbs{k}), 'YData', imag(dataSymbs{k}));
    end

    tapFig.Name = sprintf('Equalizer View | SNR %.1f dB | Frame %d/%d', ...
        cfg.snrDb, frameIdx, cfg.nSlots);
    drawnow;
    if cfg.animatePause > 0
        pause(cfg.animatePause);
    end
end
