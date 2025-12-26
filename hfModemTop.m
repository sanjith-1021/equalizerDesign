% HF modem end-to-end simulation with tap animation.
close all; clc; clear;

cfg.M = 8;
cfg.numDataSymbols = 16;
cfg.numPilotSymbols = 16;
cfg.numBlankSymbols = 16;
cfg.samplesPerSymbol = 4;
cfg.rolloff = 0.25;
cfg.filterSpan = 8;
cfg.sampRate = 9600;
cfg.numFrames = 100;
cfg.snrDb = 30;
cfg.fMax = 1;
cfg.seed = 101;
cfg.animatePause = 0.2;

generators = [1 1 1 1 0 0 1;  % G0 = 171
              1 0 1 1 0 1 1]; % G1 = 133
[next_state, out_bits] = build_trellis_rate12(generators);

rng(cfg.seed);
pilotBits = randi([0 1], cfg.numPilotSymbols * log2(cfg.M), 1);
pilotInts = bi2de(reshape(pilotBits, log2(cfg.M), []).', 'left-msb');
cfg.PilotSymbols = pskmod(pilotInts, cfg.M, pi / cfg.M);

chanLM = stdchan('iturHFLM', cfg.sampRate, cfg.fMax);
chanLM.RandomStream = 'mt19937ar with seed';
chanLM.Seed = 9999;
chanLM.PathGainsOutputPort = true;
delaySamps = round(chanLM.PathDelays(:) * cfg.sampRate);
tapCount = max(delaySamps) + 1; tapAxis = 0:(tapCount - 1);
channelModel = @(waveform) chanLM(awgn(waveform, cfg.snrDb, 'measured'));

modulator = TxModulator('Cfg', cfg);
demodLms = RxDemodulator('Cfg', cfg, 'EqualizerAlgorithm', 'LMS');
demodRls = RxDemodulator('Cfg', cfg, 'EqualizerAlgorithm', 'RLS');
demodKalman = RxDemodulator('Cfg', cfg, 'EqualizerAlgorithm', 'Kalman');

berLms = zeros(cfg.numFrames, 1);
berRls = zeros(cfg.numFrames, 1);
berKalman = zeros(cfg.numFrames, 1);
berCodedLms = zeros(cfg.numFrames, 1);
berCodedRls = zeros(cfg.numFrames, 1);
berCodedKalman = zeros(cfg.numFrames, 1);

tapFig = figure('Name', 'Channel taps vs estimates', 'NumberTitle', 'off');
tapAx = gobjects(3, 1);
tapLineChan = gobjects(3, 1);
tapLineEst = gobjects(3, 1);
algs = {'LMS', 'RLS', 'Kalman'};
for k = 1:3
    tapAx(k) = subplot(3, 1, k);
    tapLineChan(k) = plot(tapAx(k), tapAxis, zeros(1, tapCount), 'k-o', 'LineWidth', 1.1, 'MarkerSize', 3);
    hold(tapAx(k), 'on');
    tapLineEst(k) = plot(tapAx(k), tapAxis, zeros(1, tapCount), 'r.-', 'LineWidth', 1.1, 'MarkerSize', 8);
    grid(tapAx(k), 'on');
    ylabel(tapAx(k), '|tap|');
    title(tapAx(k), algs{k});
    if k == 3
        xlabel(tapAx(k), 'Tap Index');
    end
    legend(tapAx(k), {'Channel', 'Estimate'}, 'Location', 'northeast');
end

for frameIdx = 1:cfg.numFrames
    dataBits = randi([0 1], cfg.numDataSymbols * log2(cfg.M)/2, 1);
    codedBits = conv_encode_rate12(dataBits, generators);
    txWaveform = modulator(codedBits);

    [rxWaveform, pathGains] = channelModel(txWaveform);

    [rcvdBitsLms, eqSymbsLms, chanTapsLms] = demodLms(rxWaveform);
    [rcvdBitsRls, eqSymbsRls, chanTapsRls] = demodRls(rxWaveform);
    [rcvdBitsKalman, eqSymbsKalman, chanTapsKalman] = demodKalman(rxWaveform);

    rxBitsLms = viterbi_decode_rate12(rcvdBitsLms, next_state, out_bits);
    rxBitsRls = viterbi_decode_rate12(rcvdBitsRls, next_state, out_bits);
    rxBitsKalman = viterbi_decode_rate12(rcvdBitsKalman, next_state, out_bits);

    berLms(frameIdx) = mean(rcvdBitsLms ~= codedBits);
    berRls(frameIdx) = mean(rcvdBitsRls ~= codedBits);
    berKalman(frameIdx) = mean(rcvdBitsKalman ~= codedBits);
    berCodedLms(frameIdx) = mean(rxBitsLms ~= dataBits);
    berCodedRls(frameIdx) = mean(rxBitsRls ~= dataBits);
    berCodedKalman(frameIdx) = mean(rxBitsKalman ~= dataBits);

    meanGains = mean(pathGains, 1).';
    actualTaps = complex(zeros(tapCount, 1));
    for k = 1:numel(delaySamps)
        idx = delaySamps(k) + 1;
        actualTaps(idx) = actualTaps(idx) + meanGains(k);
    end

    estTaps = {chanTapsLms, chanTapsRls, chanTapsKalman};
    for k = 1:3
        est = complex(zeros(tapCount, 1));
        cnt = min(tapCount, numel(estTaps{k}));
        est(1:cnt) = estTaps{k}(1:cnt);
        set(tapLineChan(k), 'YData', abs(actualTaps));
        set(tapLineEst(k), 'YData', abs(est));
        maxVal = max([abs(actualTaps(:)); abs(est(:)); 1e-3]);
        ylim(tapAx(k), [0 maxVal * 1.1]);
    end

    tapFig.Name = sprintf('Channel taps vs estimates | SNR %.1f dB | Frame %d/%d', ...
        cfg.snrDb, frameIdx, cfg.numFrames);
    drawnow;
    if cfg.animatePause > 0
        pause(cfg.animatePause);
    end
end

berLmsMean = mean(berLms);
berRlsMean = mean(berRls);
berKalmanMean = mean(berKalman);
berCodedLmsMean = mean(berCodedLms);
berCodedRlsMean = mean(berCodedRls);
berCodedKalmanMean = mean(berCodedKalman);

figure;
semilogy(cfg.snrDb, berLmsMean,        '-o', 'LineWidth', 1.5, 'MarkerSize', 6); hold on;
semilogy(cfg.snrDb, berRlsMean,        '-s', 'LineWidth', 1.5, 'MarkerSize', 6);
semilogy(cfg.snrDb, berKalmanMean,     '-d', 'LineWidth', 1.5, 'MarkerSize', 6);
semilogy(cfg.snrDb, berCodedLmsMean,   '--^', 'LineWidth', 1.5, 'MarkerSize', 6);
semilogy(cfg.snrDb, berCodedRlsMean,   '--v', 'LineWidth', 1.5, 'MarkerSize', 6);
semilogy(cfg.snrDb, berCodedKalmanMean,'--x', 'LineWidth', 1.5, 'MarkerSize', 6);
xlabel('SNR (dB)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Bit Error Rate (BER)', 'FontSize', 12, 'FontWeight', 'bold');
title('BER at single SNR for LMS / RLS / Kalman (Coded & Uncoded)', 'FontSize', 13);
legend({'LMS', 'RLS', 'Kalman', 'Coded LMS', 'Coded RLS', 'Coded Kalman'}, ...
       'Location', 'southwest', 'FontSize', 10);
legend boxoff;
grid on;
grid minor;
