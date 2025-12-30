% HF modem end-to-end BER sweep (no animation).
close all; clc; clear;

cfg.M = 8;
cfg.numDataSymbols = 16;
cfg.numPilotSymbols = 32;
cfg.numBlankSymbols = 64;
cfg.samplesPerSymbol = 4;
cfg.rolloff = 0.25;
cfg.filterSpan = 8;
cfg.snrDb = 0:4:30;
cfg.sampRate = 9600;
cfg.numFrames = 100;
cfg.fMax = 1;
cfg.seed = 101;

generators = [1 1 1 1 0 0 1;  % G0 = 171
              1 0 1 1 0 1 1]; % G1 = 133
[next_state, out_bits] = build_trellis_rate12(generators);

rng(cfg.seed);
pilotBits = randi([0 1], cfg.numPilotSymbols * log2(cfg.M), 1);
pilotInts = bi2de(reshape(pilotBits, log2(cfg.M), []).', 'left-msb');
cfg.PilotSymbols = pskmod(pilotInts, cfg.M, pi / cfg.M);

snrPoints = cfg.snrDb(:);
numFrames = cfg.numFrames;
rndSeed = cfg.seed;

berLmsMean = zeros(numel(snrPoints),1);
berRlsMean = zeros(numel(snrPoints),1);
berKalmanMean = zeros(numel(snrPoints), 1);
berCodedLmsMean = zeros(numel(snrPoints), 1);
berCodedRlsMean = zeros(numel(snrPoints), 1);
berCodedKalmanMean = zeros(numel(snrPoints), 1);

parfor sIdx = 1:numel(snrPoints)
    snrDb = snrPoints(sIdx);

    rng(rndSeed);
    chanLM = stdchan('iturHFLM', cfg.sampRate, cfg.fMax);
    chanLM.RandomStream = 'mt19937ar with seed';
    chanLM.Seed = 9999;
    chanLM.PathGainsOutputPort = false;
    reset(chanLM);
    channelModel = @(waveform) chanLM(awgn(waveform, snrDb, 'measured'));

    % nChanTaps = 5;
    % chanCoeffs = randn(nChanTaps, 1) + 1j * randn(nChanTaps, 1);
    % chanCoeffs = chanCoeffs / norm(chanCoeffs);
    % channelModel = @(waveform) conv(awgn(waveform, snrDb, 'measured'), chanCoeffs, 'same');

    modulator = TxModulator('Cfg', cfg);
    demodLms = RxDemodulator('Cfg', cfg, 'EqualizerAlgorithm', 'LMS');
    demodRls = RxDemodulator('Cfg', cfg, 'EqualizerAlgorithm', 'RLS');
    demodKalman = RxDemodulator('Cfg', cfg, 'EqualizerAlgorithm', 'Kalman');

    berLms = zeros(1, numFrames);
    berRls = zeros(1, numFrames);
    berKalman = zeros(1, numFrames);
    berCodedLms = zeros(1, numFrames);
    berCodedRls = zeros(1, numFrames);
    berCodedKalman = zeros(1, numFrames);

    for frameIdx = 1:numFrames
        dataBits = randi([0 1], cfg.numDataSymbols * log2(cfg.M)/2, 1);
        codedBits = conv_encode_rate12(dataBits, generators);
        txWaveform = modulator(codedBits);

        rxWaveform = channelModel(txWaveform);

        [rcvdBitsLms, ~, ~] = demodLms(rxWaveform);
        [rcvdBitsRls, ~, ~] = demodRls(rxWaveform);
        [rcvdBitsKalman, ~, ~] = demodKalman(rxWaveform);

        rxBitsLms = viterbi_decode_rate12(rcvdBitsLms, next_state, out_bits);
        rxBitsRls = viterbi_decode_rate12(rcvdBitsRls, next_state, out_bits);
        rxBitsKalman = viterbi_decode_rate12(rcvdBitsKalman, next_state, out_bits);

        berLms(frameIdx) = mean(rcvdBitsLms ~= codedBits);
        berRls(frameIdx) = mean(rcvdBitsRls ~= codedBits);
        berKalman(frameIdx) = mean(rcvdBitsKalman ~= codedBits);
        berCodedLms(frameIdx) = mean(rxBitsLms ~= dataBits);
        berCodedRls(frameIdx) = mean(rxBitsRls ~= dataBits);
        berCodedKalman(frameIdx) = mean(rxBitsKalman ~= dataBits);
    end
    berLmsMean(sIdx) = mean(berLms);
    berRlsMean(sIdx) = mean(berRls);
    berKalmanMean(sIdx) = mean(berKalman);
    berCodedLmsMean(sIdx) = mean(berCodedLms);
    berCodedRlsMean(sIdx) = mean(berCodedRls);
    berCodedKalmanMean(sIdx) = mean(berCodedKalman);
end

figure;
semilogy(cfg.snrDb, berLmsMean,        '-o', 'LineWidth', 1.5, 'MarkerSize', 6); hold on;
semilogy(cfg.snrDb, berRlsMean,        '-s', 'LineWidth', 1.5, 'MarkerSize', 6);
semilogy(cfg.snrDb, berKalmanMean,     '-d', 'LineWidth', 1.5, 'MarkerSize', 6);
semilogy(cfg.snrDb, berCodedLmsMean,   '--^', 'LineWidth', 1.5, 'MarkerSize', 6);
semilogy(cfg.snrDb, berCodedRlsMean,   '--v', 'LineWidth', 1.5, 'MarkerSize', 6);
semilogy(cfg.snrDb, berCodedKalmanMean,'--x', 'LineWidth', 1.5, 'MarkerSize', 6);
xlabel('SNR (dB)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Bit Error Rate (BER)', 'FontSize', 12, 'FontWeight', 'bold');
title('BER vs SNR for LMS / RLS / Kalman (Coded & Uncoded)', 'FontSize', 13);
legend({'LMS', 'RLS', 'Kalman', 'Coded LMS', 'Coded RLS', 'Coded Kalman'}, ...
       'Location', 'southwest', 'FontSize', 10);
legend boxoff;
grid on;
grid minor;
