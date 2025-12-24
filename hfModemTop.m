% HF modem end-to-end simulation using modular System objects.
close all; clc; clear;

cfg.M = 8;
cfg.numDataSymbols =16;
cfg.numPilotSymbols = 16;
cfg.numBlankSymbols = 16;
cfg.samplesPerSymbol = 4;
cfg.rolloff = 0.25;
cfg.filterSpan = 8;
cfg.snrDb = 0:4:30;
cfg.sampRate = 9600;
cfg.numFrames = 100;


cfg.fMax = 1;
cfg.seed = 101;

generators = [1 1 1 1 0 0 1;  % G0 = 171₈
              1 0 1 1 0 1 1];  % G1 = 133₈

[next_state, out_bits] = build_trellis_rate12(generators);


rng(cfg.seed);
pilotBits = randi([0 1], cfg.numPilotSymbols * log2(cfg.M), 1);
pilotInts = bi2de(reshape(pilotBits, log2(cfg.M), []).', 'left-msb');
cfg.PilotSymbols = pskmod(pilotInts, cfg.M, pi / cfg.M);

% nChanTaps = 4;
% chanCoeffs = randn(nChanTaps, 1) + 1j * randn(nChanTaps, 1);
% chanCoeffs = chanCoeffs / norm(chanCoeffs);
% channelModel = @(waveform, snrDb) conv(awgn(waveform, snrDb, 'measured'), chanCoeffs, 'same');

chanLM = stdchan("iturHFLM", cfg.sampRate, cfg.fMax);
chanLM.RandomStream = "mt19937ar with seed";
chanLM.Seed = 9999;
chanLM.PathGainsOutputPort = false;
channelModel = @(waveform, snrDb) chanLM(awgn(waveform, snrDb, 'measured'));

modulator = TxModulator('Cfg', cfg);

demodLms = RxDemodulator('Cfg', cfg, 'EqualizerAlgorithm', 'LMS');
demodRls = RxDemodulator('Cfg', cfg, 'EqualizerAlgorithm', 'RLS');
demodKalman = RxDemodulator('Cfg', cfg, 'EqualizerAlgorithm', 'Kalman');

snrPoints = cfg.snrDb(:);
berLms = zeros(numel(snrPoints), cfg.numFrames);
berRls = zeros(numel(snrPoints), cfg.numFrames);
berKalman = zeros(numel(snrPoints), cfg.numFrames);
berLmsCoded = zeros(numel(snrPoints), cfg.numFrames);
berRlsCoded = zeros(numel(snrPoints), cfg.numFrames);

for sIdx = 1:numel(snrPoints)
    snrDb = snrPoints(sIdx);
    parfor frameIdx = 1:cfg.numFrames
        dataBits = randi([0 1], cfg.numDataSymbols * log2(cfg.M)/2, 1);
        codedBits = conv_encode_rate12(dataBits, generators);
        [txWaveform, frameSymbols] = modulator(codedBits);

        rxWaveform = channelModel(txWaveform, snrDb);

        [rcvdBitsLms, ~, errHistLms] = demodLms(rxWaveform);
        [rcvdBitsRls, ~, errHistRls] = demodRls(rxWaveform);
        [rcvdBitsKalman, ~, errHistKalman] = demodKalman(rxWaveform);

        rxBitsLms = viterbi_decode_rate12(rcvdBitsLms, next_state, out_bits);
        rxBitsRls = viterbi_decode_rate12(rcvdBitsRls, next_state, out_bits);
        rxBitsKalman = viterbi_decode_rate12(rcvdBitsKalman, next_state, out_bits);

        berLms(sIdx, frameIdx) = mean(rcvdBitsLms ~= codedBits);
        berRls(sIdx, frameIdx) = mean(rcvdBitsRls ~= codedBits);
        berKalman(sIdx, frameIdx) = mean(rcvdBitsKalman ~= codedBits);
        berCodedLms(sIdx, frameIdx) = mean(rxBitsLms ~= dataBits);
        berCodedRls(sIdx, frameIdx) = mean(rxBitsRls ~= dataBits);
        berCodedKalman(sIdx, frameIdx) = mean(rxBitsKalman ~= dataBits);
    end
end

berLmsMean = mean(berLms, 2);
berRlsMean = mean(berRls, 2);
berKalmanMean = mean(berKalman, 2);
berCodedLmsMean = mean(berCodedLms,2);
berCodedRlsMean = mean(berCodedRls,2);
berCodedKalmanMean = mean(berCodedKalman,2);

figure;
semilogy(cfg.snrDb, berLmsMean,       '-o', 'LineWidth', 1.5, 'MarkerSize', 6); hold on;
semilogy(cfg.snrDb, berRlsMean,       '-s', 'LineWidth', 1.5, 'MarkerSize', 6);
semilogy(cfg.snrDb, berKalmanMean,    '-d', 'LineWidth', 1.5, 'MarkerSize', 6);
semilogy(cfg.snrDb, berCodedLmsMean,  '--^','LineWidth', 1.5, 'MarkerSize', 6);
semilogy(cfg.snrDb, berCodedRlsMean,  '--v','LineWidth', 1.5, 'MarkerSize', 6);
semilogy(cfg.snrDb, berCodedKalmanMean,'--x','LineWidth', 1.5, 'MarkerSize', 6);
xlabel('SNR (dB)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Bit Error Rate (BER)', 'FontSize', 12, 'FontWeight', 'bold');
title('BER vs SNR for LMS / RLS / Kalman (Coded & Uncoded)', 'FontSize', 13);
legend({'LMS', 'RLS', 'Kalman', 'Coded LMS', 'Coded RLS', 'Coded Kalman'}, ...
       'Location', 'southwest', 'FontSize', 10);
legend boxoff;
grid on;
grid minor;