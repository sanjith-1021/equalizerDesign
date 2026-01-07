% HF modem end-to-end BER sweep (no animation).
close all; clc; clear;

scriptDir = fileparts(mfilename('fullpath'));
addpath(fullfile(scriptDir, 'config'));
cfg = presets('ber-sweep');

[next_state, out_bits] = build_trellis_rate12(cfg.generators);

rng(cfg.seed);
trngBits01 = randi([0 1], cfg.nTrngSymbs01 * log2(cfg.M), 1);
trngBits02 = randi([0 1], cfg.nTrngSymbs02 * log2(cfg.M), 1);
trngInts01 = bi2de(reshape(trngBits01, log2(cfg.M), []).', 'left-msb');
trngInts02 = bi2de(reshape(trngBits02, log2(cfg.M), []).', 'left-msb');
cfg.TrngSeq01 = pskmod(trngInts01, cfg.M, pi / cfg.M);
cfg.TrngSeq02 = pskmod(trngInts02, cfg.M, pi / cfg.M);

% nChanTaps = 4;
% chanCoeffs = randn(nChanTaps, 1) + 1j * randn(nChanTaps, 1);
% chanCoeffs = chanCoeffs / norm(chanCoeffs);
% channelModel = @(waveform, snrDb) conv(awgn(waveform, snrDb, 'measured'), chanCoeffs, 'same');

% chanLM = stdchan("iturHFLM", cfg.sampRate, cfg.fMax);
% chanLM.RandomStream = "mt19937ar with seed";
% chanLM = stdchan('iturHFLM', cfg.sampRate, cfg.fMax);
% chanLM.RandomStream = 'mt19937ar with seed';
% chanLM.Seed = 9999;
% chanLM.PathGainsOutputPort = false;
% channelModel = @(waveform, snrDb) chanLM(awgn(waveform, snrDb, 'measured'));

sigmaNorm = (0.5*sqrt(2))/3;
doppSpec = doppler('Gaussian', sigmaNorm);
chanLM = comm.RayleighChannel( ...
    'SampleRate', cfg.sampRate, ...
    'PathDelays', [0 1e-3], ...
    'AveragePathGains', [0 0], ...
    'NormalizePathGains', true, ...
    'FadingTechnique', 'Filtered Gaussian noise',...
    'MaximumDopplerShift', cfg.fMax, ...
    'DopplerSpectrum', doppSpec, ...
    'RandomStream', 'mt19937ar with seed', ...
    'Seed', 9999 , ...
    'PathGainsOutputPort', false ...
    );
channelModel = @(waveform, snrDb) awgn(chanLM(waveform), snrDb, 'measured');


snrPoints = cfg.snrDb(:);
nSlots = cfg.nSlots;
rndSeed = cfg.seed;

berLms = zeros(numel(snrPoints),nSlots);
berRls = zeros(numel(snrPoints),nSlots);
berKalman = zeros(numel(snrPoints), nSlots);
berCodedLms = zeros(numel(snrPoints), nSlots);
berCodedRls = zeros(numel(snrPoints), nSlots);
berCodedKalman = zeros(numel(snrPoints), nSlots);

for snrIdx = 1:numel(snrPoints)
    snrDb = snrPoints(snrIdx);

    parfor slotIdx = 1:nSlots
        
        modulator = TxModulator('Cfg', cfg);
        demodLms = RxDemodulator('Cfg', cfg, 'EqualizerAlgorithm', 'LMS');
        demodRls = RxDemodulator('Cfg', cfg, 'EqualizerAlgorithm', 'RLS');
        demodKalman = RxDemodulator('Cfg', cfg, 'EqualizerAlgorithm', 'Kalman');

        dataBits = randi([0 1], cfg.nHops * cfg.nDataSymbs * log2(cfg.M)/2, 1);
        codedBits = conv_encode_rate12(dataBits, cfg.generators);
        txWaveform = modulator(codedBits);

        rxWaveform = channelModel(txWaveform, snrDb);

        [rcvdBitsLms, ~, ~] = demodLms(rxWaveform);
        [rcvdBitsRls, ~, ~] = demodRls(rxWaveform);
        [rcvdBitsKalman, ~, ~] = demodKalman(rxWaveform);

        rxBitsLms = viterbi_decode_rate12(rcvdBitsLms, next_state, out_bits);
        rxBitsRls = viterbi_decode_rate12(rcvdBitsRls, next_state, out_bits);
        rxBitsKalman = viterbi_decode_rate12(rcvdBitsKalman, next_state, out_bits);

        berLms(snrIdx, slotIdx) = mean(rcvdBitsLms ~= codedBits);
        berRls(snrIdx, slotIdx) = mean(rcvdBitsRls ~= codedBits);
        berKalman(snrIdx, slotIdx) = mean(rcvdBitsKalman ~= codedBits);
        berCodedLms(snrIdx, slotIdx) = mean(rxBitsLms ~= dataBits);
        berCodedRls(snrIdx, slotIdx) = mean(rxBitsRls ~= dataBits);
        berCodedKalman(snrIdx, slotIdx) = mean(rxBitsKalman ~= dataBits);
    end
end

berLmsMean = mean(berLms, 2);
berRlsMean = mean(berRls, 2);
berKalmanMean = mean(berKalman, 2);
berCodedLmsMean = mean(berCodedLms, 2);
berCodedRlsMean = mean(berCodedRls, 2);
berCodedKalmanMean = mean(berCodedKalman, 2);

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
