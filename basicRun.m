% Simple transceiver run for one slot with scatter, BER, EVM, and error history.
close all; clc; clear;

scriptDir = fileparts(mfilename('fullpath'));
addpath(fullfile(scriptDir, 'config'));
cfg = presets('basic-run');

snrDb = cfg.snrDb;
% nChanTaps = 5;
% chanCoeffs = randn(nChanTaps, 1) + 1j * randn(nChanTaps, 1);
% chanCoeffs = chanCoeffs / norm(chanCoeffs);
% channelModel = @(waveform) conv(awgn(waveform, snrDb, 'measured'), chanCoeffs, 'same');

chanLM = stdchan("iturHFLM", cfg.sampRate, cfg.fMax);
chanLM.RandomStream = "mt19937ar with seed";
chanLM = stdchan('iturHFLM', cfg.sampRate, cfg.fMax);
chanLM.RandomStream = 'mt19937ar with seed';
chanLM.Seed = 9999;
chanLM.PathGainsOutputPort = false;
channelModel = @(waveform) chanLM(awgn(waveform, snrDb, 'measured'));


rng(cfg.seed);
trngBits01 = randi([0 1], cfg.nTrngSymbs01 * log2(cfg.M), 1);
trngBits02 = randi([0 1], cfg.nTrngSymbs02 * log2(cfg.M), 1);
trngInts01 = bi2de(reshape(trngBits01, log2(cfg.M), []).', 'left-msb');
trngInts02 = bi2de(reshape(trngBits02, log2(cfg.M), []).', 'left-msb');
cfg.TrngSeq01 = pskmod(trngInts01, cfg.M, pi / cfg.M);
cfg.TrngSeq02 = pskmod(trngInts02, cfg.M, pi / cfg.M);


mod = TxModulator('Cfg', cfg);
dem = RxDemodulator('Cfg', cfg, 'EqualizerAlgorithm', 'Kalman');

dataBits = randi([0 1], cfg.nHops * cfg.nDataSymbs * log2(cfg.M), 1);
[txWaveform, txSymbs] = mod(dataBits);
rxWaveform = channelModel(txWaveform);
[rxBits, eqSymbs, errHist] = dem(rxWaveform);

ber = mean(rxBits ~= dataBits);

txDataSymbs = txSymbs(mod.SlotSymbType == 3);
rxDataSymbs = eqSymbs(dem.slotSymbType == 3);
evmRms = sqrt(mean(abs(rxDataSymbs - txDataSymbs).^2) / mean(abs(txDataSymbs).^2));
evmPct = 100 * evmRms;

fprintf('BER: %.4e\n', ber);
fprintf('EVM RMS: %.2f %%\n', evmPct);

figure('Name', 'Scatter Plots', 'NumberTitle', 'off');
subplot(1, 2, 1);
plot(real(txDataSymbs), imag(txDataSymbs), 'bo');
axis equal; grid on;
title('Tx Data Symbs');
xlabel('I'); ylabel('Q');

subplot(1, 2, 2);
plot(real(rxDataSymbs), imag(rxDataSymbs), 'r.');
axis equal; grid on;
title('Rx Equalized Data Symbs');
xlabel('I'); ylabel('Q');

figure('Name', 'Error History', 'NumberTitle', 'off');
plot(abs(errHist), 'k.-');
grid on;
title('Equalizer Error History');
xlabel('Symb Index');
ylabel('|err|');
