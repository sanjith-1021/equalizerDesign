% Simple transceiver run for one slot with scatter, BER, EVM, and error history.
close all; clc; clear;

cfg = presets('basic-run');

snrDb = cfg.snrDb;
% nChanTaps = 5;
% chanCoeffs = randn(nChanTaps, 1) + 1j * randn(nChanTaps, 1);
% chanCoeffs = chanCoeffs / norm(chanCoeffs);
% channelModel = @(waveform) conv(awgn(waveform, snrDb, 'measured'), chanCoeffs, 'same');

% fdMax = 1;
% chanLM = stdchan('iturHFLM', cfg.sampRate, fdMax);
% chanLM.RandomStream = 'mt19937ar with seed';
% chanLM.Seed = 9999;
% chanLM.PathGainsOutputPort = false;
% channelModel = @(waveform) chanLM(awgn(waveform, snrDb, 'measured'));

doppSpec = doppler('Gaussian', (0.5*sqrt(2))/3);
chanLM = comm.RayleighChannel( ...
    'SampleRate', cfg.sampRate, ...
    'PathDelays', [0 1e-3], ...
    'AveragePathGains', [0 0], ...
    'NormalizePathGains', true, ...
    'FadingTechnique', 'Filtered Gaussian noise',...
    'MaximumDopplerShift', 3, ...
    'DopplerSpectrum', doppSpec, ...
    'RandomStream', 'mt19937ar with seed', ...
    'Seed', 9999 , ...
    'PathGainsOutputPort', false ...
    );
channelModel = @(waveform) chanLM(awgn(waveform, snrDb, 'measured'));

rng(cfg.seed);
trngBits01 = pnBits(cfg.nTrngSymbs01 * log2(cfg.M), [10 7], cfg.seed);
trngBits02 = randi([0 1], cfg.nTrngSymbs02 * log2(cfg.M), 1);
trngInts01 = bi2de(reshape(trngBits01, log2(cfg.M), []).', 'left-msb');
trngInts02 = bi2de(reshape(trngBits02, log2(cfg.M), []).', 'left-msb');
cfg.TrngSeq01 = pskmod(trngInts01, cfg.M, pi / cfg.M);
cfg.TrngSeq02 = pskmod(trngInts02, cfg.M, pi / cfg.M);


mod = TxModulator('Cfg', cfg);
algs = {'LMS', 'RLS', 'Kalman'};

dataBits = randi([0 1], cfg.nHops * cfg.nDataSymbs * log2(cfg.M), 1);
[txWaveform, txSymbs] = mod(dataBits);
rxWaveform = channelModel(txWaveform);
txDataSymbs = txSymbs(mod.SlotSymbType == 3);

figure('Name', 'Scatter Plots', 'NumberTitle', 'off');
figure('Name', 'Error History', 'NumberTitle', 'off');

figure(1);
subplot(2, 2, 1);
plot(real(txDataSymbs), imag(txDataSymbs), 'bo');
axis equal; grid on;
title('Tx Data Symbs');
xlabel('I'); ylabel('Q');

for algIdx = 1:numel(algs)
    dem = RxDemodulator('Cfg', cfg, 'EqualizerAlgorithm', algs{algIdx});
    [rxBits, eqSymbs, errHist] = dem(rxWaveform);

    ber = mean(rxBits ~= dataBits);
    rxDataSymbs = eqSymbs(dem.slotSymbType == 3);
    evmRms = sqrt(mean(abs(rxDataSymbs - txDataSymbs).^2) / mean(abs(txDataSymbs).^2));
    evmDb = 20 * log10(evmRms);

    fprintf('%s | BER: %.4e | EVM RMS: %.2f dB\n', algs{algIdx}, ber, evmDb);

    figure(1);
    subplot(2, 2, algIdx + 1);
    plot(real(rxDataSymbs), imag(rxDataSymbs), 'r.');
    axis equal; grid on;
    title(sprintf('Rx Equalized Data Symbs (%s)', algs{algIdx}));
    xlabel('I'); ylabel('Q');

    figure(2);
    subplot(numel(algs), 1, algIdx);
    plot(abs(errHist), 'k.-');
    grid on;
    title(sprintf('Error History (%s)', algs{algIdx}));
    xlabel('Symb Index');
    ylabel('|err|');
end
