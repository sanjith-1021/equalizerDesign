% Simple transceiver run for one slot with scatter, BER, EVM, and error history.
close all; clc; clear;

cfg = presets('intrlvd-run');

chan = channelModel(cfg, 2);

rng(cfg.seed);
trngBits01 = pnBits(cfg.nTrngSymbs01 * log2(cfg.M), [10 7], cfg.seed);
trngBits02 = randi([0 1], cfg.nTrngSymbs02 * log2(cfg.M), 1);
trngInts01 = bi2de(reshape(trngBits01, log2(cfg.M), []).', 'left-msb');
trngInts02 = bi2de(reshape(trngBits02, log2(cfg.M), []).', 'left-msb');
cfg.TrngSeq01 = pskmod(trngInts01, cfg.M, 0, 'gray');
cfg.TrngSeq02 = pskmod(trngInts02, cfg.M, 0, 'gray');


mod = TxModulator('Cfg', cfg);
algs = {'Kalman'};      % 

% intlinprog_depth = 72;
intlinprog_depth = 576;
intl_read_order = int_indx_gen(intlinprog_depth);% 72 or 576 , default depth is 40.
deintl_read_order(intl_read_order) = 1:numel(intl_read_order);

[dataBits, coded_bits, ~, data_sym_idx] = fec_conv_encode(cfg.nDataSymbs*cfg.nHops, cfg);
intlvdBits = coded_bits(intl_read_order);
[txWaveform, txSymbs] = mod(intlvdBits);
[rxWaveform, chanObj] = chan(txWaveform);
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
