% Simple transceiver run with scatter, BER, EVM, and error history.
close all; clc; clear;
addpath(genpath(pwd))

cfg = presets();
bitsPerSymb = log2(cfg.M);
nSlots = cfg.nSlots;

chan = channelModel(cfg, 1);

rng(cfg.seed);
trngBits1 = pnBits(cfg.nTrngSymbs01 * bitsPerSymb, [10 7], cfg.seed);
trngBits2 = randi([0 1], cfg.nTrngSymbs02 * bitsPerSymb, 1);
trngInts1 = bi2de(reshape(trngBits1, bitsPerSymb, []).', 'left-msb');
trngInts2 = bi2de(reshape(trngBits2, bitsPerSymb, []).', 'left-msb');
cfg.TrngSeq01 = pskmod(trngInts1, cfg.M, 0, 'gray');
cfg.TrngSeq02 = pskmod(trngInts2, cfg.M, 0, 'gray');


txMod = TxModulator('Cfg', cfg);

totalDataSymbs = nSlots * cfg.nHops * cfg.nDataSymbs;
[dataBits, codedBits] = fec_conv_encode(totalDataSymbs, cfg);

blockBits = cfg.nHops * cfg.nDataSymbs * bitsPerSymb;
interleaverDepth = blockBits / 40;
interleaverOrder = int_indx_gen(interleaverDepth);
interleavedBits = zeros(size(codedBits));
for slotIdx = 1:nSlots
    bitStart = (slotIdx - 1) * blockBits + 1;
    bitEnd = bitStart + blockBits - 1;
    block = codedBits(bitStart:bitEnd);
    interleavedBits(bitStart:bitEnd) = block(interleaverOrder);
end

slotSymbType = txMod.SlotSymbType;
slotSymbCount = numel(slotSymbType);
slotSampCount = slotSymbCount * cfg.sampsPerSymb;
dataSymbMask = slotSymbType == 3;
dataSymbCount = sum(dataSymbMask);
txWaveform = complex(zeros(slotSampCount * nSlots, 1));
txDataSymbs = complex(zeros(dataSymbCount * nSlots, 1));
for slotIdx = 1:nSlots
    bitStart = (slotIdx - 1) * blockBits + 1;
    bitEnd = bitStart + blockBits - 1;
    [slotWaveform, slotSymbs] = txMod(interleavedBits(bitStart:bitEnd));
    sampStart = (slotIdx - 1) * slotSampCount + 1;
    sampEnd = sampStart + slotSampCount - 1;
    dataStart = (slotIdx - 1) * dataSymbCount + 1;
    dataEnd = dataStart + dataSymbCount - 1;

    txWaveform(sampStart:sampEnd) = slotWaveform;
    txDataSymbs(dataStart:dataEnd) = slotSymbs(dataSymbMask);
end

[rxWaveform, ~] = chan(txWaveform);

scatterFig = figure('Name', 'Scatter Plots', 'NumberTitle', 'off');
errFig = figure('Name', 'Error History', 'NumberTitle', 'off');

figure(scatterFig);
subplot(1, 2, 1);
plot(real(txDataSymbs), imag(txDataSymbs), 'bo');
axis equal; grid on;
title('Tx Data Symbs');
xlabel('I'); ylabel('Q');

dem = RxDemodulator('Cfg', cfg);
[rxBits, eqSymbs, errHist] = dem(rxWaveform);

ber = mean(rxBits ~= dataBits);
rxDataSymbs = eqSymbs(dem.slotSymbType == 3);
evmRms = sqrt(mean(abs(rxDataSymbs - txDataSymbs).^2) / mean(abs(txDataSymbs).^2));
evmDb = 20 * log10(evmRms);

fprintf('Kalman | BER: %.4e | EVM RMS: %.2f dB\n', ber, evmDb);

figure(scatterFig);
subplot(1, 2, 2);
plot(real(rxDataSymbs), imag(rxDataSymbs), 'r.');
axis equal; grid on;
title('Rx Equalized Data Symbs (Kalman)');
xlabel('I'); ylabel('Q');

figure(errFig);
plot(abs(errHist), 'k.-');
grid on;
title('Error History (Kalman)');
xlabel('Symb Index');
ylabel('|err|');
