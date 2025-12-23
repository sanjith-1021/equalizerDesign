close all; clear; clc;
% HF modem loop with ITU low-latitude Rayleigh (stdchan) channel
M = 8;
bitsPerSymbol = log2(M);
numBlankSymbols = 32;
numDataSymbols = 32;
numPilotSymbols = 2048;
samplesPerSymbol = 4;
rolloff = 0.25;
filterSpan = 8;
snrDb = 100;
symbolRate = 2400;
sampleRate = symbolRate * samplesPerSymbol;
numFrames = 1;

rrcFilter = rcosdesign(rolloff, filterSpan, samplesPerSymbol, "sqrt");
filtDelay = filterSpan * samplesPerSymbol / 2;

rng(11);
pilotBits = randi([0 1], numPilotSymbols * bitsPerSymbol, 1);
pilotInts = bi2de(reshape(pilotBits, bitsPerSymbol, []).', "left-msb");
pilotSymbols = pskmod(pilotInts, M, pi/M);

blankSymbols = zeros(numBlankSymbols,1);

fd = 1;
chanLM = stdchan("iturHFLM", sampleRate, fd);
chanLM.RandomStream = "mt19937ar with seed";
chanLM.Seed = 9999;
chanLM.PathGainsOutputPort = true;
% chanLM.Visualization = "Impulse response";

dopplerSpec = doppler('Gaussian', .2);

fdMax = 10;
chan2 = comm.RayleighChannel( ...
    'SampleRate', sampleRate, ...
    'PathDelays', 0, ...
    'AveragePathGains', 0, ...
    'MaximumDopplerShift', fdMax, ...
    'DopplerSpectrum', dopplerSpec, ...
    'RandomStream', 'mt19937ar with seed', ...
    'Seed', 99, ...
    'PathGainsOutputPort', true);


chanDelaySamples = ceil(max(chanLM.PathDelays) * sampleRate);

rxPilotSymbols = cell(numFrames, 1);
rxDataSymbols = cell(numFrames, 1);
rxDataBits = cell(numFrames, 1);
rxPathGains = cell(numFrames, 1);
berData = zeros(numFrames, 1);

for frameIdx = 1:numFrames
    dataBits = randi([0 1], numDataSymbols * bitsPerSymbol, 1);
    dataInts = bi2de(reshape(dataBits, bitsPerSymbol, []).', "left-msb");
    dataSymbols = pskmod(dataInts, M, pi/M);

    frameSymbols = [blankSymbols;pilotSymbols; dataSymbols;blankSymbols];

    upsampled = upsample(frameSymbols, samplesPerSymbol);
    txWaveform = conv(upsampled, rrcFilter, "same");
    [rxWaveform, pathGains] = chanLM(txWaveform);
    % [rxWaveform, pathGains] = chan2(txWaveform);
    rxWaveform = awgn(rxWaveform, snrDb, "measured");

    matched = conv(rxWaveform, rrcFilter, "same");
    rxSymbols = matched(1 : samplesPerSymbol : end);

    rxInts = pskdemod(rxSymbols, M, pi/M);
    rxBits = de2bi(rxInts, bitsPerSymbol, "left-msb").';
    rxBits = rxBits(:);

    rxPilotSymbols{frameIdx} = rxSymbols(numBlankSymbols +1 :numPilotSymbols+ numBlankSymbols);
    rxDataSymbols{frameIdx} = rxSymbols(numBlankSymbols+numPilotSymbols + 1:end-numBlankSymbols);
end


freqErr =  rxPilotSymbols{1}.*(pilotSymbols.')';
scatterplot(freqErr);




