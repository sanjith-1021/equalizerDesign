% Basic loop for HF modem frame processing
M = 4;
bitsPerSymbol = log2(M);
numDataSymbols = 32;
numPilotSymbols = 16;
samplesPerSymbol = 4;
rolloff = 0.25;
filterSpan = 8;
snrDb = 10;
numFrames = 5;

rrcFilter = rcosdesign(rolloff, filterSpan, samplesPerSymbol, "sqrt");
filtDelay = filterSpan * samplesPerSymbol / 2;

rng(7);
pilotBits = randi([0 1], numPilotSymbols * bitsPerSymbol, 1);
pilotInts = bi2de(reshape(pilotBits, bitsPerSymbol, []).', "left-msb");
pilotSymbols = pskmod(pilotInts, M, pi/M);

rxPilotSymbols = cell(numFrames, 1);
rxDataSymbols = cell(numFrames, 1);
rxDataBits = cell(numFrames, 1);
berData = zeros(numFrames, 1);

for frameIdx = 1:numFrames
    dataBits = randi([0 1], numDataSymbols * bitsPerSymbol, 1);
    dataInts = bi2de(reshape(dataBits, bitsPerSymbol, []).', "left-msb");
    dataSymbols = pskmod(dataInts, M, pi/M);

    frameSymbols = [pilotSymbols; dataSymbols];

    upsampled = upsample(frameSymbols, samplesPerSymbol);
    txWaveform = conv(upsampled, rrcFilter, "same");
    rxWaveform = awgn(txWaveform, snrDb, "measured");

    matched = conv(rxWaveform, rrcFilter, "same");
    rxSymbols = matched(1 : samplesPerSymbol : end);

    rxInts = pskdemod(rxSymbols, M, pi/M);
    rxBits = de2bi(rxInts, bitsPerSymbol, "left-msb").';
    rxBits = rxBits(:);

    rxPilotSymbols{frameIdx} = rxSymbols(1:numPilotSymbols);
    rxDataSymbols{frameIdx} = rxSymbols(numPilotSymbols + 1:end);
    rxDataBits{frameIdx} = rxBits(numPilotSymbols * bitsPerSymbol + 1:end);
    berData(frameIdx) = mean(rxDataBits{frameIdx} ~= dataBits);
end

disp(berData);
