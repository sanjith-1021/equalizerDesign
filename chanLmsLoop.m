close all; clc; clear;

% HF modem loop with manual channel and fractionally spaced LMS equalizer (4 samples/symbol)
cfg = struct();
cfg.M = 8;
cfg.bitsPerSymbol = log2(cfg.M);
cfg.numDataSymbols = 32;
cfg.numPilotSymbols = 256;
cfg.numBlankSymbols = 16;
cfg.samplesPerSymbol = 4;
cfg.rolloff = 0.25;
cfg.filterSpan = 8;
cfg.snrDb = 30;
cfg.symbolRate = 2400;
cfg.sampleRate = cfg.symbolRate * cfg.samplesPerSymbol;
cfg.numFrames = 5;

cfg.numEqTaps = 5;
cfg.eqStep = 0.1;

rng(101);
cfg.chanNumTaps = 2;
chanCoeffs = randn(cfg.chanNumTaps, 1) + 1j * randn(cfg.chanNumTaps, 1);
chanCoeffs = chanCoeffs / norm(chanCoeffs);


rrcFilter = rcosdesign(cfg.rolloff, cfg.filterSpan, cfg.samplesPerSymbol, "sqrt");
filtDelay = cfg.filterSpan * cfg.samplesPerSymbol / 2;

frameSymbType = [zeros(cfg.numBlankSymbols,1); ones(cfg.numPilotSymbols,1); 2*ones(cfg.numDataSymbols,1); zeros(cfg.numBlankSymbols,1)];

blankSymbols = zeros(cfg.numBlankSymbols,1);

rng(101);
pilotBits = randi([0 1], cfg.numPilotSymbols * cfg.bitsPerSymbol, 1);
pilotInts = bi2de(reshape(pilotBits, cfg.bitsPerSymbol, []).', "left-msb");
pilotSymbols = pskmod(pilotInts, cfg.M, pi/cfg.M);

eqWeights = cell(cfg.numFrames, 1);
rxDataBits = cell(cfg.numFrames, 1);
berData = zeros(cfg.numFrames, 1);

totalSymbols = cfg.numPilotSymbols + cfg.numDataSymbols + 2*cfg.numBlankSymbols;
totalSamples = totalSymbols * cfg.samplesPerSymbol;

for frameIdx = 1:cfg.numFrames
    %% Tx Side
    dataBits = randi([0 1], cfg.numDataSymbols * cfg.bitsPerSymbol, 1);
    dataInts = bi2de(reshape(dataBits, cfg.bitsPerSymbol, []).', "left-msb");
    dataSymbols = pskmod(dataInts, cfg.M, pi/cfg.M);

    frameSymbols = [blankSymbols; pilotSymbols; dataSymbols; blankSymbols];

    upsampled = upsample(frameSymbols, cfg.samplesPerSymbol);

    txWaveform = conv(upsampled, rrcFilter, 'same');
 

    %% Channel 
    rxWaveform = awgn(txWaveform, cfg.snrDb, "measured");
    rxWaveform = conv(rxWaveform, chanCoeffs, "same");



    %% Rx Side
    matched = conv(rxWaveform, rrcFilter, 'same');

    eqTaps = zeros(cfg.numEqTaps, 1);
    mIdx =ceil(cfg.numEqTaps / 2); eqTaps(mIdx) = 1;
    DeLine = complex(zeros(cfg.numEqTaps, 1));
    eqOut = complex(zeros(totalSamples, 1));
    errHist = zeros(totalSamples);

    for dIdx = 1:totalSamples
        DeLine = [matched(dIdx); DeLine(1:end-1)];
        yi = eqTaps' * DeLine;

        yIdx = dIdx - mIdx + 1;
        symbIdx = floor((yIdx - 1)/4) + 1;
        pilotIdx = symbIdx - cfg.numBlankSymbols;

        if(yIdx > 0) && (mod(yIdx-1,4) == 0)
            if (frameSymbType(symbIdx) == 1)
                hdSymb = pilotSymbols(pilotIdx);
            elseif (frameSymbType(symbIdx) == 2)
                dInt = pskdemod(yi, cfg.M, pi/cfg.M);
                hdSymb = pskmod(dInt, cfg.M, pi/cfg.M);
            else 
                hdSymb = 0;
            end
            symbErr = hdSymb - yi;

            if (frameSymbType(symbIdx) > 0), eqTaps = eqTaps + cfg.eqStep * conj(symbErr) * DeLine; end

            errHist(yIdx) = symbErr;
        end
        if (yIdx > 0), eqOut(yIdx) = yi; end 
    end
    
    eqSymbols = eqOut(1:cfg.samplesPerSymbol:end);
    dataEq = eqSymbols(frameSymbType == 2);


    figure; plot(abs(errHist)); scatterplot(dataEq);


    rxInts = pskdemod(dataEq, cfg.M, pi/cfg.M);
    rxBits = de2bi(rxInts, cfg.bitsPerSymbol, "left-msb").';
    rxBits = rxBits(:);

    rxDataBits{frameIdx} = rxBits;
    berData(frameIdx) = mean(rxDataBits{frameIdx} ~= dataBits);
end

disp(berData);