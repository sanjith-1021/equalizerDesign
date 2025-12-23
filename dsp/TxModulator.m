classdef TxModulator < matlab.System
    %TXMODULATOR Builds a framed PSK waveform with RRC pulse shaping.

    properties
        Cfg
        PilotSymbols = [];
    end

    properties (Dependent)
        BitsPerSymbol
        FrameSymbolType
    end

    properties (Access = private)
        RrcFilter
    end

    methods
        function obj = TxModulator(varargin)
            if nargin > 0
                setProperties(obj, nargin, varargin{:});
            end
        end

        function val = get.BitsPerSymbol(obj)
            val = log2(obj.Cfg.M);
        end

        function val = get.FrameSymbolType(obj)
            c = obj.Cfg;
            val = [zeros(c.numBlankSymbols, 1); ...
                   ones(c.numPilotSymbols, 1); ...
                   2 * ones(c.numDataSymbols, 1); ...
                   zeros(c.numBlankSymbols, 1)];
        end
    end

    methods (Access = protected)
        function setupImpl(obj, ~)
            c = obj.Cfg;
            obj.RrcFilter = rcosdesign(c.rolloff, c.filterSpan, c.samplesPerSymbol, 'sqrt');
            obj.PilotSymbols = c.PilotSymbols;
        end

        function [txWaveform, frameSymbols] = stepImpl(obj, dataBits)
            c = obj.Cfg;
            dataBits = dataBits(:);
            dataInts = bi2de(reshape(dataBits, obj.BitsPerSymbol, []).', 'left-msb');
            dataSymbols = pskmod(dataInts, c.M, pi / c.M);

            frameSymbols = [zeros(c.numBlankSymbols, 1); obj.PilotSymbols; dataSymbols; zeros(c.numBlankSymbols, 1)];
            upsampled = upsample(frameSymbols, c.samplesPerSymbol);
            txWaveform = conv(upsampled, obj.RrcFilter, 'same');
        end

        function num = getNumOutputsImpl(~)
            num = 2;
        end
    end
end
