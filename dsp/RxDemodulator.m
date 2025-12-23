classdef RxDemodulator < matlab.System
    %RXDEMODULATOR Matched filter, equalization, and PSK demodulation.

    properties
        Cfg
        EqualizerAlgorithm = 'LMS';
        PilotSymbols = [];
        FrameSymbolType = [];
        BitsPerSymbol
    end


    properties (Access = private)
        RrcFilter
        EqObj
    end

    methods
        function obj = RxDemodulator(varargin)
            if nargin > 0
                setProperties(obj, nargin, varargin{:});
            end
        end

        function val = get.BitsPerSymbol(obj)
            val = log2(obj.Cfg.M);
        end

    end

    methods (Access = protected)
        function setupImpl(obj, ~)
            c = obj.Cfg;
            obj.RrcFilter = rcosdesign(c.rolloff, c.filterSpan, c.samplesPerSymbol, 'sqrt');
            obj.FrameSymbolType = [ zeros(c.numBlankSymbols, 1);ones(c.numPilotSymbols, 1);2 * ones(c.numDataSymbols, 1); zeros(c.numBlankSymbols, 1) ];
            obj.PilotSymbols = c.PilotSymbols;

            obj.EqObj = ChannelEqualizer( ...
                'algorithm', obj.EqualizerAlgorithm, ...
                'nSampsPerSymb', c.samplesPerSymbol, ...
                'modOrder', c.M);
        end

        function [rxBits, eqSymbols, errHist] = stepImpl(obj, rxWaveform)
            matched = conv(rxWaveform, obj.RrcFilter, 'same');

            reset(obj.EqObj);
            [eqSymbols, errHist] = obj.EqObj(matched, obj.PilotSymbols, obj.FrameSymbolType);

            dataEq = eqSymbols(obj.FrameSymbolType == 2);
            rxInts = pskdemod(dataEq, obj.Cfg.M, pi / obj.Cfg.M);
            rxBits = de2bi(rxInts, obj.BitsPerSymbol, 'left-msb').';
            rxBits = rxBits(:);
        end

        function num = getNumInputsImpl(~)
            num = 1;
        end

        function num = getNumOutputsImpl(~)
            num = 3;
        end
    end
end
