classdef RxDemodulator < matlab.System
    %RXDEMODULATOR Matched filter, equalization, and PSK demodulation.

    properties
        Cfg
        EqualizerAlgorithm = 'LMS';
        KnownSymbs = [];
        slotSymbType = [];
        BitsPerSymb
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

        function val = get.BitsPerSymb(obj)
            val = log2(obj.Cfg.M);
        end

    end

    methods (Access = protected)
        function setupImpl(obj, ~)
            c = obj.Cfg;
            obj.RrcFilter = rcosdesign(c.rolloff, c.filterSpan, c.sampsPerSymb, 'sqrt');
            hopTypes = repmat([2 * ones(c.nTrngSymbs02, 1); ...
                               3 * ones(c.nDataSymbs, 1)], c.nHops, 1);
            obj.slotSymbType = [zeros(c.nBlankSymbs, 1); ...
                                   ones(c.nTrngSymbs01, 1); ...
                                   hopTypes; ...
                                   zeros(c.nBlankSymbs, 1)];
            obj.KnownSymbs = [c.TrngSeq01; repmat(c.TrngSeq02, c.nHops, 1)];

            obj.EqObj = ChannelEqualizer( ...
                'algorithm', obj.EqualizerAlgorithm, ...
                'nSampsPerSymb', c.sampsPerSymb, ...
                'modOrder', c.M);
        end

        function [rxBits, eqSymbs, errHist] = stepImpl(obj, rxWaveform)
            matched = conv(rxWaveform, obj.RrcFilter, 'same');

            reset(obj.EqObj);
            [eqSymbs, errHist] = obj.EqObj(matched, obj.KnownSymbs, obj.slotSymbType);
            
            dataEq = eqSymbs(obj.slotSymbType == 3);
            rxInts = pskdemod(dataEq, obj.Cfg.M, pi / obj.Cfg.M);
            rxBits = de2bi(rxInts, obj.BitsPerSymb, 'left-msb').';
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
