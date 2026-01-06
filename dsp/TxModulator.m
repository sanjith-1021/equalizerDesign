classdef TxModulator < matlab.System
    %TXMODULATOR Builds a framed PSK waveform with RRC pulse shaping.

    properties
        Cfg
        TrngSeq01 = [];
        TrngSeq02 = [];
    end

    properties (Dependent)
        BitsPerSymb
        SlotSymbType
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

        function val = get.BitsPerSymb(obj)
            val = log2(obj.Cfg.M);
        end

        function val = get.SlotSymbType(obj)
            c = obj.Cfg;
            hopTypes = repmat([2 * ones(c.nTrngSymbs02, 1); ...
                               3 * ones(c.nDataSymbs, 1)], c.nHops, 1);
            val = [zeros(c.nBlankSymbs, 1); ...
                   ones(c.nTrngSymbs01, 1); ...
                   hopTypes; ...
                   zeros(c.nBlankSymbs, 1)];
        end
    end

    methods (Access = protected)
        function setupImpl(obj, ~)
            c = obj.Cfg;
            obj.RrcFilter = rcosdesign(c.rolloff, c.filterSpan, c.sampsPerSymb, 'sqrt');
            obj.TrngSeq01 = c.TrngSeq01;
            obj.TrngSeq02 = c.TrngSeq02;
        end

        function [txWaveform, frameSymbs] = stepImpl(obj, dataBits)
            c = obj.Cfg;
            dataBits = dataBits(:);
            dataInts = bi2de(reshape(dataBits, obj.BitsPerSymb, []).', 'left-msb');
            dataSymbs = pskmod(dataInts, c.M, pi / c.M);

            nDataSymbs = c.nHops * c.nDataSymbs;
            if numel(dataSymbs) ~= nDataSymbs
                error('Expected %d data symbols, got %d.', nDataSymbs, numel(dataSymbs));
            end
            dataSymbs = reshape(dataSymbs, c.nDataSymbs, c.nHops);
            hopSymbs = [repmat(obj.TrngSeq02, 1, c.nHops); dataSymbs];
            frameSymbs = [zeros(c.nBlankSymbs, 1); ...
                            obj.TrngSeq01; ...
                            hopSymbs(:); ...
                            zeros(c.nBlankSymbs, 1)];
            upsampled = upsample(frameSymbs, c.sampsPerSymb);
            txWaveform = conv(upsampled, obj.RrcFilter, 'same');
        end

        function num = getNumOutputsImpl(~)
            num = 2;
        end
    end
end
