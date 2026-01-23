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

        function [txWaveform, frameSymbs, dataBits] = stepImpl(obj)
            c = obj.Cfg;
            nSlots = 1;
            if isfield(c, 'nSlots') && ~isempty(c.nSlots)
                nSlots = c.nSlots;
            end
            totalDataSymbs = nSlots * c.nHops * c.nDataSymbs;
            [dataBits, codedBits] = fec_conv_encode(totalDataSymbs, c);
            dataBits = dataBits(:);

            blockBits = c.nHops * c.nDataSymbs * obj.BitsPerSymb;
            if mod(blockBits, 40) ~= 0
                error('Interleaver block size must be divisible by 40.');
            end
            interleaverDepth = blockBits / 40;
            interleaverOrder = int_indx_gen(interleaverDepth);
            interleavedBits = zeros(size(codedBits));
            for slotIdx = 1:nSlots
                bitStart = (slotIdx - 1) * blockBits + 1;
                bitEnd = bitStart + blockBits - 1;
                block = codedBits(bitStart:bitEnd);
                interleavedBits(bitStart:bitEnd) = block(interleaverOrder);
            end

            dataInts = bi2de(reshape(interleavedBits, obj.BitsPerSymb, []).', 'left-msb');
            dataSymbs = pskmod(dataInts, c.M, 0, 'gray');

            nDataSymbs = c.nHops * c.nDataSymbs;
            if mod(numel(dataSymbs), nDataSymbs) ~= 0
                error('Data symbols must be a multiple of %d.', nDataSymbs);
            end
            if numel(dataSymbs) ~= nDataSymbs * nSlots
                error('Expected %d data symbols, got %d.', nDataSymbs * nSlots, numel(dataSymbs));
            end
            dataSymbs = reshape(dataSymbs, c.nDataSymbs, c.nHops * nSlots);

            slotSymbCount = numel(obj.SlotSymbType);
            frameSymbs = complex(zeros(slotSymbCount * nSlots, 1));
            for slotIdx = 1:nSlots
                hopStart = (slotIdx - 1) * c.nHops + 1;
                hopEnd = hopStart + c.nHops - 1;
                slotData = dataSymbs(:, hopStart:hopEnd);
                hopSymbs = [repmat(obj.TrngSeq02, 1, c.nHops); slotData];
                slotFrame = [zeros(c.nBlankSymbs, 1); ...
                             obj.TrngSeq01; ...
                             hopSymbs(:); ...
                             zeros(c.nBlankSymbs, 1)];
                symbStart = (slotIdx - 1) * slotSymbCount + 1;
                symbEnd = symbStart + slotSymbCount - 1;
                frameSymbs(symbStart:symbEnd) = slotFrame;
            end

            upsampled = upsample(frameSymbs, c.sampsPerSymb);
            txWaveform = conv(upsampled, obj.RrcFilter, 'same');
        end

        function num = getNumOutputsImpl(~)
            num = 3;
        end
    end
end
