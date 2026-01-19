classdef RxDemodulator < matlab.System
    %RXDEMODULATOR Matched filter, equalization, and PSK demodulation.

    properties
        Cfg
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
            nSlots = 1;
            if isfield(c, 'nSlots') && ~isempty(c.nSlots)
                nSlots = c.nSlots;
            end
            obj.RrcFilter = rcosdesign(c.rolloff, c.filterSpan, c.sampsPerSymb, 'sqrt');
            hopTypes = repmat([2 * ones(c.nTrngSymbs02, 1); ...
                               3 * ones(c.nDataSymbs, 1)], c.nHops, 1);
            slotSymbType = [zeros(c.nBlankSymbs, 1); ...
                               ones(c.nTrngSymbs01, 1); ...
                               hopTypes; ...
                               zeros(c.nBlankSymbs, 1)];
            knownSymbs = [c.TrngSeq01; repmat(c.TrngSeq02, c.nHops, 1)];
            obj.slotSymbType = repmat(slotSymbType, nSlots, 1);
            obj.KnownSymbs = repmat(knownSymbs, nSlots, 1);

            eqCfg = c.eq;
            eqArgs = {'nSampsPerSymb', c.sampsPerSymb, ...
                      'modOrder', c.M, ...
                      'nFeedforwardTaps', eqCfg.nFeedforwardTaps, ...
                      'nFeedbackTaps', eqCfg.nFeedbackTaps, ...
                      'trainEq', eqCfg.train, ...
                      'dataEq', eqCfg.data};
            obj.EqObj = ChannelEqualizer(eqArgs{:});
        end

        function [rxBits, eqSymbs, errHist] = stepImpl(obj, rxWaveform)
            matched = conv(rxWaveform, obj.RrcFilter, 'same');

            reset(obj.EqObj);
            [eqSymbs, errHist] = obj.EqObj(matched, obj.KnownSymbs, obj.slotSymbType);
        %     eqSymbs = matched(1:obj.Cfg.sampsPerSymb:end);
        %    errHist = zeros(length(eqSymbs),1);
            
            dataEq = eqSymbs(obj.slotSymbType == 3);


            EsNo = 10^(30/10);
            nVar = 1/EsNo;
            llr = pskdemod(dataEq, obj.Cfg.M, 0, 'gray', ...
                'OutputType','llr', 'NoiseVariance', nVar);
            L = llr(:);
            Lclip = max(min(L, obj.Cfg.fec.llrClip), -obj.Cfg.fec.llrClip);
            p1 = 1 ./ (1 + exp(Lclip));
            demap_out = round(p1 * (2^obj.Cfg.fec.nSoft - 1));

            blockBits = obj.Cfg.nHops * obj.Cfg.nDataSymbs * obj.BitsPerSymb;
            if mod(blockBits, 40) ~= 0
                error('Interleaver block size must be divisible by 40.');
            end
            interleaverDepth = blockBits / 40;
            interleaverOrder = int_indx_gen(interleaverDepth);
            deinterleaverOrder = zeros(size(interleaverOrder));
            deinterleaverOrder(interleaverOrder) = 1:numel(interleaverOrder);

            if mod(numel(demap_out), blockBits) ~= 0
                error('Demapper output length does not align with interleaver block size.');
            end
            nBlocks = numel(demap_out) / blockBits;
            decodedLen = floor(blockBits * obj.Cfg.fec.rate);
            rxBits = zeros(nBlocks * decodedLen, 1);
            for blkIdx = 1:nBlocks
                blkStart = (blkIdx - 1) * blockBits + 1;
                blkEnd = blkStart + blockBits - 1;
                block = demap_out(blkStart:blkEnd);
                deintlvd_out = block(deinterleaverOrder);
                dec_bits = vitdec(deintlvd_out, obj.Cfg.fec.trellis, ...
                    obj.Cfg.fec.tblen, 'trunc', 'soft', obj.Cfg.fec.nSoft);
                outStart = (blkIdx - 1) * decodedLen + 1;
                outEnd = outStart + decodedLen - 1;
                rxBits(outStart:outEnd) = dec_bits(:);
            end
        end

        function num = getNumInputsImpl(~)
            num = 1;
        end

        function num = getNumOutputsImpl(~)
            num = 3;
        end
    end
end
