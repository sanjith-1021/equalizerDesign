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
        slotSymbTypeBase
        knownSymbsBase
        nSlots
        slotSymbCount
        slotSampCount
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
            obj.nSlots = c.nSlots;
            obj.RrcFilter = rcosdesign(c.rolloff, c.filterSpan, c.sampsPerSymb, 'sqrt');
            hopTypes = repmat([2 * ones(c.nTrngSymbs02, 1); ...
                               3 * ones(c.nDataSymbs, 1)], c.nHops, 1);
            obj.slotSymbType = [zeros(c.nBlankSymbs, 1); ...
                               ones(c.nTrngSymbs01, 1); ...
                               hopTypes; ...
                               zeros(c.nBlankSymbs, 1)];
            knownSymbs = [c.TrngSeq01; repmat(c.TrngSeq02, c.nHops, 1)];
            obj.KnownSymbs = quantizeComplex12(knownSymbs);
            obj.slotSymbCount = numel(obj.slotSymbType);
            obj.slotSampCount = obj.slotSymbCount * c.sampsPerSymb;

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
            matched = quantizeComplex12(matched);
            totalSymbs = obj.slotSymbCount * obj.nSlots;
            eqSymbs = complex(zeros(totalSymbs, 1));
            errHist = complex(zeros(totalSymbs, 1));

            blockBits = obj.Cfg.nHops * obj.Cfg.nDataSymbs * obj.BitsPerSymb;
            interleaverDepth = blockBits / 40;
            interleaverOrder = int_indx_gen(interleaverDepth);
            deinterleaverOrder = zeros(size(interleaverOrder));
            deinterleaverOrder(interleaverOrder) = 1:numel(interleaverOrder);
            deintlvdAll = zeros(blockBits * obj.nSlots, 1);

            EsNo = 10^(30/10);
            nVar = 1 / EsNo;

            for slotIdx = 1:obj.nSlots
                sampStart = (slotIdx - 1) * obj.slotSampCount + 1;
                sampEnd = sampStart + obj.slotSampCount - 1;
                symbStart = (slotIdx - 1) * obj.slotSymbCount + 1;
                symbEnd = symbStart + obj.slotSymbCount - 1;

                slotMatched = matched(sampStart:sampEnd);
                % reset(obj.EqObj);
                [slotEqSymbs, slotErr] = obj.EqObj(slotMatched, obj.KnownSymbs, obj.slotSymbType);
                eqSymbs(symbStart:symbEnd) = slotEqSymbs;
                errHist(symbStart:symbEnd) = slotErr;

                dataEq = slotEqSymbs(obj.slotSymbType == 3);
                dataEq = double(dataEq) / double(2^11 - 1);
                llr = pskdemod(dataEq, obj.Cfg.M, 0, 'gray', ...
                    'OutputType','llr', 'NoiseVariance', nVar);
                L = llr(:);
                Lclip = max(min(L, obj.Cfg.fec.llrClip), -obj.Cfg.fec.llrClip);
                p1 = 1 ./ (1 + exp(Lclip));
                demap_out = round(p1 * (2^obj.Cfg.fec.nSoft - 1));

                deintlvd_out = demap_out(deinterleaverOrder);
                bitStart = (slotIdx - 1) * blockBits + 1;
                bitEnd = bitStart + blockBits - 1;
                deintlvdAll(bitStart:bitEnd) = deintlvd_out;
            end

            rxBits = vitdec(deintlvdAll, obj.Cfg.fec.trellis, ...
                obj.Cfg.fec.tblen, 'trunc', 'soft', obj.Cfg.fec.nSoft);
            rxBits = rxBits(:);
        %     eqSymbs = matched(1:obj.Cfg.sampsPerSymb:end);
        %    errHist = zeros(length(eqSymbs),1);
        end

        function num = getNumInputsImpl(~)
            num = 1;
        end

        function num = getNumOutputsImpl(~)
            num = 3;
        end
    end
end
