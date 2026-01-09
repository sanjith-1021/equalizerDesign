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

            eqArgs = {'algorithm', obj.EqualizerAlgorithm, ...
                      'nSampsPerSymb', c.sampsPerSymb, ...
                      'modOrder', c.M};
            eqArgs = [eqArgs, obj.eqStructToArgs(c.eq)]; 
            obj.EqObj = ChannelEqualizer(eqArgs{:});
        end

        function [rxBits, eqSymbs, errHist] = stepImpl(obj, rxWaveform)
            matched = conv(rxWaveform, obj.RrcFilter, 'same');

            reset(obj.EqObj);
            [eqSymbs, errHist] = obj.EqObj(matched, obj.KnownSymbs, obj.slotSymbType);
            % eqSymbs = matched(1:obj.Cfg.sampsPerSymb:end);
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


            intlinprog_depth = 576;
            intl_read_order = int_indx_gen(intlinprog_depth);% 72 or 576 , default depth is 40.
            deintl_read_order(intl_read_order) = 1:numel(intl_read_order);
            dintlvd_out = demap_out(deintl_read_order);
            dec_bits = vitdec(dintlvd_out, obj.Cfg.fec.trellis, obj.Cfg.fec.tblen, 'trunc', 'soft', 3);
            rxBits = dec_bits(:);
        end

        function num = getNumInputsImpl(~)
            num = 1;
        end

        function num = getNumOutputsImpl(~)
            num = 3;
        end
    end

    methods (Access = private)
        function args = eqStructToArgs(~, eqCfg)
            names = fieldnames(eqCfg);
            args = {};
            for k = 1:numel(names)
                name = names{k};
                if any(strcmpi(name, {'algorithm', 'nSampsPerSymb', 'modOrder'}))
                    continue;
                end
                if strcmpi(name, 'train')
                    args = [args, {'trainEq', eqCfg.train}]; %#ok<AGROW>
                    continue;
                end
                if strcmpi(name, 'data')
                    args = [args, {'dataEq', eqCfg.data}]; %#ok<AGROW>
                    continue;
                end
                args = [args, {name, eqCfg.(name)}]; %#ok<AGROW>
            end
        end
    end
end
