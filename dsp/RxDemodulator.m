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
