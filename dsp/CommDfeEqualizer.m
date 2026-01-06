classdef CommDfeEqualizer < matlab.System
    %COMMD F EQUALIZER DFE wrapper around comm.DecisionFeedbackEqualizer.

    properties
        Algorithm = 'LMS';
        NumFeedforwardTaps = 20;
        NumFeedbackTaps = 20;
        StepSize = 5e-2;
        Lambda = 0.99;
        Delta = 1e-1;
        SampsPerSymb = 4;
        SampOffset = 1;
        ModOrder = 8;
    end

    properties (Access = private)
        DfeObj
        Constellation
        FfMidTapIndex
    end

    methods
        function obj = CommDfeEqualizer(varargin)
            if nargin > 0
                setProperties(obj, nargin, varargin{:});
            end
        end
    end

    methods (Access = protected)
        function setupImpl(obj, ~, ~, ~)
            obj.Constellation = pskmod(0:obj.ModOrder - 1, obj.ModOrder, pi / obj.ModOrder);
            obj.FfMidTapIndex = ceil(obj.NumFeedforwardTaps / 2);

            obj.DfeObj = comm.DecisionFeedbackEqualizer( ...
                'Algorithm', obj.Algorithm, ...
                'NumForwardTaps', obj.NumFeedforwardTaps, ...
                'NumFeedbackTaps', obj.NumFeedbackTaps, ...
                'ReferenceTap', obj.FfMidTapIndex, ...
                'Constellation', obj.Constellation, ...
                'TrainingFlagInputPort', true, ...
                'TrainingInputPort', true);

            if obj.isRls()
                obj.DfeObj.ForgettingFactor = obj.Lambda;
                obj.DfeObj.InitialInverseCovariance = (1 / obj.Delta) * eye(obj.NumFeedforwardTaps + obj.NumFeedbackTaps);
            else
                obj.DfeObj.StepSize = obj.StepSize;
            end

            reset(obj.DfeObj);
        end

        function resetImpl(obj)
            if ~isempty(obj.DfeObj)
                reset(obj.DfeObj);
            end
        end

        function [eqSymbs, errHist] = stepImpl(obj, matched, pilotSymbs, frameSymbType)
            totalSymbs = numel(frameSymbType);
            totalSamps = totalSymbs * obj.SampsPerSymb;

            if obj.SampOffset < 1 || obj.SampOffset > obj.SampsPerSymb
                error('CommDfeEqualizer:SampOffset', 'SampOffset must be between 1 and SampsPerSymb.');
            end

            symSamps = matched(obj.SampOffset:obj.SampsPerSymb:end);
            symSamps = symSamps(1:totalSymbs);

            trngFlag = false(totalSymbs, 1);
            trngSymbs = complex(zeros(totalSymbs, 1));
            pilotIdx = 0;

            for symIdx = 1:totalSymbs
                switch frameSymbType(symIdx)
                    case 1
                        pilotIdx = pilotIdx + 1;
                        trngFlag(symIdx) = true;
                        trngSymbs(symIdx) = pilotSymbs(pilotIdx);
                    otherwise
                        trngFlag(symIdx) = false;
                end
            end

            [eqSymbs, errVec] = obj.DfeObj(symSamps, trngSymbs, trngFlag);

            errHist = complex(zeros(totalSamps, 1));
            errHist(obj.SampOffset:obj.SampsPerSymb:end) = errVec;
        end

        function num = getNumInputsImpl(~)
            num = 3;
        end

        function num = getNumOutputsImpl(~)
            num = 2;
        end
    end

    methods (Access = private)
        function tf = isRls(obj)
            tf = strcmpi(obj.Algorithm, 'RLS');
        end
    end
end
