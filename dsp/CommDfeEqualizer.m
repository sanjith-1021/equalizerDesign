classdef CommDfeEqualizer < matlab.System
    %COMMD F EQUALIZER DFE wrapper around comm.DecisionFeedbackEqualizer.

    properties
        Algorithm = 'LMS';
        NumFeedforwardTaps = 20;
        NumFeedbackTaps = 20;
        StepSize = 5e-2;
        Lambda = 0.99;
        Delta = 1e-1;
        SamplesPerSymbol = 4;
        SampleOffset = 1;
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

        function [eqSymbols, errHist] = stepImpl(obj, matched, pilotSymbols, frameSymbType)
            totalSymbols = numel(frameSymbType);
            totalSamples = totalSymbols * obj.SamplesPerSymbol;

            if obj.SampleOffset < 1 || obj.SampleOffset > obj.SamplesPerSymbol
                error('CommDfeEqualizer:SampleOffset', 'SampleOffset must be between 1 and SamplesPerSymbol.');
            end

            symSamples = matched(obj.SampleOffset:obj.SamplesPerSymbol:end);
            symSamples = symSamples(1:totalSymbols);

            trainingFlag = false(totalSymbols, 1);
            trainingSymbols = complex(zeros(totalSymbols, 1));
            pilotIdx = 0;

            for symIdx = 1:totalSymbols
                switch frameSymbType(symIdx)
                    case 1
                        pilotIdx = pilotIdx + 1;
                        trainingFlag(symIdx) = true;
                        trainingSymbols(symIdx) = pilotSymbols(pilotIdx);
                    otherwise
                        trainingFlag(symIdx) = false;
                end
            end

            [eqSymbols, errVec] = obj.DfeObj(symSamples, trainingSymbols, trainingFlag);

            errHist = complex(zeros(totalSamples, 1));
            errHist(obj.SampleOffset:obj.SamplesPerSymbol:end) = errVec;
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
