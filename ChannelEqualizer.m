classdef ChannelEqualizer < matlab.System

    properties
        Algorithm = 'LMS';
        NumTaps = 20;
        StepSize = 5e-2;
        Lambda = 0.99;
        Delta = 1e-1;
        SamplesPerSymbol = 4;
        ModOrder = 8;
    end

    properties (Access = private)
        Weights
        DelayLine
        P
        MidTapIndex
    end

    methods
        function obj = ChannelEqualizer(varargin)
            if nargin > 0
                setProperties(obj, nargin, varargin{:});
            end
        end
    end

    methods (Access = protected)
        function setupImpl(obj, ~, ~, ~)
            reset(obj);
        end

        function resetImpl(obj)
            obj.MidTapIndex = ceil(obj.NumTaps / 2);
            obj.Weights = zeros(obj.NumTaps, 1);
            obj.Weights(obj.MidTapIndex) = 1;
            obj.DelayLine = complex(zeros(obj.NumTaps, 1));

            if obj.isRls()
                obj.P = (1 / obj.Delta) * eye(obj.NumTaps);
            else
                obj.P = [];
            end
        end

        function [eqSymbols, errHist] = stepImpl(obj, matched, pilotSymbols, frameSymbType)
            totalSymbols = numel(frameSymbType);
            totalSamples = totalSymbols * obj.SamplesPerSymbol;

            eqOut = complex(zeros(totalSamples, 1));
            errHist = complex(zeros(totalSamples, 1));

            pilotCount = 0;
            

            for dIdx = 1:totalSamples
                obj.DelayLine = [matched(dIdx); obj.DelayLine(1:end-1)];
                yi = obj.Weights' * obj.DelayLine;

                yIdx = dIdx - obj.MidTapIndex + 1;

                if yIdx > 0
                    eqOut(yIdx) = yi;

                    if mod(yIdx - 1, obj.SamplesPerSymbol) == 0
                        symbIdx = floor((yIdx - 1) / obj.SamplesPerSymbol) + 1;
                        symType = frameSymbType(symbIdx);

                        switch symType
                            case 1
                                pilotCount = pilotCount + 1;
                                hdSymb = pilotSymbols(pilotCount);
                            case 2
                                dInt = pskdemod(yi, obj.ModOrder, pi / obj.ModOrder);
                                hdSymb = pskmod(dInt, obj.ModOrder, pi / obj.ModOrder);
                            otherwise
                                hdSymb = 0;
                        end

                        if symType ~= 0
                            symErr = hdSymb - yi;
                            obj.updateWeights(symErr);
                            errHist(yIdx) = symErr;
                        end
                    end
                end
            end
            eqSymbols = eqOut(1:obj.SamplesPerSymbol:end);
        end

        function num = getNumInputsImpl(~)
            num = 3;
        end

        function num = getNumOutputsImpl(~)
            num = 2;
        end
    end

    methods (Access = private)
        function updateWeights(obj, symErr)
            if obj.isRls()
                u = obj.DelayLine;
                Pu = obj.P * u;
                denom = obj.Lambda + (u' * Pu);
                K = Pu / real(denom);

                obj.Weights = obj.Weights + K * conj(symErr);
                obj.P = (obj.P - K * (u' * obj.P)) / obj.Lambda;
            else
                obj.Weights = obj.Weights + obj.StepSize * conj(symErr) * obj.DelayLine;
            end
        end

        function tf = isRls(obj)
            tf = strcmpi(obj.Algorithm, 'RLS');
        end
    end
end
