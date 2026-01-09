classdef KalmanEqualizer < matlab.System
    properties
        nFeedforwardTaps = 32;
        nFeedbackTaps = 12;
        nSampsPerSymb = 4;
        modOrder = 8;
        trainEq = struct( ...
            'initInvCorr', 1e-2, ...
            'processNoise', 1e-2, ...
            'measurementNoise', 1e-4);
        dataEq = struct( ...
            'processNoise', 1e-4, ...
            'measurementNoise', 5e-1);
    end

    properties (Access = private)
        ffWeights
        fbWeights
        ffDelayLine
        fbDelayLine
        P
        midTapIndex
    end

    methods
        function obj = KalmanEqualizer(varargin)
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
            obj.midTapIndex = obj.nFeedforwardTaps;
            obj.ffWeights = zeros(obj.nFeedforwardTaps, 1);
            obj.ffWeights(obj.midTapIndex) = 1;
            obj.fbWeights = zeros(obj.nFeedbackTaps, 1);
            obj.ffDelayLine = complex(zeros(obj.nFeedforwardTaps, 1));
            obj.fbDelayLine = complex(zeros(obj.nFeedbackTaps, 1));

            totalTaps = obj.nFeedbackTaps + obj.nFeedforwardTaps;
            obj.P = (1 / obj.trainEq.initInvCorr) * eye(totalTaps);
        end

        function [eqSymbs, errHist] = stepImpl(obj, rxSymbs, trngSymbs, frameSymbType)
            totalSymbs = numel(frameSymbType);
            totalSamps = totalSymbs * obj.nSampsPerSymb;

            eqSamps = complex(zeros(totalSamps, 1));
            errHist = complex(zeros(totalSymbs, 1));

            trngSymbCnt = 0;

            for dIdx = 1: totalSamps
                obj.ffDelayLine = [rxSymbs(dIdx); obj.ffDelayLine(1:end-1)];
                yi = obj.ffWeights' * obj.ffDelayLine - obj.fbWeights' * obj.fbDelayLine;
                yIdx = dIdx - obj.midTapIndex + 1;

                if yIdx > 0
                    eqSamps(yIdx) = yi;

                    if mod(yIdx - 1, obj.nSampsPerSymb) == 0
                        symbIdx = floor((yIdx - 1) / obj.nSampsPerSymb) + 1;
                        symbType = frameSymbType(symbIdx);

                        switch symbType
                            case {1, 2}
                                trngSymbCnt = trngSymbCnt + 1;
                                hdSymb = trngSymbs(trngSymbCnt);
                                eqParams = obj.trainEq;
                            case 3
                                dInt = pskdemod(yi, obj.modOrder, pi / obj.modOrder);
                                hdSymb = pskmod(dInt, obj.modOrder, pi / obj.modOrder);
                                eqParams = obj.dataEq;
                            otherwise
                                hdSymb = 0;
                                eqParams = obj.trainEq;
                        end

                        if symbType ~= 0
                            symbErr = hdSymb - yi;
                            errHist(symbIdx) = symbErr;
                            obj.updateWeights(symbErr, eqParams);
                            obj.fbDelayLine = [hdSymb; obj.fbDelayLine(1:end-1)];
                        end
                    end
                end
            end
            eqSymbs = eqSamps(1:obj.nSampsPerSymb:end);
        end

        function num = getNumInputsImpl(~)
            num = 3;
        end

        function num = getNumOutputsImpl(~)
            num = 2;
        end
    end

    methods (Access = private)
        function updateWeights(obj, symbErr, eqParams)
            u = [obj.ffDelayLine; -obj.fbDelayLine];
            totalTaps = obj.nFeedforwardTaps + obj.nFeedbackTaps;

            predP = obj.P + eqParams.processNoise * eye(totalTaps);
            denom = real(u' * predP * u) + eqParams.measurementNoise;
            K = predP * u / denom;

            w = [obj.ffWeights; obj.fbWeights] + K * conj(symbErr);
            obj.ffWeights = w(1:obj.nFeedforwardTaps);
            obj.fbWeights = w(obj.nFeedforwardTaps + 1: end);

            obj.P = predP - K * (u' * predP);
        end
    end
end
