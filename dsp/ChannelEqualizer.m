classdef ChannelEqualizer < matlab.System
    properties
        algorithm = 'LMS';
        nFeedforwardTaps = 20;
        nFeedbackTaps = 10;
        stepSize = 5e-2;
        forgetFactor = .99
        initInvCorr = 1e-2
        nSampsPerSymb = 4;
        modOrder = 8;
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
            obj.midTapIndex = obj.nFeedforwardTaps;                     % NOTE: One is suppose to be only on the last tap
            obj.ffWeights = zeros(obj.nFeedforwardTaps,1);
            obj.ffWeights(obj.midTapIndex) = 1;
            obj.fbWeights = zeros(obj.nFeedbackTaps,1);
            obj.ffDelayLine = complex(zeros(obj.nFeedforwardTaps,1));
            obj.fbDelayLine = complex(zeros(obj.nFeedbackTaps,1));

            if obj.isRls()
                totalTaps = obj.nFeedbackTaps + obj.nFeedforwardTaps;
                obj.P = (1/obj.initInvCorr)*eye(totalTaps);
            else
                obj.P = [];
            end
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
                        symbIdx = floor((yIdx - 1)/obj.nSampsPerSymb) + 1;
                        symbType = frameSymbType(symbIdx);

                        switch symbType
                            case 1 
                                trngSymbCnt = trngSymbCnt + 1;
                                hdSymb = trngSymbs(trngSymbCnt);
                            case 2
                                dInt = pskdemod(yi, obj.modOrder, pi/obj.modOrder);
                                hdSymb = pskmod(dInt, obj.modOrder, pi/obj.modOrder);
                            otherwise
                                hdSymb = 0;
                        end

                        if symbType ~= 0 
                            symbErr = hdSymb - yi;
                            obj.updateWeights(symbErr);
                            errHist(symbIdx)=symbErr;
                            obj.fbDelayLine = [hdSymb;obj.fbDelayLine(1: end-1)];
                        end
                    end
                end
            end
            eqSymbs = eqSamps(1:obj.nSampsPerSymb:end);
            % chanTaps = [obj.fbWeights;obj.ffWeights;]         % NOTE: Consider passing out 
        end

        function num = getNumInputsImpl(~)
            num = 3;
        end

        function num = getNumOutputsImpl(~)
            num = 2;
        end 
    end

    methods (Access = private)
        function updateWeights(obj, symbErr)
            if obj.isRls()
                u = [obj.ffDelayLine; -obj.fbDelayLine];
                Pu = obj.P * u;
                denom = obj.forgetFactor + (u'*Pu);
                K = Pu/real(denom);

                w = [obj.ffWeights; obj.fbWeights] + K*conj(symbErr);
                obj.ffWeights = w(1:obj.nFeedforwardTaps);
                obj.fbWeights = w(obj.nFeedforwardTaps + 1: end);

                obj.P = (obj.P - K * (u'*obj.P))/obj.forgetFactor;
            else
                obj.ffWeights = obj.ffWeights + obj.stepSize * conj(symbErr) * obj.ffDelayLine;
                obj.fbWeights = obj.fbWeights - obj.stepSize * conj(symbErr) * obj.fbDelayLine;
            end
        end

        function tf = isRls(obj)
            tf = strcmpi(obj.algorithm, 'RLS'); 
        end
    end
end