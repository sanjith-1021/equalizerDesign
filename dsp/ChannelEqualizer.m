classdef ChannelEqualizer < matlab.System
    properties
        algorithm = 'LMS';
        nFeedforwardTaps = 20;
        nFeedbackTaps = 8;
        nSampsPerSymb = 4;
        modOrder = 8;
        trainEq = struct();
        dataEq = struct();
    end

    properties (Access = private)
        ffWeights 
        fbWeights
        ffDelayLine
        fbDelayLine
        P
        Q
        R
        midTapIndex

        measuredNoise = []
        initialCovariance
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
            obj.midTapIndex = obj.nFeedforwardTaps;
            obj.ffWeights = zeros(obj.nFeedforwardTaps,1);
            obj.ffWeights(obj.midTapIndex) = 1;
            obj.fbWeights = zeros(obj.nFeedbackTaps,1);
            obj.ffDelayLine = complex(zeros(obj.nFeedforwardTaps,1));
            obj.fbDelayLine = complex(zeros(obj.nFeedbackTaps,1));

            totalTaps = obj.nFeedbackTaps + obj.nFeedforwardTaps;
            if obj.isRls()
                obj.P = (1/obj.trainEq.initInvCorr)*eye(totalTaps);
            elseif obj.isKalman()
                obj.P = (1/obj.trainEq.initInvCorr)*eye(totalTaps);
                obj.Q = obj.trainEq.processNoise * eye(totalTaps);
                obj.R = obj.trainEq.measurementNoise;
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
                            case {1, 2} 
                                trngSymbCnt = trngSymbCnt + 1;
                                hdSymb = trngSymbs(trngSymbCnt);
                                eqParams = obj.trainEq;
                            case 3
                                dInt = pskdemod(yi, obj.modOrder, 0, 'gray');
                                hdSymb = pskmod(dInt, obj.modOrder, 0, 'gray');
                                eqParams = obj.dataEq;
                            otherwise
                                hdSymb = 0;
                                eqParams = obj.trainEq;
                        end

                        if symbType ~= 0 
                            symbErr = hdSymb - yi;
                            errHist(symbIdx)=symbErr;
                            if symbType~=3, obj.updateWeights(symbErr, eqParams); end
                            obj.fbDelayLine = [hdSymb;obj.fbDelayLine(1: end-1)];
                        end
                    end
                end
            end
            eqSymbs = eqSamps(1:obj.nSampsPerSymb:end);
            % chanTaps = [obj.fbWeights;obj.ffWeights;] ;        
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
            if obj.isRls()
                u = [obj.ffDelayLine; -obj.fbDelayLine];
                Pu = obj.P * u;
                denom = eqParams.forgetFactor + (u'*Pu);
                K = Pu/real(denom);

                w = [obj.ffWeights; obj.fbWeights] + K*conj(symbErr);
                obj.ffWeights = w(1:obj.nFeedforwardTaps);
                obj.fbWeights = w(obj.nFeedforwardTaps + 1: end);

                obj.P = (obj.P - K * (u'*obj.P))/eqParams.forgetFactor;

            elseif obj.isKalman()
                u = [obj.ffDelayLine; -obj.fbDelayLine];
                totalTaps = obj.nFeedforwardTaps + obj.nFeedbackTaps;
                predP = obj.P + eqParams.processNoise * eye(totalTaps);
                denom = real(u' * predP * u) + eqParams.measurementNoise;
                K = predP * u / denom;

                w = [obj.ffWeights; obj.fbWeights] + K*conj(symbErr);
                obj.ffWeights = w(1:obj.nFeedforwardTaps);
                obj.fbWeights = w(obj.nFeedforwardTaps + 1: end);

                obj.P = predP - K * (u' * predP);

                % Josepth Stabilized Form
                % iMat = eye(obj.nFeedforwardTaps + obj.nFeedbackTaps);
                % newP = (iMat - K * u') * predP * (iMat - K * u')' + K * obj.R * K';
                % obj.P = (newP + newP')/2; 
            else
                obj.ffWeights = obj.ffWeights + eqParams.stepSize * conj(symbErr) * obj.ffDelayLine;
                obj.fbWeights = obj.fbWeights - eqParams.stepSize * conj(symbErr) * obj.fbDelayLine;
            end
        end

        function tf = isRls(obj)
            tf = strcmpi(obj.algorithm, 'RLS'); 
        end

        function tf = isKalman(obj)
            tf = strcmpi(obj.algorithm, 'Kalman');
        end
    end
end
