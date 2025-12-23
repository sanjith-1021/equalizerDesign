classdef ChannelEqualizer < matlab.System
    properties
        algorithm = 'LMS';
        nFeedforwardTaps = 20;
        nFeedbackTaps = 10;
        stepSize = 5e-2;
        forgetFactor = 0.99;
        initInvCorr = 1e-2;
        nSampsPerSymb = 4;
        modOrder = 8;

        % Kalman filter parameters (tap-state model)
        %
        % For the scalar measurement model used here:
        %   d[n] = u[n]^H w[n] + v[n],   v ~ CN(0, R)
        % and random-walk tap evolution:
        %   w[n] = w[n-1] + q[n-1],      q ~ CN(0, Q)
        %
        % measurementNoise / measureNoise: observation noise variance R (scalar)
        % processNoise: process noise Q (scalar, vector (diag), or matrix)
        measurementNoise = 1e-3;
        measureNoise = [];            % optional alias; if set, overrides measurementNoise
        processNoise = 1e-5;
        initialCovariance = [];       % optional override for P0 scale (scalar)
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
            obj.midTapIndex = obj.nFeedforwardTaps;
            obj.ffWeights = zeros(obj.nFeedforwardTaps,1);
            obj.ffWeights(obj.midTapIndex) = 1;
            obj.fbWeights = zeros(obj.nFeedbackTaps,1);
            obj.ffDelayLine = complex(zeros(obj.nFeedforwardTaps,1));
            obj.fbDelayLine = complex(zeros(obj.nFeedbackTaps,1));

            if obj.isRls() || obj.isKalman()
                totalTaps = obj.nFeedbackTaps + obj.nFeedforwardTaps;
                p0Scale = 1/obj.initInvCorr;
                if ~isempty(obj.initialCovariance)
                    p0Scale = obj.initialCovariance;
                end
                obj.P = p0Scale * eye(totalTaps);
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

            elseif obj.isKalman()
                u = [obj.ffDelayLine; -obj.fbDelayLine];
                totalTaps = numel(u);

                % Build Q (process noise covariance)
                qSpec = obj.processNoise;
                if isempty(qSpec)
                    qMat = zeros(totalTaps);
                elseif isscalar(qSpec)
                    qMat = real(qSpec) * eye(totalTaps);
                elseif isvector(qSpec) && numel(qSpec) == totalTaps
                    qMat = diag(real(qSpec(:)));
                elseif isequal(size(qSpec), [totalTaps totalTaps])
                    qMat = real(qSpec);
                else
                    qMat = real(qSpec(1)) * eye(totalTaps);
                end

                % Measurement noise variance R (allow alias measureNoise)
                rVar = obj.measurementNoise;
                if ~isempty(obj.measureNoise)
                    rVar = obj.measureNoise;
                end
                rVar = real(rVar);
                if rVar <= 0
                    rVar = eps;
                end

                % Predict (random-walk tap model: F = I)
                pPred = obj.P + qMat;

                % Update (scalar measurement)
                s = real((u' * pPred * u) + rVar);
                if s <= 0
                    s = eps;
                end
                K = (pPred * u) / s;

                w = [obj.ffWeights; obj.fbWeights] + K * conj(symbErr);
                obj.ffWeights = w(1:obj.nFeedforwardTaps);
                obj.fbWeights = w(obj.nFeedforwardTaps + 1:end);

                % Joseph stabilized covariance update (keeps P PSD)
                iMat = eye(totalTaps);
                pNew = (iMat - K * u') * pPred * (iMat - K * u')' + (K * rVar * K');
                obj.P = (pNew + pNew')/2;
            else
                obj.ffWeights = obj.ffWeights + obj.stepSize * conj(symbErr) * obj.ffDelayLine;
                obj.fbWeights = obj.fbWeights - obj.stepSize * conj(symbErr) * obj.fbDelayLine;
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