classdef ChannelEqualizer < matlab.system
    properties
        nFeedforwardTaps = 20;
        nFeedbackTaps = 8;
        nSampsPerSymb = 4;
        modOrder = 8;
        trainEq = struct(...
            'initInvCorr', 1e-2, ...
            'processNoise', 5e-3,...
            'measurementNoise', 1e-5);
        dataEq = struct(...
            'processNoise', 1e-3, ...
            'measurementNoise', 1e-1);

        pBits = 20;
        measBits = 50;
        kShift = 20;
        signalScale = 12;
    end

    properties (Access = private)
        ffWeightsRe
        ffWeightsIm
        fbWeightsRe
        fbWeightsIm
        ffDelayRe
        ffDelayIm
        fbDelayRe
        fbDelayIm
        P
        midTapIndex
        trainEqFp
        dataEqFp
        constRe
        constIm
    end

    methods
        function obj = ChannelEqualizer(varargin)
            if nargin > 0, setProperties(obj, nargin, varargin{:}); end
        end
    end


    methods (Access = protected)
        function setupImpl(obj, ~, ~, ~)
            reset(obj);
        end

        function resetImpl(obj)
            obj.midTapIndex = obj.nFeedforwardTaps;
            obj.trainEqFp = obj.toEqFixed(obj.trainEq);
            obj.dataEqFp = obj.toEqFixed(obj.dataEq);

            obj.ffWeightsRe = zeros(obj.nFeedforwardTaps, 1, 'int16');
            obj.ffWeightsIm = zeros(obj.nFeedforwardTaps, 1, 'int16');
            obj.ffWeightsRe(obj.midTapIndex) = obj.signalScale;
            obj.fbWeightsRe = zeros(obj.nFeedbackTaps, 1, 'int16');
            obj.fbWeightsIm = zeros(obj.nFeedforwardTaps, 1, 'int16');

            obj.ffDelayRe = zeros(obj.nFeedforwardTaps, 1, 'int16');
            obj.ffDelayIm = zeros(obj.nFeedforwardTaps, 1, 'int16');
            obj.fbDelayRe = zeros(obj.nFeedbackTaps, 1, 'int16');
            obj.fbDelayIm = zeros(obj.nFeedbackTpas, 1, 'int16');

            totalTaps = obj.nFeedbackTaps + obj.nFeedforwardTaps;
            pInit = round((1/obj.trainEq.intiInvCorr) * 2^obj.pBits);
            obj.P = int32(pInit) * eye(totalTaps, 'int32');

            [obj.constRe, obj.constIm] = obj.buildConstellation();
        end

        function [eqSymbs, errHist] = stepImpl(obj, rxSymbs, trngSymbs, frameSymbType)
            totalSymbs = numel(frameSymbType);
            totalSamps = totalSymbs*obj.nSampsPerSymb;

            rxRe = real(rxSymbs); rxIm = imag(rxSymbs);
            trngRe = real(trngSymbs), trngIm = imag(trngSymbs);

            eqSampsRe = zeros(totalSamps, 1, 'int16');
            eqSampsIm = zeros(totalSamps, 1, 'int16');
            errRe = zeros(totalSymbs, 1, 'int16'); 
            errIm = zeros(totalSymbs, 1, 'int16');

            trngSymbCnt = 0;

            for dIdx = 1:totalSamps
                obj.ffDelayRe = [rxRe(dIdx); obj.ffDelayRe(1:end-1)];
                obj.ffDelayIm = [rxIm(dIdx); obj.ffDelayIm(1:end-1)];
                [yiRe, yiIm] = obj.fitlerOutput();
                yIdx = dIdx - obj.midTapIndex + 1;

                if yIdx > 0
                    eqSampRe(yIdx) = yiRe; eqSampIm(yIdx) = yiIm;

                    if mod(yIdx-1, obj.nSampsPerSymb) == 0
                        symbIdx = floor((yIdx - 1)/obj.nSampsPerSymb) + 1;
                        symbType = frameSymbType(symbIdx);

                        switch symbType
                            case {1, 2}
                                trngSymbCnt = trngSymbCnt + 1;
                                hdRe = trngRe(trngSymbCnt); hdIm = trngIm(trngsymbCnt);
                                eqParams = obj.trainEqFp;
                            case 3
                                [hdRe, hdIm] = obj.hardDecision(yiRe, yiIm);
                                eqParams = obj.dataEqFp;
                            otherwise
                                hdRe = int16(0); hdIm = int16(0); eqParams = obj.trainEqFp;
                        end

                        if symbType ~= 0 
                            symbErrRe = int32(hdRe) - int32(yiRe);
                            symbErrIm = int32(hdIm) - int32(yiIm);
                            errRe(symbIdx) = satInt16(symbErrRe); errIm(symbIdx) = satInt16(symbErrIm);
                            obj.updateWeightsFixed(symbErrRe, symbErrIm, eqParams);
                            obj.fbDelayRe = [hdRe; obj.fbDelayRe(1:end-1)];
                            obj.fbDelayIm = [hdIm; obj.fbDelayIm(1:end-1)];
                        end
                    end
                end

            end
            
            eqSymbsRe = eqSampRe(1:obj.nSampsPerSymb:end);
            eqSymbsIm = eqSampIm(1:obj.nSampsPerSymb:end);
            eqSymbs = complex(eqSymbsRe, eqSymbsIm);
            errHist = complex(errRe, errIm);

        end

        function num = getNumInputsImpl(~)
            num = 3;
        end

        function num = getNumOutputsImpl(~)
            num = 2;
        end 
    end


    methods (Access = private)

        function updateWeightsFixed(obj, errRe, errIm, eqParams)
            nFf = obj.nFeedforwardTaps;
            nFb = obj.nFeedbackTaps;
            totalTaps = nFf + nFb;

            uRe = [obj.ffDelayRe; -obj.fbDelayRe];
            uIm = [obj.ffDelayIm; -obj.fbDelayIm];

            predP = obj.P;
            for idx = 1:totalTaps
                predP(idx, idx) = predP(idx, idx) + eqParams.processNoise;
            end

            vRe = zeros(totalTaps, 1, 'int64');
            vIm = zeros(totalTaps, 1, 'int64');
            for row = 1:totalTaps
                accRe = int64(0);
                accIm = int64(0);
                for col = 1:totalTaps
                    pij = predP(row, col);
                    accRe = accRe + pij * int64(uRe(col));
                    accIm = accIm + pij * int64(uIm(col));
                end
                vRe(row) = accRe;
                vIm(row) = accIm;
            end

            denom = int64(0);
            for idx = 1:totalTaps
                denom = denom + int64(uRe(idx)) * vRe(idx) + int64(uIm(idx)) * vIm(idx);
            end
            denom = denom + eqParams.measurementNoise;
            if denom == 0
                denom = int64(1);
            end

            kFrac = obj.kShift - obj.signalShift;
            KRe = idivide(bitshift(vRe, obj.kShift), denom, 'fix');
            KIm = idivide(bitshift(vIm, obj.kShift), denom, 'fix');

            errRe64 = int64(errRe);
            errIm64 = int64(errIm);
            deltaRe = idivide(KRe .* errRe64 + KIm .* errIm64, int64(2^kFrac), 'fix');
            deltaIm = idivide(KRe .* errIm64 - KIm .* errRe64, int64(2^kFrac), 'fix');

            obj.ffWeightsRe = obj.satWeights(obj.ffWeightsRe + int32(deltaRe(1:nFf)));
            obj.ffWeightsIm = obj.satWeights(obj.ffWeightsIm + int32(deltaIm(1:nFf)));
            obj.fbWeightsRe = obj.satWeights(obj.fbWeightsRe + int32(deltaRe(nFf+1:end)));
            obj.fbWeightsIm = obj.satWeights(obj.fbWeightsIm + int32(deltaIm(nFf+1:end)));

            rowRe = zeros(1, totalTaps, 'int64');
            rowIm = zeros(1, totalTaps, 'int64');
            for col = 1:totalTaps
                accRe = int64(0);
                accIm = int64(0);
                for row = 1:totalTaps
                    pij = predP(row, col);
                    accRe = accRe + int64(uRe(row)) * pij;
                    accIm = accIm + int64(uIm(row)) * pij;
                end
                rowRe(col) = accRe;
                rowIm(col) = -accIm;
            end

            for row = 1:totalTaps
                kRe = KRe(row);
                kIm = KIm(row);
                for col = 1:totalTaps
                    realPart = kRe * rowRe(col) - kIm * rowIm(col);
                    update = idivide(realPart, int64(2^obj.pBits), 'fix');
                    predP(row, col) = predP(row, col) - update;
                end
            end
            obj.P = predP;
        end


        function [yiRe, yiIm] = filterOutput(obj)
            accRe = int32(0); accIm = int32(0);
            for idx = 1:obj.nFeedforwardTaps
                wr = int64(obj.ffWeightsRe(idx));
                wi = int64(obj.ffWeightsIm(idx));
                xr = int64(obj.ffDelayRe(idx));
                xi = int64(obj.ffDelayIm(idx));
                accRe = accRe + wr*xr - wi*xi;
                accIm = accIm + wr*xi + wi*xr;
            end
            for idx = 1:obj.nFeedbackTaps
                wr = int64(obj.fbWeightsRe(idx));
                wi = int64(obj.ffWeightsIm(idx));
                xr = int64(obj.ffDelayRe(idx));
                xi = int64(obj.ffDelayIm(idx));
                accRe = accRe + wr*xr - wi*xi;
                accIm = accIm + wr*xi - wi*xr;
            end
            yiRe = obj.satInt16(idivide(accRe, int64(2^obj.signalShift), 'fix'));
            yiIm = obj.satInt16(idivide(accIm, int64(2^obj.signalShift), 'fix'));
        end

        function y = satInt16(~, x)
            lo = int32(-32768);
            hi = int32(32767);
            x = int32(x);
            x = min(max(x, lo), hi);
            y = int16(x);
        end

        function [re, im] = buildConstellation(obj)
            k = (0:obj.modOrder-1).';
            symbs = pskmod(k, obj.modOrder, 0, 'gray');
            [re, im] = quantizeComplex(symbs);
        end

        function eqFp = toEqFixed(obj, eq)
            eqFp = struct();
            eqFp.processNoise = int32(round(eq.processNoise * 2^obj.pBits));
            eqFp.measurementNoise = int32(round(eq.measurementNoise*2^obj.measBits));
        end
    end
end