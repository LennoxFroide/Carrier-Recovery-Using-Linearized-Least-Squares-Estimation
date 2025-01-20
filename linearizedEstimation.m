%% CARRIER RECOVERY IN A COMMUNICATION SYSTEM IS A KEY PROBLEM IN COHERENT DEMODULATION
%% IN THIS ALGORITHM WE RECOVER A SIGNAL WITH EMBEDDED GAUSSIAN WHITE NOISE PROCESS
%% USING A LINEAR LEAST SQUARES ESTIMATION (LSE).
% theta = [0.25, pi/10];
% Generating the variances for the white noise gaussian random processes
function noiseVar = genNoiseVar(snwValue)
    % Since snw = signal var/noise var and signal has no source
    % of randomness, we can generate the noise variance using inverses
    noiseVar = (snwValue).^-1;
end

% Generating noise-free carries 
function noiseFreeSig = genNoiseFreeCarrier(N, theta)
    nValues = -N:1:N;
    noiseSigVector = NaN(1,size(nValues,2));
    for noiseVectorIdx = 1:1:size(nValues,2)
        % Undelying fomular: S(n) = cos(2pif(o)N + phi) n = -N,...,0,...N
        noiseFreeSig = cos(2*pi*theta(1)*nValues(noiseVectorIdx) + theta(2));
        noiseSigVector(1, noiseVectorIdx) = noiseFreeSig;
    end
    noiseFreeSig = noiseSigVector;
end

% Generating noisy process and noisy signal
function noisySig = genNoisyCarrier(noiseVar, mean, noiseFreeSig, totalTime)
    % Underlying fomular: x(n) = s(n) + x'(n)
    % where x'(n) is a white gaussian random process
    % Array to hold the sequence of iid gaussian random variables
    sequence = NaN(1,1);
    % Using 1 millisecond sampling observations
    times = 0.001:0.001:totalTime;
    % Getting the standard deviation for our random process
    standDev = sqrt(noiseVar);
    % Pointer to add generated rv to the array
    index = 1;
    % Generating the actual sequence
    for time = [times(1), diff(times)]
        pause(time);
        % Gerating 1 RV at eacha millisecond
        randomVar = normrnd(mean, standDev);
        % Adding this rv to the array
        sequence(1, index) = randomVar;
        %  incrementing our pointer
        index = index + 1;
    end
    % Generating the random process using cumsum() function and getting 
    % a realization
    % noisySig = noiseFreeSig + cumsum(sequence);
    noisySig = noiseFreeSig + sequence(1,500);
end

% The test estimator and a function to update it between runs
testN = 1;
testEstimator =[(8*(pi.^2)* testN.^3)/3, 0; 0, testN];
function updTestEst = updateTestEstimator(currTestEst, newN)
    % We only need to update the value at indices:
    % testEst(1,1) and testEst(2,2)
    currTestEst(1,1) = (8*(pi.^2)* newN.^3)/3;
    currTestEst(2,2) = newN;
    % This is already in the format H(theta)^T(H(theta)
    % This will now be the new test Estimator
    updTestEst = currTestEst;
end

% The ideal estimator and a function to update it between runs
idealN = 1;
% idealEstimator = Nan(2,2); % A 2*2 matrix
% theta = [f(o), phi]
theta = [0.25,pi/10];
% Elements of the ideal estimator matrix
idealR11 = 4*pi.^2*(sin(2*pi*theta(1) + theta(2))).^2 + ...
    16*pi.^2*(sin(4*pi*theta(1) + theta(2))).^2;

idealR12 = 2*pi*(sin(2*pi*theta(1) + theta(2))).^2 + ...
    4*pi*(sin(4*pi*theta(1) + theta(2))).^2;

idealR21 = 2*pi*(sin(2*pi*theta(1) + theta(2))).^2 + ...
    4*pi*(sin(4*pi*theta(1) + theta(2))).^2;

idealR22 = sin(2*pi*theta(1) + theta(2)).^2 + ...
    sin(4*pi*theta(1) + theta(2)).^2;
% Building the actual estimator
idealEstimator = [idealR11, idealR12; idealR21, idealR22];

% This function will update the ideal estimator based on past realizations of frequency and phase
function updIdealEst = updateIdealEstimator(currIdealEst, newTheta)
    % We have to update each of the four elements in the ideal estimator
    % matrix
    % This is already in the format H(theta)^T(H(theta)
    currIdealEst(1,1) = 4*pi.^2*(sin(2*pi*newTheta(1) + newTheta(2))).^2 + ...
    16*pi.^2*(sin(4*pi*newTheta(1) + newTheta(2))).^2;

    currIdealEst(1,2) = 2*pi*(sin(2*pi*newTheta(1) + newTheta(2))).^2 + ...
    4*pi*(sin(4*pi*newTheta(1) + newTheta(2))).^2;

    currIdealEst(2,1) = 2*pi*(sin(2*pi*newTheta(1) + newTheta(2))).^2 + ...
    4*pi*(sin(4*pi*newTheta(1) + newTheta(2))).^2;

    currIdealEst(2,2) = sin(2*pi*newTheta(1) + newTheta(2)).^2 + ...
    sin(4*pi*newTheta(1) + newTheta(2)).^2;

    % Returning this updated estimator
    updIdealEst = currIdealEst;
end

hMatrix11 = -2*pi*(sin(2*pi*theta(1) + theta(2)));
hMatrix12 = -sin(2*pi*theta(1) + theta(2));
% hMatrix21 = -4*pi*(sin(4*pi*theta(1) + theta(2)));
% hMatrix22 = -sin(4*pi*(theta(1) + theta(2)));

% hMatrix = [hMatrix11, hMatrix12; hMatrix21, hMatrix22];
hMatrix = [hMatrix11, hMatrix12];

% This function updates the standard matrix based on past realizations of frequency and phase
function updHMatrix = updateHMatrix(currHMat, newTheta)
    % We have to update each of the four elements in the matrix
    currHMat(1,1) = -2*pi*(sin(2*pi*newTheta(1) + newTheta(2)));
    currHMat(1,2) = -sin(2*pi*newTheta(1) + newTheta(2));
    % currHMat(2,1) = -4*pi*(sin(4*pi*newTheta(1) + newTheta(2)));
    % currHMat(2,2) = -sin(4*pi*(newTheta(1) + newTheta(2)));
    
    updHMatrix = currHMat;
end

% This function generates the array of N values used for iterations
function nArray = genNArray(nValue)
    % n = -nValue,...,0,...,nvalue
    nGen = 1:nValue;
    nArray = nGen;
end

% This function will take in the estimators and return an estimate for
% theta
function thetaEstimate = estimateThetaValue(estimator, hMatrix, noisySig, noiseFreeSig, prevTheta)
    % The formular for iteration is given us:
    % theta(k+1) = theta(k) + (estimator.^-1) * hMatrix.' *...
    % (noisySig - noiseFreeSig)
    thetaEstimate = prevTheta.' + (estimator.^-1 * hMatrix.') .* median(noisySig - noiseFreeSig);
    thetaEstimate = thetaEstimate.';
end

%% COLLECTING UNIFORM SAMPLES FOR 1 SECOND WITH RANDOM SEQUENCE OF LENGTH N=50, AND SIGNAL TO NOISE RATIOS
%% DECAYING IN THE ORDER: 30dB,20dB,10dB to 0dB THE TRUE FREQUENCY f(0) = 0.25Hz AND phi = pi/10
% Initial frequency f(0)Hz and phi
%{
thetaQ1 = [0.25, pi/10];
idealThetaQ1 = thetaQ1;
testThetaQ1 = thetaQ1;
hMatrixThetaQ1 = hMatrix;
%}

meanQ1 = 0;
timeInterval = 1; % In seconds
% Estimation window
q1EstimationWindow = 50;
% Signal-to-noise ratios in decibels
signalToNoiseRatios = [30, 20, 10, 0];
% Our estimators
q1IdealEstimator = idealEstimator;
q1TestEstimator = testEstimator;
% HMatrices
q1IdealHMatrix = hMatrix;
q1TestHMatrix = hMatrix;
% Matrices to store the generated estimates
idealOutputsMatrixQ1 = cell(size(signalToNoiseRatios,2),q1EstimationWindow);
testOutputsMatrixQ1 = cell(size(signalToNoiseRatios,2),q1EstimationWindow);
% Generating windows
windowArray = genNArray(q1EstimationWindow);
% We will use both estimators and different snr to evaluate the impact of
% the variance of the white gaussian random process on the estimates
for snrIdnex = 1:1:size(signalToNoiseRatios, 2)
    % Getting a signal-to-noise ratio value
    currentSNR = signalToNoiseRatios(1, snrIdnex);
    % Getting the variance from this signal-t-noise ratio value
    currentVariance = genNoiseVar(currentSNR);
    % Index to add outputs from estimations to a matrix
    outputsIndex = 1;
    thetaQ1 = [0.25, pi/10];
    idealThetaQ1 = thetaQ1;
    testThetaQ1 = thetaQ1;
    hMatrixThetaQ1 = hMatrix;
    % Then we will use all of the windows to run the models
    for windowIndex = 1:1:size(windowArray,2)
        % Getting a window value
        currentWindowValue = windowArray(1, windowIndex);
        % Generating the noise free carrier for this windowValue
        noiseFreeCarrier = genNoiseFreeCarrier(currentWindowValue, thetaQ1);
        % Generating the noisy signal carrier for this window value
        noisySignalCarrier = genNoisyCarrier(currentVariance, meanQ1, noiseFreeCarrier, timeInterval);
        % Updating the two estimators
        q1IdealEstimator = updateIdealEstimator(q1IdealEstimator, idealThetaQ1);
        q1IdealHMatrix = updateHMatrix(q1IdealHMatrix, idealThetaQ1);
        q1TestEstimator = updateTestEstimator(q1TestEstimator, currentWindowValue);
        q1TestHMatrix = updateHMatrix(q1TestHMatrix, testThetaQ1);
        % Getting the estimated values of theta
        idealThetaQ1 = estimateThetaValue(q1IdealEstimator, q1IdealHMatrix, noisySignalCarrier, noiseFreeCarrier, idealThetaQ1);
        %q1IdealHMatrix = updateHMatrix(q1IdealHMatrix, hMatrixThetaQ2);
        testThetaQ1 = estimateThetaValue(q1TestEstimator, q1TestHMatrix, noisySignalCarrier, noiseFreeCarrier, testThetaQ1);
        % We need to store these values from one iteration to the next
        idealOutputsMatrixQ1{snrIdnex, outputsIndex} = idealThetaQ1;
        testOutputsMatrixQ1{snrIdnex, outputsIndex} = testThetaQ1;
        % Increasing the index to store the next value
        outputsIndex = outputsIndex + 1;
    end
end


%% COLLECTING UNIFORM SAMPLES FOR 1 SECOND WITH RANDOM SEQUENCE OF LENGTH N=100, AND SIGNAL TO NOISE RATIOS
%% DECAYING IN THE ORDER: 30dB,20dB,10dB to 0dB THE TRUE FREQUENCY f(0) = 0.25Hz AND phi = pi/10
% We will use the same theta and time interval but change the estimation
% window to 100
% Part (1.)
%{
thetaQ2A = [0.25, pi/10];
idealThetaQ2A = thetaQ2A;
testThetaQ2A = thetaQ2A;
hMatrixThetaQ2 = hMatrix;
%}

meanQ1 = 0;
timeInterval = 1; % In seconds
% Estimation window
q2EstimationWindow = 100;
% Signal-to-noise ratios in decibels
signalToNoiseRatios = [30, 20, 10, 0];
% Our estimators
q2IdealEstimator = idealEstimator;
q2TestEstimator = testEstimator;
% HMatrices
q2IdealHMatrix = hMatrix;
q2TestHMatrix = hMatrix;
% Matrices to store the generated estimates
idealOutputsMatrixQ2A = cell(size(signalToNoiseRatios,2),q2EstimationWindow);
testOutputsMatrixQ2A = cell(size(signalToNoiseRatios,2),q2EstimationWindow);
% Generating windows
windowArray = genNArray(q2EstimationWindow);
% We will use both estimators and different snr to evaluate the impact of
% the variance of the white gaussian random process on the estimates
for snrIdnex = 1:1:size(signalToNoiseRatios, 2)
    % Getting a signal-to-noise ratio value
    currentSNR = signalToNoiseRatios(1, snrIdnex);
    % Getting the variance from this signal-t-noise ratio value
    currentVariance = genNoiseVar(currentSNR);
    % Index to add outputs from estimations to a matrix
    outputsIndex = 1;
    % Matrices for estimation
    thetaQ2A = [0.25, pi/10];
    idealThetaQ2A = thetaQ2A;
    testThetaQ2A = thetaQ2A;
    hMatrixThetaQ2 = hMatrix;
    % Then we will use all of the windows to run the models
    for windowIndex = 1:1:size(windowArray,2)
        % Getting a window value
        currentWindowValue = windowArray(1, windowIndex);
        % Generating the noise free carrier for this windowValue
        noiseFreeCarrier = genNoiseFreeCarrier(currentWindowValue, thetaQ2A);
        % Generating the noisy signal carrier for this window value
        noisySignalCarrier = genNoisyCarrier(currentVariance, meanQ1, noiseFreeCarrier, timeInterval);
        % Updating the two estimators
        q2IdealEstimator = updateIdealEstimator(q2IdealEstimator, idealThetaQ2A);
        q2IdealHMatrix = updateHMatrix(q2IdealHMatrix, idealThetaQ2A);
        q2TestEstimator = updateTestEstimator(q2TestEstimator, currentWindowValue);
        q2TestHMatrix = updateHMatrix(q2TestHMatrix, testThetaQ2A);
        % Getting the estimated values of theta
        idealThetaQ2A = estimateThetaValue(q2IdealEstimator, q2IdealHMatrix, noisySignalCarrier, noiseFreeCarrier, idealThetaQ2A);
        testThetaQ2A = estimateThetaValue(q2TestEstimator, q2TestHMatrix, noisySignalCarrier, noiseFreeCarrier, testThetaQ2A);
        % We need to store these values from one iteration to the next
        idealOutputsMatrixQ2A{snrIdnex, outputsIndex} = idealThetaQ2A;
        testOutputsMatrixQ2A{snrIdnex, outputsIndex} = testThetaQ2A;
        % Increasing the index to store the next value
        outputsIndex = outputsIndex + 1;
    end
end


%% COLLECTING UNIFORM SAMPLES FOR 2 SECONDS WITH RANDOM SEQUENCE OF LENGTH N=100, AND SIGNAL TO NOISE RATIOS
%% DECAYING IN THE ORDER: 30dB,20dB,10dB to 0dB THE TRUE FREQUENCY f(0) = 0.25Hz AND phi = pi/10
% We will use the same theta and time interval but change the estimation
% window to 100 and change the time interval to 2 seconds as well
%{
thetaQ2B = [0.25, pi/10];
idealThetaQ2B = thetaQ2B;
testThetaQ2B = thetaQ2B;
hMatrixThetaQ2 = hMatrix;
%}

meanQ1 = 0;
timeInterval = 2; % In seconds
% Estimation window
q2EstimationWindow = 100;
% Signal-to-noise ratios in decibels
signalToNoiseRatios = [30, 20, 10, 0];
% Our estimators
q2IdealEstimator = idealEstimator;
q2TestEstimator = testEstimator;
% HMatrices
q2IdealHMatrix = hMatrix;
q2TestHMatrix = hMatrix;
% Matrices to store the generated estimates
idealOutputsMatrixQ2B = cell(size(signalToNoiseRatios,2),q2EstimationWindow);
testOutputsMatrixQ2B = cell(size(signalToNoiseRatios,2),q2EstimationWindow);
% Generating windows
windowArray = genNArray(q2EstimationWindow);
% We will use both estimators and different snr to evaluate the impact of
% the variance of the white gaussian random process on the estimates
for snrIdnex = 1:1:size(signalToNoiseRatios, 2)
    % Getting a signal-to-noise ratio value
    currentSNR = signalToNoiseRatios(1, snrIdnex);
    % Getting the variance from this signal-t-noise ratio value
    currentVariance = genNoiseVar(currentSNR);
    % Index to add outputs from estimations to a matrix
    outputsIndex = 1;
    thetaQ2B = [0.25, pi/10];
    idealThetaQ2B = thetaQ2B;
    testThetaQ2B = thetaQ2B;
    hMatrixThetaQ2 = hMatrix;
    % Then we will use all of the windows to run the models
    for windowIndex = 1:1:size(windowArray,2)
        % Getting a window value
        currentWindowValue = windowArray(1, windowIndex);
        % Generating the noise free carrier for this windowValue
        noiseFreeCarrier = genNoiseFreeCarrier(currentWindowValue, thetaQ2B);
        % Generating the noisy signal carrier for this window value
        noisySignalCarrier = genNoisyCarrier(currentVariance, meanQ1, noiseFreeCarrier, timeInterval);
        % Updating the two estimators
        q2IdealEstimator = updateIdealEstimator(q2IdealEstimator, idealThetaQ2B);
        q2IdealHMatrix = updateHMatrix(q2IdealHMatrix, idealThetaQ2B);
        q2TestEstimator = updateTestEstimator(q2TestEstimator, currentWindowValue);
        q2TestHMatrix = updateHMatrix(q2TestHMatrix, testThetaQ2B);
        % Getting the estimated values of theta
        idealThetaQ2B = estimateThetaValue(q2IdealEstimator, q2IdealHMatrix, noisySignalCarrier, noiseFreeCarrier, idealThetaQ2B);
        testThetaQ2B = estimateThetaValue(q2TestEstimator, q2TestHMatrix, noisySignalCarrier, noiseFreeCarrier, testThetaQ2B);
        % We need to store these values from one iteration to the next
        idealOutputsMatrixQ2B{snrIdnex, outputsIndex} = idealThetaQ2B;
        testOutputsMatrixQ2B{snrIdnex, outputsIndex} = testThetaQ2B;
        % Increasing the index to store the next value
        outputsIndex = outputsIndex + 1;
    end
end


%% COLLECTING RANDOM SAMPLES FOR 2 SECONDS WITH RANDOM SEQUENCE OF LENGTH N=50, AND SIGNAL TO NOISE RATIOS
%% DECAYING IN THE ORDER: 30dB,20dB,10dB to 0dB THE TRUE FREQUENCY f(0) = 0.25Hz AND phi = pi/10
% For this question instead of getting our samples uniformly within the 
% observation window we will grab our samples randomly
% Generating noisy process and noisy signal
% To do this we will modify the function that gathers the samples to 
% randomize this process.

%{
function noisySig = genNoisyCarrierRandom(noiseVar, mean, noiseFreeSig, totalTime)
    % Underlying fomular: x(n) = s(n) + x'(n)
    % where x'(n) is a white gaussian random process
    % Array to hold the sequence of iid gaussian random variables
    sequence = NaN(1,1);
    % Using 1 millisecond sampling observations
    times = 0.001:0.001:totalTime;
    % Getting the standard deviation for our random process
    standDev = sqrt(noiseVar);
    % Pointer to add generated rv to the array
    index = 1;
    % Generating the actual sequence
    for time = [times(1), diff(times)]
        pause(time);
        % Gerating 1 RV at eacha millisecond
        randomVar = normrnd(mean, standDev);
        % Adding this rv to the array
        sequence(1, index) = randomVar;
        %  incrementing our pointer
        index = index + 1;
    end
    % Generating the random process using cumsum() function and getting 
    % a realization
    noisySig = noiseFreeSig + cumsum(sequence);    
end
%}

% Variables to be udes to run the simulations.
%{
thetaQ3 = [0.25, pi/10];
idealThetaQ3 = thetaQ3;
testThetaQ3 = thetaQ3;
hMatrixThetaQ3 = thetaQ3;
%}

meanQ1 = 0;
timeInterval = 2; % In seconds
% Estimation window
q3EstimationWindow = 50;
% Signal-to-noise ratios in decibels
signalToNoiseRatios = [30, 20, 10, 0];
% Our estimators
q3IdealEstimator = idealEstimator;
q3TestEstimator = testEstimator;
% HMatrices
q3IdealHMatrix = hMatrix;
q3TestHMatrix = hMatrix;
% Matrices to store the generated estimates
idealOutputsMatrixQ3 = cell(size(signalToNoiseRatios,2),q3EstimationWindow);
testOutputsMatrixQ3 = cell(size(signalToNoiseRatios,2),q3EstimationWindow);
% Generating windows
windowArray = genNArray(q3EstimationWindow);
% We will use both estimators and different snr to evaluate the impact of
% the variance of the white gaussian random process on the estimates
for snrIdnex = 1:1:size(signalToNoiseRatios, 2)
    % Getting a signal-to-noise ratio value
    currentSNR = signalToNoiseRatios(1, snrIdnex);
    % Getting the variance from this signal-t-noise ratio value
    currentVariance = genNoiseVar(currentSNR);
    % Index to add outputs from estimations to a matrix
    outputsIndex = 1;
    thetaQ3 = [0.25, pi/10];
    idealThetaQ3 = thetaQ3;
    testThetaQ3 = thetaQ3;
    hMatrixThetaQ3 = thetaQ3;
    % Then we will use all of the windows to run the models
    for windowIndex = 1:1:size(windowArray,2)
        % Getting a window value
        currentWindowValue = windowArray(1, windowIndex);
        % Generating the noise free carrier for this windowValue
        noiseFreeCarrier = genNoiseFreeCarrier(currentWindowValue, thetaQ3);
        % Generating the noisy signal carrier for this window value
        noisySignalCarrier = genNoisyCarrier(currentVariance, meanQ1, noiseFreeCarrier, timeInterval);
        % Updating the two estimators
        q3IdealEstimator = updateIdealEstimator(q3IdealEstimator, idealThetaQ3);
        q3IdealHMatrix = updateHMatrix(q3IdealHMatrix, idealThetaQ3);
        q3TestEstimator = updateTestEstimator(q3TestEstimator, currentWindowValue);
        q3TestHMatrix = updateHMatrix(q3TestHMatrix, testThetaQ3);
        % Getting the estimated values of theta
        idealThetaQ3 = estimateThetaValue(q3IdealEstimator, q3IdealHMatrix, noisySignalCarrier, noiseFreeCarrier, idealThetaQ3);
        testThetaQ3 = estimateThetaValue(q3TestEstimator, q3TestHMatrix, noisySignalCarrier, noiseFreeCarrier, testThetaQ3);
        % We need to store these values from one iteration to the next
        idealOutputsMatrixQ3{snrIdnex, outputsIndex} = idealThetaQ3;
        testOutputsMatrixQ3{snrIdnex, outputsIndex} = testThetaQ3;
        % Increasing the index to store the next value
        outputsIndex = outputsIndex + 1;
    end
end
%}