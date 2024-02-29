% Parameters for the signals
fs = 360; % Sampling frequency (Hz)
t = 0:1/fs:10-1/fs; % Time vector for 10 seconds
f_sin = 1; % Frequency of sinusoidal signal
f_square = 50; % Frequency of square signal
f_triangular = 9; % Frequency of triangular signal
amplitude = max(clean_ecg)-min(clean_ecg); % Amplitude of all signals

% Generate individual signals
sinusoidal_signal = amplitude * sin(2 * pi * f_sin * t);
square_signal = amplitude * square(2 * pi * f_square * t);
triangular_signal = amplitude * sawtooth(2 * pi * f_triangular * t, 0.5);

% Combine signals
combined_reference_ecg = (sinusoidal_signal + square_signal + triangular_signal)+0.2;

SNR_VAL = 20;
ecg_with_combined_noise=awgn(clean_ecg,SNR_VAL,"measured");

% Initialize variables for LMS adaptive filtering
min_error = Inf; % Initialize minimum error as infinity
optimal_order = 16; % Start with the initial filter order
optimal_W = zeros(optimal_order, 1); % Initialize optimal filter coefficients
mu = 0.01; % LMS step size

% Define the maximum order to test
max_order = 96;
mse_array = zeros(1,max_order-optimal_order);
artificial_reference_ecg =  combined_reference_ecg;

% Loop over different filter orders to find the optimal one
for test_order = optimal_order:max_order
    % Initialize filter coefficients for this test order
    W_test = zeros(test_order, 1);
    
    % Apply LMS algorithm for the current order using the artificial reference
    for n = test_order+1:length(ecg_with_combined_noise)
        x = artificial_reference_ecg(n-test_order+1:n)'; % Input signal window from the artificial reference
        d = ecg_with_combined_noise(n); % Desired output (noisy ECG)
                                                                                                                                                                                                                           
        % Calculate predicted output
        y = W_test' * x;
        
        % Error signal
        e = d - y;
        
        % Update filter coefficients using LMS algorithm
        W_test = W_test + mu * e * x;
    end
    
    % Filter the signal using the current filter coefficients
    filtered_ecg_test = filter(W_test, 1, ecg_with_combined_noise);
    
    % Calculate the mean squared error (MSE) for the current filter
    mse_test = immse(artificial_reference_ecg, filtered_ecg_test);
    mse_array(test_order)=mse_test;

    % Check if the current MSE is lower than the minimum error found so far
    if mse_test < min_error
        min_error = mse_test; % Update minimum error
        optimal_order = test_order; % Update optimal filter order
        optimal_W = W_test; % Update optimal filter coefficients
    end
    
end

% Filter the entire signal using the obtained optimal filter coefficients
filtered_ecg_optimal = filter(optimal_W, 1, ecg_with_combined_noise);
MSEwithCleanEcg = immse(clean_ecg,filtered_ecg_optimal);


% Mean Square Error with Clean Ecg -----------------------------------
MSEwithCleanEcg = immse(clean_ecg,filtered_ecg_optimal);

% Display the optimal filter order and coefficients
fprintf('Optimal Filter Order: %d\n', optimal_order);
disp('Optimal Filter Coefficients:');
disp(optimal_W);
disp('MSE: ');
disp(MSEwithCleanEcg);


%Plots and Figures -------------------------------------------------
figure(3)
subplot(411)
plot(t,clean_ecg);
xlabel('Time (second)');
ylabel('Amplitude (V)');
title('Clean ECG Signal');
grid on;

subplot(412)
plot(t,ecg_with_combined_noise);
xlabel('Time (second)');
ylabel('Amplitude (V)');
title(['ECG with Combined Noise, SNR(dB): ', num2str(SNR_VAL)]);
grid on;

subplot(413)
plot(t,filtered_ecg_optimal);
xlabel('Time (second)');
ylabel('Amplitude (V)');
title(['Adaptive Filtered Signal, MSE: ', num2str(MSEwithCleanEcg)]);
grid on;

subplot(414)
plot(t,artificial_reference_ecg);
xlabel('Time (second)');
ylabel('Amplitude (V)');
title(['Artificial Reference ECG, Amplitude:',num2str(amplitude) 'V, Frequencies: Square:',num2str(f_square),'Hz, Sine: ',num2str(f_sin),'Hz, Triangular: ',num2str(f_triangular),'Hz']);
grid on;
%Ploting for MSE versus SNR
SNR_Array = [0,5,10,20,30,50];
COM_REF_1950 = [0.0015828, 0.0011091, 0.00064659, 0.00037343, 0.00033705, 0.00033502];
COM_REF_1912 = [0.0017463, 0.001189, 0.00070117, 0.0003636, 0.00034162, 0.0003345];
SIN_REF_9 = [0.072813, 0.046682, 0.053859, 0.04701, 0.050125, 0.050156];
SIN_REF_1 = [0.14277, 0.15842, 0.14234, 0.14295, 0.14041, 0.1404];
TRI_REF_9 = [0.051394, 0.05368, 0.044411, 0.050053, 0.048272, 0.048432];
REC_REF_50 = [0.0047658, 0.0020459, 0.0012371, 0.00075284, 0.00064636, 0.00063074];
CLEAN_REF = [0.0016421, 0.00095514, 0.000637, 0.00038091, 0.00033778, 0.00033398];

%plot graphs 

figure(1);
plot(SNR_Array,COM_REF_1950,LineWidth=1,Color="b",Marker="o");
xlabel("SNR(dB)");
ylabel('Mean Square Error (MSE)');
title("MSE of Combined Reference");
grid on;
hold on
plot(SNR_Array,COM_REF_1912,LineWidth=1,Color="r",Marker="o");
xlabel("SNR(dB)");
ylabel('Mean Square Error (MSE)');
legend('Freq Sin:1Hz Freq Tri: 9Hz Freq Rec:50Hz','Freq Sin:1Hz Freq Tri: 9Hz Freq Rec:12Hz');
hold off

figure(2);
hold on
plot(SNR_Array,REC_REF_50,LineWidth=1,Color="b",Marker="o");
plot(SNR_Array,SIN_REF_9,LineWidth=1,Color="r",Marker="o");
plot(SNR_Array,TRI_REF_9,LineWidth=1,Color="g",Marker="o");
hold off
title('MSE Comparison of Different Reference Signals');
xlabel("SNR(dB)");
ylabel('Mean Square Error (MSE)');
legend('Rectengular Reference, Frequency: 50Hz','Sinusoidal Reference, Frequency: 9Hz','Triangular Reference, Frequency: 9Hz');
grid on;

figure(4);
hold on
plot(SNR_Array,REC_REF_50,LineWidth=1,Color="b",Marker="o");
plot(SNR_Array,SIN_REF_9,LineWidth=1,Color="r",Marker="o");
plot(SNR_Array,TRI_REF_9,LineWidth=1,Color="g",Marker="o");
plot(SNR_Array,COM_REF_1950,LineWidth=1,Color="m",Marker="o");
hold off
title('MSE Comparison of Different Reference Signals');
xlabel("SNR(dB)");
ylabel('Mean Square Error (MSE)');
legend('Rectengular Reference, Frequency: 50Hz','Sinusoidal Reference, Frequency: 9Hz','Triangular Reference, Frequency: 9Hz','Combined Reference, Freq Sin:1Hz Freq Tri: 9Hz Freq Rec:50Hz');
grid on;

figure(5);
hold on;
plot(SNR_Array,REC_REF_50,LineWidth=1,Color="b",Marker="o");
plot(SNR_Array,CLEAN_REF,LineWidth=1,Color="r",Marker="o");
hold off;
title('Clean ECG as Reference Signal and Combined Reference Signal');
xlabel("SNR(dB)");
ylabel('Mean Square Error (MSE)');
legend('Combined Reference, Freq Sin:1Hz Freq Tri: 9Hz Freq Rec:50Hz', 'Clean ECG Itself for Reference Signal');
grid on;



