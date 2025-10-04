%% Stochastic Aerospace Systems Practical: Henry Adam

%% Part 1: Import Equations of Motion and Conduct Stability Analysis
%Using Asymetric Equation of Motion from Flight Dynamics Notes

%First Load Cesna Variables, Including A, B
CesnaVariables; 
stateNames = ["\beta" "\phi" "p" "r" "u_{g}" "u_{g}^*" "\alpha_{g}" "\alpha_{g}^*" "\beta_{g}" "\beta_{g}^*" "a_{y}"]; 
units = ["[deg]" "[deg]" "[deg/s]" "[deg/s]" "[deg]" "[-]" "[deg]" "[-]" "[deg]" "[-]" "[m/s^2]"]; 
fileNames = ["beta" "phi" "pb_2V" "rb_2V" "u_g" "u_g_star" "alpha_g" "alpha_g_star" "beta_g" "beta_g_star" "a_y"]; 
% 

%To Include in Report: A and B values 
%NOTE: % 
%x = [beta phi pb/2V rb/2V ug_ u_g* alpha_g alpha_g* beta_g betag*]'
%u = [delta_a delta_r w1 w3 w2]'.
A
B

%Checking Stability by looking at Eigenvalues 
eig(A) 
% To Include in Report: Note that there is a positive eigenvalue, so spiral
% mode is not stable
C = eye(size(A, 1)); 
%Create Controller for phi and p to stabilize the unstable system
% State feedback control
reoptimize = false; 
if reoptimize
    damping_ratio = 0; 
    K_phi_optimal = 0;
    K_p_optimal = 0;
    
    for K_phi = linspace(-5,5,1000)
        for K_p = linspace(-5,5,1000)
            K    = [0 K_phi K_p 0  0 0  0 0  0 0];
            A1   = A-B(:,1)*K;
            [~,zeta,poles] = damp(A1);
            obj = min(zeta); 
            if (obj > damping_ratio) && (all(real(poles) < 0))
                damping_ratio = min(zeta);
                K_phi_optimal = K_phi;
                K_p_optimal = K_p;
            end
        end
    end
    disp('Optimal K_phi:'); disp(K_phi_optimal); disp('');
    disp('Optimal K_p:'); disp(K_p_optimal); disp('');
    disp('State feedback damping ratio:'); disp(damping_ratio); disp('');
end


Kphi =  -0.4154; 
Kp = 0.1151; 
K = [0 Kphi Kp 0  0 0  0 0  0 0];

%Using negative feedback to stabilize the 
A_stable = A - B(:, 1)*K; 

%Check the Eigenvalues to see if system is stable or not
eig(A_stable);

%New system is stable, so include the new augmented A matrix in the report
A_stable;
%damp(A_stable)

%% Part 2: Time Domain Analysis

%First Augment the A matrix with the lateral velocity term
%NOTE: CHANGE IN STATES
%x = [beta phi pb/2V rb/2V ug_ u_g* alpha_g alpha_g* beta_g betag* a_y]'
%u = [delta_a delta_r w1 w3 w2]'.
g = 9.81; 

disp(A_stable)
%B = [B; V * (B(1, :) + B(2, :))]; 
%Set time 
dt = 0.01; totTime = 120; t = linspace(0, totTime, totTime/dt); N = length(t);

% Describe the turbulence as white noise and scale by sqrt(dt) to adjust for lsim
%COMMENTED OUT TO MAKE SURE I DON"T RESIMULATE THE SYSTEM AND CHANGE THE
%DRIVING NOISE
u_g = randn(N, 1)/sqrt(dt); %Longitudinal Gust Velocity(Horizontal Gust Velocity)
w_g = randn(N, 1)/sqrt(dt); %Vertical Gust Velocity

%For Filling in Input Matrices
zer = zeros(N,1); 

%Create Inputs: Make sure to simulate these separately
u1 = [zer zer u_g zer  zer]; %ug input    
u2 = [zer zer zer  w_g zer]; %wg input   

%Create Output Matrices [V * (r + beta_dot)] 
augmented_line = V * (A(1, :) + [0, 0, 0, 1, 0,0,0,0,0,0]);

A_stable = [[A_stable; augmented_line], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]']; 
B = [B; [0,0,0,0,0]]; 
C = eye(size(A_stable, 1)); % In order to get all states
rate_scaler = V * 2 / b * 180/pi; 
angle_scaler = 180/pi; 
C(3,3) = rate_scaler; C(4,4) = rate_scaler; %In order to get actual rates in degrees/s instead of normalized rates
C(1,1) = angle_scaler; C(2,2) = angle_scaler; C(5,5) = angle_scaler; C(7,7) = angle_scaler; C(9,9) = angle_scaler; 
D = zeros(size(C, 1), size(B, 2)); % Causal System
%Simulate System
% u_g
y1 = lsim(A_stable,B,C,D,u1,t);
% w_g
y2 = lsim(A_stable,B,C,D,u2,t);

% NOTE: Make sure to full-screen these plots before saving them! 
%Plotting For u_g response
f = figure(1);
for i=1:size(y1, 2)
    if i == 11
        subplot(3, 4, [10, 11])
        plot(t, y1(:, i))
        xlabel("time [s]", FontSize=18)
        ylabel(stateNames(i) + " " + units(i), FontSize=18)
        %title("Response of " + stateNames(i) +" for w_{g}")
        ylim([-1,1])
    elseif i <= 4
        subplot(3, 4, [2 * (i-1) + 1, 2 * (i-1) + 2])
        plot(t, y1(:, i))
        xlabel("time [s]", FontSize=18)
        ylabel(stateNames(i) + " " + units(i), FontSize=18)
    else
        continue
    end
    sgtitle("Response of all Non-Gust States to Horizontal Input u_{g}", FontSize=16, FontWeight='bold')

end
exportgraphics(f, "horizontal_time_reponse.pdf") 

%Plotting For w_g response
f = figure(2);
for i=1:size(y1, 2)
    if i == 11
        subplot(3, 4, [10, 11])
        plot(t, y2(:, i))
        xlabel("time [s]", FontSize=18)
        ylabel(stateNames(i) + " " + units(i), FontSize=18)
        %title("Response of " + stateNames(i) +" for w_{g}")
        ylim([-1,1])
    elseif i <= 4
        subplot(3, 4, [2 * (i-1) + 1, 2 * (i-1) + 2])
        plot(t, y2(:, i))
        xlabel("time [s]", FontSize=18)
        ylabel(stateNames(i) + " " + units(i), FontSize=18)
    else
        continue
    %title("Response of " + stateNames(i) +" for w_{g}")
    end
    sgtitle("Response of all Non-Gust States to Vertical Input u_{g}", FontSize=16, FontWeight='bold')

end


exportgraphics(f, "vertical_time_reponse.pdf") 

%% Part 3a: Setup for Frequency Domain Analysis
% Not for any of the the gust states, so the C matrix will be: 
% Going from: ["beta" "phi" "pb_2V" "rb_2V" "u_g" "u_g_star" "alpha_g" "alpha_g_star" "beta_g" "beta_g_star"]; 
%         To: ["beta" "phi" "pb_2V" "rb_2V"]; (all other states are gust-related) 
% If these need to change, all you need to change is the C matrix

%pre-define omega for plotting
omega = logspace(-1.5, 2, 300); 

%% Part 3b: Spectral Analysis for Horizontal Response 

%i) Analytical PSD(Can be obtained using the bode() function)
%Get each power-spectral density individually
%Use third input to isolate horizontal response
inter = bode(A_stable, B, C(1,:), D(1,:), 3, omega); S_beta_h = inter .* conj(inter);
inter = bode(A_stable, B, C(2,:), D(2,:), 3, omega); S_phi_h = inter .* conj(inter);
inter = bode(A_stable, B, C(3,:), D(3,:), 3, omega); S_pb_h = inter .* conj(inter);
inter = bode(A_stable, B, C(4,:), D(4,:), 3, omega);  S_rb_h = inter .* conj(inter);
inter = bode(A_stable, B, C(11,:), D(11,:), 3, omega);  S_vy_h = inter .* conj(inter);

%ii) Experimental PSD with lsim() and fft() functions
%First simulate the system with only horizontal response(u1)
y_horizontal = lsim(A_stable,B,C,D, u1 ,t)'; 

%Pull state variable responses from sim
beta = y_horizontal(1, :);
phi= y_horizontal(2, :);
pb_2V = y_horizontal(3, :);
rb_2V = y_horizontal(4, :);
lateral_acceleration = y_horizontal(11, :);

%Compute periodogram with fft(), remembering to multiply by dt to adjust for lsim
%scaling. Then get PSD w/ equation from notes Sxx = (1/T) * |I_x|^2
I_beta = dt * fft(beta); S_beta_exp_h = real((1/totTime) * I_beta .* conj(I_beta)); 
I_phi =  dt * fft(phi); S_phi_exp_h = real((1/totTime) * I_phi .* conj(I_phi)); 
I_pb_2V =  dt * fft(pb_2V); S_pb_exp_h = real((1/totTime) * I_pb_2V .* conj(I_pb_2V));  
I_rb_2V =  dt * fft(rb_2V); S_rb_exp_h = real((1/totTime) * I_rb_2V .* conj(I_rb_2V)); 
I_vy = dt * fft(lateral_acceleration); S_vy_exp_h = real((1/totTime) * I_vy .* conj(I_vy)); 

%getting the omega scale for comparative plotting                              
omegaExp_h = 2*pi*(1/dt)*(0:(N/2)-1)/N;
[~, idx] = min(abs(omegaExp_h - 1e2));
omegaExp_h = omegaExp_h(1:idx); 

%iii) Experimental PSD w/ Smoothing Function
%Smoothing function is just a convolution of the periodogram
conv_ker = [0.25, 0.5, 0.25]; %convolutional kernel

%Smooth periodograms and re-calculate the PSD
I_beta_smooth = conv(I_beta, conv_ker, "same"); S_beta_smooth_h = real((1/totTime) * I_beta_smooth .* conj(I_beta_smooth)); 
I_phi_smooth =  conv(I_phi, conv_ker, "same"); S_phi_smooth_h = real((1/totTime) * I_phi_smooth .* conj(I_phi_smooth)); 
I_pb_2V_smooth =  conv(I_pb_2V, conv_ker, "same"); S_pb_2V_smooth_h = real((1/totTime) * I_pb_2V_smooth .* conj(I_pb_2V_smooth)); 
I_rb_2V_smooth =  conv(I_rb_2V, conv_ker, "same"); S_rb_2V_smooth_h = real((1/totTime) * I_rb_2V_smooth .* conj(I_rb_2V_smooth)); 
I_vy_smooth =  conv(I_vy, conv_ker, "same"); S_vy_smooth_h = real((1/totTime) * I_vy_smooth .* conj(I_vy_smooth)); 

% Graph all of the methods for generating the PSD: 
%Plot 1: Beta
figure()
loglog(omega, S_beta_h, ":", 'LineWidth',2)
hold on 
loglog(omegaExp_h, S_beta_exp_h(1:idx))
loglog(omegaExp_h, S_beta_smooth_h(1:idx))
hold off
xlim([0.05, max(omega)])
ylim([10^-14, 10^2])
ylabel("S_{\beta \beta}[deg^2*s]", FontSize=18)
xlabel("\omega [rad/s]", FontSize=18)
title("Auto-PSD for side-slip angle \beta due to Horizontal Gusts", FontSize=16, FontWeight='bold')
legend("Analytical", "Experimental: FFT", "Experimental: Smoothing", FontSize=16, Location='southwest')
exportgraphics(gca, "horizontal_Sbetabeta_comparison.pdf") 

%Plot 2: phi
figure()
loglog(omega, S_phi_h, ":", 'LineWidth',2)
hold on 
loglog(omegaExp_h, S_phi_exp_h(1:idx))
loglog(omegaExp_h, S_phi_smooth_h(1:idx))
hold off
ylabel("S_{\phi \phi} [deg^2*s]", FontSize=18)
xlabel("\omega [rad/s]", FontSize=18)
xlim([0.05, max(omega)])
ylim([10^-14, 10^2])
title("Auto-PSD for roll angle \phi due to Horizontal Gusts", FontSize=16, FontWeight='bold')
legend("Analytical", "Experimental: FFT", "Experimental: Smoothing", FontSize=16, Location='southwest')
exportgraphics(gca, "horizontal_Sphiphi_comparison.pdf") 

%Plot 3: pb
figure()
loglog(omega, S_pb_h, ":", 'LineWidth',2)
hold on 
loglog(omegaExp_h, S_pb_exp_h(1:idx))
loglog(omegaExp_h, S_pb_2V_smooth_h(1:idx))
hold off
xlim([0.05, max(omega)])
ylim([10^-14, 10^2])
ylabel("S_{pp} [deg^2/s]", FontSize=18)
xlabel("\omega [rad/s]", FontSize=18)
title("Auto-PSD for normalized roll rate pb/2V due to Horizontal Gusts", FontSize=16, FontWeight='bold')
legend("Analytical", "Experimental: FFT", "Experimental: Smoothing", FontSize=16, Location='southwest')
exportgraphics(gca, "horizontal_Spbpb_comparison.pdf") 

%Plot 4: rb
figure()
loglog(omega, S_rb_h, ":", 'LineWidth',2)
hold on 
loglog(omegaExp_h, S_rb_exp_h(1:idx))
loglog(omegaExp_h, S_rb_2V_smooth_h(1:idx))
hold off
ylabel("S_{rr} [deg^2/s]", FontSize=18)
xlabel("\omega [rad/s]", FontSize=18)
xlim([0.05, max(omega)])
ylim([10^-14, 10^2])
title("Auto-PSD for normalized yaw rate rb/2V due to Horizontal Gusts", FontSize=16, FontWeight='bold')
legend("Analytical", "Experimental: FFT", "Experimental: Smoothing", FontSize=16, Location='southwest')
exportgraphics(gca, "horizontal_Srbrb_comparison.pdf") 

%Plot 4: a_y
figure()
loglog(omega, S_vy_h, ":", 'LineWidth',2)
hold on 
loglog(omegaExp_h, S_vy_exp_h(1:idx))
loglog(omegaExp_h, S_vy_smooth_h(1:idx))
hold off
xlim([0.05, max(omega)])
ylim([10^-14, 10^2])
ylabel("S_{a_{y}a_{y}} [m^2/s]", FontSize=18)
xlabel("\omega [rad/s]", FontSize=18)
title("Auto-PSD for lateral acceleration due to Horizontal Gusts", FontSize=16, FontWeight='bold')
legend("Analytical", "Experimental: FFT", "Experimental: Smoothing", FontSize=16, Location='southwest')
exportgraphics(gca, "horizontal_Svyvy_comparison.pdf") 

%% Part 3c: Spectral Analysis for Vertical Response 

%i) Analytical PSD(Can be obtained using the bode() function)
%Get each power-spectral density individually(use 4th input to isolate veritcal terbulence)
inter = bode(A_stable, B, C(1,:), D(1,:), 4, omega); S_beta_v = inter .* conj(inter);
inter = bode(A_stable, B, C(2,:), D(2,:), 4, omega); S_phi_v = inter .* conj(inter);
inter = bode(A_stable, B, C(3,:), D(3,:), 4, omega); S_pb_v = inter .* conj(inter);
inter = bode(A_stable, B, C(4,:), D(4,:), 4, omega);  S_rb_v = inter .* conj(inter);
inter = bode(A_stable, B, C(11,:), D(5,:), 4, omega);  S_vy_v = inter .* conj(inter);

%ii) Experimental PSD with lsim() and fft() functions
%First simulate the system with only vertical response(u2) 
y_vertical = lsim(A_stable,B,C,D,u2,t)'; 

%Pull state variable responses from sim
beta = y_vertical(1, :);
phi= y_vertical(2, :);
pb_2V = y_vertical(3, :);
rb_2V = y_vertical(4, :);
lateral_acceleration = y_vertical(11, :);

%Compute periodogram with fft(), remembering to multiply by dt to adjust for lsim
%scaling. Then get PSD w/ equation from notes Sxx = (1/T) * |I_x|^2
I_beta = dt * fft(beta); S_beta_exp_v = real((1/totTime) * I_beta .* conj(I_beta)); 
I_phi =  dt * fft(phi); S_phi_exp_v = real((1/totTime) * I_phi .* conj(I_phi)); 
I_pb_2V =  dt * fft(pb_2V); S_pb_exp_v = real((1/totTime) * I_pb_2V .* conj(I_pb_2V));  
I_rb_2V =  dt * fft(rb_2V); S_rb_exp_v = real((1/totTime) * I_rb_2V .* conj(I_rb_2V)); 
I_vy = dt * fft(lateral_acceleration); S_vy_exp_v = real((1/totTime) * I_vy .* conj(I_vy)); 

%getting the omega scale for comparative plotting                              
omegaExp_v = 2*pi*(1/dt)*(0:(N/2)-1)/N; 
[~, idx] = min(abs(omegaExp_v - 1e2));
omegaExp_v = omegaExp_v(1:idx); 

%iii) Experimental PSD w/ Smoothing Function
%Smoothing function is just a convolution of the periodogram
conv_ker = [0.25, 0.5, 0.25]; %convolutional kernel

%Smooth periodograms and re-calculate the PSD
I_beta_smooth = conv(I_beta, conv_ker, "same"); S_beta_smooth_v = real((1/totTime) * I_beta_smooth .* conj(I_beta_smooth)); 
I_phi_smooth =  conv(I_phi, conv_ker, "same"); S_phi_smooth_v = real((1/totTime) * I_phi_smooth .* conj(I_phi_smooth)); 
I_pb_2V_smooth =  conv(I_pb_2V, conv_ker, "same"); S_pb_2V_smooth_v = real((1/totTime) * I_pb_2V_smooth .* conj(I_pb_2V_smooth)); 
I_rb_2V_smooth =  conv(I_rb_2V, conv_ker, "same"); S_rb_2V_smooth_v = real((1/totTime) * I_rb_2V_smooth .* conj(I_rb_2V_smooth)); 
I_vy_smooth =  conv(I_vy, conv_ker, "same"); S_vy_smooth_v = real((1/totTime) * I_vy_smooth .* conj(I_vy_smooth)); 

% Graph all of the methods for generating the PSD: 
%Plot 1: Beta
figure()
loglog(omega, S_beta_v, ":", 'LineWidth',2)
hold on 
loglog(omegaExp_v, S_beta_exp_v(1:idx))
loglog(omegaExp_v, S_beta_smooth_v(1:idx))
hold off
xlim([0.05, max(omega)])
ylim([10^-14, 10^2])
ylabel("S_{\beta \beta}  [deg^2*s]", FontSize=18)
xlabel("\omega [rad/s]", FontSize=18)
title("Auto-PSD for side-slip angle \beta due to Vertical Gusts", FontSize=16, FontWeight='bold')
legend("Analytical", "Experimental: FFT", "Experimental: Smoothing", FontSize=16, Location='southwest')
exportgraphics(gca, "vertical_Sbetabeta_comparison.pdf") 

%Plot 2: phi
figure()
loglog(omega, S_phi_v, ":", 'LineWidth',2)
hold on 
loglog(omegaExp_v, S_phi_exp_v(1:idx))
loglog(omegaExp_v, S_phi_smooth_v(1:idx))
hold off
xlim([0.05, max(omega)])
ylim([10^-14, 10^2])
ylabel("S_{\phi \phi}  [deg^2*s]", FontSize=18)
xlabel("\omega [rad/s]", FontSize=18)
title("Auto-PSD for roll angle \phi due to Vertical Gusts", FontSize=16, FontWeight='bold')
legend("Analytical", "Experimental: FFT", "Experimental: Smoothing", FontSize=16, Location='southwest')
exportgraphics(gca, "vertical_Sphiphi_comparison.pdf") 

%Plot 3: pb
figure()
loglog(omega, S_pb_v, ":", 'LineWidth',2)
hold on 
loglog(omegaExp_v, S_pb_exp_v(1:idx))
loglog(omegaExp_v, S_pb_2V_smooth_v(1:idx))
hold off
xlim([0.05, max(omega)])
ylim([10^-14, 10^2])
ylabel("S_{pp} [deg^2/s]", FontSize=18)
xlabel("\omega [rad/s]", FontSize=18)
title("Auto-PSD for normalized roll rate pb/2V due to Vertical Gusts", FontSize=16, FontWeight='bold')
legend("Analytical", "Experimental: FFT", "Experimental: Smoothing", FontSize=16, Location='southwest')
exportgraphics(gca, "vertical_Spbpb_comparison.pdf") 

%Plot 4: rb
figure()
loglog(omega, S_rb_v, ":", 'LineWidth',2)
hold on 
loglog(omegaExp_v, S_rb_exp_v(1:idx))
loglog(omegaExp_v, S_rb_2V_smooth_v(1:idx))
hold off
xlim([0.05, max(omega)])
ylim([10^-14, 10^2])
ylabel("S_{rr} [deg^2/s]", FontSize=18)
xlabel("\omega [rad/s]", FontSize=18)
title("Auto-PSD for normalized yaw rate rb/2V due to Vertical Gusts", FontSize=16, FontWeight='bold')
legend("Analytical", "Experimental: FFT", "Experimental: Smoothing", FontSize=16, Location='southwest')
exportgraphics(gca, "vertical_Srbrb_comparison.pdf") 

%Plot 5: a_y
figure()
loglog(omega, S_vy_v, ":", 'LineWidth',2)
hold on 
loglog(omegaExp_v, S_vy_exp_v(1:idx))
loglog(omegaExp_v, S_vy_smooth_v(1:idx))
hold off
xlim([0.05, max(omega)])
ylim([10^-14, 10^2])
ylabel("S_{a_{y}a_{y}} [m^2/s]", FontSize=18)
xlabel("\omega [rad/s]", FontSize=18)
title("Auto-PSD for lateral acceleration due to Vertical Gusts", FontSize=16, FontWeight='bold')
legend("Analytical", "Experimental: FFT", "Experimental: Smoothing", FontSize=16, Location='southwest')
exportgraphics(gca, "vertical_Svyvy_comparison.pdf") 

%% Finding the bounds of the functions: 
max_b1_rad = max(abs(y_horizontal(1, :))); 
max_p1_rad = max(abs(y_horizontal(2, :))); 
max_pb1_weird = max(abs(y_horizontal(3, :)));
max_rb1_weird = max(abs(y_horizontal(4, :))); 
max_ay1 = max(abs(y_horizontal(11, :)));
trad_units_horiz = [max_b1_rad, max_p1_rad, max_pb1_weird, max_rb1_weird, max_ay1]'; 

max_b2_rad = max(abs(y_vertical(1, :))); 
max_p2_rad = max(abs(y_vertical(2, :))); 
max_pb2_weird = max(abs(y_vertical(3, :))); 
max_rb2_weird = max(abs(y_vertical(4, :))); 
max_ay2 = max(abs(y_vertical(11, :)));

trad_units_vert = [max_b2_rad, max_p2_rad, max_pb2_weird, max_rb2_weird, max_ay2]'; 

max_vals_table = table(["beta" "phi" "p" "r" "a_y"]', trad_units_horiz, trad_units_vert)

%% Part 4a: Setup for Variances 
%Comparing Various Methods of obtaining the variance of the system
%get the differences to numerically integrate
dw = diff(omega); 

%get state names for making tables
SV_Name = ["beta" "phi" "pb_2V" "rb_2V" "a_y" "u_g" "u_g_star" "alpha_g" "alpha_g_star" "beta_g" "beta_g_star"]'; 

%% Part 4b: Horizontal Analysis of Variances
% FOR STATES WE'VE ALREADY BEEN USING:
%Method 1: Using the analytical PSD(integrate across psd using trapezoidal integration) 
var_beta_analytical_h = trapz(omega, S_beta_h) / (pi); 
var_phi_analytical_h = trapz(omega, S_phi_h) / (pi);
var_pb_analytical_h = trapz(omega, S_pb_h) / (pi);
var_rb_analytical_h = trapz(omega, S_rb_h) / (pi);
var_vy_analytical_h = trapz(omega, S_vy_h) / (pi);

%Method 2: Using the derrived PSDs
%First do auto-psds found with fft() 
var_beta_fft_h = trapz(omegaExp_h, S_beta_exp_h(1:idx)) / (pi); 
var_phi_fft_h = trapz(omegaExp_h, S_phi_exp_h(1:idx)) / (pi);
var_pb_fft_h = trapz(omegaExp_h, S_pb_exp_h(1:idx)) / (pi);
var_rb_fft_h = trapz(omegaExp_h, S_rb_exp_h(1:idx)) / (pi);
var_vy_fft_h = trapz(omegaExp_h, S_vy_exp_h(1:idx)) / (pi);

%Next do smoothed auto-psds
var_beta_smooth_h = trapz(omegaExp_h, S_beta_smooth_h(1:idx)) / (pi); 
var_phi_smooth_h = trapz(omegaExp_h, S_phi_smooth_h(1:idx)) / (pi);
var_pb_smooth_h = trapz(omegaExp_h, S_pb_2V_smooth_h(1:idx)) / (pi);
var_rb_smooth_h = trapz(omegaExp_h, S_rb_2V_smooth_h(1:idx)) / (pi);
var_vy_smooth_h = trapz(omegaExp_h, S_vy_smooth_h(1:idx)) / (pi);

%Method 3: Using the var() function on the simulated system
var_beta_function_h = var(y_horizontal(1, :)); 
var_phi_function_h = var(y_horizontal(2, :)); 
var_pb_function_h = var(y_horizontal(3, :)); 
var_rb_function_h = var(y_horizontal(4, :)); 
var_vy_function_h = var(y_horizontal(11, :)); 

%Put what we HAVE been using all together into a table
var_analytical_h = [var_beta_analytical_h var_phi_analytical_h var_pb_analytical_h var_rb_analytical_h var_vy_analytical_h 0 0 0 0 0 0]';
var_fft_h = [var_beta_fft_h var_phi_fft_h var_pb_fft_h var_rb_fft_h var_vy_fft_h 0 0 0 0 0 0]';
var_smooth_h = [var_beta_smooth_h var_phi_smooth_h var_pb_smooth_h var_rb_smooth_h var_vy_smooth_h 0 0 0 0 0 0]';
var_function_h = [var_beta_function_h var_phi_function_h var_pb_function_h var_rb_function_h var_vy_function_h 0 0 0 0 0 0]';


% FOR STATES WE HAVEN'T BEEN USING: 
for i = 5:10
    %Get analytical values for table
    inter = bode(A_stable, B, C(i,:), D(i,:), 3, omega);  S_i_h = inter .* conj(inter);
    var_i_analytical_h = trapz(omega, S_i_h) / pi;
    var_analytical_h(i+1) = var_i_analytical_h; 

    %Get fft variance 
    i_time_response = y_horizontal(i, :);
    I_i = dt * fft(i_time_response); S_i_exp_h = real((1/totTime) * I_i .* conj(I_i)); 
    var_exp_i = trapz(omegaExp_h, S_i_exp_h(1:idx)) / (pi); 
    var_fft_h(i+1) = var_exp_i; 

    %Get smoothed variance
    I_i_smooth = conv(I_i, conv_ker, "same"); S_i_smooth_h = real((1/totTime) * I_i_smooth .* conj(I_i_smooth)); 
    var_smooth_i = trapz(omegaExp_h, S_i_smooth_h (1:idx)) / (pi);
    var_smooth_h(i+1) = var_smooth_i;

    %Get from fucntion
    var_function_h(i+1) = var(i_time_response);
end

horizontal_variance_table = table(SV_Name, var_analytical_h, var_fft_h, var_smooth_h, var_function_h)
%% Part 4c: Vertical Analysis of Variances
%Method 1: Using the analytical PSD(integrate across psd using trapezoidal integration) 
var_beta_analytical_v = trapz(omega, S_beta_v) / (pi); 
var_phi_analytical_v = trapz(omega, S_phi_v) / (pi);
var_pb_analytical_v = trapz(omega, S_pb_v) / (pi);
var_rb_analytical_v = trapz(omega, S_rb_v) / (pi);
var_vy_analytical_v = trapz(omega, S_vy_v) / (pi);

%Method 2: Using the derrived PSDs
%First do fft() auto-psds 
var_beta_fft_v = trapz(omegaExp_v, S_beta_exp_v(1:idx)) / (pi); 
var_phi_fft_v = trapz(omegaExp_v, S_phi_exp_v(1:idx)) / (pi);
var_pb_fft_v = trapz(omegaExp_v, S_pb_exp_v(1:idx)) / (pi);
var_rb_fft_v = trapz(omegaExp_v, S_rb_exp_v(1:idx)) / (pi);
var_vy_fft_v = trapz(omegaExp_v, S_vy_exp_v(1:idx)) / (pi);

%Next do smoothed auto-psds
var_beta_smooth_v = trapz(omegaExp_v, S_beta_smooth_v(1:idx)) / (pi); 
var_phi_smooth_v = trapz(omegaExp_v, S_phi_smooth_v(1:idx)) / (pi);
var_pb_smooth_v = trapz(omegaExp_v, S_pb_2V_smooth_v(1:idx)) / (pi);
var_rb_smooth_v = trapz(omegaExp_v, S_rb_2V_smooth_v(1:idx)) / (pi);
var_vy_smooth_v = trapz(omegaExp_v, S_vy_smooth_v(1:idx)) / (pi);

%Method 3: Using the var() function on the simulated system
var_beta_function_v = var(y_vertical(1, :)); 
var_phi_function_v = var(y_vertical(2, :)); 
var_pb_function_v = var(y_vertical(3, :)); 
var_rb_function_v = var(y_vertical(4, :)); 
var_vy_function_v = var(y_vertical(11, :)); 

%Put it all together into a table
var_analytical_v = [var_beta_analytical_v var_phi_analytical_v var_pb_analytical_v var_rb_analytical_v var_vy_analytical_v 0 0 0 0 0 0]';
var_fft_v = [var_beta_fft_v var_phi_fft_v var_pb_fft_v var_rb_fft_v var_vy_fft_v 0 0 0 0 0 0]';
var_smooth_v = [var_beta_smooth_v var_phi_smooth_v var_pb_smooth_v var_rb_smooth_v var_vy_smooth_v 0 0 0 0 0 0]';
var_function_v = [var_beta_function_v var_phi_function_v var_pb_function_v var_rb_function_v var_vy_function_v 0 0 0 0 0 0]';

for i = 5:10
    %Get analytical values for table
    inter = bode(A_stable, B, C(i,:), D(i,:), 4, omega);  S_i_v = inter .* conj(inter);
    var_i_analytical_v = trapz(omega, S_i_v) / pi;
    var_analytical_v(i+1) = var_i_analytical_v; 

    %Get fft variance 
    i_time_response = y_vertical(i, :);
    I_i = dt * fft(i_time_response); S_i_exp_v = real((1/totTime) * I_i .* conj(I_i)); 
    var_exp_i = trapz(omegaExp_v, S_i_exp_v(1:idx)) / (pi); 
    var_fft_v(i+1) = var_exp_i; 

    %Get smoothed variance
    I_i_smooth = conv(I_i, conv_ker, "same"); S_i_smooth_v = real((1/totTime) * I_i_smooth .* conj(I_i_smooth)); 
    var_smooth_i = trapz(omegaExp_v, S_i_smooth_v (1:idx)) / (pi);
    var_smooth_v(i+1) = var_smooth_i;

    %Get from fucntion
    var_function_v(i+1) = var(i_time_response);
end
vertical_variance_table = table(SV_Name, var_analytical_v, var_fft_v, var_smooth_v, var_function_v)

%save("final_report")

vertical_variance_table.var_function_v ./ vertical_variance_table.var_analytical_v
horizontal_variance_table.var_function_h ./ horizontal_variance_table.var_analytical_h