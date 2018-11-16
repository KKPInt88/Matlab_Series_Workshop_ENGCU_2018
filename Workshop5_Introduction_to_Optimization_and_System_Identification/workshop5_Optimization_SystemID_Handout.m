%% Workshop 5: Optimization in Design and Intro to System Identification
% Kan Kanjanapas (Ph.D.)

clc;
close all;
clear all;


%% Optimization Part 1: Linear Programming

% linprog
% min f'*x such that  A_ineq*x <= b_ineq, A_eq*x = b_eq, lb <= x <= ub
%  x

% 
% objective function -x1 + 2*x2 --> f = [-1 2]'
% inequality constraint:  x1 +   x2 <= 10
%                         x1 - 2*x2 <= 15
%                      -3*x1 - 5*x2 <= 20
%                      -4*x1 +   x2 <= 14

% A_ineq = [ 1  1]
%          [ 1 -2]
%          [-3 -5]
%          [-4 +1]

% B_ineq = [10 15 20 14]'

% objective function
f = [-1 2];

% Inequality constraints
A_ineq = [];
B_ineq = []';

% Equality constraints
A_eq = [];
B_eq = [];

% Boundary conditions: lower/upper bounds
lb = [];
ub = [];

% options



% Solver
% x = linprog(f, A_ineq, B_ineq)


% Double check solution
x1_range = [-10 : 0.25 : 15]';
x2_range = [-10 : 0.25 : 12]';
% inequality constraint:  x1 +   x2 <= 10
%                         x1 - 2*x2 <= 15
%                      -3*x1 - 5*x2 <= 20
%                      -4*x1 +   x2 <= 14

% % % x2_line1 = 1/+1*(-1*x1_range + 10);
% % % x2_line2 = 1/-2*(-1*x1_range + 15);
% % % x2_line3 = 1/-5*(+3*x1_range + 20);
% % % x2_line4 = 1/+1*(+4*x1_range + 14);



% % % figure;
% % % set(gcf,'Position',[0 0 1280 1280]);
% % % for ii = 1:1
% % %     
% % %     subplot(2,1,1);
% % %     surf(x1_mesh, x2_mesh, z, 'EdgeColor', 'None');
% % %     hold on; 
% % %     for ii = 1:length(x1_range)
% % %         for jj = 1:length(x2_range)
% % %             
% % %             x1_ij = x1_range(ii);
% % %             x2_ij = x2_range(jj);
% % %             
% % %             ind_1 = ( 1*x1_ij + 1*x2_ij <= 10);
% % %             ind_2 = ( 1*x1_ij - 2*x2_ij <= 15);
% % %             ind_3 = (-3*x1_ij - 5*x2_ij <= 20);
% % %             ind_4 = (-4*x1_ij + 1*x2_ij <= 14);
% % %             
% % %             if (ind_1*ind_2*ind_3*ind_4 == 1) % Satisfy all inequality constraints
% % %                 plot3( x1_ij, x2_ij, f_obj(x1_ij, x2_ij), 'LineStyle', 'None', ...
% % %                      'Marker', 'o', 'MarkerSize', 5, 'Color', 'c', 'MarkerFaceColor', 'c');
% % %             end
% % %             
% % %         end % for jj = 1:length(x2_range)
% % %     end % for ii = 1:length(x1_range)  
% % %     plot3( x1_opt, x2_opt, f_obj(x1_opt, x2_opt), 'LineStyle', 'None', ...
% % %                      'Marker', 'o', 'MarkerSize', 15, 'Color', 'r', 'MarkerFaceColor', 'r');
% % %     xlabel('x1');
% % %     ylabel('x2');
% % %     zlabel('Objective Function Value');
% % %     title('Ex1: Linear Programming');
% % %     grid on;
% % %     set(gca, 'FontSize', 16);
% % %     
% % %     
% % %     
% % %     
% % %     subplot(2,1,2);
% % %     plot(x1_range, x2_line1, 'LineWidth', 2, 'Color', 'b');
% % %     hold on;
% % %     plot(x1_range, x2_line2, 'LineWidth', 2, 'Color', 'r');
% % %     plot(x1_range, x2_line3, 'LineWidth', 2, 'Color', [0 0.5 0]);
% % %     plot(x1_range, x2_line4, 'LineWidth', 2, 'Color', [1 0.5 0]);
% % %     
% % %     for ii = 1:length(x1_range)
% % %         for jj = 1:length(x2_range)
% % %             
% % %             x1_ij = x1_range(ii);
% % %             x2_ij = x2_range(jj);
% % %             
% % %             ind_1 = ( 1*x1_ij + 1*x2_ij <= 10);
% % %             ind_2 = ( 1*x1_ij - 2*x2_ij <= 15);
% % %             ind_3 = (-3*x1_ij - 5*x2_ij <= 20);
% % %             ind_4 = (-4*x1_ij + 1*x2_ij <= 14);
% % %             
% % %             if (ind_1*ind_2*ind_3*ind_4 == 1) % Satisfy all inequality constraints
% % %                 plot( x1_ij, x2_ij, 'LineStyle', 'None', ...
% % %                      'Marker', 'o', 'MarkerSize', 5, 'Color', 'c', 'MarkerFaceColor', 'c');
% % %             end
% % %             
% % %         end % for jj = 1:length(x2_range)
% % %     end % for ii = 1:length(x1_range)
% % %     
% % %     plot( x1_opt, x2_opt, 'LineStyle', 'None', ...
% % %                      'Marker', 'o', 'MarkerSize', 15, 'Color', 'r', 'MarkerFaceColor', 'r');
% % %     
% % %     
% % %     plot(x1_range, x2_line1, 'LineWidth', 2, 'Color', 'b');
% % %     plot(x1_range, x2_line2, 'LineWidth', 2, 'Color', 'r');
% % %     plot(x1_range, x2_line3, 'LineWidth', 2, 'Color', [0 0.5 0]);
% % %     plot(x1_range, x2_line4, 'LineWidth', 2, 'Color', [1 0.5 0]);
% % %     
% % %     xlabel('x1');
% % %     ylabel('x2');
% % %     grid on;
% % %     h_legend = legend('x1 + x2 <= 10', ...
% % %                       'x1 - 2*x2 <= 15', ...
% % %                       '-3*x1 - 5*x2 <= 20', ...
% % %                       '-4*x1 + x2 <= 14');
% % %     set(h_legend, 'Location', 'NorthEast', 'Color', [1 0.9 1]);
% % %     axis([min(x1_range) max(x1_range) min(x2_range) max(x2_range)]);
% % %     set(gca, 'FontSize', 16);
% % %     
% % % end


% Further read: Matlab Help --> Linear Programming Algorithms

%% Optimization Part 2: Nonlinear Optimization - Unconstrained Optimization

% Tutorial for the Optimization Toolbox? **********************************


% fminsearch : Find minimum of unconstrained multivariable function using derivative-free method
% fminsearch uses the simplex search method of Lagarias et al. [1]. 
% This is a direct search method that does not use numerical or analytic gradients as in fminunc. 
% The algorithm is described in detail in fminsearch Algorithm. 
% The algorithm is not guaranteed to converge to a local minimum.
% 
% [1] Lagarias, J. C., J. A. Reeds, M. H. Wright, and P. E. Wright. 
% ?Convergence Properties of the Nelder-Mead Simplex Method in Low Dimensions.? 
% SIAM Journal of Optimization. Vol. 9, Number 1, 1998, pp. 112?147.
%
%
% min f(x) 
%  x


% fminunc    : Find minimum of unconstrained multivariable function


% objective function
% f(x) = 20*(-x1^2 + 2*x2)^2 + (x1 - 1)^2 + (x2 - 2)^2

% Gradient
% df/dx1 = 40*(-x1^2 + 2*x2)*(-2*x1) + 2*(x1 - 1)*(1) + 0
% df/dx2 = 40*(-x1^2 + 2*x2)*(2)     + 0              + 2*(x2 - 2)*(1)


% fminsearch -------------------------------------------------------------------------------------------
fun_fminsearch = [];

% Initial condition
% X0 = [0;0];

% Set options


% use fminsearch: using derivative-free method




% fminunc -------------------------------------------------------------------------------------------



% use fminunc: 




% Display Result
% % % disp(sprintf('x_opt_fminsearch (x1,x2) = (%.4f, %.4f)', x_opt_fminsearch(1), x_opt_fminsearch(2)) );
% % % disp(sprintf('x_opt_fminunc (x1,x2) = (%.4f, %.4f)', x_opt_fminunc(1), x_opt_fminunc(2)) );

% Further Read: Matlab Help --> Local vs. Global Optima
% Unconstrained Nonlinear Optimization Algorithms



%% Optimization Part 3: Nonlinear Optimization - Constrained Optimization
clc;
close all;

% fminbnd	Find minimum of single-variable function on fixed interval
% fmincon	Find minimum of constrained nonlinear multivariable function
% fseminf	Find minimum of semi-infinitely constrained multivariable nonlinear function


% fmincon **************

% min f(x) such that 
%  x
% 
% c(x)     <= 0
% c_eq(x)   = 0
% A_ineq*x <= b_ineq
% A_eq*x    = b_eq
% lb <= x <= ub

%[x,fval,exitflag,output,lambda,grad,hessian] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options)


% peaks;

% % % x_range = [-3:0.1:3]';
% % % y_range = [-3:0.1:3]';
% % % [x, y] = meshgrid(x_range, y_range);
% % % z = 3*(1-x).^2.*exp(-(x.^2) - (y+1).^2) ... 
% % %     - 10*(x/5 - x.^3 - y.^5).*exp(-x.^2-y.^2) ... 
% % %     - 1/3*exp(-(x+1).^2 - y.^2);

% figure;
% surf(x, y, z);
% xlabel('x'); ylabel('y'); zlabel('z');


% c(x) = []
% c_eq(x) = x^2 + y^2 <= 3^2



% % % 
% % % theta = [0:0.1:2*pi];
% % % figure;
% % % set(gcf, 'Position', [0 0 1280 1280]);
% % % for ii = 1:1
% % %     surfc(x, y, z);
% % %     xlabel('x'); ylabel('y'); zlabel('z');
% % %     hold on;
% % %     line(3*cos(theta), 3*sin(theta), -8*ones(size(theta)), 'LineWidth', 2, 'Color', 'm');
% % %     plot3( X0(1), X0(2), -8, 'Marker', 'x', 'MarkerSize', 10, 'Color', 'k', 'LineWidth', 2);
% % %     plot3( x_opt_fmincon(1), x_opt_fmincon(2), -8, 'Marker', 'x', 'MarkerSize', 10, 'Color', 'r', 'LineWidth', 2);
% % %     h_legend = legend('Surface', 'Contour', 'Nonlin Constraint', 'Start', 'Final');
% % %     xlabel('x');
% % %     ylabel('y');
% % %     zlabel('z');
% % %     set(gca, 'FontSize', 16);
% % %     grid on;
% % % end



% Try to avoid local minima 
% Randomly pick IC
x_opt_fmincon2     = [];
fval_fmincon2      = []; 
exitflag_fmincon2  = [];
output_fmincon2    = [];
x_opt_fmincon2_all = [];

for ii = 1:10
    
    
    
% % %     figure;
% % %     set(gcf, 'Position', [0 0 1280 1280]);
% % %     for qq = 1:1
% % %         
% % %         surfc(x, y, z);
% % %         xlabel('x'); ylabel('y'); zlabel('z');
% % %         hold on;
% % %         
% % %         line(3*cos(theta), 3*sin(theta), -8*ones(size(theta)), 'LineWidth', 2, 'Color', 'm');
% % %         
% % %         plot3( X0(1), X0(2), -8, 'Marker', 'x', 'MarkerSize', 16, 'Color', 'r', 'LineWidth', 2);
% % %         plot3( x_opt_fmincon2{ii}(1), x_opt_fmincon2{ii}(2), -8, 'Marker', 'x', 'MarkerSize', 16, 'Color', 'g', 'LineWidth', 2);
% % %         
% % %         plot3( X0(1), X0(2), z_peak_eval(X0), 'Marker', 'o', 'MarkerSize', 16, 'Color', 'r', 'LineWidth', 2, 'MarkerFaceColor', 'r');
% % %         plot3( x_opt_fmincon2{ii}(1), x_opt_fmincon2{ii}(2), z_peak_eval(x_opt_fmincon2{ii}), 'Marker', 'o', 'MarkerSize', 16, 'Color', 'g', 'LineWidth', 2, 'MarkerFaceColor', 'g');
% % %         
% % %         h_legend = legend('Surface', 'Contour', 'Nonlin Constraint', 'Start', 'Final');
% % %         xlabel('x');
% % %         ylabel('y');
% % %         zlabel('z');
% % %         title(sprintf('Iteration = %.0d, X opt = %.4f, %.4f', ii, x_opt_fmincon2{ii}(1), x_opt_fmincon2{ii}(1)));
% % %         set(gca, 'FontSize', 16);
% % %         grid on;
% % %     end

end



% Further Read
% https://www.mathworks.com/examples/matlab/community/35112-optimization-tips-tricks
% Matlab Global Optimization


%% Part 4: Optimization in Design

% This example demonstrates the integration of optimization + simmechanics
% for parametric design optimization purpose. 

Part2_6_Four_Bar_Main_Script;



%% Part 5: Intro to System Identification

% Step
% 1) Create model for this problem (In reality, we want to find this unknown model)
% 2) Design stimulus to perturb the system u: u (input) --> sys --> y (output)
% 3) Run experiment, given u, measure y !
% 4) Create frequency response plot (bode plot) from u(t), y(t) data
% 5) Estimate transfer function (tfest)/state space (ssest) 
% 6) Model Validaton (Compare estimated model Time/FRF to Experimental data)
% 7) 

clc;
close all;
clear all;


%% Step 1: Assume the unknown system is given by:

f1      = 5;
f2      = 200;
xi1     = 0.01;
xi2     = 0.1;
wn1     = 2*pi*f1;
wn2     = 2*pi*f2;
lambda1 = 10;
lambda2 = 0.2;

% % % s = tf('s');
% % % G =   lambda1*(wn1^2/(s^2 + 2*xi1*wn1*s + wn1^2)) + ...
% % %     1*lambda2*(wn2^2/(s^2 + 2*xi2*wn2*s + wn2^2));
% % % 
% % % figure,bode(G, {2*pi*0.1,2*pi*1000});



%% Step 2: Design Stimulus: Chirp, Sinusoidal Waveform, Square Waveform, Random Noise, ...

fs = 1000; % Hz
Ts = 1/fs; % Sampling Time [s];
t  = [0:Ts:100]'; % ******************** Time duration must long enough

f_start = 0.1;
f_end   = 400;



% % % figure;
% % % subplot(2,1,1);
% % % plot(t,u_chirp); xlabel('Time [s]'); ylabel('Chirp Stimulus'); grid on;
% % % subplot(2,1,2);
% % % spectrogram(u_chirp,256,250,256,1e3,'yaxis');

Amp_mod = ones(size(u_chirp));
delta_f = (f_end - f_start)/t(end);
f_chirp = f_start + delta_f*t;
for ii = 1:length(f_chirp)
    if (f_chirp(ii) >= 180) && (f_chirp(ii) <= 220)
        Amp_mod(ii) = 1;
    end
end

% % % lpf  = 2*pi*1/(s + 2*pi*1);
% % % lpfz = c2d(lpf, Ts);
% % % Amp_mod_lpf = filtfilt(lpfz.num{:}, lpfz.den{:}, Amp_mod);

figure; plot(Amp_mod); hold on; plot(Amp_mod_lpf);



% % % figure,plot(u_chirp2); 
% % % xlabel('Time [s]');
% % % ylabel('Modified Chirp');
% % % grid on;


%% Step 3: Response of the system given the designed stimulus

% lsim(sys,u,t,x0) ;




% % % figure;
% % % set(gcf, 'Position', [0 0 2560 1280]);
% % % ax1 = subplot(4,1,1); plot(t, u_chirp2); xlabel('Time [s]'); ylabel('Chirp'); grid on;
% % % ax2 = subplot(4,1,2); plot(t, f_chirp);  xlabel('Time [s]'); ylabel('F Chirp'); grid on;
% % % ax3 = subplot(4,1,3); plot(t, y_lsim);   xlabel('Time [s]'); ylabel('Response of the system'); grid on;
% % % ax4 = subplot(4,1,4); plot(t, y_noise);  xlabel('Time [s]'); ylabel('Response of the system + Noise'); grid on;
% % % linkaxes([ax1,ax2,ax3,ax4],'x')


%% Step 4: Compute Frequency Response (FFT)

% Impose linear system as assumption

% Fit in Time Domain or Fit in Frequency Response Domain

% Method 1: Fit in Frequency Response Domain
% Calculate FFT --> Basically find frequency response of the linear system 
% G(jw) = Output(jw)/Input(jw)



% % % figure,plot(f_vec, abs(Ujw));
% % % figure,plot(f_vec, abs(Yjw));






% Plot Frequency Response of Transfer Function
% % % figure(101);
% % % set(gcf, 'Position', [0 0 2560 1280]);
% % % for ii = 1:1
% % %     
% % %     
% % %     ax1 = subplot(2,1,1);
% % %     semilogx(f_vec, Mag_Hjw, 'LineWidth', 2, 'Color', 'b');
% % %     hold on;
% % %     %semilogx(f_vec, filtfilt(LPFz.num{:}, LPFz.den{:}, Mag_Hjw), 'LineWidth', 2, 'Color', 'r');
% % %     xlabel('Frequency [Hz]');
% % %     ylabel('Magnitude H(jw) in db');
% % %     h_legend = legend('|H(jw)|');
% % %     set(h_legend, 'Location', 'NorthEast', 'Color', [1 1 0.9]);
% % %     set(gca, 'FontSize', 16);
% % %     grid on;
% % %     xlim([0.1 max(f_vec)]);
% % %     
% % %     
% % %     ax2 = subplot(2,1,2);
% % %     semilogx(f_vec, Phs_Hjw, 'LineWidth', 2, 'Color', 'b');
% % %     hold on;
% % %     xlabel('Frequency [Hz]');
% % %     ylabel('Phase H(jw) in deg');
% % %     h_legend = legend('\angle H(jw)');
% % %     set(h_legend, 'Location', 'NorthEast', 'Color', [1 1 0.9]);
% % %     set(gca, 'FontSize', 16);
% % %     grid on;
% % %     xlim([0.1 max(f_vec)]);
% % %     linkaxes([ax1,ax2],'x');
% % %     
% % % end





%% Step 5: Estimate Transfer Function

% tfest, ssest
% sys = ssest(data,nx,Name,Value)

%
% Frequency Response Data Object
% h = idfrd(Response,Freq,Ts)





% Step 6: Model Validation

% Estimated Model at 4th order


% True Model G


% % % figure(102);
% % % set(gcf, 'Position', [0 0 2560 1280]);
% % % for ii = 1:1
% % %     
% % %     ax1 = subplot(2,1,1);
% % %     semilogx(f_vec, mag_ss_G, 'LineWidth', 2, 'Color', 'g');
% % %     hold on;
% % %     semilogx(f_vec, Mag_Hjw, 'LineWidth', 2, 'Color', 'b');
% % %     semilogx(f_vec, mag_ss_est4, 'LineWidth', 2, 'Color', 'r', 'LineStyle', '--');
% % %     xlabel('Frequency [Hz]');
% % %     ylabel('Magnitude H(jw) in db');
% % %     h_legend = legend('True Model', '|H(jw)| Experimental Data', 'Fit: Estimated Model (SS)');
% % %     set(h_legend, 'Location', 'NorthEast', 'Color', [1 1 0.9]);
% % %     set(gca, 'FontSize', 16);
% % %     grid on;
% % %     xlim([0.1 max(f_vec)]);
% % %     
% % %     
% % %     ax2 = subplot(2,1,2);
% % %     semilogx(f_vec, phs_ss_G, 'LineWidth', 2, 'Color', 'g');
% % %     hold on;
% % %     semilogx(f_vec, Phs_Hjw, 'LineWidth', 2, 'Color', 'b');
% % %     semilogx(f_vec, phs_ss_est4, 'LineWidth', 2, 'Color', 'r', 'LineStyle', '--');
% % %     xlabel('Frequency [Hz]');
% % %     ylabel('Phase H(jw) in deg');
% % %     set(gca, 'FontSize', 16);
% % %     grid on;
% % %     xlim([0.1 max(f_vec)]);
% % %     linkaxes([ax1,ax2],'x');
% % %     
% % % end



% Report
% % % ss_4.Report
% % % ss_4.Report.Fit


%% Time Domain Based System Identification

% % % % Input
% % % u_chirp2;
% % % 
% % % % Output
% % % y_noise2 = y_lsim + 1*0.40*randn(size(y_lsim)); %*****************************************
% % % 
% % % 
% % % figure,plot(t,u_chirp2, t, y_noise2); grid on;
% % % 
% % % iddata_7 = iddata(y_noise2, u_chirp2, Ts);
% % % 
% % % np = 4;
% % % nz = 4;
% % % opt_tfest = tfestOptions('EnforceStability', true, ...
% % %                          'Display', 'on');
% % %                      
% % % sys_tfest = tfest(iddata_7, np, nz, opt_tfest);
% % % 
% % % 
% % % % Compare data vs estimated linear model in time/frequency domain (in this case)
% % % compare(iddata_7, sys_tfest);




% Further Read
% Frequency Domain Identification: Estimating Models Using Frequency Domain Data
% tfestimate 