% Detrended Fluctuation Analysis on time series (EEG)
% Copyright 2019-22, Maitreyee Wairagkar 
%
% ** References ** 
% Wairagkar et al., "Dynamics of Long-Range Temporal Correlations in Broadband EEG During Different Motor Execution and Imagery Tasks",  
% Frontiers in Systems Neuroscience, 2021  https://doi.org/10.3389/fnins.2021.660032
%
% Wairagkar et al., "Modeling the Ongoing Dynamics of Short and Long-Range Temporal Correlations in Broadband EEG During Movement",  
% Frontiers in Neuroscience, 2019   https://doi.org/10.3389/fnsys.2019.00066  
%
% ** Input **
% sig        : single time series 
% min_scale  : minimum windows or timescales for DFA, should atleast be 10
% max_scale  : maximum windows or timescales for DFA, should atleast be len(sig)/4
% log_base   : 2 or 10
% N1         : number of timescales n over which you want to compute RMSE fluctuations, should be atleast 4
%
% ** Output **
% Alpha1     : scaling exponent of DFA - slope of log-log fluctuation plot
% n          : time windows
% F_n        : the output RMSE fluctuation at each time window n 
% c          : intercept of the log-log fluctuation plot
% r2         : rsquared of fitted line to log-log fluctuations
% fitted_line: line fitted to the log-log fluctuation plot


function [Alpha1,n,F_n,c,r2,fitted_line]=DFA_MW(sig,min_scale,max_scale, log_base, N1)

% set range i.e. timescales for computing DFA
if log_base == 2
    n_=linspace(log2(min_scale),log2(max_scale),N1);  % log2 
else 
    n_=linspace(log10(min_scale),log10(max_scale),N1);% log10 
end
n= round(log_base.^n_); 

F_n=zeros(N1,1);

% step1 - subtract mean
mean_sig = sum(sig)/length(sig);
sig = sig - mean_sig;

% step2 - Integrate the signal
y = cumsum(sig); 

% step3 - Find F(n) - Forward direction
for i=1:N1   % for each window, find F(n)
 f_n_forward(i) = DFA(y,n(i),1);
end

% Find F(n) - Backward direction
y_ = cumsum(fliplr(sig)); 
for i=1:N1   
 f_n_backward(i) = DFA(y_,n(i),1);
end

% Final F(n) is mean in forward and backward direction
F_n = mean([f_n_forward;f_n_backward]);

% step4 - find scaling exponent alpha by fitting line to log-log plot
if log_base == 2
    log_F_n = log2(F_n);
    log_n = log2(n);
else 
    log_F_n = log10(F_n);
    log_n = log10(n);
end

[A,~,mu]=polyfit(log_n,log_F_n,1); % optimised version - polyfit_MW

Alpha1=A(1);  % slope
c = A(2);     % intercept
        
fitted_line = polyval(A,log_n,[],mu); % optimised version - polyval_MW
r2 = rsquared(log_F_n, fitted_line);

%FD = 3-Alpha1; %Fractle Dimension - not used
end
 
 function F_n = DFA(y,win_len,order)
       
   n=floor(length(y)/win_len); % find number of blocks of given win_len in the signal
   N=n*win_len; % length of signal such it it can be divided into n windows
   y=y(1:N);    
   
   % for each block fit the line
    for i=1:n 
        [P,~,MU] = polyfit_MNW(1:win_len,y(((i-1)*win_len+1):i*win_len),order);
        Yn(((i-1)*win_len+1):i*win_len)= polyval_MNW(P,1:win_len,[],MU);
    end

    F_n = sqrt ( sum((y-Yn).^2)/N ); % RMSE
   
 end
 
 function r2 = rsquared(y,f)
    r2 = 1 - (sum((y(:)-f(:)).^2)/sum((y(:)-mean(y(:))).^2));
 end
 