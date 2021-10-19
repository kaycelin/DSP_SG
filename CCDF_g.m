%% A1, 2021-10-19, add fnum_legend

function [PAR] = CCDF_g(x, Nsamps, fnum, fnum_legend)
% CCDF: Complementary Cumulative Distribution function

if ~exist('fnum_legend','var')||isempty(fnum_legend)
    flag_fnum_legendend = 0;
else
    flag_fnum_legendend = 1;
end
if ~exist('Nsamps','var')||isempty(Nsamps)
    Nsamps = 1000;
end
Nbranch = size(x,2);

% % column transformer
% x=x(:);

% Power calculation
P_x = x.*conj(x);
if Nbranch>1
    P_x_rms = mean((sqrt(mean(P_x))).^2); % mean for Nbranch
else
    P_x_rms = sqrt(mean(P_x)).^2; % mean for Nbranch
end
%% xaxis, P_vector
% generate pwr vector for xaxis
P_vector = 0:(max(P_x)/Nsamps):max(P_x);

%% yaxis
% histc the P_x(k) in range of P_vcetor(k)
P_x_histc = histc(P_x, P_vector);
% Probility Distribution Function (PDF)
PDF_P_x = P_x_histc/sum(P_x_histc);
% Cumulative Distrbution Function (CDF)
CDF_P_x = cumsum(PDF_P_x);
% Complemetary Cmulative Distrbution Function (CCDF)
CCDF_P_x = 1-CDF_P_x;

%% x and y lim, PAPR_vector
PAPR_vector = 10*log10(P_vector/P_x_rms);
PAPR_vector = PAPR_vector(2:end); % PAPR_vector min ~=0
ind_xaxis_min = find(PAPR_vector>=0);
xaxis_min = floor(PAPR_vector(ind_xaxis_min(1)));
xaxis_max = ceil(max(PAPR_vector));

CCDF_P_x = CCDF_P_x(1:end-1); % CCDF max ~=1
yaxis_min = min(CCDF_P_x);
yaxis_max = 1;

PAR = max(PAPR_vector);

%% plot
if exist('flag_debug_plot','var') || ~isempty(fnum)
    if flag_fnum_legendend
        dispLegend = [fnum_legend,', PAR: ' num2str(fix(PAR*10)/10), 'dB'];
    else
        dispLegend = ['PAR: ' num2str(fix(PAR*10)/10), 'dB'];
    end
    
    figure(fnum)
    semilogy(PAPR_vector,CCDF_P_x, 'DisplayName' ,dispLegend),hold on, grid on, legend;
    xlabel('dB above average power'),ylabel('probability')
    title('CCDF'),axis([xaxis_min xaxis_max yaxis_min yaxis_max])
    % text(2.8,.25,'sine'),text(7.2,.025,'AWGN')
end

end