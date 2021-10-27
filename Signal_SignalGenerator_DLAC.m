%% PA1, 2021-08-26, E/// radon interpolation function: [IPOL2-0, SRC, IPOL2-1, IPOL2-2, IPOL2-3, IPOL2-4, IPOL2-5, IPOL2-6]
%% PA2, 2021-10-04, Signal Generator based on Radon_interpolation_SRC_g
%% A1, 2021-10-25, Add AC signal, flag_sigType: DL/AC/DLAC
%% A2, 2021-10-26, Update for multibranches, Nbr
%% A3, 2021-10-27, Equal subcarrier for DLAC

clear all
close all
%% example: /dl/nrTdd200_23040/iFilter char "sfir_nrTdd200_23040 ipol2_0_nrTdd200_23040 src_nrTdd200_23040 ipol2_1_nrTdd200_23040 ipol2_2_nrTdd200_23040 bypass bypass bypass bypass "
flag_sigIn = 'gen'
flag_sigType = 'DLAC' % DL/AC/DLAC
fnum = 1004

%% input: signal
switch flag_sigIn
    case 'gen'
        if strcmpi(flag_sigType, 'DL')||strcmpi(flag_sigType, 'DLAC')
            flag_fs_extend = 1;
            if flag_fs_extend
                sigConfig.bw_Channel = '20MHz'
                sigConfig.fs = 122.88e6;
                [sigDL, ~, Config] = OFDM_Mod_g(sigConfig,'64QAM','off','NR',fnum,[]);
                save('waveform_OFDM_20MHz_122p88MHz.mat','sigDL')
                sigDL_cell = [{sigDL}, {Config}];
                save('waveform_OFDM_20MHz_122p88MHz_cell.mat','sigDL_cell')
            else
                [sigInDL, ~, Config] = OFDM_Mod_g('20MHz','64QAM','off','NR',fnum,[]);
                save('waveform_OFDM_20MHz_23p04MHz.mat','sigInDL')
                sigAC_cell = [{sigInDL}, {Config}];
                save('waveform_OFDM_20MHz_23p04MHz_cell.mat','sigAC_cell')
            end
            sigIn = sigInDL;
        end
        if strcmpi(flag_sigType, 'AC')||strcmpi(flag_sigType, 'DLAC') %% A1, 2021-10-25, Add AC signal, flag_sigType: DL/AC/DLAC
            flag_fs_extend = 1;
            if flag_fs_extend
                fsAC = [23.04e6, 122.88e6]
            else
                fsAC = [23.04e6]*ones(1,2)
            end
            df = 100
            Nsamps = fsAC(2)/df
            Nbr = 4
            NsampsAC = [Nsamps, Nbr]
            dACSubcarrier = 240e3;
            bwChannelAC = {'20MHz',fsAC,dACSubcarrier}
            [sigAC, ACconfig] = AntCal_genACSource_g100(bwChannelAC, [], [NsampsAC], [], [], fnum ,[]);
            %             save(['waveform_AC_20MHz_23p04MHz_4Brs.mat'],'sigAC')
            save(['waveform_AC_20MHz_122p88MHz_4Brs.mat'],'sigAC')
            sigAC_cell = [{sigAC}, {ACconfig}];
            %             save('waveform_AC_20MHz_23p04MHz_4Brs_cell.mat','sigAC_cell')
            save('waveform_AC_20MHz_122p88MHz_4Brs_cell.mat','sigAC_cell')
            sigIn = sigAC;
            Config = ACconfig;
        end
    case 'lod'
        load('waveform_OFDM_20MHz_23p04MHz.mat')
        fs = 23.04e6
        bwInband = 9e6*[-1 1]
        bwChannel = 20e6
end

%% output: carrier configuration
fs = Config.fs
bwInband = (Config.bwCarrier/2+0.5e6)*[-1 1]
bwChannel = Config.bwChannel
Nsamps = size(sigIn,1)
Nbr = size(sigIn,2)
df = fs/Nsamps
PLOT_FFT_dB_g(sigIn, fs, length(sigIn), ['sigIn'], 'df', 'full', 'pwr', [fnum+1]);

%% input: channel filter
FIRc1_Wtype = "kaiser"
FIRc1_Ftype = "LPF"
FIRc1_Order = NaN
FIRc1_fTolerance = -0.1e6
FIRc1_K_AttdB = 60
FIRc1_K_fdelta = 0.5e6
FIRc1_fcutoffL = bwChannel/2
FIRc1_fcutoffH = 0
FIRc1_Export = fs
NCarriers = 1
b_ch = SYM_FIRApp(FIRc1_Wtype,FIRc1_Ftype,FIRc1_Order,FIRc1_K_AttdB,FIRc1_K_fdelta,FIRc1_fTolerance,FIRc1_fcutoffL,FIRc1_fcutoffH,df,bwInband,NCarriers,FIRc1_Export,fnum+2);
b_ch = b_ch{:};

flag_ch_DB = 0
if flag_ch_DB
    sfir_nrTdd200_23040 = [-22,-1,21,-41,55,-59,47,-18,-23,71,-113,138,-133,94,-22,-71,167,-242,274,-246,153,-7,-167,333,-449,479,-402,217,47,-341,601,-761,769,-600,269,170,-630,1009,-1207,1157,-835,283,405,-1087,1607,-1827,1660,-1097,218,815,-1786,2470,-2677,2304,-1370,17,1496,-2850,3721,-3856,3140,-1633,-423,2631,-4517,5619,-5592,4300,-1866,-1316,4639,-7381,8848,-8520,6191,-2050,-3299,8893,-13545,16035,-15319,10738,-2166,-9916,24446,-39902,54509,-66499,74370,447210];
    b0 = [sfir_nrTdd200_23040, sfir_nrTdd200_23040(end-1:-1:1)];
    PLOT_FFT_dB_g(b0/2^19*Nsamps, fs, Nsamps, ['ch , Length:',num2str(numel(b0))], 'df', 'full', 'pwr', [fnum+2]);
end

%% output: channel filter output
if strcmpi(flag_sigType, 'DLAC')
    for ibr = 1:size(sigDL,2)
        sigIn_ch_DL(:,ibr) = conv(sigDL(:,ibr), b_ch, 'same');
    end
    for ibr = 1:size(sigAC,2)
        sigIn_ch_AC(:,ibr) = conv(sigAC(:,ibr), b_ch, 'same');
    end
    PLOT_FFT_dB_g(sigIn_ch_DL(:,:), fs, Nsamps*1, ['sigDL + channel filter'], 'df', 'full', 'pwr', [fnum+2]);
    PLOT_FFT_dB_g(sigIn_ch_AC(:,:), fs, Nsamps*1, ['sigAC + channel filter'], 'df', 'full', 'pwr', [fnum+2]);
else
    for ibr = 1:Nbr
        sigIn_ch(:,ibr) = conv(sigIn(:,ibr), b_ch, 'same');
    end
    PLOT_FFT_dB_g(sigIn_ch(:,:), fs, Nsamps*1, ['sigIn + channel filter'], 'df', 'full', 'pwr', [fnum+2]);
    for ibr = 1:Nbr
        [evm(ibr), delay(ibr)] = dsp_evm_timexcorr_inband_g(sigIn(:,ibr), circshift(sigIn_ch(:,ibr),57), fs, bwInband, 1, 2)
    end
end

%% input: digital gain/pOut assignment
if strcmpi(flag_sigType, 'DLAC')
    PdB_expected_DL = -15
    PdB_expected_AC = -15-20*log10(240/30) %% 2021-10-27, Equal subcarrier for DLAC
    v_expected_DL = sqrt(10.^((PdB_expected_DL)/10));
    v_expected_AC = sqrt(10.^((PdB_expected_AC)/10));
    v_signal_DL = mean( sqrt(mean(abs(sigIn_ch_DL).^2)) );
    v_signal_AC = mean( sqrt(mean(abs(sigIn_ch_AC).^2)) );
    sigIn_ch_DL_gain = sigIn_ch_DL.*v_expected_DL/v_signal_DL;
    sigIn_ch_AC_gain = sigIn_ch_AC.*v_expected_AC/v_signal_AC;
    PLOT_FFT_dB_g(sigIn_ch_DL_gain, fs, Nsamps, ['sigIn ch + gain, pwr output dB:',num2str(PdB_expected)], 'df', 'full', 'pwr', [fnum+2]);
    PLOT_FFT_dB_g(sigIn_ch_AC_gain, fs, Nsamps, ['sigIn ch + gain, pwr output dB:',num2str(PdB_expected)], 'df', 'full', 'pwr', [fnum+2+1]);
else
    PdB_expected = -15
    v_expected = sqrt(10.^((PdB_expected)/10));
    v_signal = mean( sqrt(mean(abs(sigIn_ch).^2)) );
    sigIn_ch_gain = sigIn_ch.*v_expected/v_signal;
    PLOT_FFT_dB_g(sigIn_ch_gain, fs, Nsamps, ['sigIn ch + gain, pwr output dB:',num2str(PdB_expected)], 'df', 'full', 'pwr', [fnum+2]);
end

%% input: interpolation filter coefficient and fsOut
if ~strcmpi(flag_sigType, 'DLAC')
    flag_ipo_filter = 'generate'
    switch flag_ipo_filter
        case 'database'
            sfir_nrTdd200_23040 = [-22,-1,21,-41,55,-59,47,-18,-23,71,-113,138,-133,94,-22,-71,167,-242,274,-246,153,-7,-167,333,-449,479,-402,217,47,-341,601,-761,769,-600,269,170,-630,1009,-1207,1157,-835,283,405,-1087,1607,-1827,1660,-1097,218,815,-1786,2470,-2677,2304,-1370,17,1496,-2850,3721,-3856,3140,-1633,-423,2631,-4517,5619,-5592,4300,-1866,-1316,4639,-7381,8848,-8520,6191,-2050,-3299,8893,-13545,16035,-15319,10738,-2166,-9916,24446,-39902,54509,-66499,74370,447210];
            ipol2_0_nrTdd200_23040 = [4,-10,22,-44,80,-134,214,-330,492,-714,1006,-1390,1886,-2516,3312,-4316,5580,-7188,9262,-12022,15878,-21710,31824,-54684,166568];
            src_nrTdd200_23040 = [-9,92,-452,1485,-3834,9974,61806,-4363,998,-173,11,2,-19,183,-912,3136,-8918,30194,49388,-9958,3149,-854,162,-16,-16,162,-854,3149,-9958,49388,30194,-8918,3136,-912,183,-19,2,11,-173,998,-4363,61806,9974,-3834,1485,-452,92,-9];
            ipol2_1_nrTdd200_23040 = [-1018,7912,-34060,158236];
            ipol2_2_nrTdd200_23040 = [-358,3324,-16012,78582];
            b0 = sfir_nrTdd200_23040;
            b1 = ipol2_0_nrTdd200_23040;
            b2 = src_nrTdd200_23040;
            b3 = ipol2_1_nrTdd200_23040;
            b4 = ipol2_2_nrTdd200_23040;
            bCell = {b1 ,b2, b3, b4};
            
            %         fs_ipo2_0 = fs*2
            %         fs_ipo4 = fs_ipo2_0*4/1
            %         bPP_src_tmp = reshape(src_nrTdd200_23040, [], 4).';
            %         bPP_src = bPP_src_tmp(:);
            %         Nbits_src = 18;
            %         PLOT_FFT_dB_g(bPP_src/2^Nbits_src*Nsamps, fs_ipo4, Nsamps, ['ipo4 filter for SRC, database'], 'df', 'full', 'pwr', [fnum+3, 2, 1, 1]);
            %         subplot(2,1,2)
            %         plot(bPP_src/2^Nbits_src, 'DisplayName', 'time domain of polyphase filter, database'), hold on, legend
            
        case 'generate'
            %% input/output: generate HB FIR and SRC FIR
            % ipol2-0
            ipo_ratio = [2, 4/3, 2]; % input
            fs_ipo2_0 = fs*ipo_ratio(1);
            %         bHB_ipo2_0 = FIR_HB_WIN([], 'Remez', [0.001, 100], 1*[0.45,0.55]*(fs_ipo2_0/2), [], fs_ipo2_0, 0, [fnum+3, numel(ipo_ratio), 1, 1], 'ipol2-0');
            bHB_ipo2_0= FIR_HB_WIN([], 'Remez', [0.001, 100], 1*[0.888,1]*(fs_ipo2_0/2), [], fs_ipo2_0, 'HB', [fnum+3, numel(ipo_ratio), 1, 1], 'ipol2-0');
            
            b1_tmp = bHB_ipo2_0*2^0;
            b1_tmp(2:2:end) = [];
            b1_tmp = b1_tmp(1:numel(b1_tmp)/2);
            
            % src
            if ipo_ratio(2)==4/3
                Nups_src = 4;
                Ndns_src = 3;
            end
            fs_ipo4 = fs_ipo2_0*Nups_src;
            Ntaps_src = Nups_src*12;
            Ntaps_src = 48 % input
            if mod(Ntaps_src, Nups_src)~=0
                error('Ntaps of src should be Nups_src*N !')
            end
            bPP_src_cell = SYM_FIRApp("HAN","LPF",[Ntaps_src],[],[],5e6,20e6,0,df,[],NCarriers,fs_ipo4,[fnum+3, numel(ipo_ratio), 1, 2], 'ipol4-src');
            bPP_src = cell2mat(bPP_src_cell);
            b2_tmp = reshape(bPP_src, Nups_src, []).';
            b2_tmp = b2_tmp(:).';
            
            % ipol2-1
            fs_src = fs_ipo2_0*ipo_ratio(2);
            fs_ipo2_1 = fs_src*ipo_ratio(3);
            f_toler = 10e6 % input
            Ntaps_ipol21 = 45 % input
            bHB_ipo2_1= FIR_HB_WIN([15], 'Remez', [0.001, 20], 1*[0.2,1]*(fs_ipo2_1/2), [], fs_ipo2_1, 'HB', [fnum+3, numel(ipo_ratio), 1, 3], 'ipol2-1');
            
            b3_tmp = bHB_ipo2_1*2^0;
            b3_tmp(2:2:end) = [];
            b3_tmp = b3_tmp(1:numel(b3_tmp)/2);
            
            Nbits = [19 18 19 18 18 18 18 18];
            b1 = b1_tmp*2^Nbits(1);
            b2 = b2_tmp*2^Nbits(2);
            b3 = b3_tmp*2^Nbits(3);
            b4 = [];
            bCell = {b1 ,b2, b3};
    end
    
    flag_debug_plt_fir = 0
    if flag_debug_plt_fir
        ipo_ratio = [2, 4/3, 2];
        Nbits = [19 18 19 18 18 18 18 18];
        fs_ipo2_0 = fs*ipo_ratio(1);
        bHB_ipo2_0 = DB_Taps_HB_Symmtric(b1, Nbits(1), fs_ipo2_0, [fnum+4, numel(ipo_ratio), 1, 1], ['ipol2-0, Length:',num2str(numel(b1)*ipo_ratio(1)*2-1)] );
        
        fs_ipo4 = fs_ipo2_0*4/1;
        bPP_src_tmp = reshape(b2, [], 4).';
        bPP_src = bPP_src_tmp(:);
        PLOT_FFT_dB_g(bPP_src/2^Nbits(2)*Nsamps, fs_ipo4, Nsamps, ['ipol4-src, Length:',num2str(numel(bPP_src))], 'df', 'full', 'pwr', [fnum+4, numel(ipo_ratio), 1, 2]);
        fs_src = fs_ipo2_0*ipo_ratio(2);
        
        fs_ipo2_1 = fs_src*ipo_ratio(3);
        bHB_ipo2_1=DB_Taps_HB_Symmtric(b3, Nbits(3), fs_ipo2_1, [fnum+4, numel(ipo_ratio), 1, 3], ['ipol2-1, Length:',num2str(numel(b3)*ipo_ratio(3)*2-1)] );
        
    end
    
    %% input: fsOut
    fsOut = 122.88e6
    
    %% output: signal interpolation with Fixed points
    for ibr = 1:Nbr
        [sigOut(:,ibr), DelayFilt, LengthFilt, DelayFiltSum, sigOutCell(ibr,:)] = radon_interpolation_g1(sigIn_ch_gain(:,ibr), fs, fsOut, bCell, 0, fnum+5);
    end
    % [sigOut, DelayFilt, LengthFilt, DelayFiltSum, sigOutCell] = radon_interpolation_g1(sigIn_ch_gain, fs, fsOut, bCell, 1, fnum+5);
    PLOT_FFT_dB_g(sigOut, fsOut, numel(sigOut), ['sigOut'], 'df', 'full', 'pwr', [fnum+6]);
    
    %% output: remove delay
    flag_removeDelay = 1
    if flag_removeDelay
        for ibr = 1:Nbr
            sigOut_o = sigOut(:,ibr);
            delay = DelayFiltSum(1)+1 % Add for fractional delay
            length_sigOut = numel(sigIn(:,ibr))*ipo_ratio(1)*ipo_ratio(2)*ipo_ratio(3)
            sigOut_rmd(:,ibr) = sigOut_o([1:length_sigOut]+delay);
            PLOT_FFT_dB_g(sigOut_rmd(:,ibr), fsOut, numel(sigOut_rmd(:,ibr)), ['sigOut remove delay'], 'df', 'full', 'pwr', [fnum+6]);
        end
    end
    
    %% check: EVM
    % Nsamps
    Nsamps_sigIn = length(sigIn)
    Nsamps_sigIn == (floor( ((length(sigOut)-DelayFilt(3,3))/2-DelayFilt(2,3))*3/4 ) - DelayFilt(1,3))/2
    
    % Synchronization
    if ~flag_removeDelay
        sigOut_DNS1 = sigOut(1:2:end,:);
    else
        sigOut_DNS1 = sigOut_rmd(1:2:end,:);
    end
    
    sigOut_DNS1_UPS = interp1( [1:3:3*length(sigOut_DNS1)], sigOut_DNS1, [1:1:3*length(sigOut_DNS1)], 'linear');
    sigOut_DNS1_UPS = interp1( [1:3:3*length(sigOut_DNS1)], sigOut_DNS1, [1:1:3*length(sigOut_DNS1)], 'spline');
    sigOut_DNS1_UPS = interp1( [1:3:3*length(sigOut_DNS1)], sigOut_DNS1, [1:1:3*length(sigOut_DNS1)], 'v5cubic');
    sigOut_DNS1_UPS = interp1( [1:3:3*length(sigOut_DNS1)], sigOut_DNS1, [1:1:3*length(sigOut_DNS1)], 'makima');
    
    sigOut_DNS1_UPS_DNS2 = sigOut_DNS1_UPS(1:4:end,:);
    sigOut_DNS1_UPS_DNS2_DNS3 = sigOut_DNS1_UPS_DNS2(1:2:end,:);
    for ibr = 1:Nbr
        [sigOut_DNS1_UPS_DNS2_DNS3_shift, DelayInt, DelayFrac] = dsp_Delay_EstCor_g2(sigIn(:,ibr), sigOut_DNS1_UPS_DNS2_DNS3(:,ibr), fnum+7);
        sigOut_DNS1_UPS_DNS2_DNS3_shift_trunk = sigOut_DNS1_UPS_DNS2_DNS3_shift(1:Nsamps_sigIn);
        [evm(ibr)] = dsp_evm_g(sigIn(:,ibr), sigOut_DNS1_UPS_DNS2_DNS3_shift_trunk)
    end
    
    %% input: NCO / output: carrier sum
    fNCO = [-20e6 20e6]
    flag_SUMCarrier = 1
    sigIn_1C = sigOut_rmd;
    disp_title = '2C'
    for ibr = 1:Nbr
        [~, sigOut_2C(:,ibr), ~, ~, ~] = SYM_NCOApp(sigIn_1C(:,ibr), fNCO, fsOut, bwInband, flag_SUMCarrier, fnum+8, disp_title, 0);
    end
else
    fsOut = fs;
    fNCO = [-20e6 0 20e6]
    flag_SUMCarrier = 1
    sigIn_DLAC_cell = [{sigIn_ch_DL_gain}; {sigIn_ch_AC_gain}; {sigIn_ch_DL_gain}]; % row
    disp_title = 'DL+AC+DL'
    [~, sigOut_DLAC_cell, ~, ~, ~] = SYM_NCOApp(sigIn_DLAC_cell, fNCO, fsOut, bwInband, flag_SUMCarrier, fnum+8, disp_title, 0);
    sigOut_DLAC = cell2mat(sigOut_DLAC_cell);
end

%% Export
if strcmpi(flag_sigType, 'DLAC')
    PLOT_FFT_dB_g(sigOut_DLAC, fsOut, length(sigOut_DLAC), ['sigOut DLAC:DL+AC+DL'], 'df', 'full', 'pwr', [fnum+6]);
    signal = sigOut_DLAC;
    save('waveform_DLAC_3x20MHz_122p88MHz.mat','signal')
else
    PLOT_FFT_dB_g(sigOut_rmd, fsOut, length(sigOut_rmd), ['sigOut remove delay'], 'df', 'full', 'pwr', [fnum+6]);
    signal = sigOut_rmd;
    save('waveform_OFDM_20MHz_122p88MHz.mat','signal')
    
    PLOT_FFT_dB_g(sigOut_2C, fsOut, length(sigOut_2C), ['sigOut 2C'], 'df', 'full', 'pwr', [fnum+9]);
    signal = sigOut_2C;
    save('waveform_OFDM_20MHz_2C_122p88MHz.mat','signal')
end

%% PAR comparsion
[PAR_DLAC] = CCDF_g(sigOut_DLAC, length(sigOut_DLAC), fnum+10, 'DLAC')
[PAR_1C] = CCDF_g(sigOut_rmd, length(sigOut_rmd), fnum+10, '1C')
[PAR_2C] = CCDF_g(sigOut_2C, numel(sigOut_2C), fnum+10, '2C')

%% function
function IQ = bin2IQ(fid, typeFormat, bitsFormat, IQFormat, fs)
if ~exist('typeFormat','var') || isempty(typeFormat)
    typeFormat = 'int32'
end
if ~exist('bitsFormat','var') || isempty(bitsFormat)
    bitsFormat = [1 14 3]
end
if ~exist('IQFormat','var') || isempty(IQFormat)
    IQFormat = 1
end

fip = fopen(fid,'rb');
[Data, ~] = fread(fip,inf,typeFormat);
fclose(fip)
Data = Data/2^(bitsFormat(3));
if IQFormat
    I=Data(1:2:end);
    Q=Data(2:2:end);
else
    I=Data(1:length(Data)/2);
    Q=Data(length(Data)/2:end);
end
IQ = I+1i*Q;
if exist('fs','var') && ~isempty(fs)
    fnum = 101
    Nfft = 2048
    PLOT_FFT_dB_g(IQ, fs, Nfft, ['IQ'], 'df', 'full', 'pwr', [fnum]);
end
end

function [meaCor, DelayInt, DelayFrac] = dsp_Delay_EstCor_g2(ref, mea, fnum)
if exist('fnum','var') && ~isempty(fnum)
    flag_fnum = 1;
else
    flag_fnum = 0;
end

if isrow(ref)||isrow(mea)
    ref = ref(:);
    mea = mea(:);
end
if length(mea)>length(ref)
    ref = [ref; zeros(abs(length(mea)-length(ref)),1)];
else
    mea = [mea; zeros(abs(length(ref)-length(mea)),1)];
end
method_FFTdiff = 0;
if method_FFTdiff
    FFT_diff = fft(mea).*conj(fft(ref)); %!!!
else
    FFT_diff = conj(fft(mea)).*fft(ref); %!!!
end
Xcor = ifft(FFT_diff); % Cross-correlation
ix = find(abs(Xcor) == max(abs(Xcor)));
ix = ix(1);
DelayInt = ix - 1;

if flag_fnum
    figure(fnum), hold on;
    plot(abs(Xcor), 'b');
    h = plot(DelayInt+1, abs(Xcor(DelayInt+1)), 'r+')
    set(h, 'lineWidth', 3);
    legend('Xcor', ['DelayInteger: ',num2str(DelayInt)]);
    title('Integer Delay Estimation');
    
    figure(fnum*10)
    plot(real(ref),'DisplayName',['ref']), legend, hold on
    if method_FFTdiff
        plot(real(circshift(mea,-DelayInt)),'DisplayName',['mea shift:-',num2str(DelayInt)]), legend, hold on
    else
        plot(real(circshift(mea,DelayInt)),'DisplayName',['mea shift:',num2str(DelayInt)]), legend, hold on
    end
end

Nsamps = numel(mea);
nbins = (mod([[0:Nsamps-1]+floor(Nsamps/2)], Nsamps) - floor(Nsamps/2)) / Nsamps; % [0:Nsamps/2-1]=>phs pos. and [Nsamps/2:end]=>phs neg.
nbins = nbins(:);
Ang_UnitDelay = 2*pi*nbins;
if nargout>=2
    FFT_DelayCor = FFT_diff.*exp(1i*Ang_UnitDelay*DelayInt);
    weight = abs(FFT_DelayCor);
    AngWeight_UnitDelay = Ang_UnitDelay.* weight;
    AngWeight_DiffFFTCor = -angle(FFT_DelayCor).* weight; % Negative Phase corresponse to Time Delay, and Positive Phase correponse to Time Advance
    
    DelayFrac = AngWeight_UnitDelay \ AngWeight_DiffFFTCor
else
    DelayFrac = 0
end
Delay = DelayInt + DelayFrac; % Delay always > 0
if method_FFTdiff
    MEACor = fft(mea).* exp(1i*Ang_UnitDelay.*Delay); % Phase<0:Time Delay, and Phase>0:Time Advance
else
    MEACor = fft(mea).* exp(1i*Ang_UnitDelay.*-Delay); % Phase<0:Time Delay, and Phase>0:Time Advance
end
meaCor = ifft(MEACor);

if flag_fnum
    figure(fnum*11)
    plot(real(ref),'DisplayName',['ref']), legend, hold on
    plot(real(meaCor),'DisplayName',['meaCor']), legend, hold on
end
end

function [sigOutFilt, DelayFilt, lengthFilt, DelayFiltSum, sigOutFiltCell] = radon_interpolation_g1(sigIn, fsIn, fsOut, bCell, flag_removeDelay, fnum)
if ~exist('flag_removeDelay','var')||isempty(flag_removeDelay) % 2021-08-24
    flag_removeDelay = 0;
end

% init.
NbitsFrac = 3;
Nstages = 8;
if numel(bCell)<Nstages
    Nstages =  numel(bCell);
end
Ratio = ones(1,Nstages);
delayFilt = ones(1,Nstages);
lengthFilt = ones(1,Nstages);
% NUps = 1;
% NDns = 1;
Nbits = [19 18 19 18 18 18 18 18];
flag_gain_compensation = 1;

DelayFilt = zeros(8,3);
DelayFiltSum = zeros(1,3);

if size(sigIn,1)<size(2)
    flag_row = 1;
    sigIn = sigIn.';
else
    flag_row = 0;
end
for i = 1:Nstages
    %     Ratio(i) = [NUps/NDns];
    %     fsIn = fsIn*Ratio(i);
    
    if (i~=2) && fsIn>fsOut/2
        NUps = 1;
        NDns = 1;
        delayFilt(i) = 0;
    elseif (i~=2) && fsIn<=fsOut
        NUps = 2;
        NDns = 1;
        bHB = DB_Taps_HB_Symmtric(cell2mat(bCell(i)), Nbits(i), fsIn*NUps/NDns, [fnum, Nstages, 1, i]);
        lengthFilt(i) = length(bHB);
        delayFilt(i) = floor(lengthFilt(i)/2);
        sigInIpo2 = INTPOL_ZERO_g(sigIn, 2, flag_gain_compensation);
        if flag_removeDelay
            sigOutFilt = conv(sigInIpo2, bHB, 'same');
        else
            sigOutFilt = conv(sigInIpo2, bHB);
        end
        PLOT_FFT_dB_g(sigOutFilt, fsIn*NUps, length(sigOutFilt), ['sig ipo2'], 'df', 'full', 'pwr', [fnum*10, Nstages, 1, i]);
        ylim([-200 10])
        
    elseif (i==2) && fsIn<fsOut
        if mod(fsOut/(fsIn*2),1)==0
            NUps = 2;
            NDns = 1;
        elseif mod(fsOut/(fsIn*4/3),1)==0
            NUps = 4;
            NDns = 3;
        elseif mod(fsOut/(fsIn*8/5),1)==0
            NUps = 8;
            NDns = 5;
        elseif mod(fsOut/(fsIn*16/15),1)==0
            NUps = 16;
            NDns = 15;
        end
        bPP = cell2mat(bCell(i))/2^Nbits(i);
        lengthFilt(i) = length(bPP)/NUps;
        delayFilt(i) = ceil( (lengthFilt(i)-1)*NUps/NDns/2 );
        sigOutFilt = FIR_PolyPhase_g(sigIn, bPP, [NUps, NDns], [], flag_gain_compensation, flag_removeDelay, 1);
        %         sigOutFiltTest = FIR_PolyPhase_g(sigIn(1+49:end-49), bPP, [NUps, NDns], [], flag_gain_compensation, flag_removeDelay, 1);
        PLOT_FFT_dB_g(sigOutFilt, fsIn*NUps/NDns, length(sigOutFilt), ['sig src'], 'df', 'full', 'pwr', [fnum*10, Nstages, 1, i]);
        ylim([-200 10])
        if 0
            LengthAdd = ceil( (length(bPP)/NUps - 1)*NUps/NDns )
            Delay = ceil(LengthAdd/2)
            Ni = mod( (numel(sigIn) + length(bPP)/NUps - 1)*NDns, NUps)
            if Ni==NUps-1
                sigOut_rmDelay = sigOutFilt(Delay+1:end-Delay);
            elseif Ni<NUps-1
                sigOut_rmDelay = sigOutFilt(Delay+1:end-(Delay-1) );
            else
                error('test?!')
            end
            sigOut_rmDelay_intp = interp1(1:NDns:NDns*numel(sigOut_rmDelay), sigOut_rmDelay, 1:1:NDns*numel(sigOut_rmDelay),'spline');
            fs = 23.04e6*2;
            PLOT_FFT_dB_g(sigOut_rmDelay_intp, fs*NUps/NDns*NDns, length(sigOut_rmDelay_intp), ['sigOut rmDelay intp'], 'df', 'full', 'pwr', [1012]);
            
            sigOut_rmDelay_intp_dmc = sigOut_rmDelay_intp(1:NUps:end);
            [evm] = dsp_evm_g(sigIn, sigOut_rmDelay_intp_dmc)
            [meaCor, DelayInt, DelayFrac] = dsp_Delay_EstCor_g2(sigIn, sigOut_rmDelay_intp_dmc, []);
        end
        
        flag_debug = 0;
        if flag_debug
            [~,~,plt1] = PLOT_FFT_dB_g(sigIn, fsIn, length(sigIn), ['1. sigIn'], 'df', 'full', 'pwr', [1004]);
            plt1.LineWidth = 3;
            sigIn_NUps = INTPOL_ZERO_g(sigIn, NUps,1);
            fs_ipo = fsIn*NUps;
            PLOT_FFT_dB_g(sigIn_NUps, fs_ipo, length(sigIn_NUps), ['2. sigIn Ipo',num2str(NUps)], 'df', 'full', 'pwr', [1004]);
            bPP_src_tmp = reshape(bPP, [], NUps).';
            bPP_src = bPP_src_tmp(:);
            yyaxis right
            PLOT_FFT_dB_g(bPP_src*length(sigIn), fs_ipo, length(sigIn), ['3. src filter'], 'df', 'full', 'pwr', [1004]);
            sigOutSRC = conv(sigIn_NUps, bPP_src, 'same');
            yyaxis left
            [~,~,plt2] = PLOT_FFT_dB_g(sigOutSRC, fs_ipo, length(sigOutSRC), ['4. sigOut polyphase filter'], 'df', 'full', 'pwr', [1004]);
            plt2.Color = [ 0.4940    0.1840    0.5560];
            sigOutDMC = sigOutSRC(1:NDns:end);
            fs_src = fsIn*NUps/NDns;
            [~,~,plt3] = PLOT_FFT_dB_g(sigOutDMC, fs_src, length(sigOutDMC), ['5. sigOut SRC:',num2str(NUps),'/',num2str(NDns)], 'df', 'full', 'pwr', [1004]);
            plt3.LineWidth = 0.5;
            plt3.LineStyle = '-';
            plt3.Color = [0.4660    0.6740    0.1880];
        end
    end
    
    flag_LengthAdd = 1;
    if nargout>=2 || flag_LengthAdd
        if mod(lengthFilt(i),2)==1 % odd
            DelayFilt(i,3) = delayFilt(i)+delayFilt(i);
            DelayFilt(i,1) = delayFilt(i);
            DelayFilt(i,2) = delayFilt(i);
        else % even
            DelayFilt(i,3) = delayFilt(i)+(delayFilt(i)-1);
            DelayFilt(i,1) = delayFilt(i);
            DelayFilt(i,2) = delayFilt(i)-1;
        end
        if flag_removeDelay
            DelayFilt = 0*DelayFilt;
        end
        
        if i~=2
            DelayFiltSum = DelayFiltSum*NUps/NDns + DelayFilt(i,:);
        elseif i==2 % SRC
            DelayFiltSum = ceil((DelayFiltSum + DelayFilt(i,:))*NUps/NDns);
        end
        if i~=2 && NUps~=1
            LengthSigOut = length(sigIn)*NUps+DelayFilt(i,3);
        elseif i==2 && NUps~=1 % SRC
            %             DelayFilt(i,:) = fix(DelayFilt(i,:)*NUps/NDns);
            LengthSigOut = ceil( (numel(sigIn)+lengthFilt(i)-1)*NUps/NDns );
        end
    end
    
    Ratio(i) = [NUps/NDns];
    fsIn = fsIn*Ratio(i);
    sigOutFilt = round(sigOutFilt,NbitsFrac);
    sigIn = sigOutFilt;
    
    if flag_row
        sigOutFiltCell{i} = sigOutFilt.';
    else
        sigOutFiltCell{i} = sigOutFilt;
    end
end
end



function bHB = DB_Taps_HB_Symmtric(bHalf_DB, Nbits, fs, fnum, fnum_legend)

%         bHalf_DB = [4,-10,22,-44,80,-134,214,-330,492,-714,1006,-1390,1886,-2516,3312,-4316,5580,-7188,9262,-12022,15878,-21710,31824,-54684,166568];
Nhalf = length(bHalf_DB);
%         Nbits = 19;

if exist('fnum','var')&&~isempty(fnum)
    flag_fnum=1;
else
    flag_fnum=0;
end
if exist('fnum_legend','var')&&~isempty(fnum_legend)
    flag_fnum_legend=1;
else
    flag_fnum_legend=0;
end
if ~exist('Nbits','var')||isempty(Nbits)
    Nbits=0;
end

bHB = zeros(2*2*Nhalf-1,1);
bHB(1:2:2*Nhalf) = bHalf_DB;
bHB(2*Nhalf+0) = 2^(Nbits-1);
bHB((2*Nhalf+1):2:end) = bHalf_DB(end:-1:1);
bHB = bHB/2^Nbits;

if flag_fnum
    Nfft = 2048;
    %     fs = 23.04e6
    %     PLOT2_FFT_dB(bHB*Nfft, fs, Nfft, fnum),hold on
    [~,~,plt] = PLOT_FFT_dB_g(bHB*Nfft, fs, Nfft, ['FIR'], 'df', 'full', 'pwr', [fnum]);
    if flag_fnum_legend
        plt.DisplayName = fnum_legend;
    end
end
end