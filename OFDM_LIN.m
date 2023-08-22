clear all
clc
snr_dB = -5:1:20;
snr = zeros(1,length(snr_dB));
for q = 1:1:length(snr_dB)
    snr(q) = 10^(snr_dB(q)/10); 
end

 % rx_grid = zeros(M,N);
 BER = zeros(1,length(snr));
    M  = 128;
    l = 10;
    index = 128;
    k = log2(index);
    n = 1280;
    cp = l;
for j = 1:1:length(snr)
    avg_ber = zeros(1,100);
    for loop= 1:1:100
        x = randi([0,1],n*k,1); %data
        in_symbols = bit2int(x,k);
        x_mod = qammod(in_symbols,index,'bin');
        N = length(x_mod)/M; % time instants
        rx_grid = zeros(M,N); 
        h = sqrt(1/2)*(randn(1,l)+1i*randn(1,l));
        H = fft(h,M);
        %h = [h zeros(1,M-cp)];
        grid = reshape(x_mod,[M,N]);
        noise = sqrt(1/(2*snr(j)))*(randn(1,M)+1i*randn(1,M));
        for nn = 1:1:N
            tx_sym = grid(:,nn);
            tx_sym_ifft = ifft(tx_sym,M);
            tx_sym_ifft = tx_sym_ifft.';
            tx_cp = [tx_sym_ifft(M-cp+1:M) tx_sym_ifft];
            rx = conv(h,tx_cp);
            %rx_with_noise = awgn(rx,snr(j),"measured");
            rx_no_cp = rx(cp+1:cp+M)+noise;
            rx_final = fft(rx_no_cp,M);%change
            rx_equ = rx_final./H;
            rx_grid(:,nn) = rx_equ;
        end
    y_demod = reshape(rx_grid,[M*N,1]);
    outsymbols = qamdemod(y_demod,index,'bin');
    y = int2bit(outsymbols,k);
    [numErrors,ber] = biterr(x,y);
    avg_ber(loop) = ber; 
    end
    BER(j) = double(mean(avg_ber));
end
%hold on;
semilogy(snr_dB,BER,'r');
%grid on;
xlabel('snr [dB]');
ylabel('Bit Error Rate');
%xlim([-5 20]);
%ylim([10^(-5) 0.45]);
title('OFDM');