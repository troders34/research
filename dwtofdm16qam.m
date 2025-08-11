function dwt_ofdm_rs_full
% End-to-end: RS -> QAM (BPSK/QPSK/16QAM) -> IDWT -> AWGN -> DWT -> Demod -> RS
% Catatan:
% - IDWT menggantikan IFFT (wavelet-OFDM).
% - DWT MATLAB real -> proses I/Q terpisah lalu digabung kompleks.
% - qammod/qamdemod dipakai dengan InputType/OutputType 'bit'.

%% ====== Parameter Utama ======
modName   = "16QAM";     % "BPSK", "QPSK", "16QAM"
M         = modM(modName);         % 2 / 4 / 16
k         = log2(M);               % bit per simbol
EbN0dB    = 12;                    % uji SNR
wname     = 'db4';                 % wavelet: 'haar','db2','db4',...
Lev       = 5;                     % level dekomposisi
Nw        = 1024;                  % panjang blok koefisien/time-domain
useCP     = true;  cpLen = 64;     % CP opsional (mulailah pakai CP di multipath)
N_RSblk   = 60;                    % jumlah blok RS untuk uji

% RS(255,239) byte-wise (GF(2^8))
Nrs = 255; Krs = 239; B = 8;       % bit per byte

%% ====== Data & Reed-Solomon Encode ======
bitsIn = randi([0 1], N_RSblk*Krs*B, 1, 'logical');   % bit mentah

[bi2u8, u82bi] = bit_packers(B);
enc = comm.RSEncoder('CodewordLength',Nrs,'MessageLength',Krs, ...
                     'GeneratorPolynomialSource','Auto','BitInput',false);

% Kemasi ke byte, bentuk blok RS, lalu encode
u8_in     = bi2u8(bitsIn);
u8_in_blk = reshape(u8_in, Krs, []);      % Krs bytes per codeword
u8_enc_blk = zeros(Nrs, size(u8_in_blk,2), 'uint8');
for i=1:size(u8_in_blk,2)
    u8_enc_blk(:,i) = enc(u8_in_blk(:,i));
end
u8_enc = u8_enc_blk(:);
bitsEnc_RS = u82bi(u8_enc);

%% ====== Modulasi QAM (bit-I/O) ======
% Pastikan panjang bit kelipatan k dan per-frame (Nw simbol)
bitsPerFrame = Nw*k;
needPadBits  = mod(-numel(bitsEnc_RS), k);         % kelipatan k
bitsModIn    = [bitsEnc_RS; false(needPadBits,1)]; % padding 0 untuk modulasi

% Pastikan kelipatan Nw*k agar frame penuh
needPadFrame = mod(-numel(bitsModIn), bitsPerFrame);
bitsModIn    = [bitsModIn; false(needPadFrame,1)];

bitsModIn = uint8(bitsModIn);
symTxAll  = qammod(bitsModIn, M, 'InputType', 'bit', ...
    'UnitAveragePower', true);

% Bagi ke frame DWT-OFDM
nFrames = numel(symTxAll)/Nw;
symTxAll = reshape(symTxAll, Nw, nFrames);

%% ====== Bookkeeping DWT ======
dwtmode('per','nodisp');                % periodization -> panjang & energi rapi
[~, book] = wavedec(zeros(Nw,1), Lev, wname);

%% ====== Transmitter (IDWT / waverec) ======
txSig = [];
for f = 1:nFrames
    sblk = symTxAll(:,f);
    x = dwt_ofdm_tx_frame(sblk, book, wname);     % IDWT I/Q -> time-domain
    if useCP && cpLen>0, x = [x(end-cpLen+1:end); x]; end
    % Headroom agar aman ke perangkat
    x = x / max(abs(x)+1e-12) * 0.85;
    txSig = [txSig; x]; %#ok<AGROW>
end

%% ====== Kanal ======
% AWGN saja (mulai dari sini). Tambah multipath Rayleigh setelah lolos AWGN.
rxSig = awgn(txSig, EbN0dB + 10*log10(k), 'measured');

%% ====== Receiver (DWT / wavedec) ======
symRxAll = complex(zeros(Nw,nFrames));
ptr = 1;
for f = 1:nFrames
    if useCP && cpLen>0
        seg = rxSig(ptr+cpLen : ptr+cpLen+Nw-1);
        ptr = ptr + cpLen + Nw;
    else
        seg = rxSig(ptr : ptr+Nw-1);
        ptr = ptr + Nw;
    end
    symRxAll(:,f) = dwt_ofdm_rx_frame(seg, Lev, wname); % koefisien (simbol)
end

symRx = symRxAll(:);

%% ====== Demodulasi QAM (bit) ======
bitsRxEnc_all = qamdemod(symRx, M, 'OutputType','bit', ...
                         'UnitAveragePower', true);

% Buang padding modulasi & padding frame
bitsRxEnc = bitsRxEnc_all(1:numel(bitsModIn)-needPadFrame);
bitsRxEnc = bitsRxEnc(1:numel(bitsEnc_RS));  % kembali ke panjang RS-encoded

%% ====== RS Decode ======
dec = comm.RSDecoder('CodewordLength',Nrs,'MessageLength',Krs, ...
                     'GeneratorPolynomialSource','Auto','BitInput',false, ...
                     'NumCorrectedErrorsOutputPort',true);

u8_rx_enc = bi2u8(bitsRxEnc);
nCW = floor(numel(u8_rx_enc)/Nrs);
u8_rx_enc = reshape(u8_rx_enc(1:nCW*Nrs), Nrs, []);

u8_dec_blk = zeros(Krs, nCW, 'uint8');
uncorr     = false(1,nCW);   % indikator ada error tak terkoreksi
for i=1:nCW
    [u8_dec_blk(:,i), uncorr] = dec(u8_rx_enc(:,i));
    uncorr(i) = (uncorr < 0);
end
u8_dec = u8_dec_blk(:);
bitsOut = u82bi(u8_dec);

% Potong ke panjang data asli
bitsOut = bitsOut(1:numel(bitsIn));

%% ====== Evaluasi ======
[bitErr, BER] = biterr(bitsIn, bitsOut);
FER = mean(uncorr);   % fraksi codeword dengan error tak terkoreksi

fprintf('Modulasi     : %s (M=%d)\n', modName, M);
fprintf('Wavelet/Lev  : %s / %d\n', wname, Lev);
fprintf('Nw / CP      : %d / %d\n', Nw, (useCP*cpLen));
fprintf('Eb/N0 (dB)   : %.1f\n', EbN0dB);
fprintf('BER akhir    : %.3e  (errors=%d dari %d bits)\n', BER, bitErr, numel(bitsIn));
fprintf('RS FER       : %.6f\n', FER);

% Sanity check rekonstruksi ideal (tanpa noise, tanpa CP)
s0 = symTxAll(:);
x0 = dwt_ofdm_tx_frame(s0(1:Nw), book, wname);
y0 = dwt_ofdm_rx_frame(x0, Lev, wname);
fprintf('Sanity ||s-y||/||s|| : %.3e\n', norm(s0(1:Nw)-y0)/norm(s0(1:Nw)));

end % ====== end main ======

%% ====== Helper: modulasi M dari nama ======
function M = modM(name)
    switch upper(string(name))
        case "BPSK",  M = 2;
        case "QPSK",  M = 4;
        case "16QAM", M = 16;
        otherwise, error('Pilihan modulasi tidak dikenal.');
    end
end

%% ====== Helper: pack/unpack bit <-> uint8 ======
function [b2u8, u82b] = bit_packers(B)
    b2u8 = @(bits) uint8(bi2de(reshape(bits, B, []).','left-msb'));
    u82b = @(u8) reshape(de2bi(double(u8(:)), B, 'left-msb').', [], 1)>0;
end

%% ====== TX frame: IDWT I/Q ======
function x = dwt_ofdm_tx_frame(symBlk, book, wname)
    % symBlk: Nw x 1 kompleks (koefisien / simbol)
    cI = real(symBlk(:));
    cQ = imag(symBlk(:));
    xI = waverec(cI, book, wname);
    xQ = waverec(cQ, book, wname);
    x  = complex(xI, xQ);
    rx = rms(x);
    if rx>0, x = x/rx; end
end

%% ====== RX frame: DWT I/Q ======
function s_hat = dwt_ofdm_rx_frame(x, Lev, wname)
    rI = real(x(:));  rQ = imag(x(:));
    cI = wavedec(rI, Lev, wname);
    cQ = wavedec(rQ, Lev, wname);
    % Jumlah koefisien = panjang sinyal
    L = numel(rI);
    s_hat = complex(cI(1:L), cQ(1:L));
end
