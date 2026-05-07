% PROJET DE TELECOMUNICATIONS
%Etude d'un système de transmission en bande passante

clear all
close all

%==========================================================================
% Paramètres
%==========================================================================

M = 4;          % 4-ASK donc 4 symboles
Fs = 12000;     % fréquence d'échantillonnage donnée
Rb = 6000;      % débit binaire donné
Rs = Rb/log2(M); % débit de symbole, étant donné qu'on a 2 bits par symbole
Ts = 1/Rs;      % durée d'un symbole
Ns = Fs/Rs;     % nombre d'échantillons par symbole
fc = 2000;      % fréquence porteuse donnée
alpha = 0.35;   % facteur de roll-off donné
N_bit_block = 2000; % nombre de bits générés par paquets dans la simulation pour ne pas saturer
tab_EbN0_dB = 0:8;  % valeurs de SNR à tester
span = 10;      % durée filtre en symboles

% --- FILTRES ---
h = rcosdesign(alpha, span, Ns, 'sqrt');
hr = fliplr(h);

%==========================================================================
% QUESTION 3 - 4-ASK
%==========================================================================

% =========================================================
% TEST SANS BRUIT : vérification que le TEB = 0
% =========================================================
bits_test = randi([0 1], 1, N_bit_block);

% -- Mapping 4-ASK --
% On regroupe les bits 2 par 2 pour former les symboles
% Mapping Gray : 00->-3, 01->-1, 11->+1, 10->+3
bits_matrix_test = reshape(bits_test, 2, []);  % 2 bits par symbole
b1_test = bits_matrix_test(1,:);               % bit de poids fort
b2_test = bits_matrix_test(2,:);               % bit de poids faible

% - Mapping Gray sur {-3, -1, +1, +3} -
symboles_test = zeros(1, size(bits_matrix_test, 2));
for i = 1:length(symboles_test)
    idx = b1_test(i)*2 + b2_test(i);
    switch idx
        case 0, symboles_test(i) = -3; % 00 -> -3
        case 1, symboles_test(i) = -1; % 01 -> -1
        case 3, symboles_test(i) = +1; % 11 -> +1
        case 2, symboles_test(i) = +3; % 10 -> +3
    end
end

% -- Modulation --
% La 4-ASK est une modulation purement réelle (pas de composante Q)
diracs_test = kron(symboles_test, [1 zeros(1, Ns-1)]);
xe_test = filter(h, 1, diracs_test);
t_test = (0:length(xe_test)-1)/Fs;
x_test = real(xe_test) .* cos(2*pi*fc*t_test);

% -- Démodulation sans bruit --
% On multiplie par 2*cos pour compenser le facteur 1/2 issu du produit cos*cos
r_bb_test = 2 * x_test .* cos(2*pi*fc*t_test);
sf_test = filter(hr, 1, r_bb_test);
t0 = span*Ns + 1; % +1 pour éviter le résidu à 2fc (fc=2000Hz proche de Rs)
se_test = real(sf_test(t0:Ns:end));

% -- Décision avec seuils à -2, 0, +2 --
sym_decides_test = zeros(1, length(se_test));
for i = 1:length(se_test)
    if     se_test(i) < -2, sym_decides_test(i) = -3;
    elseif se_test(i) < 0,  sym_decides_test(i) = -1;
    elseif se_test(i) < 2,  sym_decides_test(i) = +1;
    else,                    sym_decides_test(i) = +3;
    end
end

% -- Démapping Gray : {-3,-1,+1,+3} -> bits --
bits_decides_test = zeros(1, 2*length(sym_decides_test));
for i = 1:length(sym_decides_test)
    switch sym_decides_test(i)
        case -3, b1=0; b2=0; % 00
        case -1, b1=0; b2=1; % 01
        case +1, b1=1; b2=1; % 11
        case +3, b1=1; b2=0; % 10
    end
    bits_decides_test(2*i-1) = b1;
    bits_decides_test(2*i)   = b2;
end

% -- Vérification --
n_cmp = min(length(bits_test), length(bits_decides_test));
TEB_sans_bruit = sum(bits_test(1:n_cmp) ~= bits_decides_test(1:n_cmp)) / n_cmp;
fprintf('TEB sans bruit = %e (doit être 0)\n', TEB_sans_bruit);


% =========================================================
% PARTIE AVEC BRUIT
% =========================================================

% ---- Boucle 1 sur les valeurs de Eb/N0 pour tester les différents SNR
TEB_simu = zeros(1, length(tab_EbN0_dB));
for indice = 0:length(tab_EbN0_dB)-1

    EbN0_dB = tab_EbN0_dB(indice+1) % valeur testée qu'on affiche dans la console

    % - initialisation des compteurs
    nb_erreurs = 0;
    compteur   = 0;
    TEB_cumul  = 0;


    % ---- Boucle 2 sur le nombre d'erreurs

    while nb_erreurs < 1000

        % 1 --- MODULATEUR ---

        % La 4-ASK encode 2 bits par symbole sur l'axe réel uniquement.
        % On module en amplitude sur la porteuse cosinus.

        % 1.1 -- Mapping 4-ASK --
        bits = randi([0 1], 1, N_bit_block); % génération de bits aléatoires

        % - Regroupement des bits 2 par 2 -
        bits_matrix = reshape(bits, 2, []);  % 2 bits par symbole
        b1 = bits_matrix(1,:);               % bit de poids fort
        b2 = bits_matrix(2,:);               % bit de poids faible

        % - Mapping Gray sur {-3, -1, +1, +3} -
        symboles = zeros(1, size(bits_matrix, 2));
        for i = 1:length(symboles)
            idx = b1(i)*2 + b2(i);
            switch idx
                case 0, symboles(i) = -3; % 00 -> -3
                case 1, symboles(i) = -1; % 01 -> -1
                case 3, symboles(i) = +1; % 11 -> +1
                case 2, symboles(i) = +3; % 10 -> +3
            end
        end

        % 1.2 -- Suréchantillonnage --
        diracs = kron(symboles, [1 zeros(1, Ns-1)]); % génération de notre train de symboles avec impulsion Dirac

        % - Filtrage de mise en forme -
        xe = filter(h, 1, diracs);

        % 1.3 -- Translation sur la fréquence porteuse --
        % - Axe pour le temps -
        t = (0:length(xe)-1)/Fs;

        % - Modulation sur l'onde porteuse (4-ASK : axe réel uniquement) -
        x = real(xe) .* cos(2*pi*fc*t);

        % 1.4 -- Canal bruité AWGN --
        Px = mean(abs(x).^2);
        Eb_sur_N0 = 10^(EbN0_dB/10); % conversion
        sigma2 = (Px*Ns) / (2*log2(M)*Eb_sur_N0);
        bruit = sqrt(sigma2)*randn(1, length(x)); % notre bruit gaussien

        % - le signal bruité de notre canal -
        r = x + bruit;


        % 2 --- DEMODULATEUR ---

        % 2.1 -- Retour en bande de base --
        % On multiplie par 2*cos pour compenser le facteur 1/2 issu du produit cos*cos
        r_bb = 2 * r .* cos(2*pi*fc*t);

        % 2.2 -- Filtrage adapté --
        signal_filtre = filter(hr, 1, r_bb);
        t0 = span*Ns + 1; % +1 pour éviter le résidu à 2fc
        signal_echantillone = real(signal_filtre(t0:Ns:end)); % on récupère les symboles transmis

        % 2.3 -- Décision avec seuils à -2, 0, +2 --
        sym_decides = zeros(1, length(signal_echantillone));
        for i = 1:length(signal_echantillone)
            if     signal_echantillone(i) < -2, sym_decides(i) = -3;
            elseif signal_echantillone(i) < 0,  sym_decides(i) = -1;
            elseif signal_echantillone(i) < 2,  sym_decides(i) = +1;
            else,                                sym_decides(i) = +3;
            end
        end

        % -- Démapping Gray : {-3,-1,+1,+3} -> bits --
        bits_decides = zeros(1, 2*length(sym_decides));
        for i = 1:length(sym_decides)
            switch sym_decides(i)
                case -3, b1_d=0; b2_d=0; % 00
                case -1, b1_d=0; b2_d=1; % 01
                case +1, b1_d=1; b2_d=1; % 11
                case +3, b1_d=1; b2_d=0; % 10
            end
            bits_decides(2*i-1) = b1_d;
            bits_decides(2*i)   = b2_d;
        end

        % 3 --- COMPTAGE DES ERREURS ---
        n = min(length(bits), length(bits_decides));
        err = sum(bits(1:n) ~= bits_decides(1:n));
        nb_erreurs = nb_erreurs + err;
        TEB_cumul  = TEB_cumul + err/N_bit_block;
        compteur   = compteur + 1;

    end

    TEB_simu(indice+1) = TEB_cumul / compteur;

end

% 4 --- TEB THEORIQUE POUR 4-ASK ---
% Formule du cours (eq. 264) :
% Pb = 2*(M-1)/(M*log2(M)) * Q( sqrt( 3*log2(M)/(2*M^2-3*M+1) * Eb/N0 ) )
% Pour M=4 : 2*(4-1)/(4*2) = 3/4  et  3*log2(4)/(2*16-12+1) = 6/21
EbN0_lin = 10.^(tab_EbN0_dB/10);
TEB_th = (3/4) * qfunc(sqrt((6/21) * EbN0_lin));

% 5 --- PLOTS ---

% 1 — Signal après pulse shaping
figure
plot(real(xe(1:200)))
title('Signal xe(t) après pulse shaping (4-ASK)')
xlabel('Echantillons'); ylabel('Amplitude')

% 2 — Signal bandpass x(t)
figure
plot(t(1:500), x(1:500))
title('Signal x(t) sur porteuse fc=2000 Hz (4-ASK)')
xlabel('Temps (s)'); ylabel('Amplitude')

% 3 — DSP de x(t)
figure
pwelch(x, [], [], [], Fs)
title('DSP de x(t) - 4-ASK')

% 4 — Courbe BER
figure
semilogy(tab_EbN0_dB, TEB_th, 'r-o')
hold on
semilogy(tab_EbN0_dB, TEB_simu, 'b-+')
xlabel('Eb/N0 (dB)')
ylabel('BER')
legend('Théorique', 'Simulé')
grid on
title('BER 4-ASK')

