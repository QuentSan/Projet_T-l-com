% PROJET DE TELECOMUNICATIONS
%Etude d'un système de transmission en bande passante

clear all
close all

%==========================================================================
% Paramètres
%==========================================================================

M = 4 ; % QPSK donc 4 symboles
Fs = 10000; %fréquence d'échantillonage donnée
Rb = 2000; % débit binaire donné
Rs = Rb/log2(M); % débit de symbole, étant donné qu'on a 2 bits par symbole
Ts = 1/Rs; %durée d'un symbole
Ns = Fs/Rs; % nombre d'échantillons par symbole
fc = 2000; %fréquence porteuse donnée
alpha = 0.35; %facteur de roll-off donné
N_bit_block = 1000; % nombre de bits générer par paquets dans la simulation pour ne pas satuerer
tab_EbN0_dB = 0:8; % valeurs de SNR à tester
span = 10; %durée filtre en symboles

% --- FILTRES ---
h = rcosdesign(alpha, span, Ns, 'sqrt');
hr = fliplr(h);

%==========================================================================
% QUESTION 1 
%==========================================================================

% =========================================================
% TEST SANS BRUIT : vérification que le TEB = 0
% =========================================================
bits_test = randi([0 1], 1, N_bit_block);

% -- Mapping --
bits_I_test_tx = bits_test(1:2:end);  % On envoie les impairs sur les réels
bits_Q_test_tx = bits_test(2:2:end);  % On envoie les pairs sur les imaginaires

% - Mapping sur -1 et 1 (inversion de Q pour respecter le code Gray) -
I_test = 2*bits_I_test_tx - 1;
Q_test = 1 - 2*bits_Q_test_tx;  % inversion : 0→+1, 1→-1 pour Gray

% - Mapping symbole complexe -
symboles_test = I_test + 1j*Q_test;

% -- Modulation --
diracs_test = kron(symboles_test, [1 zeros(1, Ns-1)]);
xe_test = filter(h, 1, diracs_test);
xI_test = real(xe_test);
xQ_test = imag(xe_test);
t_test = (0:length(xe_test)-1)/Fs;
x_test = xI_test .* cos(2*pi*fc*t_test) - xQ_test .* sin(2*pi*fc*t_test);

% -- Démodulation sans bruit --
rI_test = x_test .* cos(2*pi*fc*t_test);
rQ_test = x_test .* (-sin(2*pi*fc*t_test));
r_bb_test = rI_test + 1j*rQ_test;
sf_test = filter(hr, 1, r_bb_test);
t0 = span * Ns;
se_test = sf_test(t0:Ns:end);

% -- Décision et démapping --
I_dec_test = sign(real(se_test));  % Pour repasser sur des -1 et 1 pour la partie réelle des valeurs
Q_dec_test = sign(imag(se_test));  % Pour repasser sur des -1 et 1 pour la partie imaginaire des valeurs

bits_I_decides_test = (I_dec_test + 1)/2;   % On repasse sur des bits : +1→1, -1→0
bits_Q_decides_test = (1 - Q_dec_test)/2;   % Inversion miroir du mapping : +1→0, -1→1

% - Recombinaison des bits pour reconstituer le signal original
n_sym_test = length(bits_I_decides_test);
bits_decides_test = zeros(1, 2*n_sym_test);  % On prépare notre vecteur pour accueillir notre signal recomposé
bits_decides_test(1:2:end) = bits_I_decides_test;  % On replace les impairs
bits_decides_test(2:2:end) = bits_Q_decides_test;  % On replace les pairs

% -- Vérification --
n_cmp = min(length(bits_test), length(bits_decides_test));
TEB_sans_bruit = sum(bits_test(1:n_cmp) ~= bits_decides_test(1:n_cmp)) / n_cmp;
fprintf('TEB sans bruit = %e (doit être 0)\n', TEB_sans_bruit);


% =========================================================
% PARTIE AVEC BRUIT
% =========================================================

% ---- Boucle 1 sur les valeurs de de Eb/N0 pour tester les différents SNR
TEB_simu = zeros(1, length(tab_EbN0_dB));
for indice = 0: length(tab_EbN0_dB)-1

    EbN0_dB = tab_EbN0_dB(indice +1) %valeur testée qu'on affiche dans la console
  

    % - initialisation des compteurs
    nb_erreurs = 0;
    compteur = 0;
    TEB_cumul =0;

    

    % ---- Boucle 2 sur le nombre d'erreurs 

    while nb_erreurs < 1000

        % 1 --- MODULATEUR ---

        % Le principe de la QPSK étant de regrouper les bits 2 par deux, on va donc
        % avoir des symboles complexes. Pour cela on utilise un cosinus pour
        % transmettre sur les réels et un sinus pour transmettre sur les
        % imaginaires. Ceci nous permet ainsi d'envoyer dans un repère orthogonal
        % et donc de transmettre nos bits deux à deux correctements.
        
        % 1.1 -- Mapping QPSK --
        bits = randi([0 1], 1 , N_bit_block);  % génération de bits aléatoires
        
        % - Sépartion réels et imaginaires -
        
        bits_I = bits(1:2:end);  % On envoie les impairs sur les réels
        bits_Q = bits(2:2:end);  % On envoie les pairs sur les imaginaires
        
        % - Mapping sur -1 et 1 -

        I = 2*bits_I -1;
        Q = 1 -2*bits_Q  ;

        % - Maping  symbole complexe -
        symboles = I + 1j*Q;

        % 1.2 --Suréchantillonnage --

        diracs = kron(symboles, [1 zeros(1, Ns-1)]); % génération de notre train de symboles avec impulsion Dirac
        
        % - filtrage de mise en forme -
        xe = filter(h,1,diracs);

        % - Séparation I et Q après mise en forme -
        xI = real(xe); % pour les réels
        xQ = imag(xe); % pour les imaginaires
        
        % 1.3 -- Translation sur les fréquences --
        % - Axe pour le temps -
        t = (0:length(xe)-1)/Fs;

        %- Modulation sur notre onde porteuse -

        x= xI .* cos(2*pi*fc*t) - xQ .* sin(2*pi*fc*t);

        % 1.4 -- Canal bruité AWGN --

        Px = mean(abs(x).^2);
        Eb_sur_N0 = 10^(EbN0_dB/10); %conversion
        sigma2 = (Px*Ns) / (2*log2(M)*Eb_sur_N0);
        bruit = sqrt(sigma2)*randn(1, length(x)); % notre bruit gaussien

        % - le signal bruité de notre canal -
        r = x + bruit;
        

        % 2 --- DEMODULATEUR ---

        
        % 2.1 -- Retour en bande de base --

        rI = r.* cos(2*pi*fc*t); % On récupère la partie réelle
        rQ = r.* (-sin(2*pi*fc*t)); % On récupère la partie imaginaire
        r_bb = rI + 1j*rQ; % On recompose le signal complet dans la bande de base


        % 2.2 -- filtrage --

        signal_filtre = filter(hr, 1, r_bb); % On applique le filtre inverse sur notre signal en bande de base
        t0 = span * Ns;
        signal_echantillone = signal_filtre(t0:Ns:end); % on récupère les symboles transmis

        %2.3 -- Décisiion et démapping --

        I_decides = sign(real(signal_echantillone)); % Pour repasser sur des -1 et 1 pour la partie réelle des valeur
        Q_decides = sign(imag(signal_echantillone)); % Pour repasser sur des -1 et 1 pour la partie imaginaire des valeur

        bits_I_decides = (I_decides +1)/2; % On repasse sur des bits
        bits_Q_decides = (1-Q_decides )/2; % On repasse sur des bits

        % - Etape de recombinaison des bits conetenus dans la partie réelle
        % et imaginaire pour reconstituer le signal original
        n_sym = length(bits_I_decides);
        bits_decides = zeros(1, 2*n_sym); % On prépare notre vecteur pour acceuillir notre sigan recomposé
        %bits_decides = zeros(1, N_bit_block); % On prépare notre vecteur pour acceuillir notre sigan recomposé
        bits_decides(1:2:end) = bits_I_decides; %On replace les impairs
        bits_decides(2:2:end) = bits_Q_decides; %On replace les pairs

        % 3 --- COMPTAGE DES ERREURS ---

        n = min(length(bits), length(bits_decides));
        err = sum(bits(1:n) ~= bits_decides(1:n));
        nb_erreurs = nb_erreurs + err;
        TEB_cumul  = TEB_cumul + err/N_bit_block;
        compteur   = compteur + 1;

    end
        
    TEB_simu(indice +1) = TEB_cumul / compteur; 
        
end

% 4 --- TEB THEORIQUE POUR QPSK ---

EbN0_lin = 10.^(tab_EbN0_dB/10);
TEB_th = 0.5*erfc(sqrt(EbN0_lin));

% 5 --- PLOTS ---
% 1 — Signaux I et Q après pulse shaping
figure
subplot(2,1,1)
plot(xI(1:200))
title('Voie I après pulse shaping')
xlabel('Echantillons'); ylabel('Amplitude')

subplot(2,1,2)
plot(xQ(1:200))
title('Voie Q après pulse shaping')
xlabel('Echantillons'); ylabel('Amplitude')

% 2 — Signal bandpass x(t)
figure
plot(t(1:500), x(1:500))
title('Signal x(t) sur porteuse fc=2000 Hz')
xlabel('Temps (s)'); ylabel('Amplitude')

% 3 — DSP de x(t)
figure
pwelch(x, [], [], [], Fs)
title('DSP de x(t)')

% 4 — Courbe BER
figure
semilogy(tab_EbN0_dB, TEB_th, 'r-o')
hold on
semilogy(tab_EbN0_dB, TEB_simu, 'b-+')
xlabel('Eb/N0 (dB)')
ylabel('BER')
legend('Théorique', 'Simulé')
grid on
title('BER QPSK')