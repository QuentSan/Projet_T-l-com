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
        Q = 2*bits_Q -1 ;

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

        Px = mean(abs(xe).^2)/2;
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
        bits_Q_decides = (Q_decides +1)/2; % On repasse sur des bits

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