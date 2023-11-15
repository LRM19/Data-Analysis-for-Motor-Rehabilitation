%% INIZIALIZZAZIONE

close all
clear all
clc

% Richieste prliminari all'utente, che decide il metodo di risoluzione
% della cinematica ed impone i criteri di selezione dei dati

prompt = {'Inserisci il tuo codice persona: ', ...
    sprintf(['Premi R se vuoi risolvere utilizzando il Teorema di Rivals\n' ...
    'Premi C se vuoi risolvere utilizzando la Chiusura in Posizione']), ...
    'Premi L se vuoi visualizzare i diagrammi della FDT in scala logaritmica'};
temp = inputdlg(prompt, 'Input per la risoluzione');

CODICEPERSONA = round(str2double(temp{1}(length(temp{1}))));
P = round(str2double(temp{1}(length(temp{1})-1)));
scelta = temp{2};
key = temp{3};

% IMPORTO DAL FILE 'data_stud.mat' LA TABELLA CORRETTA
switch CODICEPERSONA
    case 1  
        temp = 'uno';
    case 2
        temp = 'due';
    case 3
        temp = 'tre';
    case 4
        temp = 'quattro';
    case 5
        temp = 'cinque';
    case 6
        temp = 'sei';
    case 7
        temp = 'sette';
    case 8
        temp = 'otto';
    case 9 
        temp = 'nove';
    case 0
        temp = 'zero';
end

DATI = load("data_stud.mat").(temp);
clear temp;

dt = 1/20; %[s] , frequenza di raccolta dati da parte dei marker
tempo = 0:dt:((size(DATI,1)-1)*dt); %[s]

L12_vert = 88; % [mm] M2 si trova 88mm più in basso rispetto ad M1 

%% 3. CINEMATICA

%% STEP 3.1 : CALCOLO DELLA POSIZIONE DI M2 (x,y)
for ii = 1:length(DATI)
    M2X(ii) = DATI(ii,1);
    M2Y(ii) = DATI(ii,2) - L12_vert;
end
M2X = M2X';
M2Y = M2Y';

% Ricostruisco la tabella dei dati aggiungendo le posizioni di M2
DATI = [DATI(:,1), DATI(:,2), M2X, M2Y, DATI(:,3), DATI(:,4), DATI(:,5), DATI(:,6), DATI(:,7), DATI(:,8)];

%% STEP 3.2 : CALCOLO DELLE LUNGHEZZE M1-M2, M1-M3, M2-M4, M3-4

for jj = 1:length(DATI)
    M1M2(jj) = sqrt((DATI(jj,3)-DATI(jj,1))^2 + (DATI(jj,4)-DATI(jj,2))^2);
    M1M3(jj) = sqrt((DATI(jj,5)-DATI(jj,1))^2 + (DATI(jj,6)-DATI(jj,2))^2);
    M2M4(jj) = sqrt((DATI(jj,7)-DATI(jj,3))^2 + (DATI(jj,8)-DATI(jj,2))^2);
    M3M4(jj) = sqrt((DATI(jj,7)-DATI(jj,5))^2 + (DATI(jj,8)-DATI(jj,6))^2);
    M4M5(jj) = sqrt((DATI(jj,9)-DATI(jj,7))^2 + (DATI(jj,10)-DATI(jj,8))^2);
end

M1M2 = mean(M1M2);
M1M3 = mean(M1M3);
M2M4 = mean(M2M4);
M3M4 = mean(M3M4);
M4M5 = mean(M4M5);

%% STEP 3.4.1 : TABELLE ANTROPOMETRICHE PER BARICENTRO G1

h = 1670; %[mm]

% G1 = baricentro del braccio rispetto alla spalla
L_BRACCIO = 0.186*h; %[mm]
G1 = 0.436*L_BRACCIO; %[mm]

%% STEP 3.4.2 : TABELLE ANTROPOMETRICHE PER BARICENTRO G2

% G2 = baricentro dell'avambraccio rispetto al polso
L_AVAMBRACCIO = 0.146*h; %[mm]
G2 = 0.570*L_AVAMBRACCIO; %[mm]

%% STEP 3.4.3.1 : RAPPRESENTAZIONE TRAIETTORIE G1 G2 NEL TEMPO

video_handle = figure('Name', 'SCHEMA VIDEO DEL MOVIMENTO','NumberTitle','off');
writeObj=VideoWriter('G1G2vid.avi');
writeObj.FrameRate=10; 
open(writeObj)
for kk = 1:length(tempo)
    try 
        delete(A1);
        delete(A2);
        delete(A3);
        delete(A4);
        delete(A5);
        delete(G1_time);
        delete(G2_time);
        delete(G1anchor);
        delete(G2anchor);
        delete(M1_time);
        delete(M2_time);
        delete(M3_time);
        delete(M4_time);
        delete(M5_time);
        delete(M1anchor);
        delete(M2anchor);
        delete(M3anchor);
        delete(M4anchor);
        delete(M5anchor);
    end

    % Asse dell'esoscheletro
    ASSEY_exo = line([mean(DATI(:,1)), mean(DATI(:,3))],[-175, 550]);
    set(ASSEY_exo, 'color', 'k', 'linewidth', 1.5);

    %Linea di terra (piano su cui sta seduto il soggetto)
    ASSEX_exo = line([mean(DATI(:,1))-25, max(DATI(:,9))-75],[-150, -150]);
    set(ASSEX_exo, 'color', 'k', 'linewidth', 1.5);

    % A1 = asta M1-M2
    A1 = line([DATI(kk,1), DATI(kk,3)],[DATI(kk,2), DATI(kk,4)]);
    set(A1, 'color', 'r', 'linewidth', 0.5);
    
    % A2 = asta M1-M3
    A2 = line([DATI(kk,1), DATI(kk,5)],[DATI(kk,2), DATI(kk,6)]);
    set(A2, 'color', 'r', 'linewidth', 0.5);
    
    % A3 = asta M2-M4
    A3 = line([DATI(kk,3), DATI(kk,7)],[DATI(kk,4), DATI(kk,8)]);
    set(A3, 'color', 'r', 'linewidth', 0.5);
    
    % A4 = asta M3-M4
    A4 = line([DATI(kk,5), DATI(kk,7)],[DATI(kk,6), DATI(kk,8)]);
    set(A4, 'color', 'r', 'linewidth', 0.5);
    
    % A5 = asta M4-M5
    A5 = line([DATI(kk,7), DATI(kk,9)],[DATI(kk,8), DATI(kk,10)]);
    set(A5, 'color', 'r', 'linewidth', 0.5);

    % Calcolo angolo M2M4 (braccio, baricentro G1)
    angM2M4(kk) = atan(((DATI(kk,4))-(DATI(kk,8)))/((DATI(kk,7)-DATI(kk,3)))); %[rad]

    G1x(kk) = DATI(kk,3)+G1*cos(angM2M4(kk)); %[mm]
    G1y(kk) = DATI(kk,4)-G1*sin(angM2M4(kk)); %[mm]

    % Calcolo angolo M4M5 (avambraccio, baricentro G2)
    angM4M5(kk) = atand((DATI(kk,10)-DATI(kk,8))/(DATI(kk,9)-DATI(kk,7))); %[°]

    G2x(kk) = DATI(kk,9)-G2*cosd(angM4M5(kk)); %[mm]
    G2y(kk) = DATI(kk,10)-G2*sind(angM4M5(kk)); %[mm]
    
    % Plot di G1 e G2 nel tempo
    G1_time = viscircles([G1x(kk), G1y(kk)], 5, 'color', 'k');
    G2_time = viscircles([G2x(kk), G2y(kk)], 5, 'color', 'k');
    
    G1anchor = text(G1x(kk), G1y(kk)+35, 'G1');
    G2anchor = text(G2x(kk)+15, G2y(kk)-25, 'G2');

    % Plot dei marker M1 M2 M3 M4 M5
    M1_time = viscircles ([DATI(kk,1), DATI(kk,2)], 5, 'color', 'r');
    M2_time = viscircles ([DATI(kk,3), DATI(kk,4)], 5, 'color', 'r');
    M3_time = viscircles ([DATI(kk,5), DATI(kk,6)], 5, 'color', 'r');
    M4_time = viscircles ([DATI(kk,7), DATI(kk,8)], 5, 'color', 'r');
    M5_time = viscircles ([DATI(kk,9), DATI(kk,10)], 5, 'color', 'r');

    M1anchor = text(DATI(kk,1)+20, DATI(kk,2)+20, 'M1', 'Color','r');
    M2anchor = text(DATI(kk,3)-50, DATI(kk,4)-20, 'M2', 'Color','r');
    M3anchor = text(DATI(kk,5), DATI(kk,6)+30, 'M3', 'Color','r');
    M4anchor = text(DATI(kk,7), DATI(kk,8)-30, 'M4', 'Color','r');
    M5anchor = text(DATI(kk,9)+10, DATI(kk,10)-20, 'M5', 'Color','r');

    axis equal
    axis off
    set(gca,'xticklabel',[])
    set(gca,'yticklabel',[])
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    frame=getframe(video_handle);
    writeVideo(writeObj,frame);
    pause(0.1)
end

close(writeObj)
close(video_handle)

% Calcolo la posizione di G2 rispetto ad M5 (coordinate di M5 sono origine
% di un sistema di riferimento) e di G1 rispetto ad M2 (coordinate di M3
% sono origine di un sistema di riferimento)
for kk = 1:length(DATI)
    G2x_RM5(kk) = -G2*cosd(angM4M5(kk));
    G2y_RM5(kk) = -G2*sind(angM4M5(kk));
    G1x_RM2(kk) = +G1*cos(angM2M4(kk)); %[mm]
    G1y_RM2(kk) = -G1*sin(angM2M4(kk)); %[mm]
end

%% STEP 3.5 : CALCOLARE IL ROM DELLA SPALLA

ROM_spalla = max(rad2deg(angM2M4)) - min(rad2deg(angM2M4)); %[°]
fprintf('Il ROM della spalla è %.3f gradi \n\n\n', ROM_spalla);

%% STEP 3.6 : VELOCITÀ ED ACCELERAZIONE DEL PUNTO M4

% Per scrivere le velocità nella loro unità di misura corretta, a partire
% da POS_M4_norm converto tutto in metri [m]

for ii = 1:length(DATI)
    POS_M4_norm(ii) = sqrt((DATI(ii,7)/1000)^2+(DATI(ii,8)/1000)^2);
end

% Utlizzo la definizione di velocità ed accelerazione come derivata prima e
% seconda dello spostamento, dunque approssimo le derivate

for ii = 1:(length(DATI)-1)
    vel_M4S(ii) = sqrt((((DATI(ii+1,7)/1000)-(DATI(ii,7)/1000))^2)+(((DATI(ii+1,8)/1000)-(DATI(ii,8)/1000))^2))/dt;
    vel_M4x(ii) = ((DATI(ii+1,7)/1000)-(DATI(ii,7)/1000))/dt;
    vel_M4y(ii) = ((DATI(ii+1,8)/1000)-(DATI(ii,8)/1000))/dt;
end
for ii = 1:(length(vel_M4S)-1)
    acc_M4S(ii) = (vel_M4S(ii+1)-vel_M4S(ii))/dt;
    acc_M4x(ii) = (vel_M4x(ii+1)-vel_M4x(ii))/dt;
    acc_M4y(ii) = (vel_M4y(ii+1)-vel_M4y(ii))/dt;
end

%% STEP 3.7 : VELOCITÀ ED ACCELERAZIONE ANGOLARE DI M3 ED M5
% Utilizzo i dati ricavati in precedenza sul punto M4
% Condizioni iniziali, devo determinare alfa, gamma e delta

% Calcolo alfa e gamma

for jj = 1:length(DATI)
    angM1M3(jj) = atan(((DATI(jj,2)/1000)-(DATI(jj,6)/1000))/((DATI(jj,5)/1000)-(DATI(jj,1)/1000)));
    alfa(jj) = 2*pi-angM1M3(jj); % [rad]
    gamma(jj) = 2*pi-angM2M4(jj); % [rad]
    angM3M4(jj) = atan(((DATI(jj,6)/1000)-(DATI(jj,8)/1000))/((DATI(jj,5)/1000)-(DATI(jj,7)/1000)));
    delta(jj) = angM3M4(jj);
end

% Ricavo le velocità angolari di M3 ed M5
if (strcmp(scelta,'R') || strcmp(scelta,'r'))
    % uso la prima componene di ciascun vettore poiché le misurazioni lungo
    % l'asse x della posizione subiscono variazioni meno repentine rispetto a
    % quelle lungo l'asse y, duque si ottengono misure più uniformi e
    % verosimili

    for ii = 1:(length(DATI)-1)
        vel_M4 = [vel_M4x(ii);vel_M4y(ii);0];
        M2M4vec = [((DATI(ii,7)/1000)-(DATI(ii,3)/1000));((DATI(ii,4)/1000)-(DATI(ii,8)/1000));0];
        gammap(ii) = vel_M4(2)/(M2M4vec(1));
        
        deltap(ii) = (delta(ii+1)-delta(ii))/dt;

        deltapvec = [0;0;deltap(ii)];
        M3M4vec = [((DATI(ii,5)/1000)-(DATI(ii,7)/1000));((DATI(ii,6)/1000)-(DATI(ii,8)/1000));0];
        vel_M3 = vel_M4 + cross(deltapvec,M3M4vec);

        M1M3vec = [((DATI(ii,5)/1000)-(DATI(ii,1)/1000));((DATI(ii,2)/1000)-(DATI(ii,6)/1000));0];
        alfap(ii) = vel_M3(2)/(M1M3vec(1));
    end
elseif (strcmp(scelta,'C') || strcmp(scelta,'c'))
    for ii = 1:(length(DATI)-1)
        % Utilizzo una prima chiusura cinematica per ricavare gammap,
        % formata da M2M4 e la proiezione di M4 sull'asse di M2
        % Ricavo dalla proiezione lungo y poiché è la direzione
        % preponderante dello spostamento del punto M4
        gammap(ii) = (-vel_M4y(ii)*sin((3/2)*pi)+vel_M4x(ii)*sin(0))/((M2M4/1000)*sin(gamma(ii)+pi/2));

        % Utilizzo ora il quadrilatero formato da M1M2M3M4 e calcolo
        % deltap ed alfap
        COEFF = [(M1M3/1000)*cos(alfa(ii)), -(M3M4/1000)*cos(delta(ii));
                 (M1M3/1000)*sin(alfa(ii)), -(M3M4/1000)*sin(delta(ii))];
        NOTO_W = [(M2M4/1000)*gammap(ii)*cos(gamma(ii));
                  (M2M4/1000)*gammap(ii)*sin(gamma(ii));];
        sol1 = linsolve(COEFF,NOTO_W);
        alfap(ii) = sol1(1);
        deltap(ii) = sol1(2);
    end
end

% Ricavo le accelerazioni angolari di M3 ed M5
if (strcmp(scelta,'R') || strcmp(scelta,'r'))
    % uso la prima componene di ciascun vettore poiché le misurazioni lungo
    % l'asse x della posizione subiscono variazioni meno repentine rispetto a
    % quelle lungo l'asse y, duque si ottengono misure più uniformi e
    % verosimili
    for ii = 1:(length(DATI)-2)
        acc_M4 = [acc_M4x(ii);acc_M4y(ii);0];
        M2M4vec = [((DATI(ii,7)/1000)-(DATI(ii,3)/1000));((DATI(ii,4)/1000)-(DATI(ii,8)/1000));0];
        gammapp(ii) = (((gammap(ii)^2)*(M2M4vec(2)))+acc_M4(2))/(M2M4vec(1));
        
        deltapp(ii) = (deltap(ii+1)-deltap(ii))/dt;

        deltappvec = [0;0;deltapp(ii)];
        M3M4vec = [((DATI(ii,5)/1000)-(DATI(ii,7)/1000));((DATI(ii,6)/1000)-(DATI(ii,8)/1000));0];
        acc_M3 = acc_M4 + cross(deltappvec,M3M4vec);

        M1M3vec = [((DATI(ii,5)/1000)-(DATI(ii,1)/1000));((DATI(ii,2)/1000)-(DATI(ii,6)/1000));0];
        alfapp(ii) = (((alfap(ii)^2)*(M1M3vec(2)))+acc_M3(2))/(M1M3vec(1));
    end
elseif (strcmp(scelta,'C') || strcmp(scelta,'c'))
    for ii = 1:(length(DATI)-2)
        % Seguo il procedimento usato per calcolare le velocità angolari
        % utilizzando le stesse due chiusure

        % Ricavo dalla proiezione lungo y poiché è la direzione
        % preponderante dello spostamento del punto M4
        gammapp(ii) = (-acc_M4y(ii)*sin((3/2)*pi)+acc_M4x(ii)*sin(0)+(M2M4/1000)*((gammap(ii))^2)*sin(gamma(ii)+pi));

        COEFF = [(M1M3/1000)*cos(alfa(ii)), -(M3M4/1000)*cos(delta(ii));
                 (M1M3/1000)*sin(alfa(ii)), -(M3M4/1000)*sin(delta(ii))];
        NOTO_WP = [(M3M4/1000)*((deltap(ii))^2)*cos(delta(ii)+pi/2)-(M1M3/1000)*((alfap(ii))^2)*cos(alfa(ii)+pi/2)+(M2M4/1000)*gammapp(ii)*cos(gamma(ii))+(M2M4/1000)*((gammap(ii))^2)*cos(gamma(ii)+pi/2);
                   (M3M4/1000)*((deltap(ii))^2)*sin(delta(ii)+pi/2)-(M1M3/1000)*((alfap(ii))^2)*sin(alfa(ii)+pi/2)+(M2M4/1000)*gammapp(ii)*sin(gamma(ii))+(M2M4/1000)*((gammap(ii))^2)*sin(gamma(ii)+pi/2);];
        sol2 = linsolve(COEFF,NOTO_WP);
        alfapp(ii) = sol2(1);
        deltapp(ii) = sol2(2);
    end
end

%% STEP 3.8 : CONFRONTO FRA VELOCITÀ ED ACCELERAZIONE DI M5 TEORICHE E SPERIMENTALI
% Calcolo le velocità ed accelerazioni teoriche di M5 da gammap e gammapp

% V M5 (Teorica) - Rivals
for ii = 1:length(gammap)
    vel_M4 = [vel_M4x(ii);vel_M4y(ii);0]; %[m/s]
    gammapvec = [0;0;gammap(ii)]; %[rad/s]
    M4M5vec = [((DATI(ii,9)/1000)-(DATI(ii,7)/1000));((DATI(ii,10)/1000)-(DATI(ii,8)/1000));0]; %[m]

    vel_M5T_vec = vel_M4 + cross(gammapvec,M4M5vec); %[m/s]
    vel_M5Tx(ii) = vel_M5T_vec(1); 
    vel_M5Ty(ii) = vel_M5T_vec(2);
    
    if vel_M4S(ii) < 0
        vel_M5T(ii) = -sqrt((vel_M5Tx(ii))^2+(vel_M5Ty(ii))^2);
    else
        vel_M5T(ii) = sqrt((vel_M5Tx(ii))^2+(vel_M5Ty(ii))^2);
    end
end

% A M5 (Teorica) - Rivals
for ii = 1:length(gammapp)
    acc_M4 = [acc_M4x(ii);acc_M4y(ii);0]; %[m/s^2]
    gammappvec = [0;0;gammapp(ii)]; %[rad/s^2]
    M4M5vec = [((DATI(ii,9)/1000)-(DATI(ii,7)/1000));((DATI(ii,10)/1000)-(DATI(ii,8)/1000));0]; %[m]
    
    acc_M5T_vec = acc_M4 + cross(gammappvec,M4M5vec)-(gammap(ii))^2.*(M4M5vec); %[m/s^2]
    acc_M5Tx(ii) = acc_M5T_vec(1);
    acc_M5Ty(ii) = acc_M5T_vec(2);
    
    if acc_M4S(ii) < 0
        acc_M5T(ii) = -sqrt((acc_M5Tx(ii))^2+(acc_M5Ty(ii))^2);
    else
        acc_M5T(ii) = sqrt((acc_M5Tx(ii))^2+(acc_M5Ty(ii))^2);
    end
end

% V M5 (Sperimentale) - Differenze finite
for jj = 1:(length(DATI)-1)
    vel_M5S(jj) = sqrt((((DATI(jj+1,9)/1000)-(DATI(jj,9)/1000))^2)+(((DATI(jj+1,10)/1000)-(DATI(jj,10)/1000))^2))/dt;
    vel_M5Sx(jj) = ((DATI(jj+1,9)/1000)-(DATI(jj,9)/1000))/dt;
    vel_M5Sy(jj) = ((DATI(jj+1,10)/1000)-(DATI(jj,10)/1000))/dt;
end

% A M5 (Sperimentale) - Differenze finite
for jj = 1:(length(vel_M5S)-1)
    acc_M5S(jj) = (vel_M5S(jj+1)-vel_M5S(jj))/dt;
    acc_M5Sx(jj) = (vel_M5Sx(jj+1)-vel_M5Sx(jj))/dt;
    acc_M5Sy(jj) = (vel_M5Sy(jj+1)-vel_M5Sy(jj))/dt;
end

%% 4. DINAMICA
% Identificare atto di moto in cui M2M4 è orizzontale (APPROSSIMATO) 
% e scelgo atto di moto fra quelli calcolati

% Seleziono gli atti di moto per il quale M2M4 è orizzontale (APPROSSIMATO)
atti_di_moto = find(abs(angM2M4)<deg2rad(0.5));

% Scelgo un atto di moto fra quelli selezionati
fprintf('Gli atti di moto in cui \nl''asta M2M4 è orrizontale sono: \n\n');
disp(atti_di_moto);
atto = input('Per proseguire indica la posizione \ndi quello che ti interessa (1,2,ecc...): ');
fprintf('\n\n\n');
atto = atti_di_moto(atto);

% calcolo le variabili che descrivono il moto di M4 nell'atto scelto
vel_M4vec_DIN = [vel_M4x(atto); vel_M4y(atto);0]; %[m/s]
acc_M4vec_DIN = [acc_M4x(atto); acc_M4y(atto);0]; %[m/s]

% Calcolo le variabili descrittive del sistema in corrispondenza dell'atto
% di moto individuato
alfa_DIN = alfa(atto); %[rad]
gamma_DIN = gamma(atto); %[rad]
delta_DIN = delta(atto); %[rad]

M2M4vec_DIN = [((DATI(atto,7)/1000)-(DATI(atto,3)/1000));((DATI(atto,4)/1000)-(DATI(atto,8)/1000));0]; %[m]

gammap_DIN = vel_M4(2)/(M2M4vec_DIN(1)); %[rad/s]
gammapvec_DIN = [0;0;gammap_DIN];

deltapvec_DIN = [0;0;deltap(atto)]; %[rad/s]
M3M4vec_DIN = [((DATI(atto,5)/1000)-(DATI(atto,7)/1000));((DATI(atto,6)/1000)-(DATI(atto,8)/1000));0]; %[m]
vel_M3vec_DIN = vel_M4vec_DIN + cross(deltapvec_DIN,M3M4vec_DIN); %[m/s]

M1M3vec_DIN = [((DATI(atto,5)/1000)-(DATI(atto,1)/1000));((DATI(atto,2)/1000)-(DATI(atto,6)/1000));0]; %[m]

alfap_DIN = vel_M3vec_DIN(2)/(M1M3vec_DIN(1)); %[rad/s]
alfapvec_DIN = [0;0;alfap_DIN];

gammapp_DIN = (((gammap_DIN^2)*(M2M4vec_DIN(2)))+acc_M4vec_DIN(2))/(M2M4vec_DIN(1));
gammappvec_DIN = [0;0;gammapp_DIN];

deltappvec_DIN = [0;0;deltapp(atto)];
acc_M3vec_DIN = acc_M4vec_DIN + cross(deltappvec_DIN,M3M4vec_DIN);

alfapp_DIN = (((alfap_DIN^2)*(M1M3vec_DIN(2)))+acc_M3vec_DIN(2))/(M1M3vec_DIN(1));
alfappvec_DIN = [0;0;alfapp_DIN];

angM4M5 = mean(angM4M5);
M4M5_vec = [M4M5/1000*cos(deg2rad(angM4M5)); M4M5/1000*sin(deg2rad(angM4M5)); 0];
vel_M5vec_DIN = vel_M4vec_DIN + cross(gammapvec_DIN,M4M5_vec);
acc_M5vec_DIN = acc_M4vec_DIN + cross(gammappvec_DIN,M4M5_vec)-(gammapvec_DIN(3)^2).*M4M5_vec;

if gammappvec_DIN(3) < 0
    fprintf(['Durante l''atto di moto considerato (atto %d; t = %.2f s)\n' ...
        'si sta compiendo un movimento di ESTENSIONE (rotazione ORARIA)\n\n\n'], atto, atto*dt);
else
    fprintf(['Durante l''atto di moto considerato (atto %d; t = %.2f s)\n' ...
        'si sta compiendo un movimento di FLESSIONE (rotazione ANTIORARIA)\n\n\n'], atto, atto*dt);
end

%% STEP 4.1-4.2 : DETERMINARE E RAPPRESENTARE I VETTORI VEL ED ACC DI G1,G2

angM2G1 = mean(angM2M4); %[rad]
angM4G2 = deg2rad(angM4M5); %[rad]

M2G1_1 = (G1x(atto)/1000)-(DATI(atto,3)/1000); %[m]
M2G1_2 = (DATI(atto,4)/1000)-(G1y(atto)/1000); %[m]
M4G2_1 = (G2x(atto)/1000)-(DATI(atto,7)/1000); %[m]
M4G2_2 = (G2y(atto)/1000)-(DATI(atto,8)/1000); %[m]

M2G1 = [M2G1_1; M2G1_2;0];
M4G2 = [M4G2_1; M4G2_2;0];

% Rivals velocità
vel_G1_DIN = cross(gammapvec_DIN,M2G1); %[m/s]
vel_G2_DIN = vel_M4vec_DIN + cross(gammapvec_DIN,M4G2); %[m/s]

% Rivals accelerazioni
acc_G1_DIN = cross(gammappvec_DIN,M2G1)-(gammap_DIN)^2.*(M2G1); %[m/s^2]
acc_G2_DIN = acc_M4vec_DIN + cross(gammappvec_DIN,M4G2)-(gammap_DIN)^2.*(M4G2); %[m/s^2]

%% ALTRI DATI DA TABELLE ANTROPOMETRICHE (MASSE E MOMENTI D'INERZIA)

% Posizione dei baricentri di braccio ed avambraccio+mano in metri
G1 = 0.436*L_BRACCIO/1000; %[m] (rispetto alla spalla)
G2 = 0.430*L_AVAMBRACCIO/1000; %[m] (rispetto al gomito)

% Massa del soggetto utilizzato per acquisire i dati
m = 60 + CODICEPERSONA; %[kg]

% Massa di avambraccio e braccio
m_BRACCIO = 0.028*m; %[kg]
m_AVAMBRACCIO = 0.016*m; %[kg]
m_MANO = 0.006*m; %[kg]

% Momenti d'inerzia di braccio e avambraccio rispetto ai baricentri
JG1_BRACCIO = m_BRACCIO*(0.322*L_BRACCIO/1000)^2; %[kgm^2]
JG2_AVAMBRACCIO = m_AVAMBRACCIO*(0.303*L_AVAMBRACCIO/1000)^2; %[kgm^2]

%% STEP 4.3 : CALCOLO DI COPPIE DI ESTENSIONE E FLESSIONE (+SUPPORTO ANTIGRAVITARIO)

g = [0;-9.81;0]; %[m/s^2]

% Bilancio di potenze
% Considero positivo il verso antiorario
Ctot_ANTIG = (m_BRACCIO*(dot(vel_G1_DIN,acc_G1_DIN))+ ...
              m_AVAMBRACCIO*(dot(vel_G2_DIN,acc_G2_DIN))+ ...
              m_MANO*(dot(vel_M5vec_DIN,acc_M5vec_DIN)-dot(vel_M5vec_DIN,g))+ ...
              JG1_BRACCIO*dot(gammapvec_DIN,gammappvec_DIN)+JG2_AVAMBRACCIO*dot(gammapvec_DIN,gammappvec_DIN))/gammapvec_DIN(3);

% Individuo i valori della coppia di estensione e flessione, a seconda che
% il movimento sia verso il basso o verso l'alto
if Ctot_ANTIG > 0
    if gammappvec_DIN(3) > 0
       Cf_ANTIG = Ctot_ANTIG; 
       Ce_ANTIG = -0.1*Cf_ANTIG;
       
       fprintf('CON SUPPORTO ANTIGRAVITARIO - I muscoli del braccio sviluppano: \n');
       fprintf(['Per permettere il movimento (muscoli FLESSORI): %.3f Nm\n' ...
            'Per stabilizzare l''articolazione (muscoli ESTENSORI): %.3f Nm \n\n\n'], Cf_ANTIG, Ce_ANTIG);
    elseif gammappvec_DIN(3) < 0
       Ce_ANTIG = -Ctot_ANTIG; 
       Cf_ANTIG = -0.1*Ce_ANTIG;
        
       fprintf('CON SUPPORTO ANTIGRAVITARIO - I muscoli del braccio sviluppano: \n');
       fprintf(['Per permettere il movimento (muscoli ESTENSORI): %.3f Nm\n' ...
            'Per stabilizzare l''articolazione (muscoli FLESSORI): %.3f Nm \n\n\n'], Ce_ANTIG, Cf_ANTIG);
    end
elseif Ctot_ANTIG < 0
    if gammappvec_DIN(3) < 0
       Cf_ANTIG = Ctot_ANTIG; 
       Ce_ANTIG = -0.1*Cf_ANTIG;
        
       fprintf('CON SUPPORTO ANTIGRAVITARIO - I muscoli del braccio sviluppano: \n');
       fprintf(['Per permettere il movimento (muscoli FLESSORI): %.3f Nm\n' ...
            'Per stabilizzare l''articolazione (muscoli ESTENSORI): %.3f Nm \n\n\n'], Cf_ANTIG, Ce_ANTIG);
    elseif gammappvec_DIN(3) > 0
       Ce_ANTIG = -Ctot_ANTIG; 
       Cf_ANTIG = -0.1*Ce_ANTIG;
        
       fprintf('CON SUPPORTO ANTIGRAVITARIO - I muscoli del braccio sviluppano: \n');
       fprintf(['Per permettere il movimento (muscoli ESTENSORI): %.3f Nm\n' ...
            'Per stabilizzare l''articolazione (muscoli FLESSORI): %.3f Nm \n\n\n'], Ce_ANTIG, Cf_ANTIG);
    end
end

%% STEP 4.4 : CALCOLO DI COPPIE DI ESTENSIONE E FLESSIONE (NO SUPPORTO ANTIGRAVITARIO)

% Seguo il procedimento dello step precedente; la forza peso di braccio ed 
% avambraccio non viene compensata quindi concorre nell'equazione
Ctot_NOSUPP = (m_BRACCIO*(dot(vel_G1_DIN,acc_G1_DIN)-dot(vel_G1_DIN,g))+ ...
               m_AVAMBRACCIO*(dot(vel_G2_DIN,acc_G2_DIN)-dot(vel_G2_DIN,g))+ ...
               m_MANO*(dot(vel_M5vec_DIN,acc_M5vec_DIN)-dot(vel_M5vec_DIN,g))+ ...
               JG1_BRACCIO*dot(gammapvec_DIN,gammappvec_DIN)+JG2_AVAMBRACCIO*dot(gammapvec_DIN,gammappvec_DIN))/gammapvec_DIN(3);

% Individuo i valori della coppia di estensione e flessione, a seconda che
% il movimento sia verso il basso o verso l'alto
if Ctot_NOSUPP > 0
    if gammappvec_DIN(3) > 0
       Cf_NOSUPP = Ctot_NOSUPP; 
       Ce_NOSUPP = -0.1*Cf_NOSUPP;
       
       fprintf('SENZA SUPPORTO ANTIGRAVITARIO - I muscoli del braccio sviluppano: \n');
       fprintf(['Per permettere il movimento (muscoli FLESSORI): %.3f Nm\n' ...
            'Per stabilizzare l''articolazione (muscoli ESTENSORI): %.3f Nm \n\n\n'], Cf_NOSUPP, Ce_NOSUPP);
    elseif gammappvec_DIN(3) < 0
       Ce_NOSUPP = -Ctot_NOSUPP; 
       Cf_NOSUPP = -0.1*Ce_NOSUPP;
        
       fprintf('SENZA SUPPORTO ANTIGRAVITARIO - I muscoli del braccio sviluppano: \n');
       fprintf(['Per permettere il movimento (muscoli ESTENSORI): %.3f Nm\n' ...
            'Per stabilizzare l''articolazione (muscoli FLESSORI): %.3f Nm \n\n\n'], Ce_NOSUPP, Cf_NOSUPP);
    end
elseif Ctot_NOSUPP < 0
    if gammappvec_DIN(3) < 0
       Cf_NOSUPP = Ctot_NOSUPP; 
       Ce_NOSUPP = -0.1*Cf_NOSUPP;
        
       fprintf('SENZA SUPPORTO ANTIGRAVITARIO - I muscoli del braccio sviluppano: \n');
       fprintf(['Per permettere il movimento (muscoli FLESSORI): %.3f Nm\n' ...
            'Per stabilizzare l''articolazione (muscoli ESTENSORI): %.3f Nm \n\n\n'], Cf_NOSUPP, Ce_NOSUPP);
    elseif gammappvec_DIN(3) > 0
       Ce_NOSUPP = -Ctot_NOSUPP; 
       Cf_NOSUPP = -0.1*Ce_NOSUPP;
        
       fprintf('CON SUPPORTO ANTIGRAVITARIO - I muscoli del braccio sviluppano: \n');
       fprintf(['Per permettere il movimento (muscoli ESTENSORI): %.3f Nm\n' ...
            'Per stabilizzare l''articolazione (muscoli FLESSORI): %.3f Nm \n\n\n'], Ce_NOSUPP, Cf_NOSUPP);
    end
end

% Raggruppo le coppie di estensione e flessione per 
% confrontarle graficamente
Ce_CONFRONTO = [Ce_NOSUPP, Ce_ANTIG];
Cf_CONFRONTO = [Cf_NOSUPP, Cf_ANTIG];
C_ANTIG_CONFRONTO = [Cf_ANTIG, Ce_ANTIG];
C_NOSUPP_CONFRONTO = [Cf_NOSUPP, Ce_NOSUPP];

%% STEP 4.5 : CALCOLO DELLE REAZIONI VINCOLARI SULL'AGGANCIO AL TELAIO
L_CONG = M3M4/1000;
L_BRACCIO = L_BRACCIO/1000;
L_AVAMBRACCIO = L_AVAMBRACCIO/1000;

angM1M3 = mean(angM1M3);
angM3M4 = mean(angM3M4);

ang_ACC_G1= atan(acc_G1_DIN(2)/acc_G1_DIN(1)); % [rad] angolo del vettore accelerazione di G1
ang_ACC_G2= atan(acc_G2_DIN(2)/acc_G2_DIN(1)); % [rad] angolo del vettore accelerazione di G2

% Calcolo e nuove masse con +20% di tolleranza per l'esoscheletro
m_BRACCIO = 1.2*m_BRACCIO; %[kg]
m_AVAMBRACCIO = 1.2*m_AVAMBRACCIO; %[kg]
JG1_BRACCIO = m_BRACCIO*(0.322*L_BRACCIO)^2; %[kgm^2]
JG2_AVAMBRACCIO = m_AVAMBRACCIO*(0.303*L_AVAMBRACCIO)^2; %[kgm^2]
g = 9.81;

% Risolvo ora il sistema ricavando le seguenti equazioni

K = [0 0 0 1 0 0 0 -1;
     0 0 -1 0 0 0 1 0;
     0 0 0 0 0 0 L_BRACCIO*cos(angM2G1) -L_BRACCIO*sin(angM2G1);
     0 0 0 0 0 -1 0 1;
     0 0 0 0 1 0 -1 0;
     0 0 0 0 0 0 L_CONG*cos(angM3M4) L_CONG*sin(angM3M4);
     0 -1 0 0 0 1 0 0;
     1 0 0 0 -1 0 0 0;];

F = [m_BRACCIO*acc_G1_DIN(1)+m_AVAMBRACCIO*acc_G2_DIN(1)+m_MANO*acc_M5vec_DIN(1);
     m_BRACCIO*g+m_BRACCIO*acc_G1_DIN(2)+m_AVAMBRACCIO*g+m_AVAMBRACCIO*acc_G2_DIN(2)+m_MANO*g+m_MANO*acc_M5vec_DIN(2);
     (JG1_BRACCIO+JG2_AVAMBRACCIO)*gammapp_DIN+m_BRACCIO*acc_G1_DIN(2)*G1*cos(angM2G1)+m_BRACCIO*g*G1*cos(angM2G1)+m_BRACCIO*acc_G1_DIN(1)*G1*sin(angM2G1)+m_AVAMBRACCIO*acc_G2_DIN(1)*(L_BRACCIO*sin(angM2G1)-G2*sin(angM4G2))+m_AVAMBRACCIO*g*(L_BRACCIO*cos(angM2G1)+G2*cos(angM4G2))+m_AVAMBRACCIO*acc_G2_DIN(2)*(L_BRACCIO*cos(angM2G1)+G2*cos(angM4G2))+m_MANO*acc_M5vec_DIN(1)*(L_BRACCIO*sin(angM2G1)-L_AVAMBRACCIO*sin(angM4G2))+m_MANO*g*(L_BRACCIO*cos(angM2G1)+L_AVAMBRACCIO*cos(angM4G2))+m_MANO*acc_M5vec_DIN(2)*(L_BRACCIO*cos(angM2G1)+L_AVAMBRACCIO*cos(angM4G2));
     0;
     0;
     0;
     0;
     0;];

RV = linsolve(K,F);
V1_REAL = RV(1);
H1_REAL = RV(2);
V2_REAL = RV(3);
H2_REAL = RV(4);

fprintf(['Le reazioni vincolari calcolate sul telaio valgono:\n (Per un soggetto di 167cm X 6%dkg) \nIN M1:\n ' ...
    'H_M1 = %.3f N verso sinistra \n V_M1 = %.3f N verso l''alto\n' ...
    'IN M2:\n H_M2 = %.3f N verso destra \n V_M2 = %.3f N verso il basso\n\n\n'], CODICEPERSONA, H1_REAL, V1_REAL, H2_REAL, V2_REAL);

%% STEP 4.6 : REAZIONI VINCOLARI CALCOLATE PER UN PAZIENTE DI 90KG
L_BRACCIO_FAT = 0.186*1.9;
L_AVAMBRACCIO_FAT = 0.146*1.9;

m_FAT = 90; %[kg]
m_BRACCIO_FAT = 1.2*0.028*m_FAT; %[kg]
m_AVAMBRACCIO_FAT = 1.2*0.016*m_FAT; %[kg]
m_MANO_FAT = 0.006*m_FAT; %[kg]
G1_FAT = 0.436*L_BRACCIO_FAT; %[m] (rispetto alla spalla)
G2_FAT = 0.430*L_AVAMBRACCIO_FAT; %[m] (rispetto al gomito)

JG1_BRACCIO_FAT = m_BRACCIO_FAT*(0.322*L_BRACCIO_FAT)^2; %[kgm^2]
JG2_AVAMBRACCIO_FAT = m_AVAMBRACCIO_FAT*(0.303*L_AVAMBRACCIO_FAT)^2; %[kgm^2]

% Risolvo ora il sistema ricavando le seguenti equazioni

K_FAT = [0 0 0 1 0 0 0 -1;
     0 0 -1 0 0 0 1 0;
     0 0 0 0 0 0 L_BRACCIO_FAT*cos(angM2G1) -L_BRACCIO_FAT*sin(angM2G1);
     0 0 0 0 0 -1 0 1;
     0 0 0 0 1 0 -1 0;
     0 0 0 0 0 0 L_CONG*cos(angM3M4) L_CONG*sin(angM3M4);
     0 -1 0 0 0 1 0 0;
     1 0 0 0 -1 0 0 0;];

F_FAT = [m_BRACCIO_FAT*acc_G1_DIN(1)+m_AVAMBRACCIO_FAT*acc_G2_DIN(1)+m_MANO_FAT*acc_M5vec_DIN(1);
     m_BRACCIO_FAT*g+m_BRACCIO_FAT*acc_G1_DIN(2)+m_AVAMBRACCIO_FAT*g+m_AVAMBRACCIO_FAT*acc_G2_DIN(2)+m_MANO_FAT*g-m_MANO_FAT*acc_M5vec_DIN(2);
     (JG1_BRACCIO_FAT+JG2_AVAMBRACCIO_FAT)*gammapp_DIN+m_BRACCIO_FAT*acc_G1_DIN(2)*G1_FAT*cos(angM2G1)+m_BRACCIO_FAT*g*G1_FAT*cos(angM2G1)+m_BRACCIO_FAT*acc_G1_DIN(1)*G1_FAT*sin(angM2G1)+m_AVAMBRACCIO_FAT*acc_G2_DIN(1)*(L_BRACCIO_FAT*sin(angM2G1)-G2_FAT*sin(angM4G2))+m_AVAMBRACCIO_FAT*g*(L_BRACCIO_FAT*cos(angM2G1)+G2_FAT*cos(angM4G2))+m_AVAMBRACCIO_FAT*acc_G2_DIN(2)*(L_BRACCIO_FAT*cos(angM2G1)+G2_FAT*cos(angM4G2))+m_MANO_FAT*acc_M5vec_DIN(1)*(L_BRACCIO_FAT*sin(angM2G1)-L_AVAMBRACCIO_FAT*sin(angM4G2))+m_MANO_FAT*g*(L_BRACCIO_FAT*cos(angM2G1)+L_AVAMBRACCIO_FAT*cos(angM4G2))+m_MANO_FAT*acc_M5vec_DIN(2)*(L_BRACCIO_FAT*cos(angM2G1)+L_AVAMBRACCIO_FAT*cos(angM4G2));
     0;
     0;
     0;
     0;
     0;];

RV_FAT = linsolve(K_FAT,F_FAT);
V1_FAT = RV_FAT(1);
H1_FAT = RV_FAT(2);
V2_FAT = RV_FAT(3);
H2_FAT = RV_FAT(4);

fprintf(['Le reazioni vincolari calcolate sul telaio valgono:\n (Per un soggetto di 190cm X 90kg) \nIN M1:\n ' ...
    'H_M1 = %.3f N verso sinistra \n V_M1 = %.3f N verso l''alto\n ' ...
    'IN M2:\n H_M2 = %.3f N verso destra \n V_M2 = %.3f N verso il basso\n\n\n'], H1_FAT, V1_FAT, H2_FAT, V2_FAT);

%% 5. VIBRAZIONI

% Definisco i parametri di calcolo per il sistema vibrante
m_MANO; %[kg]
k = 1e7/(P+1); %[N/m]
r = 500; %[Ns/m]    0
A = 5/1000; %[m]
OO = P*1000; %[rad/s]

%% STEP 5.1 : DEFINIZIONE DI FREQUENZA PROPRIA E SMORZAMENTO ADIMENSIONALE DEL SISTEMA

w = sqrt(k/m_MANO); %[rad/s]
f = w/(2*pi); %[Hz]

rc = 2*m_MANO*w; %[Ns/m]
h = r/rc; %[-]

fprintf(['La pulsazione propria del sistema vibrante risulta: %.3f rad/s \n' ...
    'La frequenza propria del sistema vibrante risulta: %.3f Hz\n\n' ...
    'Smorzamento adimensionale: h = %.3f\n'], w, f, h);

if h>1
    fprintf('Il sistema vibrante è sovrasmorzato\n\n\n');
elseif h==1
    fprintf('Il sistema vibrante è smorzato criticamente\n\n\n');
else
    fprintf('Il sistema vibrante è sottosmorzato\n\n\n');
end

% Nella risoluzione considero x0 = A ed xp0 = 0
x0 = A;
xp0 = 0;

% Calcolo e rappresentazione della risposta libera non smorzata del sistema
tempo_sym = linspace(0,0.05,1000);
syms C1 C2 t real
xlibns_sym = real(C1*exp(1i*w*t)+C2*exp(-1i*w*t));
xplibns_sym = diff(xlibns_sym,t);
xlibns_sym0 = vpa(subs(xlibns_sym,t,0)) == x0;
xplibns_sym0 = vpa(subs(xplibns_sym,t,0)) == 0;

[C1,C2] = solve([xlibns_sym0 xplibns_sym0], [C1 C2]);
C1num = double(C1);
C2num = double(C2);

xlibns_subnum = real(C1num*exp(1i*w.*tempo_sym)+C2num*exp(-1i*w.*tempo_sym));

% Calcolo e rappresentazione della risposta libera smorzata del sistema
alpha = r/(2*m_MANO);
% si può calcolare come alpha = h*w
ws = imag((sqrt(r^2-4*m_MANO*k))/(2*m_MANO));
% si può calcolare come ws = w*sqrt(1-h^2)
syms C phi t real
xlibs_sym = exp(-alpha*t)*(C*cos(ws*t+phi));
xplibs_sym = diff(xlibs_sym,t);
xlibs_sym0 = vpa(subs(xlibs_sym,t,0)) == x0;
xplibs_sym0 = vpa(subs(xplibs_sym,t,0)) == 0;

[C,phi] = solve([xlibs_sym0 xplibs_sym0], [C phi]);
Cnum = double(C(1));
phinum = rad2deg(double(phi(1)));

xlibs_subnum = exp(-alpha*tempo_sym).*(Cnum.*cos(ws*tempo_sym+phinum));


%% STEP 5.2 : CALCOLO E RAPPRESENTAZIONE DI FDT E RISPOSTA FORZATA
syms t X0 W
Y = A*exp(1i*W*t);
YP = diff(Y,t);
XIP = X0*exp(1i*W*t);
XIP_P = diff(XIP,t);
XIP_PP = diff(XIP_P,t);

eq_moto = collect((m_MANO*XIP_PP+r*XIP_P+k*XIP),X0) == collect((k*Y+r*YP),A);
FDT = simplify((solve(eq_moto,X0))/A);

if (strcmp(key,'L') || strcmp(key,'l'))
    MOD_FDT_sym = 20*log10(norm(FDT));
else
    MOD_FDT_sym = norm(FDT);
end
    FASE_FDT_sym = angle(FDT);

OMEGAvec = 0:100:(10000-1);
for ii = 1:length(OMEGAvec)
    MOD_FDT_subnum(ii) = double(vpa(subs(MOD_FDT_sym,W,OMEGAvec(ii))));
    FASE_FDT_subnum(ii) = rad2deg(double(vpa(subs(FASE_FDT_sym,W,OMEGAvec(ii)))));
end

MOD_FDT_OO = double(vpa(subs(MOD_FDT_sym,W,OO)));
FASE_FDT_OO = rad2deg(double(vpa(subs(FASE_FDT_sym,W,OO))));

X0_mf = A*MOD_FDT_OO;
phi0_mf = deg2rad(FASE_FDT_OO);
int_part_mf = real((X0_mf)*exp(1i*phi0_mf)*exp(1i*OO.*tempo_sym));
S_completa_mf = xlibs_subnum + int_part_mf;

%% GRAFICI E FIGURE CINEMATICA

figure('Name','TRAIETTORIE DEI MARKER (rispetto agli assi x e y)','NumberTitle','off')
tiledlayout(5,2);

nexttile
plot(tempo, DATI(:,1), 'color', 'blue');
grid on
ylabel('M1 $x$ [mm]', 'Interpreter','latex');
xlim([0 tempo(end)]);
nexttile
plot(tempo, DATI(:,2), 'color', 'red');
grid on
ylabel('M1 $y$ [mm]', 'Interpreter','latex');
xlim([0 tempo(end)]);

nexttile
plot(tempo, DATI(:,3), 'color', 'blue');
grid on
ylabel('M2 $x$ [mm]', 'Interpreter','latex');
xlim([0 tempo(end)]);
nexttile
plot(tempo, DATI(:,4), 'color', 'red');
grid on
ylabel('M2 $y$ [mm]', 'Interpreter','latex');
xlim([0 tempo(end)]);

nexttile
plot(tempo, DATI(:,5), 'color', 'blue');
grid on
ylabel('M3 $x$ [mm]', 'Interpreter','latex');
xlim([0 tempo(end)]);
nexttile
plot(tempo, DATI(:,6), 'color', 'red');
grid on
ylabel('M3 $y$ [mm]', 'Interpreter','latex');
xlim([0 tempo(end)]);

nexttile
plot(tempo, DATI(:,7), 'color', 'blue');
grid on
ylabel('M4 $x$ [mm]', 'Interpreter','latex');
xlim([0 tempo(end)]);
nexttile
plot(tempo, DATI(:,8), 'color', 'red');
grid on
ylabel('M4 $y$ [mm]', 'Interpreter','latex');
xlim([0 tempo(end)]);

nexttile
plot(tempo, DATI(:,9), 'color', 'blue');
grid on
ylabel('M5 $x$ [mm]', 'Interpreter','latex');
xlabel('t [s]', 'Interpreter','latex');
xlim([0 tempo(end)]);
nexttile
plot(tempo, DATI(:,10), 'color', 'red');
grid on
ylabel('M5 $y$ [mm]', 'Interpreter','latex');
xlabel('t [s]', 'Interpreter','latex');
xlim([0 tempo(end)]);

figure('Name','MOVIMENTO DEI BARICENTRI G1 (R/spalla) E G2 (R/mano)','NumberTitle','off')
subplot(121)
plot(G1x_RM2,G1y_RM2, 'Color', 'red');
grid on
xlabel('x [mm]', 'Interpreter', 'latex');
ylabel('y [mm]', 'Interpreter', 'latex');
title('Movimento di G1 (R/ spalla)');

subplot(122)
plot(G2x_RM5,G2y_RM5, 'Color', 'blue');
grid on
xlabel('x [mm]', 'Interpreter', 'latex');
ylabel('y [mm]', 'Interpreter', 'latex');
title('Movimento di G2 (R/ mano)');

figure('Name','VELOCITÀ ED ACCELERAZIONE DEL PUNTO M4','NumberTitle','off')
subplot(1,2,1)
plot(tempo(1:(end-1)), vel_M4S, 'Color', 'r');
grid on
title('Velocità del punto M4');
ylabel('v [m/s]', 'Interpreter', 'latex');
xlabel('t [s]', 'Interpreter', 'latex');
xlim([0 tempo((end)-1)]);

subplot(1,2,2)
plot(tempo(1:(end-2)),acc_M4S, 'Color', 'b');
grid on
title('Accelerazione del punto M4');
ylabel('a [m/$\mathrm{s}^2$]', 'Interpreter', 'latex');
xlabel('t [s]', 'Interpreter', 'latex');
xlim([0 tempo((end)-2)]);

figure('Name','ANGOLI ALFA (M1M3) E GAMMA(M2M4) - rispetto al semiasse positivo delle x','NumberTitle','off')
hold on
plot(tempo, alfa, 'Color', 'red');
plot(tempo, gamma, 'Color', 'blue');
grid on
ylabel('Ampiezza [rad]', 'Interpreter','latex');
xlabel('t [s]', 'Interpreter', 'latex');
legend('$\alpha$', '$\gamma$', 'Interpreter', 'latex', 'Location', 'northwest');
xlim([0 tempo(end)]);

figure('Name','ALFAP E GAMMAP','NumberTitle','off')
hold on
plot(tempo(1:(end-1)), alfap);
plot(tempo(1:(end-1)), gammap);
grid on
xlabel('t [s]', 'Interpreter', 'latex');
legend('$\dot{\alpha}$ ($\omega_{M1M3}$)', '$\dot{\gamma}$ ($\omega_{M2M4}$)', ...
    'Interpreter', 'latex', 'Location', 'northwest')
ylabel('$\omega$ [rad/s]', 'Interpreter', 'latex');
xlim([0 tempo((end)-1)]);

figure('Name','ALFAPP E GAMMAPP','NumberTitle','off')
hold on
plot(tempo(1:(end-2)), alfapp);
plot(tempo(1:(end-2)), gammapp);
grid on
xlabel('t [s]', 'Interpreter', 'latex');
ylabel('$\dot{\omega}$ [rad/$\mathrm{s}^2$]', 'Interpreter', 'latex');
legend('$\ddot{\alpha}$ ($\dot{\omega}_{M1M3}$)', '$\ddot{\gamma}$ ($\dot{\omega}_{M2M4}$)', ...
    'Interpreter', 'latex', 'Location', 'northwest');
xlim([0 tempo((end)-2)]);

figure('Name','DELTA, DELTAP E DELTAPP','NumberTitle','off')
subplot(311)
plot(tempo(1:(end-2)), delta(1:(end-2)), 'Color', 'k');
grid on
ylabel('$\delta$ [rad]', 'Interpreter', 'latex');
xlim([0 tempo(end)]);

subplot(312)
plot(tempo(1:(end-2)), deltap(1:(end-1)), 'Color', 'r');
grid on
ylabel('$\dot{\delta}$  [rad/s]', 'Interpreter', 'latex');
xlim([0 tempo((end)-1)]);

subplot(313)
plot(tempo(1:(end-2)), deltapp, 'Color', 'b');
grid on
ylabel('$\ddot{\delta}$ [rad/$\mathrm{s}^2$]', 'Interpreter', 'latex')
xlabel('t [s]', 'Interpreter', 'latex');
xlim([0 tempo((end)-2)]);

figure('Name','CONFRONTO FRA VELOCITÀ ED ACCELERAZIONE DI M5 (lungo l''asse x)','NumberTitle','off')
subplot (121)
hold on
plot(tempo(1:(end-1)), vel_M5Tx);
plot(tempo(1:(end-1)), vel_M5Sx);
title('$v_{M5}$ teorica e sperimentale (asse x)', 'Interpreter', 'latex');
grid on
legend('Teorico', 'Sperimentale', 'Location', 'northwest');
ylabel('$v_x$ [m/s]', 'Interpreter', 'latex');
xlabel('t [s]', 'Interpreter', 'latex');
xlim([0 tempo((end)-1)]);

subplot(122)
hold on
plot(tempo(1:(end-2)), acc_M5Tx);
plot(tempo(1:(end-2)), acc_M5Sx);
title('$a_{M5}$ teorica e sperimentale (asse x)', 'Interpreter', 'latex');
grid on
legend('Teorico', 'Sperimentale', 'Location', 'northwest');
ylabel('$a_x$ [m/$\mathrm{s}^2$]', 'Interpreter', 'latex');
xlabel('t [s]', 'Interpreter', 'latex');
xlim([0 tempo((end)-2)]);

figure('Name','CONFRONTO FRA VELOCITÀ ED ACCELERAZIONE DI M5 (lungo l''asse y)','NumberTitle','off')
subplot (121)
hold on
plot(tempo(1:(end-1)), vel_M5Ty);
plot(tempo(1:(end-1)), vel_M5Sy);
title('$v_{M5}$ teorica e sperimentale (asse y)', 'Interpreter', 'latex');
grid on
legend('Teorico', 'Sperimentale', 'Location', 'northwest');
ylabel('$v_y$ [m/s]', 'Interpreter', 'latex');
xlabel('t [s]', 'Interpreter', 'latex');
xlim([0 tempo((end)-1)]);

subplot(122)
hold on
plot(tempo(1:(end-2)), acc_M5Ty);
plot(tempo(1:(end-2)), acc_M5Sy);
title('$a_{M5}$ teorica e sperimentale (asse y)', 'Interpreter', 'latex');
grid on
legend('Teorico', 'Sperimentale', 'Location', 'northwest');
ylabel('$a_y$ [m/$\mathrm{s}^2$]', 'Interpreter', 'latex');
xlabel('t [s]', 'Interpreter', 'latex');
xlim([0 tempo((end)-2)]);

figure('Name','CONFRONTO FRA VELOCITÀ ED ACCELERAZIONE DI M5 (modulo)','NumberTitle','off')
subplot(121)
hold on
plot(tempo(1:(end-1)), vel_M5T);
plot(tempo(1:(end-1)), vel_M5S);
title('$v_{M5}$ teorica e sperimentale (modulo)', 'Interpreter', 'latex');
grid on
legend('Teorico', 'Sperimentale', 'Location', 'northwest');
ylabel('v [m/s]', 'Interpreter', 'latex');
xlabel('t [s]', 'Interpreter', 'latex');
xlim([0 tempo((end)-1)]);

subplot(122)
hold on
plot(tempo(1:(end-2)), acc_M5T);
plot(tempo(1:(end-2)), acc_M5S);
title('$a_{M5}$ teorica e sperimentale (modulo)', 'Interpreter', 'latex');
grid on
legend('Teorico', 'Sperimentale', 'Location', 'northwest');
ylabel('a [m/$\mathrm{s}^2$]', 'Interpreter', 'latex');
xlabel('t [s]', 'Interpreter', 'latex');
xlim([0 tempo((end)-2)]);

%% GRAFICI E FIGURE DINAMICA

figure('Name','VETTORI VELOCITÀ DI G1 E G2','NumberTitle','off')
hold on
quiver(0,0,vel_G1_DIN(1),vel_G1_DIN(2));
quiver(0,0,vel_G2_DIN(1),vel_G2_DIN(2));
title('Vettori velocità di G1 e G2 (scala 1:1)');
legend('$v_{G1}$', '$v_{G2}$', 'Interpreter', 'latex', ...
    'Location', 'northwest');
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$y$', 'Interpreter', 'latex');

figure('Name','VETTORI ACCELERAZIONE DI G1 E G2','NumberTitle','off')
hold on
quiver(0,0,acc_G1_DIN(1),acc_G1_DIN(2));
quiver(0,0,acc_G2_DIN(1),acc_G2_DIN(2));
title('Vettori accelerazione di G1 e G2 (scala 1:1)');
legend('$a_{G1}$', '$a_{G2}$', 'Interpreter', 'latex', ...
    'Location', 'northwest');
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$y$', 'Interpreter', 'latex');

figure('Name','VALORI DI COPPIA ESERCITATI DAI MUSCOLI','NumberTitle','off')
subplot(121)
cfrF = bar(categorical({'C_f no supp', 'C_f supp'}), Cf_CONFRONTO, 0.4);
title('Coppia muscoli flessori');
xtips1F = cfrF(1).XEndPoints;
ytips1F = cfrF(1).YEndPoints;
labels1F = string(cfrF(1).YData);
cfrF.FaceColor = 'flat';
cfrF.CData(2,:) = [1 0 0];
text(xtips1F,ytips1F,labels1F,'HorizontalAlignment','center',...
    'VerticalAlignment','top')
ylabel('[$\mathrm{Nm}$]', 'Interpreter', 'latex');

subplot(122)
cfrE = bar(categorical({'C_e no supp', 'C_e supp'}), Ce_CONFRONTO, 0.4);
title('Coppia muscoli estensori');
xtips1E = cfrE(1).XEndPoints;
ytips1E = cfrE(1).YEndPoints;
labels1E = string(cfrE(1).YData);
cfrE.FaceColor = 'flat';
cfrE.CData(2,:) = [1 0 0];
text(xtips1E,ytips1E,labels1E,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
ylabel('[$\mathrm{Nm}$]', 'Interpreter', 'latex');

figure('Name','VALORI DI COPPIA ESERCITATI CON E SENZA SUPPORTO','NumberTitle','off')
subplot(121)
cfrC_ANTIG = bar(categorical({'C_f supp', 'C_e supp'}), C_ANTIG_CONFRONTO, 0.4);
title('Coppie con supporto antigravitario');
xtips1C_ANTIG = cfrC_ANTIG(1).XEndPoints;
ytips1C_ANTIG = cfrC_ANTIG(1).YEndPoints;
labels1C_ANTIG = string(cfrC_ANTIG(1).YData);
cfrC_ANTIG.FaceColor = 'flat';
cfrC_ANTIG.CData(2,:) = [1 0 0];
text(xtips1C_ANTIG,ytips1C_ANTIG,labels1C_ANTIG,'HorizontalAlignment','center',...
    'VerticalAlignment','top')
ylabel('[$\mathrm{Nm}$]', 'Interpreter', 'latex');

subplot(122)
cfrCNOSUPP = bar(categorical({'C_f no supp', 'C_e no supp'}), C_NOSUPP_CONFRONTO, 0.4);
title('Coppie senza supporto antigravitario');
xtips1CNOSUPP = cfrCNOSUPP(1).XEndPoints;
ytips1CNOSUPP = cfrCNOSUPP(1).YEndPoints;
labels1CNOSUPP = string(cfrCNOSUPP(1).YData);
cfrCNOSUPP.FaceColor = 'flat';
cfrCNOSUPP.CData(2,:) = [1 0 0];
text(xtips1CNOSUPP,ytips1CNOSUPP,labels1CNOSUPP,'HorizontalAlignment','center',...
    'VerticalAlignment','top')
ylabel('[$\mathrm{Nm}$]', 'Interpreter', 'latex');

%% GRAFICI E FIGURE VIBRAZIONI
figure('Name','RISPOSTA LIBERA NON SMORZATA DEL SISTEMA VIBRANTE','NumberTitle','off')
plot(tempo_sym,xlibns_subnum)
grid on
title('Sistema NON SMORZATO')
xlabel('t [s]','Interpreter','latex');
ylabel('x(t) [m]', 'Interpreter','latex');
xlim([0 tempo_sym(end)]);
ylim([1.1*min(xlibns_subnum) 1.1*max(xlibns_subnum)]);

figure('Name','RISPOSTA LIBERA SMORZATA DEL SISTEMA VIBRANTE','NumberTitle','off')
plot(tempo_sym,xlibs_subnum)
grid on
title('Sistema SMORZATO')
xlabel('t [s]','Interpreter','latex');
ylabel('x(t) [m]', 'Interpreter','latex');
xlim([0 tempo_sym(end)]);

figure('Name','DIAGRAMMI DEL MODULO E DELLA FASE DELLA FDT','NumberTitle','off')
if (strcmp(key,'L') || strcmp(key,'l'))
    subplot(211)
    semilogx(OMEGAvec,MOD_FDT_subnum)
    hold on
    semilogx(OO,MOD_FDT_OO,'r*')
    grid on
    xlabel('$\omega$ [rad/s]','Interpreter','latex');
    ylabel('$|F(\omega)|$ [dB]','Interpreter','latex');
    legend('','$|F(\omega = \Omega)|$','Interpreter','latex');
    xlim([0 OMEGAvec(end)]);
    subplot(212)
    semilogx(OMEGAvec,FASE_FDT_subnum)
    hold on
    semilogx(OO,FASE_FDT_OO,'r*')
    grid on
    xlabel('$\omega$ [rad/s]','Interpreter','latex');
    ylabel('$\mathrm{arg}(F(\omega)) [^{\mathrm{o}}]$','Interpreter','latex');
    xlim([0 OMEGAvec(end)]);
    ylim([1.1*min(FASE_FDT_subnum) 5]);
    legend('','$\mathrm{arg}(F(\omega = \Omega))$','Interpreter','latex');
else
    subplot(211)
    plot(OMEGAvec,MOD_FDT_subnum)
    hold on
    plot(OO,MOD_FDT_OO,'r*')
    grid on
    xlabel('$\omega$ [rad/s]','Interpreter','latex');
    ylabel('$|F(\omega)|$ [m]','Interpreter','latex');
    legend('','$|F(\omega = \Omega)|$','Interpreter','latex');
    ylim([0 1.1*max(MOD_FDT_subnum)]);
    xlim([0 OMEGAvec(end)]);
    subplot(212)
    plot(OMEGAvec,FASE_FDT_subnum)
    hold on
    plot(OO,FASE_FDT_OO,'r*')
    hold on
    grid on
    xlabel('$\omega$ [rad/s]','Interpreter','latex');
    ylabel('$\mathrm{arg}(F(\omega)) [^{\mathrm{o}}]$','Interpreter','latex');
    xlim([0 OMEGAvec(end)]);
    ylim([1.1*min(FASE_FDT_subnum) 5]);
    legend('','$\mathrm{arg}(F(\omega = \Omega))$','Interpreter','latex');
end

figure('Name','RISPOSTA FORZATA DA y(t) DEL SISTEMA VIBRANTE A REGIME','NumberTitle','off')
plot(tempo_sym,S_completa_mf)
title('RISPOSTA FORZATA DEL SISTEMA');
grid on
xlabel('t [s]','Interpreter','latex');
ylabel('x(t) [m]', 'Interpreter','latex');
xlim([0 tempo_sym(end)]);

%% STEP 5.3 : COMMENTO A POSSIBILI ALTERNATIVE DI FORZANTI CON PULSAZIONE DIVERSA
prompt = {sprintf(['La pulsazione di eccitazione Omega della forzante è direttamente correlata alla frequenza, e quindi alla precisione\n' ...
    'del feedback dato al paziente dal grip module, in quanto più essa è alta maggiore è il numero di oscillazioni nel tempo, che quindi\n' ...
    'vengono avvertite con minore ritardo e si interrompono con minore ritardo, oltre che con una maggiore probabilità statistica di trovare il sistema vibrante\n' ...
    'in una posizione prossima a quella di equilibrio x(t) = 0, che dunque richiede di non esservici riportato successivamente ed in modo artificiale.\n' ...
    'Una pulsazione di eccitazione minore richiede alla forzante un tempo maggiore per far compiere al sistema un''oscillazione completa,\n' ...
    'dunque influenza negativamente la precisione del feedback del grip module.\n' ...
    'L''ampiezza della forzante influenza invece l''intensità del feedback di vibrazione, dunque aumentando lo rende più intenso,\n' ...
    'viceversa diminuendo. Vanno imposti dei limiti all''ampiezza della forzante poiché se essa è troppo grande aumenta\n' ...
    'la complessità della gestione del sistema, in quanto si ha a che fare con grandezze caratteristiche di dimensioni maggiori che richiedono un\n' ...
    'maggiore controllo per essere gestite correttamente (ad esempio per riportare il sistema all''equilibrio)\n\n'])};
msgbox(prompt,'influenza di pulsazione ed ampiezza');
clear prompt
