TMsteluta = TM2+TTohm;
KMsteluta = K2*KTohm;

%% Proiectarea primului regulator cu Guillman-Truxal

sigma = 0.13; % Ales de mine
tita = abs(log(sigma)/sqrt(pi^2+(log(sigma)^2)))
tr = 1; % Ales de mine

wn = 4/tita/tr
Wb = wn*sqrt(1-2*tita*tita+sqrt(2-4*tita*tita+4*tita^4))
Estv = 2*tita/wn

H02 = tf(wn^2,[1 2*tita*wn wn^2])

% Pentru verificare vom face step(H02)
step(H02);
title("Intrare de tip treapta");
figure;
t = 0:0.1:30;
lsim(H02,t,t);
title("Intrare de tip rampa");


%%
Hf = tf(Kap*K1*KMsteluta,conv([TAUap 1],conv([TM1 1],[TMsteluta 1])));
Hd = tf(wn^2,[1 2*tita*wn 0]);


Hajutor = tf(wn/2/tita,[1/2/tita/wn 1 0]);
HR1 = Hajutor/Hf;

%% d.1)
tita2 = 0.51;
wn2 = 7.1;
HajutorPrim = tf(wn/2/tita,[1/2/tita2/wn2 1 0]);
HR1prim = HajutorPrim/Hf;

HR1primSimplificat = tf(0.0034139*conv([1 100],[1 1.049]),[1 0]);
H02prim = feedback(series(HR1primSimplificat,Hf),1)

%% d.2)

TMstelutasteluta = TMsteluta + TAUap;
HR1secund = tf((wn/2/tita) * conv([TM1 1],[TMstelutasteluta 1]) , (Kap*K1*KMsteluta) * [1/2/tita/wn 1 0]);

H02secund = feedback(series(HR1secund,Hf),1);
% Verificare performante dupa simplificare
step(H02secund)

%% e)
% Se scot din HR1primSimplificat
Vr = 0.345;
TAUi = 0.963;
TAUd = 0.009895;
TN = 0.1 * TAUi;

% Comparatiile
%La treapta
step(H02prim,H02secund,H02);
legend;
title("Raspunsul la treapta")
figure;
t = 0:0.1:30;
lsim(H02,t,t);
hold on;
lsim(H02prim,t,t);
hold on;
lsim(H02secund,t,t);
legend;
title("Intrare de tip rampa");

%% Proiectare al doilea regulator Guillman-Truxal cu corectie
sigmastelat = 0.1;
tita = abs(log(sigmastelat))/sqrt(pi^2+log(sigmastelat)^2);
trstelat = 0.9;
wn = 4/tita/trstelat;
Ho = tf(wn^2,[1 2*tita*wn wn^2]);
%step(Ho);
Wb = wn*sqrt(1-2*tita^2+sqrt(2-4*tita^2+4*tita^4));
% Nu respecta Estv si trecem la corectie

deltasigma = 0.01; % 5%
sigma2 = sigmastelat - deltasigma;
tita2 = abs(log(sigma2))/sqrt(pi^2+log(sigma2)^2);
wn2 = Wb/sqrt(1-2*tita2^2+sqrt(2-4*tita2^2+4*tita2^4));
Estvstelat = 0.05;
Pc = deltasigma/(2*(tita2/wn2)-Estvstelat);
Zc = Pc / (1+deltasigma);

H0c = tf(wn2^2*Pc*[1 Zc],Zc*conv([1 2*tita2*wn2 wn2^2],[1 Pc]));
step(H0c);
title('Raspunsul la treapta')
figure;
lsim(H0c,t,t)
title('Raspunsul la rampa')
HR2 = minreal(zpk(H0c/((1-H0c)*Hf)));



%% d.1) 

HR2prim = tf(0.9847*conv([0.01 1],conv([0.1386 1],conv([0.8468 1],[10.9733 1]))),conv([1 0],[34.7102 1]));
%HR2prim = tf(0.9847*conv([0.01 1],conv([0.1386 1],conv([0.9533 1],[10.8668 1]))),conv([1 0],[34.7102 1]));
H0cprim = feedback(HR2prim*Hf,1);
step(H0cprim);

%% d.2) timp de raspuns bun nu si suprareglaj
HR2secund = tf(0.9846*conv([0.1386 1],conv([0.9533 1],[10.9833 1])),conv([1 0],conv([0.1065 1],[34.7102 1])));
H0csecund = feedback(HR2secund*Hf,1);
step(H0csecund);

%% d.3) 
HR2primSecund = tf(3.6602*conv([0.01 1],[0.9533 1]),[1 0]);
H0cprimSecund = feedback(HR2primSecund*Hf,1);
step(H0cprimSecund)

%% d) Parametrii de acord a regulatorului
Vr = 0.9484;
TAUi = 0.9633;
TAUd = 0.0099;
TN = 0.1 * TAUi;

%% Raspunsurile la treapta 
step(H0c,H0cprim);
title('Raspunsul la rampa');
legend;
figure;
step(H0c,H0csecund);
title('Raspunsul la rampa');
legend;
figure;
step(H0c,H0cprimSecund);
title('Raspunsul la rampa');
legend;

%% Raspunsurile la rampa
t = 0:0.1:30;
lsim(H0c,t,t);
hold on;
lsim(H0cprim,t,t);
hold on;
lsim(H0csecund,t,t);
hold on;
lsim(H0cprimSecund,t,t);
title("Intrare de tip rampa");
legend