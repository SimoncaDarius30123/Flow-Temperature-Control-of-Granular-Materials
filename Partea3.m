%% Metode frecventiale.Asigurarea unei margini de faza impusa

%% Regulator PI
Hf = minreal(tf(KceKv*Kc*Kthetac,conv(conv([Tv 1],[Tc 1]),[Tthetac 1]),'iodelay',Tauc));
bode(Hf);
wtstelat = 0.00496;

Ti = 4/wtstelat;
Kp = sqrt(wtstelat^2+ 0.0625) * sqrt(wtstelat^2+0.0039) * sqrt(wtstelat^2+1.4811e-6) / 1.3723e-5;

Hpi = tf(Kp*[Ti 1],[Ti 0]);
Hd = series(Hf,Hpi);
bode(Hd)
%% Teste PD
Hf = minreal(tf(KceKv*Kc*Kthetac,conv(conv([Tv 1],[Tc 1]),[Tthetac 1]),'iodelay',Tauc));
bode(Hf);
beta = 0.12;
w0 = 0.0108;
modulHfjw0 = 1.3723e-5/sqrt(w0^2+0.25^2)/sqrt(w0^2+0.0625^2)/sqrt(w0^2+0.001217^2);
VR = sqrt(beta)/modulHfjw0;
Td = 1/w0/sqrt(beta);
TN = beta*Td;

HR = tf(VR*[Td 1],[TN 1]);
Hd = series(HR,Hf);
bode(Hd)

%% Regulator PID

Hf = minreal(tf(KceKv*Kc*Kthetac,conv(conv([Tv 1],[Tc 1]),[Tthetac 1]),'iodelay',Tauc));
bode(Hf);
wtstelat = 0.00389;
w0 = 0.0108;
modulHfjwt = 1.3723e-5/sqrt(wtstelat^2+0.25^2)/sqrt(wtstelat^2+0.0625^2)/sqrt(wtstelat^2+0.001217^2);
modulHfjw0 = 1.3723e-5/sqrt(w0^2+0.25^2)/sqrt(w0^2+0.0625^2)/sqrt(w0^2+0.001217^2);
VR = 0.228 / modulHfjw0;
T0 = 2*pi/w0;
Ti = 1.2*T0;
Td = 0.5*T0;
beta = 0.12;

HR = tf(VR*conv([Td 1],[Ti 1]),conv([beta*Td 1],[Ti 0]));
Hd = series(Hf,HR);
bode(Hd)

%% Metode cvasioptim
Hf = tf(Kap*K1*KMsteluta,conv([TAUap 1],conv([TM1 1],[TMsteluta 1])));

%% 1)Metoda modulului
Tsigma = TAUap;
Hdstelat = tf(1,conv([2*Tsigma 0],[Tsigma 1]));
HR = tf(conv([TM1 1],[TMsteluta 1]),[2*Tsigma*Kap*K1*KMsteluta 0]);

%% Calculul unui regulator PI apare la o singură constantă de timp preponderentă
Hf = tf(Kap*K1*KMsteluta,conv([TAUap 1],[(TMsteluta+TM1) 1]));
Tsigma = TAUap;
Hdstelat = tf(1,conv([2*Tsigma 0],[Tsigma 1]));
HR = Hdstelat / Hf;
Ho = feedback(series(Hf,HR),1)
step(Ho)
figure;
t=0:0.1:30;
lsim(Ho,t,t)
%% 2)Metoda simetriei
Tsigma = TAUap;
Hdstelat = tf([4*Tsigma 1],conv([8*Tsigma^2 0 0],[Tsigma 1]));
Hf = tf(Kap*K1*KMsteluta,conv([TAUap 1],[(TMsteluta+TM1) 1]));

HR = Hdstelat / Hf;
Ho = feedback(series(Hf,HR),1);

Hoprim = tf([TauT 1],KTohm*[2*Tsigma^2 2*Tsigma 1]);

step(Ho);
figure;
t=0:0.1:30;
lsim(Ho,t,t)

%%
subplot(212);
step(Hoprim);

%% Reglare in cascada
KTMsteluta = KTm/KTohm;
TTmsteluta = TTm-TTohm;
Ksteluta = K*Kg;
Tsteluta = TM1 + Taum + Tg;

Kf = KceKv*Kc*Kthetac;
Tf = Tv + Tc + Tthetac;

KthetaTstelat = KtetaT/Kthetac;
TthetaTstelat = TtetaT - Tthetac;

HRohm = tf(conv([TM1 1],[TMsteluta 1]),[2*TAUap*Kap*K1*KMsteluta 0]);
Hohm0 = tf(1,[2*TAUap 1]);

TTMstelutasteluta = TTmsteluta + 2*TAUap;
TsigmaQ = TTMstelutasteluta;

HRQ = tf([0.797 0.2539],[3.1386 0]);
tau1 = TM1;
tau2 = TMsteluta;
VRohm = 1/2/TAUap/Kap/K1/KMsteluta;
VRohm = Tsteluta/2/TTMstelutasteluta/KTMsteluta/Ksteluta

HRQ = tf(VRohm*conv([tau1 1],[tau2 1]),[1 0]);
tauiQ = Tsteluta;

HRQ = tf(VRohm*[tauiQ 1],[tauiQ 0]);

Mr = 0.1;
uohm = tf(KMsteluta * Mr * [2*Tsigma*TsigmaQ 2*Tsigma 0],conv([TMsteluta 1 0],[2*TAUap^2 2*TAUap 1]))
step(uohm,0.1)

%%
Hf = tf(Kap*K1,conv([TAUap 1],[TM1 1]));
bode(Hf);
wstelat = 73;
modulHfjwt = 15538 / sqrt(wstelat^2 + 100^2) / sqrt(wstelat^2 + 7.215^2);
VR = 1/modulHfjwt
Ti = 4/wstelat;
HR = tf(VR*[Ti 1],[Ti 0]);
bode(series(HR,Hf));

Hf1 = zpk(minreal(tf(Kap*K1,conv([TAUap 1],[TM1 1]))));
Hf2 = zpk(minreal(tf(KMsteluta,[TMsteluta 1])));
Hf3 = zpk(minreal(tf(KTMsteluta*Ksteluta,conv([TTmsteluta 1],[Tsteluta 1]))));

UQ1 = zpk(minreal(Hf2/(1+HR*Hf1*Hf2*Hf3)));
 
t = 0:0.001:0.1;
subplot(211);
step(uohm,t);
title('uohm')
subplot(212);
step(UQ1,t);
title('UQ1')

%% 5.3)
VRC = 0.9 * Tauc / Tf / Kf;
tauic = 3.3*Tauc;

H0c = zpk(minreal(tf(VRC*Kf*[2.3*Tauc 1],[3.3*Tauc*conv([Tf 1],[Tauc 1]) VRC*Kf*[2.3*Tauc 1]])));
B = 0.01257 * 0.003393;
T1 = 1/0.01257;
T2 = 1/0.003393;

K = VRC*Kf/B;
T0c =T1+T2-2.3*Tauc;

H0c = zpk(minreal(tf(K,[T0c 1])));
Hajutor = tf(KtetaT,[TtetaT 1],'iodelay',TauT);

Hf = series(H0c,Hajutor);
bode(Hf);
w0 = 0.669;
T0 = 2*pi/w0;
mk = db2mag(-9.83);
VRM = 0.75/mk;

HR = VRM;

%% 6)
T = 0.1 * (TtetaT + Tthetam);
T1steluta = T + Tg;
Hf = zpk(minreal(tf(KTMsteluta*Ksteluta,conv([2*TAUap 1],conv([TTmsteluta 1],[T1steluta 1])),'iodelay',Taum)))

%% 6.2)
Tmin = min(2*TAUap+TTmsteluta,T1steluta);
Tr = 2 * Tmin;
TB = 2 * Tr;

Bs = tf(1,conv([4*Tmin 0],[Tmin 1]));
Hfprim = zpk(minreal(tf(1311.9,conv([1 50],conv([1 0.2004],[1 0.04704])))));
HR1 = zpk(minreal(Bs/Hfprim));
Hozz = feedback(series(HR1,Hf),1);

HR1 = zpk(minreal(tf(T1steluta*[T1steluta 1],4*KTMsteluta*Ksteluta*(2*TAUap+TTmsteluta)*[T1steluta 0])));

Ti = T1steluta
VR = T1steluta / 4/KTMsteluta/Ksteluta / (2*TAUap + TTmsteluta);

Hozz = zpk(minreal(tf(1,[4*Tmin^2 4*Tmin 1],'iodelay',Taum)));
step(Hozz);

%%
Hdelay = tf(1,'iodelay',1);
Hfprim = series(tf(1311.9,conv([1 50],conv([1 0.2004],[1 0.04704]))),Hdelay);
wtstelat = 0.115;
[mag,~] = bode(Hfprim,wtstelat);
VR = 1 / mag;
Ti = 4/wtstelat;

HR = tf(VR*[Ti 1],[Ti 0]);
Hdelay = tf(1,'iodelay',1);
Hf = Hfprim;

Ho = feedback(series(HR,Hf),1);


Hoz = c2d(Ho,0.01,'tustin');

step(Hoz,Hozz);
legend("Regulatorul prin impunerea marginii de faza","Regulatorul cu predictie")
