Tf = TAUap + TM1;
Kf = Kap * K1;
T = TMsteluta / KMsteluta;

%% 2.2 Regulatorul P

Hf = tf(Kf/T,[Tf 1 0]);
sigma = 0.12;
tita = abs(log(sigma)/sqrt(pi^2+log(sigma)^2));
modulA = 1/4/sqrt(2)/tita/tita;
modulAdB = 20*log10(modulA);
F = 7.0011;
FN = -(F-modulAdB);
VR = 10^(FN/20);

Hc = VR;
Hd = Hf*Hc;

wn = sqrt(36.18);

Ho = zpk(feedback(Hd,1))

step(Ho)
title('Raspunsul la treapta')

figure
t = 0:0.1:30;
lsim(Ho,t,t);
title('Raspunsul la rampa')

%% 2.3 Regulator PI

Hf = tf(Kf/T,[Tf 1 0]);
sigma = 0.075;
tita = abs(log(sigma)/sqrt(pi^2+log(sigma)^2));
modulA = 1/4/tita/tita/sqrt(2);
modulAdB = 20*log10(modulA);

F = 7.001;
FN = -(F-modulAdB);
VR = 10^(FN/20);

Hc = VR;
Hd = series(Hf,Hc);

bode(Hd)
Ho = zpk(feedback(Hd,1));
step(Ho)
t = 0:0.1:20;
x = lsim(Ho,t,t);

wt = 3.65; 
wz = 0.1 * wt;
Tz = (1/wz);
cvstelat = 10;
cv = 1/(20-max(x));
wp = wz * (cv/cvstelat);
Tp = (1/wp);
VRpi = VR * cvstelat / cv;
Hc = tf(VRpi*[Tz 1],[Tp 1]);
Hdes = series(Hc,Hf);
Ho = feedback(Hdes,1);
Ho = tf(26,[1 6.49 26]);
step(Ho)
figure
lsim(Ho,t,t)



%% Regulator PD
Hf = tf(Kf/T,[Tf 1 0]);
sigma = 0.05;
tita = abs(log(sigma)/sqrt(pi^2+log(sigma)^2));
modulA = 1/4/tita/tita/sqrt(2);
modulAdB = 20*log10(modulA);
wf = 6.73;

F = 7.0011;
FN = -(F-modulAdB);
VR = 10^(FN/20);

Hc = VR;
Hd = Hf*Hc;

wt1 = 3.19;
tr = 2/tita/tita/wt1;
trstelat = 0.5;
wt2 = wt1 * (tr/trstelat);

Kpd = wt2/wt1;

Td = 0.1486;
TN = Td * (trstelat/tr);

HPD = tf(VR*Kpd*[Td 1],[TN 1]);
Hd = series(HPD,Hf);

Ho = feedback(Hd,1);
step(Ho)

figure
t = 0:0.1:30;
lsim(Ho,t,t);
title('Raspunsul la rampa')

%% Regulator PID 

Hf = tf(Kf/T,[Tf 1 0]);
sigma = 0.075;
tita = abs(log(sigma)/sqrt(pi^2+log(sigma)^2));
modulA = 1/4/tita/tita/sqrt(2);
modulAdB = 20*log10(modulA);
wf = 6.73;

F = 7.0011;
FN = -(F-modulAdB);
VR = 10^(FN/20);

Hc = VR;
Hd = Hf*Hc;

Ho = feedback(Hd,1);
step(Ho)
t = 0:0.1:30;
x = lsim(Ho,t,t);

wt1 = 3.65;
tsstelat = 0.5;
ts = 1.2;
wt2 = wt1 * (ts/tsstelat);
Td = Tf;
TN = Td * (tsstelat/ts);
cv = 1/(30-max(x));
cvstelat = 12;
VRpid = VR * (wt2/wt1) * (cvstelat/cv);
wz = 0.1*wt2;
Tz = 1/wz;
wp = wz * (cv/cvstelat);
Tp = 1/wp;

num = conv([Td 1],[Tz 1]);
den = conv([TN 1],[Tp 1]);
Hc = tf(VRpid * num, den);

Hdes = series(Hc,Hf);
Ho = feedback(Hdes,1);
Ho = tf(151,[1 15.52 151])
step(Ho);

figure;
lsim(Ho,t,t)