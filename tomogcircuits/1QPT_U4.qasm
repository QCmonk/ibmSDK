OPENQASM 2.0;

include "qelib1.inc";
qreg q[5];
creg c[5];

//PROCESS

u3(0.194657*pi, 0.0408352*pi, 0.703968*pi) q[0];
u3(0.436393*pi, 0.622953*pi, 0.674483*pi) q[0];
u3(0.90592*pi, 0.407216*pi, 0.921755*pi) q[0];
u3(0.184731*pi, 0.867235*pi, 0.209647*pi) q[0];
