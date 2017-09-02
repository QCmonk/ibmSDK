include "qelib1.inc";
qreg q[5];
creg c[5];
gate icu3 (theta, phi, lambda) c,t {
x c;
u1 ((lambda-phi)/2) t;
cx c,t;
u3 (-theta/2,0,-(phi+lambda)/2) t;
cx c,t;
u3 (theta/2,phi,0) t;
x c;
}

x q[0];
cu3(pi,0,pi) q[0],q[1];
measure q[0] -> c[0];
measure q[1] -> c[1];