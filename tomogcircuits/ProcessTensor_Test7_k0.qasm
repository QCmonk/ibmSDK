OPENQASM 2.0;

include "qelib1.inc";
qreg q[5];
creg c[5];

h q[1];
h q[4];
cx q[1],q[2];
cx q[4],q[3];





cx q[0],q[1];
h q[0];
h q[1];
cx q[0],q[1];
h q[0];
h q[1];
cx q[0],q[1];

u3(0.123*pi, 0.167*pi, 0.099*pi) q[0];

