OPENQASM 2.0;

include "qelib1.inc";
qreg q[5];
creg c[5];

h q[0];
h q[3];
cx q[0], q[1];
cx q[3], q[2];

id q[0];

cx q[0],q[2];
h q[0];
h q[2];
cx q[0], q[2];
h q[0];
h q[2];
cx q[0], q[2];

u3(0.5*pi, 0.25*pi, 0.75*pi) q[0];
u3(0.25*pi, 0.5*pi, 0.25*pi) q[1];
