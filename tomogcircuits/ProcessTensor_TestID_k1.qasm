OPENQASM 2.0;

include "qelib1.inc";
qreg q[5];
creg c[5];


h q[0];
h q[3];
cx q[0],q[1];
cx q[3],q[2];

x q[0];

cx q[0],q[2];
h q[0];
h q[2];
cx q[0],q[2];
h q[0];
h q[2];
cx q[0],q[2];

u3(0.123*pi, 0.167*pi, 0.099*pi) q[0];
u3(0.401*pi, 0.322*pi, 0.222*pi) q[2];
