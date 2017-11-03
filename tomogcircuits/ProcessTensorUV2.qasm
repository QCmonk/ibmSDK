OPENQASM 2.0;

include "qelib1.inc";
qreg q[5];
creg c[5];

//PROCESS

h q[1];
h q[3];
cx q[1], q[0];
cx q[3], q[2];

u3(0.22*pi, 0.34*pi, 0.78*pi) q[0];

barrier q;

cx q[2],q[0];
h q[0];
h q[2];
cx q[2], q[0];
h q[0];
h q[2];
cx q[2], q[0];

barrier q;

z q[0];
x q[0];