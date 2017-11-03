OPENQASM 2.0;

include "qelib1.inc";
qreg q[5];
creg c[5];

//PROCESS

h q[1];
h q[3];
cx q[1], q[0];
cx q[3], q[2];