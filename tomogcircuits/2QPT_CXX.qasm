OPENQASM 2.0;

include "qelib1.inc";
qreg q[5];
creg c[5];

//PROCESS

cx q[1],q[0];
x q[1];
