OPENQASM 2.0;

include "qelib1.inc";
qreg q[5];
creg c[5];

//PROCESS



cx q[2],q[0];
h q[2];
h q[0];
cx q[2], q[0];
h q[2];
h q[0];
cx q[2], q[0];
