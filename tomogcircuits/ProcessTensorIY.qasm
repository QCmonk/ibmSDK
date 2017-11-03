OPENQASM 2.0;

include "qelib1.inc";
qreg q[5];
creg c[5];

//PROCESS

h q[1];
h q[3];
cx q[1], q[0];
cx q[3], q[2];

id q[0];

barrier q;

cx q[2],q[0];
h q[0];
h q[2];
cx q[2], q[0];
h q[0];
h q[2];
cx q[2], q[0];

barrier q;

y q[0];