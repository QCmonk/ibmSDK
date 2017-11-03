OPENQASM 2.0;

include "qelib1.inc";
qreg q[5];
creg c[5];

h q[1];
h q[3];
cx q[1], q[0];
cx q[3], q[2];

barrier q;

s q[1];
barrier q;

cx q[2],q[1];
h q[1];
h q[2];
cx q[2], q[1];
h q[1];
h q[2];
cx q[2], q[1];

barrier q;

t q[1];