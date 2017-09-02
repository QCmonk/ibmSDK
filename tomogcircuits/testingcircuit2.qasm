OPENQASM 2.0;

include "qelib1.inc";
qreg q[5];
creg c[5];
u1((0.225-0)/2) q[1];
cx q[0],q[2];


