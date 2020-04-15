1. IKE test case, there is one inconsistency in the node number.

ADCIRC: node number start from 1
SSWM: node number start from 0

In IKE, the node number is named after ADCIRC.

A). When plot density plots, I pick the node number of SSWM correctly, e.g. SSWM=3701 whilst ADCIRC=3702
B). When plot ADCIRC line mean plots, I pick the node number of SSWM wrong, e.g. SSWM=3702 whilst ADCIRC=3702, which is not the same spatial point.


2. In Harvey test case. I did it right.

In HARVEY, the node number is named after SSWM. Everything is one less than ADCIRC.
