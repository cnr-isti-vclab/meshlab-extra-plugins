# Oriented bounding box filter

This filter calculates the rotation matrix and the dimensions of the
optimal oriented bounding box (OBB) for a set of vertices. It creates
two new meshes, one with the convex hull for the vertices (this is a
needed initial step), then another one with the optimal OBB.

There is a known exact algorithm [1] for solving this problem, however
it is O(n^3) and it is difficult to implement. The algorithm this
filter implements instead is [2], which their authors call HYBBRID and
that is based on combining genetic and Nelder-Mead algorithms. The
matlab authors code [3] was used as reference for this implementation.

Author: Alfonso Sanchez-Beato (alfonsosanchezbeato_at_yahoo.es)

[1] O'Rourke, J. (1985). Finding minimal enclosing
boxes. International journal of computer & information sciences,
14(3), 183-199.

[2] Chang, C. T., Gorissen, B., & Melchior, S. (2011). Fast oriented
bounding box optimization on the rotation group SO(3,R).  ACM
Transactions on Graphics (TOG), 30(5), 1-16.

[3] https://github.com/chadogome/OptimalOBB
