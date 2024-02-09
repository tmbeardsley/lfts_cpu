# lfts_cpu
Langevin Field-Theoretic Simulations (L-FTS) of Diblock Copolymers on CPUs

See https://tbeardsley.com/projects/lfts/fts_gpu for a detailed discussion of this project.
A GPU-accelerated version of this code can be found at: https://github.com/tmbeardsley/lfts_gpu.

Input file format:<br>
Line1: N NA XN C Ndt isXeN<br>
Line2: mx my mz Lx Ly Lz<br>
Line3: n_eq n_st n_smpl loadType<br>
Lines 4->(M+3): W-(r)<br>
Lines (M+4)->(2M+3): w+(r)<br>

Parameter descriptions:<br>
N is the number of monomers in a single polymer chain (integer).<br>
NA is the number of monomers in the A-block of a polymer chain (integer).<br>
XN is the interaction strength between A and B-type monomers (double).<br>
C is the square root of the invariant polymerisation index, Nbar (double).<br>
Ndt is the size of the time step in the Langevin update of W-(r) (double).<br>
isXeN instructs the program whether the parameter XN is in terms of bare (isXeN=0) or effective (isXeN=1) chi (integer).<br>
mx, my, mz are the number of mesh points in the x, y, and z dimensions of the simulation box (integers).<br>
Lx, Ly, Lz are the dimensions of the simulation box (in units of the polymer end-to-end length, R0) in the x, y, and z dimensions (doubles).<br>
n_eq is the number of langevin steps performed to equilibrate the system (integer).<br>
n_st is the number of langevin steps performed after equilibration has ended, during which statistics are sampled (integer).<br>
n_smpl is the number of steps between samples taken in the statistics period (integer).<br>
loadType instructs the program whether to load the W-(r) and w+(r) fields from the proceeding file lines (loadType=1), start from a disordered state (loadType=0) or start from a (300) lamellar phase (loadType=2).<br><br>
M = (mx x my x mz) is the total number of mesh points, such that the proceeding 2*M lines of the file can hold W-(r) and w+(r) fields to load.
