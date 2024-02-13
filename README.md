# lfts_cpu
Langevin Field-Theoretic Simulation of Diblock Copolymers on CPUs

See https://www.tbeardsley.com/projects/lfts/fts_gpu for a detailed discussion of this project.<br>

<b>comp.sh:</b><br>
An example compile script

<b>Running the program:</b><br>
./<program_name> <input_file_name><br>

<b>input:</b><br>
An example input file.

<b>Input file format:</b><br>
Line 1: <em>N NA XN C Ndt isXeN</em><br>
Line 2: <em>mx my mz Lx Ly Lz</em><br>
Line 3: <em>n_eq n_st n_smpl loadType</em><br>
Lines 4->(M+3): W-(r)<br>
Lines (M+4)->(2M+3): w+(r)<br>

Note: A real-space position r = (x,y,z) corresponds to a mesh point position r_m = (i,j,k), where i=0->mx-1, j=0->my-1 and k=0->mz-1 are integers. The elements of the fields, W-(r) and w+(r), are then written in ascending order of the row-major index: p = mx\*(i\*my+j)+k.

<b>Parameter descriptions:</b><br>
<em>N</em> is the number of monomers in a single polymer chain (integer).<br>
<em>NA</em> is the number of monomers in the A-block of a polymer chain (integer).<br>
<em>XN</em> is the interaction strength between A and B-type monomers (double).<br>
<em>C</em> is the square root of the invariant polymerisation index, Nbar (double).<br>
<em>Ndt</em> is the size of the time step in the Langevin update of W-(r) (double).<br>
<em>isXeN</em> instructs the program whether the parameter XN is in terms of bare (isXeN=0) or effective (isXeN=1) chi (integer).<br>
<em>mx, my, mz</em> are the number of mesh points in the x, y, and z dimensions of the simulation box (integers).<br>
<em>Lx, Ly, Lz</em> are the dimensions of the simulation box (in units of the polymer end-to-end length, R0) in the x, y, and z dimensions (doubles).<br>
<em>n_eq</em> is the number of langevin steps performed to equilibrate the system (integer).<br>
<em>n_st</em> is the number of langevin steps performed after equilibration has ended, during which statistics are sampled (integer).<br>
<em>n_smpl</em> is the number of steps between samples being taken in the statistics period (integer).<br>
<em>loadType</em> instructs the program whether to load the W-(r) and w+(r) fields from the proceeding file lines (loadType=1), start from a disordered state (loadType=0) or start from a (300) lamellar phase (loadType=2).<br><br>
M = (mx\*my\*mz) is the total number of mesh points, such that the proceeding 2*M lines of the file can hold W-(r) and w+(r) fields to load.

<b>Output files:</b><br>
"w_eq_<step_number>": The state of the W-(r) and w+(r) fields at simulation step number <step_number> during the equilibration period. First three lines are simulation parameters so it can be used as an input file.<br>
"w_st_<step_number>": The state of the W-(r) and w+(r) fields at simulation step number <step_number> during the statistics gathering period. First three lines are simulation parameters so it can be used as an input file.<br>
"phi_eq_<step_number>": The state of the phi-(r) and phi+(r) fields at simulation step number <step_number> during the equilibration period.<br>
"phi_eq_<step_number>": The state of the phi-(r) and phi+(r) fields at simulation step number <step_number> during the statistics gathering period.<br>

