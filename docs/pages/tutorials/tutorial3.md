
# Tutorial: End-to-end analysis of the Shallow Water Equations with LSPG 

@m_class{m-block m-info}

@par
This tutorial works through an end-to-end analysis where we use the LSPG approach to accelerate a forward model of the shallow water equations (SWEs). In this tutorial, we will:
1. Interface an SWE solver, written with Eigen data structures, to Pressio.
2. Use Pressio's time marching schemes to execute solves of the SWEs to construct training data.
3. Use supporting Python scripts to analyze the training data and construct ROM basis vectors.
4. Construct and run a standard LSPG ROM for novel parameter training instances.
5. Construct and run a hyper-reduced LSPG ROM for novel parameter training instances.

@image html swetut_f1.gif width=90%

# The shallow water equations
We consider the shallow water equations on the spatial domain @f$\Omega = [-\frac{L}{2},\frac{L}{2}] \times  [-\frac{L}{2},\frac{L}{2}]@f$, which comprises the PDE system
@f[
\begin{split}
&\frac{\partial h}{\partial t} + \frac{\partial}{\partial x }(  h u) + \frac{\partial}{\partial y }( h v) = 0\\
&\frac{\partial h u}{\partial t} + \frac{\partial}{\partial x} (h u^2 + \frac{1}{2} \mu_1 h^2) + \frac{\partial}{\partial y }( h u v) = \mu_3 hv\\
&\frac{\partial h v}{\partial t} + \frac{\partial}{\partial x} (h u v) + \frac{\partial}{\partial y }( h v^2 +  \frac{1}{2} \mu_1 h^2) = \mu_3 hu.
\end{split}
@f]
In the above, @f$h : \Omega \rightarrow \mathbb{R}@f$ is the height of the water surface,  @f$u : \Omega \rightarrow \mathbb{R}@f$ is the x-velocity, and @f$v : \Omega \rightarrow \mathbb{R}@f$ is the y-velocity. The system has three parameters: 
@f$\mu_1@f$ is the gravity parameter, @f$\mu_2@f$ controls the magnitude of the initial pulse, and @f$\mu_3@f$ controls the magnitude of the Coriolis forcing. 

# Offline phase
In the offline phase, we
1. Interface with and run the full-order model
2. Extract training data and perform proper orthogonal decomposition to find a basis for the ROM
3. Select indices for hyper-reduction by Q-sampling

We now walk through this outline phase. To run the commands provided below, set the $SWE2D_DIR as 

```bash
export $SWE2D_DIR="build_location"/tutorials/swe2d
```

## Running the full-order model 
A full-order model that can be used to solve the SWEs is located in the packages/apps section of the Pressio repo (see [here](https://github.com/Pressio/pressio/tree/swe2d/packages/apps/src/swe2d) ). The full-order model employs a first-order finite volume 
discretization with the Rusanov flux scheme at the cell interfaces.  

The first step in our analysis is to run the full-order model for training parameter instances. To do this, we can write the driver file *run_fom_for_training_params.cc*. This file couples Pressio to the application, and uses Pressio's time marching schemes to solve the model. See [here](./md_pages_tutorials_tutorial3_fom.html) for a step-by-step walk through of constructing this driver file. In summary, the driver file executes the full-order model for 9 parameter instances on the grid @f$\mu_1 \times \mu_3 = [3,6,9]\times [0.05,0.15,0.25]@f$, and saves the solutions to file.

To run the driver file, move to the offline_phase directory and run the script: 
```bash
cd $SWE2D_DIR/offline_phase
./run_fom_for_training_params
```
This will take some time to run, approximately 30 minutes. If successful, a series of *solution#.bin* files should have been written. These solution files contain the FOM solutions at every time step for each of the nine training parameter instances. To view the results of one sample simulation, we can go to the supporting_python_scripts directory and run the *viewSolutionAndMakePlots.py* script
```bash
cd $SWE2D_DIR/offline_phase/supporting_python_scripts
python viewSolutionAndMakePlots.py
```
This script will bring up a live animation of the solution for the first parameter instance.

## Extracting the bases and building the sample mesh
We now need to construct the basis vectors used in the ROM. To do this, we again move to the supportingPythonFiles directory and run the *makeBasisAndHyperReducedBasis.py* script
```bash
cd $SWE2D_DIR/offline_phase/supporting_python_scripts/
python makeBasisAndHyperReducedBasis.py
```
This script loads in the snapshots and performs POD to obtain the ROM basis. Additionally, this script selects cells for the sample mesh employed in hyper-reduction, and saves the relevant information of this sample mesh to file. Specifically, it makes the following files:
1. *info.txt* This file contains information on the size of the ROM, the size of the sample mesh, and the size of the sample mesh and stencil mesh. 
2. *basis.txt* This file contains the basis vectors for the ROM on the global mesh 
3. *sample_mesh_gids.txt* This file contains the global IDs of the indices used for the sample mesh
4. *sample_mesh_plus_stencil_gids.txt* This file contains the global IDs of the indices used for the sample *and* stencil mesh
5. *PhiSamplePlusStencil.txt* This file contains the ROM basis, but only at the sample mesh plus stencil mesh
Additionally, this script will create a file, *samplemesh.png*, depicting the sample and stencil mesh. Cells in black are the sample mesh, while cells in red are on the stencil mesh.

@image html samplemesh.png width=50%

# Online phase
With the offline stage complete, we can now run our ROMs for novel parameter instances. We will first run a standard LSPG ROM without hyper-reduction, followed by an LSPG ROM with hyper-reduction. To set a novel parameter instance, we switch to the online directory and look at the *novel_params.txt* file
```bash
cd $SWE2D_DIR/online_phase/
vim novel_params.txt
```
By default, we have the novel parameter instance set to be @f$\mu_1 = 7.5, \; \mu_2=0.125, \; \mu_3 = 0.2@f$. The rest of this tutorial will present results for this parameter instance, but the user is encouraged to play around with different parameters and see how it impacts the results. Before we run the ROM, we first run a FOM for this new parameter instance so we can assess the accuracy of our ROM (of course, in a practical scenario we would not do this step!). We do not provide a detailed explanation on this driver script, since it closely follows that written previously. To run the FOM for our new parameter instance, we do the following: 
 ```bash
cd $SWE2D_DIR/online_phase/fom
./run_fom
```
This will run the FOM and save the solution to file. The FOM was tested on a 2.7 GHz 12-Core Intel Xeon E5 core, and took 152 seconds to run. 

##LSPG ROM 
To run an LSPG ROM, we write a driver file, called [run_lspg.cc](https://github.com/Pressio/pressio-tutorials/blob/swe2d_tutorial/tutorials/swe2d/online_phase/lspg_rom/run_lspg.cc). See [here](./md_pages_tutorials_tutorial3_lspg.html) for a step-by-step walk-through of constructing this driver file. In summary, this script couples the application to Pressio, loads in the basis information we generated in the offline phase, and couples to Pressio's ROM capabilities to run an LSPG ROM. 

To run the LSPG ROM, we move to the *lspg_rom* directory, copy our ROM basis, and run the ROM, 

```bash
cd $SWE2D_DIR/online_phase/lspg_rom
cp ../../offline_phase/supporting_python_scripts/basis.txt .
./lspg_rom
python viewSolutionAndMakePlots.py
```
This process saves the generalized coordinates of the ROM to the *solution.bin* file, and *viewSolutionAndMakePlots.py* plots the height of the water surface for a given spatial location as a function of time, and saves the plot to *result.png*. This plot looks as follows:
@image html result_lspg.png width=50%
The ROM was tested on a 2.7 GHz 12-Core Intel Xeon E5 core, and took 179 seconds to run. We immediately note that our *ROM is slower than the FOM!* This, of course, is due to the well known bottleneck associated with nonlinear systems. To gain computational speedups, we need hyper-reduction. We now detail this. 

##Hyperreduced LSPG ROM 
We now run construct and run a hyper-reduced LSPG ROM. To do this, we again need to write a driver file, which here we call [run_lspg_with_hyperreduction.cc](https://github.com/Pressio/pressio-tutorials/blob/swe2d_tutorial/tutorials/swe2d/online_phase/lspg_hyperReducedRom/run_lspg_with_hyperreduction.cc). A step-by-step tutorial for what is entailed in constructing this driver file is provided [here](./md_pages_tutorials_tutorial3_lspg_hyper.html). In summary, this file loads the basis on the *stencil mesh*, loads in information about the *sample mesh* and *stencil mesh*, and then constructs and runs an LSPG ROM employing the collocation hyper-reduction technique.
  
To run the LSPG ROM with hyper-reduction, we move to the *lspg_hyperReducedRom* directory, copy over our basis and sample mesh information, and then run our ROM.
```bash
cd $SWE2D_DIR/online_phase/lspg_hyperReducedRom
cp ../../offline_phase/supporting_python_scripts/*.txt .
./run_lspg_with_hyperreduction 
python viewSolutionAndMakePlots.py
```
If successful, the following plot will be generated. 
@image html result_lspgHyper.png width=50%
Running the ROM on the same 2.7 GHz 12-Core Intel Xeon E5 core machine took 20 seconds, which is about a 7.5x speedup over the FOM!

This completes our tutorial on ROMs for the shallow water equations. 
