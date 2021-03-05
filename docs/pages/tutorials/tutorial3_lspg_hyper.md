
# Tutorial: End-to-end analysis of the Shallow Water Equations with LSPG: Writing the hyperreduced LSPG ROM 
This page walks through writing a driver file for an LSPG ROM with collocation-based hyper-reduction of the shallow water equations. The full code for this coupler is available in the [pressio-tutorials repo](https://github.com/Pressio/pressio-tutorials/blob/swe2d_tutorial/tutorials/swe2d/online_phase/lspg_hyperReducedRom/run_lspg_with_hyperreduction.cc). The first step in writing our driver is to include the relevant headers from the Pressio library. Here, we include [pressio_apps.hpp](https://github.com/Pressio/pressio/blob/master/packages/pressio_apps.hpp), and [pressio_rom_lspg.hpp](https://github.com/Pressio/pressio/blob/master/packages/pressio_lspg.hpp), which enable us to access our Shallow Water application and LSPG capabilities, respectively. 

```cpp
@codesnippet
../../../../pressio-tutorials/tutorials/swe2d/online_phase/lspg_hyperReducedRom/run_lspg_with_hyperreduction.cc
1:4
```

The next step for us is to start out int main file, and initialize the types for our FOM, as well as the scalar type 
```cpp
@codesnippet
../../../../pressio-tutorials/tutorials/swe2d/online_phase/lspg_hyperReducedRom/run_lspg_with_hyperreduction.cc
61:66
```
Note that, here we use the [swe2d_hyper](.) application, opposed to the [swe2d](.) application. These two applications are the same, with the exception that the *swe2d_hyper* application loads in indices for the sample and stencil mesh, and only computes the velocity/Jacobians at these points. 

The next thing that we will do is read in information on our ROM, and information for hyper-reduction. First, we read in *info_file.txt*. We wrote this file in the [training phase](.), and it contains information on the size of the ROM, sample mesh, and sample plus stencil mesh, respectively:
```cpp
@codesnippet
../../../../pressio-tutorials/tutorials/swe2d/online_phase/lspg_hyperReducedRom/run_lspg_with_hyperreduction.cc
68:84
```
Next we create data structures that will contain our sample mesh information: 
```cpp
@codesnippet
../../../../pressio-tutorials/tutorials/swe2d/online_phase/lspg_hyperReducedRom/run_lspg_with_hyperreduction.cc
85:90
```
In the above, we created
1. An *std::vector<int> sm_gids* that will contain the global x,y indices of the sample mesh indices (this is required by the app). Note that this vector only contains information on the x,y cells used in the sample mesh, but not the number of conserved variables. 
2. An *std::vector<int> smps_gids* that will contain the global IDs of the stencil mesh (this is also required by the app). Note that this vector only contains information on the x,y cells used in the stencil mesh, but not the number of conserved variables.  
3. An *Eigen::VectorXd sm_rel_gids* that will contain the IDs of the sample mesh *relative* to the IDs of the stencil mesh. Note that this vector will contain information about not only the x,y cells, but also the conserved variables, and hence it has *3sampleMeshSize* entries. This data is required by Pressio, and is discussed more below.

In point 3, we create an Eigen vector that contains information about the sample mesh indices, relative to the stencil mesh. What do we mean by this? For example, let's say our FOM has the (zero-based) indices 
@f[ \mathcal{I}_{\text{FOM}} = \begin{Bmatrix} 0& 1 & 2 & 3 & 4 & 5 & 6 & 7 & 8 & 9  \end{Bmatrix}@f] 
Let's say our sample mesh uses indices 3, 5, and 7,
 @f[ \mathcal{I}_{\text{sample}} = \begin{Bmatrix}  3 &  5  & 7 \end{Bmatrix}@f]
Let's suppose we use a three-point stencil in our numerical method. In this case, the stencil mesh would be
  @f[ \mathcal{I}_{\text{stencil}} = \begin{Bmatrix}  2& 3 & 4& 5 &6  & 7 & 8\end{Bmatrix}@f]
The indices of the *sample mesh* relative to the indices of the *stencil mesh* are thus
@f[ \mathcal{I}_{\text{rel}} =\begin{Bmatrix} 1 & 3 & 5 \end{Bmatrix}@f] 
Why do we need this information? In Pressio, when we perform hyperreduction, we need to load in the basis on the *stencil mesh* so that we can reconstruct the state at these points, and eventually evaluate the velocity/residual. However, we only then evaluate the residual at the *sample mesh*, and thus we need to know what indices in the *stencil mesh* correspond to those that are on the *sample mesh*. Makes sense? Great!

We load in our sample mesh, stencil mesh, and relative sample mesh information into these vectors. Recall that we created this information in the offline stage. We load  this information in by
```cpp
@codesnippet
../../../../pressio-tutorials/tutorials/swe2d/online_phase/lspg_hyperReducedRom/run_lspg_with_hyperreduction.cc
90:126
```
In our final step for the sample mesh, we create a wrap our sample mesh relative IDs with a Pressio type,
```cpp
@codesnippet
../../../../pressio-tutorials/tutorials/swe2d/online_phase/lspg_hyperReducedRom/run_lspg_with_hyperreduction.cc
127:128
```

With these steps complete, we can now proceed with constructing our application. We define the number of grid points, domain size, and load in our testing parameter from file,
```cpp
@codesnippet
../../../../pressio-tutorials/tutorials/swe2d/online_phase/lspg_hyperReducedRom/run_lspg_with_hyperreduction.cc
129:143
```
We note that, in the last line where we construct the application, the construct additionally requires the sample mesh information that we had loaded in.

Next, we read in the basis. Note that for hyper-reduction, we only read in the basis on the *stencil mesh*. The basis on the full mesh is not required.
```cpp
@codesnippet
../../../../pressio-tutorials/tutorials/swe2d/online_phase/lspg_hyperReducedRom/run_lspg_with_hyperreduction.cc
144:151
```

We then use this basis information to construct the decoder. For hyper-reduced LSPG the decoder, in essence, is used to compute the product @f$\boldsymbol \Phi_{\text{stencil}} \hat{\boldsymbol x} @f$, where @f$\boldsymbol \Phi_{\text{stencil}}@f$ is the basis on the stencil mesh. 
```cpp
@codesnippet
../../../../pressio-tutorials/tutorials/swe2d/online_phase/lspg_hyperReducedRom/run_lspg_with_hyperreduction.cc
152:158
```

Next we create our reference state. For hyper-reduction, Pressio requires the reference state on the *stencil_mesh*. If we want to view the entire solution field, however, we will also need the reference state on the full mesh. Here, we create both entities:
```cpp
@codesnippet
../../../../pressio-tutorials/tutorials/swe2d/online_phase/lspg_hyperReducedRom/run_lspg_with_hyperreduction.cc
160:164
```

We can now create the LSPG problem. First, we create the ROM state vector and initialize it to zero (again, the initial conditions are the reference state),
```cpp
@codesnippet
../../../../pressio-tutorials/tutorials/swe2d/online_phase/lspg_hyperReducedRom/run_lspg_with_hyperreduction.cc
165:174
```

Next, we define our LSPG problem. We again use the Crank Nicolson time marching scheme, and create our LSPG problem with *pressio::rom::lspg::createHyperReducedProblemUnsteady*:
```cpp
@codesnippet
../../../../pressio-tutorials/tutorials/swe2d/online_phase/lspg_hyperReducedRom/run_lspg_with_hyperreduction.cc
165:174
```
Note that in the constructor we pass the relative sample mesh indices array.

The remainder of our driver follows closely that which was written for LSPG. We define our linear and nonlinear solvers
```cpp
@codesnippet
../../../../pressio-tutorials/tutorials/swe2d/online_phase/lspg_hyperReducedRom/run_lspg_with_hyperreduction.cc
179:189
```

Then we finally define our time steps sizes, create our observer, and solve the ROM
```cpp
@codesnippet
../../../../pressio-tutorials/tutorials/swe2d/online_phase/lspg_hyperReducedRom/run_lspg_with_hyperreduction.cc
192:209
```


This completes our description of writing the LSPG hyper-reduction coupler. Click [here](./md_pages_tutorials_tutorial3.html) to return to the SWE tutorial.


