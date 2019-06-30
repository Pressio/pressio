#include<assert.h>
#include<math.h>
#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include"cblas.h"

#include"pde.h"

int N_GRID=0, N_CHUNK=0, N_STEP=0, N_MODE=1;
int OBJ_IND=0, GRAD_IND=0;
int PROJ=0;
int SHAD_PROJ=0, MEM=0;
double S_CONST[5]={ 0.0 };
double DT_STEP=0;
double JBAR[5]={ 0.0 };
double T_TOTAL=0;
double *** SOLN_U=0;
double *** SOLN_V=0;
double *** SOLN_W=0;
double * S_FIELD=0;

// Runge Kutta coefficients
double RK[3][3] = {{ 1./2,  0,    0 },
                   {-1./6,  1./3, 0 },
                   { 0,    -2./3, 1.}};


// Linear Convection Operator (1st spatial derivative)
void
du_dx(const double * u, double * dudt)
{
    double dx = 128. / (N_GRID + 1);
    for (int i = 0; i < N_GRID; ++i)
    {
        double up = (i < N_GRID - 1) ? u[i+1] : 0;//u[0];
        double um = (i > 0) ? u[i-1] : 0; //u[N_GRID-1];
        double ux = (up - um) / (2 * dx);
        dudt[i] = ux;
    }
}

// Non-Linear Convection Operator
void
du2_dx(const double * u, double * dudt)
{
    double dx = 128. / (N_GRID + 1);
    for (int i = 0; i < N_GRID; ++i)
    {
        double up = (i < N_GRID - 1) ? u[i+1] : 0;//u[0];
        double um = (i > 0) ? u[i-1] : 0;//u[N_GRID-1];
        double u2x = (up*up - um*um) / (2 * dx);
        dudt[i] = u2x;
    }
}

// 2nd order Diffusion operator
void
d2u_dx2(const double * u, double * dudt)
{
    double dx = 128. / (N_GRID + 1);
    for (int i = 0; i < N_GRID; ++i)
    {
        double up = (i < N_GRID - 1) ? u[i+1] : 0;//u[0];
        double um = (i > 0) ? u[i-1] : 0;//u[N_GRID-1];
        double uxx = (up + um - 2 * u[i]) / (dx * dx);
        dudt[i] = uxx;
    }
}

// 4th order anti-diffusion operator
void
d4u_dx4(const double * u, double * dudt)
{
    double dx = 128. / (N_GRID + 1);
    for (int i = 0; i < N_GRID; ++i)
    {
        double up = (i < N_GRID - 1) ? u[i+1] : 0;//u[0];
        double upp = (i < N_GRID - 2) ? u[i+2] : u[i];//u[0];
        double um = (i > 0) ? u[i-1] : 0;//u[N_GRID-1];
        double umm = (i > 1) ? u[i-2] : u[i];//u[N_GRID-1];
        double uxx = (up + um - 2 * u[i]) / (dx * dx);
        double uxxp = (upp + u[i] - 2 * up) / (dx * dx);
        double uxxm = (umm + u[i] - 2 * um) / (dx * dx);
        double uxxxx = (uxxp + uxxm - 2 * uxx) / (dx * dx);
        dudt[i] = uxxxx;
    }
}

// pde right hand side, assigns f(u) to dudt
void
ddt(const double * u, double * dudt)
{
    double dx = 128. / (N_GRID + 1);
	double ux[N_GRID];
    double u2x[N_GRID];
    double uxx[N_GRID];
    double uxxxx[N_GRID];
    du_dx(u,ux);
	du2_dx(u, u2x);
    d2u_dx2(u,uxx);
    d4u_dx4(u,uxxxx);
    
    //double x = 0.0;
    for (int i = 0; i < N_GRID; ++i)
    {
        dudt[i] = -S_FIELD[i] * ux[i] - 0.5 * u2x[i] - uxx[i] - uxxxx[i];
        //x = ((i + 1) * dx - 64.) / 64.;
        //dudt[i] = -(S_CONST[0] + S_CONST[1] * x) * ux[i] - 0.5 * u2x[i] - uxx[i] - uxxxx[i];
    }
}


// TODO: Update for multiple objectives
// Objective function J 
double
Obj(const double * u, int obj_ind)
{
    double J = 0;
    if(obj_ind == 1)
    // time averaged u(x=32)
    {
        J = u[N_GRID/4];
    }
    else if(obj_ind == 2)
    // time averaged u(x=96)
    {
        J = u[3*N_GRID/4];
    }
    else if(obj_ind == 3)
    // time averaged u(x=32)^2
    {
        J = pow(u[N_GRID/4],2);
    }
    else if(obj_ind == 4)
    // time averaged u(x=96)^2
    {
        J = pow(u[3*N_GRID/4],2);
    }
    else
    {
    // time and space averaged u
    	for (int i = 0; i < N_GRID; ++i) {
            J = J + u[i] / (N_GRID + 2.0);
        }
    }

	return J;
}

// Derivative of objective function J with respect to the primal u
// Needed for computing forcing in adjoint solver
void
ddObj(const double * u, double * dJdu, int obj_ind)
{
    if(obj_ind == 1)
    // Sensitivity of time averaged u(x=32)
    {
        for (int i = 0; i < N_GRID; ++i) {
            dJdu[i] = 0.;
        }
        dJdu[N_GRID/4] = 1.;
    }
    else if(obj_ind == 2)
    // Sensitivity of time averaged u(x=96)
    {
        for (int i = 0; i < N_GRID; ++i) {
            dJdu[i] = 0.;
        }
        dJdu[3*N_GRID/4] = 1.;
    }
    else if(obj_ind == 3)
    // Sensitivity of time averaged u(x=32)^2
    {
        for (int i = 0; i < N_GRID; ++i) {
            dJdu[i] = 0.;
        }
        dJdu[N_GRID/4] = 2. * u[N_GRID/4];
    }
    else if(obj_ind == 4)
    // Sensitivity of time averaged u(x=96)^2
    {
        for (int i = 0; i < N_GRID; ++i) {
            dJdu[i] = 0.;
        }
        dJdu[3*N_GRID/4] = 2. * u[3*N_GRID/4];
    }
    else
    // Sensitivity of time and space averaged u
    {
        for (int i = 0; i < N_GRID; ++i) {
            dJdu[i] = 1. / (N_GRID + 2.0);
        }
    }

}

// Derivative of f(u) with respect to s, df/ds.  
// where f is the KS equation right hand side (see the ddt function)
// and s is the parameter S_CONST in the ddt function
void
ddtds(const double * u, double * dfds, int grad_ind)
{
    double dx = 128. / (N_GRID + 1);
    double x = 0.0;
    for (int i = 0; i < N_GRID; ++i)
    {
        double up = (i < N_GRID - 1) ? u[i+1] : 0;// u[0];
        double um = (i > 0) ? u[i-1] : 0;//u[N_GRID-1];
        x = ((i + 1) * dx - 64.) / 64.;
        //x = ((i + 1) * dx * 3.14159265359) / 64.;

        if(grad_ind == 1)
        //Sensitivity to 1st mode of s(x) 
        {
            dfds[i] = -x * (up - um) / (2 * dx);
            //dfds[i] = - cos(x) * (up - um) / (2 * dx);

        }
        else if(grad_ind == 2)
        //Sensitivity to 2nd mode of s(x)
        {
            dfds[i] = -0.5 * (3.0 * pow(x,2) - 1.0) * (up - um) / (2 * dx);
            //dfds[i] = - cos(2*x) * (up - um) / (2 * dx);
        }
        else if(grad_ind == 3)
        //Sensitivity to 3rd mode of s(x)
        {
            dfds[i] = -0.5 * (5.0 * pow(x,3) - 3 * x) * (up - um) / (2 * dx);
            //dfds[i] = - cos(3*x) * (up - um) / (2 * dx);

        }
        else if(grad_ind == 4)
        //Sensitivity to 4th mode of s(x)
        {
            dfds[i] = -0.125 * (35.0 * pow(x,4) - 30 * pow(x,2) + 3.0) * (up - um) / (2 * dx);
            //dfds[i] = - cos(4*x) * (up - um) / (2 * dx);
            //dfds[i] = - uxxxx[i];
        }
        else
        //Sensitivity to 0th mode of s(x)
        {
            dfds[i] = - (up - um) / (2 * dx);
        }


    }
}




// Linearized Convection Operator
void
udv_dx(const double * u, const double * v, double * dvdt)
{
    double dx = 128. / (N_GRID + 1);
    for (int i = 0; i < N_GRID; ++i)
    {
        double up = (i < N_GRID - 1) ? u[i+1] : 0;//u[0];
        double um = (i > 0) ? u[i-1] : 0; //u[N_GRID-1];
        double vp = (i < N_GRID - 1) ? v[i+1] : 0; //u[0];
        double vm = (i > 0) ? v[i-1] : 0; //u[N_GRID-1];
        double v2x = 2 * (up*vp - um*vm) / (2 * dx);
        dvdt[i] = v2x;
    }
}


// Tangent equation rhs, assigns (df/du) * v + inhomo * (df/ds) to dvdt
// where f is the pde right hand side (see the ddt function)
// and s is the parameter S_CONST in the ddt function
void
ddtTan(const double * u, const double * v, double * dvdt, int inhomo)
{
    double dx = 128. / (N_GRID + 1);
    double vx[N_GRID];
    double uvx[N_GRID];
    double vxx[N_GRID];
    double vxxxx[N_GRID];
    du_dx(v,vx);
    udv_dx(u,v,uvx);
    d2u_dx2(v,vxx);
    d4u_dx4(v,vxxxx);
	double dfds[N_GRID];
	ddtds(u, dfds, GRAD_IND);
    
    for (int i = 0; i < N_GRID; ++i)
    {
        dvdt[i] = -S_FIELD[i] * vx[i] - 0.5 * uvx[i] - vxx[i] - vxxxx[i] + inhomo * dfds[i];
        //dvdt[i] = -(S_CONST[0] + S_CONST[1]*x) * vx[i] - 0.5 * uvx[i] - vxx[i] - vxxxx[i] + inhomo * dfds[i];
    }
}


// Adjoint equation rhs, assigns -(df/du)' * w to dwdt
// where f is the pde equation right hand side (see the ddt function)
void
ddtAdj(const double * u, const double * w, double * dwdt)
{
    double dx = 128. / (N_GRID + 1);
    double wx[N_GRID];
    double lwx[N_GRID];
    double wxx[N_GRID]; 
    double wxxxx[N_GRID];
    du_dx(w,wx);
    udv_dx(S_FIELD,w,lwx);
    d2u_dx2(w,wxx);
    d4u_dx4(w,wxxxx);
    for (int i = 0; i < N_GRID; ++i) {
        dwdt[i] = - 0.5 * lwx[i] - u[i] * wx[i] + wxx[i] + wxxxx[i];
        //dwdt[i] = -(S_CONST[0] + S_CONST[1]*x + u[i]) * wx[i] + wxx[i] + wxxxx[i];
        //dwdt[i] = -S_CONST * wx[i] + wxx[i] + wxxxx[i];
    }
}


// dual consistent explicit RK (See Shan Yang's thesis)
void
stepPrimal(const double * u0, double * u, double dt)
{
    double dudt0[N_GRID], dudt1[N_GRID], *dudt2 = dudt0;

    memmove(u, u0, sizeof(double) * N_GRID);

    ddt(u, dudt0);
    cblas_daxpy(N_GRID, dt * RK[0][0], dudt0, 1, u, 1);

    ddt(u, dudt1);
    cblas_daxpy(N_GRID, dt * RK[1][0], dudt0, 1, u, 1);
    cblas_daxpy(N_GRID, dt * RK[1][1], dudt1, 1, u, 1);

    ddt(u, dudt2);
    cblas_daxpy(N_GRID, dt * RK[2][1], dudt1, 1, u, 1);
    cblas_daxpy(N_GRID, dt * RK[2][2], dudt2, 1, u, 1);
}

// Tangent of stepPrimal
void
stepTangent(const double * u0, const double * v0, double * v, double dt,
            int inhomo)
{
    double u1[N_GRID], u2[N_GRID], dudt0[N_GRID], dudt1[N_GRID];
    double * dvdt0 = dudt0, * dvdt1 = dudt1, * dvdt2 = dudt0;

    // Tangent forcing
    //double dfds[N_GRID];
	//ddtds(u0, dfds, GRAD_IND);
    //cblas_daxpy(N_GRID, inhomo * dt, dfds, 1, v, 1);

    // Primal
    ddt(u0, dudt0);
    memmove(u1, u0, sizeof(double) * N_GRID);
    cblas_daxpy(N_GRID, dt * RK[0][0], dudt0, 1, u1, 1);

    ddt(u1, dudt1);
    memmove(u2, u1, sizeof(double) * N_GRID);
    cblas_daxpy(N_GRID, dt * RK[1][0], dudt0, 1, u2, 1);
    cblas_daxpy(N_GRID, dt * RK[1][1], dudt1, 1, u2, 1);

    // Tangent
    memmove(v, v0, sizeof(double) * N_GRID);
    ddtTan(u0, v, dvdt0, inhomo);
    cblas_daxpy(N_GRID, dt * RK[0][0], dvdt0, 1, v, 1);

    ddtTan(u1, v, dvdt1, inhomo);
    cblas_daxpy(N_GRID, dt * RK[1][0], dvdt0, 1, v, 1); 
    cblas_daxpy(N_GRID, dt * RK[1][1], dvdt1, 1, v, 1);

    ddtTan(u2, v, dvdt2, inhomo);
    cblas_daxpy(N_GRID, dt * RK[2][1], dvdt1, 1, v, 1);
    cblas_daxpy(N_GRID, dt * RK[2][2], dvdt2, 1, v, 1);

}

// Adjoint of stepPrimal, forced with strength * v0 as rhs
// Note that the rhs is added per time-step instead of per rk-step for
// consistency with discrete objective function evaluation
void
stepAdjoint(const double * u0, const double * w0, double * w, double dt, int inhomo)
{
    double u1[N_GRID], u2[N_GRID], dudt0[N_GRID], dudt1[N_GRID];
    double dwdt_u2w3[N_GRID], dwdt_u1w2[N_GRID], dwdt_u1w3[N_GRID],
           dwdt_u0w1[N_GRID], dwdt_u0w2[N_GRID];

    // Primal
    ddt(u0, dudt0);
    memmove(u1, u0, sizeof(double) * N_GRID);
    cblas_daxpy(N_GRID, dt * RK[0][0], dudt0, 1, u1, 1);

    ddt(u1, dudt1);
    memmove(u2, u1, sizeof(double) * N_GRID);
    cblas_daxpy(N_GRID, dt * RK[1][0], dudt0, 1, u2, 1);
    cblas_daxpy(N_GRID, dt * RK[1][1], dudt1, 1, u2, 1);

    // Adjoint -- w is w3
    memmove(w, w0, sizeof(double) * N_GRID);
    ddtAdj(u2, w, dwdt_u2w3);
    ddtAdj(u1, w, dwdt_u1w3);
    cblas_daxpy(N_GRID, -dt * RK[2][2], dwdt_u2w3, 1, w, 1);

    // w is now w2
    ddtAdj(u1, w, dwdt_u1w2);
    ddtAdj(u0, w, dwdt_u0w2);
    cblas_daxpy(N_GRID, -dt * RK[2][1], dwdt_u1w3, 1, w, 1); 
    cblas_daxpy(N_GRID, -dt * RK[1][1], dwdt_u1w2, 1, w, 1);

    // w is now w1
    ddtAdj(u0, w, dwdt_u0w1);
    cblas_daxpy(N_GRID, -dt * RK[1][0], dwdt_u0w2, 1, w, 1);
    cblas_daxpy(N_GRID, -dt * RK[0][0], dwdt_u0w1, 1, w, 1);

	
	// Adjoint -- source term from forcing
	double dJdu[N_GRID];
	ddObj(u0, dJdu, OBJ_IND);
    cblas_daxpy(N_GRID, inhomo * dt / T_TOTAL , dJdu, 1, w, 1);	
}


// projects an in-place vector v to orthogonal direction of dudt at
// the i_chunks'th time chunk and i_step'th time step
double
project_ddt(int i_chunk, int i_step, double * v)
{
    assert (i_chunk >= 0 && i_chunk <= N_CHUNK);
    assert (i_step >= 0 && i_step <= N_STEP);
    assert (i_step == 0 || i_chunk < N_CHUNK);
    if (PROJ)
    {
        double dudt[N_GRID];
        ddt(SOLN_U[i_chunk][i_step], dudt);
    
        int i = (i_chunk == 0 && i_step == 0) ? 0 : i_step - 1;
        for (int j = 0; j < N_GRID; ++ j)
            dudt[j] = (SOLN_U[i_chunk][i+1][j] - SOLN_U[i_chunk][i][j]) / DT_STEP;
    
        double vDotUt = cblas_ddot(N_GRID, v, 1, dudt, 1);
        double utDotUt = cblas_ddot(N_GRID, dudt, 1, dudt, 1);
    
        cblas_daxpy(N_GRID, -vDotUt / utDotUt, dudt, 1, v, 1);
        return vDotUt / utDotUt;
    } else {
        return 0.0;
    }
}


// It's beneficial to make sure that SOLN_U and SOLN_V
// are both in a contiguous chunk of physical memory
void
alloc_space_for_big_arrays_U_and_V()
{
    MEM = 1;
    // allocate space for big arrays
    SOLN_U = (double ***) malloc(sizeof(double **) * (N_CHUNK + 1));
    SOLN_V = (double ***) malloc(sizeof(double **) * N_CHUNK);
    SOLN_W = (double ***) malloc(sizeof(double **) * N_CHUNK);
	assert (SOLN_U != 0 && SOLN_V != 0);
	assert (SOLN_W != 0);

    SOLN_U[0] = (double **) malloc(sizeof(double*) * (N_CHUNK*N_STEP + 1));
    SOLN_V[0] = (double **) malloc(sizeof(double*) * N_CHUNK*N_STEP);
    SOLN_W[0] = (double **) malloc(sizeof(double*) * N_CHUNK*N_STEP);
	assert (SOLN_U[0] != 0 && SOLN_V[0] != 0);
	assert (SOLN_W[0] != 0);

    SOLN_U[0][0] = (double *) malloc(sizeof(double) * 
                                            (N_GRID + N_CHUNK*N_STEP*N_GRID));
    SOLN_V[0][0] = (double *) malloc(sizeof(double) * N_CHUNK*N_STEP*N_GRID);
	SOLN_W[0][0] = (double *) malloc(sizeof(double) * N_CHUNK*N_STEP*N_GRID);
    assert (SOLN_U[0][0] != 0 && SOLN_V[0][0] != 0);
	assert (SOLN_W[0][0] != 0);

    for (int i = 0; i < N_CHUNK; ++ i)
    {
        SOLN_U[i] = SOLN_U[0] + i * N_STEP;
        SOLN_V[i] = SOLN_V[0] + i * N_STEP;
		SOLN_W[i] = SOLN_W[0] + i * N_STEP;

        for (int j = 0; j < N_STEP; ++ j)
        {
            SOLN_U[i][j] = SOLN_U[0][0] + (i * N_STEP + j) * N_GRID;
            SOLN_V[i][j] = SOLN_V[0][0] + (i * N_STEP + j) * N_GRID;
			SOLN_W[i][j] = SOLN_W[0][0] + (i * N_STEP + j) * N_GRID;
        }
    }

    // This stores the last dangling solution
    SOLN_U[N_CHUNK] = SOLN_U[N_CHUNK - 1] + N_STEP;
    SOLN_U[N_CHUNK][0] = SOLN_U[N_CHUNK - 1][N_STEP - 1] + N_GRID;
}

// s field
void
s_field()
{
    S_FIELD = (double *) malloc(sizeof(double) * N_GRID);
    
    double dx = 128. / (N_GRID + 1);
    double x = 0.0;
    for (int i = 0; i < N_GRID; ++i)
    {
        x = ((i + 1) * dx - 64.) / 64.;
        S_FIELD[i] = S_CONST[0]*1.0 
             + S_CONST[1]*x 
             + S_CONST[2]*0.5*(3.0*pow(x,2) - 1.0)
             + S_CONST[3]*0.5*(5.0*pow(x,3) - 3*x)
             + S_CONST[4]*0.125*(35.0*pow(x,4) - 30*pow(x,2) + 3.0);
        //x = ((i + 1) * dx * 3.14159265359) / 64.;
        //S_FIELD[i] = S_CONST[0]*1.0
        //    + S_CONST[1]*cos(x)
        //    + S_CONST[2]*cos(2*x)
        //    + S_CONST[3]*cos(3*x)
        //    + S_CONST[4]*cos(4*x);
    }
}

void
run_up_to_T0(double * u, double T0, double dt_max)
{
    assert(T0 >= 0);
    assert(dt_max > 0);
    int n0_steps = (int) ceil(T0 / dt_max);
    double dt0 = T0 / n0_steps;

    assert(u != 0);
    for (int i = 0; i < n0_steps; ++ i)
    {
        stepPrimal(u, u, dt0);
    }
}


// This function initializes this module, must be called from Python
// before using any other functionality.
void
init(double * s, int n_grid, int n_mode, int proj, int n_proc)
{
    S_CONST[0] = s[0];
    S_CONST[1] = s[1];
    S_CONST[2] = s[2];
    S_CONST[3] = s[3];
    S_CONST[4] = s[4];

    assert(n_grid > 0);
    N_GRID = n_grid;

    assert(n_mode > 0);
    N_MODE = n_mode;

    PROJ = proj;

    s_field();

}

// Solve n_tan homogeneous + 1 inhomogenous tangent equations simultaneously
void
tangents(double * u, double * v, double * vc, 
        double T0,double t_chunk, double dt_max, int grad_ind)
{

    assert(t_chunk > 0);
    N_STEP = (int) ceil(t_chunk / dt_max);
    DT_STEP = t_chunk / N_STEP;
    
    GRAD_IND = grad_ind;

    run_up_to_T0(u, T0, dt_max);

    double vj[N_GRID]; 

    for (int i = 0; i < N_STEP; ++ i)
    {
        // time steps
 
        // homogeneous tangents
        for (int j = 0; j < N_MODE; ++ j)
        {
            memmove(vj, v + j * N_GRID, sizeof(double) * N_GRID); 
            stepTangent(u, vj, vj, DT_STEP, 0);
            memmove(v + j * N_GRID, vj, sizeof(double) * N_GRID);
        }
       
        // inhomogeneous tangent
        stepTangent(u, vc, vc, DT_STEP, 1);
        stepPrimal(u, u, DT_STEP);
    }

}    

// Solve n_tan homogeneous + 1 inhomogenous tangent equations simultaneously
// generate integrated quantities needed for tangent nilss with checkpoint minimization

void
nilss_tan(double * u, double * v, double *vc, double * g, int n_chunk, double t_chunk, double dt_max)
{
    SHAD_PROJ = 1; // flag to indicate that shad_proj_tan or shad_proj_adj has been run
    
    assert(n_chunk > 0);
    N_CHUNK = n_chunk;
    assert(t_chunk > 0);
    N_STEP = (int) ceil(t_chunk / dt_max);
    DT_STEP = t_chunk / N_STEP;

    T_TOTAL = N_CHUNK * N_STEP * DT_STEP; 

    double Jbar[5] = { 0.0 };

    double vj[N_GRID]; 
	double dJdu[N_GRID];

    for (int i = 0; i < N_STEP; ++ i)
    {
       
        // time steps
 
        // homogeneous tangents
        for (int j = 0; j < N_MODE; ++ j)
        {
            memmove(vj, v + j * N_GRID, sizeof(double) * N_GRID); 
            stepTangent(u, vj, vj, DT_STEP, 0);
            memmove(v + j * N_GRID, vj, sizeof(double) * N_GRID);
        }
       
        // inhomogeneous tangent
        stepTangent(u, vc, vc, DT_STEP, 1);
        
        stepPrimal(u, u, DT_STEP);

        // integrate W^T dJdu
        ddObj(u, dJdu, OBJ_IND);
        for (int j = 0; j < N_MODE; ++ j)
        {
            memmove(vj, v + j * N_GRID, sizeof(double) * N_GRID);
            g[j] = g[j] +  (DT_STEP/T_TOTAL) * cblas_ddot(N_GRID, vj, 1, dJdu, 1);
        }


        // objective function
        for (int j = 0; j < 5; ++ j)
        {
            Jbar[j] = Jbar[j] + (DT_STEP/T_TOTAL) * Obj(u,j);
        }
    }

    // add to total average:
    for (int i = 0; i < 5; ++ i)
    {
	    JBAR[i] = JBAR[i] + Jbar[i];
    }
    
    // time dilation
    if (PROJ)
        {
           double dudt[N_GRID];
           ddt(u , dudt);
           double vDotUt = 0.0;
           double utDotUt = cblas_ddot(N_GRID, dudt, 1, dudt, 1);

           // time dilation for homogeneous tangent
           for (int j = 0; j < N_MODE; ++ j)
               {
                   memmove(vj, v + j * N_GRID, sizeof(double) * N_GRID); 
                   vDotUt = cblas_ddot(N_GRID, vj, 1, dudt, 1);
                   cblas_daxpy(N_GRID, -vDotUt / utDotUt, dudt, 1, vj, 1);
                   memmove(v + j * N_GRID, vj, sizeof(double) * N_GRID);
               }

           // time dilation for inhomogeneous tangent
           vDotUt = cblas_ddot(N_GRID, vc, 1, dudt, 1);
           cblas_daxpy(N_GRID, -vDotUt / utDotUt, dudt, 1, vc, 1);

        } 


}

// Solve n_tan homogeneous + 1 inhomogenous tangent equations simultaneously
// generate integrated quantities needed for adjoint nilss with checkpoint minimization

double
nilss_adj(double * u, double * v, double * g1, double * g2, int n_chunk, double t_chunk, double dt_max)
{
    SHAD_PROJ = 1; // flag to indicate that shad_proj_tan or shad_proj_adj has been run
    
    assert(n_chunk > 0);
    N_CHUNK = n_chunk;
    assert(t_chunk > 0);
    N_STEP = (int) ceil(t_chunk / dt_max);
    DT_STEP = t_chunk / N_STEP;

    T_TOTAL = N_CHUNK * N_STEP * DT_STEP; 

    double Jbar[5] = { 0.0 };

    double vj[N_GRID]; 
	double dJdu[N_GRID];

    for (int i = 0; i < N_STEP; ++ i)
    {
       
        // time steps
 
        // homogeneous tangents
        for (int j = 0; j < N_MODE; ++ j)
        {
            memmove(vj, v + j * N_GRID, sizeof(double) * N_GRID); 
            stepTangent(u, vj, vj, DT_STEP, 0);
            memmove(v + j * N_GRID, vj, sizeof(double) * N_GRID);
        }
       
        stepPrimal(u, u, DT_STEP);

        // integrate W^T dJdu
        ddObj(u, dJdu, OBJ_IND);
        for (int j = 0; j < N_MODE; ++ j)
        {
            memmove(vj, v + j * N_GRID, sizeof(double) * N_GRID);
            g1[j] += cblas_ddot(N_GRID, vj, 1, dJdu, 1) * (DT_STEP/T_TOTAL);
        }

        // objective function
        for (int j = 0; j < 5; ++ j)
        {
            Jbar[j] = Jbar[j] + (DT_STEP/T_TOTAL) * Obj(u,j);
        }
    }

    // add to total average:
    for (int i = 0; i < 5; ++ i)
    {
	    JBAR[i] = JBAR[i] + Jbar[i];
    }

    // time dilation
    if (PROJ)
        {
           double dudt[N_GRID];
           ddt(u , dudt);
           double vDotUt = 0.0;
           double utDotUt = cblas_ddot(N_GRID, dudt, 1, dudt, 1);

           // time dilation for homogeneous tangent
           for (int j = 0; j < N_MODE; ++ j)
               {
                   memmove(vj, v + j * N_GRID, sizeof(double) * N_GRID); 
                   vDotUt = cblas_ddot(N_GRID, vj, 1, dudt, 1);
                   g2[j] = vDotUt / utDotUt;
                   cblas_daxpy(N_GRID, -vDotUt / utDotUt, dudt, 1, vj, 1);
                   memmove(v + j * N_GRID, vj, sizeof(double) * N_GRID);
               }

        } 
    return Obj(u,OBJ_IND);
}

// Assign value to global variable JBAR, the long time averaged objective function J
void
assignJBAR(double * Jbar)
{
    for (int i = 0; i < 5; ++ i)
    {
	    JBAR[i] = Jbar[i];
    }
}

// Solve for primal over time horizon initiated in shad_proj_tan or shad_proj_adj
// Allocate memory for tangent and adjoint solutions
// must be run AFTER shad_proj_tan or shad_proj_adj
void
generate_solution(double * u0)
{
    //assert (SHAD_PROJ == 1);

    alloc_space_for_big_arrays_U_and_V();

    memmove(SOLN_U[0][0], u0, sizeof(double) * N_GRID);

    // run in each time chunk, compute time averaged objective fcn Jbar
	double Jbar[5] = { 0.0 };
    for (int i_chunk = 0; i_chunk < N_CHUNK; ++ i_chunk)
    {
        double ** u = SOLN_U[i_chunk];
        for (int i = 0; i < N_STEP; ++ i)
        {
            // when i == N_STEP - 1, u[i+1] is SOLN_U{i_chunk + 1][0]
            // there's one more memory slot at the last i_chunk, so don't worry
            stepPrimal(u[i], u[i+1], DT_STEP);
		    for (int j = 0; j < 5; ++ j)
            {
                Jbar[j] = Jbar[j] + (DT_STEP/T_TOTAL) * Obj(u[i],j);
            }
        }
    }
	assignJBAR(Jbar);
    memmove(u0, SOLN_U[N_CHUNK][0], sizeof(double) * N_GRID);
}

// Solve the tangent equation in the i_chunk'th chunk with v0 as the
// initial condition.  This can be called from Python.
// The solution at the end of the time chunk then overwrites v0.
// See ddtTan for the argument inhomo.
void
tangent(int i_chunk, double * v0, int inhomo)
{
    assert (SHAD_PROJ == 1);
    assert (MEM == 1);
    assert (i_chunk >= 0 && i_chunk < N_CHUNK);
    double ** u = SOLN_U[i_chunk];
    double ** v = SOLN_V[i_chunk];
	double dx = 128. / (N_GRID + 1);
    memmove(v[0], v0, sizeof(double) * N_GRID);

    for (int i = 0; i < N_STEP - 1; ++ i)
    {
        stepTangent(u[i], v[i], v[i+1], DT_STEP, inhomo);
    }
    stepTangent(u[N_STEP - 1], v[N_STEP - 1], v0, DT_STEP, inhomo);

}


// Solve the adjoint equation in the i_chunk'th chunk with w0 as the
// terminal condiiton.  This can be called from Python.
// The solution at the beginning of the time chunk then overwrites v0.
// See stepAdjoint for the argument forcing.
void
adjoint(int i_chunk, double * w0, double forcing, int inhomo)
{
    assert (SHAD_PROJ == 1);
    assert (MEM == 1);

    assert (i_chunk >= 0 && i_chunk < N_CHUNK);
    double ** u = SOLN_U[i_chunk];
	double ** w = SOLN_W[i_chunk];


	// Terminal condition
	double dudtT[N_GRID];
	ddt(u[N_STEP],dudtT);  
    double utDotUt = cblas_ddot(N_GRID, dudtT, 1, dudtT, 1);
    //TODO: uncomment for projections!
    if (PROJ)
        {
        cblas_daxpy(N_GRID, inhomo * ((JBAR[OBJ_IND] - Obj(u[N_STEP],OBJ_IND)) / (utDotUt * T_TOTAL)), dudtT, 1, w0, 1);
        }

	stepAdjoint(u[N_STEP - 1], w0, w[N_STEP - 1], DT_STEP, inhomo);
    for (int i = N_STEP - 2; i >= 0; -- i)
    {
        stepAdjoint(u[i], w[i+1], w[i], DT_STEP, inhomo);
    }
	stepAdjoint(u[0], w[1], w0, DT_STEP, inhomo);

/*	for (int i = N_STEP - 1; i >= 0; -- i)
    {
        stepAdjoint(u[i], v[i], 1, w0, w0, DT_STEP);
    }
*/

}

// Compute the gradient with respect to s (Tangent LSS)
void
grad(double * dJbar_ds)
{

    double dJdu[N_GRID];


    for (int i = 0; i < 5; ++i)
    {
        dJbar_ds[i] = 0;
    }

    for (int i_chunk = 0; i_chunk < N_CHUNK; ++ i_chunk)
    {
        double ** v = SOLN_V[i_chunk];
        double ** u = SOLN_U[i_chunk];
        for (int j = 0; j < 5; ++ j)
        {
            for (int i = 0; i < N_STEP; ++ i)
            {
                ddObj(u[i], dJdu, j);
                dJbar_ds[j] = dJbar_ds[j] + (DT_STEP/T_TOTAL) * 
                     cblas_ddot(N_GRID, dJdu, 1, v[i], 1);
            }
            if (PROJ)
                {
    		    dJbar_ds[j] = dJbar_ds[j] + project_ddt(i_chunk, N_STEP, v[N_STEP-1]) 
                 * (JBAR[j] - Obj(u[N_STEP],j)) / T_TOTAL;
                }
            
        }
    }

}



// Compute the gradient with respect to s (Adjoint LSS)
void 
gradAdj(double * dJbar_ds)
{

	double dfds[N_GRID];

    for (int i = 0; i < 5; ++i)
    {
        dJbar_ds[i] = 0;
    }


    for (int i_chunk = 0; i_chunk < N_CHUNK; ++ i_chunk)
    {
        double ** w = SOLN_W[i_chunk];
		double ** u = SOLN_U[i_chunk];
        for (int i = 0; i < N_STEP; ++ i)
        {
            for (int j = 0; j < 5; ++j)
            {    
                ddtds(u[i], dfds, j);
			    dJbar_ds[j] = dJbar_ds[j] + DT_STEP * cblas_ddot(N_GRID, dfds, 1, w[i], 1);
            }
        }
    }

}


// This function runs the primal and computes time averaged statistics over the time horizon [T0, T0 + T1]. The statistics are stored in global variable JBAR and can be accessed in python using pde.c_Jbar 
void
pdeTavg(double * s, double * u0, int n_grid, double T0, double T1, double dt_max)
{
    S_CONST[0] = s[0];
    S_CONST[1] = s[1];
    S_CONST[2] = s[2];
    S_CONST[3] = s[3];
    S_CONST[4] = s[4];

    assert(n_grid > 0);
    N_GRID = n_grid;

    assert(T1 > 0);
    N_STEP = (int) ceil(T1 / dt_max);
    DT_STEP = T1 / N_STEP;

	T_TOTAL = T1;

    s_field();

    run_up_to_T0(u0, T0, dt_max);

    // run in each time chunk, compute time averaged objective fcn Jbar
	double ubar[N_GRID];
    double Jbar[5] = { 0.0 };

    for (int k = 0; k < N_GRID; ++k)
    {
        ubar[k] = 0.0;
    }


    for (int i = 0; i < N_STEP; ++ i)
    {
	    for (int j = 0; j < 5; ++ j)
        {
            Jbar[j] = Jbar[j] + (DT_STEP/T_TOTAL) * Obj(u0,j);
        }
        stepPrimal(u0, u0, DT_STEP);
        
        for (int k = 0; k < N_GRID; ++k)
        {
            ubar[k] = ubar[k] + (DT_STEP/T_TOTAL) * u0[k];
        }
    }
	assignJBAR(Jbar);
    
    
    for (int k = 0; k < N_GRID; ++k)
    {
        u0[k] = ubar[k];
    }
}


