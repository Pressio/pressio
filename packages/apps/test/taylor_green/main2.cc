
// #include "taygreen.hpp"
#include "apps_taylor_green.hpp"
#include "apps_shear_layer.hpp"


int main(int argc, char *argv[])
{
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);

  //apps::TaylorGreen2d obj(64, &Comm);
  apps::ShearLayer obj(270, &Comm);
  obj.setup();
  double dt = 0.001;
  obj.run(dt, static_cast<int>(dt/dt) );
    
  MPI_Finalize();
  return 0;
}
