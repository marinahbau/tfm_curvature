/// \file
#include <iostream>
#include <igl/opengl/glfw/Viewer.h>

 #include <igl/gaussian_curvature.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>
#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>

int main(int argc, char *argv[])
{
  //Leer fichero 
  std::string filename = (argc > 1) ? argv[1] : "bunny.off"; //Malla por defecto

  using namespace Eigen;
  using namespace std;
  MatrixXd V;
  MatrixXi F;
  igl::readOFF(filename,V,F);

  VectorXd K;
  // Compute integral of Gaussian curvature
  igl::gaussian_curvature(V,F,K);
  // Compute mass matrix
  SparseMatrix<double> M,Minv;
  igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_DEFAULT,M);
  igl::invert_diag(M,Minv);
  // Divide by area to get integral average
  K = (Minv*K).eval();
   

  //Crear archivo de valores de curvatura
  ofstream s("ficheros_curvaturas/curvatura.txt");

  // Loop over K
  for(int i = 0;i<(int)K.rows();i++)
  {
    for(int j = 0;j<(int)K.cols();++j)
    {
      s << K(i,j) << ", ";
    }
  }
  s.close();


  // Plot the mesh with pseudocolors
  /*igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V, F);
  viewer.data().set_data(K);
  viewer.launch();*/
}
