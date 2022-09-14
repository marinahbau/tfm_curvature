/// \file
#include <iostream>
#include <igl/opengl/glfw/Viewer.h>
 #include <igl/gaussian_curvature.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>
#include <igl/readOFF.h>
#include <igl/avg_edge_length.h>
#include <igl/cotmatrix.h>
#include <igl/parula.h>
#include <igl/per_corner_normals.h>
#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/principal_curvature.h>
#include <igl/read_triangle_mesh.h>

Eigen::MatrixXd V;
Eigen::MatrixXi F;

std::string LOCALNAME = "";
int main(int argc, char *argv[])
{
   using namespace Eigen;
  std::string filename = (argc > 1) ? argv[1] : "../tfm/data/SFcoBorja.ply"; //Malla por defecto
  int variacion = std::atoi(argv[2]);
  LOCALNAME = filename;
  LOCALNAME.erase(LOCALNAME.begin(), LOCALNAME.begin() + variacion); //Para ejecutar solo este archivo hay que poner +12
  LOCALNAME.erase(LOCALNAME.end() - 4, LOCALNAME.end());
  // Load a mesh in OFF format
  igl::read_triangle_mesh(filename, V, F);
  if(argc>1)
  {
    filename = argv[1];
  }

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
  std::ofstream s("../tfm/curvaturas/" + LOCALNAME + "_gaussian_curvatura.txt");

  // Loop over K
  for(int i = 0;i<(int)K.rows();i++)
  {
    for(int j = 0;j<(int)K.cols();++j)
    {
      s << K(i,j) << ",";
    }
  }
  s.close();

}
