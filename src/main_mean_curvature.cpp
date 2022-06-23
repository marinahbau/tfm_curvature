#include <igl/avg_edge_length.h>
#include <igl/cotmatrix.h>
#include <igl/invert_diag.h>
#include <igl/massmatrix.h>
#include <igl/parula.h>
#include <igl/per_corner_normals.h>
#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/principal_curvature.h>
#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>

Eigen::MatrixXd V;
Eigen::MatrixXi F;
std::string LOCALNAME = "";

int main(int argc, char *argv[])
{
  using namespace Eigen;
  std::string filename = (argc > 1) ? argv[1] : "../tfm/data/sfcoborjacienmilnormalizado.ply"; //Malla por defecto
  LOCALNAME = filename;
  LOCALNAME.erase(LOCALNAME.begin(), LOCALNAME.begin() + 12);
  LOCALNAME.erase(LOCALNAME.end() - 4, LOCALNAME.end());

  if(argc>1)
  {
    filename = argv[1];
  }

  std::cout << "Iniciando TFM curvatura para la malla " << LOCALNAME << std::endl;
  // Load a mesh in OFF format
  igl::read_triangle_mesh(filename, V, F);

  // Alternative discrete mean curvature
  MatrixXd HN;
  SparseMatrix<double> L,M,Minv;
  igl::cotmatrix(V,F,L);
  igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_VORONOI,M);
  igl::invert_diag(M,Minv);
  // Laplace-Beltrami of position
  HN = -Minv*(L*V);
  // Extract magnitude as mean curvature
  VectorXd H = HN.rowwise().norm();

  // Compute curvature directions via quadric fitting
  MatrixXd PD1,PD2;
  VectorXd PV1,PV2;
  igl::principal_curvature(V,F,PD1,PD2,PV1,PV2);
  // mean curvature
  H = 0.5*(PV1+PV2);

  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V, F);

  viewer.data().set_data(H);

  //Crear archivo de valores de curvatura
  std::ofstream s("../tfm/curvaturas/" + LOCALNAME + "_mean_curvatura.txt");

  // Loop over H
  for(int i = 0;i<(int)H.rows();i++)
  {
    for(int j = 0;j<(int)H.cols();++j)
    {
      s << H(i,j) << ",";
    }
  }
  s.close();

}