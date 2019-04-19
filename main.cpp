#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>
#include <igl/arap.h>
#include <igl/boundary_loop.h>
#include <igl/harmonic.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/readOFF.h>
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include<map>
#include<algorithm>
#include "tutorial_shared_path.h"
using namespace std;
bool readFeaturePoints(const string fpname, igl::opengl::ViewerData &data);
bool findSymmetryAxis(igl::opengl::ViewerData &data);
bool dealRotate(igl::opengl::ViewerData &data);
bool dealEdit(igl::opengl::ViewerData &data);
double calSEA(igl::opengl::ViewerData &data, const string fpname);
double calSEB(igl::opengl::ViewerData &data, const string fpname);
bool SEAAdjust(const string fpname);
bool SEBAdjust(const string fpname);
bool makeFeatureChange(igl::opengl::ViewerData &data, igl::ARAPData &arap_data);
bool show_uv = false;
map<int, int> SymPairs;
//Edit Vars
int dealingindex = 0;
enum Orientation { Up = 0, Down, Left, Right };
Orientation dir = Up;
bool boolSymEdit = true;
Eigen::MatrixXd V;
Eigen::MatrixXi F;
Eigen::MatrixXd V_uv;
Eigen::MatrixXd UV_Adjust;
Eigen::MatrixXd initial_guess;
Eigen::MatrixXd AxisData(2, 1);
int main(int argc, char *argv[])
{
  // Load a mesh in OFF format
  igl::readOFF(TUTORIAL_SHARED_PATH "/meanshape.off", V, F);
  igl::opengl::glfw::Viewer viewer;
  readFeaturePoints(TUTORIAL_SHARED_PATH "/MPEG4_FDP_face05.fp", viewer.data()); //Reading feature points

  Eigen::VectorXi bnd;
  igl::boundary_loop(F, bnd);
  Eigen::MatrixXd bnd_uv;
  igl::map_vertices_to_circle(V, bnd, bnd_uv);
  igl::harmonic(V, F, bnd, bnd_uv, 1, initial_guess);//First map the mesh to circle with harmonic
													 // Add dynamic regularization to avoid to specify boundary conditions
  igl::ARAPData arap_data;
  arap_data.with_dynamics = false;
  Eigen::VectorXi b = Eigen::VectorXi::Zero(0);
  Eigen::MatrixXd bc = Eigen::MatrixXd::Zero(0, 0);
  // Initialize ARAP
  arap_data.max_iter = 10;
  // 2 means that we're going to *solve* in 2d
  arap_precomputation(V, F, 2, b, arap_data);
  // Solve arap using the harmonic map as initial guess
  V_uv = initial_guess;
  arap_solve(bc, arap_data, V_uv);
  // Scale UV to make the texture more clear
  V_uv *= 20;
  // Plot the mesh
  viewer.data().set_mesh(V, F);
  viewer.data().set_uv(V_uv);
  findSymmetryAxis(viewer.data());
  dealRotate(viewer.data());
  // Attach a menu plugin
  UV_Adjust = V_uv;
  igl::opengl::glfw::imgui::ImGuiMenu menu;
  viewer.plugins.push_back(&menu);

  // Customize the menu
  double doubleVariable = 0.1f; // Shared between two menus
  // Add content to the default menu window
  menu.callback_draw_viewer_menu = [&]()
  {
    // Draw parent menu content
    menu.draw_viewer_menu();

    // Add new group
    if (ImGui::CollapsingHeader("New Group", ImGuiTreeNodeFlags_DefaultOpen))
    {
      // Expose variable directly ...
      ImGui::InputDouble("double", &doubleVariable, 0, 0, "%.4f");

      // ... or using a custom callback
      static bool boolShowUV = false;
      if (ImGui::Checkbox("Show Parametric results", &boolShowUV))
      {
        // do something
		  if (boolShowUV)
		  {
			  viewer.data().set_mesh(V_uv, F);
			  viewer.core.align_camera_center(V_uv, F);
		  }
		  else
		  {
			  viewer.data().set_mesh(V, F);
			  viewer.core.align_camera_center(V, F);
		  }
		  viewer.data().compute_normals();
      }
      // Expose an enumeration type
      ImGui::Combo("Direction", (int *)(&dir), "Up\0Down\0Left\0Right\0\0");
	  if (ImGui::Checkbox("Symmetrical Edit", &boolSymEdit))
		  cout << boolSymEdit << endl;
	  if (ImGui::InputInt("Edit Index", &dealingindex))
		  cout << dealingindex << endl;
	  if (ImGui::Button("Deal Edit", ImVec2(-1, 0)))
	  {
		  //makeFeatureChange(viewer.data(), arap_data);
		  dealRotate(viewer.data());
	  }
      // We can also use a vector<string> defined dynamically
	  /*
	  static int num_choices = 3;
      static vector<string> choices;
      static int idx_choice = 0;
      if (ImGui::InputInt("Num letters", &num_choices))
      {
        num_choices = max(1, min(26, num_choices));
      }
      if (num_choices != (int) choices.size())
      {
        choices.resize(num_choices);
        for (int i = 0; i < num_choices; ++i)
          choices[i] = string(1, 'A' + i);
        if (idx_choice >= num_choices)
          idx_choice = num_choices - 1;
      }
      ImGui::Combo("Letter", &idx_choice, choices);
	  */
      // Add a button
	  /*
	  if (ImGui::Button("Testing Axis", ImVec2(-1,0)))
      {
		  findSymmetryAxis(viewer.data());
      }
	  if (ImGui::Button("Testing Rotate", ImVec2(-1, 0)))
	  {
		  //makeFeatureChange(viewer.data(), arap_data);
		  dealRotate(viewer.data());
	  }
	  */
	  if (ImGui::Button("Calculate SE A", ImVec2(-1, 0)))
		  //Axis points
		  calSEA(viewer.data(), TUTORIAL_SHARED_PATH "/SymPoints.SP");
	  if (ImGui::Button("Calculate SE B", ImVec2(-1, 0)))
		  //Axis points
		  calSEB(viewer.data(), TUTORIAL_SHARED_PATH "/SymPoints.SP");
	  if (ImGui::Button("Adjust SE A", ImVec2(-1, 0)))
		  SEAAdjust(TUTORIAL_SHARED_PATH "/SymPoints.SP");
	  if (ImGui::Button("Adjust SE B", ImVec2(-1, 0)))
		  SEBAdjust(TUTORIAL_SHARED_PATH "/SymPoints.SP");
	  if (ImGui::Button("ARAP Solve", ImVec2(-1, 0)))
	  {
		  makeFeatureChange(viewer.data(), arap_data);
	  }
    }
  };
  viewer.launch();
}
bool findSymmetryAxis(igl::opengl::ViewerData &data)
{
	//Symmetry Axis point : 40777 8309 8319 8334 8344 8366 8354 8374
	Eigen::MatrixXd x(8, 2);
	for (int i = 0; i<8; ++i) x(i, 1) = 1;
	x(0, 0) = data.V_uv(40777, 0);
	x(1, 0) = data.V_uv(8309, 0);
	x(2, 0) = data.V_uv(8319, 0);
	x(3, 0) = data.V_uv(8334, 0);
	x(4, 0) = data.V_uv(8344, 0);
	x(5, 0) = data.V_uv(8366, 0);
	x(6, 0) = data.V_uv(8354, 0);
	x(7, 0) = data.V_uv(8374, 0);
	Eigen::MatrixXd xt = x.transpose();
	Eigen::MatrixXd y(8, 1);
	y(0, 0) = data.V_uv(40777, 1);
	y(1, 0) = data.V_uv(8309, 1);
	y(2, 0) = data.V_uv(8319, 1);
	y(3, 0) = data.V_uv(8334, 1);
	y(4, 0) = data.V_uv(8344, 1);
	y(5, 0) = data.V_uv(8366, 1);
	y(6, 0) = data.V_uv(8354, 1);
	y(7, 0) = data.V_uv(8374, 1);
	Eigen::MatrixXd b(2, 1);
	b = (xt * x).inverse() * xt * y; //主要公式
	cout<< "Y = " << b(0, 0) << "X +" << b(1, 0) << endl;
	cout << data.V_uv.row(40777) << endl;
	cout << data.V_uv.row(8309) << endl;
	cout << data.V_uv.row(8319) << endl;
	cout << data.V_uv.row(8334) << endl;
	cout << data.V_uv.row(8344) << endl;
	cout << data.V_uv.row(8366) << endl;
	cout << data.V_uv.row(8354) << endl;
	cout << data.V_uv.row(8374) << endl;
	//Eigen::MatrixXd P1(1, 3);
	//P1(0, 0) = data.V_uv.colwise().minCoeff()(0); P1(0, 1) = data.V_uv.colwise().minCoeff()(0)*b(0, 0) + b(1, 0); P1(0, 2) = 0;
	//Eigen::MatrixXd P2(1, 3);
	//P2(0, 0) = data.V_uv.colwise().maxCoeff()(1); P2(0, 1) = data.V_uv.colwise().maxCoeff()(1)*b(0, 0) + b(1, 0); P2(0, 2) = 0;
	//data.add_edges(P1,P2,Eigen::RowVector3d(1, 0, 0));
	AxisData = b;
	return true;
}

bool dealRotate(igl::opengl::ViewerData &data)
{
	double sinA = sqrt(pow(AxisData(0, 0),2) + 1)/ AxisData(0, 0);
	double cosA = sqrt(pow(AxisData(0, 0), 2) + 1);
	//cout << sinA << endl;
	//cout << cosA << endl;
	//x = (x1 - x2)cosθ - (y1 - y2)sinθ + x2
	//y = (y1 - y2)cosθ + (x1 - x2)sinθ + y2
	//8334
	//cout << data.V_uv.row(8334) << endl;
	double x2 = data.V_uv(8334, 0);
	double y2 = data.V_uv(8334, 1);
	cout << data.V_uv(40777) << endl;
	for (int i = 0; i < data.V_uv.rows(); i++)
	{
		//cout << i << data.V_uv(i, 0) << data.V_uv(i, 1) << endl;
		double nx = (data.V_uv(i, 0) - x2)*cosA - (data.V_uv(i, 1) - y2)*sinA;
		double ny = -1*((data.V_uv(i, 1) - y2)*cosA + (data.V_uv(i, 0) - x2)*sinA) + y2;
		//cout << i << nx << ny << endl;
		data.V_uv(i, 0) = nx;
		data.V_uv(i, 1) = ny;
	}
	data.set_mesh(data.V_uv, data.F);
	V_uv = data.V_uv;
	//cout << data.V_uv(40777) << endl;
	return true;
}
bool dealEdit(igl::opengl::ViewerData &data)
{
	if (dir == Up)
	{
		V_uv(dealingindex, 1) += 1;
		data.set_uv(V_uv);
		data.set_mesh(V_uv, data.F);
		UV_Adjust = V_uv;
	}
	if (boolSymEdit == true)
	{

	}
	return true;
}
bool makeFeatureChange(igl::opengl::ViewerData &data, igl::ARAPData &arap_data)
{
	Eigen::MatrixXd V_temp = data.V_uv;
	arap_data.b = data.FPs;
	Eigen::MatrixXd bc(data.FPs.size(), data.V_uv.cols());
	for (int i = 0; i < data.FPs.rows(); i++)
	{
		bc(i, 0) = UV_Adjust(data.FPs(i), 0);
		bc(i, 1) = UV_Adjust(data.FPs(i), 1);
	}
	cout << "BC:" << endl;
	for (int i = 0; i < bc.rows(); i++)
	{
		cout << arap_data.b.row(i) << " : " << bc.row(i) << endl;
	}
	//arap_data.with_dynamics = false;
	arap_data.max_iter = 5;
	arap_precomputation(data.V_uv, F, 2, arap_data.b, arap_data);
	arap_solve(bc, arap_data, V_temp);
	data.set_uv(V_temp);
	data.set_mesh(V_temp, data.F);
	V_uv = V_temp;
	UV_Adjust = V_temp;
	cout << "Results:" << endl;
	for (int i = 0; i < bc.rows(); i++)
	{
		cout << arap_data.b(i, 0) << " : " << V_temp.row(arap_data.b(i,0)) << endl;

	}
	return true;
}
double calSEA(igl::opengl::ViewerData &data, const string fpname)
{
	string line;
	ifstream fp_file(fpname);
	stringstream ss;
	double TotalE = 0;
	if (fp_file)
	{
		int p1; int p2;
		data.FaceModel = true;//Confirm face
		data.FPs.resize(54);
		int i = 0;
		while (getline(fp_file, line)) // line中不包括每行的换行符
		{
			if (line[0] == '#') continue;
			ss << line.substr(0, line.find_first_of(' '));
			ss >> p1;
			ss.clear();
			ss << line.substr(line.find_first_of(' '), line.size() - 1);
			ss >> p2;
			ss.clear();
			double V;
			if (p2 > 0) V = (data.V_uv(p1, 1) - data.V_uv(p2, 1)) * 1;
			else V = 0;
			V *= V;
			cout << p1 << " " << p2 << " " << V << endl;
			if (p2 > 0) cout << data.V_uv.row(p1) << " " << data.V_uv.row(p2) << endl;
			else cout << data.V_uv.row(p1) << endl;
			TotalE += V;
		}
	}
	else
	{
		cout << "IOError: " << fpname << " could not be opened...\n" << endl;
		return false;
	}
	cout << TotalE << endl;
	return TotalE;
}
double calSEB(igl::opengl::ViewerData &data, const string fpname)
{
	string line;
	ifstream fp_file(fpname);
	stringstream ss;
	double TotalE = 0;
	if (fp_file)
	{
		int p1; int p2;
		int i = 0;
		while (getline(fp_file, line)) // line中不包括每行的换行符
		{
			if (line[0] == '#') continue;
			ss << line.substr(0, line.find_first_of(' '));
			ss >> p1;
			ss.clear();
			ss << line.substr(line.find_first_of(' '), line.size()-1);
			ss >> p2;
			ss.clear();
			double V;
			if (p2 > 0) V = data.V_uv(p1, 0) + data.V_uv(p2, 0);
			else V = data.V_uv(p1, 0);
			V *= V;
			cout << p1 << " " << p2 << " " << V << endl;
			TotalE += V;
		}
	}
	else
	{
		cout << "IOError: " << fpname << " could not be opened...\n" << endl;
		return false;
	}
	cout << TotalE << endl;
	return TotalE;
}
bool SEAAdjust(const string fpname)
{
	string line;
	ifstream fp_file(fpname);
	stringstream ss;
	if (fp_file)
	{
		int p1; int p2;
		int i = 0;
		while (getline(fp_file, line)) // line中不包括每行的换行符
		{
			if (line[0] == '#') continue;
			ss << line.substr(0, line.find_first_of(' '));
			ss >> p1;
			ss.clear();
			ss << line.substr(line.find_first_of(' '), line.size() - 1);
			ss >> p2;
			ss.clear();
			if (p2 > 0)
			{
				cout << UV_Adjust.row(p1) << " " << UV_Adjust.row(p2) << endl;
				double k = (UV_Adjust(p2, 1) - UV_Adjust(p1, 1)) / (UV_Adjust(p2, 0) - UV_Adjust(p1, 0)); //(y2-y1)/(x2-x1)
				double b = UV_Adjust(p2, 1) - UV_Adjust(p2, 0) * k;
				UV_Adjust(p1, 1) = b;
				UV_Adjust(p2, 1) = b;
				cout << p1 << " " << p2 << " " << k << " " << b << endl;
				cout << UV_Adjust.row(p1) << " " << UV_Adjust.row(p2) << endl;
			}
		}
	}
	else
	{
		cout << "IOError: " << fpname << " could not be opened...\n" << endl;
		return false;
	}
	return true;
}
bool SEBAdjust(const string fpname)
{
	string line;
	ifstream fp_file(fpname);
	stringstream ss;
	if (fp_file)
	{
		int p1; int p2;
		int i = 0;
		while (getline(fp_file, line)) // line中不包括每行的换行符
		{
			if (line[0] == '#') continue;
			ss << line.substr(0, line.find_first_of(' '));
			ss >> p1;
			ss.clear();
			ss << line.substr(line.find_first_of(' '), line.size() - 1);
			ss >> p2;
			ss.clear();
			cout << p1 << " " << p2 << endl;
			if (p2 > 0)
			{
				cout << UV_Adjust.row(p1) << " " << UV_Adjust.row(p2) << endl;
				if (UV_Adjust(p1, 0) > 0)
				{
					double ax = UV_Adjust(p1, 0) - UV_Adjust(p2, 0);
					ax /= 2;
					UV_Adjust(p1, 0) = ax;
					UV_Adjust(p2, 0) = -1 * ax;
				}
				else
				{
					double ax = UV_Adjust(p2, 0) - UV_Adjust(p1, 0);
					ax /= 2;
					UV_Adjust(p1, 0) = -1 * ax;
					UV_Adjust(p2, 0) = ax;
				}
				cout << UV_Adjust.row(p1) << " " << UV_Adjust.row(p2) << endl;
			}
			else
			{
				cout << UV_Adjust.row(p1) << endl;
				UV_Adjust(p1, 0) = 0;
				cout << UV_Adjust.row(p1) << endl;
			}
			
		}
	}
	else
	{
		cout << "IOError: " << fpname << " could not be opened...\n" << endl;
		return false;
	}
	return true;
}
bool readFeaturePoints(const string fpname, igl::opengl::ViewerData &data)
{
	string line;
	ifstream fp_file(fpname);
	stringstream ss;
	if (fp_file)
	{
		int fpindex;
		data.FaceModel = true;//Confirm face
		data.FPs.resize(54);
		int i = 0;
		while (getline(fp_file, line)) // line中不包括每行的换行符
		{
			if (line[0] == '#') continue;
			ss << line.substr(0, line.find_first_of(' '));
			ss >> fpindex;
			ss.clear();
			data.FPs(i++) = fpindex;
		}
	}
	else
	{
		cout << "IOError: " << fpname << " could not be opened...\n" << endl;
		return false;
	}
}
