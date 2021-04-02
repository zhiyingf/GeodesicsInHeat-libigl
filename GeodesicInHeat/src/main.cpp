#include "GeodesicInHeat.h"
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/read_triangle_mesh.h>
#include <igl/png/readPNG.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>


Eigen::MatrixXd V, TC, N;
Eigen::MatrixXi F, FTC, FN;


void setTexture(std::string pngName, Eigen::VectorXd geodistance, igl::opengl::glfw::Viewer& viewer)
{
	//distance play:use texture coordinate
	TC.setZero(geodistance.size(), 2);
	Eigen::MatrixXd::Index maxRow, maxCol;
	double max = geodistance.maxCoeff(&maxRow, &maxCol);
	if (max > 1.0)
	{
		Eigen::VectorXd geodistanceNorm = geodistance * (1.0 / max);
		TC.block(0, 0, geodistance.size(), 1) = geodistanceNorm;
	}
	else
	{
		TC.block(0, 0, geodistance.size(), 1) = geodistance;
	}
	// Allocate temporary buffers
	Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> R, G, B, A;

	// Read the PNG
	igl::png::readPNG(pngName, R, G, B, A);
	viewer.data().set_uv(TC, FTC);
	viewer.data().show_texture = true;
	viewer.data().set_mesh(V, F);
	viewer.core().align_camera_center(V, F);
	viewer.data().set_texture(R, G, B);
}

int main(int argc, char* argv[])
{	
	std::string fileName = "../data/bunny.obj";//result1.obj
	igl::readOBJ(argc > 1 ? argv[1] : fileName, V, TC, N, F, FTC, FN);

	TriangleMesh triMesh(&V, &F);
	GeodesicInHeat GeoInHeat(triMesh, GeodesicInHeat::NEUMAN);

	Eigen::VectorXd geodistance;

	size_t sourceVi = 0;
	
	// Plot the mesh
	igl::opengl::glfw::Viewer viewer;
	igl::opengl::glfw::imgui::ImGuiMenu menu;
	viewer.plugins.push_back(&menu);

	//texture image
	std::string pngName = "../MyColorBar/MyColorBar2.png";

	bool down_on_mesh = false;
	//bool threadFlag = true;



	const auto update = [&]()->bool
	{
		int fid;
		Eigen::Vector3f bc;
		double x = viewer.current_mouse_x;//viewer.down_mouse_x
		double y = viewer.core().viewport(3) - viewer.current_mouse_y;//viewer.down_mouse_y
		if (igl::unproject_onto_mesh(Eigen::Vector2f(x, y), viewer.core().view, viewer.core().proj, viewer.core().viewport, V, F, fid, bc))
		{
			const Eigen::RowVector3d m3 = V.row(F(fid, 0)) * bc(0) + V.row(F(fid, 1)) * bc(1) + V.row(F(fid, 2)) * bc(2);
			int cid = 0;
			Eigen::Vector3d(
				(V.row(F(fid, 0)) - m3).squaredNorm(),
				(V.row(F(fid, 1)) - m3).squaredNorm(),
				(V.row(F(fid, 2)) - m3).squaredNorm()).minCoeff(&cid);
			sourceVi = (size_t)F.coeffRef(fid, cid);

			//compute distance
			GeoInHeat.computeDistance(sourceVi, geodistance);

			//set texture
			setTexture(pngName, geodistance, viewer);

			return true;
		}
		return false;
	};

	//MEUN----BEGIN
	// Draw additional windows
	double smooth = 1.0;
	menu.callback_draw_custom_window = [&]()
	{
		// Define next window position + size
		double w = 200.0, h = 160.0;
		ImGui::SetNextWindowPos(ImVec2(180.f * menu.menu_scaling(), 10), ImGuiCond_FirstUseEver);
		ImGui::SetNextWindowSize(ImVec2(w, h), ImGuiCond_FirstUseEver);
		ImGui::Begin(
			"Single Source:", nullptr,
			ImGuiWindowFlags_NoSavedSettings
		);

		ImGui::Text("Vertex %llu : (%.2f, %.2f, %.2f)", sourceVi, V.coeffRef(sourceVi, 0), V.coeffRef(sourceVi, 1), V.coeffRef(sourceVi, 2));
		ImGui::End();

		//
		ImGui::SetNextWindowPos(ImVec2(180.f * menu.menu_scaling(), h + 20), ImGuiCond_FirstUseEver);
		ImGui::SetNextWindowSize(ImVec2(w, h), ImGuiCond_FirstUseEver);
		ImGui::Begin(
			"Min Distance From Single Source", nullptr,
			ImGuiWindowFlags_NoSavedSettings
		);
		for (size_t i = 0; i < geodistance.rows(); ++i) {
			ImGui::Text("Vertex %llu : %.9f", i, geodistance(i));
		}
		ImGui::End();

		//
		ImGui::SetNextWindowPos(ImVec2(180.f * menu.menu_scaling(), h * 2 + 20), ImGuiCond_FirstUseEver);
		ImGui::SetNextWindowSize(ImVec2(w * 1.1, h * 0.6), ImGuiCond_FirstUseEver);
		ImGui::Begin(
			"Seting", nullptr,
			ImGuiWindowFlags_NoSavedSettings
		);
		ImGui::PushItemWidth(-80);
		if (ImGui::InputDouble("Scale time step", &smooth, 0.00000005f, 0, "%.8f")) {
			//compute distance
			GeoInHeat.computeDistance(sourceVi, geodistance, smooth);
			
			//set texture
			setTexture(pngName, geodistance, viewer);
		}
		ImGui::PopItemWidth();

		ImGui::Text("Current time step: %.16f ", GeoInHeat.getTimeStep());
		ImGui::End();
	};
	//MEUN----END

	viewer.callback_mouse_down =
		[&](igl::opengl::glfw::Viewer& viewer, int button, int modifier)->bool
	{
		if (update())
		{
			down_on_mesh = true;
			return true;
		}
		return false;
	};
	
	viewer.callback_mouse_up =
		[&down_on_mesh](igl::opengl::glfw::Viewer& viewer, int button, int modifier)->bool
	{
		down_on_mesh = false;
		return false;
	};

	viewer.data().set_mesh(V, F);
	viewer.data().show_lines = false;
	viewer.launch();

	
	////test before use ui
	////UI??? sourceVi??? view distance???
	//std::string fileName = "../data/result1.obj";
	//igl::readOBJ(fileName, V, TC, N, F, FTC, FN);
	//TriangleMesh triMesh(&V,&F);
	//GeodesicInheat GeoInHeat(triMesh);

	//size_t sourceVi = 10;
	//Eigen::VectorXd geodistance;
	//GeoInHeat.computeDistance(sourceVi, geodistance);
}

