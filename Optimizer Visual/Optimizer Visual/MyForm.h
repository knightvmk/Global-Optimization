#pragma once

#include <omp.h>
#include <msclr\marshal_cppstd.h>
#include <msclr\marshal.h>
#include <msclr\marshal_windows.h>
#include "TGlobal.h"
#include <math.h>

const int PPM_BY_R_INTERVALS = 1;
const int PPM_BY_DIVISIONS = 2;

using namespace OpenTK::Platform::Windows;
using namespace OpenTK::Graphics::OpenGL;

namespace Optimizer_Visual {

	using namespace System;
	using namespace System::ComponentModel;
	using namespace System::Collections;
	using namespace System::Windows::Forms;
	using namespace System::Data;
	using namespace System::Drawing;

	/// <summary>
	/// Summary for MyForm
	/// </summary>
	public ref class MyForm : public System::Windows::Forms::Form
	{
	public:
		MyForm(void)
		{
			InitializeComponent();
			//Init Thread Max count
			InitParams();
			//
			//TODO: Add the constructor code here
			//
			CHART->Series->Clear();
			CHART->Visible = false;
			glControl1->Visible = false;
			statusStrip1->Items[0]->Text = "Ожидание инициализации функции и границ поиска...";
		}
		int DIMENSION = 0;
		int MAX_THREADS = 1;
		bool ready_for_use = false;
		bool optimized = false;
		bool draw_graph = true;
		double move_x = 0, move_y = 0, move_z = 0;
		double rotate_grad = 0;
		std::vector<double> *x_test_points = new std::vector<double>;

		std::vector<double> *y_test_points = new std::vector<double>;

		void InitParams()
		{
			MAX_THREADS = omp_get_max_threads();
			CPU_COUNT->Text = Convert::ToString(MAX_THREADS);
		}

		bool CheckDimension()
		{
			DIMENSION = 0;
			std::string expression;
			msclr::interop::marshal_context context;
			expression = context.marshal_as<std::string>(MAIN_FUNC->Text);
			TCALC *matexp = new TCALC();
			char *expr = new char[expression.length() + 1];
			try		
			{
				strcpy(expr, expression.c_str());
				matexp->Compile(expr, 1, 1);
				matexp->Evaluate();
				delete matexp;
				delete[]expr;
			}
			catch (...)
			{
				MessageBox::Show("Плохой ввод функции. Проверьте корректность ввода. Используйте инструкцию.\n", "Ошибка ввода", MessageBoxButtons::OK, MessageBoxIcon::Error);
				delete matexp;
				delete[]expr;
				ready_for_use = false;
				return false;
			}
			int j = expression.find("X");
			int k = expression.find("Y");
			if (j != std::string::npos) DIMENSION++;
			if (k != std::string::npos) DIMENSION++;
			if (k != std::string::npos && j == std::string::npos)
			{
				{
					SET_EXPRESSION->BackColor = Color::Red;
					MessageBox::Show("Плохой ввод функции. Проверьте ввод. Используйте инструкцию.\nДля одномерной функции должен быть только X, а не Y", "Ошибка ввода", MessageBoxButtons::OK, MessageBoxIcon::Error);
					ready_for_use = false;
					return false;
				}
			}
			return true;
		}

		void Optimize(bool is_minimize, std::string expression, double *Left, double *Right, double r, double Epsilon, double N_max, 
			int N_DIMENSION, bool IS_PARALLEL, int PARALLEL_MODE, int THREADS, 
			std::vector<double> &out_x_test_points, std::vector<double> &out_y_test_points)
		{
			statusStrip1->Items[0]->Text = "Поиск глобального экстремума...";
			if (N_DIMENSION == 1)
			{ 
				//Plot1D(expression, Left[0], Right[0], 300);
				if (IS_PARALLEL)
				{
					if (PARALLEL_MODE == PPM_BY_R_INTERVALS) // ПП по характеристикам
					{
						Parallel_1D test(Left[0], Right[0], r);
						test.SetExpression(expression);
						test.Optimize(Epsilon, N_max, THREADS);
						TIME_ELAPSED->Text = Convert::ToString(test.Time()) + L" сек";
						STEPS_ELAPSED->Text = Convert::ToString(test.GetSteps());
						if (is_minimize)
						{
							Z_VAL_OPT->Text = Convert::ToString(test.GetSolutionZ());
							X_VAL_OPT->Text = Convert::ToString(test.GetSolutionX());
							out_x_test_points = std::vector<double>(test.Arguments_X.begin(), test.Arguments_X.end());
							out_y_test_points = std::vector<double>(test.Values_Z.begin(), test.Values_Z.end());
							Plot1D_extrema(test.GetSolutionX(), test.GetSolutionZ());
						}
						else
						{
							Z_VAL_OPT->Text = Convert::ToString(-1*test.GetSolutionZ());
							X_VAL_OPT->Text = Convert::ToString(test.GetSolutionX());
							out_x_test_points = std::vector<double>(test.Arguments_X.begin(), test.Arguments_X.end());
							for (auto iter = test.Values_Z.begin(); iter != test.Values_Z.end(); iter++) *iter *= -1;
							out_y_test_points = std::vector<double>(test.Values_Z.begin(), test.Values_Z.end());
							Plot1D_extrema(test.GetSolutionX(), -1*test.GetSolutionZ());
						}
					}
					else if (PARALLEL_MODE == PPM_BY_DIVISIONS)
					{
						Parallel_1D test(Left[0], Right[0], r);
						test.SetExpression(expression);
						test.Optimize(Left[0], Right[0], N_max, Epsilon, THREADS);
						TIME_ELAPSED->Text = Convert::ToString(test.Time()) + L" сек";
						STEPS_ELAPSED->Text = Convert::ToString(test.GetSteps());
						if (is_minimize)
						{
							Z_VAL_OPT->Text = Convert::ToString(test.GetSolutionZ());
							X_VAL_OPT->Text = Convert::ToString(test.GetSolutionX());
							out_x_test_points = std::vector<double>(test.Arguments_X.begin(), test.Arguments_X.end());
							out_y_test_points = std::vector<double>(test.Values_Z.begin(), test.Values_Z.end());
							Plot1D_extrema(test.GetSolutionX(), test.GetSolutionZ());
						}
						else
						{
							Z_VAL_OPT->Text = Convert::ToString(-1 * test.GetSolutionZ());
							X_VAL_OPT->Text = Convert::ToString(test.GetSolutionX());
							out_x_test_points = std::vector<double>(test.Arguments_X.begin(), test.Arguments_X.end());
							for (auto iter = test.Values_Z.begin(); iter != test.Values_Z.end(); iter++) *iter *= -1;
							out_y_test_points = std::vector<double>(test.Values_Z.begin(), test.Values_Z.end());
							Plot1D_extrema(test.GetSolutionX(), -1 * test.GetSolutionZ());
						}
					}
				}
				else
				{
					Sequental_1D test(Left[0], Right[0], r);
					test.SetExpression(expression);
					test.Optimize(Left[0], Right[0], N_max, Epsilon);
					TIME_ELAPSED->Text = Convert::ToString(test.Time()) + L" сек";
					STEPS_ELAPSED->Text = Convert::ToString(test.GetSteps());
					if (is_minimize)
					{
						Z_VAL_OPT->Text = Convert::ToString(test.GetSolutionZ());
						X_VAL_OPT->Text = Convert::ToString(test.GetSolutionX());
						out_x_test_points = std::vector<double>(test.Arguments_X.begin(), test.Arguments_X.end());
						out_y_test_points = std::vector<double>(test.Values_Z.begin(), test.Values_Z.end());
						Plot1D_extrema(test.GetSolutionX(), test.GetSolutionZ());
					}
					else
					{
						Z_VAL_OPT->Text = Convert::ToString(-1 * test.GetSolutionZ());
						X_VAL_OPT->Text = Convert::ToString(test.GetSolutionX());
						out_x_test_points = std::vector<double>(test.Arguments_X.begin(), test.Arguments_X.end());
						for (auto iter = test.Values_Z.begin(); iter != test.Values_Z.end(); iter++) *iter *= -1;
						out_y_test_points = std::vector<double>(test.Values_Z.begin(), test.Values_Z.end());
						Plot1D_extrema(test.GetSolutionX(), -1 * test.GetSolutionZ());
					}
				}
			}
			else if (N_DIMENSION == 2)
			{
				if (IS_PARALLEL)
				{
					Parallel_2D test(Left, Right, r);
					test.SetExpression(expression);
					test.Optimize(Epsilon, N_max, THREADS);
					if(is_minimize)
					Z_VAL_OPT->Text = Convert::ToString(test.GetSolutionZ());
					else Z_VAL_OPT->Text = Convert::ToString(-1*test.GetSolutionZ());
					X_VAL_OPT->Text = Convert::ToString(test.GetSolutionX());
					Y_VAL_OPT->Text = Convert::ToString(test.GetSolutionY());
					TIME_ELAPSED->Text = Convert::ToString(test.Time()) + L" сек";
					STEPS_ELAPSED->Text = Convert::ToString(test.GetSteps());
					Plot2D_extrema();
				}
				else
				{
					Sequental_2D test(Left, Right, r);
					test.SetExpression(expression);
					test.Optimize(Left, Right, N_max, Epsilon);
					if (is_minimize)
						Z_VAL_OPT->Text = Convert::ToString(test.GetSolutionZ());
					else Z_VAL_OPT->Text = Convert::ToString(-1 * test.GetSolutionZ());
					X_VAL_OPT->Text = Convert::ToString(test.GetSolutionX());
					Y_VAL_OPT->Text = Convert::ToString(test.GetSolutionY());
					TIME_ELAPSED->Text = Convert::ToString(test.Time()) + L" сек";
					STEPS_ELAPSED->Text = Convert::ToString(test.GetSteps());
					Plot2D_extrema();
				}
			}
		}

		void Plot1D(std::string expression, double left, double right, double density)
		{
			CHART->Series->Clear();
			CHART->ChartAreas[0]->AxisX->Minimum = left;
			CHART->ChartAreas[0]->AxisX->Maximum = right;
			DataVisualization::Charting::Series ^Lines = gcnew DataVisualization::Charting::Series();
			Lines->ChartArea = L"ChartArea1";
			Lines->ChartType = DataVisualization::Charting::SeriesChartType::Line;
			Lines->Color = Color::MintCream;
			//Lines->BorderWidth += 1;
			char *expr = new char[expression.length() + 1];
			strcpy(expr, expression.c_str());
			double step = abs((right - left) / density);
			double x = left;
			double y = 0;
			TCALC *matexp = new TCALC();
			while (x < right)
			{
				matexp->Compile(expr, x, 1);
				matexp->Evaluate();
				y = matexp->GetResult();
				Lines->Points->AddXY(x, y);
				x += step;
			}
			CHART->Series->Add(Lines);
			delete matexp;
			delete[]expr;
		}
		void Plot1D_extrema(double x, double y)
		{
			DataVisualization::Charting::Series ^Point = gcnew DataVisualization::Charting::Series();
			Point->ChartArea = L"ChartArea1";
			Point->ChartType = DataVisualization::Charting::SeriesChartType::Point;
			Point->Color = Color::Red;
			Point->BorderWidth += 3;
			Point->Points->AddXY(x, y);
			CHART->Series->Add(Point);
		}
		void Plot1D_test_points(std::vector<double> &x, std::vector<double> &y)
		{
			DataVisualization::Charting::Series ^Point = gcnew DataVisualization::Charting::Series();
			Point->ChartArea = L"ChartArea1";
			Point->ChartType = DataVisualization::Charting::SeriesChartType::Point;
			//Point->Color = Color::LightGreen;
			Point->Color = Color::Chartreuse;
			//Point->BorderWidth += 1;
			for (int i = 0; i < x.size(); i++)
			{
				Point->Points->AddXY(x[i], y[i]);
			}
			CHART->Series->Add(Point);
		}

		void Plot2D(std::string expression, double *Left, double *Right, int density)
		{
			double step_x = abs((Right[0] - Left[0]) / density);
			double step_y = abs((Right[1] - Left[1]) / density);
			double x = Left[0];
			double y = Left[1];
			double z = 0;

			TCALC *matexp = new TCALC();

			char *expr = new char[expression.length() + 1];
			strcpy(expr, expression.c_str());

			double init_x = x;
			double init_y = y;

			double **Matrix = new double*[density];
			for (int i = 0; i < density; i++)
				Matrix[i] = new double[density];
//#pragma ivdep
//#pragma vector always
			for (int i = 0; i < density; i++)
			{
				x = Left[0];
				for (int j = 0; j < density; j++)
				{
					matexp->Compile(expr, x, y);
					matexp->Evaluate();
					Matrix[i][j] = matexp->GetResult();
					x += step_x;
				}
				y += step_y;
			}

			GL::Clear(ClearBufferMask::ColorBufferBit | ClearBufferMask::DepthBufferBit);

			// рисуем сам график
			y = Left[1];
			for (int i = 0; i < density - 1; i++)
			{
				x = Left[0];
				for (int j = 0; j < density - 1; j++)
				{
					GL::PointSize(1);
					GL::Color3(Color::White);
					GL::Begin(BeginMode::LineStrip);

					GL::Vertex3(x, y, Matrix[i][j]);
					GL::Vertex3(x + step_x, y, Matrix[i][j + 1]);
					GL::Vertex3(x + step_x, y + step_y, Matrix[i + 1][j + 1]);
					GL::Vertex3(x, y + step_y, Matrix[i + 1][j]);
					GL::Vertex3(x, y, Matrix[i][j]);

					GL::End();
					x += step_x;
				}
				y += step_y;
			}
			
			// рисуем оси
			{
				GL::PointSize(1);

				GL::Color3(Color::CadetBlue);
				GL::Begin(BeginMode::Lines);

				GL::Vertex3(0,0,0); // OX
				GL::Vertex3(density * 2, 0, 0);

				GL::End();

				GL::Color3(Color::OrangeRed);
				GL::Begin(BeginMode::Lines);

				GL::Vertex3(0, 0, 0); // OY
				GL::Vertex3(0, density * 2, 0);

				GL::End();

				GL::Color3(Color::LightYellow);
				GL::Begin(BeginMode::Lines);

				GL::Vertex3(0, 0, 0); // OZ
				GL::Vertex3(0, 0, density * 2);

				GL::End();
			}

			// рисуем точки экстремума
			if (optimized)
			{
				double x_extrema = Convert::ToDouble(X_VAL_OPT->Text);
				double y_extrema = Convert::ToDouble(Y_VAL_OPT->Text);
				double z_extrema = Convert::ToDouble(Z_VAL_OPT->Text);

				// точкой
				GL::PointSize(10);
				GL::Color3(Color::Red);
				GL::Begin(BeginMode::Points);

				GL::Vertex3(x_extrema, y_extrema, z_extrema);

				GL::End();

				// а теперь линию прицела

				GL::PointSize(1);
				GL::Color3(Color::Red);
				GL::Begin(BeginMode::Lines);

				GL::Vertex3(x_extrema, y_extrema, z_extrema-10);
				GL::Vertex3(x_extrema, y_extrema, z_extrema+10);

				GL::End();

				GL::Color3(Color::Red);
				GL::Begin(BeginMode::Lines);

				GL::Vertex3(x_extrema, y_extrema-10, z_extrema);
				GL::Vertex3(x_extrema, y_extrema+10, z_extrema);

				GL::End();

				GL::Color3(Color::Red);
				GL::Begin(BeginMode::Lines);

				GL::Vertex3(x_extrema-10, y_extrema, z_extrema);
				GL::Vertex3(x_extrema+10, y_extrema, z_extrema);

				GL::End();
			}
	
			glControl1->SwapBuffers();
						

			delete[]expr;
			delete matexp;
			for (int i = 0; i < density; i++)
				delete Matrix[i];
			delete[]Matrix;
		}
		void Plot2D_extrema()
		{
			glControl1->Invalidate();
		}

		void Rotate(double gradus)
		{
			X_LEFT->Text->Replace('.', ',');
			Y_LEFT->Text->Replace('.', ',');
			X_RIGHT->Text->Replace('.', ',');
			Y_RIGHT->Text->Replace('.', ',');
			R_PARAM->Text->Replace('.', ',');

			double *Left = new double[2];
			double *Right = new double[2];

			Left[0] = Convert::ToDouble(X_LEFT->Text);
			Left[1] = Convert::ToDouble(Y_LEFT->Text);
			Right[0] = Convert::ToDouble(X_RIGHT->Text);
			Right[1] = Convert::ToDouble(Y_RIGHT->Text);

			double x = abs((Right[0] - Left[0]));
			double y = abs((Right[1] - Left[1]));

			GL::MatrixMode(MatrixMode::Modelview);
			GL::Rotate(gradus, 0, 0, abs(sqrt((pow(x,2)+pow(y,2))))/2);

			glControl1->Invalidate();
			//glControl1->Refresh();

			delete[]Left;
			delete[]Right;
		}
		void MoveCamera(double xstep, double ystep, double zstep)
		{
			X_LEFT->Text->Replace('.', ',');
			Y_LEFT->Text->Replace('.', ',');
			X_RIGHT->Text->Replace('.', ',');
			Y_RIGHT->Text->Replace('.', ',');
			R_PARAM->Text->Replace('.', ',');

			double *Left = new double[2];
			double *Right = new double[2];

			Left[0] = Convert::ToDouble(X_LEFT->Text);
			Left[1] = Convert::ToDouble(Y_LEFT->Text);
			Right[0] = Convert::ToDouble(X_RIGHT->Text);
			Right[1] = Convert::ToDouble(Y_RIGHT->Text);

			double x = (Right[0] - Left[0]) / 2;
			double y = (Right[1] - Left[1]) / 2;

			GL::ClearColor(Color::MidnightBlue);
			GL::Enable(EnableCap::DepthTest);

			//OpenTK::Matrix4 p = OpenTK::Matrix4::CreatePerspectiveFieldOfView((float)(50 * PI / 180), (float)(glControl1->Width / glControl1->Height), (float)1, (float)1000);
			//GL::MatrixMode(MatrixMode::Projection);
			//GL::LoadMatrix(p);

			OpenTK::Matrix4 modelview = OpenTK::Matrix4::LookAt(x*2.5 + y + xstep, y*2.5 + x + ystep, (x + y) * 3 + zstep,
				x, y, abs(sqrt((pow(x * 2, 2) + pow(y * 2, 2)))) / 2 + (x + y) / 2,
				0, 0, abs(sqrt((pow(x * 2, 2) + pow(y * 2, 2)))) / 2 + (x + y) / 2);
			GL::MatrixMode(MatrixMode::Modelview);
			GL::LoadMatrix(modelview);

			glControl1->Invalidate();
			//glControl1->Refresh();

			delete[]Left;
			delete[]Right;

		}
	protected:
		/// <summary>
		/// Clean up any resources being used.
		/// </summary>
		~MyForm()
		{
			if (components)
			{
				delete components;
			}
		}
	private: System::Windows::Forms::CheckBox^  SHOW_TEST_POINTS;
	private: System::Windows::Forms::StatusStrip^  statusStrip1;
	private: System::Windows::Forms::ToolStripStatusLabel^  toolStripStatusLabel1;
	private: OpenTK::GLControl^  glControl1;
	private: System::Windows::Forms::CheckBox^  CAN_DRAW;
	private: System::Windows::Forms::ToolStripMenuItem^  справкаToolStripMenuItem;
	private: System::Windows::Forms::ToolStripMenuItem^  оПрограммеToolStripMenuItem;
	private: System::Windows::Forms::ToolStripMenuItem^  вводФункцийToolStripMenuItem;
	private: System::Windows::Forms::Button^  LEFT;
	private: System::Windows::Forms::Button^  RIGHT;
	private: System::Windows::Forms::Button^  DOWN;
	private: System::Windows::Forms::Button^  UP;
	private: System::Windows::Forms::Button^  ROT_RIGHT;
	private: System::Windows::Forms::Button^  ROT_LEFT;
	private: System::Windows::Forms::MenuStrip^  menuStrip1;
	private: System::Windows::Forms::GroupBox^  groupBox7;
	private: System::Windows::Forms::DataVisualization::Charting::Chart^  CHART;
	private: System::Windows::Forms::ToolStripMenuItem^  файлToolStripMenuItem;
	private: System::Windows::Forms::ToolStripSeparator^  toolStripMenuItem1;
	private: System::Windows::Forms::ToolStripMenuItem^  выходToolStripMenuItem;
	private: System::Windows::Forms::GroupBox^  groupBox1;
	private: System::Windows::Forms::GroupBox^  groupBox3;
	private: System::Windows::Forms::GroupBox^  groupBox2;
	private: System::Windows::Forms::TextBox^  MAIN_FUNC;
	private: System::Windows::Forms::GroupBox^  groupBox4;
	private: System::Windows::Forms::GroupBox^  groupBox5;
	private: System::Windows::Forms::Label^  label1;
	private: System::Windows::Forms::TextBox^  CPU_COUNT;
	private: System::Windows::Forms::RadioButton^  ON_OPTIMIZE;
	private: System::Windows::Forms::RadioButton^  OFF_OPTIMIZE;
	private: System::Windows::Forms::RadioButton^  PP_BY_INTERVALS;
	private: System::Windows::Forms::RadioButton^  PP_BY_R_MAX;
	private: System::Windows::Forms::Button^  button1;
	private: System::Windows::Forms::Label^  label14;
	private: System::Windows::Forms::Label^  label15;
	private: System::Windows::Forms::Label^  label12;
	private: System::Windows::Forms::Label^  label13;
	private: System::Windows::Forms::Label^  label11;
	private: System::Windows::Forms::Label^  label10;
	private: System::Windows::Forms::TextBox^  Y_RIGHT;
	private: System::Windows::Forms::TextBox^  X_RIGHT;
	private: System::Windows::Forms::TextBox^  Y_LEFT;
	private: System::Windows::Forms::Label^  label8;
	private: System::Windows::Forms::TextBox^  X_LEFT;
	private: System::Windows::Forms::Label^  label9;
	private: System::Windows::Forms::Label^  label7;
	private: System::Windows::Forms::Label^  label6;
	private: System::Windows::Forms::Label^  label2;
	private: System::Windows::Forms::GroupBox^  groupBox6;
	private: System::Windows::Forms::TextBox^  Y_VAL_OPT;
	private: System::Windows::Forms::Label^  label5;
	private: System::Windows::Forms::TextBox^  X_VAL_OPT;
	private: System::Windows::Forms::Label^  label4;
	private: System::Windows::Forms::TextBox^  Z_VAL_OPT;
	private: System::Windows::Forms::Label^  label3;
	private: System::Windows::Forms::Button^  SET_EXPRESSION;
	private: System::Windows::Forms::Button^  GET_MIN;
	private: System::Windows::Forms::Button^  GET_MAX;
	private: System::Windows::Forms::TextBox^  TIME_ELAPSED;
	private: System::Windows::Forms::Label^  label16;
	private: System::Windows::Forms::TextBox^  STEPS_ELAPSED;
	private: System::Windows::Forms::Label^  labelnn2;
	private: System::Windows::Forms::TextBox^  EPS;
	private: System::Windows::Forms::Label^  labeln;
	private: System::Windows::Forms::TextBox^  STEPS;
	private: System::Windows::Forms::Label^  labeln1;
	private: System::Windows::Forms::TextBox^  R_PARAM;
	private: System::Windows::Forms::Label^  label17;
	private:
		/// <summary>
		/// Required designer variable.
		/// </summary>
		System::ComponentModel::Container ^components;

#pragma region Windows Form Designer generated code
		/// <summary>
		/// Required method for Designer support - do not modify
		/// the contents of this method with the code editor.
		/// </summary>
		void InitializeComponent(void)
		{
			System::ComponentModel::ComponentResourceManager^  resources = (gcnew System::ComponentModel::ComponentResourceManager(MyForm::typeid));
			System::Windows::Forms::DataVisualization::Charting::ChartArea^  chartArea1 = (gcnew System::Windows::Forms::DataVisualization::Charting::ChartArea());
			System::Windows::Forms::DataVisualization::Charting::Legend^  legend1 = (gcnew System::Windows::Forms::DataVisualization::Charting::Legend());
			System::Windows::Forms::DataVisualization::Charting::Series^  series1 = (gcnew System::Windows::Forms::DataVisualization::Charting::Series());
			this->menuStrip1 = (gcnew System::Windows::Forms::MenuStrip());
			this->файлToolStripMenuItem = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->toolStripMenuItem1 = (gcnew System::Windows::Forms::ToolStripSeparator());
			this->выходToolStripMenuItem = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->справкаToolStripMenuItem = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->вводФункцийToolStripMenuItem = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->оПрограммеToolStripMenuItem = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->groupBox1 = (gcnew System::Windows::Forms::GroupBox());
			this->EPS = (gcnew System::Windows::Forms::TextBox());
			this->labeln = (gcnew System::Windows::Forms::Label());
			this->STEPS = (gcnew System::Windows::Forms::TextBox());
			this->labeln1 = (gcnew System::Windows::Forms::Label());
			this->R_PARAM = (gcnew System::Windows::Forms::TextBox());
			this->label17 = (gcnew System::Windows::Forms::Label());
			this->groupBox3 = (gcnew System::Windows::Forms::GroupBox());
			this->button1 = (gcnew System::Windows::Forms::Button());
			this->label14 = (gcnew System::Windows::Forms::Label());
			this->label15 = (gcnew System::Windows::Forms::Label());
			this->label12 = (gcnew System::Windows::Forms::Label());
			this->label13 = (gcnew System::Windows::Forms::Label());
			this->label11 = (gcnew System::Windows::Forms::Label());
			this->label10 = (gcnew System::Windows::Forms::Label());
			this->Y_RIGHT = (gcnew System::Windows::Forms::TextBox());
			this->X_RIGHT = (gcnew System::Windows::Forms::TextBox());
			this->Y_LEFT = (gcnew System::Windows::Forms::TextBox());
			this->label8 = (gcnew System::Windows::Forms::Label());
			this->X_LEFT = (gcnew System::Windows::Forms::TextBox());
			this->label9 = (gcnew System::Windows::Forms::Label());
			this->label7 = (gcnew System::Windows::Forms::Label());
			this->label6 = (gcnew System::Windows::Forms::Label());
			this->groupBox2 = (gcnew System::Windows::Forms::GroupBox());
			this->SET_EXPRESSION = (gcnew System::Windows::Forms::Button());
			this->label2 = (gcnew System::Windows::Forms::Label());
			this->MAIN_FUNC = (gcnew System::Windows::Forms::TextBox());
			this->groupBox4 = (gcnew System::Windows::Forms::GroupBox());
			this->CAN_DRAW = (gcnew System::Windows::Forms::CheckBox());
			this->groupBox5 = (gcnew System::Windows::Forms::GroupBox());
			this->PP_BY_INTERVALS = (gcnew System::Windows::Forms::RadioButton());
			this->PP_BY_R_MAX = (gcnew System::Windows::Forms::RadioButton());
			this->label1 = (gcnew System::Windows::Forms::Label());
			this->CPU_COUNT = (gcnew System::Windows::Forms::TextBox());
			this->ON_OPTIMIZE = (gcnew System::Windows::Forms::RadioButton());
			this->OFF_OPTIMIZE = (gcnew System::Windows::Forms::RadioButton());
			this->groupBox6 = (gcnew System::Windows::Forms::GroupBox());
			this->SHOW_TEST_POINTS = (gcnew System::Windows::Forms::CheckBox());
			this->STEPS_ELAPSED = (gcnew System::Windows::Forms::TextBox());
			this->labelnn2 = (gcnew System::Windows::Forms::Label());
			this->TIME_ELAPSED = (gcnew System::Windows::Forms::TextBox());
			this->label16 = (gcnew System::Windows::Forms::Label());
			this->Y_VAL_OPT = (gcnew System::Windows::Forms::TextBox());
			this->label5 = (gcnew System::Windows::Forms::Label());
			this->X_VAL_OPT = (gcnew System::Windows::Forms::TextBox());
			this->label4 = (gcnew System::Windows::Forms::Label());
			this->Z_VAL_OPT = (gcnew System::Windows::Forms::TextBox());
			this->label3 = (gcnew System::Windows::Forms::Label());
			this->GET_MIN = (gcnew System::Windows::Forms::Button());
			this->GET_MAX = (gcnew System::Windows::Forms::Button());
			this->groupBox7 = (gcnew System::Windows::Forms::GroupBox());
			this->ROT_RIGHT = (gcnew System::Windows::Forms::Button());
			this->ROT_LEFT = (gcnew System::Windows::Forms::Button());
			this->LEFT = (gcnew System::Windows::Forms::Button());
			this->RIGHT = (gcnew System::Windows::Forms::Button());
			this->DOWN = (gcnew System::Windows::Forms::Button());
			this->UP = (gcnew System::Windows::Forms::Button());
			this->glControl1 = (gcnew OpenTK::GLControl());
			this->CHART = (gcnew System::Windows::Forms::DataVisualization::Charting::Chart());
			this->statusStrip1 = (gcnew System::Windows::Forms::StatusStrip());
			this->toolStripStatusLabel1 = (gcnew System::Windows::Forms::ToolStripStatusLabel());
			this->menuStrip1->SuspendLayout();
			this->groupBox1->SuspendLayout();
			this->groupBox3->SuspendLayout();
			this->groupBox2->SuspendLayout();
			this->groupBox4->SuspendLayout();
			this->groupBox5->SuspendLayout();
			this->groupBox6->SuspendLayout();
			this->groupBox7->SuspendLayout();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->CHART))->BeginInit();
			this->statusStrip1->SuspendLayout();
			this->SuspendLayout();
			// 
			// menuStrip1
			// 
			this->menuStrip1->ImageScalingSize = System::Drawing::Size(19, 19);
			this->menuStrip1->Items->AddRange(gcnew cli::array< System::Windows::Forms::ToolStripItem^  >(3) {
				this->файлToolStripMenuItem,
					this->справкаToolStripMenuItem, this->оПрограммеToolStripMenuItem
			});
			this->menuStrip1->Location = System::Drawing::Point(0, 0);
			this->menuStrip1->Name = L"menuStrip1";
			this->menuStrip1->Size = System::Drawing::Size(1194, 24);
			this->menuStrip1->TabIndex = 0;
			this->menuStrip1->Text = L"menuStrip1";
			// 
			// файлToolStripMenuItem
			// 
			this->файлToolStripMenuItem->DropDownItems->AddRange(gcnew cli::array< System::Windows::Forms::ToolStripItem^  >(2) {
				this->toolStripMenuItem1,
					this->выходToolStripMenuItem
			});
			this->файлToolStripMenuItem->Name = L"файлToolStripMenuItem";
			this->файлToolStripMenuItem->Size = System::Drawing::Size(48, 20);
			this->файлToolStripMenuItem->Text = L"Файл";
			// 
			// toolStripMenuItem1
			// 
			this->toolStripMenuItem1->Name = L"toolStripMenuItem1";
			this->toolStripMenuItem1->Size = System::Drawing::Size(105, 6);
			// 
			// выходToolStripMenuItem
			// 
			this->выходToolStripMenuItem->Name = L"выходToolStripMenuItem";
			this->выходToolStripMenuItem->Size = System::Drawing::Size(108, 22);
			this->выходToolStripMenuItem->Text = L"Выход";
			this->выходToolStripMenuItem->Click += gcnew System::EventHandler(this, &MyForm::выходToolStripMenuItem_Click);
			// 
			// справкаToolStripMenuItem
			// 
			this->справкаToolStripMenuItem->DropDownItems->AddRange(gcnew cli::array< System::Windows::Forms::ToolStripItem^  >(1) { this->вводФункцийToolStripMenuItem });
			this->справкаToolStripMenuItem->Name = L"справкаToolStripMenuItem";
			this->справкаToolStripMenuItem->Size = System::Drawing::Size(65, 20);
			this->справкаToolStripMenuItem->Text = L"Справка";
			// 
			// вводФункцийToolStripMenuItem
			// 
			this->вводФункцийToolStripMenuItem->Name = L"вводФункцийToolStripMenuItem";
			this->вводФункцийToolStripMenuItem->Size = System::Drawing::Size(152, 22);
			this->вводФункцийToolStripMenuItem->Text = L"Ввод функций";
			this->вводФункцийToolStripMenuItem->Click += gcnew System::EventHandler(this, &MyForm::вводФункцийToolStripMenuItem_Click);
			// 
			// оПрограммеToolStripMenuItem
			// 
			this->оПрограммеToolStripMenuItem->Name = L"оПрограммеToolStripMenuItem";
			this->оПрограммеToolStripMenuItem->Size = System::Drawing::Size(94, 20);
			this->оПрограммеToolStripMenuItem->Text = L"О программе";
			this->оПрограммеToolStripMenuItem->Click += gcnew System::EventHandler(this, &MyForm::оПрограммеToolStripMenuItem_Click);
			// 
			// groupBox1
			// 
			this->groupBox1->BackColor = System::Drawing::Color::Transparent;
			this->groupBox1->BackgroundImageLayout = System::Windows::Forms::ImageLayout::Zoom;
			this->groupBox1->Controls->Add(this->EPS);
			this->groupBox1->Controls->Add(this->labeln);
			this->groupBox1->Controls->Add(this->STEPS);
			this->groupBox1->Controls->Add(this->labeln1);
			this->groupBox1->Controls->Add(this->R_PARAM);
			this->groupBox1->Controls->Add(this->label17);
			this->groupBox1->Controls->Add(this->groupBox3);
			this->groupBox1->Controls->Add(this->groupBox2);
			this->groupBox1->ForeColor = System::Drawing::SystemColors::ControlLightLight;
			this->groupBox1->Location = System::Drawing::Point(4, 27);
			this->groupBox1->Name = L"groupBox1";
			this->groupBox1->Size = System::Drawing::Size(303, 224);
			this->groupBox1->TabIndex = 1;
			this->groupBox1->TabStop = false;
			this->groupBox1->Text = L"[ Входные данные глобального поиска ]";
			// 
			// EPS
			// 
			this->EPS->Location = System::Drawing::Point(187, 197);
			this->EPS->Name = L"EPS";
			this->EPS->Size = System::Drawing::Size(104, 20);
			this->EPS->TabIndex = 16;
			this->EPS->Text = L"0,00001";
			// 
			// labeln
			// 
			this->labeln->AutoSize = true;
			this->labeln->Location = System::Drawing::Point(123, 200);
			this->labeln->Name = L"labeln";
			this->labeln->Size = System::Drawing::Size(57, 13);
			this->labeln->TabIndex = 15;
			this->labeln->Text = L"Точность:";
			// 
			// STEPS
			// 
			this->STEPS->Location = System::Drawing::Point(187, 171);
			this->STEPS->Name = L"STEPS";
			this->STEPS->Size = System::Drawing::Size(104, 20);
			this->STEPS->TabIndex = 14;
			this->STEPS->Text = L"300";
			// 
			// labeln1
			// 
			this->labeln1->AutoSize = true;
			this->labeln1->Location = System::Drawing::Point(139, 174);
			this->labeln1->Name = L"labeln1";
			this->labeln1->Size = System::Drawing::Size(42, 13);
			this->labeln1->TabIndex = 13;
			this->labeln1->Text = L"Шагов:";
			// 
			// R_PARAM
			// 
			this->R_PARAM->Location = System::Drawing::Point(85, 171);
			this->R_PARAM->Name = L"R_PARAM";
			this->R_PARAM->Size = System::Drawing::Size(42, 20);
			this->R_PARAM->TabIndex = 10;
			this->R_PARAM->Text = L"2";
			// 
			// label17
			// 
			this->label17->AutoSize = true;
			this->label17->Location = System::Drawing::Point(12, 174);
			this->label17->Name = L"label17";
			this->label17->Size = System::Drawing::Size(67, 13);
			this->label17->TabIndex = 2;
			this->label17->Text = L"Параметр r:";
			// 
			// groupBox3
			// 
			this->groupBox3->BackColor = System::Drawing::Color::Transparent;
			this->groupBox3->Controls->Add(this->button1);
			this->groupBox3->Controls->Add(this->label14);
			this->groupBox3->Controls->Add(this->label15);
			this->groupBox3->Controls->Add(this->label12);
			this->groupBox3->Controls->Add(this->label13);
			this->groupBox3->Controls->Add(this->label11);
			this->groupBox3->Controls->Add(this->label10);
			this->groupBox3->Controls->Add(this->Y_RIGHT);
			this->groupBox3->Controls->Add(this->X_RIGHT);
			this->groupBox3->Controls->Add(this->Y_LEFT);
			this->groupBox3->Controls->Add(this->label8);
			this->groupBox3->Controls->Add(this->X_LEFT);
			this->groupBox3->Controls->Add(this->label9);
			this->groupBox3->Controls->Add(this->label7);
			this->groupBox3->Controls->Add(this->label6);
			this->groupBox3->ForeColor = System::Drawing::SystemColors::ControlLightLight;
			this->groupBox3->Location = System::Drawing::Point(9, 75);
			this->groupBox3->Name = L"groupBox3";
			this->groupBox3->Size = System::Drawing::Size(288, 92);
			this->groupBox3->TabIndex = 1;
			this->groupBox3->TabStop = false;
			this->groupBox3->Text = L"Линейные ограничения";
			// 
			// button1
			// 
			this->button1->Image = (cli::safe_cast<System::Drawing::Image^>(resources->GetObject(L"button1.Image")));
			this->button1->Location = System::Drawing::Point(244, 36);
			this->button1->Name = L"button1";
			this->button1->Size = System::Drawing::Size(38, 46);
			this->button1->TabIndex = 18;
			this->button1->UseVisualStyleBackColor = true;
			this->button1->Click += gcnew System::EventHandler(this, &MyForm::button1_Click);
			// 
			// label14
			// 
			this->label14->AutoSize = true;
			this->label14->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 9.75F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label14->Location = System::Drawing::Point(123, 63);
			this->label14->Name = L"label14";
			this->label14->Size = System::Drawing::Size(11, 16);
			this->label14->TabIndex = 17;
			this->label14->Text = L":";
			// 
			// label15
			// 
			this->label15->AutoSize = true;
			this->label15->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 9.75F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label15->Location = System::Drawing::Point(123, 37);
			this->label15->Name = L"label15";
			this->label15->Size = System::Drawing::Size(11, 16);
			this->label15->TabIndex = 16;
			this->label15->Text = L":";
			// 
			// label12
			// 
			this->label12->AutoSize = true;
			this->label12->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 9.75F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label12->Location = System::Drawing::Point(28, 63);
			this->label12->Name = L"label12";
			this->label12->Size = System::Drawing::Size(12, 16);
			this->label12->TabIndex = 15;
			this->label12->Text = L"[";
			// 
			// label13
			// 
			this->label13->AutoSize = true;
			this->label13->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 9.75F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label13->Location = System::Drawing::Point(28, 37);
			this->label13->Name = L"label13";
			this->label13->Size = System::Drawing::Size(12, 16);
			this->label13->TabIndex = 14;
			this->label13->Text = L"[";
			// 
			// label11
			// 
			this->label11->AutoSize = true;
			this->label11->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 9.75F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label11->Location = System::Drawing::Point(224, 63);
			this->label11->Name = L"label11";
			this->label11->Size = System::Drawing::Size(12, 16);
			this->label11->TabIndex = 13;
			this->label11->Text = L"]";
			// 
			// label10
			// 
			this->label10->AutoSize = true;
			this->label10->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 9.75F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label10->Location = System::Drawing::Point(224, 37);
			this->label10->Name = L"label10";
			this->label10->Size = System::Drawing::Size(12, 16);
			this->label10->TabIndex = 12;
			this->label10->Text = L"]";
			// 
			// Y_RIGHT
			// 
			this->Y_RIGHT->Location = System::Drawing::Point(137, 62);
			this->Y_RIGHT->Name = L"Y_RIGHT";
			this->Y_RIGHT->Size = System::Drawing::Size(83, 20);
			this->Y_RIGHT->TabIndex = 11;
			this->Y_RIGHT->Text = L"0";
			// 
			// X_RIGHT
			// 
			this->X_RIGHT->Location = System::Drawing::Point(137, 36);
			this->X_RIGHT->Name = L"X_RIGHT";
			this->X_RIGHT->Size = System::Drawing::Size(83, 20);
			this->X_RIGHT->TabIndex = 10;
			this->X_RIGHT->Text = L"2";
			// 
			// Y_LEFT
			// 
			this->Y_LEFT->Location = System::Drawing::Point(40, 62);
			this->Y_LEFT->Name = L"Y_LEFT";
			this->Y_LEFT->Size = System::Drawing::Size(78, 20);
			this->Y_LEFT->TabIndex = 9;
			this->Y_LEFT->Text = L"-5";
			// 
			// label8
			// 
			this->label8->AutoSize = true;
			this->label8->Location = System::Drawing::Point(6, 65);
			this->label8->Name = L"label8";
			this->label8->Size = System::Drawing::Size(20, 13);
			this->label8->TabIndex = 8;
			this->label8->Text = L"Y :";
			// 
			// X_LEFT
			// 
			this->X_LEFT->Location = System::Drawing::Point(40, 36);
			this->X_LEFT->Name = L"X_LEFT";
			this->X_LEFT->Size = System::Drawing::Size(78, 20);
			this->X_LEFT->TabIndex = 7;
			this->X_LEFT->Text = L"-4";
			// 
			// label9
			// 
			this->label9->AutoSize = true;
			this->label9->Location = System::Drawing::Point(6, 39);
			this->label9->Name = L"label9";
			this->label9->Size = System::Drawing::Size(20, 13);
			this->label9->TabIndex = 6;
			this->label9->Text = L"X :";
			// 
			// label7
			// 
			this->label7->AutoSize = true;
			this->label7->Location = System::Drawing::Point(124, 20);
			this->label7->Name = L"label7";
			this->label7->Size = System::Drawing::Size(47, 13);
			this->label7->TabIndex = 1;
			this->label7->Text = L"Справа:";
			// 
			// label6
			// 
			this->label6->AutoSize = true;
			this->label6->Location = System::Drawing::Point(37, 20);
			this->label6->Name = L"label6";
			this->label6->Size = System::Drawing::Size(41, 13);
			this->label6->TabIndex = 0;
			this->label6->Text = L"Слева:";
			// 
			// groupBox2
			// 
			this->groupBox2->Controls->Add(this->SET_EXPRESSION);
			this->groupBox2->Controls->Add(this->label2);
			this->groupBox2->Controls->Add(this->MAIN_FUNC);
			this->groupBox2->ForeColor = System::Drawing::SystemColors::ControlLightLight;
			this->groupBox2->Location = System::Drawing::Point(9, 20);
			this->groupBox2->Name = L"groupBox2";
			this->groupBox2->Size = System::Drawing::Size(288, 49);
			this->groupBox2->TabIndex = 0;
			this->groupBox2->TabStop = false;
			this->groupBox2->Text = L"Функция";
			// 
			// SET_EXPRESSION
			// 
			this->SET_EXPRESSION->BackColor = System::Drawing::Color::DarkGray;
			this->SET_EXPRESSION->BackgroundImageLayout = System::Windows::Forms::ImageLayout::Center;
			this->SET_EXPRESSION->FlatStyle = System::Windows::Forms::FlatStyle::Popup;
			this->SET_EXPRESSION->Location = System::Drawing::Point(247, 20);
			this->SET_EXPRESSION->Name = L"SET_EXPRESSION";
			this->SET_EXPRESSION->Size = System::Drawing::Size(35, 21);
			this->SET_EXPRESSION->TabIndex = 2;
			this->SET_EXPRESSION->Text = L"OK";
			this->SET_EXPRESSION->UseVisualStyleBackColor = false;
			this->SET_EXPRESSION->Click += gcnew System::EventHandler(this, &MyForm::SET_EXPRESSION_Click);
			// 
			// label2
			// 
			this->label2->AutoSize = true;
			this->label2->Location = System::Drawing::Point(6, 23);
			this->label2->Name = L"label2";
			this->label2->Size = System::Drawing::Size(23, 13);
			this->label2->TabIndex = 1;
			this->label2->Text = L"Z =";
			// 
			// MAIN_FUNC
			// 
			this->MAIN_FUNC->Location = System::Drawing::Point(31, 20);
			this->MAIN_FUNC->Name = L"MAIN_FUNC";
			this->MAIN_FUNC->Size = System::Drawing::Size(214, 20);
			this->MAIN_FUNC->TabIndex = 0;
			this->MAIN_FUNC->Text = L"2 * cos(Y) + 3 * sin(X)";
			this->MAIN_FUNC->TextChanged += gcnew System::EventHandler(this, &MyForm::MAIN_FUNC_TextChanged);
			// 
			// groupBox4
			// 
			this->groupBox4->Anchor = static_cast<System::Windows::Forms::AnchorStyles>(((System::Windows::Forms::AnchorStyles::Top | System::Windows::Forms::AnchorStyles::Bottom)
				| System::Windows::Forms::AnchorStyles::Left));
			this->groupBox4->BackColor = System::Drawing::Color::Transparent;
			this->groupBox4->Controls->Add(this->CAN_DRAW);
			this->groupBox4->Controls->Add(this->groupBox5);
			this->groupBox4->Controls->Add(this->ON_OPTIMIZE);
			this->groupBox4->Controls->Add(this->OFF_OPTIMIZE);
			this->groupBox4->ForeColor = System::Drawing::SystemColors::ControlLightLight;
			this->groupBox4->Location = System::Drawing::Point(4, 250);
			this->groupBox4->Name = L"groupBox4";
			this->groupBox4->Size = System::Drawing::Size(303, 166);
			this->groupBox4->TabIndex = 2;
			this->groupBox4->TabStop = false;
			this->groupBox4->Text = L"[ Оптимизация вычислений ]";
			// 
			// CAN_DRAW
			// 
			this->CAN_DRAW->AutoSize = true;
			this->CAN_DRAW->Checked = true;
			this->CAN_DRAW->CheckState = System::Windows::Forms::CheckState::Checked;
			this->CAN_DRAW->Location = System::Drawing::Point(201, 20);
			this->CAN_DRAW->Name = L"CAN_DRAW";
			this->CAN_DRAW->Size = System::Drawing::Size(70, 30);
			this->CAN_DRAW->TabIndex = 3;
			this->CAN_DRAW->Text = L"Строить \r\nграфик";
			this->CAN_DRAW->UseVisualStyleBackColor = true;
			this->CAN_DRAW->CheckedChanged += gcnew System::EventHandler(this, &MyForm::CAN_DRAW_CheckedChanged);
			// 
			// groupBox5
			// 
			this->groupBox5->Controls->Add(this->PP_BY_INTERVALS);
			this->groupBox5->Controls->Add(this->PP_BY_R_MAX);
			this->groupBox5->Controls->Add(this->label1);
			this->groupBox5->Controls->Add(this->CPU_COUNT);
			this->groupBox5->Enabled = false;
			this->groupBox5->ForeColor = System::Drawing::SystemColors::ControlLightLight;
			this->groupBox5->Location = System::Drawing::Point(9, 65);
			this->groupBox5->Name = L"groupBox5";
			this->groupBox5->Size = System::Drawing::Size(266, 94);
			this->groupBox5->TabIndex = 2;
			this->groupBox5->TabStop = false;
			this->groupBox5->Text = L"Тонкие настройки ускорения";
			// 
			// PP_BY_INTERVALS
			// 
			this->PP_BY_INTERVALS->AutoSize = true;
			this->PP_BY_INTERVALS->Location = System::Drawing::Point(10, 72);
			this->PP_BY_INTERVALS->Name = L"PP_BY_INTERVALS";
			this->PP_BY_INTERVALS->Size = System::Drawing::Size(234, 17);
			this->PP_BY_INTERVALS->TabIndex = 3;
			this->PP_BY_INTERVALS->Text = L"Распараллеливание по отрезкам поиска";
			this->PP_BY_INTERVALS->UseVisualStyleBackColor = true;
			// 
			// PP_BY_R_MAX
			// 
			this->PP_BY_R_MAX->AutoSize = true;
			this->PP_BY_R_MAX->Checked = true;
			this->PP_BY_R_MAX->Location = System::Drawing::Point(10, 49);
			this->PP_BY_R_MAX->Name = L"PP_BY_R_MAX";
			this->PP_BY_R_MAX->Size = System::Drawing::Size(235, 17);
			this->PP_BY_R_MAX->TabIndex = 2;
			this->PP_BY_R_MAX->TabStop = true;
			this->PP_BY_R_MAX->Text = L"Распараллеливание по характеристикам";
			this->PP_BY_R_MAX->UseVisualStyleBackColor = true;
			// 
			// label1
			// 
			this->label1->AutoSize = true;
			this->label1->Location = System::Drawing::Point(7, 26);
			this->label1->Name = L"label1";
			this->label1->Size = System::Drawing::Size(132, 13);
			this->label1->TabIndex = 1;
			this->label1->Text = L"Количество потоков ЦП:";
			// 
			// CPU_COUNT
			// 
			this->CPU_COUNT->Location = System::Drawing::Point(178, 23);
			this->CPU_COUNT->Name = L"CPU_COUNT";
			this->CPU_COUNT->Size = System::Drawing::Size(77, 20);
			this->CPU_COUNT->TabIndex = 0;
			this->CPU_COUNT->TextChanged += gcnew System::EventHandler(this, &MyForm::textBox2_TextChanged);
			// 
			// ON_OPTIMIZE
			// 
			this->ON_OPTIMIZE->AutoSize = true;
			this->ON_OPTIMIZE->Location = System::Drawing::Point(9, 42);
			this->ON_OPTIMIZE->Name = L"ON_OPTIMIZE";
			this->ON_OPTIMIZE->Size = System::Drawing::Size(75, 17);
			this->ON_OPTIMIZE->TabIndex = 1;
			this->ON_OPTIMIZE->Text = L"Включено";
			this->ON_OPTIMIZE->UseVisualStyleBackColor = true;
			this->ON_OPTIMIZE->CheckedChanged += gcnew System::EventHandler(this, &MyForm::ON_OPTIMIZE_CheckedChanged);
			// 
			// OFF_OPTIMIZE
			// 
			this->OFF_OPTIMIZE->AutoSize = true;
			this->OFF_OPTIMIZE->Checked = true;
			this->OFF_OPTIMIZE->Location = System::Drawing::Point(9, 19);
			this->OFF_OPTIMIZE->Name = L"OFF_OPTIMIZE";
			this->OFF_OPTIMIZE->Size = System::Drawing::Size(83, 17);
			this->OFF_OPTIMIZE->TabIndex = 0;
			this->OFF_OPTIMIZE->TabStop = true;
			this->OFF_OPTIMIZE->Text = L"Выключено";
			this->OFF_OPTIMIZE->UseVisualStyleBackColor = true;
			this->OFF_OPTIMIZE->CheckedChanged += gcnew System::EventHandler(this, &MyForm::radioButton1_CheckedChanged);
			// 
			// groupBox6
			// 
			this->groupBox6->Anchor = static_cast<System::Windows::Forms::AnchorStyles>((System::Windows::Forms::AnchorStyles::Bottom | System::Windows::Forms::AnchorStyles::Left));
			this->groupBox6->BackColor = System::Drawing::Color::Transparent;
			this->groupBox6->Controls->Add(this->SHOW_TEST_POINTS);
			this->groupBox6->Controls->Add(this->STEPS_ELAPSED);
			this->groupBox6->Controls->Add(this->labelnn2);
			this->groupBox6->Controls->Add(this->TIME_ELAPSED);
			this->groupBox6->Controls->Add(this->label16);
			this->groupBox6->Controls->Add(this->Y_VAL_OPT);
			this->groupBox6->Controls->Add(this->label5);
			this->groupBox6->Controls->Add(this->X_VAL_OPT);
			this->groupBox6->Controls->Add(this->label4);
			this->groupBox6->Controls->Add(this->Z_VAL_OPT);
			this->groupBox6->Controls->Add(this->label3);
			this->groupBox6->Enabled = false;
			this->groupBox6->ForeColor = System::Drawing::SystemColors::ControlLightLight;
			this->groupBox6->Location = System::Drawing::Point(4, 449);
			this->groupBox6->Name = L"groupBox6";
			this->groupBox6->Size = System::Drawing::Size(303, 141);
			this->groupBox6->TabIndex = 3;
			this->groupBox6->TabStop = false;
			this->groupBox6->Text = L"[ Результаты ]";
			// 
			// SHOW_TEST_POINTS
			// 
			this->SHOW_TEST_POINTS->AutoSize = true;
			this->SHOW_TEST_POINTS->Enabled = false;
			this->SHOW_TEST_POINTS->Location = System::Drawing::Point(201, 20);
			this->SHOW_TEST_POINTS->Name = L"SHOW_TEST_POINTS";
			this->SHOW_TEST_POINTS->Size = System::Drawing::Size(81, 30);
			this->SHOW_TEST_POINTS->TabIndex = 10;
			this->SHOW_TEST_POINTS->Text = L"Точки \r\nиспытаний";
			this->SHOW_TEST_POINTS->UseVisualStyleBackColor = true;
			this->SHOW_TEST_POINTS->CheckedChanged += gcnew System::EventHandler(this, &MyForm::SHOW_TEST_POINTS_CheckedChanged);
			// 
			// STEPS_ELAPSED
			// 
			this->STEPS_ELAPSED->Location = System::Drawing::Point(58, 114);
			this->STEPS_ELAPSED->Name = L"STEPS_ELAPSED";
			this->STEPS_ELAPSED->ReadOnly = true;
			this->STEPS_ELAPSED->Size = System::Drawing::Size(129, 20);
			this->STEPS_ELAPSED->TabIndex = 9;
			// 
			// labelnn2
			// 
			this->labelnn2->AutoSize = true;
			this->labelnn2->Location = System::Drawing::Point(9, 117);
			this->labelnn2->Name = L"labelnn2";
			this->labelnn2->Size = System::Drawing::Size(42, 13);
			this->labelnn2->TabIndex = 8;
			this->labelnn2->Text = L"Шагов:";
			// 
			// TIME_ELAPSED
			// 
			this->TIME_ELAPSED->Location = System::Drawing::Point(58, 92);
			this->TIME_ELAPSED->Name = L"TIME_ELAPSED";
			this->TIME_ELAPSED->ReadOnly = true;
			this->TIME_ELAPSED->Size = System::Drawing::Size(129, 20);
			this->TIME_ELAPSED->TabIndex = 7;
			// 
			// label16
			// 
			this->label16->AutoSize = true;
			this->label16->Location = System::Drawing::Point(9, 95);
			this->label16->Name = L"label16";
			this->label16->Size = System::Drawing::Size(43, 13);
			this->label16->TabIndex = 6;
			this->label16->Text = L"Время:";
			// 
			// Y_VAL_OPT
			// 
			this->Y_VAL_OPT->Location = System::Drawing::Point(40, 69);
			this->Y_VAL_OPT->Name = L"Y_VAL_OPT";
			this->Y_VAL_OPT->ReadOnly = true;
			this->Y_VAL_OPT->Size = System::Drawing::Size(147, 20);
			this->Y_VAL_OPT->TabIndex = 5;
			// 
			// label5
			// 
			this->label5->AutoSize = true;
			this->label5->Location = System::Drawing::Point(9, 72);
			this->label5->Name = L"label5";
			this->label5->Size = System::Drawing::Size(20, 13);
			this->label5->TabIndex = 4;
			this->label5->Text = L"Y :";
			// 
			// X_VAL_OPT
			// 
			this->X_VAL_OPT->Location = System::Drawing::Point(40, 43);
			this->X_VAL_OPT->Name = L"X_VAL_OPT";
			this->X_VAL_OPT->ReadOnly = true;
			this->X_VAL_OPT->Size = System::Drawing::Size(147, 20);
			this->X_VAL_OPT->TabIndex = 3;
			// 
			// label4
			// 
			this->label4->AutoSize = true;
			this->label4->Location = System::Drawing::Point(9, 46);
			this->label4->Name = L"label4";
			this->label4->Size = System::Drawing::Size(20, 13);
			this->label4->TabIndex = 2;
			this->label4->Text = L"X :";
			// 
			// Z_VAL_OPT
			// 
			this->Z_VAL_OPT->Location = System::Drawing::Point(40, 17);
			this->Z_VAL_OPT->Name = L"Z_VAL_OPT";
			this->Z_VAL_OPT->ReadOnly = true;
			this->Z_VAL_OPT->Size = System::Drawing::Size(147, 20);
			this->Z_VAL_OPT->TabIndex = 1;
			// 
			// label3
			// 
			this->label3->AutoSize = true;
			this->label3->Location = System::Drawing::Point(9, 20);
			this->label3->Name = L"label3";
			this->label3->Size = System::Drawing::Size(20, 13);
			this->label3->TabIndex = 0;
			this->label3->Text = L"Z :";
			// 
			// GET_MIN
			// 
			this->GET_MIN->Anchor = static_cast<System::Windows::Forms::AnchorStyles>((System::Windows::Forms::AnchorStyles::Bottom | System::Windows::Forms::AnchorStyles::Left));
			this->GET_MIN->Location = System::Drawing::Point(176, 419);
			this->GET_MIN->Name = L"GET_MIN";
			this->GET_MIN->Size = System::Drawing::Size(131, 23);
			this->GET_MIN->TabIndex = 4;
			this->GET_MIN->Text = L"Найти минимум";
			this->GET_MIN->UseVisualStyleBackColor = true;
			this->GET_MIN->Click += gcnew System::EventHandler(this, &MyForm::button2_Click);
			// 
			// GET_MAX
			// 
			this->GET_MAX->Anchor = static_cast<System::Windows::Forms::AnchorStyles>((System::Windows::Forms::AnchorStyles::Bottom | System::Windows::Forms::AnchorStyles::Left));
			this->GET_MAX->Location = System::Drawing::Point(4, 419);
			this->GET_MAX->Name = L"GET_MAX";
			this->GET_MAX->Size = System::Drawing::Size(131, 23);
			this->GET_MAX->TabIndex = 5;
			this->GET_MAX->Text = L"Найти максимум";
			this->GET_MAX->UseVisualStyleBackColor = true;
			this->GET_MAX->Click += gcnew System::EventHandler(this, &MyForm::GET_MAX_Click);
			// 
			// groupBox7
			// 
			this->groupBox7->Anchor = static_cast<System::Windows::Forms::AnchorStyles>((((System::Windows::Forms::AnchorStyles::Top | System::Windows::Forms::AnchorStyles::Bottom)
				| System::Windows::Forms::AnchorStyles::Left)
				| System::Windows::Forms::AnchorStyles::Right));
			this->groupBox7->BackColor = System::Drawing::Color::Transparent;
			this->groupBox7->Controls->Add(this->ROT_RIGHT);
			this->groupBox7->Controls->Add(this->ROT_LEFT);
			this->groupBox7->Controls->Add(this->LEFT);
			this->groupBox7->Controls->Add(this->RIGHT);
			this->groupBox7->Controls->Add(this->DOWN);
			this->groupBox7->Controls->Add(this->UP);
			this->groupBox7->Controls->Add(this->glControl1);
			this->groupBox7->Controls->Add(this->CHART);
			this->groupBox7->ForeColor = System::Drawing::SystemColors::ControlLightLight;
			this->groupBox7->Location = System::Drawing::Point(314, 28);
			this->groupBox7->Name = L"groupBox7";
			this->groupBox7->Size = System::Drawing::Size(876, 562);
			this->groupBox7->TabIndex = 6;
			this->groupBox7->TabStop = false;
			this->groupBox7->Text = L"[ Визуализация ]";
			// 
			// ROT_RIGHT
			// 
			this->ROT_RIGHT->Anchor = static_cast<System::Windows::Forms::AnchorStyles>((System::Windows::Forms::AnchorStyles::Bottom | System::Windows::Forms::AnchorStyles::Left));
			this->ROT_RIGHT->BackColor = System::Drawing::Color::SteelBlue;
			this->ROT_RIGHT->BackgroundImage = (cli::safe_cast<System::Drawing::Image^>(resources->GetObject(L"ROT_RIGHT.BackgroundImage")));
			this->ROT_RIGHT->BackgroundImageLayout = System::Windows::Forms::ImageLayout::Zoom;
			this->ROT_RIGHT->ForeColor = System::Drawing::SystemColors::ActiveCaption;
			this->ROT_RIGHT->Location = System::Drawing::Point(140, 509);
			this->ROT_RIGHT->Name = L"ROT_RIGHT";
			this->ROT_RIGHT->Size = System::Drawing::Size(30, 30);
			this->ROT_RIGHT->TabIndex = 7;
			this->ROT_RIGHT->UseVisualStyleBackColor = false;
			this->ROT_RIGHT->Visible = false;
			this->ROT_RIGHT->Click += gcnew System::EventHandler(this, &MyForm::ROT_RIGHT_Click);
			// 
			// ROT_LEFT
			// 
			this->ROT_LEFT->Anchor = static_cast<System::Windows::Forms::AnchorStyles>((System::Windows::Forms::AnchorStyles::Bottom | System::Windows::Forms::AnchorStyles::Left));
			this->ROT_LEFT->BackColor = System::Drawing::Color::SteelBlue;
			this->ROT_LEFT->BackgroundImage = (cli::safe_cast<System::Drawing::Image^>(resources->GetObject(L"ROT_LEFT.BackgroundImage")));
			this->ROT_LEFT->BackgroundImageLayout = System::Windows::Forms::ImageLayout::Zoom;
			this->ROT_LEFT->ForeColor = System::Drawing::SystemColors::ActiveCaption;
			this->ROT_LEFT->Location = System::Drawing::Point(104, 509);
			this->ROT_LEFT->Name = L"ROT_LEFT";
			this->ROT_LEFT->Size = System::Drawing::Size(30, 30);
			this->ROT_LEFT->TabIndex = 6;
			this->ROT_LEFT->UseVisualStyleBackColor = false;
			this->ROT_LEFT->Visible = false;
			this->ROT_LEFT->Click += gcnew System::EventHandler(this, &MyForm::ROT_LEFT_Click);
			// 
			// LEFT
			// 
			this->LEFT->Anchor = static_cast<System::Windows::Forms::AnchorStyles>((System::Windows::Forms::AnchorStyles::Bottom | System::Windows::Forms::AnchorStyles::Left));
			this->LEFT->BackColor = System::Drawing::Color::Bisque;
			this->LEFT->BackgroundImage = (cli::safe_cast<System::Drawing::Image^>(resources->GetObject(L"LEFT.BackgroundImage")));
			this->LEFT->BackgroundImageLayout = System::Windows::Forms::ImageLayout::Zoom;
			this->LEFT->ForeColor = System::Drawing::Color::Bisque;
			this->LEFT->Location = System::Drawing::Point(21, 509);
			this->LEFT->Name = L"LEFT";
			this->LEFT->Size = System::Drawing::Size(30, 30);
			this->LEFT->TabIndex = 5;
			this->LEFT->UseVisualStyleBackColor = false;
			this->LEFT->Visible = false;
			this->LEFT->Click += gcnew System::EventHandler(this, &MyForm::LEFT_Click);
			// 
			// RIGHT
			// 
			this->RIGHT->Anchor = static_cast<System::Windows::Forms::AnchorStyles>((System::Windows::Forms::AnchorStyles::Bottom | System::Windows::Forms::AnchorStyles::Left));
			this->RIGHT->BackColor = System::Drawing::Color::Bisque;
			this->RIGHT->BackgroundImage = (cli::safe_cast<System::Drawing::Image^>(resources->GetObject(L"RIGHT.BackgroundImage")));
			this->RIGHT->BackgroundImageLayout = System::Windows::Forms::ImageLayout::Zoom;
			this->RIGHT->ForeColor = System::Drawing::Color::Bisque;
			this->RIGHT->Location = System::Drawing::Point(21, 473);
			this->RIGHT->Name = L"RIGHT";
			this->RIGHT->Size = System::Drawing::Size(30, 30);
			this->RIGHT->TabIndex = 4;
			this->RIGHT->UseVisualStyleBackColor = false;
			this->RIGHT->Visible = false;
			this->RIGHT->Click += gcnew System::EventHandler(this, &MyForm::RIGHT_Click);
			// 
			// DOWN
			// 
			this->DOWN->Anchor = static_cast<System::Windows::Forms::AnchorStyles>((System::Windows::Forms::AnchorStyles::Bottom | System::Windows::Forms::AnchorStyles::Left));
			this->DOWN->BackColor = System::Drawing::Color::SteelBlue;
			this->DOWN->BackgroundImage = (cli::safe_cast<System::Drawing::Image^>(resources->GetObject(L"DOWN.BackgroundImage")));
			this->DOWN->BackgroundImageLayout = System::Windows::Forms::ImageLayout::Zoom;
			this->DOWN->ForeColor = System::Drawing::SystemColors::ActiveCaption;
			this->DOWN->Location = System::Drawing::Point(57, 509);
			this->DOWN->Name = L"DOWN";
			this->DOWN->Size = System::Drawing::Size(30, 30);
			this->DOWN->TabIndex = 3;
			this->DOWN->UseVisualStyleBackColor = false;
			this->DOWN->Visible = false;
			this->DOWN->Click += gcnew System::EventHandler(this, &MyForm::DOWN_Click);
			// 
			// UP
			// 
			this->UP->Anchor = static_cast<System::Windows::Forms::AnchorStyles>((System::Windows::Forms::AnchorStyles::Bottom | System::Windows::Forms::AnchorStyles::Left));
			this->UP->BackColor = System::Drawing::Color::SteelBlue;
			this->UP->BackgroundImage = (cli::safe_cast<System::Drawing::Image^>(resources->GetObject(L"UP.BackgroundImage")));
			this->UP->BackgroundImageLayout = System::Windows::Forms::ImageLayout::Zoom;
			this->UP->ForeColor = System::Drawing::SystemColors::ActiveCaption;
			this->UP->Location = System::Drawing::Point(57, 473);
			this->UP->Name = L"UP";
			this->UP->Size = System::Drawing::Size(30, 30);
			this->UP->TabIndex = 2;
			this->UP->UseVisualStyleBackColor = false;
			this->UP->Visible = false;
			this->UP->Click += gcnew System::EventHandler(this, &MyForm::UP_Click);
			// 
			// glControl1
			// 
			this->glControl1->Anchor = static_cast<System::Windows::Forms::AnchorStyles>((((System::Windows::Forms::AnchorStyles::Top | System::Windows::Forms::AnchorStyles::Bottom)
				| System::Windows::Forms::AnchorStyles::Left)
				| System::Windows::Forms::AnchorStyles::Right));
			this->glControl1->BackColor = System::Drawing::Color::Transparent;
			this->glControl1->Location = System::Drawing::Point(18, 23);
			this->glControl1->Margin = System::Windows::Forms::Padding(4, 4, 4, 4);
			this->glControl1->Name = L"glControl1";
			this->glControl1->Size = System::Drawing::Size(839, 519);
			this->glControl1->TabIndex = 1;
			this->glControl1->VSync = true;
			this->glControl1->Load += gcnew System::EventHandler(this, &MyForm::glControl1_Load);
			this->glControl1->Paint += gcnew System::Windows::Forms::PaintEventHandler(this, &MyForm::glControl1_Paint);
			this->glControl1->Resize += gcnew System::EventHandler(this, &MyForm::glControl1_Resize);
			// 
			// CHART
			// 
			this->CHART->Anchor = static_cast<System::Windows::Forms::AnchorStyles>((((System::Windows::Forms::AnchorStyles::Top | System::Windows::Forms::AnchorStyles::Bottom)
				| System::Windows::Forms::AnchorStyles::Left)
				| System::Windows::Forms::AnchorStyles::Right));
			this->CHART->BackColor = System::Drawing::Color::Transparent;
			chartArea1->Area3DStyle->LightStyle = System::Windows::Forms::DataVisualization::Charting::LightStyle::Realistic;
			chartArea1->AxisX->LabelStyle->ForeColor = System::Drawing::Color::WhiteSmoke;
			chartArea1->AxisX->LineColor = System::Drawing::Color::White;
			chartArea1->AxisX->MajorGrid->LineColor = System::Drawing::Color::Gray;
			chartArea1->AxisX->MajorTickMark->LineColor = System::Drawing::Color::Silver;
			chartArea1->AxisX->MinorGrid->LineColor = System::Drawing::Color::Silver;
			chartArea1->AxisX->MinorTickMark->LineColor = System::Drawing::Color::Silver;
			chartArea1->AxisX2->LineColor = System::Drawing::Color::Gray;
			chartArea1->AxisX2->MajorGrid->LineColor = System::Drawing::Color::DarkGray;
			chartArea1->AxisY->LabelStyle->ForeColor = System::Drawing::Color::WhiteSmoke;
			chartArea1->AxisY->LineColor = System::Drawing::Color::White;
			chartArea1->AxisY->MajorGrid->LineColor = System::Drawing::Color::Gray;
			chartArea1->AxisY->MajorTickMark->LineColor = System::Drawing::Color::Silver;
			chartArea1->AxisY->MinorGrid->LineColor = System::Drawing::Color::Silver;
			chartArea1->AxisY->MinorTickMark->LineColor = System::Drawing::Color::Silver;
			chartArea1->AxisY2->LineColor = System::Drawing::Color::DarkGray;
			chartArea1->AxisY2->MajorGrid->LineColor = System::Drawing::Color::DarkGray;
			chartArea1->BackColor = System::Drawing::Color::MidnightBlue;
			chartArea1->BackImageTransparentColor = System::Drawing::Color::Transparent;
			chartArea1->BackSecondaryColor = System::Drawing::Color::Transparent;
			chartArea1->Name = L"ChartArea1";
			this->CHART->ChartAreas->Add(chartArea1);
			legend1->Enabled = false;
			legend1->Name = L"Legend1";
			this->CHART->Legends->Add(legend1);
			this->CHART->Location = System::Drawing::Point(7, 29);
			this->CHART->Name = L"CHART";
			series1->ChartArea = L"ChartArea1";
			series1->ChartType = System::Windows::Forms::DataVisualization::Charting::SeriesChartType::Line;
			series1->Legend = L"Legend1";
			series1->Name = L"Series1";
			this->CHART->Series->Add(series1);
			this->CHART->Size = System::Drawing::Size(861, 527);
			this->CHART->TabIndex = 0;
			this->CHART->Text = L"chart1";
			// 
			// statusStrip1
			// 
			this->statusStrip1->BackColor = System::Drawing::Color::Transparent;
			this->statusStrip1->ImageScalingSize = System::Drawing::Size(19, 19);
			this->statusStrip1->Items->AddRange(gcnew cli::array< System::Windows::Forms::ToolStripItem^  >(1) { this->toolStripStatusLabel1 });
			this->statusStrip1->Location = System::Drawing::Point(0, 604);
			this->statusStrip1->Name = L"statusStrip1";
			this->statusStrip1->Size = System::Drawing::Size(1194, 22);
			this->statusStrip1->TabIndex = 7;
			this->statusStrip1->Text = L"statusStrip1";
			// 
			// toolStripStatusLabel1
			// 
			this->toolStripStatusLabel1->ForeColor = System::Drawing::SystemColors::ControlLight;
			this->toolStripStatusLabel1->Name = L"toolStripStatusLabel1";
			this->toolStripStatusLabel1->Size = System::Drawing::Size(118, 17);
			this->toolStripStatusLabel1->Text = L"toolStripStatusLabel1";
			// 
			// MyForm
			// 
			this->AutoScaleDimensions = System::Drawing::SizeF(96, 96);
			this->AutoScaleMode = System::Windows::Forms::AutoScaleMode::Dpi;
			this->BackColor = System::Drawing::Color::LightSkyBlue;
			this->BackgroundImage = (cli::safe_cast<System::Drawing::Image^>(resources->GetObject(L"$this.BackgroundImage")));
			this->BackgroundImageLayout = System::Windows::Forms::ImageLayout::Center;
			this->ClientSize = System::Drawing::Size(1194, 626);
			this->Controls->Add(this->statusStrip1);
			this->Controls->Add(this->groupBox7);
			this->Controls->Add(this->GET_MAX);
			this->Controls->Add(this->GET_MIN);
			this->Controls->Add(this->groupBox6);
			this->Controls->Add(this->groupBox4);
			this->Controls->Add(this->groupBox1);
			this->Controls->Add(this->menuStrip1);
			this->Icon = (cli::safe_cast<System::Drawing::Icon^>(resources->GetObject(L"$this.Icon")));
			this->MainMenuStrip = this->menuStrip1;
			this->MinimumSize = System::Drawing::Size(1210, 664);
			this->Name = L"MyForm";
			this->Text = L"Optimizer Visual 2.0";
			this->menuStrip1->ResumeLayout(false);
			this->menuStrip1->PerformLayout();
			this->groupBox1->ResumeLayout(false);
			this->groupBox1->PerformLayout();
			this->groupBox3->ResumeLayout(false);
			this->groupBox3->PerformLayout();
			this->groupBox2->ResumeLayout(false);
			this->groupBox2->PerformLayout();
			this->groupBox4->ResumeLayout(false);
			this->groupBox4->PerformLayout();
			this->groupBox5->ResumeLayout(false);
			this->groupBox5->PerformLayout();
			this->groupBox6->ResumeLayout(false);
			this->groupBox6->PerformLayout();
			this->groupBox7->ResumeLayout(false);
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->CHART))->EndInit();
			this->statusStrip1->ResumeLayout(false);
			this->statusStrip1->PerformLayout();
			this->ResumeLayout(false);
			this->PerformLayout();

		}
#pragma endregion

	private: System::Void выходToolStripMenuItem_Click(System::Object^  sender, System::EventArgs^  e) 
	{
		Application::Exit();
	}
			
private: System::Void textBox2_TextChanged(System::Object^  sender, System::EventArgs^  e) 
{
	int threads = 0;
	try
	{
		if (CPU_COUNT->Text == "")
		{
			threads = MAX_THREADS;
		}
		return;
		threads = Convert::ToInt32(CPU_COUNT->Text);
	}
	catch (...)
	{
		MessageBox::Show("Число потоков ЦП должно быть целочисленным", "Ошибка ввода", MessageBoxButtons::OK, MessageBoxIcon::Error);
		CPU_COUNT->Text = Convert::ToString(MAX_THREADS);
		return;
	}
	if (threads > MAX_THREADS)
	{
		MessageBox::Show("Число потоков ЦП должно быть меньше максимального = [ " + Convert::ToString(MAX_THREADS) + " ]!", "Ошибка ввода", MessageBoxButtons::OK, MessageBoxIcon::Error);
		CPU_COUNT->Text = Convert::ToString(MAX_THREADS);
		return;
	}
	if (threads <= 0)
	{
		MessageBox::Show("Число потоков ЦП должно быть больше нуля", "Ошибка ввода", MessageBoxButtons::OK, MessageBoxIcon::Error);
		CPU_COUNT->Text = Convert::ToString(MAX_THREADS);
		return;
	}

}
private: System::Void radioButton1_CheckedChanged(System::Object^  sender, System::EventArgs^  e) 
{
	if (OFF_OPTIMIZE->Checked)
	{
		groupBox5->Enabled = false;
		PP_BY_R_MAX->Checked = true;
	}
}

private: System::Void ON_OPTIMIZE_CheckedChanged(System::Object^  sender, System::EventArgs^  e) 
{
	if (ON_OPTIMIZE->Checked)
	{
		groupBox5->Enabled = true;
		PP_BY_R_MAX->Checked = true;
	}
}
private: System::Void button1_Click(System::Object^  sender, System::EventArgs^  e) 
{
	MessageBox::Show("Переменные вводятся в виде X и Y в случае двумерной функции или X в случае одномерной функции. \nАргументы элементарных функций необходимо заключать в скобки.",
		"Информация по вводу данных", MessageBoxButtons::OK, MessageBoxIcon::Information);
}
private: System::Void SET_EXPRESSION_Click(System::Object^  sender, System::EventArgs^  e) 
{
	Z_VAL_OPT->Text = "";
	X_VAL_OPT->Text = "";
	Y_VAL_OPT->Text = "";
	TIME_ELAPSED->Text = "";
	STEPS_ELAPSED->Text = "";
	move_x = move_y = move_z = 0;
	rotate_grad = 0;
	if (CheckDimension())
	{
		if (DIMENSION < 1)
		{
			MessageBox::Show("Плохой ввод функции. Проверьте ввод. Используйте инструкцию.", "Ошибка ввода", MessageBoxButtons::OK, MessageBoxIcon::Error);
			ready_for_use = false;
			return;
		}
		if (DIMENSION == 1)
		{
			label8->Enabled = false;
			Y_VAL_OPT->Enabled = false;
			Y_LEFT->Enabled = false;
			Y_RIGHT->Enabled = false;
			SHOW_TEST_POINTS->Visible = true;
			SHOW_TEST_POINTS->Enabled = false;
			PP_BY_INTERVALS->Enabled = true;
		}
		else if (DIMENSION == 2)
		{
			label8->Enabled = true;
			Y_VAL_OPT->Enabled = true;
			Y_LEFT->Enabled = true;
			Y_RIGHT->Enabled = true;
			SHOW_TEST_POINTS->Visible = false;
			SHOW_TEST_POINTS->Enabled = false;
			PP_BY_INTERVALS->Enabled = false;
		}
		else
		{
			SET_EXPRESSION->BackColor = Color::Red;
			MessageBox::Show("Плохой ввод функции. Проверьте ввод. Используйте инструкцию.", "Ошибка ввода", MessageBoxButtons::OK, MessageBoxIcon::Error);
			ready_for_use = false;
			return;
		}
		SET_EXPRESSION->BackColor = Color::LightGreen;
		ready_for_use = true;
		optimized = false;
		if (draw_graph)
		{
			std::string expression;
			msclr::interop::marshal_context context;
			expression = context.marshal_as<std::string>(MAIN_FUNC->Text);
			X_LEFT->Text->Replace('.', ',');
			Y_LEFT->Text->Replace('.', ',');
			X_RIGHT->Text->Replace('.', ',');
			Y_RIGHT->Text->Replace('.', ',');
			R_PARAM->Text->Replace('.', ',');

			double *Left = new double[2];
			double *Right = new double[2];

			Left[0] = Convert::ToDouble(X_LEFT->Text);
			Left[1] = Convert::ToDouble(Y_LEFT->Text);
			Right[0] = Convert::ToDouble(X_RIGHT->Text);
			Right[1] = Convert::ToDouble(Y_RIGHT->Text);

			if (DIMENSION == 1)
			{
				CHART->Visible = true;
				glControl1->Visible = false;
				ROT_LEFT->Visible = false;
				ROT_RIGHT->Visible = false;
				LEFT->Visible = false;
				RIGHT->Visible = false;
				UP->Visible = false;
				DOWN->Visible = false;
				if(draw_graph)
				Plot1D(expression, Left[0], Right[0], 300);
			}
			else if (DIMENSION == 2)
			{
				CHART->Visible = false;
				glControl1->Visible = true;
				ROT_LEFT->Visible = true;
				ROT_RIGHT->Visible = true;
				LEFT->Visible = true;
				RIGHT->Visible = true;
				UP->Visible = true;
				DOWN->Visible = true;
				glControl1->Invalidate();
				glControl1->Refresh();
				glControl1_Resize(sender, e);
			}
			statusStrip1->Items[0]->Text = "График функции отрисован";
		}
		
	}
}
private: System::Void button2_Click(System::Object^  sender, System::EventArgs^  e) 
{
	if (!ready_for_use || DIMENSION == 0)
	{
		MessageBox::Show("Необходима инициализация функции", "Информация", MessageBoxButtons::OK, MessageBoxIcon::Stop);
		return;
	}
	groupBox6->Enabled = false;
	std::string expression;
	msclr::interop::marshal_context context;
	expression = context.marshal_as<std::string>(MAIN_FUNC->Text);
	EPS->Text->Replace('.', ',');
	X_LEFT->Text->Replace('.', ',');
	Y_LEFT->Text->Replace('.', ',');
	X_RIGHT->Text->Replace('.', ',');
	Y_RIGHT->Text->Replace('.', ',');
	R_PARAM->Text->Replace('.', ',');
	Z_VAL_OPT->Text = "";
	X_VAL_OPT->Text = "";
	Y_VAL_OPT->Text = "";
	TIME_ELAPSED->Text = "";
	STEPS_ELAPSED->Text = "";
	double r = Convert::ToDouble(R_PARAM->Text);
	int N_max = 0;
	double epsilon = Convert::ToDouble(EPS->Text);
	if (r < 1)
	{
		MessageBox::Show("Некорректный параметр 'r' для метода.\n     r > 1", "Ошибка ввода", MessageBoxButtons::OK, MessageBoxIcon::Error);
		return;
	}
	try
	{
		N_max = Convert::ToInt32(STEPS->Text);
	}
	catch(...)
	{
		MessageBox::Show("Максимальное количество шагов должно быть целочисленным", "Ошибка ввода", MessageBoxButtons::OK, MessageBoxIcon::Error);
		return;
	}
	if (N_max < 1)
	{
		MessageBox::Show("Максимальное количество шагов должно быть больше 1", "Ошибка ввода", MessageBoxButtons::OK, MessageBoxIcon::Error);
		return;
	}
	double *Left = new double[2];
	double *Right = new double[2];
	Left[0] = Convert::ToDouble(X_LEFT->Text);
	Left[1] = Convert::ToDouble(Y_LEFT->Text);
	Right[0] = Convert::ToDouble(X_RIGHT->Text);
	Right[1] = Convert::ToDouble(Y_RIGHT->Text);

	bool is_parallel = ON_OPTIMIZE->Checked;
	int parallel_mode = 0;
	if (PP_BY_INTERVALS->Checked) parallel_mode = PPM_BY_DIVISIONS;
	else parallel_mode = PPM_BY_R_INTERVALS;

	int threads = Convert::ToInt32(CPU_COUNT->Text);
	SHOW_TEST_POINTS->Checked = false;
	while(CHART->Series->Count>1)
	CHART->Series->RemoveAt(1);
	
	Optimize(true, expression, Left, Right, r, epsilon, N_max, DIMENSION, is_parallel, parallel_mode, threads, *x_test_points, *y_test_points);
	optimized = true;
	groupBox6->Enabled = true;
	SHOW_TEST_POINTS->Enabled = true;
	omp_set_num_threads(MAX_THREADS);
	// костыльная защита
	CPU_COUNT->Text = Convert::ToString(threads);
	statusStrip1->Items[0]->Text = "Глобальный минимум найден";
}


private: System::Void SHOW_TEST_POINTS_CheckedChanged(System::Object^  sender, System::EventArgs^  e) 
{
	if (!draw_graph) return;
	if (SHOW_TEST_POINTS->Checked)
	{
		Plot1D_test_points(*x_test_points, *y_test_points);
	}
	else
	{
		CHART->Series->RemoveAt(2);
	}
	statusStrip1->Items[0]->Text = "Точки испытаний отрисованы";
}
private: System::Void MAIN_FUNC_TextChanged(System::Object^  sender, System::EventArgs^  e) 
{
	ready_for_use = false;
	SET_EXPRESSION->BackColor = Color::DarkGray;
	CHART->Visible = false;
	glControl1->Visible = false;
	statusStrip1->Items[0]->Text = "Ожидание инициализации функции...";
	optimized = false;
	X_VAL_OPT->Text = "";
	Y_VAL_OPT->Text = "";
	Z_VAL_OPT->Text = "";
	TIME_ELAPSED->Text = "";
	STEPS_ELAPSED->Text = "";
	DIMENSION = 0;

	ROT_LEFT->Visible = false;
	ROT_RIGHT->Visible = false;
	LEFT->Visible = false;
	RIGHT->Visible = false;
	UP->Visible = false;
	DOWN->Visible = false;
}
private: System::Void GET_MAX_Click(System::Object^  sender, System::EventArgs^  e) 
{
	if (!ready_for_use || DIMENSION == 0)
	{
		MessageBox::Show("Необходима инициализация функции", "Информация", MessageBoxButtons::OK, MessageBoxIcon::Stop);
		return;
	}
	groupBox6->Enabled = false;
	std::string expression;
	msclr::interop::marshal_context context;
	expression = context.marshal_as<std::string>(MAIN_FUNC->Text);
	expression.insert(0, "-1*(", strlen("-1*("));
	expression.insert(expression.end(), ')');

	EPS->Text->Replace('.', ',');
	X_LEFT->Text->Replace('.', ',');
	Y_LEFT->Text->Replace('.', ',');
	X_RIGHT->Text->Replace('.', ',');
	Y_RIGHT->Text->Replace('.', ',');
	R_PARAM->Text->Replace('.', ',');
	Z_VAL_OPT->Text = "";
	X_VAL_OPT->Text = "";
	Y_VAL_OPT->Text = "";
	TIME_ELAPSED->Text = "";
	STEPS_ELAPSED->Text = "";
	double r = Convert::ToDouble(R_PARAM->Text);
	int N_max = 0;
	double epsilon = Convert::ToDouble(EPS->Text);
	if (r < 1)
	{
		MessageBox::Show("Некорректный параметр 'r' для метода.\n     r > 1", "Ошибка ввода", MessageBoxButtons::OK, MessageBoxIcon::Error);
		return;
	}
	try
	{
		N_max = Convert::ToInt32(STEPS->Text);
	}
	catch (...)
	{
		MessageBox::Show("Максимальное количество шагов должно быть целочисленным", "Ошибка ввода", MessageBoxButtons::OK, MessageBoxIcon::Error);
		return;
	}
	if (N_max < 1)
	{
		MessageBox::Show("Максимальное количество шагов должно быть больше 1", "Ошибка ввода", MessageBoxButtons::OK, MessageBoxIcon::Error);
		return;
	}
	double *Left = new double[2];
	double *Right = new double[2];
	Left[0] = Convert::ToDouble(X_LEFT->Text);
	Left[1] = Convert::ToDouble(Y_LEFT->Text);
	Right[0] = Convert::ToDouble(X_RIGHT->Text);
	Right[1] = Convert::ToDouble(Y_RIGHT->Text);

	bool is_parallel = ON_OPTIMIZE->Checked;
	int parallel_mode = 0;
	if (PP_BY_INTERVALS->Checked) parallel_mode = PPM_BY_DIVISIONS;
	else parallel_mode = PPM_BY_R_INTERVALS;

	int threads = Convert::ToInt32(CPU_COUNT->Text);
	SHOW_TEST_POINTS->Checked = false;
	while (CHART->Series->Count>1)
		CHART->Series->RemoveAt(1);

	Optimize(false, expression, Left, Right, r, epsilon, N_max, DIMENSION, is_parallel, parallel_mode, threads, *x_test_points, *y_test_points);
	optimized = true;
	groupBox6->Enabled = true;
	SHOW_TEST_POINTS->Enabled = true;
	omp_set_num_threads(MAX_THREADS);
	// костыльная защита
	CPU_COUNT->Text = Convert::ToString(threads);
	statusStrip1->Items[0]->Text = "Глобальный максимум найден";
}
private: System::Void glControl1_Paint(System::Object^  sender, System::Windows::Forms::PaintEventArgs^  e) 
{
	if (draw_graph && ready_for_use && DIMENSION == 2)
	{
		statusStrip1->Items[0]->Text = "Отрисовка графика функции...";
		std::string expression;
		msclr::interop::marshal_context context;
		expression = context.marshal_as<std::string>(MAIN_FUNC->Text);
		X_LEFT->Text->Replace('.', ',');
		Y_LEFT->Text->Replace('.', ',');
		X_RIGHT->Text->Replace('.', ',');
		Y_RIGHT->Text->Replace('.', ',');
		R_PARAM->Text->Replace('.', ',');

		double *Left = new double[2];
		double *Right = new double[2];

		Left[0] = Convert::ToDouble(X_LEFT->Text);
		Left[1] = Convert::ToDouble(Y_LEFT->Text);
		Right[0] = Convert::ToDouble(X_RIGHT->Text);
		Right[1] = Convert::ToDouble(Y_RIGHT->Text);

		Plot2D(expression, Left, Right, Convert::ToInt32(sqrt(pow(abs(Left[0]),2) + pow(abs(Right[0]),2)) + sqrt(pow(abs(Left[1]),2) + pow(abs(Right[1]),2))) * 2);
		
		delete[]Left;
		delete[]Right;
		statusStrip1->Items[0]->Text = "Готов к использованию...";
	}
}
private: System::Void glControl1_Load(System::Object^  sender, System::EventArgs^  e) 
{
	if (!draw_graph || !ready_for_use) return;
	X_LEFT->Text->Replace('.', ',');
	Y_LEFT->Text->Replace('.', ',');
	X_RIGHT->Text->Replace('.', ',');
	Y_RIGHT->Text->Replace('.', ',');
	R_PARAM->Text->Replace('.', ',');

	double *Left = new double[2];
	double *Right = new double[2];

	Left[0] = Convert::ToDouble(X_LEFT->Text);
	Left[1] = Convert::ToDouble(Y_LEFT->Text);
	Right[0] = Convert::ToDouble(X_RIGHT->Text);
	Right[1] = Convert::ToDouble(Y_RIGHT->Text);

	double x = (Right[0] - Left[0]) / 2;
	double y = (Right[1] - Left[1]) / 2;

	GL::ClearColor(Color::MidnightBlue);
	GL::Enable(EnableCap::DepthTest);

	//OpenTK::Matrix4 p = OpenTK::Matrix4::CreatePerspectiveFieldOfView((float)(80 * PI / 180), 1, 20, 500);
	OpenTK::Matrix4 p = OpenTK::Matrix4::CreatePerspectiveFieldOfView((float)(50 * PI / 180), (float)(glControl1->Width / glControl1->Height), (float)1, (float)1000);
	GL::MatrixMode(MatrixMode::Projection);
	GL::LoadMatrix(p);

	//OpenTK::Matrix4 modelview = OpenTK::Matrix4::LookAt(70, 70, 70, 0, 0, 0, 0, 1, 0);
	//OpenTK::Matrix4 modelview = OpenTK::Matrix4::LookAt(x*2.5 + y/2, y*2.5 + x/2, (x + y) * 0.5*pow((x/2 + y/2),2), x, y, 0, 0, 0, 1);
	//OpenTK::Matrix4 modelview = OpenTK::Matrix4::LookAt(x*2.5 + y, y*2.5 + x, (x + y) * 3, x, y, (x+y), 0, 0, 1);
	OpenTK::Matrix4 modelview = OpenTK::Matrix4::LookAt(x*2.5 + y + move_x, y*2.5 + x + move_y, (x + y) * 3 + move_z, x, y, abs(sqrt((pow(x * 2, 2) + pow(y * 2, 2)))) / 2 + (x + y) / 2, 0, 0, abs(sqrt((pow(x * 2, 2) + pow(y * 2, 2)))) / 2 + (x + y) / 2);
	GL::MatrixMode(MatrixMode::Modelview);
	GL::LoadMatrix(modelview);
	Rotate(rotate_grad);

	glControl1->Invalidate();
	glControl1->Refresh();

	delete[]Left;
	delete[]Right;

}
private: System::Void glControl1_Resize(System::Object^  sender, System::EventArgs^  e) 
{
	if (!draw_graph || !ready_for_use) return;
	//if (!ready_for_use) return;
	statusStrip1->Items[0]->Text = "Отрисовка графика функции...";
	X_LEFT->Text->Replace('.', ',');
	Y_LEFT->Text->Replace('.', ',');
	X_RIGHT->Text->Replace('.', ',');
	Y_RIGHT->Text->Replace('.', ',');
	R_PARAM->Text->Replace('.', ',');

	double *Left = new double[2];
	double *Right = new double[2];

	Left[0] = Convert::ToDouble(X_LEFT->Text);
	Left[1] = Convert::ToDouble(Y_LEFT->Text);
	Right[0] = Convert::ToDouble(X_RIGHT->Text);
	Right[1] = Convert::ToDouble(Y_RIGHT->Text);

	double x = abs((Right[0] - Left[0]) / 2);
	double y = abs((Right[1] - Left[1]) / 2);

	//OpenTK::Matrix4 p = OpenTK::Matrix4::CreatePerspectiveFieldOfView((float)(80 * PI / 180), 1, 20, 500);
	OpenTK::Matrix4 p = OpenTK::Matrix4::CreatePerspectiveFieldOfView((float)(50 * PI / 180), (float)(glControl1->Width / glControl1->Height), (float)1, (float)1000);
	GL::MatrixMode(MatrixMode::Projection);
	GL::LoadMatrix(p);

	//OpenTK::Matrix4 modelview = OpenTK::Matrix4::LookAt(70, 70, 70, 0, 0, 0, 0, 1, 0);
	//OpenTK::Matrix4 modelview = OpenTK::Matrix4::LookAt(x*2.5 + y, y*2.5 + x, (x + y) * 3, x, y, (x+y), 0, 0, 1);
	OpenTK::Matrix4 modelview = OpenTK::Matrix4::LookAt(x*2.5 + y + move_x, y*2.5 + x + move_y, (x + y) * 3 + move_z, x, y, abs(sqrt((pow(x * 2, 2) + pow(y * 2, 2)))) / 2 + (x + y) / 2, 0, 0, abs(sqrt((pow(x * 2, 2) + pow(y * 2, 2)))) / 2 + (x + y) / 2);
	GL::MatrixMode(MatrixMode::Modelview);
	GL::LoadMatrix(modelview);
	Rotate(rotate_grad);

	glControl1->Invalidate();
	glControl1->Refresh();

	delete[]Left;
	delete[]Right;
	statusStrip1->Items[0]->Text = "Готов к использованию...";
}
private: System::Void CAN_DRAW_CheckedChanged(System::Object^  sender, System::EventArgs^  e) 
{
	if (CAN_DRAW->Checked)
	{
		draw_graph = true;
		if (DIMENSION == 1)
		{
			CHART->Visible = true;

			ROT_LEFT->Visible = false;
			ROT_RIGHT->Visible = false;
			LEFT->Visible = false;
			RIGHT->Visible = false;
			UP->Visible = false;
			DOWN->Visible = false;
		}
		else if (DIMENSION == 2) 
		{
			glControl1->Visible = true;

			ROT_LEFT->Visible = true;
			ROT_RIGHT->Visible = true;
			LEFT->Visible = true;
			RIGHT->Visible = true;
			UP->Visible = true;
			DOWN->Visible = true;
		}
	}
	else
	{
		draw_graph = false;
		CHART->Visible = false;
		glControl1->Visible = false;

		CHART->Visible = false;

		glControl1->Visible = false;
		ROT_LEFT->Visible = false;
		ROT_RIGHT->Visible = false;
		LEFT->Visible = false;
		RIGHT->Visible = false;
		UP->Visible = false;
		DOWN->Visible = false;
	}
}
private: System::Void оПрограммеToolStripMenuItem_Click(System::Object^  sender, System::EventArgs^  e) 
{
	MessageBox::Show("Программа для демонстрации алгоритма Стронгина Р.Г. поиска глобального минимума (АГП).\n\n" +
		"- Поиск глобального минимума для одномерных и двумерных функций методом редукции размерности через развертки Пеано плотностью m = 20,\n"
		"- Поддержка многопроцессорных систем с общей памятью,\n" + 
		"- Отрисовка 2D и 3D графиков для одномерных и двумерных функций.\n\n" +
		"Программа выполнена в рамках Курсового проекта.\n" + 
		"Программа и реализованные алгоритмы являются интеллектуальной собственостью автора.\n\n" +
		"Горб Н.Н. | 2018 | ННГУ им. Лобачевского | ИИТММ | ПМИ | каф. ТУиДС\n" + 
		"Науч.руководитель: \nСысоев А.В. доцент каф. МОСТ, к.т.н.,\n" + 
		"Городецкий С.Ю. доцент каф. ТУиДС, к.ф.-м.н.", "О программе", MessageBoxButtons::OK, MessageBoxIcon::Information);
}
private: System::Void вводФункцийToolStripMenuItem_Click(System::Object^  sender, System::EventArgs^  e) 
{
	MessageBox::Show("Функции вводятся в строковом представлении.\n\n * Список поддерживаемых элементарных функций:\n" +
		"- sin, cos, arcsin, arccos, tg, ctg, arctg, arcctg;\n" +
		"- sh, ch, th, cth;\n" +
		"- lg, ln, e^();\n\n" +
		" * Поддерживаемые математические операции:\n" +
		"- '+' '-' '*' '/';\n" +
		"- '%', 'sqrt', '^';\n\n" +
		" * Встроенные константы:\n" + 
		"- g = 9.81, pi = 3.14, e = 2.71."		, "Справка", MessageBoxButtons::OK, MessageBoxIcon::Information);

}
private: System::Void UP_Click(System::Object^  sender, System::EventArgs^  e)
{
	move_z -= 2;
	MoveCamera(move_x, move_y, move_z);
	Rotate(rotate_grad);
}
private: System::Void DOWN_Click(System::Object^  sender, System::EventArgs^  e)
{
	move_z += 2;
	MoveCamera(move_x, move_y, move_z);
	Rotate(rotate_grad);
}
private: System::Void RIGHT_Click(System::Object^  sender, System::EventArgs^  e)
{
	move_x -= 1;
	move_y -= 1;
	MoveCamera(move_x, move_y, move_z);
	Rotate(rotate_grad);
}
private: System::Void LEFT_Click(System::Object^  sender, System::EventArgs^  e)
{
	move_x += 1;
	move_y += 1;
	MoveCamera(move_x, move_y, move_z);
	Rotate(rotate_grad);
}
private: System::Void ROT_LEFT_Click(System::Object^  sender, System::EventArgs^  e)
{
	double rotate = -7;
	rotate_grad -= rotate;
	Rotate(rotate);
}
private: System::Void ROT_RIGHT_Click(System::Object^  sender, System::EventArgs^  e)
{
	double rotate = 7;
	rotate_grad += rotate;
	Rotate(rotate);
}

};//END OF FUNCTIONS
}
