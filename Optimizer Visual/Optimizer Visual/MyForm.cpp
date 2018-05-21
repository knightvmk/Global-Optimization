#include "MyForm.h"

using namespace System;
using namespace System::Windows::Forms;

//-1/(2*sin(3*X)+3*cos(5*X))
//2*sin(3*X)+3*cos(5*X)
//2 * cos(Y) + 3 * sin(X)
//array<String^>^ args
[STAThread]
void Main(array<String^>^ args)
{
	Application::EnableVisualStyles();
	Application::SetCompatibleTextRenderingDefault(false);

	Optimizer_Visual::MyForm form;
	Application::Run(%form);
}
