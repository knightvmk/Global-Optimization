#ifndef __MPARCER_H__
#define __MPARCER_H__

#define _CRT_SECURE_NO_WARNINGS
#pragma warning (disable: 4996)

#include <ios>
#include <string>
#include <sstream>
#include <iomanip>

//using namespace std;

// Defines
#define MAX_EXPR_LEN   255
#define MAX_TOKEN_LEN  80
#define MAX_DOUBLE_DECISION 16

struct TCALCNode
{
	double value;

	TCALCNode *left;

	TCALCNode *right;

	TCALCNode(double _value = 0.0, TCALCNode *_left = NULL, TCALCNode *_right = NULL)
	{
		value = _value;
		left = _left;
		right = _right;
	}
};

struct TError
{
	char *error;
	int pos;

	TError() {	};

	TError(char *_error, int _pos)
	{
		error = _error;
		pos = _pos + 1;
	}
};

class TCALC
{
private:

	TCALCNode *root;
	char *expr;
	char curToken[MAX_TOKEN_LEN];

	enum {
		CALC_PLUS, CALC_MINUS, CALC_MULTIPLY, CALC_DIVIDE, CALC_PERCENT, CALC_POWER,
		CALC_SIN, CALC_COS, CALC_TG, CALC_CTG, CALC_ARCSIN, CALC_ARCCOS, CALC_ARCTG, 
		CALC_ARCCTG, CALC_SH, CALC_CH, CALC_TH, CALC_CTH, CALC_EXP, CALC_LG, CALC_LN, 
		CALC_SQRT, CALC_X, CALC_L_BRACKET, CALC_R_BRACKET, CALC_E, CALC_PI, CALC_NUMBER, 
		CALC_END, CALC_G, CALC_EXP1, CALC_LEET
	} typToken; //используемые действия, константы

	int pos;
	const double *x;
	double result;



private:

	inline TCALCNode *CreateNode(double _value = 0.0, TCALCNode *_left = NULL, TCALCNode  *_right = NULL);

	inline TCALCNode *Expr(void);

	inline TCALCNode *Expr1(void);

	inline TCALCNode *Expr2(void);

	inline TCALCNode *Expr3(void);

	inline TCALCNode *Expr4(void);

	inline TCALCNode *Expr5(void);


	inline bool GetToken(void);

	inline bool IsDelim(void)
	{
		return (strchr("+-*/%^()[]", expr[pos]) != NULL);    //разрешенные символы 
	}

	inline bool IsLetter(void)
	{
		return ((expr[pos] >= 'a' && expr[pos] <= 'z') ||
			(expr[pos] >= 'A' && expr[pos] <= 'Z'));
	}

	inline bool IsDigit(void)
	{
		return (expr[pos] >= '0' && expr[pos] <= '9');       // разрешенные цифры
	}

	inline bool IsPoint(void)
	{
		return (expr[pos] == '.');                         //и точка
	}

	inline double CalcTree(TCALCNode *tree);

	inline void DelTree(TCALCNode *tree);

	inline void SendError(int errNum);


public:

	TCALC()
	{
		result = 0.0;
		expr = new char[2];
		x = NULL;
		root = NULL;
	}

	~TCALC()
	{
		DelTree(root);
		root = NULL;
	}

	void SetX(const double *_x)
	{
		x = _x;
	}

	bool Compile(char *expr, double _x, double _y);

	void Decompile()
	{
		DelTree(root);
		root = NULL;
	}

	double Evaluate();

	double Evaluate(double *_x)
	{
		SetX(_x);
		return Evaluate();
	}

	double Evaluate(double x, ...)
	{
		SetX(&x);
		return Evaluate();
	}

	double GetResult(void)
	{
		return result;
	}
};


#endif
