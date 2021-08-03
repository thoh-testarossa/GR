#pragma once
#include<iostream>
#include<string>
#include<vector>
using namespace std;

class Data_instance
{
public:
	int Data_object_id = -1;//object id
	double Data_instance_prob;//prob
	double Data_instance_Measure;//Measure(y值)
	double Data_instance_x;//x值
	vector<string> Standard_Dimension;//标准维度
	vector<string> Time_Dimension;//时间维度-年-月
	Data_instance() {};
	Data_instance(int id,double x,double y,double p, vector<string> &S_D, vector<string> &T_D);

};
Data_instance::Data_instance(int id, double x, double y, double p, vector<string> &S_D, vector<string> &T_D)
{
	this->Data_object_id = id;
	this->Data_instance_prob = p;
	this->Data_instance_x = x;
	this->Data_instance_Measure = y;
	this->Standard_Dimension.assign(S_D.begin(),S_D.end());
	this->Time_Dimension.assign(T_D.begin(), T_D.end());
};


struct instance_pmf
{
	int id;
	double x, y, p;
	instance_pmf(int id,double x,double y,double p) {
		this->id =id;
		this->x = x;
		this->y = y;
		this->p = p;
	}
};
struct Object_pmfs {

	int id;
	bool isTruncated_pmf = false;
	vector<instance_pmf>pmfs;
	Object_pmfs(instance_pmf pmf) {
		id = pmf.id;
		pmfs.push_back(pmf);
	}
	bool isTruncated()//判断是否是截断pmf
	{
		double sum = 0.0;
		for (auto pmf : pmfs)
		{
			sum += pmf.p;
		}
		if (sum < 1) { 
			this->isTruncated_pmf = true;
			return true; }
		else {
			this->isTruncated_pmf = false;
			return false;
		}
		
	}
	void add_new_pmfs(Object_pmfs pmf)
	{
		if (this->id != pmf.id)
		{
			cout << "wrong!!"<< endl;
			return;
		}
		this->pmfs.insert(this->pmfs.end(), pmf.pmfs.begin(), pmf.pmfs.end());
		return;
	}
	void Make_Truncated_complete()
	{
		double sum = 0.0;
		for (auto pmf : pmfs)
		{
			sum += pmf.p;
		}
		double p = 1.0 - sum;
		instance_pmf com_pmf(this->id, 0, 0, p);
		this->pmfs.push_back(com_pmf);
	
}

};
struct Sketch
{
public:
	double E_k=0.0, E_b=0.0, E_kb=0.0, Var_k=0.0, Var_b=0.0,Cov=0.0;
	Sketch() {};
	Sketch(double ek,double eb,double ekb,double vark,double varb) {
		this->E_k = ek;
		this->E_b = eb;
		this->E_kb = ekb;
		this->Var_k = vark;
		this->Var_b = varb;
		this->Cov = this->E_kb-this->E_b*this->E_k;
	};
};
struct test_struct
{
	int number = 0;
	double time = 0;
	double space = 0;
	test_struct() {};
	test_struct(int num, double t, double s)
	{
		this->number = num;
		this->time = t;
		this->space = s*3*8/1024;//大小为KB

	}








};