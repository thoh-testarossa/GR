#pragma once
#include<iostream>
#include<string>
#include<vector>
#include"data_instance.h"
#include"data_object_pmf.h"
#include<math.h>
using namespace std;
static int number =0;
class Line {
public:
	double k=0.0, b=0.0, p=0.0,k_r=0.0;
	Line(double k, double b, double p)
	{
		//int k_i = k * pow(10, number);
		//int b_i = b * pow(10, number);
		this->k = k;
		this->b = b;
		this->p = p;
	}
	Line()
	{
		this->p = -0.1;
	
	
	}
};
class Lines_pmf {
public:
	vector<Line> pmf;
	//double p[120][350] = { 0 };//【截距】【斜率】
	//double time_p[][]
	//vector<int> b_search, k_search;
	int N = 1;//代表当前的聚合个数，1的时候就是最基础cell的 k范围是(-0.6,0.5),b范围是(0,11),2的时候则是k的范围为(-0.12,0.1),b是(0,22)
	vector<int> Now;
	bool isTime = false;
	void add_new_line(Line line)
	{
		//line.k_r = atan(line.k);//获得k的弧度
		//line.k_r = int(line.k_r * pow(10, 1));
		//line.k_r = line.k_r / pow(10, 1);
		int k_i = line.k * pow(10, number);
		int b_i = line.b * pow(10, number);
		line.k = double(k_i) / pow(10, number);
		line.b = double(b_i) / pow(10, number);
		for (auto& line1 : pmf)
		{
			if (line1.k == line.k && line1.b == line.b)
			{
				line1.p += line.p;
				return;
			}
		}
		pmf.push_back(line);
	}
	void add_new_line_k(Line line)
	{
		//line.k_r = atan(line.k);//获得k的弧度
		//line.k_r = int(line.k_r * pow(10, 1));
		//line.k_r = line.k_r / pow(10, 1);
		int k_i = line.k * pow(10, 0);
		int b_i = int(line.b);
		//b_i = b_i * 2;
		line.k = double(k_i) / pow(10, 0);
		line.b = b_i;
		for (auto& line1 : pmf)
		{
			if (line1.k == line.k && line1.b == line.b)
			{
				line1.p += line.p;
				return;
			}
		}
		pmf.push_back(line);
	}
	/*void add_new_linetest(Line line)
	{
		line.k_r = atan(line.k);//获得k的弧度
		line.k_r = int(line.k_r * pow(10, 1));
		line.k_r = line.k_r / pow(10, 1);
		int k_i = line.k * pow(10, number);
		int b_i = line.b * pow(10, number);
		line.k = double(k_i) / pow(10, number);
		line.b = double(b_i) / pow(10, number);
		for (auto& line1 : pmf)
		{
			if (line1.k_r == line.k_r && line1.b == line.b)
			{
				line1.p += line.p;
				return;
			}
		}
		pmf.push_back(line);
	}*/
	void add_new_line2(Line line)
	{
		int number1 = 0;
		int k_i = line.k * pow(10, number1)+5;
		int b_i = line.b * pow(10, number1)+5;
		line.k = double(k_i) / pow(10, number1);
		line.b = double(b_i) / pow(10, number1);
		for (auto& line1 : pmf)
		{
			if (line1.k == line.k && line1.b == line.b)
			{
				line1.p += line.p;
				return;
			}
		}
		pmf.push_back(line);
	}
	void add_new_line3(Line line)
	{
		int number1 = 0;
		int k_i = line.k * pow(10, number1);
		int b_i = line.b * pow(10, number1);
		line.k = double(k_i) / pow(10, number1);
		line.b = double(b_i) / pow(10, number1);
		for (auto& line1 : pmf)
		{
			if (line1.k == line.k && line1.b == line.b)
			{
				line1.p += line.p;
				return;
			}
		}
		pmf.push_back(line);
	}
	void add_new_line_pws(Line line)
	{
		int Number_pws =1;
		int k_i = line.k * pow(10, Number_pws);
		int b_i = line.b * pow(10, Number_pws);
		line.k = double(k_i) / pow(10, Number_pws);
		line.b = double(b_i) / pow(10, Number_pws);
		for (auto& line1 : pmf)
		{
			
			if (line1.k == line.k && line1.b == line.b)
			{
				line1.p += line.p;
				return;
			}
		}
		pmf.push_back(line);
	}
	//划分了格子以后的添加方法,把格子定量！
	/*void add_line_lattice(Line line)
	{
		
		if (this->isTime == false) {
		if (this->pmf.size() < 1)//如果当前没有线，则画出所有格子
		{
			Line A;
			this->pmf = vector<Line>(10000, A);
		}
		double k_radian = atan(line.k);//获得k的弧度
		int k_search = 0;
		double k_s = k_radian + 0.6;
		k_search = int(k_s * 100);
		if (line.b < 0)
			line.b = 0;

		//寻找所需格子数 search=b*11*N+k;
		int search = int(line.b*10) * 110*this->N + k_search;
		if (search > this->pmf.size())
		{
			search -= 49999;
			return;
		}
		if (this->pmf[search].p < 0)
		{
			this->pmf[search].k = line.k;
			this->pmf[search].b = int(line.b);
			this->pmf[search].p = line.p;
			this->pmf[search].k_r = k_radian;
			//this->N = 1;
			Now.push_back(search);
		}
		else
		{
			this->pmf[search].p += line.p;
		}
	
	}
	}
	*/
	//定量划分格子，归一化到0-100以后，斜率精确到[0.01弧度],其取值是-1.58到正1.59，取小数点后两位，b取值是0-115，如果两次叠加，则b的格子放大一倍
	/*void add_line_lattice(Line line)
	{

		 
			
			double k_radian = atan(line.k);//获得k的弧度
			int k_search = 0;
			double k_s = k_radian + 1.58;
			k_search = int(k_s * 10);
			if (line.b < 0)
				line.b = 0;
			int b_search = int(line.b);
			b_search = b_search / this->N;
			if (this->p[b_search][k_search] <= 0)
			{
				this->p[b_search][k_search] = line.p;
				this->k_search.push_back(k_search);
				this->b_search.push_back(b_search);
				this->pmf.push_back(line);
			}
			else
			{
				this->p[b_search][k_search]+= line.p;
				
			}
			
		
	}
	*/

};
