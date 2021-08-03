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
	//double p[120][350] = { 0 };//���ؾࡿ��б�ʡ�
	//double time_p[][]
	//vector<int> b_search, k_search;
	int N = 1;//����ǰ�ľۺϸ�����1��ʱ����������cell�� k��Χ��(-0.6,0.5),b��Χ��(0,11),2��ʱ������k�ķ�ΧΪ(-0.12,0.1),b��(0,22)
	vector<int> Now;
	bool isTime = false;
	void add_new_line(Line line)
	{
		//line.k_r = atan(line.k);//���k�Ļ���
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
		//line.k_r = atan(line.k);//���k�Ļ���
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
		line.k_r = atan(line.k);//���k�Ļ���
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
	//�����˸����Ժ����ӷ���,�Ѹ��Ӷ�����
	/*void add_line_lattice(Line line)
	{
		
		if (this->isTime == false) {
		if (this->pmf.size() < 1)//�����ǰû���ߣ��򻭳����и���
		{
			Line A;
			this->pmf = vector<Line>(10000, A);
		}
		double k_radian = atan(line.k);//���k�Ļ���
		int k_search = 0;
		double k_s = k_radian + 0.6;
		k_search = int(k_s * 100);
		if (line.b < 0)
			line.b = 0;

		//Ѱ����������� search=b*11*N+k;
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
	//�������ָ��ӣ���һ����0-100�Ժ�б�ʾ�ȷ��[0.01����],��ȡֵ��-1.58����1.59��ȡС�������λ��bȡֵ��0-115��������ε��ӣ���b�ĸ��ӷŴ�һ��
	/*void add_line_lattice(Line line)
	{

		 
			
			double k_radian = atan(line.k);//���k�Ļ���
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
