#pragma once
#include"Lines.h"
#include"data_object_pmf.h"
//叶结点中计算fc的函数
void use_objects_to_get_pmfs(double T, Object_pmfs O_i_pmfs, Lines_pmf &F_x)
{
	//当前pmf为空
	if (O_i_pmfs.pmfs.size() == 0)
		return;
	double A_T = double(pow(T, 3) - T) / (12.0);
	double avg_x = double(T + 1) / (2.0);
	if (F_x.pmf.size() ==0)
	{//当前pmf是第一个pmf，视为{(1,y),(2,0),(T,0),p}构造第一根线
		for (auto pmf : O_i_pmfs.pmfs)
		{
			double sum_x_i_y_i = pmf.x * pmf.y;
			double avg_x_y = avg_x * pmf.y;
			double k_i = (sum_x_i_y_i - avg_x_y) / A_T;
			double b_i = pmf.y - k_i * avg_x;
			double p_i = pmf.p;
			Line l(k_i, b_i, p_i);
			//F_x.add_new_line(l);
			F_x.add_new_line_k(l);
		}
		return;
	}
	else {
		//把pmf累加到已经聚合出来的线上
		int size = O_i_pmfs.pmfs.size();//当前pmf由几部分组成
		if (size <= 0)
			return;
		//if (size > 1)
		{
			Lines_pmf re_pmf;
			//每一个instance都对现有的所有line进行修正
			for (int i = 0; i < O_i_pmfs.pmfs.size(); i++)
			{//累加公式为k=（xi*yi-avg(x)yi）/A(T)
				for (auto line : F_x.pmf)
				{
					double sum_x_i_y_i = O_i_pmfs.pmfs[i].x * O_i_pmfs.pmfs[i].y;
					double avgx_y = avg_x * O_i_pmfs.pmfs[i].y;
					double k_1 = line.k;
					double b_1 = line.b;
					double avg_y_1 = k_1 * avg_x + b_1;
					double p_2 = line.p * O_i_pmfs.pmfs[i].p;//更新概率
					double k_2 = k_1 + (sum_x_i_y_i - avgx_y) / A_T;//更新斜率
					double avg_y_2 = avg_y_1 + O_i_pmfs.pmfs[i].y / T;
					double b_2 = avg_y_2 - k_2 * avg_x;
					Line l(k_2, b_2, p_2);
					//re_pmf.add_new_line(l);
					re_pmf.add_new_line_k(l);
				}
			}
			F_x.pmf.clear();
			F_x.pmf.assign(re_pmf.pmf.begin(), re_pmf.pmf.end());

		}
		/*else if (size ==1)//当下只有一个pmf，其概率为1,则直接在fx里更新即可
		{
			for (auto line : F_x.pmf)
			{
				double sum_x_i_y_i = O_i_pmfs.pmfs[0].x * O_i_pmfs.pmfs[0].y;
				double avgx_y = avg_x * O_i_pmfs.pmfs[0].y;
				double k_1 = line.k;
				double b_1 = line.b;
				double avg_y_1 = k_1 * avg_x + b_1;
				line.p = line.p * O_i_pmfs.pmfs[0].p;//更新概率
				line.k = k_1 + (sum_x_i_y_i - avgx_y) / A_T;//更新斜率
				double avg_y_2 = avg_y_1 + O_i_pmfs.pmfs[0].y / T;
				line.b = avg_y_2 - line.k * avg_x;
			}
		}*/
		return;
	}



	return;


}
//标准维度上的卷积函数
Lines_pmf Convolution_at_Standard(Lines_pmf F_1, Lines_pmf F_2) {
	if (F_1.pmf.size() == 0)
		return F_2;
	if (F_2.pmf.size() == 0)
		return F_1;

	Lines_pmf re_pmf;
	for (auto line1:F_1.pmf)
		for (auto line2 : F_2.pmf)
		{
			double k = line1.k + line2.k;
			double b = line1.b + line2.b;
			double p = line1.p * line2.p;
			Line l(k, b, p);
			re_pmf.add_new_line(l);
		}
	return re_pmf;
}
//时间维度上的卷积函数
Lines_pmf Convolution_at_Time(Lines_pmf F_1, Lines_pmf F_2,double T_1,double T_2) {
	if (F_1.pmf.size() == 0)
		return F_2;
	if (F_2.pmf.size() == 0)
		return F_1;
	Lines_pmf re_pmf;
	double T_3 = T_1 + T_2;
	double avg_x_1 = (1 + T_1) / 2;
	double avg_x_2 = (1 + T_2) / 2;
	double avg_x_3 = (1 + T_3) / 2;

	for (auto line1 : F_1.pmf)
		for (auto line2 : F_2.pmf)
		{
			double avg_y_1 = line1.k * avg_x_1 + line1.b;
			double avg_y_2 = line2.k * avg_x_2 + line2.b;
			double avg_y_3 = (T_1 * avg_y_1 + T_2 * avg_y_2) / T_3;
			double p = line1.p * line2.p;
			double k = ((pow(T_1, 3) - T_1) / (pow(T_3, 3) - T_3)) * line1.k + ((pow(T_2, 3) - T_2) / (pow(T_3, 3) - T_3)) * line2.k;
			k += 6 * T_1 * T_2 * (avg_y_2 - avg_y_1) / ((pow(T_3, 3) - T_3));
			double b = avg_y_3 - k * (1 + T_3) / 2;
			Line l(k, b, p);
			re_pmf.add_new_line(l);
		}

	return re_pmf;
}
//标准维度skech叠加,s1+=s2;
void Sketh_on_Standard(Sketch& S1, Sketch& S2)
{
	double Ek = S1.E_k + S2.E_k;
	double Eb = S1.E_b + S2.E_b;
	double Vark = S1.Var_b + S2.Var_k;
	double Varb = S1.Var_b + S2.Var_b;
	double Ekb = (S1.E_k * S1.E_b) + (S1.E_k * S2.E_b) + (S2.E_k * S1.E_b) + (S2.E_k * S2.E_b);
	S1 = Sketch(Ek, Eb, Ekb, Vark, Varb);
}
//时间维度skech叠加,s1+=s2;
void Sketh_on_Time(Sketch& S1, Sketch& S2,double T_1, double T_2)
{
	//首先从T_1，T_3确定系数A1~4 B1~4
	//k_3 = A_1*k1+A2*b2+A3k2+A4b2;
	//b_3= B_1*k1+B2*b2+B3k2+B4b2;
	double T_3 = T_1 + T_2;
	double avg_x_1 = (1 + T_1) / 2;
	double avg_x_2 = (1 + T_2) / 2;
	double avg_x_3 = (1 + T_3) / 2;
	double A_1 = ((pow(T_1, 3) - T_1) / (pow(T_3, 3) - T_3)) - (6 * T_1 * T_2 * avg_x_1) / (pow(T_3, 3) - T_3);
	double A_2 = - (6 * T_1 * T_2) / (pow(T_3, 3) - T_3);
	double A_3 = ((pow(T_2, 3) - T_2) / (pow(T_3, 3) - T_3)) + (6 * T_1 * T_2 * avg_x_2) / (pow(T_3, 3) - T_3);
	double A_4 =(6 * T_1 * T_2) / (pow(T_3, 3) - T_3);
	double B_1 = T_1 * avg_x_1 / T_3 - (T_3 / 2) * A_1;
	double B_2 = T_1 / T_3 - (T_3 / 2) * A_2;
	double B_3 = T_2 * avg_x_2 / T_3 - (T_3 / 2) * A_3;
	double B_4 = T_2 / T_3 - (T_3 / 2) * A_4;
	double Ek = A_1 * S1.E_k + A_2 * S1.E_b + A_3 * S2.E_k + A_4 * S2.E_b;
	double Eb = B_1 * S1.E_k + B_2 * S1.E_b + B_3 * S2.E_k + B_4 * S2.E_b;
	double Cov_k_b_1 = S1.E_kb - S1.E_k * S1.E_b;
	double Cov_k_b_2 = S2.E_kb - S2.E_k * S2.E_b;
	double Vark = pow(A_1, 2) * S1.Var_k + pow(A_2, 2) * S1.Var_b + 2 * A_1 * A_2*Cov_k_b_1;
	Vark+= pow(A_3, 2) * S2.Var_k + pow(A_4, 2) * S2.Var_b + 2 * A_3 * A_4 * Cov_k_b_2;
	double Varb = pow(B_1, 2) * S1.Var_k + pow(B_2, 2) * S1.Var_b + 2 * B_1 * B_2 * Cov_k_b_1;
	Varb += pow(B_3, 2) * S2.Var_k + pow(B_4, 2) * S2.Var_b + 2 * B_3 * B_4 * Cov_k_b_2;
	//e(k_1^2)....
	double E_K1_2 = pow(S1.E_k, 2) + S1.Var_k;
	double E_B1_2 = pow(S1.E_b, 2) + S1.Var_b;
	double E_K2_2 = pow(S2.E_k, 2) + S2.Var_k;
	double E_B2_2 = pow(S2.E_b, 2) + S2.Var_b;
	//计算新的E(kb)
	double Ekb = A_1 * B_1 * E_K1_2 + (A_1 * B_2 + A_2 * B_1) * S1.E_kb + (A_1 * B_3 + A_3 * B_1) * S1.E_k * S2.E_k + (A_1 * B_4 + A_4 * B_1) * S1.E_k * S2.E_b;
	Ekb += A_2 * B_2 * E_B1_2 + (A_2 * B_3 + A_3 * B_2) * S1.E_b * S2.E_k + (A_2 * B_4 + A_4 * B_2) * S1.E_b * S2.E_b + A_3 * B_3 * E_K2_2 + (A_3 * B_4 + A_4 * B_3) * S2.E_kb + A_4 * B_4 * E_B2_2;
	S1 = Sketch(Ek, Eb, Ekb, Vark, Varb);
}

//计算一个FX的Sketch 
Sketch Get_sketch_of_Fx(Lines_pmf F_1)
{
	double Var_k, Var_b, E_kb;
	double ek = 0.0, eb = 0.0, ek2 = 0.0, eb2 = 0.0, ekb = 0.0;
	for (auto line : F_1.pmf)
	{
		ek += line.k * line.p;
		eb += line.b * line.p;
		ek2 += pow(line.k, 2) * line.p;
		eb2 += pow(line.b, 2) * line.p;
		ekb += line.k * line.b * line.p;
	}
	Var_k = ek2 - ek * ek;
	Var_b = eb2 - eb * eb;
	return Sketch(ek, eb, ekb, Var_k, Var_b);
}

//标准维度上柱状图压缩的卷积函数
void Convolution_at_Standard_Histogram(Lines_pmf F_1, Lines_pmf F_2,Lines_pmf &F_A,Lines_pmf &F_C) {
	

	for (auto line1 : F_1.pmf)
		for (auto line2 : F_2.pmf)
		{
			double k = line1.k + line2.k;
			double b = line1.b + line2.b;
			double p = line1.p * line2.p;
			Line l(k, b, p);
			F_A.add_new_line(l);
			F_C.add_new_line(l);
		}
}

//时间维度上的卷积函数
void Convolution_at_Time_Histogram(Lines_pmf F_1, Lines_pmf F_2, double T_1, double T_2,Lines_pmf& F_A, Lines_pmf& F_C) {
	
	Lines_pmf re_pmf;
	double T_3 = T_1 + T_2;
	double avg_x_1 = (1 + T_1) / 2;
	double avg_x_2 = (1 + T_2) / 2;
	double avg_x_3 = (1 + T_3) / 2;

	for (auto line1 : F_1.pmf)
		for (auto line2 : F_2.pmf)
		{
			double avg_y_1 = line1.k * avg_x_1 + line1.b;
			double avg_y_2 = line2.k * avg_x_2 + line2.b;
			double avg_y_3 = (T_1 * avg_y_1 + T_2 * avg_y_2) / T_3;
			double p = line1.p * line2.p;
			double k = ((pow(T_1, 3) - T_1) / (pow(T_3, 3) - T_3)) * line1.k + ((pow(T_2, 3) - T_2) / (pow(T_3, 3) - T_3)) * line2.k;
			k += 6 * T_1 * T_2 * (avg_y_2 - avg_y_1) / ((pow(T_3, 3) - T_3));
			double b = avg_y_3 - k * (1 + T_3) / 2;
			Line l(k, b, p);
			F_A.add_new_line(l);
			F_C.add_new_line3(l);
		}

}
void use_objects_to_get_pmfs_pws(double T, Object_pmfs O_i_pmfs, Lines_pmf& F_x)
{
	//当前pmf为空
	if (O_i_pmfs.pmfs.size() == 0)
		return;
	double A_T = double(pow(T, 3) - T) / (12.0);
	double avg_x = double(T + 1) / (2.0);
	if (F_x.pmf.size() == 0)
	{//当前pmf是第一个pmf，视为{(1,y),(2,0),(T,0),p}构造第一根线
		for (auto pmf : O_i_pmfs.pmfs)
		{
			double sum_x_i_y_i = pmf.x * pmf.y;
			double avg_x_y = avg_x * pmf.y;
			double k_i = (sum_x_i_y_i - avg_x_y) / A_T;
			double b_i = pmf.y - k_i * avg_x;
			double p_i = pmf.p;
			Line l(k_i, b_i, p_i);
			F_x.add_new_line_pws(l);
		}
		return;
	}
	else {
		//把pmf累加到已经聚合出来的线上
		int size = O_i_pmfs.pmfs.size();//当前pmf由几部分组成
		if (size <= 0)
			return;
		if (size > 1)
		{
			Lines_pmf re_pmf;
			//每一个instance都对现有的所有line进行修正
			for (int i = 0; i < O_i_pmfs.pmfs.size(); i++)
			{//累加公式为k=（xi*yi-avg(x)yi）/A(T)
				for (auto line : F_x.pmf)
				{
					double sum_x_i_y_i = O_i_pmfs.pmfs[i].x * O_i_pmfs.pmfs[i].y;
					double avgx_y = avg_x * O_i_pmfs.pmfs[i].y;
					double k_1 = line.k;
					double b_1 = line.b;
					double avg_y_1 = k_1 * avg_x + b_1;
					double p_2 = line.p * O_i_pmfs.pmfs[i].p;//更新概率
					double k_2 = k_1 + (sum_x_i_y_i - avgx_y) / A_T;//更新斜率
					double avg_y_2 = avg_y_1 + O_i_pmfs.pmfs[i].y / T;
					double b_2 = avg_y_2 - k_2 * avg_x;
					Line l(k_2, b_2, p_2);
					re_pmf.add_new_line_pws(l);
				}
			}
			F_x.pmf.clear();
			F_x.pmf.assign(re_pmf.pmf.begin(), re_pmf.pmf.end());

		}
		else if (size == 1)//当下只有一个pmf，其概率为1,则直接在fx里更新即可
		{
			for (auto line : F_x.pmf)
			{
				double sum_x_i_y_i = O_i_pmfs.pmfs[0].x * O_i_pmfs.pmfs[0].y;
				double avgx_y = avg_x * O_i_pmfs.pmfs[0].y;
				double k_1 = line.k;
				double b_1 = line.b;
				double avg_y_1 = k_1 * avg_x + b_1;
				line.p = line.p * O_i_pmfs.pmfs[0].p;//更新概率
				line.k = k_1 + (sum_x_i_y_i - avgx_y) / A_T;//更新斜率
				double avg_y_2 = avg_y_1 + O_i_pmfs.pmfs[0].y / T;
				line.b = avg_y_2 - line.k * avg_x;
			}
		}
		return;
	}



	return;


}
void use_objects_to_get_pmfs1(double T, Object_pmfs O_i_pmfs, Lines_pmf& F_x)
{
	//当前pmf为空
	if (O_i_pmfs.pmfs.size() == 0)
		return;
	double A_T = double(pow(T, 3) - T) / (12.0);
	double avg_x = double(T + 1) / (2.0);
	if (F_x.pmf.size() == 0)
	{//当前pmf是第一个pmf，视为{(1,y),(2,0),(T,0),p}构造第一根线
		for (auto pmf : O_i_pmfs.pmfs)
		{
			double sum_x_i_y_i = pmf.x * pmf.y;
			double avg_x_y = avg_x * pmf.y;
			double k_i = (sum_x_i_y_i - avg_x_y) / A_T;
			double b_i = pmf.y - k_i * avg_x;
			double p_i = pmf.p;
			Line l(k_i, b_i, p_i);
			F_x.add_new_line3(l);
		}
		return;
	}
	else {
		//把pmf累加到已经聚合出来的线上
		int size = O_i_pmfs.pmfs.size();//当前pmf由几部分组成
		if (size <= 0)
			return;
		if (size > 1)
		{
			Lines_pmf re_pmf;
			//每一个instance都对现有的所有line进行修正
			for (int i = 0; i < O_i_pmfs.pmfs.size(); i++)
			{//累加公式为k=（xi*yi-avg(x)yi）/A(T)
				for (auto line : F_x.pmf)
				{
					double sum_x_i_y_i = O_i_pmfs.pmfs[i].x * O_i_pmfs.pmfs[i].y;
					double avgx_y = avg_x * O_i_pmfs.pmfs[i].y;
					double k_1 = line.k;
					double b_1 = line.b;
					double avg_y_1 = k_1 * avg_x + b_1;
					double p_2 = line.p * O_i_pmfs.pmfs[i].p;//更新概率
					double k_2 = k_1 + (sum_x_i_y_i - avgx_y) / A_T;//更新斜率
					double avg_y_2 = avg_y_1 + O_i_pmfs.pmfs[i].y / T;
					double b_2 = avg_y_2 - k_2 * avg_x;
					Line l(k_2, b_2, p_2);
					re_pmf.add_new_line3(l);
				}
			}
			F_x.pmf.clear();
			F_x.pmf.assign(re_pmf.pmf.begin(), re_pmf.pmf.end());

		}
		else if (size == 1)//当下只有一个pmf，其概率为1,则直接在fx里更新即可
		{
			for (auto line : F_x.pmf)
			{
				double sum_x_i_y_i = O_i_pmfs.pmfs[0].x * O_i_pmfs.pmfs[0].y;
				double avgx_y = avg_x * O_i_pmfs.pmfs[0].y;
				double k_1 = line.k;
				double b_1 = line.b;
				double avg_y_1 = k_1 * avg_x + b_1;
				line.p = line.p * O_i_pmfs.pmfs[0].p;//更新概率
				line.k = k_1 + (sum_x_i_y_i - avgx_y) / A_T;//更新斜率
				double avg_y_2 = avg_y_1 + O_i_pmfs.pmfs[0].y / T;
				line.b = avg_y_2 - line.k * avg_x;
			}
		}
		return;
	}



	return;


}

void use_objects_to_get_pmfs_US(double T, Object_pmfs O_i_pmfs, Lines_pmf& F_x,bool islast=false)
{
	//当前pmf为空
	if (O_i_pmfs.pmfs.size() == 0)
		return;
	double A_T = double(pow(T, 3) - T) / (12.0);
	double avg_x = double(T + 1) / (2.0);
	if (F_x.pmf.size() == 0)
	{//当前pmf是第一个pmf，视为{(1,y),(2,0),(T,0),p}构造第一根线
		for (auto pmf : O_i_pmfs.pmfs)
		{
			double sum_x_i_y_i = pmf.x * pmf.y;
			double avg_x_y = avg_x * pmf.y;
			double  k_i = (sum_x_i_y_i - avg_x_y);
			double b_i = pmf.y;
			if (islast == true)
			{
				k_i = k_i / A_T;
				b_i = b_i/A_T - k_i * avg_x;
			}
			
			double p_i = pmf.p;
			Line l(k_i, b_i, p_i);
			F_x.add_new_line3(l);
		}
		return;
	}
	else {
		//把pmf累加到已经聚合出来的线上
		int size = O_i_pmfs.pmfs.size();//当前pmf由几部分组成
		if (size <= 0)
			return;
		if (size > 1)
		{
			Lines_pmf re_pmf;
			//每一个instance都对现有的所有line进行修正
			for (int i = 0; i < O_i_pmfs.pmfs.size(); i++)
			{//累加公式为k=（xi*yi-avg(x)yi）/A(T)
				for (auto line : F_x.pmf)
				{
					double sum_x_i_y_i = O_i_pmfs.pmfs[i].x * O_i_pmfs.pmfs[i].y;
					double avgx_y = avg_x * O_i_pmfs.pmfs[i].y;
					double k_1 = line.k;
					double b_1 = line.b;
					
					double avg_y_1 = k_1 * avg_x + b_1;
					double p_2 = line.p * O_i_pmfs.pmfs[i].p;//更新概率
					double k_2 = k_1 + (sum_x_i_y_i - avgx_y);//更新斜率
					//double avg_y_2 = avg_y_1 + O_i_pmfs.pmfs[i].y / T;
					double b_2 = b_1+ O_i_pmfs.pmfs[i].y;
					if (islast == true)
					{
						k_2 = k_2 / A_T;
						b_2 = b_2 / A_T - k_2 * avg_x;
					}
					
					
					Line l(k_2, b_2, p_2);
					re_pmf.add_new_line3(l);
				}
			}
			F_x.pmf.clear();
			F_x.pmf.assign(re_pmf.pmf.begin(), re_pmf.pmf.end());

		}
		else if (size == 1)//当下只有一个pmf，其概率为1,则直接在fx里更新即可
		{
			for (auto line : F_x.pmf)
			{
				double sum_x_i_y_i = O_i_pmfs.pmfs[0].x * O_i_pmfs.pmfs[0].y;
				double avgx_y = avg_x * O_i_pmfs.pmfs[0].y;
				double k_1 = line.k;
				double b_1 = line.b;
				double avg_y_1 = k_1 * avg_x + b_1;
				line.p = line.p * O_i_pmfs.pmfs[0].p;//更新概率
				line.k = k_1 + (sum_x_i_y_i - avgx_y);//更新斜率	
				line.b = line.b + O_i_pmfs.pmfs[0].y;
				if (islast == true)
				{
					line.k = line.k / A_T;
					line.b = line.b / A_T - line.b * avg_x;
				}

			}
		}
		return;
	}



	return;


}
/*
//叶结点中计算fc的函数
Lines_pmf use_objects_to_get_pmfs_lattice(double T, Object_pmfs O_i_pmfs, Lines_pmf& F_x,bool isfirst)
{
	Lines_pmf re_pmf;
	//当前pmf为空
	if (O_i_pmfs.pmfs.size() == 0)
		return F_x;
	double A_T = double(pow(T, 3) - T) / (12.0);
	double avg_x = double(T + 1) / (2.0);
	if (F_x.k_search.size() == 0)
	{//当前pmf是第一个pmf，视为{(1,y),(2,0),(T,0),p}构造第一根线
		for (auto pmf : O_i_pmfs.pmfs)
		{
			double sum_x_i_y_i = pmf.x * pmf.y;
			double avg_x_y = avg_x * pmf.y;
			double k_i = (sum_x_i_y_i - avg_x_y) / A_T;
			double b_i = pmf.y - k_i * avg_x;
			double p_i = pmf.p;
			Line l(k_i, b_i, p_i);
			re_pmf.add_line_lattice(l);
		}
		return re_pmf;
	}
	else {
		//把pmf累加到已经聚合出来的线上
		int size = O_i_pmfs.pmfs.size();//当前pmf由几部分组成
		/*if (size <= 0)
			return re_pmf;
			*/
			//每一个instance都对现有的所有line进行修正
			/*if (isfirst == true) {
				for (int i = 0; i < O_i_pmfs.pmfs.size(); i++)
				{//累加公式为k=（xi*yi-avg(x)yi）/A(T)
					for (auto line : F_x.pmf)
					{
						double sum_x_i_y_i = O_i_pmfs.pmfs[i].x * O_i_pmfs.pmfs[i].y;
						double avgx_y = avg_x * O_i_pmfs.pmfs[i].y;
						double k_1 = line.k;
						double b_1 = line.b;
						double avg_y_1 = k_1 * avg_x + b_1;
						double p_2 = line.p * O_i_pmfs.pmfs[i].p;//更新概率
						double k_2 = k_1 + (sum_x_i_y_i - avgx_y) / A_T;//更新斜率
						double avg_y_2 = avg_y_1 + O_i_pmfs.pmfs[i].y / T;
						double b_2 = avg_y_2 - k_2 * avg_x;
						Line l(k_2, b_2, p_2);
						re_pmf.add_line_lattice(l);
					}
				}
				F_x.pmf.clear();
				F_x.pmf.assign(re_pmf.pmf.begin(), re_pmf.pmf.end());
				F_x.Now.clear();
				F_x.Now.assign(re_pmf.Now.begin(), re_pmf.Now.end());
			}
			//else
			{
				for (int i = 0; i < O_i_pmfs.pmfs.size(); i++)
				{//累加公式为k=（xi*yi-avg(x)yi）/A(T)
					for(int j=0;j< F_x.k_search.size();j++)
					{
						//Line line = F_x.pmf[k];
						double k_1 = tan(double(F_x.k_search[j]) / 10 - 1.58);
						double b_1 = F_x.b_search[j] * F_x.N;
						double p_1 = F_x.p[F_x.b_search[j]][F_x.k_search[j]];
						F_x.p[F_x.b_search[j]][F_x.k_search[j]] = 0;
						double sum_x_i_y_i = O_i_pmfs.pmfs[i].x * O_i_pmfs.pmfs[i].y;
						double avgx_y = avg_x * O_i_pmfs.pmfs[i].y;
						//double k_1 = line.k;
						//double b_1 = line.b;
						double avg_y_1 = k_1 * avg_x + b_1;
						double p_2 =p_1 * O_i_pmfs.pmfs[i].p;//更新概率
						double k_2 = k_1 + (sum_x_i_y_i - avgx_y) / A_T;//更新斜率
						double avg_y_2 = avg_y_1 + O_i_pmfs.pmfs[i].y / T;
						double b_2 = avg_y_2 - k_2 * avg_x;
						Line l(k_2, b_2, p_2);
						re_pmf.add_line_lattice(l);
					}
				}
				//F_x.pmf.clear();
				//F_x.pmf.assign(re_pmf.pmf.begin(), re_pmf.pmf.end());
				//F_x.Now.clear();
				//F_x.Now.assign(re_pmf.Now.begin(), re_pmf.Now.end());
				
			}
			

		//}
		else if (size == 1)//当下只有一个pmf，其概率为1,则直接在fx里更新即可
		{
			for (auto line : F_x.pmf)
			{
				double sum_x_i_y_i = O_i_pmfs.pmfs[0].x * O_i_pmfs.pmfs[0].y;
				double avgx_y = avg_x * O_i_pmfs.pmfs[0].y;
				double k_1 = line.k;
				double b_1 = line.b;
				double avg_y_1 = k_1 * avg_x + b_1;
				line.p = line.p * O_i_pmfs.pmfs[0].p;//更新概率
				line.k = k_1 + (sum_x_i_y_i - avgx_y) / A_T;//更新斜率
				double avg_y_2 = avg_y_1 + O_i_pmfs.pmfs[0].y / T;
				line.b = avg_y_2 - line.k * avg_x;
			}
		}
		return re_pmf;
	}



	return re_pmf;


}


Lines_pmf Convolution_at_Standard_lattice(Lines_pmf F_1, Lines_pmf F_2) {
	if (F_1.pmf.size() == 0)
		return F_2;
	if (F_2.pmf.size() == 0)
		return F_1;

	Lines_pmf re_pmf;
	re_pmf.N = F_1.N + F_2.N;
	//Line A;
	//re_pmf.pmf = vector<Line>(re_pmf.N * re_pmf.N * 130, A);
	for(int i=0;i<F_1.k_search.size();i++)
	{
		for (int j=0; j < F_2.k_search.size(); j++)
		{
			//double k_1
			double k_1 = tan(double(F_1.k_search[i]) / 10 - 1.58);
			double b_1 = F_1.b_search[i] * F_1.N;
			double p_1 = F_1.p[F_1.b_search[i]][F_1.k_search[i]];
			double k_2 = tan(double(F_2.k_search[j]) / 10 - 1.58);
			double b_2 = F_2.b_search[j] * F_2.N;
			double p_2 = F_2.p[F_2.b_search[j]][F_2.k_search[j]];
			
			double k = k_1 + k_2;
			double b = b_1 + b_2;
			double p = p_1 * p_2;
			Line l(k, b, p);
			re_pmf.add_line_lattice(l);
		}
	}
	return re_pmf;
}

Lines_pmf Convolution_at_Time_lattice(Lines_pmf F_1, Lines_pmf F_2, double T_1, double T_2) {
	if (F_1.k_search.size() == 0)
		return F_2;
	if (F_2.k_search.size() == 0)
		return F_1;
	Lines_pmf re_pmf;
	re_pmf.isTime = true;
	double T_3 = T_1 + T_2;
	double avg_x_1 = (1 + T_1) / 2;
	double avg_x_2 = (1 + T_2) / 2;
	double avg_x_3 = (1 + T_3) / 2;

	for (auto i : F_1.Now)
	{
		for (auto j : F_2.Now)
		{
			Line line1 = F_1.pmf[i];
			Line line2 = F_2.pmf[j];
			double avg_y_1 = line1.k * avg_x_1 + line1.b;
			double avg_y_2 = line2.k * avg_x_2 + line2.b;
			double avg_y_3 = (T_1 * avg_y_1 + T_2 * avg_y_2) / T_3;
			double p = line1.p * line2.p;
			double k = ((pow(T_1, 3) - T_1) / (pow(T_3, 3) - T_3)) * line1.k + ((pow(T_2, 3) - T_2) / (pow(T_3, 3) - T_3)) * line2.k;
			k += 6 * T_1 * T_2 * (avg_y_2 - avg_y_1) / ((pow(T_3, 3) - T_3));
			double b = avg_y_3 - k * (1 + T_3) / 2;
			Line l(k, b, p);
			re_pmf.add_new_line(l);
		}
	}
	return re_pmf;
}*/