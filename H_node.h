#pragma once
#include"data_instance.h"
#include"data_object_pmf.h"
#include"Lines.h"
#include"until.h"
static int T = 30;
extern double Space_Cost_all;
class H_node
{
public:
	H_node();
	~H_node();
	double Materalize_time = 0.0;
	bool IsLeaf = false;//�ж��Ƿ�ΪҶ�ӽ��
	bool IsRoot = false;//�ж��Ƿ�Ϊ�����
	bool IsTime = false;//�жϵ�ǰά���Ƿ���ʱ��ά��
	bool isMaterialize = false;//�жϵ�ǰcuboid�Ƿ��Ѿ����ﻯ
	double Space_cost = 0.0;//�ߵ���������λ�Ǹ�
	double Space_cost_sketch = 0.0;//skech�ĸ�����
	double Space_data = 0.0;//��λ��B���ֽ�
	double Space_pws = 0.0;//�ϲ�ʹ��pws����Ŀռ俪������λ��B
	vector<string> Dimension;//ά�����б�
	string Now_Dimenson;//��ǰά��
	int T=-1;//����һ��T
	int Now_window = 0;
	vector<Object_pmfs> Quarantine_area;//������
	vector<Object_pmfs> New_complete;//�µ���������
	vector<Object_pmfs> Complete_obj;//Ҷ�ӽ���е��������������ǻ��fc
	vector<Object_pmfs> All_data_instance;//����pmfs
	vector<H_node> Children_node_list;//�ӽ��
	Lines_pmf F_C, F_A, F_T, F_N;//fc,fa
	Sketch S_C, S_A;//sc,sa
	void add_pmf(instance_pmf pmf);
	void Use_object_pmf_to_get_line_pmf_on_leaf_node();//Ҷ��㣨��ϸ�����ϻ�ûع���pmf�ĺ���
	void Use_object_pmf_to_get_line_pmf_on_leaf_node_pws();//Ҷ��㣨��ϸ�����ϻ�ûع���pmf�ĺ���
	void Materialize_Not_leaf();
	void Materialize_node();
	void Materialize_node_by_sketch();
	void Materialize_node_Histogram();
	void Materialize_Not_leaf_Histogram();
	void Materialize_node_by_pws();
	void Materialize_Not_leaf_by_sketch();
	void Materialize_Not_leaf_by_pws();
	void get_sketch_on_leaf_node();
	void Use_Fc_Quarantine_area_to_get_Fa_in_node();//�ڵ���ʹ��FC�͸���������FA�ĺ���
	void Use_Fc_New_complete_to_get_FC_in_node();//�ڵ���ʹ��FC�������ɵ������������FC�ĺ���
	void Use_SC_Quarantine_area_to_get_SA_in_node();//�ڵ���ʹ��FC�͸���������FA�ĺ���
	void Use_Sc_New_complete_to_get_SC_in_node();//�ڵ���ʹ��FC�������ɵ������������FC�ĺ���
	void query_pmf_DFS();
	void query_sketch_DFS();
	void query_sketch_BFS();
	void query_pmf_BFS();
	void Materialize_node_with_path();
	void agg_cell();
	void on_line_sketch(double eM, double varM, double tc);
	void Materialize_node_for_time();
	void on_line_sketch(Data_object oi, double tc);
private:

};

H_node::H_node()
{
}

H_node::~H_node()
{
}
void H_node::add_pmf(instance_pmf pmf) {
	for (auto &object : this->All_data_instance)
	{
		if (object.id == pmf.id)
		{
			object.pmfs.push_back(pmf);
			return;
		}
	}
	Object_pmfs o_i(pmf);
	this->All_data_instance.push_back(o_i);
}


void H_node::Use_object_pmf_to_get_line_pmf_on_leaf_node()
{
	if (this->IsLeaf == false)
	{
		cout << "Not leaf to get lines' pmf,error!" << endl;
		return;
	}
	this->T = 30;
	cout << "��ǰ�ﻯcuboidΪ:";
	for (auto s : this->Dimension)
	{
		cout << s << "  ";
	}
	/*
	//������ǰϸ���е�����pmf,ʹ�����е�����pmf�����㼯fC���ضϵķ��������
	int sum_pmf_size = 0;
	
	for (int i=0;i< this->All_data_instance.size();i++)
	{
		Object_pmfs pmf_i = this->All_data_instance[i];
		sum_pmf_size += pmf_i.pmfs.size();
		//����pmf
		if (pmf_i.isTruncated()==false)
		{
			bool islate = false;
			if (i == this->All_data_instance.size() - 1)
			{
				islate = true;
			}
			//use_objects_to_get_pmfs_US(T,pmf_i,F_C,islate);
			use_objects_to_get_pmfs(T, pmf_i, F_C);
		}
		else//���������
		{
			Quarantine_area.push_back(pmf_i);
		}
	}
	cout << "F_C�Ĵ�СΪ��" << this->F_C.pmf.size() << endl;
	this->Use_Fc_Quarantine_area_to_get_Fa_in_node();
	this->Space_cost = double(int(this->F_A.pmf.size()) + (this->F_C.pmf.size()));
	this->Space_data = double(sum_pmf_size * (4 + 3 * 8));
	isMaterialize = true;
	return;
	*/
	//������ǰϸ���е�����pmf,ʹ�����е�����pmf�����㼯fC���ضϵķ��������
	int sum_pmf_size = 0;

	for (int i = 0; i < this->All_data_instance.size(); i++)
	{
		Object_pmfs pmf_i = this->All_data_instance[i];
		//sum_pmf_size += pmf_i.pmfs.size();
		//����pmf
		if (pmf_i.isTruncated() == false){
			use_objects_to_get_pmfs(T, pmf_i, F_C);
		}
		{
			Quarantine_area.push_back(pmf_i);
		}
		/*int number = 4;
		if (i < number)
		{
			use_objects_to_get_pmfs(T, pmf_i, F_C);
		}
		else
		{
			bool isfirst = false;
			if (i == number)
			{
				isfirst = true;
			}*/
			//use_objects_to_get_pmfs_lattice(T, pmf_i, F_C, isfirst);
		//}
		{
			/*bool islate = false;
			if (i == this->All_data_instance.size() - 1)
			{
				islate = true;
			}*/
			//use_objects_to_get_pmfs_US(T,pmf_i,F_C,islate);
			//�����cell�� k��Χ��(-0.6,0.5),b��Χ��(0,11),
		//	use_objects_to_get_pmfs_lattice(T, pmf_i, F_C);
		}
		/*else//���������
		{
			Quarantine_area.push_back(pmf_i);
		}*/
	}
	cout << "F_C�Ĵ�СΪ��" << this->F_C.pmf.size() << endl;
	this->Use_Fc_Quarantine_area_to_get_Fa_in_node();
	this->Space_cost = double(int(this->F_A.pmf.size()) + (this->F_C.pmf.size()));
	this->Space_data = double(sum_pmf_size * (4 + 3 * 8));
	isMaterialize = true;
	return;

}
void H_node::Use_Fc_Quarantine_area_to_get_Fa_in_node()
{
	if (this->Quarantine_area.size() == 0)
	{
		this->F_A = this->F_C;
		return;
	}
	//double A_T = double(pow(T, 3) - T) / (12.0);
	//double avg_x = double(T + 1) / (2.0);

	Lines_pmf F_T;
	//����������
	for (auto pmf : this->Quarantine_area)
	{
		if (pmf.isTruncated_pmf == false)
			continue;
		pmf.Make_Truncated_complete();
		use_objects_to_get_pmfs(T,pmf,F_T);
	}
	//
	this->F_T = F_T;
	this->F_A = Convolution_at_Standard(this->F_C, this->F_T);
}
void  H_node::Use_Fc_New_complete_to_get_FC_in_node()
{
	if (this->New_complete.size() == 0)
		return;
	

	Lines_pmf F_N;
	//����������
	for (auto pmf : this->Quarantine_area)
	{
		if (pmf.isTruncated_pmf == true)
		{
			cout << "wrong!!at 127 H_node.h" << endl;
			continue;
		}
			
		pmf.Make_Truncated_complete();
		use_objects_to_get_pmfs(T, pmf, F_N);
	}
	//
	//Lines_pmf fc= Convolution_at_Standard(this->F_C, F_N);
	//this->F_C.pmf.clear();
	//this->F_C.pmf.assign(fc.pmf.begin(), fc.pmf.end());
	this->F_C = Convolution_at_Standard(this->F_C, F_N);
}
void H_node::Use_SC_Quarantine_area_to_get_SA_in_node()
{
	if (this->Quarantine_area.size() <= 0)
	{
		this->F_A = this->F_C;
		return;
	}
	Lines_pmf F_T;
	//����������
	for (auto pmf : this->Quarantine_area)
	{
		if (pmf.isTruncated_pmf == false)
			continue;
		pmf.Make_Truncated_complete();
		use_objects_to_get_pmfs(T, pmf, F_T);
	}
	//
	this->F_T = F_T;
	this->S_A = this->S_C;
	Sketch S_T = Get_sketch_of_Fx(this->F_T);
	Sketh_on_Standard(this->S_A, S_T);
}
void H_node::Use_Sc_New_complete_to_get_SC_in_node()
{
	if (this->New_complete.size() <= 0)
		return;


	Lines_pmf F_N;
	//����������
	for (auto pmf : this->Quarantine_area)
	{
		if (pmf.isTruncated_pmf == true)
		{
			cout << "wrong!!at 127 H_node.h" << endl;
			continue;
		}

		pmf.Make_Truncated_complete();
		use_objects_to_get_pmfs(T, pmf, F_N);
	}
	Sketch S_N = Get_sketch_of_Fx(F_N);
	Sketh_on_Standard(this->S_C, S_N);


}

void H_node::Materialize_Not_leaf()
{

	if (this->IsLeaf == true)
	{
		cout << "Leaf node,dont use this funcation" << endl;
		return;
	}
	cout << "��ǰ�ﻯcuboidΪ:";
	for (auto s : this->Dimension)
	{
		cout << s << "  ";
	}
	if (this->Children_node_list.size() == 0)
	{
		cout << "*****" << endl;
		return;
	}


		//��ǰ����һ���ӽ��
		if (this->Children_node_list.size() == 1)
		{
			this->T = this->Children_node_list[0].T;
			this->F_C.pmf.assign(this->Children_node_list[0].F_C.pmf.begin(), this->Children_node_list[0].F_C.pmf.end());
			this->F_A.pmf.assign(this->Children_node_list[0].F_A.pmf.begin(), this->Children_node_list[0].F_A.pmf.end());
			this->Quarantine_area.assign(this->Children_node_list[0].Quarantine_area.begin(), this->Children_node_list[0].Quarantine_area.end());
			this->F_C.N = this->Children_node_list[0].F_C.N;
			this->F_C.Now.assign(this->Children_node_list[0].F_C.Now.begin(), this->Children_node_list[0].F_C.Now.end());
			cout << "F_C�Ĵ�СΪ��" << this->F_C.pmf.size() << endl;
			this->Space_cost += this->Children_node_list[0].Space_cost;
			this->Space_cost += double(int(this->F_C.pmf.size()) + this->F_A.pmf.size());
			this->Space_data = this->Children_node_list[0].Space_data;
			return;
		}
		
		//�����ȸ��Ƹ�����,�ȸ��Ƶ�һ���ӽڵ�ĸ�����
		this->Quarantine_area.assign(this->Children_node_list[0].Quarantine_area.begin(), this->Children_node_list[0].Quarantine_area.end());
		//�ӵڶ�����ʼ�������û���ܺϲ���
		for (int i=1;i<this->Children_node_list.size();i++)
		{
			
			//�����������ж���pmf
			int size_now = this->Quarantine_area.size();
			if (size_now <= 0)
				continue;
			//��������ÿһ���ӽ��ĸ�����
			for (auto pmf : this->Children_node_list[i].Quarantine_area)
			{
				//�ܷ�ϲ��ı�־
				bool flag = false;
				for (int j = 0; j < size_now; j++)
				{
					//���кϲ�
					if (pmf.id == this->Quarantine_area[j].id)
					{
						this->Quarantine_area[j].add_new_pmfs(pmf);
						flag = true;
						//�����������,��������ֱ����
						if (this->Quarantine_area[j].isTruncated() == false)
						{
							this->Quarantine_area[j].isTruncated_pmf = false;
							this->New_complete.push_back(this->Quarantine_area[j]);
							
						}
						break;
					}

				}
				if (flag == false)//���뵱ǰ������
				{
					this->Quarantine_area.push_back(pmf);
				}
			}
		}
		//���ȸ��Ƶ�һ���ӽ���fc
		Lines_pmf F_C_now;
		F_C_now.pmf.assign(this->Children_node_list[0].F_C.pmf.begin(), this->Children_node_list[0].F_C.pmf.end());
		F_C_now.N = this->Children_node_list[0].F_C.N;
		F_C_now.Now.assign(this->Children_node_list[0].F_C.Now.begin(), this->Children_node_list[0].F_C.Now.end());
		this->Space_cost += this->Children_node_list[0].Space_cost;
		this->Space_data = this->Children_node_list[0].Space_data;
		if (this->IsTime==true)//ʱ��ά���ϵľۺ�
		{
			this->T = this->Children_node_list[0].T;
			for (int i = 1; i < this->Children_node_list.size(); i++)
			{
				
				this->Space_cost += this->Children_node_list[i].Space_cost;
				this->Space_data += this->Children_node_list[i].Space_data;
				F_C_now = Convolution_at_Time(F_C_now, this->Children_node_list[i].F_C,this->T, this->Children_node_list[i].T);
				this->T += this->Children_node_list[i].T;
			}
			this->F_C.pmf.clear();
			this->F_C.pmf.swap(F_C_now.pmf);

		}
		else//��׼ά���ϵľۺ�
		{
			this->T = this->Children_node_list[0].T;
			for (int i = 1; i < this->Children_node_list.size(); i++)
			{
				this->Space_cost += this->Children_node_list[i].Space_cost;
				this->Space_data += this->Children_node_list[i].Space_data;
				F_C_now =Convolution_at_Standard(F_C_now, this->Children_node_list[i].F_C);
				
			}
			this->F_C.pmf.clear();
			this->F_C.pmf.swap(F_C_now.pmf);
			//this->F_C.Now.clear();
			//this->F_C.Now.swap(F_C_now.Now);
		}
		
		//this->Use_Fc_New_complete_to_get_FC_in_node();
		//this->Use_Fc_Quarantine_area_to_get_Fa_in_node();
		this->isMaterialize = true;
		this->Space_cost += double(int(this->F_C.pmf.size()) + int(this->F_A.pmf.size()));
		cout << "F_C�Ĵ�СΪ��" << this->F_C.pmf.size() << endl;
		return;







}
void H_node::Materialize_node()
{
	if (this->IsLeaf == true)
	{
		double start = clock();
		this->Use_object_pmf_to_get_line_pmf_on_leaf_node();
		this->isMaterialize = true;
		double end = clock();
		double t=(end-start)/ CLOCKS_PER_SEC;
		this->Materalize_time = t;
		return;
	}

	for (auto &child_node : this->Children_node_list)
	{
		if (child_node.isMaterialize == false)//�ӽ��û�б��ﻯ
		{
			child_node.Materialize_node();
		}
	}
	//�����ӽ�㶼���ﻯ��
	double start = clock();
	this->Materialize_Not_leaf();
	double end = clock();
	double t = (end - start) / CLOCKS_PER_SEC;
	this->Materalize_time = t;
	return;
}
void H_node::Materialize_node_Histogram()
{
	if (this->IsLeaf == true)
	{
		this->Use_object_pmf_to_get_line_pmf_on_leaf_node();
		this->isMaterialize = true;
		return;
	}

	for (auto& child_node : this->Children_node_list)
	{
		if (child_node.isMaterialize == false)//�ӽ��û�б��ﻯ
		{
			child_node.Materialize_node_Histogram();
		}
	}
	//�����ӽ�㶼���ﻯ��
	this->Materialize_Not_leaf_Histogram();
	return;
}
void H_node::Materialize_Not_leaf_Histogram() {
	string D = "\"\"\"1992\"\"\"";
	if (this->IsTime == true && this->Now_Dimenson == "\"\"\"1992\"\"\"")
	{
		if (this->IsLeaf == true)
		{
			std::cout << "Leaf node,dont use this funcation" << endl;
			return;
		}
		std::cout << "��ǰ�ﻯcuboidΪ:";
		for (auto s : this->Dimension)
		{
			cout << s << "  ";
		}
		if (this->Children_node_list.size() == 0)
		{
			std::cout << "*****" << endl;
			return;
		}
		//��ǰ����һ���ӽ��
		if (this->Children_node_list.size() == 1)
		{
			this->T = this->Children_node_list[0].T;
			this->F_C.pmf.assign(this->Children_node_list[0].F_C.pmf.begin(), this->Children_node_list[0].F_C.pmf.end());
			this->F_A.pmf.assign(this->Children_node_list[0].F_A.pmf.begin(), this->Children_node_list[0].F_A.pmf.end());
			this->Quarantine_area.assign(this->Children_node_list[0].Quarantine_area.begin(), this->Children_node_list[0].Quarantine_area.end());
			std::cout << "F_C�Ĵ�СΪ��" << this->F_C.pmf.size() << endl;
			this->Space_cost += this->Children_node_list[0].Space_cost;
			this->Space_cost += double(int(this->F_C.pmf.size()) + this->F_A.pmf.size());
			this->Space_data = this->Children_node_list[0].Space_data;
			return;
		}
		//�����ȸ��Ƹ�����,�ȸ��Ƶ�һ���ӽڵ�ĸ�����
		this->Quarantine_area.assign(this->Children_node_list[0].Quarantine_area.begin(), this->Children_node_list[0].Quarantine_area.end());
		//�ӵڶ�����ʼ�������û���ܺϲ���
		for (int i = 1; i < this->Children_node_list.size(); i++)
		{

			//�����������ж���pmf
			int size_now = this->Quarantine_area.size();
			if (size_now <= 0)
				continue;
			//��������ÿһ���ӽ��ĸ�����
			for (auto pmf : this->Children_node_list[i].Quarantine_area)
			{
				//�ܷ�ϲ��ı�־
				bool flag = false;
				for (int j = 0; j < size_now; j++)
				{
					//���кϲ�
					if (pmf.id == this->Quarantine_area[j].id)
					{
						this->Quarantine_area[j].add_new_pmfs(pmf);
						flag = true;
						//�����������,��������ֱ����
						if (this->Quarantine_area[j].isTruncated() == false)
						{
							this->Quarantine_area[j].isTruncated_pmf = false;
							this->New_complete.push_back(this->Quarantine_area[j]);

						}
						break;
					}

				}
				if (flag == false)//���뵱ǰ������
				{
					this->Quarantine_area.push_back(pmf);
				}
			}
		}
		//���ȸ��Ƶ�һ���ӽ���fc
		Lines_pmf F_C_now;
		F_C_now.pmf.assign(this->Children_node_list[0].F_C.pmf.begin(), this->Children_node_list[0].F_C.pmf.end());
		this->Space_cost += this->Children_node_list[0].Space_cost;
		this->Space_data = this->Children_node_list[0].Space_data;
		if (this->Children_node_list[0].IsTime == true)//ʱ��ά���ϵľۺ�
		{
			this->T = this->Children_node_list[0].T;
			for (int i = 1; i < this->Children_node_list.size(); i++)
			{

				this->Space_cost += this->Children_node_list[i].Space_cost;
				this->Space_data += this->Children_node_list[i].Space_data;
				Convolution_at_Time_Histogram(F_C_now, this->Children_node_list[i].F_C, this->T, this->Children_node_list[i].T, this->F_A, this->F_C);
				this->T += this->Children_node_list[i].T;
			}
		}
			
		else//��׼ά���ϵľۺ�
		{
			this->T = this->Children_node_list[0].T;
			for (int i = 1; i < this->Children_node_list.size(); i++)
			{
				this->Space_cost += this->Children_node_list[i].Space_cost;
				this->Space_data += this->Children_node_list[i].Space_data;
				 Convolution_at_Standard_Histogram(F_C_now, this->Children_node_list[i].F_C, this->F_A, this->F_C);

			}
		}
		this->isMaterialize = true;
		this->Space_cost += double(int(this->F_C.pmf.size()) + int(this->F_A.pmf.size()));
		std::cout << "F_C�Ĵ�СΪ��" << this->F_C.pmf.size() << endl;
		return;
	}
	else
	{
	this->Materialize_Not_leaf();
	}














}
void H_node::get_sketch_on_leaf_node()
{
	double A_T = this->T;
	double Var_k=0.0, Var_b=0.0, Cov_kb =0.0;
	double ek = 0.0, eb = 0.0, ek2 = 0.0, eb2 = 0.0, ekb = 0.0;
	//this->Use_object_pmf_to_get_line_pmf_on_leaf_node();
	for (auto pmf_i : this->All_data_instance)
	{
		//����pmf
		//if (pmf_i.isTruncated() == false)
		{
			double eM=0,eM_2=0, varM=0;
			double sum = 0.0;
			for (auto pmf : pmf_i.pmfs)
			{
				eM += pmf.y * pmf.p;
				eM_2 += pmf.y * pmf.y * pmf.p;
				sum += pmf.p;
				
			}
			if (sum > 0.96)
			{
				varM = eM_2 - eM * eM;
				double tc = pmf_i.pmfs[0].x;
				ek += (double(6 * (2 * tc - A_T - 1)) / double(A_T * A_T * A_T - A_T)) * eM;
				Var_k += (double(6 * (2 * tc - A_T - 1)) / double(A_T * A_T * A_T - A_T)) * (double(6 * (2 * tc - A_T - 1)) / double(A_T * A_T * A_T - A_T)) * varM;
				eb += (double(2.0 * (2 * A_T - 3 * tc + 1)) / double(A_T * (A_T - 1))) * eM;
				Var_b += (double(2.0 * (2 * A_T - 3 * tc + 1)) / double(A_T * (A_T - 1))) * (double(2.0 * (2 * A_T - 3 * tc + 1)) / double(A_T * (A_T - 1))) * varM;
				Cov_kb += (double(6 * (2 * tc - A_T - 1)) / double(A_T * A_T * A_T - A_T)) * (double(2.0 * (2 * A_T - 3 * tc + 1)) / double(A_T * (A_T - 1))) * varM;
				ekb = Cov_kb + ek * eb;
				this->S_C = Sketch(ek, eb, ekb, Var_k, Var_b);
				this->S_A = Sketch(ek, eb, ekb, Var_k, Var_b);
				Sketch S_t = Get_sketch_of_Fx(this->F_T);
				Sketh_on_Standard(this->S_A, S_t);
				//this->Space_cost_sketch = 5;
				//cout << S_C.E_k << "  " << S_C.E_b << endl;
				this->Space_cost_sketch = 5;

			}
			else//���������
			{
				Quarantine_area.push_back(pmf_i);
			}
			}
	}
	
}
void H_node::Materialize_node_by_sketch()
{
	if (this->IsLeaf == true)
	{
		this->get_sketch_on_leaf_node();
		this->isMaterialize = true;
		return;
	}

	for (auto& child_node : this->Children_node_list)
	{
		if (child_node.isMaterialize == false)//�ӽ��û�б��ﻯ
		{
			child_node.Materialize_node_by_sketch();
		}
	}
	//�����ӽ�㶼���ﻯ��
	if (this->IsRoot == true)
	{
		this->Space_cost_sketch += 5.0;
		return;
	}
	this->Materialize_Not_leaf_by_sketch();
	return;


}
void H_node::Materialize_Not_leaf_by_sketch()
{

	/*if (this->IsLeaf == true)
	{
		std::cout << "Leaf node,dont use this funcation" << endl;
		return;
	}
	std::cout << "��ǰ�ﻯcuboidΪ:";
	for (auto s : this->Dimension)
	{
		std::cout << s << "  ";
	}
	*/
	//��ǰ����һ���ӽ��
	if (this->Children_node_list.size() == 1)
	{
		this->T = this->Children_node_list[0].T;
		this->S_C = this->Children_node_list[0].S_C;
		this->S_A = this->Children_node_list[0].S_A;
		this->Quarantine_area.assign(this->Children_node_list[0].Quarantine_area.begin(), this->Children_node_list[0].Quarantine_area.end());
		cout << "F_C�Ĵ�СΪ��" << this->F_C.pmf.size() << endl;
		cout << "EK and Eb:   " << this->S_C.E_k << "  " << this->S_C.E_b << endl;
		this->Space_cost += this->Children_node_list[0].Space_cost;
		this->Space_data = this->Children_node_list[0].Space_data;
		this->Space_cost_sketch = 5;
		this->Space_cost_sketch += this->Children_node_list[0].Space_cost_sketch;
		return;
	}

	//�����ȸ��Ƹ�����,�ȸ��Ƶ�һ���ӽڵ�ĸ�����
	this->Quarantine_area.assign(this->Children_node_list[0].Quarantine_area.begin(), this->Children_node_list[0].Quarantine_area.end());
	//�ӵڶ�����ʼ�������û���ܺϲ���
	for (int i = 1; i < this->Children_node_list.size(); i++)
	{
		//�����������ж���pmf
		int size_now = this->Quarantine_area.size();
		//��������ÿһ���ӽ��ĸ�����
		for (auto pmf : this->Children_node_list[i].Quarantine_area)
		{
			//�ܷ�ϲ��ı�־
			bool flag = false;
			for (int j = 0; j < size_now; j++)
			{
				//���кϲ�
				if (pmf.id == this->Quarantine_area[j].id)
				{
					this->Quarantine_area[j].add_new_pmfs(pmf);
					flag = true;
					//�����������,��������ֱ����
					if (this->Quarantine_area[j].isTruncated() == false)
					{
						this->New_complete.push_back(this->Quarantine_area[j]);
						this->Quarantine_area[j].isTruncated_pmf = false;
					}
					break;
				}

			}
			if (flag == false)//���뵱ǰ������
			{
				this->Quarantine_area.push_back(pmf);
			}
		}
	}
	//���ȸ��Ƶ�һ���ӽ���SC
	this->S_C=this->Children_node_list[0].S_C;
	if (this->Children_node_list[0].IsTime == true)//ʱ��ά���ϵľۺ�
	{
		this->T = this->Children_node_list[0].T;
		this->Space_cost += this->Children_node_list[0].Space_cost;
		this->Space_cost_sketch = this->Children_node_list[0].Space_cost_sketch;
		this->Space_data = this->Children_node_list[0].Space_data;
		for (int i = 1; i < this->Children_node_list.size(); i++)
		{	
			Sketh_on_Time(this->S_C, this->Children_node_list[i].S_C,this->T,this->Children_node_list[i].T);
			this->T += this->Children_node_list[i].T;
			this->Space_cost += this->Children_node_list[i].Space_cost;
			this->Space_cost_sketch += this->Children_node_list[i].Space_cost_sketch;
			//this->Space_data += this->Children_node_list[i].Space_data;
		}
	
	}
	else//��׼ά���ϵľۺ�
	{
		this->T = this->Children_node_list[0].T;
		this->Space_cost += this->Children_node_list[0].Space_cost;
		this->Space_cost_sketch = this->Children_node_list[0].Space_cost_sketch;
		this->Space_data = this->Children_node_list[0].Space_data;
		for (int i = 1; i < this->Children_node_list.size(); i++)
		{

		
			Sketh_on_Standard(this->S_C, this->Children_node_list[i].S_C);
			this->Space_cost += this->Children_node_list[i].Space_cost;
			this->Space_cost_sketch += this->Children_node_list[i].Space_cost_sketch;
			//this->Space_data += this->Children_node_list[i].Space_data;

		}
		
	}

	//std::cout << "EK and Eb:   " << this->S_C.E_k << "  " << this->S_C.E_b << endl;
	this->isMaterialize = true;
	this->Use_Sc_New_complete_to_get_SC_in_node();
	this->Use_SC_Quarantine_area_to_get_SA_in_node();
	this->Space_cost_sketch += 5;
	return;
}
void H_node::Materialize_node_by_pws()
{
	if (this->IsLeaf == true)
	{
		this->Use_object_pmf_to_get_line_pmf_on_leaf_node_pws();
		this->isMaterialize = true;
		return;
	}

	for (auto& child_node : this->Children_node_list)
	{
		if (child_node.isMaterialize == false)//�ӽ��û�б��ﻯ
		{
			child_node.Materialize_node_by_pws();
		}
	}
	//�����ӽ�㶼���ﻯ��
	if (this->IsRoot == true)
	{
		this->Space_cost += double(this->Children_node_list[0].F_C.pmf.size()) + double(this->Children_node_list[0].F_C.pmf.size());
		return;
	}
	this->Materialize_Not_leaf_by_pws();
	return;






}
void H_node::Materialize_Not_leaf_by_pws()
{
	if (this->IsLeaf == true)
	{
		cout << "Leaf node,dont use this funcation" << endl;
		return;
	}
	cout << "��ǰ�ﻯcuboidΪ:";
	for (auto s : this->Dimension)
	{
		cout << s << "  ";
	}
	if (this->Children_node_list.size() == 0)
	{
		cout << "********" << endl;
		return;
	}
	//���ݴ���
	if (this->Children_node_list[0].IsTime == true)//ʱ��ά���ϵľۺ�
	{
		for (int i = 0; i < this->Children_node_list.size(); i++)
		{
			if (i == 0)
			{
				this->T = this->Children_node_list[0].T;
			}
			else
			{
				this->T += this->Children_node_list[i].T;
			}
			this->Space_pws += this->Children_node_list[i].Space_pws;
			this->Space_cost += this->Children_node_list[i].Space_cost;
			for (auto pmf_i : this->Children_node_list[i].All_data_instance)
			{
				//��һ��ʱ���Ӻ���
				for (auto &pmf_j : pmf_i.pmfs)
				{
					{
						pmf_j.x += this->T- this->Children_node_list[i].T;
					}
					
				}
				this->All_data_instance.push_back(pmf_i);
			}
		}

	}
	else//��׼ά���ϵľۺ�
	{
		this->T = this->Children_node_list[0].T;
		
		for (int i = 0; i < this->Children_node_list.size(); i++)
		{
			this->Space_pws += this->Children_node_list[i].Space_pws;
			this->Space_cost += this->Children_node_list[i].Space_cost;
			for (auto pmf_i : this->Children_node_list[i].All_data_instance)
			{
				this->All_data_instance.push_back(pmf_i);
			}

		}

	}
int sum_pmf_size = 0;
	for (auto pmf_i : this->All_data_instance)
	{
		sum_pmf_size += pmf_i.pmfs.size();
		//����pmf
		if (pmf_i.isTruncated()==false)
		{
			use_objects_to_get_pmfs_pws(T,pmf_i,F_C);
		}
		else//���������
		{
			Quarantine_area.push_back(pmf_i);
		}
	}
	cout << "F_C�Ĵ�СΪ��" << this->F_C.pmf.size() << endl;
	this->Use_Fc_Quarantine_area_to_get_Fa_in_node();
	this->Space_cost += double(this->F_C.pmf.size()) + double(this->F_A.pmf.size());
	this->Space_pws += double(sum_pmf_size * (4 + 3 * 8));
	isMaterialize = true;
	return;
}
void H_node::query_pmf_DFS()
{
	//������ֵ
	double k_m=0, p_m = 0;
	
	
	double sum_k = 0.0, sum_p = 0.0, sum_p2 = 0.0;
		for (auto line : this->F_C.pmf)
		{
			if (line.k > k_m)
				sum_p += line.p;
			else
			{
				sum_p2+= line.p;
			}
		}
		if (sum_p >= p_m)
		{
				for (auto s : this->Dimension)
				{
					std::cout << s << "  ";
				}
				std::cout << "���쳣cuboid" << endl;
		}
		if (this->IsLeaf == true)
			return;
	for (auto child : Children_node_list)
	{
		child.query_pmf_DFS();
	}
}
void H_node::query_sketch_DFS()
{
	//������ֵ
	double e_k = -0.1;


	//double sum_k = 0.0, sum_p = 0.0;
	
	if (this->S_A.E_k > e_k)
	{
		std::cout << this->Now_Dimenson;
		std::cout << "���쳣cuboid" << endl;
	}
	if (this->IsLeaf == true)
		return;
	for (auto child : Children_node_list)
	{
		child.query_sketch_DFS();
	}
}
void H_node::query_sketch_BFS()
{
	//������ֵ
	double e_k = -0.1;
	//double sum_k = 0.0, sum_p = 0.0;

	if (this->S_A.E_k > e_k)
	{
		std::cout << this->Now_Dimenson;
		std::cout << "���쳣cuboid" << endl;
	}
	return;
}
void H_node::query_pmf_BFS()
{


	//������ֵ
	double k_m = 0, p_m = 0;


	double sum_k = 0.0, sum_p = 0.0, sum_p2 = 0.0;
	for (auto line : this->F_C.pmf)
	{
		if (line.k > k_m)
			sum_p += line.p;
		else
		{
			sum_p2 += line.p;
		}
	}
	if (sum_p >= p_m)
	{
		for (auto s : this->Dimension)
		{
			std::cout << s << "  ";
		}
		std::cout << "���쳣cuboid" << endl;
	}
	return;






}
void H_node::Use_object_pmf_to_get_line_pmf_on_leaf_node_pws()
{
	if (this->IsLeaf == false)
	{
		cout << "Not leaf to get lines' pmf,error!" << endl;
		return;
	}
	cout << "��ǰ�ﻯcuboidΪ:";
	for (auto s : this->Dimension)
	{
		cout << s << "  ";
	}

	//������ǰϸ���е�����pmf,ʹ�����е�����pmf�����㼯fC���ضϵķ��������
	int sum_pmf_size = 0;
	for (auto pmf_i : this->All_data_instance)
	{
		sum_pmf_size += pmf_i.pmfs.size();
		//����pmf
		if (pmf_i.isTruncated() == false)
		{
			use_objects_to_get_pmfs_pws(T, pmf_i, F_C);
		}
		else//���������
		{
			Quarantine_area.push_back(pmf_i);
		}
	}
	cout << "F_C�Ĵ�СΪ��" << this->F_C.pmf.size() << endl;
	this->Use_Fc_Quarantine_area_to_get_Fa_in_node();
	this->Space_cost = double(int(this->F_A.pmf.size()) + (this->F_C.pmf.size()));
	//this->Space_data = double(sum_pmf_size * (4 + 3 * 8));
	isMaterialize = true;
	return;



}
void H_node::agg_cell()
{
	if (this->IsLeaf == true)
	{
		if (this->IsLeaf == false)
		{
			cout << "Not leaf to get lines' pmf,error!" << endl;
			return;
		}
		cout << "��ǰ�ﻯcuboidΪ:";
		for (auto s : this->Dimension)
		{
			cout << s << "  ";
		}

		int k = 0;
		vector< Lines_pmf> F;
		Lines_pmf fc;
		int t = 0;
		bool isfirst = true;
		for (auto pmf_i : this->All_data_instance)
		{
			//Lines_pmf fn;
			k++;
			//fc=use_objects_to_get_pmfs_lattice(this->T, pmf_i, fc, false);
			//isfirst = true;
			/*if (k < 5)
			{
				use_objects_to_get_pmfs1(this->T, pmf_i, fc);
			}*/
			//else 
			use_objects_to_get_pmfs(this->T, pmf_i, fc);
			if (k % 30 == 0)
			{
				
				Lines_pmf fc1;
				fc1.pmf.swap(fc.pmf);
				F.push_back(fc1);
				//fc.pmf.swap(fc1.pmf);
				//fc.Now.clear();
				//fc.pmf.clear();
				//k = 0;
				//fc = fc1;
				//isfirst = false;
				//fc.Now.clear();
			}
		}
		if (fc.pmf.size() > 0)
		{
			F.push_back(fc);
		}
		//this->F_C.pmf.assign(F[0].pmf.begin(), F[0].pmf.end());
		//this->F_C.Now.assign(F[0].Now.begin(), F[0].Now.end());
		this->F_C = F[0];
		for (int i = 1; i < F.size(); i++)
		{
			this->F_C = Convolution_at_Standard(this->F_C, F[i]);
			//this->F_C=Convolution_at_Standard_lattice(this->F_C, F[i]);
		}
		//this->F_C.pmf.clear();
		//this->F_C.pmf.swap(fc.pmf);
		cout << "F_C�Ĵ�СΪ��" << this->F_C.pmf.size() << endl;
		//this->Use_Fc_Quarantine_area_to_get_Fa_in_node();
		this->Space_cost = double(int(this->F_A.pmf.size()) + (this->F_C.pmf.size()));
		//this->Space_data = double(sum_pmf_size * (4 + 3 * 8));
		isMaterialize = true;
		return;
	}
	/*if (this->IsLeaf == true)
	{
		this->Use_object_pmf_to_get_line_pmf_on_leaf_node();
		this->isMaterialize = true;
		return;
	}*/
	for (auto& child_node : this->Children_node_list)
	{
		if (child_node.isMaterialize == false)//�ӽ��û�б��ﻯ
		{
			child_node.agg_cell();
		}
	}
	//�����ӽ�㶼���ﻯ��
	if (this->IsRoot == true)
	{
		this->Space_cost = this->Children_node_list[0].Space_cost;
		return;
	}
	return;
}
void H_node::on_line_sketch(double eM,double varM,double tc)
{
	//this->T = T;
	if ((this->Now_window < tc)||(this->Now_window==1&&tc==1)) {
		this->Now_window++;
		double EK = this->S_C.E_k, EB = this->S_C.E_b, VarK = this->S_C.Var_k, VarB = this->S_C.Var_b;
		double Cov = this->S_C.Cov;
		double K = tc - 1;
		double C1 = double(K - 4) / (K + 2), C2 = -6.0 / ((K + 1) * (K + 2)), C3 = -C2;
		double D1 = 2, D2 = (K + 3) / ((K + 1)), D3 = -2 / ((K + 1));
		this->S_C.E_k = C1 * EK + C2 * EB + C3 * eM;
		this->S_C.E_b = D1 * EK + D2 * EB + D3 * eM;
		this->S_C.Var_k = C1 * C1 * VarK + C2 * C2 * VarB + 2 * C1 * C2 * Cov + C3 * C3 * varM;
		this->S_C.Var_b = D1 * D1 * VarK + D2 * D2 * VarB + 2 * D1 * D2 * Cov + D3 * D3 * varM;
		this->S_C.Cov = (C1 * D2 + C2 * D1) * Cov + C1 * D1 * VarK + C2 * D2 * VarB + C3 * D3 * varM;
	}
	else
	{

		double K = tc;
		double C_1 = 6.0 * (2 * K - this->Now_window - 1) / (pow(this->Now_window, 3) - this->Now_window);
		double C_2 = 2.0 * (2 * this->Now_window - 3 * K + 1) / (this->Now_window * (this->Now_window - 1));
		double EK = this->S_C.E_k, EB = this->S_C.E_b, VarK = this->S_C.Var_k, VarB = this->S_C.Var_b;
		double Cov = this->S_C.Cov;
		this->S_C.E_k =  EK + C_1 * eM;
		this->S_C.E_b =  C_2 * eM;
		this->S_C.Var_k =  VarK + C_1 * C_1 * eM;
		this->S_C.Var_b =  VarB + C_2 * C_2 * eM;
		this->S_C.Cov = Cov + C_1 * C_2  * varM;





	}
}
void H_node::Materialize_node_for_time()
{
	
	if (this->IsLeaf == true)
	{
		double t_start = clock();
		this->Use_object_pmf_to_get_line_pmf_on_leaf_node();
		this->isMaterialize = true;
		double t_end = clock();
		this->Materalize_time += (t_end - t_start) / CLOCKS_PER_SEC;
		return;
	}

	for (auto& child_node : this->Children_node_list)
	{
		if (child_node.isMaterialize == false)//�ӽ��û�б��ﻯ
		{
			child_node.Materialize_node_for_time();
			this->Materalize_time += child_node.Materalize_time;
		}
	}
	//�����ӽ�㶼���ﻯ��
	
	if(this->IsTime==true)
	{
		this->Materalize_time = 10000;
		return;
	}
	double t_start = clock();
	this->Materialize_Not_leaf();
	double t_end = clock();
	this->Materalize_time += (t_end - t_start) / CLOCKS_PER_SEC;
	return;






}

void H_node::on_line_sketch(Data_object o_i, double tc)
{
	double eM = 0.0, varM = 0.0, eM2 = 0.0;
	for (auto pmi : o_i.Data_instance_list)
	{
		eM += pmi.Data_instance_Measure * pmi.Data_instance_prob;
		eM2 += pmi.Data_instance_Measure * pmi.Data_instance_Measure * pmi.Data_instance_prob;
	}
	varM = eM2 - eM * eM;
	//this->T = T;
	if ((this->Now_window < tc) || (this->Now_window == 1 && tc == 1)) {
		this->Now_window++;
		double EK = this->S_C.E_k, EB = this->S_C.E_b, VarK = this->S_C.Var_k, VarB = this->S_C.Var_b;
		double Cov = this->S_C.Cov;
		double K = tc - 1;
		double C1 = double(K - 4) / (K + 2), C2 = -6.0 / ((K + 1) * (K + 2)), C3 = -C2;
		double D1 = 2, D2 = (K + 3) / ((K + 1)), D3 = -2 / ((K + 1));
		this->S_C.E_k = C1 * EK + C2 * EB + C3 * eM;
		this->S_C.E_b = D1 * EK + D2 * EB + D3 * eM;
		this->S_C.Var_k = C1 * C1 * VarK + C2 * C2 * VarB + 2 * C1 * C2 * Cov + C3 * C3 * varM;
		this->S_C.Var_b = D1 * D1 * VarK + D2 * D2 * VarB + 2 * D1 * D2 * Cov + D3 * D3 * varM;
		this->S_C.Cov = (C1 * D2 + C2 * D1) * Cov + C1 * D1 * VarK + C2 * D2 * VarB + C3 * D3 * varM;
	}
	else
	{

		double K = tc;
		double C_1 = 6.0 * (2 * K - this->Now_window - 1) / (pow(this->Now_window, 3) - this->Now_window);
		double C_2 = 2.0 * (2 * this->Now_window - 3 * K + 1) / (this->Now_window * (this->Now_window - 1));
		double EK = this->S_C.E_k, EB = this->S_C.E_b, VarK = this->S_C.Var_k, VarB = this->S_C.Var_b;
		double Cov = this->S_C.Cov;
		this->S_C.E_k = EK + C_1 * eM;
		this->S_C.E_b = C_2 * eM;
		this->S_C.Var_k = VarK + C_1 * C_1 * eM;
		this->S_C.Var_b = VarB + C_2 * C_2 * eM;
		this->S_C.Cov = Cov + C_1 * C_2 * varM;





	}
}
