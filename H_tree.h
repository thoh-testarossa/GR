#pragma once
#include"H_node.h"
#include"data_instance.h"
#include<stack>
#include"Path.h"
class H_tree {

public:
	H_node Root;//�ӽ��
	H_tree();
	vector<cuboid> Path;
	
	double Space_Cost=0.0;
	bool isMateralize = false;
	bool  isMateralize_sketch = false;
	void Read_stream_instance(Data_instance t);//����instance
	void Read_stream_instance_new_path(Data_instance t);//����instance�ĺ���
	void Materialize_Tree();//�ﻯ����
	void Materialize_Tree_H();
	void Materialize_Tree_by_sketch();//�ﻯ����,��sketch
	void Materialize_Tree_by_PWS();//����PWS�ľۺ�
	void query_pmf();
	void query_sketch();
	void quety_sketch_BFS();
	void query_pmf_BFS();
	void Bulid_path();//����popular-path�ĺ���
	vector<double>  Con_materzation_with_path();//�鿴����ﻯ��popular-path����֮���ϵ�ĺ���
	vector<double>  Sketch_materzation_with_path();////�鿴sketch�ﻯ��popular-path����֮���ϵ�ĺ���
	void  agg_con();////����ۺϵĺ���
	void test_Monitorng(Data_instance t, int mounth, int x, double EM, double VarM, int popular_length);//�����������ĺ���
	void Materalize_for_time();
	void test_MonitorngUS(Data_instance t, int mounth, int x, double EM, double VarM, int popular_length);
	void test_MonitorngClimate(Data_instance t, int mounth, int x, double EM, double VarM, int popular_length);

};
H_tree::H_tree() {
	Root.Dimension.push_back("*");
	Root.IsRoot = true;
};
void H_tree::Read_stream_instance(Data_instance t) {
	H_node *Now_node;
	Now_node= &this->Root;
	
	//��ʱ��ά��,����ʱ��ά
	int k =0 ;
	for (auto &d : t.Time_Dimension)
	{
		bool search = false;
		for(int i=0;i<Now_node->Children_node_list.size();i++)
		{
			//���ڵ�ǰά��
			if (Now_node->Children_node_list[i].Now_Dimenson == d)
			{
				Now_node = &Now_node->Children_node_list[i];
				search = true;
				break;
			}
			
		}
		if (search == false)
		{
			//�¼�һ��cell
			H_node New_node;
			New_node.Now_Dimenson = d;
			if (k == 0)
			{
				New_node.IsTime = true;
				
			}
			New_node.Dimension.assign(Now_node->Dimension.begin(), Now_node->Dimension.end());
			New_node.Dimension.push_back(d);
			Now_node->Children_node_list.push_back(New_node);
			
			Now_node = &Now_node->Children_node_list.back();
			
			
		}
		k++;
	}
	//��׼ά��
	for (auto &d : t.Standard_Dimension)
	{
		bool search = false;
		for (int i = 0; i < Now_node->Children_node_list.size(); i++)
		{
			//���ڵ�ǰά��
			if (Now_node->Children_node_list[i].Now_Dimenson == d)
			{
				Now_node = &Now_node->Children_node_list[i];
				search = true;
				break;
			}

		}
		if (search == false)
		{
			H_node New_node;
			New_node.Now_Dimenson = d;
			New_node.Dimension.assign(Now_node->Dimension.begin(), Now_node->Dimension.end());
			New_node.Dimension.push_back(d);
			Now_node->Children_node_list.push_back(New_node);
			
			Now_node = &Now_node->Children_node_list.back();
		}
		
	}
	Now_node->IsLeaf = true;
	if (Now_node->T < t.Data_instance_x)
	{
		Now_node->T = t.Data_instance_x;
	}
	//Now_node->T = 30;//Ҷ�ӽ��ȫ������Ϊ30
	instance_pmf pmf(t.Data_object_id,t.Data_instance_x,t.Data_instance_Measure,t.Data_instance_prob);
	Now_node->add_pmf(pmf);
};
void H_tree::Materialize_Tree()
{

	this->Root.Materialize_node();
	this->isMateralize = true;


}
void H_tree::Materialize_Tree_H()
{

	this->Root.Materialize_node_Histogram();
	this->isMateralize = true;


}
void H_tree::Materialize_Tree_by_sketch()
{

	this->Root.Materialize_node_by_sketch();
	this->isMateralize_sketch = true;

}
void H_tree::Materialize_Tree_by_PWS()
{

	this->Root.Materialize_node_by_pws();

}
void H_tree::query_pmf()
{
	if (this->isMateralize == false)
	{
		cout << "*****�ﻯδ���*****" << endl;
		return;
	}
	
	this->Root.query_pmf_DFS();






}
void H_tree::query_sketch()
{
	if (this->isMateralize_sketch == false)
	{
		cout << "*****�ﻯδ���*****" << endl;
		return;
	}

	this->Root.query_sketch_DFS();
}
void H_tree::quety_sketch_BFS()
{
	if (this->isMateralize_sketch == false)
	{
		cout << "*****�ﻯδ���*****" << endl;
		return;
	}
	vector<H_node> list;
	list.push_back(this->Root);
	while (list.empty()==false)
	{
		auto Now_node = list[0];
		Now_node.query_sketch_BFS();
		vector<H_node>::iterator k = list.begin();
		list.erase(k);//ɾ����һ��Ԫ��
		if (Now_node.IsLeaf == false)
		{
			for (auto child : Now_node.Children_node_list)
			{
				list.push_back(child);
			}
		}
	}
	
}
void H_tree::query_pmf_BFS()
{
	if (this->isMateralize== false)
	{
		cout << "*****�ﻯδ���*****" << endl;
		return;
	}
	vector<H_node> list;
	list.push_back(this->Root);
	while (list.empty() == false)
	{
		auto Now_node = list[0];
		Now_node.query_pmf_BFS();
		vector<H_node>::iterator k = list.begin();
		list.erase(k);//ɾ����һ��Ԫ��
		if (Now_node.IsLeaf == false)
		{
			for (auto child : Now_node.Children_node_list)
			{
				list.push_back(child);
			}
		}
	}

}
void H_tree::Bulid_path() {
	cuboid root_cuboid;
	root_cuboid.cells.push_back(&this->Root);
	this->Path.push_back(root_cuboid);
	while (true)
	{
		cuboid now_cuboid = Path.back();
		cuboid new_cuboids;
		bool isjump=false;
		for(auto Nowcells : now_cuboid.cells)
		{
			
			for(int i=0;i<Nowcells->Children_node_list.size();i++)
			{

				new_cuboids.cells.push_back(&Nowcells->Children_node_list[i]);
				isjump = Nowcells->Children_node_list[i].IsLeaf;
			}
			
			//Path.push_back(new_cuboids);
		}
		this->Path.push_back(new_cuboids);
		if (isjump == true)
			return;
	}
}
vector<double>  H_tree::Con_materzation_with_path()
{
	vector<double> time_list;
	int length = Path.size() - 1;
	double start_t = clock();
	int i = length;
	{
		for (auto &cells : Path[i].cells)
		{
			cells->Use_object_pmf_to_get_line_pmf_on_leaf_node();
			cells->isMaterialize = true;
		}
	}
	//int path_length = 1;
	
	
	for (int i = length - 1; i > 0; i--)
	{

		for (auto &cells : Path[i].cells)
		{
			cells->Materialize_Not_leaf();
		}
		double end_t = clock();   //����ʱ��
		double time = double(end_t - start_t) / CLOCKS_PER_SEC;
		time_list.push_back(time);
	}
	
	return time_list;



}
vector<double>  H_tree::Sketch_materzation_with_path()
{
	vector<double> time_list;
	int length = Path.size() - 1;
	double start_t = clock();
	int i = length;
	{
		for (auto& cells : Path[i].cells)
		{
			cells->get_sketch_on_leaf_node();
			cells->isMaterialize = true;
		}
	}
	//int path_length = 1;
	
	double start_2 = clock();
	for (int i = length - 1; i > 0; i--)
	{

		for (auto& cells : Path[i].cells)
		{
			cells->Materialize_Not_leaf_by_sketch();
		}
		double end_t = clock();   //����ʱ��
		double time = double(end_t - start_2) / CLOCKS_PER_SEC;
		time_list.push_back(time);
	}
	time_list.push_back(double(start_2 - start_t) / CLOCKS_PER_SEC);
	return time_list;



}
void H_tree::agg_con()
{
	this->Root.agg_cell();
	this->isMateralize = true;
}
/*void H_tree::Read_stream_instance(Data_instance t) {
	H_node *Now_node;
	Now_node= &this->Root;
	
	//��ʱ��ά��,����ʱ��ά
	for (auto &d : t.Time_Dimension)
	{
		bool search = false;
		for(int i=0;i<Now_node->Children_node_list.size();i++)
		{
			//���ڵ�ǰά��
			if (Now_node->Children_node_list[i].Now_Dimenson == d)
			{
				Now_node = &Now_node->Children_node_list[i];
				search = true;
				break;
			}
			
		}
		if (search == false)
		{
			//�¼�һ��cell
			H_node New_node;
			New_node.Now_Dimenson = d;
			New_node.IsTime = true;
			New_node.Dimension.assign(Now_node->Dimension.begin(), Now_node->Dimension.end());
			New_node.Dimension.push_back(d);
			Now_node->Children_node_list.push_back(New_node);

			Now_node = &Now_node->Children_node_list.back();
			
		}
		
	}
	//��׼ά��
	for (auto &d : t.Standard_Dimension)
	{
		bool search = false;
		for (int i = 0; i < Now_node->Children_node_list.size(); i++)
		{
			//���ڵ�ǰά��
			if (Now_node->Children_node_list[i].Now_Dimenson == d)
			{
				Now_node = &Now_node->Children_node_list[i];
				search = true;
				break;
			}

		}
		if (search == false)
		{
			H_node New_node;
			New_node.Now_Dimenson = d;
			New_node.Dimension.assign(Now_node->Dimension.begin(), Now_node->Dimension.end());
			New_node.Dimension.push_back(d);
			Now_node->Children_node_list.push_back(New_node);
			
			Now_node = &Now_node->Children_node_list.back();
		}
		
	}
	Now_node->IsLeaf = true;
	Now_node->T = 30;//Ҷ�ӽ��ȫ������Ϊ30
	instance_pmf pmf(t.Data_object_id,t.Data_instance_x,t.Data_instance_Measure,t.Data_instance_prob);
	Now_node->add_pmf(pmf);
};*/
void H_tree::test_Monitorng(Data_instance t,int mounth,int x,double EM,double VarM,int popular_length) {
	H_node* Now_node;
	Now_node = &this->Root;

	//��ʱ��ά��,����ʱ��ά
	for(int j=0;j<t.Time_Dimension.size();j++)
	{
		string& d = t.Time_Dimension[j];
		bool search = false;
		for (int i = 0; i < Now_node->Children_node_list.size(); i++)
		{
			//���ڵ�ǰά��
			if (Now_node->Children_node_list[i].Now_Dimenson == d)
			{
				Now_node = &Now_node->Children_node_list[i];
				search = true;
				if (popular_length == 7)//˵����ʱ1912��m
				{
					if (j == 0)//���д���
					{
						Now_node->on_line_sketch(EM, VarM, mounth * x);
					}
				}
				break;
			}

		}
		if (search == false)
		{
			//�¼�һ��cell
			H_node New_node;
			New_node.Now_Dimenson = d;
			New_node.IsTime = true;
			New_node.Dimension.assign(Now_node->Dimension.begin(), Now_node->Dimension.end());
			New_node.Dimension.push_back(d);
			Now_node->Children_node_list.push_back(New_node);
			Now_node = &Now_node->Children_node_list.back();
			if (popular_length == 7)//˵����ʱ1912��m
			{
				if (j == 0)//���д���
				{
					Now_node->on_line_sketch(EM, VarM, mounth * x);
					break;
				}
			}
		}

	}
	//��׼ά��
	for(int j=0;j< t.Standard_Dimension.size();j++)
	{
		bool search = false;
		string d = t.Standard_Dimension[j];
		for (int i = 0; i < Now_node->Children_node_list.size(); i++)
		{
			//���ڵ�ǰά��
			if (Now_node->Children_node_list[i].Now_Dimenson == d)
			{
				Now_node = &Now_node->Children_node_list[i];
				search = true;
				if(popular_length==4)
				{
					if (j == 1)
					{
						Now_node->on_line_sketch(EM, VarM, x);
					}
					
				}
				else if (popular_length == 1)
				{
					if (j == 4)
					{
						Now_node->on_line_sketch(EM, VarM, x);
					}
				}
				break;
				
			}

		}
		if (search == false)
		{
			H_node New_node;
			New_node.Now_Dimenson = d;
			New_node.Dimension.assign(Now_node->Dimension.begin(), Now_node->Dimension.end());
			New_node.Dimension.push_back(d);
			Now_node->Children_node_list.push_back(New_node);

			Now_node = &Now_node->Children_node_list.back();

			if (popular_length == 4)
			{
				if (j == 1)
				{
					Now_node->on_line_sketch(EM, VarM, x);
				}

			}
			else if (popular_length == 1)
			{
				if (j == 4)
				{
					Now_node->on_line_sketch(EM, VarM, x);
				}
			}
		}

	}
	Now_node->IsLeaf = true;
	Now_node->T = 30;//Ҷ�ӽ��ȫ������Ϊ30
	Now_node->on_line_sketch(EM, VarM, x);
	//instance_pmf pmf(t.Data_object_id, t.Data_instance_x, t.Data_instance_Measure, t.Data_instance_prob);
	//Now_node->add_pmf(pmf);
};
void H_tree::Materalize_for_time() {

	this->Root.Materialize_node_for_time();

	this->isMateralize = true;









}

void H_tree::Read_stream_instance_new_path(Data_instance t) {
	H_node* Now_node;
	Now_node = &this->Root;
	for (int i = 0; i <2; i++)
	{
		string d = t.Standard_Dimension[i];

		bool search = false;
		for (int i = 0; i < Now_node->Children_node_list.size(); i++)
		{
			//���ڵ�ǰά��
			if (Now_node->Children_node_list[i].Now_Dimenson == d)
			{
				Now_node = &Now_node->Children_node_list[i];
				search = true;
				break;
			}

		}
		if (search == false)
		{
			H_node New_node;
			New_node.Now_Dimenson = d;
			New_node.Dimension.assign(Now_node->Dimension.begin(), Now_node->Dimension.end());
			New_node.Dimension.push_back(d);
			Now_node->Children_node_list.push_back(New_node);

			Now_node = &Now_node->Children_node_list.back();
		}




	}
	
		string d = t.Time_Dimension[0];
		bool search = false;
		for (int i = 0; i < Now_node->Children_node_list.size(); i++)
		{
			//���ڵ�ǰά��
			if (Now_node->Children_node_list[i].Now_Dimenson == d)
			{
				Now_node = &Now_node->Children_node_list[i];
				search = true;
				break;
			}

		}
		if (search == false)
		{
			//�¼�һ��cell
			H_node New_node;
			New_node.Now_Dimenson = d;
			New_node.Dimension.assign(Now_node->Dimension.begin(), Now_node->Dimension.end());
			New_node.Dimension.push_back(d);
			Now_node->Children_node_list.push_back(New_node);

			Now_node = &Now_node->Children_node_list.back();

		}
		{
			d = t.Time_Dimension[1];
			bool search = false;
			for (int i = 0; i < Now_node->Children_node_list.size(); i++)
			{
				//���ڵ�ǰά��
				if (Now_node->Children_node_list[i].Now_Dimenson == d)
				{
					Now_node = &Now_node->Children_node_list[i];
					search = true;
					break;
				}

			}
			if (search == false)
			{
				//�¼�һ��cell
				H_node New_node;
				New_node.Now_Dimenson = d;
				New_node.Dimension.assign(Now_node->Dimension.begin(), Now_node->Dimension.end());
				New_node.Dimension.push_back(d);
				Now_node->Children_node_list.push_back(New_node);

				Now_node = &Now_node->Children_node_list.back();

			}

		}
	
	//��׼ά��
	for (int i = 2; i < 5; i++)
	{
		d = t.Standard_Dimension[i];
		bool search = false;
		for (int i = 0; i < Now_node->Children_node_list.size(); i++)
		{
			//���ڵ�ǰά��
			if (Now_node->Children_node_list[i].Now_Dimenson == d)
			{
				Now_node = &Now_node->Children_node_list[i];
				search = true;
				break;
			}

		}
		if (search == false)
		{
			H_node New_node;
			New_node.Now_Dimenson = d;
			New_node.Dimension.assign(Now_node->Dimension.begin(), Now_node->Dimension.end());
			New_node.Dimension.push_back(d);
			Now_node->Children_node_list.push_back(New_node);

			Now_node = &Now_node->Children_node_list.back();
		}

	}
	//��ʱ��ά��,����ʱ��ά
	//for (auto& d : t.Time_Dimension)
	
	for (int i = 5; i < t.Standard_Dimension.size(); i++)
	{
		string d = t.Standard_Dimension[i];

		bool search = false;
		for (int i = 0; i < Now_node->Children_node_list.size(); i++)
		{
			//���ڵ�ǰά��
			if (Now_node->Children_node_list[i].Now_Dimenson == d)
			{
				Now_node = &Now_node->Children_node_list[i];
				search = true;
				break;
			}

		}
		if (search == false)
		{
			H_node New_node;
			New_node.Now_Dimenson = d;
			New_node.Dimension.assign(Now_node->Dimension.begin(), Now_node->Dimension.end());
			New_node.Dimension.push_back(d);
			Now_node->Children_node_list.push_back(New_node);

			Now_node = &Now_node->Children_node_list.back();
		}




	}
	Now_node->IsLeaf = true;
	/*if (Now_node->T < t.Data_instance_x)
	{
		Now_node->T = t.Data_instance_x;
	}*/
	Now_node->T = 30;//Ҷ�ӽ��ȫ������Ϊ30
	instance_pmf pmf(t.Data_object_id, t.Data_instance_x, t.Data_instance_Measure, t.Data_instance_prob);
	Now_node->add_pmf(pmf);
};
void H_tree::test_MonitorngUS(Data_instance t, int mounth, int x, double EM, double VarM, int popular_length) {
	H_node* Now_node;
	Now_node = &this->Root;

	//��ʱ��ά��,����ʱ��ά
	for (int j = 0; j < t.Time_Dimension.size(); j++)
	{
		string& d = t.Time_Dimension[j];
		bool search = false;
		for (int i = 0; i < Now_node->Children_node_list.size(); i++)
		{
			//���ڵ�ǰά��
			if (Now_node->Children_node_list[i].Now_Dimenson == d)
			{
				Now_node = &Now_node->Children_node_list[i];
				search = true;
				if (popular_length == 4)//˵����ʱ��
				{
					if (j == 0)//���д���
					{
						Now_node->on_line_sketch(EM, VarM, mounth * x);
					}
				}
				else if (popular_length == 3)
				{
					if (j == 1)//���д���
					{
						Now_node->on_line_sketch(EM, VarM, x);
					}
				}
				break;
			}

		}
		if (search == false)
		{
			//�¼�һ��cell
			H_node New_node;
			New_node.Now_Dimenson = d;
			New_node.IsTime = true;
			New_node.Dimension.assign(Now_node->Dimension.begin(), Now_node->Dimension.end());
			New_node.Dimension.push_back(d);
			Now_node->Children_node_list.push_back(New_node);
			Now_node = &Now_node->Children_node_list.back();
			if (popular_length == 3)//˵����ʱ��
			{
				if (j == 0)//���д���
				{
					Now_node->on_line_sketch(EM, VarM, mounth * x);
				}
			}
			else if (popular_length == 2)
			{
				if (j == 1)//���д���
				{
					Now_node->on_line_sketch(EM, VarM, x);
				}
			}
		}

	}
	//��׼ά��
	for (int j = 0; j < t.Standard_Dimension.size(); j++)
	{
		bool search = false;
		string d = t.Standard_Dimension[j];
		for (int i = 0; i < Now_node->Children_node_list.size(); i++)
		{
			//���ڵ�ǰά��
			if (Now_node->Children_node_list[i].Now_Dimenson == d)
			{
				Now_node = &Now_node->Children_node_list[i];
				search = true;
				if (popular_length==1)
				{
					if(j==0)
					Now_node->on_line_sketch(EM, VarM, x);

				}
				break;

			}

		}
		if (search == false)
		{
			H_node New_node;
			New_node.Now_Dimenson = d;
			New_node.Dimension.assign(Now_node->Dimension.begin(), Now_node->Dimension.end());
			New_node.Dimension.push_back(d);
			Now_node->Children_node_list.push_back(New_node);

			Now_node = &Now_node->Children_node_list.back();
			if (popular_length == 1)
			{
				if (j == 0)
					Now_node->on_line_sketch(EM, VarM, x);

			}
		}

	}
	Now_node->IsLeaf = true;
	Now_node->T = 30;//Ҷ�ӽ��ȫ������Ϊ30
	Now_node->on_line_sketch(EM, VarM, x);
	//instance_pmf pmf(t.Data_object_id, t.Data_instance_x, t.Data_instance_Measure, t.Data_instance_prob);
	//Now_node->add_pmf(pmf);
};

void H_tree::test_MonitorngClimate(Data_instance t, int mounth, int x, double EM, double VarM, int popular_length) {
	H_node* Now_node;
	Now_node = &this->Root;

	//��ʱ��ά��,����ʱ��ά
	for (int j = 0; j < t.Time_Dimension.size(); j++)
	{
		string& d = t.Time_Dimension[j];
		bool search = false;
		for (int i = 0; i < Now_node->Children_node_list.size(); i++)
		{
			//���ڵ�ǰά��
			if (Now_node->Children_node_list[i].Now_Dimenson == d)
			{
				Now_node = &Now_node->Children_node_list[i];
				search = true;
				if (popular_length == 4)//
				{
					if (j == 0)//���д���
					{
						Now_node->on_line_sketch(EM, VarM, mounth * x);
						//Now_node->on_line_sketch(o_i, mounth * x);
						
					}
				}
				else if (popular_length == 3)
				{
					if (j == 1)//���д���
					{
						Now_node->on_line_sketch(EM, VarM, x);
						//Now_node->on_line_sketch(o_i, x);
					}
				}
				break;
			}

		}
		if (search == false)
		{
			//�¼�һ��cell
			H_node New_node;
			New_node.Now_Dimenson = d;
			New_node.IsTime = true;
			New_node.Dimension.assign(Now_node->Dimension.begin(), Now_node->Dimension.end());
			New_node.Dimension.push_back(d);
			Now_node->Children_node_list.push_back(New_node);
			Now_node = &Now_node->Children_node_list.back();
			if (popular_length == 4)//˵����ʱ��
			{
				if (j == 0)//���д���
				{
					Now_node->on_line_sketch(EM, VarM, mounth * x);
				
					//Now_node->on_line_sketch(o_i, mounth * x);
				}
			}
			else if (popular_length == 3)
			{
				if (j == 1)//���д���
				{
					Now_node->on_line_sketch(EM, VarM, x);
					//Now_node->on_line_sketch(o_i, x);
				}
			}
		}

	}
	//��׼ά��
	for (int j = 0; j < t.Standard_Dimension.size(); j++)
	{
		bool search = false;
		string d =t.Standard_Dimension[j];
		for (int i = 0; i < Now_node->Children_node_list.size(); i++)
		{
			//���ڵ�ǰά��
			if (Now_node->Children_node_list[i].Now_Dimenson == d)
			{
				Now_node = &Now_node->Children_node_list[i];
				search = true;
				if (popular_length == 2)
				{
					if (j == 0)
					{
						Now_node->on_line_sketch(EM, VarM, x);
					}
						//Now_node->on_line_sketch(EM, VarM, x);
						//Now_node->on_line_sketch(o_i, x);
						//Now_node->on_line_sketch(EM, VarM, x);

				}
				if (popular_length == 1)
				{
					if (j == 1)
					{
						Now_node->on_line_sketch(EM, VarM, x);
					}
						//Now_node->on_line_sketch(o_i, x);
						//Now_node->on_line_sketch(EM, VarM, x);

				}
				break;

			}

		}
		if (search == false)
		{
			H_node New_node;
			New_node.Now_Dimenson = d;
			New_node.Dimension.assign(Now_node->Dimension.begin(), Now_node->Dimension.end());
			New_node.Dimension.push_back(d);
			Now_node->Children_node_list.push_back(New_node);

			Now_node = &Now_node->Children_node_list.back();
			if (popular_length == 2)
			{
				if (j == 0)
				{
					Now_node->on_line_sketch(EM, VarM, x);
				}
					//Now_node->on_line_sketch(o_i, x);
					//Now_node->on_line_sketch(EM, VarM, x);

			}
			if (popular_length == 1)
			{
				if (j == 1)
				{
					Now_node->on_line_sketch(EM, VarM, x);
				}
					//Now_node->on_line_sketch(o_i, x);
					//Now_node->on_line_sketch(EM, VarM, x);

			}
		}

	}
	Now_node->IsLeaf = true;
	Now_node->T = 30;//Ҷ�ӽ��ȫ������Ϊ30
	Now_node->on_line_sketch(EM, VarM, x);
	//instance_pmf pmf(t.Data_object_id, t.Data_instance_x, t.Data_instance_Measure, t.Data_instance_prob);
	//Now_node->add_pmf(pmf);
};
