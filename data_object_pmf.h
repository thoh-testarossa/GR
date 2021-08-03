#pragma once
#include"data_instance.h"
class Data_object {

public:
	int Data_object_id = -1;
	vector<Data_instance> Data_instance_list;
	void Add_new_instance(Data_instance t);
	bool IsTruncated();
};

bool Data_object::IsTruncated() {
	int sum = 0;
	for (auto t : Data_instance_list)
	{
		sum += t.Data_instance_prob;
	}
	double sum1 = sum * 100000;
	return sum1 < 99999;


};
void Data_object::Add_new_instance(Data_instance t) {
	if (Data_object_id < 0)
		Data_object_id = t.Data_object_id;
	this->Data_instance_list.push_back(t);
};