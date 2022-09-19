#include "Snap.h"
#include <ctime>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>


const std::string delim = "\t";
const int delimlength = 1;

using namespace TSnap;

struct config {
	bool is_directed = false;
	int num_nodes;
	int num_edges;
	std::string edge_list_path;

	config() {
		is_directed = false;
		num_nodes = 0;
		num_edges = 0;
		edge_list_path = "";
	}

	bool isDirected() {
		return is_directed;
	}

	bool readPath(std::string config_path) {
		if (config_path.empty()) return false;
		
		std::string line;
		std::istream *config_file = new std::ifstream(config_path.c_str());
		while (std::getline(*config_file, line)) {
			if (line.find("-n ") != std::string::npos)
			{
				num_nodes = std::stoi(line.substr(line.find_first_of("0123456789")));
			}
			else if (line.find("-m ") != std::string::npos)
			{
				num_edges = std::stoi(line.substr(line.find_first_of("0123456789")));
			}
			else if (line.find("-d ") != std::string::npos)
			{
				is_directed = std::stoi(line.substr(line.find_first_of("0123456789")));
			}
			
		}

		edge_list_path = config_path.substr(0, config_path.size() - 11);
		edge_list_path = edge_list_path + ".txt";
	}
};

bool readEdgeList(Eigen::Matrix<long long int, Eigen::Dynamic, Eigen::Dynamic> &mat, config &config) {

	std::string line;
	std::ifstream fin(config.edge_list_path);
	mat.setZero();
	if (fin) {
		for (int i = 0; getline(fin, line); i++) {
			std::cout << line << std::endl;
			int a = line.find_first_of(delim);
			
			mat(atoi(line.substr(a + delimlength).c_str()), atoi(line.substr(0, a).c_str())) = 1;
		}
	}
}

bool readMatrixData(Eigen::Matrix<long long int, Eigen::Dynamic, Eigen::Dynamic> &mat, config &_config) {
	int value;
	std::string line;
	//std::istream *graph_file;
	std::ifstream fin(_config.edge_list_path);


	//mat.resize(_config.num_nodes, _config.num_nodes);
	if (fin) {
		for (int i = 0; getline(fin, line); i++) {
			std::stringstream input(line);
			for (int j = 0; input >> value; ++j) {
				mat(j, i) = value;
			}
			if (i % 100 == 0) std::cout << i << std::endl;
		}
	}
	return true;
}



long long countUndirectedCycle(Eigen::Matrix<long long int, Eigen::Dynamic, Eigen::Dynamic> &mat) {
	long long cycle4_part2 = mat.sum();
	
	std::cout << cycle4_part2 << std::endl;
	Eigen::Matrix<long long int, Eigen::Dynamic, Eigen::Dynamic> newmat = mat * mat;
	long long cycle4_part3 = newmat.sum() - newmat.trace();
	
	std::cout << cycle4_part3 << std::endl;
	newmat = newmat * mat;
	long long cycle3_part1 = newmat.trace();
	long long cycle5_part2 = cycle3_part1;
	long long cycle5_part3 = ((mat.colwise().sum().array() - 2).matrix() * newmat.diagonal()).sum();
	long long cycle3_num = cycle3_part1 / 6;
	
	newmat = newmat * mat;
	long long cycle4_part1 = newmat.trace();
	long long cycle4_num = (cycle4_part1 - 2 * cycle4_part2 - 2 * cycle4_part3) / 8;

	newmat = newmat * mat;
	long long cycle5_part1 = newmat.trace();
	long long cycle5_num = (cycle5_part1 - 5 * cycle5_part2 - 5 * cycle5_part3) / 10;


	std::cout << "��Ԫ����" << cycle3_num << "��" << std::endl;
	std::cout << "��Ԫ����" << cycle4_num << "��" << std::endl;
	std::cout << "��Ԫ����" << cycle5_num << "��" << std::endl;
	return cycle3_num + cycle4_num + cycle5_num;
	
}

long long countDirectedCycle(Eigen::Matrix<long long int, Eigen::Dynamic, Eigen::Dynamic> &mat) {
	Eigen::Matrix<long long int, Eigen::Dynamic, Eigen::Dynamic> unmat = mat.transpose().array() * mat.array();		
	long long cycle4_part2 = unmat.sum();
	
	Eigen::Matrix<long long int, Eigen::Dynamic, Eigen::Dynamic> newmat = mat * mat;									//A^2
	Eigen::Matrix<long long int, Eigen::Dynamic, Eigen::Dynamic> newunmat = newmat * unmat;
	long long cycle5_part3 = (newunmat.array() * unmat.array()).matrix().sum() - (newunmat.array() * unmat.array()).matrix().trace();
	newunmat.resize(0, 0);

	newmat = newmat * mat;
	std::cout << newmat.diagonal() << std::endl;
	long long cycle5_part2 = (unmat.colwise().sum() * newmat.diagonal()).sum();
	long long cycle3_part1 = newmat.trace();

	unmat = unmat * unmat;
	long long cycle4_part3 = unmat.sum() - unmat.trace();

	newmat = newmat * mat;
	long long cycle4_part1 = newmat.trace();

	newmat = newmat * mat;
	long long cycle5_part1 = newmat.trace();
	
	long long cycle3_num = cycle3_part1 / 3;
	long long cycle4_num = (cycle4_part1 - cycle4_part2 - 2 * cycle4_part3) / 4;
	long long cycle5_num = (cycle5_part1 - 5 * (cycle5_part2 - cycle5_part3)) / 5;
	


	std::cout << "��Ԫ����" << cycle3_num << "��" << std::endl;
	std::cout << "��Ԫ����" << cycle4_num << "��" << std::endl;
	std::cout << "��Ԫ����" << cycle5_num << "��" << std::endl;
	return cycle3_num + cycle4_num + cycle5_num;

}

template <class PGraph>
void DelInDeg(PGraph& Graph, const int& InDegK) {
	TIntV DelNIdV;
	for (typename PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
		if (NI.GetOutDeg() == 0 && NI.GetInDeg() == InDegK) {
			DelNIdV.Add(NI.GetId());
		}
	}
	for (int i = 0; i < DelNIdV.Len(); i++) {
		Graph->DelNode(DelNIdV[i]);
	}
}

template <class PGraph>
void DelOutDeg(PGraph& Graph, const int& OutDegK) {
	TIntV DelNIdV;
	for (typename PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
		if (NI.GetOutDeg() == OutDegK && NI.GetInDeg() == 0) {
			DelNIdV.Add(NI.GetId());
		}
	}
	for (int i = 0; i < DelNIdV.Len(); i++) {
		Graph->DelNode(DelNIdV[i]);
	}
}

template <class PGraph>
void DelDeg(PGraph& Graph, const int& DegK) {
	TIntV DelNIdV;
	for (typename PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
		if (NI.GetDeg() == DegK) {
			DelNIdV.Add(NI.GetId());
		}
	}
	for (int i = 0; i < DelNIdV.Len(); i++) {
		Graph->DelNode(DelNIdV[i]);
	}
}

template <class PGraph>
std::vector<long long int> CntSCCCycle(TCnComV &SCnComV, bool is_directed, PGraph& Graph)
{
	std::vector<long long int> cycle_num;

	if (is_directed) 
	{
		for (int i = 0; i < SCnComV.Len(); i++) 
		{
			if (SCnComV[i].Len() >= 5)
			{
				PGraph sub_graph = GetSubGraph(Graph, SCnComV[i].NIdV, true);
				Eigen::Matrix<long long int, Eigen::Dynamic, Eigen::Dynamic> mat(SCnComV[i].Len(), SCnComV[i].Len());
				mat.setZero();
				for (auto EI = sub_graph->BegEI(); EI < sub_graph->EndEI(); EI++) {
					mat(EI.GetSrcNId(), EI.GetDstNId()) = 1;
				}
				cycle_num.push_back(countDirectedCycle(mat));
			}
		}
	}
	else
	{
		for (int i = 0; i < SCnComV.Len(); i++) {
			if (SCnComV[i].Len() >= 5)
			{
				PGraph sub_graph = GetSubGraph(Graph, SCnComV[i].NIdV, true);
				Eigen::Matrix<long long int, Eigen::Dynamic, Eigen::Dynamic> mat(SCnComV[i].Len(), SCnComV[i].Len());
				mat.setZero();
				for (auto EI = sub_graph->BegEI(); EI < sub_graph->EndEI(); EI++) {
					mat(EI.GetSrcNId(), EI.GetDstNId()) = 1;
					mat(EI.GetDstNId(), EI.GetSrcNId()) = 1;
				}
				cycle_num.push_back(countDirectedCycle(mat));
			}
		}
	}
	return cycle_num;
	
}


int main() 
{

	std::string config_path = "D:\\code\\C++\\count_cycle\\cntcycle\\Data\\email-Enron_config.txt";
	std::string result_path = "D:\\code\\C++\\count_cycle\\cntcycle\\Data\\email-Enron_result.txt";

	struct config config;
	std::vector<long long> cycle_num;
	config.readPath(config_path);

	if (config.isDirected()) {
		// load edge list
		PNGraph G = LoadEdgeList<PNGraph>(config.edge_list_path.c_str());
		
		std::cout << G->GetNodes() << std::endl;
		// Removes all the node of out-degree and in-degree = 1 or 0.
		DelInDeg(G, 0);
		DelOutDeg(G, 0);
		DelInDeg(G, 1);
		DelOutDeg(G, 1);
		std::cout << G->GetNodes() << std::endl;

		// Get Scc
		TCnComV SCnComV;
		GetSccs(G, SCnComV);
		cycle_num = CntSCCCycle(SCnComV, config.isDirected(), G);
	}
	else {
		// load edge list
		PUNGraph G = LoadEdgeList<PUNGraph>(config.edge_list_path.c_str(), 0, 1);

		std::cout << G->GetNodes() << std::endl;
		// Removes all the node of out-degree and in-degree = 1 or 0.
		DelDeg(G, 1);
		DelDeg(G, 0);
		std::cout << G->GetNodes() << std::endl;

		// Get Scc
		TCnComV SCnComV;
		GetSccs(G, SCnComV);
		cycle_num = CntSCCCycle(SCnComV, config.isDirected(), G);
	}

	for (auto i : cycle_num) {
		std::cout << i << std::endl;
	}


	//Eigen::Matrix<long long int, Eigen::Dynamic, Eigen::Dynamic> mat(config.num_nodes + 1, config.num_nodes + 1);
	//readEdgeList(mat, config);
	//std::cout << mat << std::endl;
	//cout << config.num_nodes << endl;


	//clock_t start = clock();
	//Matrix<long long int, Dynamic, Dynamic> mat(config.num_nodes, config.num_nodes);
	//readMatrixData(mat, config);
	//clock_t end = clock();
	//double elapsed_secs = static_cast<double>(end - start) / CLOCKS_PER_SEC;

	//ofstream out(result_path);

	//out << "��ȡ���󻨷�:" << elapsed_secs << " s" << endl;

	//start = clock();
	//long long cnt = 0;
	//if (!config.is_directed) {
	//	cnt = countUndirectedCycle(mat);
	//}
	//else {
	//	cnt = countDirectedCycle(mat);
	//}
	//end = clock();
	//elapsed_secs = static_cast<double>(end - start) / CLOCKS_PER_SEC;
	//out << "����3��4��5Ԫ������ʱ��:" << elapsed_secs << " s" << endl;
	//out << "����3��4��5Ԫ������Ϊ:" << cnt << "����" << endl;

	return 0;
}