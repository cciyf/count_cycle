#include "Snap.h"
#include <ctime>
#include <vector>
#include <string>
#include <numeric>
#include <fstream>
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>


const std::string delim = "\t";
const int delimlength = 1;

using namespace TSnap;

/*
 * struct for cycles
 */
struct cyclenum {
	long long num_3_cycles;
	long long num_4_cycles;
	long long num_5_cycles;
	cyclenum(long long _num_3_cycles, long long _num_4_cycles, long long _num_5_cycles) {
		num_3_cycles = _num_3_cycles;
		num_4_cycles = _num_4_cycles;
		num_5_cycles = _num_5_cycles;
	}
};


/*
 * Configuration structure of graph data, 
 * including graph directivity, 
 * number of nodes, number of edges and data path.
 */
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

/**
 * Get the graph edge list file path from the config file and read it to the mat.
 * @param mat Eigen' matrix to store graph data
 * @param config The graph data configuration.
 * @return Data Reading Results
 */
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
		return true;
	}
	return false;
}

/**
 * Get the graph edge matrix file path from the config file and read it to the mat.
 * @param mat Eigen' matrix to store graph data
 * @param config The graph data configuration.
 * @return Data Reading Results
 */
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


/**
 * Get the number of undirected graph 3-cycles, 4-cycles and 5-cycles.
 * Let cn(G) be the number of n-cycles in G.
 * Let A be the adjacency matrix of G.
 * c3(G) = 1 / 6 * tr(A)
 * c4(G) = 1 / 8 * (tr(A^4) - 2 q - 2 ??a_ij(2)(i != j))
 * c5(G) = 1 / 10 * (tr(A^5) - 5 * tr(A^3) - 5 ??(??aij - 2)aii(2)))
 * @param mat Eigen' matrix to store graph data
 * @return cycles' number
 */
cyclenum countUndirectedCycle(Eigen::Matrix<long long int, Eigen::Dynamic, Eigen::Dynamic> &mat) {
	long long cycle4_part2 = mat.sum();
	
	Eigen::Matrix<long long int, Eigen::Dynamic, Eigen::Dynamic> newmat = mat * mat;
	long long cycle4_part3 = newmat.sum() - newmat.trace();
	
	newmat = newmat * mat;
	long long cycle3_part1 = newmat.trace();
	long long cycle5_part2 = cycle3_part1;
	long long cycle5_part3 = ((mat.colwise().sum().array() - 2).matrix() * newmat.diagonal()).sum();
	long long cycle3_num = cycle3_part1 / 6;
	
	newmat = newmat * mat;
	long long cycle4_part1 = newmat.trace();
	long long cycle4_num = (cycle4_part1 - cycle4_part2 - 2 * cycle4_part3) / 8;

	newmat = newmat * mat;
	long long cycle5_part1 = newmat.trace();
	long long cycle5_num = (cycle5_part1 - 5 * cycle5_part2 - 5 * cycle5_part3) / 10;


	//std::cout << "There are " << cycle3_num << " 3-cycle(s)." << std::endl;
	//std::cout << "There are " << cycle4_num << " 4-cycle(s)." << std::endl;
	//std::cout << "There are " << cycle5_num << " 5-cycle(s)." << std::endl;


	return cyclenum(cycle3_num, cycle4_num, cycle5_num);
}

/**
 * Get the number of directed graph 3-cycles, 4-cycles and 5-cycles.
 * @param mat Eigen' matrix to store graph data
 * @return cycles' number
 */
cyclenum countDirectedCycle(Eigen::Matrix<long long int, Eigen::Dynamic, Eigen::Dynamic> &mat) {
	Eigen::Matrix<long long int, Eigen::Dynamic, Eigen::Dynamic> unmat = mat.transpose().array() * mat.array();		
	long long cycle4_part2 = unmat.sum();
	
	Eigen::Matrix<long long int, Eigen::Dynamic, Eigen::Dynamic> newmat = mat * mat;									//A^2
	Eigen::Matrix<long long int, Eigen::Dynamic, Eigen::Dynamic> newunmat = newmat * unmat;
	long long cycle5_part3 = (newunmat.array() * unmat.array()).matrix().sum() - (newunmat.array() * unmat.array()).matrix().trace();
	newunmat.resize(0, 0);

	newmat = newmat * mat;
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
	


	//std::cout << "There are " << cycle3_num << " 3-cycle(s)." << std::endl;
	//std::cout << "There are " << cycle4_num << " 4-cycle(s)." << std::endl;
	//std::cout << "There are " << cycle5_num << " 5-cycle(s)." << std::endl;
	return cyclenum(cycle3_num, cycle4_num, cycle5_num);
}

/**
 * Remove k-in_degree nodes in directed graph
 * @param Graph
 * @param DegK k-degree
 * @return cycles' number vector, each element represents the number of cycles of a SCC
 */
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

/**
 * Remove k-out_degree nodes in directed graph
 * @param Graph
 * @param DegK k-degree
 * @return cycles' number vector, each element represents the number of cycles of a SCC
 */
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

/**
 * Remove k-degree nodes in undirected graph
 * @param Graph
 * @param DegK k-degree
 * @return cycles' number vector, each element represents the number of cycles of a SCC
 */
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

/**
 * Calculate the number of cycles for each strongly connected component
 * @param SCnComV, graph's scc
 * @param is_directed, indicates whether the graph is a directed graph.
 * @param Graph, graph
 * @return cycles' number vector, each element represents the number of cycles of a SCC
 */
template <class PGraph>
std::vector<cyclenum> CntSCCCycle(TCnComV &SCnComV, bool is_directed, PGraph& Graph)
{
	std::vector<cyclenum> cycle_num;

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


				cycle_num.push_back(countUndirectedCycle(mat));
				
			}
		}
	}
	return cycle_num;
	
}


int main() 
{
	// the config data path and result path;
	std::string config_path = "D:/code/C++/count_cycle/cntcycle/Data/web-Stanford_2_config.txt";
	std::string result_path = "D:/code/C++/count_cycle/cntcycle/Data/web-Stanford_2_result.txt";

	// read config data
	struct config config;
	std::vector<cyclenum> cycle_num;
	config.readPath(config_path);

	// write result
	std::ofstream out(result_path);


	if (config.isDirected()) {
		// load edge list
		clock_t start = clock();
		PNGraph G = LoadEdgeList<PNGraph>(config.edge_list_path.c_str());
		clock_t end = clock();
		double elapsed_secs = static_cast<double>(end - start) / CLOCKS_PER_SEC;

		// write result
		std::cout << "Read matrix : " << config_path << std::endl;
		std::cout << "Read matrix cost :" << elapsed_secs << " s." << std::endl;
		out << "Read matrix : " << config_path << std::endl;
		out << "Read matrix cost :" << elapsed_secs << " s." << std::endl;
		
		
		std::cout << "Before remove nodes, there are " << G->GetNodes() << std::endl;
		out << "Before remove nodes, there are " << G->GetNodes() << std::endl;

		// Removes all the node of out-degree and in-degree = 1 or 0.
		start = clock();
		DelInDeg(G, 0);
		DelOutDeg(G, 0);
		DelInDeg(G, 1);
		DelOutDeg(G, 1);
		end = clock();
		elapsed_secs = static_cast<double>(end - start) / CLOCKS_PER_SEC;

		// write result
		std::cout << "After remove nodes, there are " << G->GetNodes() << std::endl;
		out << "After remove nodes, there are " << G->GetNodes() << std::endl;
		std::cout << "Remove nodes cost : " << elapsed_secs << " s." << std::endl;
		out << "Read matrix cost :" << elapsed_secs << " s." << std::endl;

		// Get Sccs and calculate the cycles' number.
		TCnComV SCnComV;
		start = clock();
		GetSccs(G, SCnComV);
		end = clock();
		elapsed_secs = static_cast<double>(end - start) / CLOCKS_PER_SEC;
		std::cout << "There are " << SCnComV.Len() << " strongly connected components." << std::endl;
		std::cout << "There are " << G->GetNodes() << " nodes." << std::endl;
		std::cout << "There are " << G->GetEdges() << " edges." << std::endl;
		std::cout << "Calculate strongly connected components costs: " << elapsed_secs << " s." << std::endl;
		out << "There are " << SCnComV.Len() << " strongly connected components." << std::endl;
		out << "Calculate strongly connected components costs: " << elapsed_secs << " s." << std::endl;

		// Calculate the number of cycles for each strongly connected component
		start = clock();
		cycle_num = CntSCCCycle(SCnComV, config.isDirected(), G);
		end = clock();
		elapsed_secs = static_cast<double>(end - start) / CLOCKS_PER_SEC;
		std::cout << "Count all cycles cost : " << elapsed_secs << " s." << std::endl;
		out << "Count all cycles cost : " << elapsed_secs << " s." << std::endl;
		std::cout << "There are " << accumulate(cycle_num.begin(), cycle_num.end(), 0, [](long long a, cyclenum b) {return a + b.num_3_cycles + b.num_4_cycles + b.num_5_cycles; }) << " cycles in graph." << std::endl;
		out << "There are " << accumulate(cycle_num.begin(), cycle_num.end(), 0, [](long long a, cyclenum b) {return a + b.num_3_cycles + b.num_4_cycles + b.num_5_cycles; }) << " cycles in graph." << std::endl;

		out << "There are " << accumulate(cycle_num.begin(), cycle_num.end(), 0, [](long long a, cyclenum b) {return a + b.num_3_cycles; }) << " 3-cycles in graph." << std::endl;
		out << "There are " << accumulate(cycle_num.begin(), cycle_num.end(), 0, [](long long a, cyclenum b) {return a + b.num_4_cycles; }) << " 4-cycles in graph." << std::endl;
		out << "There are " << accumulate(cycle_num.begin(), cycle_num.end(), 0, [](long long a, cyclenum b) {return a + b.num_5_cycles; }) << " 5-cycles in graph." << std::endl;

		std::cout << "There are " << accumulate(cycle_num.begin(), cycle_num.end(), 0, [](long long a, cyclenum b) {return a + b.num_3_cycles; }) << " 3-cycles in graph." << std::endl;
		std::cout << "There are " << accumulate(cycle_num.begin(), cycle_num.end(), 0, [](long long a, cyclenum b) {return a + b.num_4_cycles; }) << " 4-cycles in graph." << std::endl;
		std::cout << "There are " << accumulate(cycle_num.begin(), cycle_num.end(), 0, [](long long a, cyclenum b) {return a + b.num_5_cycles; }) << " 5-cycles in graph." << std::endl;

	}
	else {
		// load edge list
		clock_t start = clock();
		PUNGraph G = LoadEdgeList<PUNGraph>(config.edge_list_path.c_str());
		clock_t end = clock();
		double elapsed_secs = static_cast<double>(end - start) / CLOCKS_PER_SEC;

		// write result
		std::cout << "Read matrix : " << config_path << std::endl;
		std::cout << "Read matrix cost :" << elapsed_secs << " s." << std::endl;
		out << "Read matrix : " << config_path << std::endl;
		out << "Read matrix cost :" << elapsed_secs << " s." << std::endl;


		std::cout << "Before remove nodes, there are " << G->GetNodes() << std::endl;
		out << "Before remove nodes, there are " << G->GetNodes() << std::endl;

		// Removes all the node of out-degree and in-degree = 1 or 0.
		start = clock();
		DelDeg(G, 1);
		DelDeg(G, 0);
		end = clock();
		elapsed_secs = static_cast<double>(end - start) / CLOCKS_PER_SEC;

		// write result
		std::cout << "After remove nodes, there are " << G->GetNodes() << std::endl;
		out << "After remove nodes, there are " << G->GetNodes() << std::endl;
		std::cout << "Remove nodes cost : " << elapsed_secs << " s." << std::endl;
		out << "Read matrix cost :" << elapsed_secs << " s." << std::endl;


		// Get Sccs and calculate the cycles' number.
		TCnComV SCnComV;
		start = clock();
		GetSccs(G, SCnComV);
		end = clock();
		elapsed_secs = static_cast<double>(end - start) / CLOCKS_PER_SEC;
		std::cout << "There are " << SCnComV.Len() << " strongly connected components." << std::endl;
		std::cout << "Calculate strongly connected components costs: " << elapsed_secs << " s." << std::endl;
		out << "There are " << SCnComV.Len() << " strongly connected components." << std::endl;
		out << "Calculate strongly connected components costs: " << elapsed_secs << " s." << std::endl;

		// Calculate the number of cycles for each strongly connected component
		start = clock();
		cycle_num = CntSCCCycle(SCnComV, config.isDirected(), G);
		end = clock();
		elapsed_secs = static_cast<double>(end - start) / CLOCKS_PER_SEC;
		std::cout << "Count all cycles cost : " << elapsed_secs << " s." << std::endl;
		out << "Count all cycles cost : " << elapsed_secs << " s." << std::endl;
		std::cout << "There are " << accumulate(cycle_num.begin(), cycle_num.end(), 0, [](long long a, cyclenum b) {return a + b.num_3_cycles + b.num_4_cycles + b.num_5_cycles; }) << " cycles in graph." << std::endl;
		out << "There are " << accumulate(cycle_num.begin(), cycle_num.end(), 0, [](long long a, cyclenum b) {return a + b.num_3_cycles + b.num_4_cycles + b.num_5_cycles; }) << " cycles in graph." << std::endl;


		out << "There are " << accumulate(cycle_num.begin(), cycle_num.end(), 0, [](long long a, cyclenum b) {return a + b.num_3_cycles; }) << " 3-cycles in graph." << std::endl;
		out << "There are " << accumulate(cycle_num.begin(), cycle_num.end(), 0, [](long long a, cyclenum b) {return a + b.num_4_cycles; }) << " 4-cycles in graph." << std::endl;
		out << "There are " << accumulate(cycle_num.begin(), cycle_num.end(), 0, [](long long a, cyclenum b) {return a + b.num_5_cycles; }) << " 5-cycles in graph." << std::endl;
		std::cout << "There are " << accumulate(cycle_num.begin(), cycle_num.end(), 0, [](long long a, cyclenum b) {return a + b.num_3_cycles; }) << " 3-cycles in graph." << std::endl;
		std::cout << "There are " << accumulate(cycle_num.begin(), cycle_num.end(), 0, [](long long a, cyclenum b) {return a + b.num_4_cycles; }) << " 4-cycles in graph." << std::endl;
		std::cout << "There are " << accumulate(cycle_num.begin(), cycle_num.end(), 0, [](long long a, cyclenum b) {return a + b.num_5_cycles; }) << " 5-cycles in graph." << std::endl;

	}

	return 0;
}