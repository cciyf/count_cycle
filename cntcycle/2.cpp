//#include "Snap.h"
//
//
//
//// Get connected components on directed graph that is not fully connected
//void GetCnCom() {
//
//	PNGraph G;
//
//	// Create benchmark graph
//	const int NNodes = 90;
//	const int NEdges = 200;
//	G = GenRndGnm<PNGraph>(NNodes, NEdges);
//
//	printf("IsConnected(G) = %s\n", (IsConnected(G) ? "True" : "False"));
//	printf("IsWeaklyConnected(G) = %s\n", (IsWeaklyConn(G) ? "True" : "False"));
//
//	TCnComV SCnComV;
//	GetSccs(G, SCnComV);
//	for (int i = 0; i < SCnComV.Len(); i++) {
//		printf("SCnComV[%d].Len() = %d\n", i, SCnComV[i].Len());
//	}
//}
//
//int main() {
//	PNGraph G = TNGraph::New();
//
//}
