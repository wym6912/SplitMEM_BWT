/*
Copyright 2016 Uwe Baier, Timo Beller, Enno Ohlebusch

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef INCLUDED_HANDLE_GRAPH
#define INCLUDED_HANDLE_GRAPH

#include <iostream>
#include <fstream>
#include <limits>

using namespace std;

template<class t_node>
void print_graph(const vector<t_node>& graph, const vector<uint64_t>& start_nodes, std::ostream& out=std::cout, std::ostream& out2=std::cout)
{
	// Print start nodes
	for(const auto& node : start_nodes)
	{
		out2 << node << endl;
	}
	uint64_t labels = 0;
	uint64_t edges = 0;
	out << "digraph G {\n";
	uint64_t node_number = 0;
	for(const auto& node : graph)
	{
		labels += graph[node_number].pos_list.size();
		edges += graph[node_number].adj_list.size();
		// Print node_number
		out << "  " << node_number << " ";
		// Print label
		out << "[label=\"";
		if(node.pos_list.size())
		{
			out << node.pos_list[node.pos_list.size()-1];
			for(uint64_t j=node.pos_list.size()-2; j<node.pos_list.size()-1; --j)
			{
				out << "," << node.pos_list[j];
			}
		}
		out << ":" << node.len << "\"]\n";
		// Print edges
		for(uint64_t j=node.adj_list.size()-1; j<node.adj_list.size(); --j)
		{
			out << "  " << node_number << " -> " << node.adj_list[j] << "\n";
		}
		++node_number;
	}
	out << "}" << endl;
	cerr << "nodes=" << graph.size() << endl;
	cerr << "label=" << labels << endl;
	cerr << "edges=" << edges << endl;
}

#endif
