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

#define DEBUG_GRAPH 0

#include <iostream>
#include <fstream>
#include <limits>
#include <algorithm>
#include <tuple>

using namespace std;

#include <sdsl/sdsl_concepts.hpp>

template<class node_t, class int_t>
void dfs(const vector<node_t>& graph, const int_t graph_id, const int_t &now_seq_color, const vector<int_t> &len_, const int_t seqs, std::ostream &out)
{
	if(graph[graph_id].color != now_seq_color) return; // has visited
	auto &this_node = graph[graph_id];
	auto list_size = this_node.pos_list.size();
	this_node.color = numeric_limits<int_t>::max();
	if(list_size >= 2)
	{
		typedef std::tuple<int_t, int_t, int_t> node_info_t;
		// save node list
		vector<node_info_t> node_list(list_size);
		// check duplicate nodes
		map<int_t, int_t> m;
		for(int_t i = 0; i < list_size; ++ i)
		{
			int_t node_number = std::lower_bound(len_.begin(), len_.end(), this_node.pos_list[i]) - len_.begin(), 
				  pre_node_sum = node_number == 0 ? 0 : len_[node_number - 1] + 1;
			node_list[i] = std::make_tuple(this_node.adj_list[i], node_number, this_node.pos_list[i] - 1 - pre_node_sum);
			m[node_number] ++;
		}
		if(m.size() == list_size) // check the nodes are same
		{
			// directly output
			out << this_node.len << ": " << this_node.pos_list.size() << '\n';
			for(auto &data: node_list)
				out << std::get<1>(data) << ' ' << std::get<2>(data) << '\n';
		}
		else
		{
			// cerr << "recount node " << graph_id << " with length = " << this_node.len << endl;
			// re count the adj node
			m.clear();
			for(auto &item: node_list) ++ m[std::get<0>(item)];
			std::sort(node_list.begin(), node_list.end());
			auto node_list_beg = node_list.begin();
			for(auto &item: m)
			{
				// cerr << item.first << " " << item.second << endl;
				if(2 <= item.second && item.second <= len_.size())
				{
#if DEBUG_GRAPH
					std::set<int_t> s; s.clear();
#endif
					while(item.first != std::get<0>(*node_list_beg)) node_list_beg ++;
					out << this_node.len << ": " << item.second << '\n';
					for(int_t i = 0; i < item.second; ++ i, ++ node_list_beg)
					{
#if DEBUG_GRAPH
						if(s.find(std::get<1>(*node_list_beg)) != s.end()) 
						{
							cerr << "Error!!" << std::get<1>(*node_list_beg) << endl;
							exit(1);
						}
#endif
						out << std::get<1>(*node_list_beg) << ' ' << std::get<2>(*node_list_beg) << '\n';
					}
				}
			}
		}
	}
	for(auto &next_id: this_node.adj_list)
	{
		dfs(graph, next_id, now_seq_color, len_, seqs, out);
	}
}

template<class t_node>
void print_graph(const vector<t_node>& graph, const vector<uint64_t>& start_nodes, const vector<uint64_t> &seq_colors, const vector<uint64_t> &cluster_ids, std::ostream& out, sdsl::cache_config& config, vector<uint64_t>& len_)
{
	// len -> prefix sum len
	// +1 means text[i] has 0
	uint64_t seqs_ = len_.size();
	for(uint64_t i = 1; i < seqs_; ++ i) len_[i] += len_[i - 1] + 1;
	// temp data in load file
	sdsl::int_vector<8> text;
	sdsl::load_from_cache(text, conf::KEY_TEXT, config);
	uint64_t node_number;
	std::vector<uint64_t> this_start; this_start.reserve(seqs_);
	for(const auto &this_cluster_id: cluster_ids)
	{
		this_start.clear();
		for(node_number = 0; node_number < seqs_; ++ node_number)
		{
			if(seq_colors[node_number] == this_cluster_id)
			{
				dfs(graph, start_nodes[node_number], this_cluster_id, len_, seqs_, out);
				this_start.emplace_back(node_number);
			}
		}
		if(this_start.size() > 1)
		{
			out << "0:0\n" << this_start.size();
			for(auto &i: this_start) out << ' ' << i;
			out << '\n';
		}
	}
	out << flush;
}

#endif
