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

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <limits>
#include <sdsl/sdsl_concepts.hpp>
#include <sdsl/csa_wt.hpp>
#include <sdsl/suffix_array_algorithm.hpp>
#include "create_datastructures.hpp"
#include "partial_lcp.hpp"
#include "handle_graph.hpp"

#if (defined(_WIN32) || defined(_WIN64))
#include "../include/getopt9/include/getopt.h"
#else
#include <getopt.h>
#endif

using namespace std;
using namespace sdsl;
using namespace std::chrono;

void version(char *name)
{
	cerr << name << " version 1.2.1.0" << endl;
	exit(0);
}

void usage(char *name)
{
	cerr << name << "\n\tUsage: " << name << " [options]\n\n";
	cerr << "Available options:\n";
	cerr << "\t--inputfile/-i    FILE  input file name (required)\n";
	cerr << "\t--outputfile/-o   FILE  output file name (required)\n";
	cerr << "\t--kmer/-k         k     K-mer size (required)\n";
	cerr << "\t--method/-m             use method to construct SA {SAIS, DIV, PARDIV} (default: DIV)\n";
	cerr << "\t--threads/-t      t     use t threads to construct SA (only effected in PARDIV)\n";
	cerr << "\t--tmpdir/-d       DIR   temp dir name (default: \".\")\n";
	cerr << "\t--nocache/-n            do not use cache in program\n";
	cerr << "\t--help/-h               print help message\n";
	cerr << "\t--version/-v            print program version\n";
	cerr << "Sample usage: " << name << " -i 1.fasta -k 10 -o 1.ans" << endl;

	exit(EXIT_FAILURE);
}

struct node
{
	uint64_t len;
	mutable uint64_t color;
	vector<uint64_t> adj_list;
	vector<uint64_t> pos_list;
	node()
	{
		len = color = 0;
	}
};

struct dsu
{
	vector<uint64_t> pa;
  	explicit dsu(uint64_t size) : pa(size) { iota(pa.begin(), pa.end(), 0); }
  	uint64_t find(uint64_t x) { return pa[x] == x ? x : pa[x] = find(pa[x]); }
	uint64_t rev_find(uint64_t x) { return pa[pa.size() - x - 1]; }
	void     merge(uint64_t x, uint64_t y) { pa[find(x)] = find(y); }
};

auto create_cdbg_with_lf(cache_config& config, uint64_t k, uint64_t seqs)
{
	// Create WT of the BWT
	typedef wt_huff<bit_vector, rank_support_v<>, select_support_scan<1>, select_support_scan<0>> wt;
	wt wt_bwt;
	construct(wt_bwt, cache_file_name(conf::KEY_BWT, config));

	// Create C-array (needed for interval_symbols)
	vector<uint64_t> carray(256, 0);
	for(uint64_t i=0, sum=0; i<256; ++i)
	{
		carray[i] = sum;
		sum += wt_bwt.rank(wt_bwt.size(), i);
	}

	// Create bit-vectors
	bit_vector bv(wt_bwt.size(), 0);
	bit_vector bv2(wt_bwt.size(), 0);
	bit_vector bv3(wt_bwt.size(), 0);
	{
		auto start = high_resolution_clock::now();
		int_vector<2> lcp_k = construct_partial_lcp<wt>(wt_bwt, carray, k);
		auto stop = high_resolution_clock::now();
		cerr << std::setw(10) << duration_cast<milliseconds>(stop-start).count() << "ms for partial LCP construction (sdsl-bblaca)" << endl;
		start = high_resolution_clock::now();
		int_vector_buffer<8> bwt(cache_file_name(conf::KEY_BWT, config));
		bool open=false;
		uint64_t kvalue=0;
		uint64_t lb=0;
		uint64_t last_change=0;
		vector<uint64_t> occ(256, 0);
		vector<uint8_t> occ_list;
		occ_list.reserve(256);
		vector<uint64_t> lf = carray;
		for(uint64_t i=1; i<lcp_k.size(); ++i)
		{
			++lf[bwt[i-1]];
			if(lcp_k[i] == gt_k or lcp_k[i] == eq_k)
			{
				open = true;
				if(lcp_k[i] == eq_k)
				{
					kvalue = i;
				}
			}
			else
			{
				if(open)
				{
					if(kvalue > lb)
					{
						bv[lb] = true;
						bv[i-1] = true;
					}
					if(last_change > lb)
					{
						for(uint64_t j=lb; j<=i-1; ++j)
						{
							bv2[j] = true;
							uint8_t c = bwt[j];
							if(occ[c]==0)
							{
								occ[c] = j;
								occ_list.emplace_back(c);
							}
						}
						for(uint64_t j=0; j<occ_list.size(); ++j)
						{
							bv3[lf[occ_list[j]]-1] = true;
							occ[occ_list[j]] = 0;
						}
						occ_list.resize(0);
					}
					open = false;
				}
				lb = i;
			}
			if(bwt[i] != bwt[i-1] or bwt[i] <= 1)
			{
				last_change = i;
			}
		}
		if(open)
		{
			++lf[bwt[lcp_k.size()-1]];
			if(kvalue > lb)
			{
				bv[lb] = true;
				bv[lcp_k.size()-1] = true;
			}
			if(last_change > lb)
			{
				for(uint64_t j=lb; j<=lcp_k.size()-1; ++j)
				{
					bv2[j] = true;
					uint8_t c = bwt[j];
					if(occ[c]==0)
					{
						occ[c] = j;
						occ_list.emplace_back(c);
					}
				}
				for(const auto& c : occ_list)
				{
					bv3[lf[c]-1] = true;
					occ[c] = 0;
				}
				occ_list.resize(0);
			}
		}
		bv3[0] = 1;
		for(uint64_t i=1; i<carray[2]; ++i)
		{
			bv3[i] = 0;
		}
		open = false;
		for(uint64_t i=0; i<bv.size(); ++i)
		{
			if(open)
			{
				bv3[i] = 0;
				if(bv[i])
				{
					open = false;
				}
			}
			else if(bv[i])
			{
				open = true;
				bv3[i] = 0;
			}
		}
		stop = high_resolution_clock::now();
		cerr << std::setw(10) << duration_cast<milliseconds>(stop-start).count() << "ms for create bv and initial nodes" << endl;
	}

	// Init rank support and graph
	bit_vector::rank_1_type bv_rank, bv3_rank;
	util::init_support(bv_rank, &bv);
	util::init_support(bv3_rank, &bv3);
	uint64_t number_right_max_nodes = bv_rank(bv.size())/2, now_seq;
	vector<node> graph(number_right_max_nodes+bv3_rank(bv3.size()));
	vector<uint64_t> start_nodes;
	start_nodes.reserve(seqs);
	dsu seq_colors(seqs);
	set<uint64_t> s;
	// Create compressed de bruijn graph
	{
		auto start = high_resolution_clock::now();
		uint64_t lb = 0;
		uint64_t cur_node = number_right_max_nodes;
		graph[cur_node].len = 1;
		now_seq = 0;
		for(uint64_t i=wt_bwt.size()-1; i>0; --i)
		{
			// LF
			auto res = wt_bwt.inverse_select(lb);
			uint64_t c = res.second;
			uint64_t lb_new = carray[c] + res.first;

			// Right maximal
			uint64_t ones = bv_rank(lb_new+1);
			uint64_t next_node = numeric_limits<uint64_t>::max();
			if(ones % 2 == 1 or bv[lb_new] == 1)
			{
				next_node = (ones-1)/2;
			}

			if(c <= 1) // c == sentinal
			{
				graph[cur_node].pos_list.emplace_back(i+1);
				start_nodes.emplace_back(cur_node);
				cur_node = graph.size();
				graph.emplace_back(node());
				graph[cur_node].len = 1;
				// change now_seq
				graph[cur_node].color = ++ now_seq;
			}
			else if(next_node != numeric_limits<uint64_t>::max() or bv2[lb]) // Next node is right max or cur node is left maximal => split
			{
				if(next_node == numeric_limits<uint64_t>::max())
				{
					next_node = number_right_max_nodes + bv3_rank(lb_new);
				}
				graph[cur_node].pos_list.emplace_back(i+1);
				graph[next_node].adj_list.emplace_back(cur_node);
				graph[next_node].len = k;

				seq_colors.merge(graph[next_node].color, graph[cur_node].color);
				graph[next_node].color = graph[cur_node].color = seq_colors.find(graph[cur_node].color);

				cur_node = next_node;
			}
			else
			{
				++graph[cur_node].len;
			}
			lb = lb_new;
		}
		graph[cur_node].pos_list.emplace_back(1);
		start_nodes.emplace_back(cur_node);
		reverse(begin(start_nodes), end(start_nodes));
		// update colors
		for(auto &node: graph) node.color = seq_colors.find(node.color);
		// for(uint64_t i = 0; i < seqs; ++ i) seq_colors.find(i);
		reverse(begin(seq_colors.pa), end(seq_colors.pa));
		s.clear();
		for(uint64_t i = 0; i < seqs; ++ i) s.insert(seq_colors.rev_find(i));
		auto stop = high_resolution_clock::now();
		cerr << std::setw(10) << duration_cast<milliseconds>(stop-start).count() << "ms for creating graph" << endl;
	}
	return make_tuple(std::move(graph), std::move(start_nodes), std::move(seq_colors.pa), std::move(vector<uint64_t>(s.begin(), s.end())));
}

int main(int argc, char *argv[])
{
	ios::sync_with_stdio(false);
	cin.tie(nullptr);
	// Get parameters
	string inputfile, outputfile, tmp_dir = ".";
	uint64_t k = 0;
	byte_sa_algo_type sa_method = LIBDIVSUFSORT;
	bool use_cache = true;
	int t = 1;
	while (1)
	{
		//int this_option_optind = optind ? optind : 1;
		int option_index = 0, c;
		static struct option long_options[] =
		{
			{"inputfile",  required_argument, 0, 'i'},
			{"outputfile", required_argument, 0, 'o'},
			{"help",       no_argument,       0, 'h'},
			{"kmer",       required_argument, 0, 'k'},
			{"method",     required_argument, 0, 'm'},
			{"threads",    required_argument, 0, 't'},
			{"nocache",    required_argument, 0, 'n'},
			{"tmpdir",     required_argument, 0, 'd'},
			{"version",    no_argument,       0, 'v'},
			{0,            0,                 0,  0 }
		};
		c = getopt_long(argc, argv, "i:o:hnk:d:m:v", long_options, &option_index);
		if (c == -1) break;
		switch (c)
		{
			case 0:
				cerr << "option " << long_options[option_index].name << ' ';
				if (optarg) cerr << "with arg %s ";
				cerr << "not supported" << endl;
				break;
			case 'h': usage(argv[0]); break;
			case 'k':
				k = stoull(argv[optind - 1]);
				cerr << "Info: k = " << k << '\n';
				break;
			case 'i':
				inputfile = argv[optind - 1];
				cerr << "Input file = " << inputfile << '\n';
				break;
			case 'o':
				outputfile = argv[optind - 1];
				cerr << "Output file = " << outputfile << '\n';
				break;
			case 'm':
				// {SAIS, DIV, PARDIV}
				if(! strcmp("SAIS", argv[optind - 1]))
				{
					cerr << "SA construct method = SA-IS\n";
					sa_method = SE_SAIS;
				}
				break;
			case 't':
				t = atoi(argv[optind - 1]);
				setWorkers(t);
				cerr << "threads = " << getWorkers() << '\n';
				break;
			case 'n':
				use_cache = false;
				cerr << "Force generate BWT file" << '\n';
				break;
			case 'd':
				tmp_dir = argv[optind - 1];
				break;
			case '?':
				usage(argv[0]);
				break;
			case 'v':
				version(argv[0]);
				break;
		}
	}
	if(k == 0)
	{
		cerr << "Error: k is not defined. Check your arugments and try again." << endl;
		exit(EXIT_FAILURE);
	}
	else if(inputfile == "" || outputfile == "")
	{
		cerr << "Error: not determined input/output file. Check your arugments and try again." << endl;
		exit(EXIT_FAILURE);
	}

	// Create datastructures
	vector <uint64_t> seq_len;
	cache_config config(true, tmp_dir, "tmp");
	uint64_t errors = create_datastructures(config, inputfile, k, sa_method, t, seq_len, use_cache, outputfile);
	if(errors)
	{
		cerr << "Error: Error on create BWT. Program will exit." << endl;
		exit(EXIT_FAILURE);
	}

	// Read k-value
	auto [graph, start_nodes, seq_colors, cluster_colors] = create_cdbg_with_lf(config, k, seq_len.size());
#if 0
	cerr << "start_nodes\tseq_colors\tcluster_colors" << '\n';
	for(uint64_t i = 0; i < seq_len.size(); ++ i)
	{
		cerr << start_nodes[i] << "\t\t" << seq_colors[i];
		if(i < cluster_colors.size()) cerr << "\t\t" << cluster_colors[i];
		cerr << '\n';
	}
	cerr << flush;
#endif
	// Print graph
	auto start = high_resolution_clock::now();
	ofstream output(outputfile);
	print_graph(graph, start_nodes, seq_colors, cluster_colors, output, config, seq_len);
	// print_graph(graph, start_nodes, output, output_sum);
	auto stop = high_resolution_clock::now();
	cerr << std::setw(10) << duration_cast<milliseconds>(stop-start).count() << "ms for printing graph" << endl;

	// Delete files
	if(config.delete_files)
	{
		util::delete_all_files(config.file_map);
	}

	return errors;
}
