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

#ifndef INCLUDED_CREATE_DATASTRUCTURES
#define INCLUDED_CREATE_DATASTRUCTURES

#include <iostream>
#include <fstream>
#include <limits>
#include <climits>
#include <sdsl/sdsl_concepts.hpp>
#include <sdsl/suffix_array_algorithm.hpp>

using namespace std;
using namespace sdsl;
using namespace std::chrono;

uint64_t create_datastructures(cache_config& config, const string& inputfile, const uint64_t& k, const byte_sa_algo_type& sa_method, const int &threads, vector<uint64_t> &sequences, bool use_cache, const string& outputfile)
{
	auto start = high_resolution_clock::now();
	uint64_t min_len = ULLONG_MAX;

	string method;
	// (1) Check, if the text is cached
	if ( !use_cache || !cache_file_exists(conf::KEY_TEXT, config) )
	{
		method = "sdsl-format";
		int_vector<8> text;
		uint64_t num_bytes = 1;
		load_vector_from_file(text, inputfile, num_bytes);
		bool skip = false;
		uint64_t target=0;
		uint64_t sequence_length = 0;
		for(uint64_t i=0; i<text.size(); ++i)
		{
			if(text[i] == '>')
			{
				skip = true;
				if(i>0)
				{
					text[target] = 1;
					++target;
					sequences.emplace_back(sequence_length);
					min_len = std::min(min_len, sequence_length);
					sequence_length = 0;
				}
			}
			else if(text[i] == '\n')
			{
				skip = false;
			}
			else if(!skip)
			{
				text[target] = text[i];
				++target;
				++sequence_length;
			}
		}
		sequences.emplace_back(sequence_length);
		text.resize(target);
		if(contains_no_zero_symbol(text, inputfile))
		{
			append_zero_symbol(text);
			store_to_cache(text, conf::KEY_TEXT, config);
		}
	}
	else
	{
		method = "register only";
		config.delete_files = false;
		int_vector<8> text;
		load_from_cache(text, conf::KEY_TEXT, config);
		for(uint64_t i=0, len=0; i<text.size(); ++i)
		{
			if(text[i] <= 1)
			{
				sequences.emplace_back(len);
				min_len = std::min(min_len, len);
				len = 0;
			}
			else
			{
				++len;
			}
		}
	}
	register_cache_file(conf::KEY_TEXT, config);
	auto stop = high_resolution_clock::now();
	cerr << std::setw(10) << duration_cast<milliseconds>(stop-start).count() << "ms for Text-Transformation (" << method << ")" << endl;

	// Check input
	if(!sequences.size() || ULLONG_MAX == min_len)
	{
		cerr << "Could not find a sequence." << endl;
		return 1;
	}

	if (min_len < k)
	{
		cerr << "k = " << k << " must smaller than sequence length, but in input file '" << inputfile << "' there is a sequence with length " << min_len << ".\n";
		cerr << "Will create empty file for remaining process." << endl;
		ofstream g(outputfile); g << endl; g.close();
		exit(EXIT_SUCCESS); // ok on empty file
	}

	// (2) Check, if the suffix array is cached
	start = high_resolution_clock::now();
	if ( !use_cache || !cache_file_exists(conf::KEY_SA, config) )
	{
		if(sa_method == LIBDIVSUFSORT)
		{
			method = "sdsl-divsufsort";
		}
		else if(sa_method == SE_SAIS)
		{
			method = "sdsl-sesais";
		}
		else
		{
			construct_config::threads = threads;
			method = "sdsl-divsufsort (" + std::to_string(construct_config::threads) + " threads)";
		}
		construct_config::byte_algo_sa = sa_method;
		construct_sa<8>(config);
	}
	else
	{
		method = "register only";
		config.delete_files = false;
	}
	register_cache_file(conf::KEY_SA, config);
	stop = high_resolution_clock::now();
	cerr << std::setw(10) << duration_cast<milliseconds>(stop-start).count() << "ms for SA construction (" << method << ")" << endl;

	// (3) Check, if bwt is cached
	start = high_resolution_clock::now();
	if ( !use_cache || !cache_file_exists(conf::KEY_BWT, config) )
	{
		method = "sdsl";
		construct_bwt<8>(config);
	}
	else
	{
		method = "register only";
		config.delete_files = false;
	}
	register_cache_file(conf::KEY_BWT, config);
	stop = high_resolution_clock::now();
	cerr << std::setw(10) << duration_cast<milliseconds>(stop-start).count() << "ms for BWT construction (" << method << ")" << endl;

	return 0;
}

#endif
