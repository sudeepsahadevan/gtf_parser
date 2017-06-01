#include <iostream>
#include <fstream>
#include <regex>
#include <stdio.h>
#include <getopt.h>
#include <unordered_set>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/zlib.hpp>

/*
  	parse_gtf.cpp
  	GTF/GFF file parser
  	Created on: May 29, 2017
  	Copyright (C) 2017 Sudeep Sahadevan

 	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>
*/

using namespace std;

class GTF_parser {
private:
	unordered_set<string> idset;
	bool remove_version;
//	regular expressions
	regex gz_re = regex("^.*\\.gz$",regex_constants::icase);
	regex version_re = regex("\\.\\d{1,}$");
	regex gtf_re = regex("^.*\\.gtf.*$",regex_constants::icase);
	regex gff_re = regex("^.*\\.gff.*$",regex_constants::icase);
//	functions
	void ids_fh(string id_file);
	void parse_ids(string ids);
	void gff_fh(string feature_file);
	void gff_fh(string feature_file, string feature);
	void filter_gff(string line,regex feat_re,regex attrib_re);
	void filter_gff_noversion(string line,regex feat_re,regex attrib_re);
	void filter_gff(string line,regex attrib_re);
	void filter_gff_noversion(string line,regex attrib_re);
	regex get_attrib_regex(string feature_file);

public:
	GTF_parser(string id_file,string feature_file,bool remove_version);
	GTF_parser(string id_file,string feature_file,string feature,bool remove_version);

};

GTF_parser::GTF_parser(string id_file, string feature_file, bool remove_version){
	this->remove_version = remove_version;
	ids_fh(id_file);
	gff_fh(feature_file);
}

GTF_parser::GTF_parser(string id_file, string feature_file, string feature, bool remove_version){
	this->remove_version = remove_version;
	ids_fh(id_file);
	gff_fh(feature_file,feature);
}

void GTF_parser::ids_fh(string id_file){
	/*
	 * stream .gz files line by line
	 * 	source: https://stackoverflow.com/questions/29574044/how-to-use-boostiostreamsmapped-file-source-with-a-gzipped-input-file
	 */
	cerr << "Parsing name/id list from file: "<<id_file<<endl;
	string ids;
	if(regex_match(id_file,gz_re)){
		ifstream infh(id_file, ios_base::in | ios_base::binary);
		boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
		inbuf.push(boost::iostreams::gzip_decompressor());
		inbuf.push(infh);
		istream instream(&inbuf);
		while(getline(instream, ids)) {
			parse_ids(ids);
		}
		infh.close();
	}else{
		ifstream infh;
		infh.open(id_file,fstream::in);
		if(infh.is_open()){
			while(getline(infh,ids)){
				parse_ids(ids);
			}
			infh.close();
		}
	}
	if(idset.size()==0){
		cerr << "ERROR! cannot parse feature names/ids from file "+id_file+". Check your input" << endl;
		cerr << "exiting..." <<  endl;
		exit(1);
	}else{
		cerr<< "Omit version: "<<boolalpha<<remove_version << endl;
		cerr << "Found " << idset.size()<< " unique names/ids" << endl;
	}
}

void GTF_parser::parse_ids(string ids){
	transform(ids.begin(), ids.end(), ids.begin(), ::toupper);
	idset.insert(ids);
	if(remove_version){
		idset.insert(regex_replace(ids, version_re,""));
	}
}

void GTF_parser::gff_fh(string feature_file){
	cerr << "Parsing feature file: "<<feature_file<<endl;
	regex attrib_re = get_attrib_regex(feature_file);
	string line;
	if(regex_match(feature_file,gz_re) && remove_version){
		ifstream infh(feature_file, ios_base::in | ios_base::binary);
		boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
		inbuf.push(boost::iostreams::gzip_decompressor());
		inbuf.push(infh);
		istream instream(&inbuf);
		while(getline(instream, line)) {
			filter_gff_noversion(line, attrib_re);
		}
		infh.close();
	}else if(regex_match(feature_file,gz_re) && (!remove_version)){
		ifstream infh(feature_file, ios_base::in | ios_base::binary);
		boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
		inbuf.push(boost::iostreams::gzip_decompressor());
		inbuf.push(infh);
		istream instream(&inbuf);
		while(getline(instream, line)) {
			filter_gff(line, attrib_re);
		}
		infh.close();
	}else if((!regex_match(feature_file,gz_re)) && remove_version){
		ifstream infh;
		infh.open(feature_file,fstream::in);
		if(infh.is_open()){
			while(getline(infh,line)){
				filter_gff_noversion(line, attrib_re);
			}
			infh.close();
		}
	}else if((!regex_match(feature_file,gz_re)) && (!remove_version)){
		ifstream infh;
		infh.open(feature_file,fstream::in);
		if(infh.is_open()){
			while(getline(infh,line)){
				filter_gff(line, attrib_re);
			}
			infh.close();
		}
	}
}

void GTF_parser::gff_fh(string feature_file, string feature){
	cerr << "Parsing feature file: "<<feature_file<<" for feature: "<<feature<<endl;
	regex attrib_re = get_attrib_regex(feature_file);
	regex feat_re = regex("^.*\\t"+feature+"\\t.*$");
	string line;
	if(regex_match(feature_file,gz_re) && remove_version){
		ifstream infh(feature_file, ios_base::in | ios_base::binary);
		boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
		inbuf.push(boost::iostreams::gzip_decompressor());
		inbuf.push(infh);
		istream instream(&inbuf);
		while(getline(instream, line)) {
			filter_gff_noversion(line, feat_re, attrib_re);
		}
		infh.close();
	}else if(regex_match(feature_file,gz_re) && (!remove_version)){
		ifstream infh(feature_file, ios_base::in | ios_base::binary);
		boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
		inbuf.push(boost::iostreams::gzip_decompressor());
		inbuf.push(infh);
		istream instream(&inbuf);
		while(getline(instream, line)) {
			filter_gff(line, feat_re, attrib_re);
		}
		infh.close();
	}else if((!regex_match(feature_file,gz_re)) && remove_version){
		ifstream infh;
		infh.open(feature_file,fstream::in);
		if(infh.is_open()){
			while(getline(infh,line)){
				filter_gff_noversion(line, feat_re, attrib_re);
			}
			infh.close();
		}
	}else if((!regex_match(feature_file,gz_re)) && (!remove_version)){
		ifstream infh;
		infh.open(feature_file,fstream::in);
		if(infh.is_open()){
			while(getline(infh,line)){
				filter_gff(line, feat_re, attrib_re);
			}
			infh.close();
		}
	}
}

regex GTF_parser::get_attrib_regex(string feature_file){
	regex attrib_re;
	if(regex_match(feature_file,gtf_re)){
		attrib_re = regex("(?:id|name)\\s{1,}\"([^\"]+)\"",regex_constants::icase);
	}else if(regex_match(feature_file,gff_re)){
		attrib_re = regex("(?:id|name)\\=([^=]+)(?:\\;|$)",regex_constants::icase);
	}else{
		cerr << "ERROR! cannot determine file type for file "+feature_file+". Check your input" << endl;
		cerr << "exiting..." <<  endl;
		exit(1);
	}
	return attrib_re;
}

void GTF_parser::filter_gff(string line, regex feat_re, regex attrib_re){
	if(regex_match(line,feat_re)){
		sregex_token_iterator swb(line.begin(),line.end(),attrib_re,1);
		sregex_token_iterator swe;
		while(swb!=swe){
			string attrib_id =*swb++;
			transform(attrib_id.begin(), attrib_id.end(), attrib_id.begin(), ::toupper);
			if(idset.find(attrib_id)!=idset.end()){
				cout<<line<<endl;
				break;
			}
		}
	}
}

void GTF_parser::filter_gff_noversion(string line, regex feat_re, regex attrib_re){
	if(regex_match(line,feat_re)){
		sregex_token_iterator swb(line.begin(),line.end(),attrib_re,1);
		sregex_token_iterator swe;
		while(swb!=swe){
			string attrib_id =*swb++;
			transform(attrib_id.begin(), attrib_id.end(), attrib_id.begin(), ::toupper);
			if(idset.find(attrib_id)!=idset.end()){
				cout<<line<<endl;
				break;
			}else if(idset.find(regex_replace(attrib_id,version_re,""))!=idset.end()){
				cout<<line<<endl;
				break;
			}
		}
	}
}

void GTF_parser::filter_gff(string line, regex attrib_re){
	sregex_token_iterator swb(line.begin(),line.end(),attrib_re,1);
	sregex_token_iterator swe;
	while(swb!=swe){
		string attrib_id =*swb++;
		transform(attrib_id.begin(), attrib_id.end(), attrib_id.begin(), ::toupper);
		if(idset.find(attrib_id)!=idset.end()){
			cout<<line<<endl;
			break;
		}
	}
}

void GTF_parser::filter_gff_noversion(string line, regex attrib_re){
	sregex_token_iterator swb(line.begin(),line.end(),attrib_re,1);
	sregex_token_iterator swe;
	while(swb!=swe){
		string attrib_id =*swb++;
		transform(attrib_id.begin(), attrib_id.end(), attrib_id.begin(), ::toupper);
		if(idset.find(attrib_id)!=idset.end()){
			cout<<line<<endl;
			break;
		}else if(idset.find(regex_replace(attrib_id,version_re,""))!=idset.end()){
			cout<<line<<endl;
			break;
		}
	}
}

void show_help(string prog){
	cout << "Usage: "<<prog<<" [-argument] <Val>"<<endl;
//	cout << "Copyright (C) 2017  Sudeep Sahadevan" << endl;
	cout << "Given a list of gene names/ids in a file and a gtf/gff file, filter the gtf/gff file for features annotated with input names/ids" << endl;
	cout << "--help, -h	Show this message and exit"<<endl;
	cout << "--version, -v	Show version information and exit"<<endl;
	cout << "--ids, -i	<path>	Path to file with gene names/ids. One entry per line, supports gzipped (.gz) files"<<endl;
	cout << "--gff, -g	<path>	Path to GTF/GFF formatted file. Supports gzipped (.gz) files"<<endl;
	cout << "--feature, -f	<string> [optional]	Feature to parse, from the 3rd column in GTF/GFF formatted files" << endl;
	cout << "--noversion, -n	[optional]	This flag omits version number from gene names/ids before filtering" << endl;
	cout << "example use: " << endl;
	cout << "	"+prog+" --ids /path/to/ids --gff /path/to/gtf --feature exon --noversion > gff_filtered" << endl;
}

int main(int argc, char *argv[]){
	/*
	 * getopt parse optional arguments:
	 * 	source : https://stackoverflow.com/questions/1052746/getopt-does-not-parse-optional-arguments-to-parameters
	 */
	int in_arg;
	int opt_ind;
	string ids_file;
	string feature_file;
	string feature;
	bool remove_version = false;
	regex prog_re("^.*\\/");
	static struct option long_options[] = {
			{"help",no_argument,NULL,'h'},
			{"version",no_argument,NULL,'v'},
			{"noversion",no_argument,NULL,'n'},
			{"ids",required_argument,NULL,'i'},
			{"gff",required_argument,NULL,'g'},
			{"feature",optional_argument,NULL,'f'},
			{NULL,0,NULL,0}
	};
	if(argc == 1){
			show_help(regex_replace(argv[0],prog_re,""));
			exit(1);
		}
	while((in_arg=getopt_long(argc,argv,"hvni:g:f::",long_options,&opt_ind))!=-1){
		const char *tmp_optarg = optarg;
		switch(in_arg){
		case 'h':
			show_help(regex_replace(argv[0],prog_re,""));
			exit(0);
			break;
		case 'v':
			cout << "version 0.0.1.1" << endl;
			exit(0);
			break;
		case 'n':
			remove_version = true;
			break;
		case 'i':
			ids_file = optarg;
			break;
		case 'g':
			feature_file = optarg;
			break;
		case 'f':
			if( !optarg && NULL != argv[optind] && '-' != argv[optind][0] && optind < argc ){
				 tmp_optarg = argv[optind++];
			}
			if(tmp_optarg){
				feature = tmp_optarg;
			}
			break;
		default:
			show_help(regex_replace(argv[0],prog_re,""));
			exit(1);
			break;
		}
	}
	if(ids_file.empty()){
		cout << "ERROR! required argument --ids|-i is missing" << endl;
		show_help(regex_replace(argv[0],prog_re,""));
		exit(1);
	}else if(feature_file.empty()){
		cout << "ERROR! required argument --gff|-g is missing" << endl;
		show_help(regex_replace(argv[0],prog_re,""));
		exit(1);
	}
	if(feature.empty()){
		GTF_parser gtfparser(ids_file, feature_file, remove_version);
	}else{
		GTF_parser gtfparser(ids_file, feature_file, feature, remove_version);
	}
	return 0;
}
