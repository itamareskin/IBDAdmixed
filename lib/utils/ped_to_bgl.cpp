#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "ped_to_bgl.h"

using namespace std;

const char MIS='0';
const char HET='h';

int ped2bgl(char *ped_file, char *map_file, char *bgl_file, char *markers_file)
{
	string chrom, rsid, dist, pos, fam, id;
	int nr_snp;

	ifstream file_ped(ped_file);
	ifstream file_map(map_file);
	ofstream file_bgl(bgl_file);
	ofstream file_markers(markers_file);
	if(!file_ped || !file_map || !file_bgl || !markers_file) { cerr << "file could not be opened" << endl; return 0; }

	// load marker names
	vector<string> rsids;
	vector<string> positions;
	while(!file_map.eof())
	{
		chrom = rsid = dist = pos = "";
		file_map >> chrom >> rsid >> dist >> pos;
		if(!file_ped.good() || chrom == "" || rsid == "" || dist == "" || pos == "") continue;
		rsids.push_back(rsid);
		positions.push_back(pos);
	}
	nr_snp = rsids.size();
	cerr << nr_snp << " SNPs" << endl;

	// load individuals
	vector<string> sample;
	vector< vector<char> > haplotype[2];
	int ctr = 0;
	string discard, line;
	while(!file_ped.eof())
	{
		fam = id = discard = "";
		file_ped >> fam >> id >> discard >> discard >> discard >> discard;
		if(fam == "" || id == "" || discard == "") continue;
		getline(file_ped,line);
		sample.push_back(id);
		vector<char> seq[2];
		for(int i=0;i<nr_snp;i++)
		{
			seq[0].push_back(line.at(i*4+1));
			seq[1].push_back(line.at(i*4+3));
		}
		haplotype[0].push_back(seq[0]);
		haplotype[1].push_back(seq[1]);
		ctr++;
	}
	int nr_samples = sample.size();
	cerr << nr_samples << " samples" << endl;
	file_map.close();
	file_ped.close();

	// print header
	file_bgl << "# sampleID";
	for(int i=0;i<nr_samples;i++)
	{
		file_bgl << '\t' << sample[i] << '\t' << sample[i];
	}
	// print sequences
	for(int h=0;h<nr_snp;h++)
	{
		file_bgl << endl << "M " << rsids[h];
		for(int i=0;i<nr_samples;i++)
		{
			file_bgl << '\t' << haplotype[0][i][h] << '\t' << haplotype[1][i][h];
		}
	}
	file_bgl.close();

	// print markers file
	for(int h=0;h<nr_snp;h++)
	{
		file_markers << rsids[h] << '\t' << positions[h] << '\t' << "1" << '\t' << "2" << endl;
	}
	file_markers.close();

	return 1;
}