#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "ped_to_lamp.h"

using namespace std;

const char MIS='0';
const char HET='h';

int ped2lamp(char *ped_file, char *map_file, char *lamp_file, bool phased)
{
	string chrom, rsid, dist, pos, fam, id;
	int nr_snp;

	ifstream file_ped(ped_file);
	ifstream file_map(map_file);
	ofstream file_lamp(lamp_file);
	if(!file_ped || !file_map || !file_lamp) { cerr << "file could not be opened" << endl; return 0; }

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

	// print sequences

	for(int i=0;i<nr_samples;i++)
	{
	    if (phased) {
	        for(int a=0;a<2;a++) {
                for(int h=0;h<nr_snp;h++) {
                    if ((haplotype[a][i][h] != '0') && (haplotype[a][i][h] != '1') && (haplotype[a][i][h] != '2')) {
                        cerr << "alleles should be coded 0/1/2" << endl;
                        return 0;
                    }
                    int allele1 = haplotype[a][i][h] - '0';
                    if (allele1 == 0) {
                        file_lamp << '?';
                    } else {
                        file_lamp << allele1 - 1;
                    }
                }
                file_lamp << endl;
            }
	    } else {
            for(int h=0;h<nr_snp;h++)
            {
                if (((haplotype[0][i][h] != '0') &&
                     (haplotype[0][i][h] != '1') &&
                     (haplotype[0][i][h] != '2')) ||
                    ((haplotype[1][i][h] != '0') &&
                     (haplotype[1][i][h] != '1') &&
                     (haplotype[1][i][h] != '2'))) {
                    cerr << "alleles should be coded 0/1/2" << endl;
                    return 0;
                }
                int allele1 = haplotype[0][i][h] - '0';
                int allele2 = haplotype[1][i][h] - '0';
                if ((allele1 == 0) || (allele2 == 0)) {
                    file_lamp << '?';
                } else {
                    file_lamp << allele1 + allele2 - 2;
                }
            }
            file_lamp << endl;
        }
	}
	file_ped.close();

	return 1;
}