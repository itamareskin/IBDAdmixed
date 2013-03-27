#include <math.h>

int c_isfinite(double x);

typedef struct struct_state
{
	double prob_em[2];
	int out_trans_num;
	int in_trans_num;
} state;

typedef struct struct_transition
{
	double prob;
	int next_node;
} transition;

typedef struct struct_gen_map_entry
{
	int position;
	double recomb_rate;
	double genetic_dist;
} gen_map_entry;

state create_state(double allele_0_prob, double allele_1_prob, int out_trans_num, int in_trans_num);
transition create_transition(double prob, int next_node);
gen_map_entry create_gen_map_entry(int position, double recomb_rate, double genetic_dist);
char get_likely_allele(state s);
