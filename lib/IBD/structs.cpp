#include <stdlib.h>
#include <stdio.h>
#include "structs.h"

int c_isfinite(double x)
{
	return isfinite(x);
}

state create_state(double allele_0_prob, double allele_1_prob, int out_trans_num, int in_trans_num)
{
	//state *s = malloc(sizeof(state));
	state s;
    s.prob_em[0] = (double)allele_0_prob;
    s.prob_em[1] = (double)allele_1_prob;
    s.out_trans_num = out_trans_num;
    s.in_trans_num = in_trans_num;
    if (allele_0_prob > 0.5)
    	s.likely_allele = 0;
    else
    	s.likely_allele = 1;
    return s;
}

transition create_transition(double prob, int next_node)
{
	transition t;
	t.prob = prob;
	t.next_node = next_node;
	return t;
}

gen_map_entry create_gen_map_entry(int position, double recomb_rate, double genetic_dist)
{
	gen_map_entry g;
	g.position = position;
	g.recomb_rate = recomb_rate;
	g.genetic_dist = genetic_dist;
	return g;
}


bool get_likely_allele(state s)
{
	if (s.prob_em[0] > 0.5)
		return 0;
	else
        return 1;
}

