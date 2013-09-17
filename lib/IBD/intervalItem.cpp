#include "intervalItem.h"

interval_item create_interval_item(int start, int end, float value, float length)
{
	//state *s = malloc(sizeof(state));
	interval_item s;
    s.start = (int)start;
    s.end = (int)end;
    s.value = (float)value;
    s.length = (float)length;
    return s;
}
