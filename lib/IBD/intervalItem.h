typedef struct struct_interval_item
{
	int start;
	int end;
	float value;
} interval_item;

interval_item create_interval_item(int start, int end, float value);
