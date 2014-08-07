#ifndef DEFINES_GERMLINE_0001
#define DEFINES_GERMLINE_0001

extern float MIN_MATCH_LEN;
extern int MARKER_SET_SIZE;
extern bool PRINT_MATCH_HAPS;
extern bool ROI;
extern bool HAP_EXT;
extern bool WIN_EXT;
extern bool ALLOW_HOM;
extern bool HOM_ONLY;
extern bool HAPLOID;
extern bool SILENT;
extern bool DEBUG;
extern bool BINARY_OUT;
extern int MAX_ERR_HOM;
extern int MAX_ERR_HET;

int dummy(int x);
int germline_main(int argc, char* argv[]);

#endif //DEFINES_GERMLINE_0001