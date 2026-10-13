#include "v10/format.h"
int main(int argc, char **argv) {
 if (argc != 2) return 2;
 hic10::Header h; h.genome="review";
 h.chromosomes={{"chrA",10000000,{}},{"chrB",10000000,{}}};
 h.resolutions[0]={{10000},{10000000}};
 hic10::Options o; o.includeAll=false; o.verifyDerived=false; o.derived={{10000000,10000}};
 hic10::Writer w(argv[1],h,o);
 w.matrix(0,1,[](uint8_t,uint32_t){hic10::Matrix m; m.scores=true; m.cells={{0,0,hic10::bits(0x1p60f)},{501,0,hic10::bits(-0x1p60f)},{0,1,hic10::bits(1.0f)}}; return m;});
 w.finish({});
}
