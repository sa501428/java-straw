package javastraw.tools;

import javastraw.reader.Dataset;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.datastructures.ListOfDoubleArrays;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationType;

public class Tester {

    public static void main(String[] args) {


        Dataset ds = HiCFileTools.extractDatasetForCLT("/Users/mshamim/Desktop/hicfiles/gm12878_rh14_30.hic",
                true, false);
        ChromosomeHandler chromosomeHandler = ds.getChromosomeHandler();
        NormalizationType norm = ds.getNormalizationHandler().getNormTypeFromString("KR");

        ListOfDoubleArrays d1 = ds.getExpectedValues(new HiCZoom(HiCZoom.HiCUnit.BP, 50000),
                norm, false).getExpectedValuesWithNormalization(1);
        ListOfDoubleArrays d2 = ds.getExpectedValues(new HiCZoom(HiCZoom.HiCUnit.BP, 50000),
                norm, true).getExpectedValuesWithNormalization(1);

        MatrixTools.saveMatrixTextNumpy("/Users/mshamim/Desktop/exp1.npy", d1.getValues().get(0));
        MatrixTools.saveMatrixTextNumpy("/Users/mshamim/Desktop/exp2.npy", d2.getValues().get(0));


        /*
        float[][] temp = new float[5][6];
        String outpath = "/Users/mshamim/Desktop/temp.npy";
        MatrixTools.saveMatrixTextNumpy(outpath, temp);
        System.exit(0);

        /*
        Dataset ds = HiCFileTools.extractDatasetForCLT("/Users/mshamim/Desktop/hicfiles/gm12878_rh14_30.hic",
                true, false);
        ChromosomeHandler chromosomeHandler = ds.getChromosomeHandler();
        List<String> givenChromosomes = new ArrayList<>();
        givenChromosomes.add("1");
        givenChromosomes.add("14");
        givenChromosomes.add("19");

        for(String s : givenChromosomes){
            System.out.println(s);
        }

        chromosomeHandler = HiCFileTools.stringToChromosomes(givenChromosomes, chromosomeHandler);
        */
    }
}
