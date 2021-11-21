package javastraw.tools;

import javastraw.reader.Dataset;
import javastraw.reader.datastructures.ListOfDoubleArrays;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationType;

public class Tester {

    public static void main(String[] args) {


        Dataset ds = HiCFileTools.extractDatasetForCLT("/Users/mshamim/Desktop/hicfiles/gm12878_rh14_30.hic",
                true, false);
        NormalizationType norm = ds.getNormalizationHandler().getNormTypeFromString("KR");

        long t1 = System.nanoTime();
        ListOfDoubleArrays d1 = ds.getExpectedValues(new HiCZoom(HiCZoom.HiCUnit.BP, 5000),
                norm, false).getExpectedValuesWithNormalization(1);
        long t2 = System.nanoTime();
        ListOfDoubleArrays d2 = ds.getExpectedValues(new HiCZoom(HiCZoom.HiCUnit.BP, 5000),
                norm, true).getExpectedValuesWithNormalization(1);
        long t3 = System.nanoTime();
        MatrixTools.saveMatrixTextNumpy("/Users/mshamim/Desktop/exp1_5k.npy", d1.getValues().get(0));
        MatrixTools.saveMatrixTextNumpy("/Users/mshamim/Desktop/exp2_5k.npy", d2.getValues().get(0));

        System.out.println((t2 - t1) * 1e-9);
        System.out.println((t3 - t2) * 1e-9);


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
