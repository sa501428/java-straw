package javastraw.tools;

import javastraw.reader.Dataset;
import javastraw.reader.datastructures.ListOfDoubleArrays;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationType;

public class Tester {

    public static void main(String[] args) {

        String[] files = new String[]{
                "/Volumes/AidenLabWD6/hicfiles/GM12878_insitu_Mnase_8.15.20.hic",
                "/Volumes/AidenLabWD6/hicfiles/GM12878_insitu_quadRE_8.15.20_30.hic",
                "/Volumes/AidenLabWD6/hicfiles/GM12878_intact_0.1Mnase_8.15.20.hic",
                "/Volumes/AidenLabWD6/hicfiles/GM12878_intact_quadRE_8.15.20_30.hic"
        };

        String[] stems = new String[]{
                "InSituMNase", "InSitu4RE", "IntactMNase", "Intact4RE"
        };
        String[] norms = new String[]{
                "SCALE", "SCALE", "SCALE", "SCALE"
        };

        for (int res : new int[]{1000, 100, 50, 10}) {
            for (int q = 0; q < files.length; q++) {
                Dataset ds = HiCFileTools.extractDatasetForCLT(files[q], true, false);
                NormalizationType norm = ds.getNormalizationHandler().getNormTypeFromString(norms[q]);

                long t2 = System.nanoTime();
                ListOfDoubleArrays d2 = ds.getCustomCorrectedExpectedValues(new HiCZoom(HiCZoom.HiCUnit.BP, res),
                        norm, 50000).getExpectedValuesWithNormalization(1);
                long t3 = System.nanoTime();

                MatrixTools.saveMatrixTextNumpy("/Users/mshamim/Desktop/radius_subcompartment_explore/batch3/50k_smooth_" + stems[q] + "_expected_" + res + "BP.npy", d2.getValues().get(0));
                System.out.println((t3 - t2) * 1e-9);
            }
        }


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
