package javastraw.tools;

import javastraw.reader.Dataset;

public class Tester {

    public static void main(String[] args) {
        String newList = ("https://s3.us-central-1.wasabisys.com/aiden-encode-hic-mirror/mapq30/S5_IMR90.hic");
        Dataset ds = HiCFileTools.extractDatasetForCLT(newList, false, false, true);
        //System.out.println(ds.getStatistics());
        //System.out.println("\n\n\n\n");
        //System.out.println(ds.getGraphs());
    }

    /*
    public static void main33(String[] args) throws IOException {
        long s1 = System.nanoTime();
        IGVSeekableStreamFactory streamFactory = IGVSeekableStreamFactory.getInstance();
        SeekableStream stream = streamFactory.getStreamFor("https://s3.us-central-1.wasabisys.com/aiden-encode-hic-mirror/mapq30/S5_IMR90.hic");
        long s2 = System.nanoTime();
        System.out.println((s2-s1)*1e-9);
    }
    */

    public static void main22(String[] args) {
        String path = "/Users/mshamim/Desktop/chr10.hic";
        long start = System.currentTimeMillis();
        //for(int k = 0; k < 1; k++) {
        IterationUtils.processDataset(path, 500, 1); // 0 new
        //}
        long end = System.currentTimeMillis();
        long elapsed = (end - start) / 1000;
        System.out.println("Runtime " + elapsed + " (s)");
    }
}
