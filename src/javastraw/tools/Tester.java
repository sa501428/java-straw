package javastraw.tools;

public class Tester {

    public static void main(String[] args) {
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
