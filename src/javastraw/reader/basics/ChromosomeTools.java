package javastraw.reader.basics;

import javastraw.StrawGlobals;
import javastraw.reader.basics.chrom.sizes.ChromosomeSizes;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Pattern;

public class ChromosomeTools {
    /**
     * Load the list of chromosomes based on given genome id or file
     *
     * @param idOrFile string
     * @return list of chromosomes
     */
    public static ChromosomeHandler loadChromosomes(String idOrFile) {

        InputStream is = null;

        try {
            // Note: to get this to work, had to edit Intellij settings
            // so that "?*.sizes" are considered sources to be copied to class path
            is = ChromosomeSizes.class.getResourceAsStream(idOrFile + ".chrom.sizes");

            if (is == null) {
                // Not an ID,  see if its a file
                File file = new File(idOrFile);

                try {
                    if (file.exists()) {
                        is = new FileInputStream(file);
                    } else {
                        System.err.println("Could not find chromosome sizes file for: " + idOrFile);
                        System.exit(35);
                    }
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }

            List<Chromosome> chromosomes = new ArrayList<>();
            chromosomes.add(0, null);   // Index 0 reserved for "whole genome" pseudo-chromosome

            Pattern pattern = Pattern.compile("\\s+");
            BufferedReader reader = new BufferedReader(new InputStreamReader(is), StrawGlobals.bufferSize);
            String nextLine;
            int idx = 1;

            try {
                while ((nextLine = reader.readLine()) != null) {
                    String[] tokens = pattern.split(nextLine);
                    if (tokens.length == 2) {
                        String name = tokens[0];
                        int length = Integer.parseInt(tokens[1]);
                        chromosomes.add(idx, new Chromosome(idx, name, length));
                        idx++;
                    } else {
                        System.out.println("Skipping " + nextLine);
                    }
                }
            } catch (Exception e) {
                e.printStackTrace();
            }

            // "pseudo-chromosome" All taken care of by by chromosome handler
            return new ChromosomeHandler(chromosomes, idOrFile, false);
        } finally {
            if (is != null) {
                try {
                    is.close();
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        }
    }
}
