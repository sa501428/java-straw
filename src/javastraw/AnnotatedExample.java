package javastraw;

import javastraw.reader.Dataset;
import javastraw.reader.Matrix;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.block.Block;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.norm.NormalizationPicker;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;

import java.util.Iterator;
import java.util.List;

public class AnnotatedExample {
    public void run() {
        // do you want to cache portions of the file?
        // this uses more RAM, but if you want to repeatedly
        // query nearby regions, it can improve speed by a lot
        boolean useCache = false;
        String filename = "file.hic";

        // create a hic dataset object
        Dataset ds = HiCFileTools.extractDatasetForCLT(filename, false, useCache, false);

        // pick the normalization we would like
        // this line will check multiple possible norms
        // and pick whichever is available (in order of preference)
        NormalizationType norm = NormalizationPicker.getFirstValidNormInThisOrder(ds, new String[]{"KR", "SCALE", "VC", "VC_SQRT", "NONE"});
        System.out.println("Norm being used: " + norm.getLabel());

        // let's set our resolution
        int resolution = 5000;

        // let's grab the chromosomes
        Chromosome[] chromosomes = ds.getChromosomeHandler().getChromosomeArrayWithoutAllByAll();

        // now let's iterate on every chromosome (only intra-chromosomal regions for now)
        for (Chromosome chromosome : chromosomes) {
            Matrix matrix = ds.getMatrix(chromosome, chromosome);
            if (matrix == null) continue;
            MatrixZoomData zd = matrix.getZoomData(new HiCZoom(resolution));
            if (zd == null) continue;

            // zd is now a data structure that contains pointers to the data
            // *** Let's show 2 different ways to access data ***

            // OPTION 1
            // iterate on all the data for the whole chromosome in sparse format
            Iterator<ContactRecord> iterator = zd.getNormalizedIterator(norm);
            while (iterator.hasNext()) {
                ContactRecord record = iterator.next();
                // now do whatever you want with the contact record
                int binX = record.getBinX();
                int binY = record.getBinY();
                float counts = record.getCounts();

                // binX and binY are in BIN coordinates, not genome coordinates
                // to switch, we can just multiply by the resolution
                int genomeX = binX * resolution;
                int genomeY = binY * resolution;

                if (counts > 0) { // will skip NaNs

                    // do task

                    // the iterator only iterates above the diagonal
                    // to also fill in data below the diagonal, flip it
                    if (binX != binY) {
                        binX = record.getBinY();
                        binY = record.getBinX();
                        counts = record.getCounts();

                        // do task
                    }
                }
            }

            // OPTION 2
            // just grab sparse data for a specific region

            // choose your setting for when the diagonal is in the region
            boolean getDataUnderTheDiagonal = true;

            // our bounds will be binXStart, binYStart, binXEnd, binYEnd
            // these are in BIN coordinates, not genome coordinates
            int binXStart = 500, binYStart = 600, binXEnd = 1000, binYEnd = 1200;
            List<Block> blocks = zd.getNormalizedBlocksOverlapping(binXStart, binYStart, binXEnd, binYEnd, norm, getDataUnderTheDiagonal);
            for (Block b : blocks) {
                if (b != null) {
                    for (ContactRecord rec : b.getContactRecords()) {
                        if (rec.getCounts() > 0) { // will skip NaNs
                            // can choose to use the BIN coordinates
                            int binX = rec.getBinX();
                            int binY = rec.getBinY();

                            // you could choose to use relative coordinates for the box given
                            int relativeX = rec.getBinX() - binXStart;
                            int relativeY = rec.getBinY() - binYStart;

                            float counts = rec.getCounts();
                        }
                    }
                }
            }
        }

        // to iterate over the whole genome
        for (int i = 0; i < chromosomes.length; i++) {
            for (int j = i; i < chromosomes.length; i++) {
                Matrix matrix = ds.getMatrix(chromosomes[i], chromosomes[j]);
                if (matrix == null) continue;
                MatrixZoomData zd = matrix.getZoomData(new HiCZoom(resolution));
                if (zd == null) continue;

                if (i == j) {
                    // intra-chromosomal region
                } else {
                    // inter-chromosomal region
                }
            }
        }
    }
}
