package javastraw.reader.norm;

import javastraw.reader.Dataset;
import javastraw.reader.type.NormalizationType;

import java.util.Arrays;
import java.util.Map;

public class NormalizationHandler {
    public static NormalizationType getFirstValidNormInThisOrder(Dataset ds, String[] norms) {
        Map<String, NormalizationType> normsForDataset = ds.getNormalizationTypesMap();

        for (String key : norms) {
            if (normsForDataset.containsKey(key)) {
                return normsForDataset.get(key);
            }
        }

        System.err.println("None of these normalizations were found: " + Arrays.toString(norms));
        System.exit(11);
        return null;
    }
}
