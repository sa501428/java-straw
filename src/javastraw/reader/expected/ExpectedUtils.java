package javastraw.reader.expected;

import javastraw.reader.datastructures.ListOfDoubleArrays;

public class ExpectedUtils {
    public static ListOfDoubleArrays smooth(ListOfDoubleArrays expectedValues) {
        ListOfDoubleArrays smoothArray = new ListOfDoubleArrays(expectedValues.getLength());
        smoothArray.setDataAfterSmoothing(expectedValues);
        return smoothArray;
    }
}
