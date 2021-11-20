package javastraw.reader.expected;

import java.util.List;

public class QuickMedian {

    public static float fastMedian(List<Double> list) {
        int size = list.size();
        if (size < 1) {
            return Float.NaN;
        } else if (size == 1) {
            return list.get(0).floatValue();
        }

        float[] arr = new float[size];
        for (int k = 0; k < arr.length; k++) {
            arr[k] = list.get(k).floatValue();
        }

        return fastMedian(arr);
    }

    public static float fastMedian(float[] arr) {
        int len = arr.length;
        if (len % 2 == 1) {
            return kSelection(arr, 0, len - 1, len / 2);
        } else {
            float a = kSelection(arr, 0, len - 1, len / 2);
            float b = kSelection(arr, 0, len - 1, len / 2 - 1);
            return (a + b) / 2;
        }
    }

    public static float kSelection(float[] arr, int low, int high, int k) {
        int localLow = low;
        int localHigh = high;

        int partitionSortingValue = partition(arr, localLow, localHigh);
        while (partitionSortingValue != k) {
            if (partitionSortingValue < k) {
                localLow = partitionSortingValue + 1;
            } else {
                localHigh = partitionSortingValue - 1;
            }
            partitionSortingValue = partition(arr, localLow, localHigh);
        }
        return arr[partitionSortingValue];
    }

    static int partition(float[] arr, int low, int high) {
        float pivot = arr[high];
        int z = (low - 1);
        for (int j = low; j < high; j++) {
            if (arr[j] < pivot) {
                z++;
                float temp = arr[z];
                arr[z] = arr[j];
                arr[j] = temp;
            }
        }
        float temp = arr[z + 1];
        arr[z + 1] = arr[high];
        arr[high] = temp;
        return z + 1;
    }
}
