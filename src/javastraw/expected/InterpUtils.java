package javastraw.expected;

import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;

import java.util.ArrayList;
import java.util.List;

public class InterpUtils {

    public static PolynomialSplineFunction cleanInterp(double[] x, double[] y) {
        List<Double> x2 = new ArrayList<>();
        List<Double> y2 = new ArrayList<>();
        addValidPoints(x, y, x2, y2);
        double[] x3 = convert(x2);
        double[] y3 = convert(y2);
        return new SplineInterpolator().interpolate(x3, y3);
    }

    public static void addValidPoints(double[] x, double[] y, List<Double> x2, List<Double> y2) {
        x2.add(x[0]);
        y2.add(y[0]);
        for (int k = 1; k < x.length; k++) {
            if (x[k] > x[k - 1]) {
                x2.add(x[k]);
                y2.add(y[k]);
            }
        }
    }

    public static double[] convert(List<Double> input) {
        double[] vec = new double[input.size()];
        int counter = 0;
        for (Double value : input) {
            vec[counter] = value;
            counter++;
        }
        return vec;
    }
}
