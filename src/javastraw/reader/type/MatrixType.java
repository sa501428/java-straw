/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2011-2020 Broad Institute, Aiden Lab, Rice University, Baylor College of Medicine
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 *  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 *  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 *  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 *  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 *  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 *  THE SOFTWARE.
 */


package javastraw.reader.type;

public enum MatrixType {
    OBSERVED("Observed"),
    EXPECTED("Expected"),
    OE("Observed/Expected"),
    PROB("Observed Probability"),
    PROB2("Observed Probability Inf"),
    OEV2("Log(Observed/Expected)"),
    OEP1("(Observed+1)/(Expected+1)"),
    OEP1V2("Log((Observed+1)/(Expected+1))"),
    OME("Observed-Expected"),
    PEARSON("Observed Pearson"),
    LOG("Log(Observed+1)"),
    LOGEO("Log(Observed+1)/Log(Expected+1)"),
    METALOGEO("(Log(Observed+1)+1)/(Log(Expected+1)+1)"),
    EXPLOGEO("e^(Log(Observed+1)/Log(Expected+1))"),
    EXPLOGEOINV("e^(Log(Expected+1)/Log(Observed+1))"),
    NORM2("Observed Norm^2"),
    CONTROL("Control"),
    OECTRL("Control/ExpectedC"),
    PROBCTRL("Control Probability"),
    PROB2CTRL("Control Probability Inf"),
    OECTRLV2("Log(Control/ExpectedC)"),
    OECTRLP1("(Control+1)/(ExpectedC+1)"),
    OECTRLP1V2("Log((Control+1)/(ExpectedC+1))"),
    CME("Control-ExpectedC"),
    PEARSONCTRL("Control Pearson"),
    LOGC("Log(Control+1)"),
    LOGCEO("Log[Control+1]/Log[ExpectedC+1]"),
    EXPLOGCEO("e^(Log[Control+1]/Log[ExpectedC+1])"),
    NORM2CTRL("Control Norm^2"),
    RATIO("Observed/Control * (AvgC/AvgO)"),
    RATIOV2("Log(Observed/Control * (AvgC/AvgO))"),
    RATIOP1("(Observed+1)/(Control+1) * (AvgC+1)/(AvgO+1)"),
    RATIOP1V2("Log((Observed+1)/(Control+1) * (AvgC+1)/(AvgO+1))"),
    RATIO0("Observed/Control * (ExpC0/Exp0)"),
    RATIO0V2("Log(Observed/Control * (ExpC0/Exp0))"),
    RATIO0P1("(Observed+1)/(Control+1) * (ExpC0+1)/(Exp0+1)"),
    RATIO0P1V2("Log((Observed+1)/(Control+1) * (ExpC0+1)/(Exp0+1))"),
    VS("Observed vs Control"),
    OEVS("Observed/Expected vs Control/ExpectedC"),
    PROBVS("Probability Observed vs Control"),
    PROB2VS("Probability Inf Observed vs Control"),
    OEVSV2("Log(Observed/Expected) vs Log(Control/ExpectedC)"),
    OEVSP1("(Observed+1)/(Expected+1) vs (Control+1)/(ExpectedC+1)"),
    OEVSP1V2("Log((Observed+1)/(Expected+1)) vs Log((Control+1)/(ExpectedC+1))"),
    OERATIO("(Observed/Expected) / (Control/ExpectedC)"),
    OERATIOV2("Log((Observed/Expected) / (Control/ExpectedC))"),
    OERATIOP1("((Observed+1)/(Expected+1)) / ((Control+1)/(ExpectedC+1))"),
    OERATIOP1V2("Log(((Observed+1)/(Expected+1)) / ((Control+1)/(ExpectedC+1)))"),
    OERATIOMINUS("(Observed/Expected) - (Control/ExpectedC)"),
    OERATIOMINUSP1("(Observed+1)/(Expected+1) - (Control+1)/(ExpectedC+1)"),
    OCMEVS("Observed-Expected vs Control-Expected"),
    PEARSONVS("Observed Pearson vs Control Pearson"),
    LOGVS("Log(Observed/AvgO+1) vs Log(Control/AvgC+1)"),
    LOGEOVS("Log(Observed+1)/Log(Expected+1) vs Log(Control+1)/Log(ExpectedC+1)"),
    LOGRATIO("Log(Observed/AvgO+1)/Log(Control/AvgC+1)"),
    LOGRATIOV2("Log(Log(Observed/AvgO+1)/Log(Control/AvgC+1))"),
    LOGEORATIO("(Log(Observed+1)/Log(Expected+1)) / (Log(Control+1)/Log(ExpectedC+1))"),
    LOGEORATIOV2("Log((Log(Observed+1)/Log(Expected+1)) / (Log(Control+1)/Log(ExpectedC+1)))"),
    DIFF("Observed-Control"),
    NORM("Norm"),
    EIGENVECTOR("Eigenvector"),
    NORM2OBSVSCTRL("Observed Norm^2 vs Control Norm^2");

    public static final MatrixType[] enabledMatrixTypesWithControl = new MatrixType[]{
            OBSERVED, CONTROL, VS,
            EXPECTED,
            OE, OECTRL, OEVS,
            OEV2, OECTRLV2, OEVSV2,
            PROB, PROBCTRL, PROBVS,
            PROB2, PROB2CTRL, PROB2VS,
            RATIO, RATIOV2,
            PEARSON, PEARSONCTRL, PEARSONVS,
            LOG, LOGC, LOGVS,
            RATIO0, RATIO0V2,
            OERATIO, OERATIOV2,
            EXPLOGEO, EXPLOGCEO
    };

    public static final MatrixType[] enabledMatrixTypesNoControl =
            new MatrixType[]{OBSERVED, EXPECTED, OE, OEV2, PROB, PROB2, PEARSON, LOG, EXPLOGEO};

    private final String value;

    MatrixType(String value) {
        this.value = value;
    }

    public static MatrixType enumValueFromString(String text) {
        if (text != null) {
            for (MatrixType matrix : MatrixType.values()) {
                if (text.equalsIgnoreCase(matrix.value)) {
                    return matrix;
                }
            }
            if (text.equalsIgnoreCase("oe")) {
                return OE;
            }
        }
        return null;
    }

    /**
     * @param option
     * @return true is the option is generally available all maps or resolutions
     */
    public static boolean isObservedOrControl(MatrixType option) {
        return option == OBSERVED || option == CONTROL;
    }

    /**
     * @param option
     * @return true is the option can be manipulated by the color range slider
     */
    public static boolean isColorScaleType(MatrixType option) {
        return !isPearsonType(option);
    }

    /**
     * @param option
     * @return true if the option should allowed in genome-wide view
     */
    public static boolean isValidGenomeWideOption(MatrixType option) {
        return !(option.toString().toLowerCase().contains("expected") ||
                option.toString().toLowerCase().contains("probability"));
    }

    /**
     * @param option
     * @return true if the option requires control map but not expected vector
     */
    public static boolean isSimpleControlType(MatrixType option) {
        return option.toString().toLowerCase().contains("control") &&
                !option.toString().toLowerCase().contains("expected");
    }


    /**
     * @param option
     * @return true if the option involves comparison/divis (but not pearsons)
     */
    public static boolean isOEColorScaleType(MatrixType option) {
        return option == OEV2 || option == OEP1V2 || option == RATIOV2 || option == RATIOP1V2 || option == RATIO0V2
                || option == RATIO0P1V2 || option == OECTRLV2 || option == OECTRLP1V2
                || option == OEVSV2 || option == OEVSP1V2 || option == OERATIOV2 || option == OERATIOP1V2
                || option == LOGRATIOV2 || option == LOGEORATIOV2
                || isSubtactType(option);
    }

    public static boolean isSubtactType(MatrixType option) {
        return option.toString().toLowerCase().contains("-");
    }

    /**
     * @param option
     * @return true if the option only works for intrachromosomal, not interchromosomal (genomewide may still be allowed)
     */
    public static boolean isOnlyIntrachromosomalType(MatrixType option) {
        return isPearsonType(option) || isVSTypeDisplay(option) ||
                option.toString().toLowerCase().contains("probability");
    }

    /**
     * @param option
     * @return true if the option requires the expected vector
     */
    public static boolean isExpectedValueType(MatrixType option) {
        return option.toString().toLowerCase().contains("expected") ||
                option.toString().toLowerCase().contains("probability");
    }

    /**
     * @param option
     * @return true if the option uses pearson's
     */
    public static boolean isPearsonType(MatrixType option) {
        return option.toString().toLowerCase().contains("pearson");
    }

    /**
     * @param option
     * @return true if the option is dumped (clt) as a vector
     */
    public static boolean isDumpVectorType(MatrixType option) {
        return option == NORM || option == EXPECTED;
    }

    /**
     * @param option
     * @return true if the option is dumped (clt) as a matrix
     */
    public static boolean isDumpMatrixType(MatrixType option) {
        return option == OE || option == OBSERVED;
    }

    public static boolean isVSTypeDisplay(MatrixType option) {
        return option.toString().toLowerCase().contains("vs");
    }

    public static boolean isControlPearsonType(MatrixType option) {
        return option.toString().toLowerCase().contains("control") && option.toString().toLowerCase().contains("pearson");
    }

    public String toString() {
        return value;
    }
}
