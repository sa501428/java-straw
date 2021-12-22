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


package javastraw.feature2D;

import java.awt.*;
import java.util.*;


/**
 * reflection only used for plotting, should not be used by CLTs
 */
public class Feature2D implements Comparable<Feature2D> {

    public static final String genericHeader = "#chr1\tx1\tx2\tchr2\ty1\ty2\tname\tscore\tstrand1\tstrand2\tcolor";
    private static final String BEDPE_SPACER = "\t.\t.\t.\t.";
    public static int tolerance = 0;
    protected final FeatureType featureType;
    protected final Map<String, String> attributes;
    protected final String chr1;
    protected final String chr2;
    protected final long start1;
    protected final long start2;
    protected long end1;
    protected long end2;
    protected Feature2D reflection = null;
    protected Color color;

    public Feature2D(FeatureType featureType, String chr1, long start1, long end1, String chr2, long start2, long end2, Color c,
                     Map<String, String> attributes) {
        this.featureType = featureType;
        this.chr1 = chr1;
        this.start1 = start1;
        this.end1 = end1;
        this.chr2 = chr2;
        this.start2 = start2;
        this.end2 = end2;
        this.color = (c == null ? Color.black : c);
        this.attributes = attributes;
    }

    public static String getDefaultOutputFileHeader() {
        return genericHeader;
    }

    public FeatureType getFeatureType() {
        return this.featureType;
    }

    public String getFeatureName() {
        switch (featureType) {
            case PEAK:
                return "Peak";
            case DOMAIN:
                return "Contact Domain";
            case GENERIC:
            case NONE:
            default:
                return "Feature";
        }
    }

    public String getChr1() {
        return chr1;
    }

    public String getChr2() {
        return chr2;
    }

    public long getStart1() {
        return start1;
    }

    public long getStart2() {
        return start2;
    }

    public long getEnd1() {
        return end1;
    }

    public void setEnd1(int end1) {
        this.end1 = end1;
        if (reflection != null)
            reflection.end2 = end1;
    }

    public long getEnd2() {
        return end2;
    }

    public void setEnd2(int end2) {
        this.end2 = end2;
        if (reflection != null)
            reflection.end1 = end2;
    }

    public long getWidth1() {
        return end1 - start1;
    }

    public long getWidth2() {
        return end2 - start2;
    }

    public long getMidPt1() {
        return midPoint(start1, end1);
    }

    public long getMidPt2() {
        return midPoint(start2, end2);
    }

    public long midPoint(long start, long end) {
        return (long) (start + (end - start) / 2.0);
    }

    public Color getColor() {
        return color;
    }


    public void setColor(Color color) {
        this.color = color;
        if (reflection != null)
            reflection.color = color;
    }

    public String getOutputFileHeader() {
        StringBuilder output = new StringBuilder(getDefaultOutputFileHeader());

        ArrayList<String> keys = new ArrayList<>(attributes.keySet());
        Collections.sort(keys);

        for (String key : keys) {
            output.append("\t").append(key);
        }

        return output.toString();
    }

    public String simpleString() {
        return chr1 + "\t" + start1 + "\t" + end1 + "\t" + chr2 + "\t" + start2 + "\t" + end2;
    }

    public String justColorString() {
        return "\t" + color.getRed() + "," + color.getGreen() + "," + color.getBlue();
    }

    public String simpleStringWithColor() {
        return simpleString() + BEDPE_SPACER + justColorString();
    }

    @Override
    public String toString() {
        StringBuilder output = new StringBuilder(simpleStringWithColor());

        ArrayList<String> keys = new ArrayList<>(attributes.keySet());
        Collections.sort(keys);
        for (String key : keys) {
            output.append("\t").append(attributes.get(key));
        }

        return output.toString();
    }

    public ArrayList<String> getAttributeKeys() {
        ArrayList<String> keys = new ArrayList<>(attributes.keySet());
        Collections.sort(keys);
        return keys;
    }

    public boolean hasAttributeKey(String key) {
        return attributes.containsKey(key);
    }

    public String getAttribute(String key) {
        return attributes.get(key);
    }

    public Map<String, String> getAttributes() {
        return attributes;
    }

    public void setAttribute(String key, String newVal) {
        attributes.put(key, newVal);
        // attribute directly shared between reflections
        if (reflection != null)
            reflection.attributes.put(key, newVal);
    }

    public float getFloatAttribute(String key) {
        return Float.parseFloat(attributes.get(key));
    }

    public void addIntAttribute(String key, int value) {
        attributes.put(key, "" + value);
    }

    public void addFloatAttribute(String key, Float value) {
        attributes.put(key, "" + value);
    }

    public void addStringAttribute(String key, String value) {
        attributes.put(key, value);
    }

    /**
     * @param otherFeature
     * @return
     */
    public boolean overlapsWith(Feature2D otherFeature) {

        float window1 = (otherFeature.getEnd1() - otherFeature.getStart1()) / 2;
        float window2 = (otherFeature.getEnd2() - otherFeature.getStart2()) / 2;

        long midOther1 = otherFeature.getMidPt1();
        long midOther2 = otherFeature.getMidPt2();

        return midOther1 >= (this.start1 - window1) && midOther1 <= (this.end1 + window1) && midOther2 >= (this.start2 - window2) && midOther2 <= (this.end2 + window2);
    }

    @Override
    public int compareTo(Feature2D o) {
        // highest observed point ordering needed for hiccups sorting
        long[] comparisons = new long[]{chr1.compareTo(o.chr1), chr2.compareTo(o.chr2), start1 - o.start1,
                start2 - o.start2, end1 - o.end1, end2 - o.end2};
        for (long i : comparisons) {
            if (i != 0)
                if (i > 0) {
                    return 1;
                } else {
                    return -1;
                }
        }
        return 0;
    }

    public boolean isOnDiagonal() {
        return chr1.equals(chr2) && start1 == start2 && end1 == end2;
    }

    public Feature2D reflectionAcrossDiagonal() {
        if (reflection == null) {
            reflection = new Feature2D(featureType, chr2, start2, end2, chr1, start1, end1, color, attributes);
            reflection.reflection = this;
        }
        return reflection;
    }

    public boolean isInLowerLeft() {
        return chr1.equals(chr2) && start2 > start1;
    }

    public boolean isInUpperRight() {
        return chr1.equals(chr2) && start2 < start1;
    }

    public boolean doesNotContainAttributeKey(String attribute) {
        return !attributes.containsKey(attribute);
    }

    public boolean containsAttributeValue(String attribute) {
        return attributes.containsValue(attribute);
    }

    public String getLocationKey() {
        return start1 + "_" + start2;
    }

    @Override
    public boolean equals(Object obj) {
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }
        if (this == obj) {
            return true;
        }

        final Feature2D other = (Feature2D) obj;
        if (chr1.equals(other.chr1)) {
            if (chr2.equals(other.chr2)) {
                if (Math.abs(start1 - other.start1) <= tolerance) {
                    if (Math.abs(start2 - other.start2) <= tolerance) {
                        if (Math.abs(end1 - other.end1) <= tolerance) {
                            return Math.abs(end2 - other.end2) <= tolerance;
                        }
                    }
                }
            }
        }

        return false;
    }

    @Override
    public int hashCode() {
        return Objects.hash(chr1, end1, start1, chr2, end2, start2);
    }

    public void clearAttributes() {
        attributes.clear();
    }

    public Feature2D deepCopy() {
        Map<String, String> attrClone = new HashMap<>();
        for (String key : attributes.keySet()) {
            attrClone.put(key, attributes.get(key));
        }
        return new Feature2D(featureType, chr1, start1, end1, chr2, start2, end2, color, attrClone);
    }

    public boolean containsPoint(float x, float y) {
        return start1 <= x && x <= end1 && start2 <= y && y <= end2;
    }

    public enum FeatureType {
        NONE, PEAK, DOMAIN, GENERIC, SCAFFOLD, SUPERSCAFFOLD, SELECTED_GROUP
    }
}
