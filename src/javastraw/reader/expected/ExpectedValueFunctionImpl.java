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

package javastraw.reader.expected;

import javastraw.reader.DatasetReader;
import javastraw.reader.datastructures.ListOfDoubleArrays;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationType;

import java.io.IOException;
import java.util.Map;

/**
 * Utility holder for Density calculation, for O/E maps.
 */
public class ExpectedValueFunctionImpl implements ExpectedValueFunction {

	private final DatasetReader reader;
	private final int binSize;
	private final NormalizationType type;
	private final HiCZoom.HiCUnit unit;
	private final Map<Integer, Double> normFactors;
	private final long nValues;
	private final long filePosition;
	private final int streamSize = 500000;
	private int version = -1;
	private ListOfDoubleArrays expectedValues;
	private long streamBound1 = 0;
	private long streamBound2 = 0;

	public ExpectedValueFunctionImpl(NormalizationType type, HiCZoom.HiCUnit unit, int binSize,
									 ListOfDoubleArrays expectedValues, Map<Integer, Double> normFactors) {
		this.type = type;
		this.unit = unit;
		this.binSize = binSize;
		this.normFactors = normFactors;
		this.expectedValues = expectedValues;
		this.nValues = expectedValues.getLength();
		this.filePosition = 0;
		this.reader = null;
	}

	public ExpectedValueFunctionImpl(NormalizationType type, HiCZoom.HiCUnit unit, int binSize, long nValues,
									 long filePosition, Map<Integer, Double> normFactors, DatasetReader reader) {
		this.type = type;
		this.unit = unit;
		this.binSize = binSize;
		this.normFactors = normFactors;
		this.expectedValues = null;
		this.nValues = nValues;
		this.filePosition = filePosition;
		this.reader = reader;
		this.version = reader.getVersion();
	}

	// This is exposed for testing, should not use directly
	public Map<Integer, Double> getNormFactors() {
		return normFactors;
	}

	/**
	 * Expected value vector.  No chromosome normalization
	 *
	 * @return Genome-wide expected value vector
	 */
	@Override
	public ListOfDoubleArrays getExpectedValuesNoNormalization() {
		streamExpectedVectorFromFileIfNeeded();
		return expectedValues;
	}

	private void streamExpectedVectorFromFileIfNeeded() {
		if (expectedValues == null && (streamBound1 == 0 && streamBound2 == 0)) {
			try {
				expectedValues = reader.readExpectedVectorPart(filePosition, nValues);
				streamBound1 = 0;
				streamBound2 = 0;
			} catch (IOException e) {
				System.err.println("Error reading expected vector");
				e.printStackTrace();
			}
		}
	}

	/**
	 * Gets the expected value, distance and coverage normalized, chromosome-length normalized
	 *
	 * @param chrIdx   Chromosome index
	 * @param distance Distance from diagonal in bins
	 * @return Expected value, distance and coverage normalized
	 */
	@Override
	public double getExpectedValue(int chrIdx, long distance) {
		double normFactor = 1.0;
		long streamIndex;
		if (normFactors != null && normFactors.containsKey(chrIdx)) {
			normFactor = normFactors.get(chrIdx);
		}
		if (expectedValues == null) {
			streamPortionOfExpectedVector(distance);
		}

		if (streamBound1 == 0 && streamBound2 == 0) {
			if (expectedValues.getLength() > 0) {
				if (distance >= expectedValues.getLength()) {
					return expectedValues.getLastValue() / normFactor;
				} else {
					return expectedValues.get(distance) / normFactor;
				}
			} else {
				System.err.println("Expected values array is empty");
				return -1;
			}
		} else if (distance >= streamBound1 && distance < streamBound2) {
			streamIndex = distance - streamBound1;
			return expectedValues.get(streamIndex) / normFactor;
		} else {
			streamPortionOfExpectedVector(distance);
			streamIndex = distance - streamBound1;
			return expectedValues.get(streamIndex) / normFactor;
		}
	}

	private void streamPortionOfExpectedVector(long distance) {
		streamBound1 = Math.max(0, distance - (streamSize / 2));
		streamBound2 = streamBound1 + streamSize;
		long position = version > 8 ? filePosition + (streamBound1 * 4) : filePosition + (streamBound1 * 8);
		try {
			expectedValues = reader.readExpectedVectorPart(position, streamSize);
		} catch (IOException e) {
			System.err.println("Error reading expected vector");
			e.printStackTrace();
		}
	}

	@Override
	public ListOfDoubleArrays getExpectedValuesWithNormalization(int chrIdx) {
		streamExpectedVectorFromFileIfNeeded();
		double normFactor = 1.0;
		if (normFactors != null && normFactors.containsKey(chrIdx)) {
			normFactor = normFactors.get(chrIdx);
		}

		if (expectedValues.getLength() > 0) {
			ListOfDoubleArrays normedExpectedValues = expectedValues.deepClone();
			normedExpectedValues.multiplyEverythingBy(1.0 / normFactor);
			return normedExpectedValues;
		} else {
			System.err.println("Expected values array is empty");
			return null;
		}
	}

	@Override
	public ExpectedValueFunction getCorrectedVersion(int window) {
		streamExpectedVectorFromFileIfNeeded();
		ListOfDoubleArrays smoothVector = expectedValues.deepClone();
		smoothVector.doRollingMedian(window / binSize);
		return new ExpectedValueFunctionImpl(type, unit, binSize, smoothVector, normFactors);
	}

	@Override
	public long getLength() {
		return nValues;
	}

	@Override
	public NormalizationType getNormalizationType() {
		return type;
	}

	@Override
	public HiCZoom.HiCUnit getUnit() {
		return unit;
	}

	@Override
	public int getBinSize() {
		return binSize;
	}
}
