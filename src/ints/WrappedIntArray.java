/*
 * Copyright 2020 Brian L. Browning
 *
 * This file is part of the ibd-ends program.
 *
 * Licensed under the Apache License, Version 2.0 (the License);
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an AS IS BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package ints;

import java.util.Arrays;
import java.util.stream.IntStream;

/**
 * <p>Class {@code WrappedIntArray} represents an immutable
 * {@code int[]} array.
 * </p>
 * Instances of {@code WrappedIntArray} are immutable.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class WrappedIntArray implements IntArray {

    private final int[] ia;

    /**
     * Constructs a new {@code WrappedIntArray} instance.
     * @param ia an array of integers
     * @throws NullPointerException if {@code ia == null}
     */
    public WrappedIntArray(int[] ia) {
        this.ia = ia.clone();
    }

    /**
     * Constructs an {@code WrappedIntrray} object with the specified
     * values.
     * @param values a stream of integer values
     * @throws NullPointerException if {@code values == null}
     * @throws IllegalStateException if the stream has previously been used
     */
    public WrappedIntArray(IntStream values) {
        this.ia = values.toArray();
    }

    /**
     * Constructs a new {@code WrappedIntArray} instance.
     * @param ia an array of integers
     * @param valueSize the exclusive end of the range of non-negative
     * array values
     * @throws IllegalArgumentException if
     * {@code (ia[j] < 0 || ia[j] > valueSize)} for any index {@code j}
     * satisfying  {@code (j >= 0 && j < ia.length)}
     * @throws NullPointerException if {@code ia == null}
     */
    public WrappedIntArray(int[] ia, int valueSize) {
        this.ia = new int[ia.length];
        for (int j=0; j<ia.length; ++j) {
            if (ia[j]<0 || ia[j]>=valueSize) {
                throw new IllegalArgumentException(String.valueOf(ia[j]));
            }
            this.ia[j] = ia[j];
        }
    }

    /**
     * Constructs a new {@code WrappedIntArray} instance.
     * @param il a list of integers
     * @throws NullPointerException if {@code il == null}
     */
    public WrappedIntArray(IntList il) {
        this.ia = il.toArray();
    }

    /**
     * Constructs a new {@code WrappedIntArray} instance.
     * @param il a list of integers
     * @param valueSize the exclusive end of the range of non-negative
     * array values
     * @throws IllegalArgumentException if
     * {@code (il[j] < 0 || il[j] > valueSize)} for any index {@code j}
     * satisfying  {@code (j >= 0 && j < il.length)}
     * @throws NullPointerException if {@code il == null}
     */
    public WrappedIntArray(IntList il, int valueSize) {
        this.ia = new int[il.size()];
        for (int j=0; j<ia.length; ++j) {
            int value = il.get(j);
            if (value<0 || value>=valueSize) {
                throw new IllegalArgumentException(String.valueOf(value));
            }
            ia[j] = value;
        }
    }

    @Override
    public int size() {
        return ia.length;
    }

    @Override
    public int get(int index) {
        return ia[index];
    }

    @Override
    public String toString() {
        return Arrays.toString(ia);
    }

    /**
     * Returns the list of integers.
     * @return the list of integers
     */
    public int[] toArray() {
        return ia.clone();
    }

     /**
     * Searches {@code this} for the specified value using the binary search
     * algorithm. This list must be sorted (as by the
     * {@code java.util.Arrays.sort(int[])} method) prior to making
     * this call. If it is not sorted, the results are undefined.
     * If the list contains multiple elements with the specified
     * value, there is no guarantee which one will be found.
     *
     * @param key the value to be searched for
     *
     * @return index of the search key, if it is contained in the list;
     * otherwise, {@code (-(insertion point) - 1)}. The insertion point is
     * defined as the point at which the key would be inserted into the list:
     * the index of the first element greater than the key, or
     * {@code this.size()} if all elements in the list are less than the
     * specified key. Note that this guarantees that the return value will
     * be {@code >= 0} if and only if the key is found.
     */
    public int binarySearch(int key) {
        return Arrays.binarySearch(ia, key);
    }

    /**
     * Searches the specified range of {@code this} for the specified value
     * using the binary search algorithm. This range must be sorted (as by the
     * {@code java.util.Arrays.sort(int[])} method) prior to making
     * this call. If it is not sorted, the results are undefined.
     * If the range contains multiple elements with the specified
     * value, there is no guarantee which one will be found. This method
     * considers all NaN values to be equivalent and equal.
     *
     * @param fromIndex the index of the first element (inclusive) to be searched
     * @param toIndex the index of the last element (exclusive) to be searched
     * @param key the value to be searched for
     *
     * @return index of the search key, if it is contained in the list;
     * otherwise, {@code (-(insertion point) - 1)}. The insertion point is
     * defined as the point at which the key would be inserted into the list:
     * the index of the first element greater than the key, or
     * {@code this.size()} if all elements in the list are less than the
     * specified key. Note that this guarantees that the return value will
     * be {@code >= 0} if and only if the key is found.
     *
     * @throws IllegalArgumentException if {@code fromIndex > toIndex}
     * @throws ArrayIndexOutOfBoundsException if
     * {@code fromIndex < 0 || toIndex > this.size()}
     */
    public int binarySearch(int fromIndex, int toIndex, int key) {
        return Arrays.binarySearch(ia, fromIndex, toIndex, key);
    }
}
