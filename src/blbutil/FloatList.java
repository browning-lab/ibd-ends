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
package blbutil;

import java.util.Arrays;

/**
 * Class {@code FloatList} represents a list of floats.
 * Class {@code FloatList} supports a {@code clear()} method, but does not
 * support a {@code remove()} method.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class FloatList {

    /**
     * The default initial capacity of an {@code FloatList}, which is 10.
     */
    public static final int DEFAULT_INIT_CAPACITY = 10;

    private int size;
    private float[] values;

    /**
     * Constructs an {@code FloatList} object with the default
     * initial capacity.
     *
     * @see #DEFAULT_INIT_CAPACITY
     */
    public FloatList() {
        this(DEFAULT_INIT_CAPACITY);
    }

    /**
     * Constructs an {@code FloatList} object with the specified
     * initial capacity.
     *
     * @param initCapacity the initial capacity of this list
     * @throws IllegalArgumentException if {@code initCapacity<0}.
     */
    public FloatList(int initCapacity) {
        if (initCapacity < 0) {
            throw new IllegalArgumentException(String.valueOf(initCapacity));
        }
        this.size = 0;
        this.values = new float[initCapacity];
    }

    /**
     * Adds the specified integer to the end of this list.
     *
     * @param element the value to be added to the end of this list.
     */
    public void add(float element) {
        if (size==values.length) {
            int newCapacity = (values.length * 3)/2 + 1;
            this.values = Arrays.copyOf(this.values, newCapacity);
        }
        this.values[size++] = element;
    }

    /**
     * Adds the specified value to the specified element.
     *
     * @param index the index of the element to which the specified value
     * will be added
     * @param value the to be added
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 || index >= this.size()}
     */
    public void addToElement(int index, float value) {
        if (index >= size) {
            throw new IndexOutOfBoundsException(String.valueOf(index));
        }
        this.values[index] += value;
    }

    /**
     * Returns the float at the specified position in this list.
     * @param index the index of the returned float.
     * @return the float at the specified position in this list.
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 || index >= this.size()}
     */
    public float get(int index) {
        if (index >= size) {
            throw new IndexOutOfBoundsException(String.valueOf(index));
        }
        return values[index];
    }

    /**
     * Replaces the element at the specified position in this list with the
     * specified element.
     * @param index the index of the element to be replaced
     * @param value the value to be stored at the specified position
     * in this list
     * @return the previous value at the specified position in this list
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 || index >= this.size()}
     */
    public float set(int index, float value) {
        if (index >= size) {
            throw new IndexOutOfBoundsException(String.valueOf(index));
        }
        float oldValue = values[index];
        values[index] = value;
        return oldValue;
    }

    /**
     * Returns the number of elements in this list.
     * @return the number of elements in this list.
     */
    public int size() {
        return size;
    }

    /**
     * Returns {@code true} if this list has no elements, and returns
     * {@code false} otherwise.
     * @return {@code true} if this list has no elements, and returns
     * {@code false} otherwise.
     */
    public boolean isEmpty() {
        return size==0;
    }

    /**
     * Returns an integer array containing the sequence of elements in this
     * list.
     * @return an integer array containing the sequence of elements in this
     * list.
     */
    public float[] toArray() {
        return Arrays.copyOf(values, size);
    }

    /**
     * Removes all elements from this list.
     */
    public void clear() {
        this.size = 0;
    }

    /**
     * Returns a string representation of this list.  The
     * exact details of the representation are unspecified and
     * subject to change.
     *
     * @return a string representation of this list.
     */
    @Override
    public String toString() {
        return Arrays.toString(Arrays.copyOf(values, size));
    }
}
