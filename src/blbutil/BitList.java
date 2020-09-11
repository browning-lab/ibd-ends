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
 * <p>Interface {@code BitList} represents a mutable sequence of bits
 * with a fixed length.</p>
 *
 * <p>Instances of {@code BitList} are not thread-safe.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class BitList {

    private static final int LOG2_BITS_PER_WORD = 6;
    private static final long WORD_MASK = 0xffffffffffffffffL;
    private static final int BIT_INDEX_MASK = (1 << LOG2_BITS_PER_WORD) - 1;

    private final long[] words;
    private final int size;

    /**
     * Constructs a {@code BitList} instance with the specified {@code size}
     * and having all bits set to 0 (unset).
     * @param size the number of bits
     * @throws IllegalArgumentException if {@code size < 0}
     */
    public BitList(int size) {
        if (size<0) {
            throw new IllegalArgumentException(String.valueOf(size));
        }
        int nWords = (size + Long.SIZE - 1) / Long.SIZE;
        this.words = new long[nWords];
        this.size = size;
    }

    /**
     * Constructs a {@code BitList} instance from the specified values.
     * @param values a sequence of bits
     * @param size the number of bits
     * @throws IllegalArgumentException if {@code size < 0}
     * @throws IllegalArgumentException if
     * {@code values.length != (size + Long.SIZE - 1) / Long.SIZE}
     * @throws NullPointerException if {@code values == null}
     */
    public BitList(long[] values, int size) {
        if (size<0) {
            throw new IllegalArgumentException(String.valueOf(size));
        }
        int nWords = (size + Long.SIZE - 1) / Long.SIZE;
        if (values.length!=nWords) {
            throw new IllegalArgumentException(String.valueOf(values.length));
        }
        this.words = Arrays.copyOf(values, nWords);
        this.size = size;
    }

    /**
     * Constructs a new {@code BitList} instance with the same sequence of
     * bits as the specified {@code BitList}.
     * @param bitList a sequence of bits to be copied
     * @throws IllegalArgumentException if {@code size < 0}
     * @throws NullPointerException if {@code bitList == null}
     */
    public BitList(BitList bitList) {
        this.words = bitList.words.clone();
        this.size = bitList.size;
    }

    /**
     * Returns the size of this {@code BitList}.
     * @return the size of this {@code BitList}
     */
    public int size() {
        return size;
    }

    /**
     * Returns the specified bit as a {@code boolean} value.
     * @param index a bit index
     * @return the specified bit as a {@code boolean} value.
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 ||  index >= this.size()}
     */
    public boolean get(int index) {
        if (index>size) {
            throw new IndexOutOfBoundsException(String.valueOf(index));
        }
        int wordIndex =  index >> LOG2_BITS_PER_WORD;
        return (words[wordIndex] & (1L << index)) != 0L;
    }

    /**
     * Sets the specified bit.
     * @param index a bit index
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 ||  index >= this.size()}
     */
    public void set(int index) {
        if (index>size) {
            throw new IndexOutOfBoundsException(String.valueOf(index));
        }
        int wordIndex =  index >> LOG2_BITS_PER_WORD;
        words[wordIndex] |= (1L << index);
    }

    /**
     * Clears the specified bit.
     * @param index a bit index
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 ||  index >= this.size()}
     */
    public void clear(int index) {
        if (index>size) {
            throw new IndexOutOfBoundsException(String.valueOf(index));
        }
        int wordIndex =  index >> LOG2_BITS_PER_WORD;
        words[wordIndex] &= ~(1L << index);    // need to check this
    }

    /**
     * Clears all bits.
     */
    public void clear() {
        Arrays.fill(words, 0L);
    }

    /**
     * Returns a new {@code BitList} of size {@code (from - to)} that is a
     * copy of the specified bit indices of this {@code BitList}.
     * @param from the first bit to be copied (inclusive)
     * @param to the last bit to be copied (exclusive)
     * @return a new {@code BitList} of size {@code (from - to)} that is a
     * copy of the specified bit indices of this {@code BitList}
     * @throws IndexOutOfBoundsException if
     * {@code from < 0 || from > to || from >= this.size}
     * @throws NullPointerException if {@code source == null}
     */
    public BitList restrict(int from, int to) {
        if (from>to || from>=size) {
            throw new IndexOutOfBoundsException(String.valueOf(from));
        }
        if (from==to) {
            return  new BitList(0);
        }
        BitList result = new BitList(to-from);
        int nResultWordsM1 = result.words.length - 1;
        final boolean isWordAligned = ((from & BIT_INDEX_MASK) == 0);
        int srcWord = from >> LOG2_BITS_PER_WORD;
        for (int w=0; w<nResultWordsM1; w++, srcWord++) {
            result.words[w] = isWordAligned ? words[srcWord]
                    : (words[srcWord] >>> from) | (words[srcWord+1] << -from);
        }
        long endWordMask = WORD_MASK >>> -to;
        result.words[nResultWordsM1] =
                ((to-1) & BIT_INDEX_MASK) < (from & BIT_INDEX_MASK)
                ? ((words[srcWord] >>> from) | (words[srcWord+1] & endWordMask) << -from)
                : (words[srcWord] & endWordMask) >>> from;

        return result;
    }

    /**
     * Replaced the specified bits in this {@code Bitlist} with the corresponding
     * bits in the specified @code BitList}.
     * @param src the {@code BitList} to be copied from
     * @param from the first bit to be copied (inclusive)
     * @param to the last bit to be copied (exclusive)
     * @throws IndexOutOfBoundsException if
     * {@code from < 0 || from > to || from >= this.size || from >= source.size()}
     * @throws NullPointerException if {@code source == null}
     */
    public void copyFrom(BitList src, int from, int to) {
        if (from>to || from>=size || from >= src.size()) {
            throw new IndexOutOfBoundsException(String.valueOf(from));
        }
        if (from==to) {
            return;
        }
        int startWord = from >> LOG2_BITS_PER_WORD;
        int endWord = (to-1) >> LOG2_BITS_PER_WORD;
        long startWordMask = WORD_MASK << from;
        long endWordMask = WORD_MASK >>> -to;
        if (startWord==endWord) {
            long mask = (startWordMask & endWordMask);
            words[startWord] ^= ((words[startWord] ^ src.words[startWord]) & mask);
        }
        else {
            words[startWord] ^= ((words[startWord] ^ src.words[startWord]) & startWordMask);
            for (int j=startWord+1; j<endWord; ++j) {
                words[j] = src.words[j];
            }
            words[endWord] ^= ((words[endWord] ^ src.words[endWord]) & endWordMask);
        }
    }

    /**
     * Returns a hash code for the specified bits in this {@code Bitlist}
     * @param from the first bit (inclusive)
     * @param to the last bit (exclusive)
     * @return a hash code for the specified bits in this {@code Bitlist}
     * @throws IndexOutOfBoundsException if
     * {@code from < 0 || from > to || from >= this.size}
     */
    public int hash(int from, int to) {
        if (from<0 || from>to || from>=size) {
            throw new IndexOutOfBoundsException(String.valueOf(from));
        }
        if (from==to) {
            return 0;
        }
        int startWord = from >> LOG2_BITS_PER_WORD;
        int endWord = (to-1) >> LOG2_BITS_PER_WORD;
        long startWordMask = WORD_MASK << from;
        long endWordMask = WORD_MASK >>> -to;
        if (startWord==endWord) {
            long mask = (startWordMask & endWordMask);
            return Long.hashCode(words[startWord] & mask);
        }
        else {
            long longHash = (words[startWord] & startWordMask);
            for (int j=startWord+1; j<endWord; ++j) {
                longHash ^= words[j];
            }
            longHash ^= (words[endWord] & endWordMask);
            return Long.hashCode(longHash);
        }
    }

    /**
     * Swaps the specified bits of the two specified {@code Bitlist} objects.
     * @param a the first {@code BitList}
     * @param b the second {@code BitList}
     * @param from the first bit to be copied (inclusive)
     * @param to the last bit to be copied (exclusive)
     * @throws IllegalArgumentException if
     * {@code s.size() != b.size()}
     * @throws IndexOutOfBoundsException if
     * {@code from < 0 || from > to || from >= a.size()}
     * @throws NullPointerException if {@code a == null || b == null}
     */
    public static void swapBits(BitList a, BitList b, int from, int to) {
        if (from>to || from >= a.size()) {
            throw new IndexOutOfBoundsException(String.valueOf(from));
        }
        if (a.size()!=b.size()) {
            throw new IllegalArgumentException("inconsistent data");
        }
        if (from==to) {
            return;
        }
        int startWord = from >> LOG2_BITS_PER_WORD;
        int endWord = (to-1) >> LOG2_BITS_PER_WORD;
        long startWordMask = WORD_MASK << from;
        long endWordMask = WORD_MASK >>> -to;
        if (startWord==endWord) {
            long mask = (startWordMask & endWordMask);
            a.words[startWord] ^= (b.words[startWord] & mask);
            b.words[startWord] ^= (a.words[startWord] & mask);
            a.words[startWord] ^= (b.words[startWord] & mask);
        }
        else {
            a.words[startWord] ^= (b.words[startWord] & startWordMask);
            b.words[startWord] ^= (a.words[startWord] & startWordMask);
            a.words[startWord] ^= (b.words[startWord] & startWordMask);
            for (int j=startWord+1; j<endWord; ++j) {
                a.words[j] ^= b.words[j];
                b.words[j] ^= a.words[j];
                a.words[j] ^= b.words[j];
            }
            a.words[endWord] ^= (b.words[endWord] & endWordMask);
            b.words[endWord] ^= (a.words[endWord] & endWordMask);
            a.words[endWord] ^= (b.words[endWord] & endWordMask);
        }
    }

    /**
     * Returns {@code true} if this {@code Bitlist} and the specified
     * {@code BitList} have identical sequences of bits for the specified
     * indices, and returns {@code false} otherwise.
     * @param other the {@code BitList} to be compared with {@code this} for
     * equality.
     * @param from the first bit to be compared (inclusive)
     * @param to the last bit to be compared (exclusive)
     * @return {@code true} if this {@code Bitlist} and the specified
     * {@code BitList} have identical sequences of bits for the specified
     * indices.
     * @throws IndexOutOfBoundsException if
     * {@code from < 0 || from > to || from >= this.size || from >= other.size()}
     * @throws NullPointerException if {@code other == null}
     */
    public boolean equal(BitList other, int from, int to) {
        if (from>to || from>=size || from>=other.size()) {
            throw new IndexOutOfBoundsException(String.valueOf(from));
        }
        int startWord = from >> LOG2_BITS_PER_WORD;
        int endWord = (to-1) >> LOG2_BITS_PER_WORD;
        long startWordMask = WORD_MASK << from;
        long endWordMask = WORD_MASK >>> -to;
        if (startWord==endWord) {
            long mask = (startWordMask & endWordMask);
            return ((words[startWord] ^ other.words[startWord]) & mask) == 0L;
        }
        else {
            boolean areEqual = true;
            areEqual &= (((words[startWord] ^ other.words[startWord]) & startWordMask) == 0L);
            for (int j=startWord+1; j<endWord; ++j) {
                areEqual &= words[j]==other.words[j];
            }
            areEqual &= (((words[endWord] ^ other.words[endWord]) & endWordMask) == 0L);

            return areEqual;
        }
    }

    /**
     * Returns this {@code BitList} as a {@code long} array.
     * @return this {@code BitList} as a {@code long} array
     */
    public long[] toLongArray() {
        return words.clone();
    }

    /**
     * Returns a string representation of this {@code BitList}.
     *
     * @return a string representation of this {@code BitList}.
     */
    public String asString() {
        StringBuilder sb = new StringBuilder(size);
        for (int j=0; j<size; ++j) {
            sb.append(get(j) ? '1' : '0');
        }
        return sb.toString();
    }

    /**
     * Returns {@code true} if the specified {@code BitList} objects
     * represent the same sequence of long values, and returns {@code false}
     * otherwise.
     * @param a a sequence of long values
     * @param b a sequence of long values
     * @return {@code true} if the specified {@code BitList} objects
     * represent the same sequence of long values
     */
    public static boolean equals(BitList a, BitList b) {
        return Arrays.equals(a.words, b.words);
    }
}
