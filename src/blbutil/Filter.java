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

import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

/**
 * <p>A filter for accepting or rejecting objects.
 * </p>
 * Instances of class {@code Filter} are required to be immutable.
 *
 * @param <E> the type of object that is filtered.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public interface Filter<E> {

    /**
     * Returns a filter that accepts all non-null objects.
     * @param <E> the type of object that is filtered
     * @return a filter that accepts all non-null objects
     */
    static <E> Filter<E> acceptAllFilter() {
        return (E e) -> {
            if (e == null) {
                throw new NullPointerException("e==null");
            }
            return true;
        };
    }

    /**
     * Returns a filter that accepts all non-null objects that are
     * contained in the specified collection.
     * @param <E> the type of object that is filtered
     * @param include the collection of objects that will be accepted by
     * the filter
     * @return a filter that accepts all non-null objects that are
     * contained in the specified collection
     * @throws NullPointerException if {@code include == null}
     */
    static <E> Filter<E> includeFilter(Collection<E> include) {
        final Set<E> includeSet = new HashSet<>(include);
        return (E e) -> {
            if (e == null) {
                throw new NullPointerException("e==null");
            }
            return includeSet.contains(e);
        };
    }

    /**
     * Returns a filter that accepts all non-null objects that are not
     * contained in the specified collection.
     * @param <E> the type of object that is filtered
     * @param exclude the collection of objects that will be rejected
     * by the filter
     * @return a filter that accepts all non-null objects that are not
     * contained in the specified collection
     * @throws NullPointerException if {@code exclude == null}
     */
    static <E> Filter<E> excludeFilter(Collection<E> exclude) {
        final Set<E> includeSet = new HashSet<>(exclude);
        return (E e) -> {
            if (e == null) {
                throw new NullPointerException("e==null");
            }
            return !includeSet.contains(e);
        };
    }

    /**
     * Returns {@code true} if the specified object is
     * accepted and returns {@code false} if the specified object
     * is rejected.
     * @param e the object to be filtered
     * @return {@code true} if the specified object is
     * accepted
     * @throws NullPointerException if {@code e==null}
     */
    boolean accept(E e);
}
