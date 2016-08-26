package org.nextprot.api.commons.bio.variation.impl;

import org.nextprot.api.commons.bio.variation.SequenceChange;

/**
 * A simple deletion with no associated value.
 *
 * Created by fnikitin on 10/07/15.
 */
public class Deletion implements SequenceChange<Object> {

    private static Deletion INSTANCE = new Deletion();

    private Deletion() { }

    public static Deletion getInstance() {

        return INSTANCE;
    }

    /**
     * No value is associated with a deletion
     * @return null
     */
    @Override
    public Object getValue() {
        return null;
    }

    @Override
    public Type getType() {

        return Type.DELETION;
    }
}