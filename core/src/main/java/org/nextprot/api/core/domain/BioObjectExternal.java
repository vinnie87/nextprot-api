package org.nextprot.api.core.domain;

/**
 * A biological domain object external to nextProt
 *
 * Created by fnikitin on 26/08/15.
 */
public class BioObjectExternal extends BioObject<DbXref> {

    private static final long serialVersionUID = 0L;

    public BioObjectExternal(BioType bioType) {

        super(bioType, ResourceType.EXTERNAL);
    }
}