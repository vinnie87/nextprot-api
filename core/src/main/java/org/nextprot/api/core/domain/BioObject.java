package org.nextprot.api.core.domain;

import com.google.common.base.Preconditions;

import java.io.Serializable;

/**
 * A wrapper over biological domain object
 *
 * Created by fnikitin on 26/08/15.
 */
public abstract class BioObject<T> implements Serializable {

    private static final long serialVersionUID = 1L;

    protected static final String NEXTPROT = "neXtProt";

    public enum BioType { CHEMICAL, PROTEIN, PROTEIN_ISOFORM, COMPLEX, GROUP}
    public enum ResourceType { INTERNAL, EXTERNAL, MIXED }

    private long id;
    private String accession;
    private String database;
    private final BioType bioType;
    private final ResourceType resourceType;
    transient private T content;

    protected BioObject(BioType bioType, ResourceType resourceType) {

        Preconditions.checkNotNull(bioType);

        this.bioType = bioType;
        this.resourceType = resourceType;
    }

    public long getId() {
        return id;
    }

    public void setId(long id) {
        this.id = id;
    }

    public String getAccession() {
        return accession;
    }

    public void setAccession(String accession) {
        this.accession = accession;
    }

    public String getDatabase() {
        return database;
    }

    public void setDatabase(String database) {
        this.database = database;
    }

    public BioType getBioType() {
        return bioType;
    }

    public ResourceType getResourceType() {
        return resourceType;
    }

    public int size() {
        return 1;
    }

    public T getContent() {
        return content;
    }

    public void setContent(T content) {
        this.content = content;
    }
}
