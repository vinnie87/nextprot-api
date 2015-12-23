package org.nextprot.api.core.utils.dbxref.conv;

import org.nextprot.api.commons.constants.IdentifierOffset;
import org.nextprot.api.core.domain.CvDatabasePreferredLink;
import org.nextprot.api.core.domain.DbXref;

import java.util.ArrayList;
import java.util.List;

class EmblDbXrefConverter implements DbXrefPropertyToXrefConverter {

    @Override
    public List<DbXref> convert(DbXref xref) {

        List<DbXref> list = new ArrayList<>(3);

        for (DbXref.DbXrefProperty property : xref.getProperties()) {

            if ("genomic sequence ID".equals(property.getName())) {

                list.add(createEmblDbXrefFromDbXrefProperty(xref, property, CvDatabasePreferredLink.EMBL_GENE));
            }
            else if ("protein sequence ID".equals(property.getName())) {

                list.add(createEmblDbXrefFromDbXrefProperty(xref, property, CvDatabasePreferredLink.EMBL_PROTEIN));
            }
        }

        return list;
    }

    private DbXref createEmblDbXrefFromDbXrefProperty(DbXref xref, DbXref.DbXrefProperty property, CvDatabasePreferredLink preferredLink) {

        DbXref dbXRef = new DbXref();

        dbXRef.setDbXrefId(IdentifierOffset.XREF_PROPERTY_OFFSET +property.getPropertyId());
        dbXRef.setAccession(property.getValue());
        dbXRef.setDatabaseCategory("Sequence databases");
        dbXRef.setDatabaseName(CvDatabasePreferredLink.EMBL_GENE.getDbName());
        dbXRef.setUrl(xref.getUrl());
        dbXRef.setLinkUrl(preferredLink.getLink());
        dbXRef.setProperties(new ArrayList<DbXref.DbXrefProperty>());

        return dbXRef;
    }
}