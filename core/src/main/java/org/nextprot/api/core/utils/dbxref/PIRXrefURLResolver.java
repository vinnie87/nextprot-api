package org.nextprot.api.core.utils.dbxref;

import org.nextprot.api.core.domain.DbXref;


class PirXrefURLResolver extends DbXrefURLBaseResolver {

    @Override
    protected String getAccessionNumber(DbXref xref) {

        DbXref.DbXrefProperty property = xref.getPropertyByName("entry name");

        if (property != null)
            return property.getValue();

        throw new UnresolvedXrefURLException("missing accession number: could not find 'entry name' property");
    }
}