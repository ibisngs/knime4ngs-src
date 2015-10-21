package de.helmholtz_muenchen.ibis.utils.datatypes.file;

import org.knime.core.data.DataValue;

public interface FileValue extends DataValue {

    /**
     * @return A String value.
     */
    String getFilePath();

}
