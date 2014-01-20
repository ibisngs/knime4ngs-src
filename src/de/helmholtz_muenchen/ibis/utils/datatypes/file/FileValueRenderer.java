package de.helmholtz_muenchen.ibis.utils.datatypes.file;


import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.renderer.AbstractDataValueRendererFactory;
import org.knime.core.data.renderer.DataValueRenderer;
import org.knime.core.data.renderer.MultiLineStringValueRenderer;

@SuppressWarnings("serial")
public class FileValueRenderer extends MultiLineStringValueRenderer {
	   /**
     * Factory for {@link CMLValueRenderer}.
     */
    public static final class Factory extends AbstractDataValueRendererFactory {
        /**
         * {@inheritDoc}
         */
        @Override
        public String getDescription() {
            return "File String";
        }

        /**
         * {@inheritDoc}
         */
        @Override
        public DataValueRenderer createRenderer(final DataColumnSpec colSpec) {
            return new FileValueRenderer(getDescription());
        }
    }

    /**
     * Constructor.
     *
     * @param description a description for the renderer
     */
    FileValueRenderer(final String description) {
        super(description);
    }

    /** {@inheritDoc} */
    @Override
    protected void setValue(final Object value) {
        if (value instanceof FileValue) {
            super.setValue(((FileValue)value).getFilePath());
        } else {
            super.setValue(value);
        }
    }
}
