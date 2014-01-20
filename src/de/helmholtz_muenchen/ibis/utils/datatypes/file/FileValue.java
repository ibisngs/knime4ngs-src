package de.helmholtz_muenchen.ibis.utils.datatypes.file;



import javax.swing.Icon;

import org.knime.core.data.DataValue;
import org.knime.core.data.DataValueComparator;
import org.knime.core.data.ExtensibleUtilityFactory;
import org.knime.core.data.StringValueComparator;

public interface FileValue extends DataValue {
    /**
     * Meta information to this value type.
     * @see DataValue#UTILITY
     */
	//FileUtilityFactory UTILITY = new FileUtilityFactory();

    /**
     * @return A String value.
     */
    String getFilePath();

   
    /** Implementations of the meta information of this value class. */
    class FileUtilityFactory extends ExtensibleUtilityFactory {
        /** Singleton icon to be used to display this cell type. */
        private static final Icon ICON =
            loadIcon(FileValue.class, "/file.png");

        private static final StringValueComparator STRING_COMPARATOR =
            new StringValueComparator();

        /** Only subclasses are allowed to instantiate this class. */
        protected FileUtilityFactory() {
            super(FileValue.class);
        }

        /**
         * {@inheritDoc}
         */
        @Override
        public Icon getIcon() {
            return ICON;
        }

        /**
         * {@inheritDoc}
         */
        @Override
        protected DataValueComparator getComparator() {
            return STRING_COMPARATOR;
        }

        /**
         * {@inheritDoc}
         */
        @Override
        public String getName() {
            return "File";
        }

        /**
         * {@inheritDoc}
         */
        @Override
        public String getGroupName() {
            return "Basic";
        }
    }
}
