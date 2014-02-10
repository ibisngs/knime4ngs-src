package de.helmholtz_muenchen.ibis.utils.datatypes.file;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.eclipse.jface.preference.IPreferenceStore;
import org.eclipse.jface.util.IPropertyChangeListener;
import org.eclipse.jface.util.PropertyChangeEvent;
import org.eclipse.ui.plugin.AbstractUIPlugin;
import org.knime.core.data.DataValue;
import org.osgi.framework.BundleContext;

public class DataTypesNodePlugin extends AbstractUIPlugin{
    // The shared instance.
    private static DataTypesNodePlugin plugin;

    private List<DataTypeTuple> m_typesList;

    /**
     * The constructor.
     */
    public DataTypesNodePlugin() {
        super();
        plugin = this;
    }

    /**
     * This method is called upon plug-in activation.
     *
     * @param context the OSGI bundle context
     * @throws Exception if this plugin could not be started
     */
    @Override
    public void start(final BundleContext context) throws Exception {
        super.start(context);
        final IPreferenceStore pStore = getPreferenceStore();
        pStore.addPropertyChangeListener(new IPropertyChangeListener() {
           /** {@inheritDoc} */
            @SuppressWarnings("deprecation")
			@Override
            public void propertyChange(final PropertyChangeEvent event) {
                for (DataTypeTuple t : getCustomizableTypeList()) {
                    if (t.getPreferenceKey().equals(event.getProperty())) {
                        Object n = event.getNewValue();
                        String newValue = n != null && n.toString().length() > 0
                            ? n.toString() : null;
                        t.getUtilityFactory().setPreferredRendererClassName(
                                newValue);
                    }
                }
            }
        });
    }

    /**
     * This method is called when the plug-in is stopped.
     *
     * @param context the OSGI bundle context
     * @throws Exception if this plugin could not be stopped
     */
    @Override
    public void stop(final BundleContext context) throws Exception {
        super.stop(context);
        plugin = null;
    }

    /**
     * Lookups the preferred renderer for a given molecular type.
     *
     * @param valueClass The DataValue class of interest, one of {@link SdfValue}, {@link Mol2Value},
     *            {@link SmilesValue}
     * @return The class name of the stored preferred renderer or <code>null</code> if none has been saved.
     *
     * @deprecated the preferred method now is to make use of the {@link ExtensibleUtilityFactory} of the specific
     *             {@link DataValue}s
     */
    @Deprecated
    public String getPreferredRendererClassName(
            final Class<? extends DataValue> valueClass) {
        final IPreferenceStore pStore = getPreferenceStore();
        String preferenceIdentifier = getPreferenceIdentifier(valueClass);
        String resultClassName = pStore.getString(preferenceIdentifier);
        if (resultClassName == null || resultClassName.length() == 0) {
            return null;
        }
        return resultClassName;
    }

    /** Get the list of types which can be customized in the preference page.
     * This currently contains tuples for SDF, Smiles, and Mol2.
     * @return An unmodifiable list of such types.
     *
     * @deprecated every data value type is now customizable
     */
    @Deprecated
    public List<DataTypeTuple> getCustomizableTypeList() {
        if (m_typesList == null) {
            m_typesList = new ArrayList<DataTypeTuple>();
          //  m_typesList.add(new DataTypeTuple("File", FileValue.UTILITY,
          //          getPreferenceIdentifier(FileValue.class)));
        }
        return Collections.unmodifiableList(m_typesList);
    }

    /**
     * Returns the shared instance.
     *
     * @return singleton instance of the Plugin
     */
    public static DataTypesNodePlugin getDefault() {
        return plugin;
    }

    /** Get preference name for a given data value class. */
    private static String getPreferenceIdentifier(
            final Class<? extends DataValue> valueClass) {
        return "prefRenderer_" + valueClass.getName();
    }

    /** Utility class that summarizes different information to a chemical type.
     * This includes the label (as shown in the preference page), the
     * accompanying utility class and the key for the preference store.
     *
     * @deprecated not needed any more
     */
    @Deprecated
    public static final class DataTypeTuple {
        private final String m_label;
        private final MyUtilityFactory m_utilityFactory;
        private final String m_preferenceKey;

        /**
         *
         */
        private DataTypeTuple(final String label,
                final MyUtilityFactory utilityFactory,
                final String preferenceKey) {
            m_label = label;
            m_utilityFactory = utilityFactory;
            m_preferenceKey = preferenceKey;
        }

        /** @return the label */
        public String getLabel() {
            return m_label;
        }

        /** @return the utilityFactory */
        public MyUtilityFactory getUtilityFactory() {
            return m_utilityFactory;
        }

        /** @return the preferenceKey */
        public String getPreferenceKey() {
            return m_preferenceKey;
        }
    }
}
