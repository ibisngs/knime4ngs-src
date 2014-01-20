package de.helmholtz_muenchen.ibis.utils.datatypes.file;

import java.util.ArrayList;
import java.util.List;

import javax.swing.Icon;
import javax.swing.ImageIcon;

import org.eclipse.core.runtime.preferences.IEclipsePreferences;
import org.eclipse.core.runtime.preferences.InstanceScope;
import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DataValue;
import org.knime.core.data.DataValueComparator;
import org.knime.core.data.ExtensibleUtilityFactory;
import org.knime.core.data.StringValueComparator;
import org.knime.core.data.renderer.DataValueRenderer;
import org.knime.core.data.renderer.DataValueRendererFactory;
import org.knime.core.data.renderer.StringValueRenderer;
import org.osgi.framework.FrameworkUtil;

public class MyUtilityFactory extends ExtensibleUtilityFactory{
	   private boolean m_preferencesRead;


	    /**
	     * Converts a renderer instance into a factory that always returns the instance. This is only for compatibility
	     * reasons and can cause all kinds of problems because the renderers are not thread safe.
	     */
	    private static final class RendererFactoryAdapter implements DataValueRendererFactory {
	        private final DataValueRenderer m_renderer;

	        RendererFactoryAdapter(final DataValueRenderer renderer) {
	            m_renderer = renderer;
	        }

	        /**
	         * {@inheritDoc}
	         */
	        @Override
	        public String getDescription() {
	            return m_renderer.getDescription();
	        }

	        /**
	         * {@inheritDoc}
	         */
	        @Override
	        public String getId() {
	            return m_renderer.getClass().getName();
	        }

	        /**
	         * {@inheritDoc}
	         */
	        @Override
	        public DataValueRenderer createRenderer(final DataColumnSpec colSpec) {
	            return m_renderer;
	        }
	    }

	    /** Singleton icon to be used to display this cell type. */
	    private final Icon m_icon;

	    /**
	     * If you need to the functionality of this class, your class either gets one of the KNIME core chem types or your
	     * need to copy code.
	     *
	     * @param iconName the name of the icon
	     * @param valueClass the value class this class is the utility object for
	     */
	    protected MyUtilityFactory(final String iconName, final Class<? extends DataValue> valueClass) {
	        super(valueClass);
	        ImageIcon icon;
	        try {
	            ClassLoader loader = getClass().getClassLoader();
	            String path = getClass().getPackage().getName().replace('.', '/');
	            icon = new ImageIcon(loader.getResource(path + iconName));
	        } catch (Exception e) {
	            m_logger.info("Could not find icon for " + valueClass.getName(), e);
	            icon = null;
	        }
	        m_icon = icon;
	    }

	    /**
	     * {@inheritDoc}
	     */
	    @Override
	    public Icon getIcon() {
	        return m_icon;
	    }

	    private static final StringValueComparator STRING_COMPARATOR = new StringValueComparator();

	    /**
	     * {@inheritDoc}
	     */
	    @Override
	    protected DataValueComparator getComparator() {
	        return STRING_COMPARATOR;
	    }

	    /**
	     * Give others the possibility to register new renders for this value class. Note: This renderer instance may be
	     * used in different columns (or tables) at the same time. It must not rely on column specific information.
	     *
	     * <p>
	     * This method calls {@link #addRenderer(DataValueRenderer, boolean)}, whereby the flag is set to <code>true</code>
	     * (the renderer is used as default if possible).
	     *
	     * @param renderer a new renderer instance for this value class
	     * @throws IllegalArgumentException if argument is <code>null</code>
	     * @deprecated register renderers via the extension point
	     */
	    @Deprecated
	    public void addRenderer(final DataValueRenderer renderer) {
	        addRenderer(renderer, true);
	    }

	    /**
	     * Give others the possibility to register new renders for this value class. Note: This renderer instance may be
	     * used in different columns (or tables) at the same time. It must not rely on column specific information.
	     *
	     * @param renderer a new renderer instance for this value class
	     * @param suggestAsDefault Whether to use the renderer as default renderer for new workspaces (that may or may not
	     *            succeed - depending on other plug-ins also calling this method).
	     * @throws IllegalArgumentException if <code>renderer</code> is <code>null</code>
	     * @deprecated register renderers via the extension point
	     */
	    @Deprecated
	    public void addRenderer(final DataValueRenderer renderer, final boolean suggestAsDefault) {
	        if (renderer == null) {
	            throw new IllegalArgumentException("Renderer must not be null");
	        }
	        addRendererFactory(new RendererFactoryAdapter(renderer), suggestAsDefault);
	    }

	    /**
	     * Get list of all registered renders. This list is not ordered (that is, the preferred renderer as set in the
	     * preference page is not necessarily at first position.
	     *
	     * @return a list of available renderers
	     *
	     * @deprecated use {@link #getAvailableRenderers()} instead
	     */
	    @Deprecated
	    public List<DataValueRenderer> getAllAvailableRenderer() {
	        List<DataValueRenderer> renderers = new ArrayList<DataValueRenderer>();
	        for (DataValueRendererFactory fac : getAvailableRenderers()) {
	            try {
	                renderers.add(fac.createRenderer(null));
	            } catch (NullPointerException ex) {
	                // seems like the renderer does not like null specs
	                m_logger.debug("Could not create renderer from " + fac.getClass().getName() + " without column spec",
	                    ex);
	            }
	        }

	        return renderers;
	    }

	    /**
	     * Set the class name of the preferred renderer or <code>null</code> if to unset any preference ordering, which has
	     * been set previously. If the argument is invalid, this method does nothing except for printing a warning to the
	     * log.
	     *
	     * @param rendererClassName the class name of the renderer or <code>null</code>
	     * @deprecated preferred renderers should only be set by the user via the preference page
	     */
	    @Deprecated
	    public void setPreferredRendererClassName(final String rendererClassName) {
	        if (isPreferredRendererClassNameValid(rendererClassName)) {
	            if (rendererClassName == null) {
	                setPreferredRenderer((String)null);
	            } else {
	                for (DataValueRendererFactory fac : getAvailableRenderers()) {
	                    try {
	                        DataValueRenderer r = fac.createRenderer(null);
	                        if (r.getClass().getName().equals(rendererClassName)) {
	                            setPreferredRenderer(fac.getId());
	                            break;
	                        }
	                    } catch (NullPointerException ex) {
	                        // seems like the renderer does not like null specs
	                        m_logger.debug("Could not create renderer from " + fac.getClass().getName()
	                            + " without column spec", ex);
	                    }
	                }
	            }
	        } else {
	            m_logger.warn("Unknown preferred renderer class " + rendererClassName + " - ignoring");
	        }
	    }

	    /**
	     * @return class name of preferred renderer, <code>null</code> if none is set
	     *
	     * @deprecated use {@link #getPreferredRenderer()} instead.
	     */
	    @Deprecated
	    public String getPreferredRendererClassName() {
	        try {
	            return getPreferredRenderer().createRenderer(null).getClass().getName();
	        } catch (NullPointerException ex) {
	            // seems like the renderer does not like null specs
	            m_logger.debug("Could not get class name of preferred renderer for " + m_valueClass, ex);
	            return null;
	        }
	    }

	    /**
	     * Method to override and return a default "fall back" renderer. This default implementation returns
	     * a {@link StringValueRenderer} but derived classes may return a different renderer here (for instance
	     * molecular types such as, e.g. Mol2, that have line breaks in their string representation return a
	     * {@link org.knime.core.data.renderer.MultiLineStringValueRenderer} here.
	     *
	     * @return a fall back renderer, should not be <code>null</code>
	     * @deprecated register at least one renderer via the extension point instead
	     */
	    @Deprecated
	    protected DataValueRenderer getFallBackRenderer() {
	        return new StringValueRenderer();
	    }

	    private boolean isPreferredRendererClassNameValid(final String className) {
	        if ((className == null) || (className.length() == 0)) {
	            return true;
	        }
	        for (DataValueRenderer r : getAllAvailableRenderer()) {
	            if (className.equals(r.getClass().getName())) {
	                return true;
	            }
	        }
	        return false;
	    }

	    /**
	     * {@inheritDoc}
	     */
	    @Override
	    public String getName() {
	        return m_valueClass.getSimpleName();
	    }

	    /**
	     * {@inheritDoc}
	     */
	    @Override
	    protected void readPreferredRendererFromPreferences() {
	        if (m_preferencesRead) {
	            return;
	        }

	        // only for backwards compatibility, it reads the old chem-only preferred renderers
	        String rendererClassName = DataTypesNodePlugin.getDefault().getPreferredRendererClassName(m_valueClass);
	        if (rendererClassName != null) {
	            if (rendererClassName.matches("^org\\.knime\\.chem\\.types\\..+UtilityFactory\\$\\d$")) {
	                // the old fallback-renderer, replace it with the default renderer
	                setPreferredRenderer(getDefaultRenderer().getId());
	            } else {
	                setPreferredRendererClassName(rendererClassName);
	            }
	            m_preferencesRead = true;

	            if (getPreferredRenderer() != null) {
	                // write to new preferences
	                IEclipsePreferences corePrefs =
	                    InstanceScope.INSTANCE.getNode(FrameworkUtil.getBundle(DataValueRendererFactory.class)
	                        .getSymbolicName());
	                corePrefs.put(ExtensibleUtilityFactory.getPreferenceKey(m_valueClass), getPreferredRenderer().getId());
	            }
	        }

	        // new style preferred renderers override old-style
	        super.readPreferredRendererFromPreferences();
	    }

	    /**
	     * {@inheritDoc}
	     */
	    @Override
	    public String getGroupName() {
	        return "Basic";
	    }
}
